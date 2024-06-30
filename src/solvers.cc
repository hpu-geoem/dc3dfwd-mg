#include <deal.II/base/quadrature_lib.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/multigrid.h>

#include "dc_problem.h"

/*
 * This function sets up the multigrid data structures.
*/
void DCProblem::setup_multigrid() {
  // Distribute the degrees of freedom on the multigrid levels
  dof_handler.distribute_mg_dofs();

  // Initialize the constrained dofs
  mg_constrained_dofs.clear();
  mg_constrained_dofs.initialize(dof_handler);

  const unsigned int n_levels = tria.n_global_levels();

  mg_matrix.resize(0, n_levels - 1);
  mg_matrix.clear_elements();
  mg_interface_in.resize(0, n_levels - 1);
  mg_interface_in.clear_elements();

  // Create the sparsity patterns for the matrices on each level
  for (unsigned int level = 0; level < n_levels; ++level) {
    const dealii::IndexSet dof_set = dealii::DoFTools::extract_locally_relevant_level_dofs(dof_handler, level);

    {
      dealii::DynamicSparsityPattern dsp(dof_set);
      dealii::MGTools::make_sparsity_pattern(dof_handler, dsp, level);
      dealii::SparsityTools::distribute_sparsity_pattern(dsp, dof_handler.locally_owned_mg_dofs(level), MPI_COMM_WORLD, dof_set);
      dsp.compress();

      mg_matrix[level].reinit(dof_handler.locally_owned_mg_dofs(level), dof_handler.locally_owned_mg_dofs(level), dsp,
                              MPI_COMM_WORLD);
    }

    {
      dealii::DynamicSparsityPattern dsp(dof_set);
      dealii::MGTools::make_interface_sparsity_pattern(dof_handler, mg_constrained_dofs, dsp, level);
      dealii::SparsityTools::distribute_sparsity_pattern(dsp, dof_handler.locally_owned_mg_dofs(level), MPI_COMM_WORLD, dof_set);
      dsp.compress();

      mg_interface_in[level].reinit(dof_handler.locally_owned_mg_dofs(level), dof_handler.locally_owned_mg_dofs(level), dsp,
                                    MPI_COMM_WORLD);
    }
  }
}

/*
 * This function assembles the multigrid matrices.
*/
void DCProblem::assemble_multigrid() {
  // Initialize the quadrature formula and the FEValues object
  dealii::QGauss<3> quadrature_formula(fe.degree + 1);
  dealii::FEValues<3> fe_values(fe, quadrature_formula,
                                dealii::update_values | dealii::update_gradients | dealii::update_quadrature_points |
                                    dealii::update_JxW_values);

  // Initialize the quadrature formula and the FEValues object for the faces
  dealii::QGauss<2> face_quadrature_formula(fe.degree + 1);
  dealii::FEFaceValues<3> fe_face_values(fe, face_quadrature_formula,
                                         dealii::update_values | dealii::update_quadrature_points |
                                             dealii::update_normal_vectors | dealii::update_JxW_values);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

  dealii::FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

  std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);

  // Initialize the boundary constraints
  std::vector<dealii::AffineConstraints<double>> boundary_constraints(tria.n_global_levels());
  for (unsigned int level = 0; level < tria.n_global_levels(); ++level) {
    const dealii::IndexSet dof_set = dealii::DoFTools::extract_locally_relevant_level_dofs(dof_handler, level);
    boundary_constraints[level].reinit(dof_set);
    boundary_constraints[level].add_lines(mg_constrained_dofs.get_refinement_edge_indices(level));
    boundary_constraints[level].close();
  }

  for (const auto &cell : dof_handler.cell_iterators()) {
    // Assemble the matrix only for the locally owned cells
    if (cell->level_subdomain_id() == tria.locally_owned_subdomain()) {
      cell_matrix = 0;
      fe_values.reinit(cell);

      // Compute the conductivity on the coarse cells
      // We use the volume weighted average of the conductivities of the active children cells of the coarse cell
      double sigma = 0.0, volume = 0.0;
      auto active_child_cells = dealii::GridTools::get_active_child_cells<dealii::DoFHandler<3>>(cell);
      if (active_child_cells.size() > 0) {
        for (auto active_cell : active_child_cells) {
          sigma += (1.0 / rho[active_cell->material_id()]) * active_cell->measure();
          volume += active_cell->measure();
        }
        sigma /= volume;
      } else {
        sigma = 1.0 / rho[cell->material_id()];
      }

      // Assemble the volume integral: \int_{\Omega} \sigma \nabla u \cdot \nabla v dx
      for (const unsigned int q_index : fe_values.quadrature_point_indices()) {
        for (const unsigned int i : fe_values.dof_indices())
          for (const unsigned int j : fe_values.dof_indices())
            cell_matrix(i, j) +=
                sigma * (fe_values.shape_grad(i, q_index) * fe_values.shape_grad(j, q_index) * fe_values.JxW(q_index));
      }
      cell->get_mg_dof_indices(local_dof_indices);

      // Assemble the face integrals
      for (unsigned int f : cell->face_indices()) {
        auto face = cell->face(f);
        if (!face->at_boundary()) {
          continue;
        }

        // Skip the surface boundary, i.e., faces with normal vector pointing upwards
        fe_face_values.reinit(cell, f);
        if (fe_face_values.get_normal_vectors()[0][2] < 0) {
          continue;
        }

        // Compute the face integral
        dealii::Tensor<1, 3> r = face->center() - source; // vector from the source to the face center
        for (unsigned int q_index : fe_face_values.quadrature_point_indices()) {
          dealii::Tensor<1, 3> n = fe_face_values.normal_vector(q_index);                       // outward normal vector
          double foo = sigma * dealii::scalar_product(r, n) / (r.norm() * n.norm() * r.norm()); // \sigma \frac{\cos(r, n)}{|r|}
          for (unsigned int i : fe_face_values.dof_indices()) {
            for (unsigned int j : fe_face_values.dof_indices()) {
              cell_matrix(i, j) += (foo * fe_face_values.shape_value(i, q_index) * fe_face_values.shape_value(j, q_index) *
                                    fe_face_values.JxW(q_index));
            }
          }
        }
      }
      cell->get_mg_dof_indices(local_dof_indices);

      boundary_constraints[cell->level()].distribute_local_to_global(cell_matrix, local_dof_indices, mg_matrix[cell->level()]);

      // Take care of the interface matrix
      for (unsigned int i = 0; i < dofs_per_cell; ++i) {
        for (unsigned int j = 0; j < dofs_per_cell; ++j) {
          if (mg_constrained_dofs.is_interface_matrix_entry(cell->level(), local_dof_indices[i], local_dof_indices[j])) {
            mg_interface_in[cell->level()].add(local_dof_indices[i], local_dof_indices[j], cell_matrix(i, j));
          }
        }
      }
    }
  }

  // Finalize the assembly
  for (unsigned int i = 0; i < tria.n_global_levels(); ++i) {
    mg_matrix[i].compress(dealii::VectorOperation::add);
    mg_interface_in[i].compress(dealii::VectorOperation::add);
  }
}

/*
 * This function solves the linear system using the selected solver.
*/
void DCProblem::solve() {
  if (settings.solver == Settings::amg) {
    solve_amg();
  } else if (settings.solver == Settings::gmg) {
    solve_gmg();
  }
}

/*
 * This function solves the linear system using the geometric multigrid preconditioner.
*/
void DCProblem::solve_gmg() {
  dealii::TimerOutput::Scope timing(timer, "Solve");

  // Setup the multigrid data structures
  timer.enter_subsection("Solve: setup multigrid");
  setup_multigrid();
  timer.leave_subsection("Solve: setup multigrid");

  // Assemble the multigrid matrices
  timer.enter_subsection("Solve: assemble multigrid");
  assemble_multigrid();
  timer.leave_subsection("Solve: assemble multigrid");

  timer.enter_subsection("Solve: setup preconditioner");

  // Setup the multigrid transfer operators
  dealii::MGTransferPrebuilt<PETScVector> mg_transfer(mg_constrained_dofs);
  mg_transfer.build(dof_handler);

  // Setup the coarse grid solver
  // We use CG solver and BoomerAMG as the preconditioner
  dealii::PETScWrappers::PreconditionBoomerAMG coarse_preconditioner;
  dealii::PETScWrappers::PreconditionBoomerAMG::AdditionalData coarse_amg_data;
  coarse_amg_data.symmetric_operator = true;
  coarse_amg_data.output_details = false;
  coarse_preconditioner.initialize(mg_matrix[0], coarse_amg_data);

  dealii::SolverControl coarse_solver_control(settings.max_iter);
  dealii::SolverCG<PETScVector> coarse_solver(coarse_solver_control);

  dealii::MGCoarseGridIterativeSolver<PETScVector, dealii::SolverCG<PETScVector>, PETScSparseMatrix,
                                      dealii::PETScWrappers::PreconditionBoomerAMG>
      coarse_grid_solver(coarse_solver, mg_matrix[0], coarse_preconditioner);

  // Setup the smoother: we use block Jacobi preconditioner
  using Smoother = dealii::PETScWrappers::PreconditionBlockJacobi;
  dealii::MGSmootherPrecondition<PETScSparseMatrix, Smoother, PETScVector> smoother;
  smoother.initialize(mg_matrix);
  smoother.set_steps(settings.smoother_steps);

  // Setup the multigrid object
  dealii::mg::Matrix<PETScVector> mg_m(mg_matrix);
  dealii::mg::Matrix<PETScVector> mg_in(mg_interface_in);
  dealii::mg::Matrix<PETScVector> mg_out(mg_interface_in);

  dealii::Multigrid<PETScVector> mg(mg_m, coarse_grid_solver, mg_transfer, smoother, smoother);
  mg.set_edge_matrices(mg_out, mg_in);

  // Setup the iterative solver
  dealii::SolverControl solver_control(settings.max_iter);
  dealii::SolverCG<PETScVector> solver(solver_control);
  dealii::PreconditionMG<3, PETScVector, dealii::MGTransferPrebuilt<PETScVector>> preconditioner(dof_handler, mg, mg_transfer);

  timer.leave_subsection("Solve: setup preconditioner");

  // Solve the primary problem
  {
    timer.enter_subsection("Solve: primary problem");
    solver_control.set_tolerance(settings.rtol * system_rhs.l2_norm());
    solver.solve(system_matrix, solution, system_rhs, preconditioner);
    pcout << "  Number of Iterations: " << solver_control.last_step() << std::endl;
    timer.leave_subsection("Solve: primary problem");
  }
  constraints.distribute(solution);

  // Solve the dual problem
  {
    timer.enter_subsection("Solve: dual problem");
    solver_control.set_tolerance(settings.rtol * rhs_dual.l2_norm());
    solver.solve(system_matrix, solution_dual, rhs_dual, preconditioner);
    pcout << "  Number of Iterations: " << solver_control.last_step() << std::endl;
    timer.leave_subsection("Solve: dual problem");
  }
  constraints.distribute(solution_dual);
}

/*
 * This function solves the linear system using the algebraic multigrid preconditioner.
*/
void DCProblem::solve_amg() {
  dealii::TimerOutput::Scope timing(timer, "Solve");

  dealii::SolverControl solver_control(settings.max_iter);
  dealii::SolverCG<PETScVector> solver(solver_control);
  dealii::PETScWrappers::PreconditionBoomerAMG preconditioner;
  dealii::PETScWrappers::PreconditionBoomerAMG::AdditionalData amg_data;

  amg_data.symmetric_operator = true;
  amg_data.output_details = false;

  // Setup the algebraic multigrid preconditioner
  {
    timer.enter_subsection("Solve: setup amg");
    preconditioner.initialize(system_matrix, amg_data);
    preconditioner.setup();
    timer.leave_subsection("Solve: setup amg");
  }

  // Solve the primary problem
  {
    timer.enter_subsection("Solve: primary problem");
    solver_control.set_tolerance(settings.rtol * system_rhs.l2_norm());
    solver.solve(system_matrix, solution, system_rhs, preconditioner);
    pcout << "  Number of Iterations: " << solver_control.last_step() << std::endl;
    timer.leave_subsection("Solve: primary problem");
  }
  constraints.distribute(solution);

  // Solve the dual problem
  {
    timer.enter_subsection("Solve: dual problem");
    solver_control.set_tolerance(settings.rtol * rhs_dual.l2_norm());
    solver.solve(system_matrix, solution_dual, rhs_dual, preconditioner);
    pcout << "  Number of Iterations: " << solver_control.last_step() << std::endl;
    timer.leave_subsection("Solve: dual problem");
  }
  constraints.distribute(solution_dual);
}
