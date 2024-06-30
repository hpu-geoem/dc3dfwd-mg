#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>

#include "dc_problem.h"
#include "settings.h"
#include "tria_finder.h"

/**
 * Initialize various members of the class, such as the finite element, the triangulation,
 * the DoF handler, the TriaFinder, the output stream, and the timer.
 * Note that the `limite_level_difference_at_vertices` and `construct_multigrid_hierarchy` flags
 * are required to enable the geometric multigrid solver.
 *
 * @param s The parameters for the DCProblem.
 */
DCProblem::DCProblem(Settings s)
    : settings(s), fe(settings.fe_degree), tria(MPI_COMM_WORLD, dealii::Triangulation<3>::limit_level_difference_at_vertices,
                                                true, dealii::parallel::shared::Triangulation<3>::construct_multigrid_hierarchy),
      dof_handler(tria), finder(tria), pcout(std::cout, dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0),
      timer(MPI_COMM_WORLD, pcout, dealii::TimerOutput::never, dealii::TimerOutput::wall_times) {}

/**
 * Reads the model and data files for the DCProblem.
 *
 * The model file contains the mesh and the resistivity values, while the data files contain
 * the source and the receiver locations.
 */
void DCProblem::read_model_and_data() {
  // The file name of the mesh file is the input prefix with the extension ".tria".
  std::ifstream ifs_tria(settings.iprefix + ".tria");

  // If the file does not exist, abort the program.
  if (!ifs_tria.good()) {
    abort();
  }
  // Otherwise, create a binary input archive
  boost::archive::binary_iarchive ia(ifs_tria);

  // Load the triangulation from the archive.
  tria.clear();
  tria.load(ia, 0);

  // Global refinement of the mesh.
  tria.refine_global(settings.n_global_refinements);

  // The file name of the resistivity file is the input prefix with the extension ".rho".
  // Open the file for reading.
  std::ifstream ifs_rho(settings.iprefix + ".rho");
  int nrhos;
  ifs_rho >> nrhos;
  for (int i = 0; i < nrhos; ++i) {
    double r;
    ifs_rho >> r;
    rho.push_back(r);
  }

  // The file name of the data file is the input prefix with the extension ".emd".
  std::ifstream ifs_emd(settings.iprefix + ".emd");

  // Read the source location from the file.
  double x, y, z;
  ifs_emd >> x >> y >> z;
  source = dealii::Point<3>(x, y, z);

  // Read the receiver locations from the file.
  int nsites;
  ifs_emd >> nsites;
  for (int i = 0; i < nsites; ++i) {
    ifs_emd >> x >> y >> z;
    sites.push_back(dealii::Point<3>(x, y, z));
  }
}

/**
 * Sets up the system for solving the DC problem.
 */
void DCProblem::setup_system() {
  dealii::TimerOutput::Scope timing(timer, "Setup system");

  // Distribute the degrees of freedom.
  dof_handler.distribute_dofs(fe);

  pcout << "  Number of active cells: " << tria.n_active_cells() << std::endl;
  pcout << "  Number of DoFs: " << dof_handler.n_dofs() << std::endl;

  // Extract the locally relevant and owned degrees of freedom.
  locally_relevant_dofs = dealii::DoFTools::extract_locally_relevant_dofs(dof_handler);
  locally_owned_dofs = dof_handler.locally_owned_dofs();

  // Initialize the constraints
  constraints.reinit(locally_relevant_dofs);
  dealii::DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  // Initialize the sparsity pattern of the system matrix.
  dealii::DynamicSparsityPattern dsp(locally_relevant_dofs);
  dealii::DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints);
  dealii::SparsityTools::distribute_sparsity_pattern(dsp, locally_owned_dofs, MPI_COMM_WORLD, locally_relevant_dofs);
  dsp.compress();

  // Initialize the system matrix, the solution, and the right-hand side vector.
  system_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, MPI_COMM_WORLD);

  solution.reinit(locally_owned_dofs, MPI_COMM_WORLD);
  system_rhs.reinit(locally_owned_dofs, MPI_COMM_WORLD);

  solution_dual.reinit(locally_owned_dofs, MPI_COMM_WORLD);
  rhs_dual.reinit(locally_owned_dofs, MPI_COMM_WORLD);
}

/**
 * Assembles the system matrix for the DC problem.
 *
 * This function calculates the element-wise contributions to the system matrix
 * by looping over all locally owned cells and faces, and distributes the local
 * contributions to the global system matrix.
 */
void DCProblem::assemble_system() {
  dealii::TimerOutput::Scope timing(timer, "Assemble system matrix");

  // Initialize the quadrature formula and the FEValues object
  dealii::QGauss<3> quadrature_formula(fe.degree + 1);
  dealii::FEValues<3> fe_values(fe, quadrature_formula,
                                dealii::update_values | dealii::update_gradients | dealii::update_JxW_values);

  // Initialize the quadrature formula and the FEValues object for the faces
  dealii::QGauss<2> face_quadrature_formula(fe.degree + 1);
  dealii::FEFaceValues<3> fe_face_values(fe, face_quadrature_formula,
                                         dealii::update_values | dealii::update_quadrature_points |
                                             dealii::update_normal_vectors | dealii::update_JxW_values);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

  dealii::FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);

  std::vector<dealii::types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators()) {
    // Skip the cell if it is not locally owned.
    if (!cell->is_locally_owned()) {
      continue;
    }

    double sigma = 1.0 / rho[cell->material_id()];

    fe_values.reinit(cell);

    // Assemble the volume integral: \int_{\Omega} \sigma \nabla u \cdot \nabla v dx
    cell_matrix = 0;
    for (const unsigned int q_index : fe_values.quadrature_point_indices()) {
      for (const unsigned int i : fe_values.dof_indices())
        for (const unsigned int j : fe_values.dof_indices())
          cell_matrix(i, j) +=
              sigma * (fe_values.shape_grad(i, q_index) * fe_values.shape_grad(j, q_index) * fe_values.JxW(q_index));
    }
    cell->get_dof_indices(local_dof_indices);

    // Assemble the face integrals
    for (unsigned int f : cell->face_indices()) {
      auto face = cell->face(f);
      if (!face->at_boundary()) {
        continue;
      }

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

    constraints.distribute_local_to_global(cell_matrix, local_dof_indices, system_matrix);
  }

  // Finalize the assembly of the system matrix
  system_matrix.compress(dealii::VectorOperation::add);
}

/**
 * Assembles the right-hand side vector for the primal problem.
 *
 * This function calculates the right-hand side vector caused by the source, i.e., an electrode placed on the surface.
 */
void DCProblem::assemble_rhs() {
  dealii::TimerOutput::Scope timing(timer, "Assemble rhs");

  // Find the active cell containing the source point
  auto cell_point = finder.find_active_cell_around_point(source);

  // Initialize the right-hand side vector
  system_rhs = 0.0;

  if (cell_point.first.state() == dealii::IteratorState::valid && cell_point.first->is_locally_owned()) {
    auto cell_dh = dealii::DoFHandler<3>::active_cell_iterator(*(cell_point.first), &dof_handler);

    // Initialize the quadrature formula using the reference point
    dealii::Quadrature<3> quadrature(dealii::GeometryInfo<3>::project_to_unit_cell(cell_point.second));
    dealii::FEValues<3> fe_values(fe, quadrature, dealii::update_values);

    fe_values.reinit(cell_dh);

    // Compute the integral of the source term
    dealii::Vector<double> local_rhs(fe.dofs_per_cell);
    for (unsigned int i = 0; i < fe.dofs_per_cell; ++i) {
      local_rhs(i) = fe_values.shape_value(i, 0);
    }

    // Distribute the local contributions to the global right-hand side vector
    std::vector<dealii::types::global_dof_index> dof_indices(fe.dofs_per_cell);
    cell_dh->get_dof_indices(dof_indices);

    constraints.distribute_local_to_global(local_rhs, dof_indices, system_rhs);
  }

  // Finalize the assembly of the right-hand side vector
  system_rhs.compress(dealii::VectorOperation::add);
}

/**
 * Assembles the right-hand side vector for the dual problem.
 *
 * The right-hand side vector for the dual problem is assembled by looping over all receivers and
 * treat each receiver as a source. The contribution of each receiver are then assembled into the global vector.
 */
void DCProblem::assemble_rhs_dual() {
  dealii::TimerOutput::Scope timing(timer, "Assemble dual rhs");

  rhs_dual = 0.0;

  // Loop over all receivers
  dealii::Vector<double> local_rhs(fe.dofs_per_cell);
  for (unsigned int i = 0; i < sites.size(); ++i) {
    // Find the active cell containing the receiver
    auto cell_point = finder.find_active_cell_around_point(sites[i]);

    // Skip the receiver if it is not locally owned
    if (cell_point.first.state() != dealii::IteratorState::valid || !cell_point.first->is_locally_owned()) {
      continue;
    }

    auto cell_dh = dealii::DoFHandler<3>::active_cell_iterator(*(cell_point.first), &dof_handler);

    // Initialize the quadrature formula using the reference point
    dealii::Quadrature<3> quadrature(dealii::GeometryInfo<3>::project_to_unit_cell(cell_point.second));
    dealii::FEValues<3> fe_values(fe, quadrature, dealii::update_values);

    fe_values.reinit(cell_dh);

    // Compute the integral of the source term
    local_rhs = 0.0;
    for (unsigned int j = 0; j < fe.dofs_per_cell; ++j) {
      local_rhs(j) = fe_values.shape_value(j, 0);
    }

    // Distribute the local contributions to the global right-hand side vector
    std::vector<dealii::types::global_dof_index> dof_indices(fe.dofs_per_cell);
    cell_dh->get_dof_indices(dof_indices);

    constraints.distribute_local_to_global(local_rhs, dof_indices, rhs_dual);
  }

  // Finalize the assembly of the right-hand side vector
  rhs_dual.compress(dealii::VectorOperation::add);
}

/**
 * Estimates the error for the primal and dual problems.
 */
void DCProblem::estimate_error() {
  dealii::TimerOutput::Scope timing(timer, "Estimate error");

  // Create vectors for the primal and dual solutions with ghost values
  PETScVector solution_ghost(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
  solution_ghost = solution;

  PETScVector solution_dual_ghost(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
  solution_dual_ghost = solution_dual;

  error.reinit(tria.n_active_cells());
  error_dual.reinit(tria.n_active_cells());
  error_go.reinit(tria.n_active_cells());

  std::vector<const dealii::ReadVector<double> *> solutions(2);
  solutions[0] = &solution_ghost;
  solutions[1] = &solution_dual_ghost;

  std::vector<dealii::Vector<float> *> errors_ptr{ &error, &error_dual };
  dealii::ArrayView<dealii::Vector<float> *> errors_view(errors_ptr);

  // Estimate the error
  dealii::KellyErrorEstimator<3>::estimate(dof_handler, dealii::QGauss<2>(fe.degree + 1), {}, dealii::make_array_view(solutions),
                                           errors_view, dealii::ComponentMask(), nullptr, 1, tria.locally_owned_subdomain());

  // Compute the go-oriented error by multiplying the primal and dual errors
  for (unsigned int i = 0; i < tria.n_active_cells(); ++i) {
    error_go[i] = error[i] * error_dual[i];
  }

  // Sum the error contributions from all processors since the errors are only locally available
  MPI_Allreduce(MPI_IN_PLACE, &error_go[0], error_go.size(), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
}

/**
 * Outputs the results for a given refinement step.
 *
 * This function writes the solution, the resistivity, and the error to a VTU file.
 * It also computes the potential the receiver locations and writes the results to a file.
 *
 * @param step The adaptive refinement step.
 */
void DCProblem::output_results(int step) {
  dealii::TimerOutput::Scope timing(timer, "Output results");

  // Create solution vector with ghost values
  PETScVector solution_ghost(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
  solution_ghost = solution;

  // Collect resistivity values for each cell
  dealii::Vector<double> rho(tria.n_active_cells());
  for (auto cell : tria.active_cell_iterators()) {
    rho[cell->active_cell_index()] = this->rho[cell->material_id()];
  }

  // Output the results to a VTU file
  dealii::DataOut<3> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "u");
  data_out.add_data_vector(rho, "rho", dealii::DataOut<3>::type_cell_data);
  data_out.add_data_vector(error, "error", dealii::DataOut<3>::type_cell_data);
  data_out.add_data_vector(error_dual, "error_dual", dealii::DataOut<3>::type_cell_data);
  data_out.add_data_vector(error_go, "error_go", dealii::DataOut<3>::type_cell_data);
  data_out.build_patches();

  data_out.write_vtu_in_parallel((settings.oprefix + "-" + dealii::Utilities::int_to_string(step, 2) + ".vtu").c_str(),
                                 MPI_COMM_WORLD);

  // Compute the potential at the receiver locations
  std::vector<double> u(sites.size(), 0.0);
  for (int i = 0; i < (int)sites.size(); ++i) {
    auto cell_point = finder.find_active_cell_around_point(sites[i]);

    // Only process the receiver if it is located in a locally owned cell
    if (cell_point.first.state() != dealii::IteratorState::valid || !cell_point.first->is_locally_owned()) {
      continue;
    }

    auto cell_dh = dealii::DoFHandler<3>::active_cell_iterator(*(cell_point.first), &dof_handler);

    dealii::Quadrature<3> quadrature(cell_point.second);
    dealii::FEValues<3> fe_values(fe, quadrature, dealii::update_values);

    fe_values.reinit(cell_dh);

    std::vector<dealii::types::global_dof_index> dof_indices(fe.dofs_per_cell);
    cell_dh->get_dof_indices(dof_indices);

    // Compute the potential: u = \sum_i u_i \phi_i
    for (unsigned int j = 0; j < fe.dofs_per_cell; ++j) {
      u[i] += fe_values.shape_value(j, 0) * solution_ghost[dof_indices[j]];
    }
  }

  // Reduce the potential contributions from all processors
  MPI_Allreduce(MPI_IN_PLACE, &u[0], u.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // Write the potential to a file: <output-prefix>-<step>.rsp
  // File format: <receiver-x> <receiver-y> <receiver-z> <potential>
  if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) {
    std::ofstream ofs_rsp(settings.oprefix + "-" + dealii::Utilities::int_to_string(step, 2) + ".rsp");
    for (int i = 0; i < (int)sites.size(); ++i) {
      ofs_rsp << sites[i] << " " << u[i] << std::endl;
    }
  }
}

/**
 * Refine and coarsen the mesh according to the error indicators.
 *
 * Cells with top `refine_fraction` error are refined, and cells with bottom `coarsen_fraction` error are coarsened.
 */
void DCProblem::refine_mesh() {
  dealii::TimerOutput::Scope timing(timer, "Refine mesh");

  dealii::GridRefinement::refine_and_coarsen_fixed_number(tria, error_go, settings.refine_fraction, settings.coarsen_fraction);
  tria.execute_coarsening_and_refinement();
}

/**
 * Runs the DC forward modeling problem.
 */
void DCProblem::run() {
  // Read the model and data files
  read_model_and_data();

  // Loop over the adaptive refinement steps
  for (unsigned int step = 0; step < settings.n_adaptive_refinements + 1; ++step) {
    pcout << "Cycle: " << step << std::endl;
    timer.reset();

    setup_system();       // Set up the system
    assemble_system();    // Assemble the system matrix
    assemble_rhs();       // Assemble the right-hand side vector for the primal problem
    assemble_rhs_dual();  // Assemble the right-hand side vector for the dual problem
    solve();              // Solve the linear system
    estimate_error();     // Estimate the error
    output_results(step); // Output the results
    refine_mesh();        // Refine and coarsen the mesh

    timer.print_summary();
  }
}
