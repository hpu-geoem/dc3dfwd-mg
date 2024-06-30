#include "tria_finder.h"
#ifndef _DC_PROBLEM_H_
#define _DC_PROBLEM_H_ 1

#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/timer.h>

#include <deal.II/grid/tria.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>

#include "settings.h"
#include "tria_finder.h"

/**
 * Class representing a DC problem.
 *
 * This class encapsulates the functionality for solving a DC problem.
 * It provides methods for reading the model and data, setting up the system,
 * assembling the system matrix and right-hand side, solving the problem,
 * estimating the error, and outputting the results.
 *
 * The documentation for the methods is provided in the implementation file.
 */
class DCProblem {
  using PETScSparseMatrix = dealii::PETScWrappers::MPI::SparseMatrix;
  using PETScVector = dealii::PETScWrappers::MPI::Vector;

public:
  DCProblem(Settings s);

  void run();

private:
  void read_model_and_data();
  void setup_system();
  void setup_multigrid();
  void assemble_system();
  void assemble_multigrid();
  void assemble_rhs();
  void assemble_rhs_dual();
  void solve();
  void solve_amg();
  void solve_gmg();
  void estimate_error();
  void output_results(int step);
  void refine_mesh();

  Settings settings; // Parameters for the forward modeling process.

  dealii::Point<3> source;             // The source location.
  std::vector<dealii::Point<3>> sites; // The receiver locations.

  std::vector<double> rho; // The resistivity values.

  dealii::FE_Q<3> fe;                              // The finite element.
  dealii::parallel::shared::Triangulation<3> tria; // The triangulation.
  dealii::DoFHandler<3> dof_handler;               // The DoF handler.

  TriaFinder finder; // The triangulation finder.

  dealii::AffineConstraints<double> constraints; // The constraints object for the hanging nodes.
  PETScSparseMatrix system_matrix;               // The system matrix.

  PETScVector solution, system_rhs; // The solution and right-hand side vectors for the primal problem.

  PETScVector solution_dual; // The solution vector for the dual problem.
  PETScVector rhs_dual;      // The right-hand side vector for the dual problem.

  dealii::Vector<float> error, error_dual, error_go; // The error vectors.

  dealii::IndexSet locally_owned_dofs, locally_relevant_dofs; // Index sets for the locally owned and relevant DoFs.

  dealii::MGLevelObject<PETScSparseMatrix> mg_matrix;       // The system matrix for the multigrid hierarchy.
  dealii::MGLevelObject<PETScSparseMatrix> mg_interface_in; // The interface matrix for the multigrid hierarchy.
  dealii::MGConstrainedDoFs mg_constrained_dofs;            // The constrained DoFs for the multigrid hierarchy.

  dealii::ConditionalOStream pcout; // Conditional output stream for parallel programs.
  dealii::TimerOutput timer;        // Timer for the various parts of the program.
};

#endif
