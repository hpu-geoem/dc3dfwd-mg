#include <iostream>

#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>

#include "settings.h"

/**
 * Parse the parameters from a specified file.
 *
 * This function tries to parse the parameters from a specified file.
 * If the parsing is successful, the parameters are stored in the current object.
 * Otherwise, the function returns false.
 *
 * @param prm_filename The filename of the parameter file.
 * @return True if the parsing was successful, false otherwise.
 */
bool Settings::try_parse(const std::string &prm_filename) {
  // Declare the parameters.
  dealii::ParameterHandler prm;
  prm.declare_entry("fe_degree", "1", dealii::Patterns::Integer(1), "Degree of the Finite Element.");
  prm.declare_entry("n_adaptive_refinements", "0", dealii::Patterns::Integer(0, 30), "Number of adaptive refinement steps.");
  prm.declare_entry("n_global_refinements", "0", dealii::Patterns::Integer(0, 10), "Number of global refinement steps.");
  prm.declare_entry("refine_fraction", "0.1", dealii::Patterns::Double(0.0, 1.0), "The fraction of cells to be refined.");
  prm.declare_entry("coarsen_fraction", "0.0", dealii::Patterns::Double(0.0, 1.0), "The fraction of cells to be coarsened.");
  prm.declare_entry("solver", "GMG", dealii::Patterns::Selection("GMG|AMG"), "Switch between GMG and AMG.");
  prm.declare_entry("smoother_steps", "1", dealii::Patterns::Integer(1), "Number of smoother steps.");
  prm.declare_entry("max_iter", "100", dealii::Patterns::Integer(1),
                    "Maximum number of iterations for iterative solvers.");
  prm.declare_entry("rtol", "1E-8", dealii::Patterns::Double(1E-15, 1E-4), "Relative residual tolerance.");

  prm.declare_entry("iprefix", "", dealii::Patterns::FileName(), "Prefix of the input files");
  prm.declare_entry("oprefix", "", dealii::Patterns::FileName(dealii::Patterns::FileName::output), "Prefix of the output files");

  // Check if the input file is provided.
  if (prm_filename.size() == 0) {
    std::cout << "****  Error: No input file provided!\n"
              << "****  Error: Call this program as './dc3dfwd input.prm\n"
              << '\n';
    if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      prm.print_parameters(std::cout, dealii::ParameterHandler::Text);
    return false;
  }

  // Parse the input file.
  try {
    prm.parse_input(prm_filename);
  } catch (std::exception &e) {
    if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      std::cerr << e.what() << std::endl;
    return false;
  }

  // Store the parameters in the current object.
  if (prm.get("solver") == "GMG")
    this->solver = gmg;
  else if (prm.get("solver") == "AMG")
    this->solver = amg;
  else
    AssertThrow(false, dealii::ExcNotImplemented());

  this->fe_degree = prm.get_integer("fe_degree");
  this->smoother_steps = prm.get_integer("smoother_steps");
  this->n_adaptive_refinements = prm.get_integer("n_adaptive_refinements");
  this->n_global_refinements = prm.get_integer("n_global_refinements");
  this->refine_fraction = prm.get_double("refine_fraction");
  this->coarsen_fraction = prm.get_double("coarsen_fraction");
  this->max_iter = prm.get_integer("max_iter");
  this->rtol = prm.get_double("rtol");
  this->iprefix = prm.get("iprefix");
  this->oprefix = prm.get("oprefix");

  return true;
}
