#include <deal.II/base/mpi.h>

#include "dc_problem.h"

int main(int argc, char **argv) {
  // Initialize MPI and finalize it at the end of the program.
  dealii::Utilities::MPI::MPI_InitFinalize init(argc, argv, 1);

  // Parse the Parameters from the input file.
  Settings settings;
  if (!settings.try_parse((argc > 1) ? (argv[1]) : ""))
    return 0;

  // Create a DC problem and run it.
  DCProblem dc(settings);
  dc.run();

  return 0;
}
