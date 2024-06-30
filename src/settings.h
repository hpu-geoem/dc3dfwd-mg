#ifndef _SETTINGS_H_
#define _SETTINGS_H_ 1

#include <string>

/**
 * Various configuration parameters for the forward modeling process.
 */
struct Settings {
  /**
   * Tries to parse the settings from a given file.
   *
   * @param prm_filename The filename of the settings file. Note that the file must has extension .prm.
   * @return True if the parsing was successful, false otherwise.
   */
  bool try_parse(const std::string &prm_filename);

  /**
   * Represents the type of solver to be used. The options are geometric multigrid (gmg) and algebraic multigrid (amg).
   */
  enum SolverType { gmg, amg };

  SolverType solver;                   /* The type of solver. */
  unsigned int fe_degree;              /* The degree of the finite element. */
  unsigned int smoother_steps;         /* The number of smoother steps. */
  unsigned int n_adaptive_refinements; /* The number of adaptive refinements. */
  unsigned int n_global_refinements;   /* The number of global refinements. */
  double refine_fraction;              /* The fraction used for refinement. */
  double coarsen_fraction;             /* The fraction used for coarsening. */
  unsigned int max_iter;               /* The maximum number of iterations. */
  double rtol;                         /* The relative tolerance. */
  std::string iprefix, oprefix;        /* The input and output prefixes. */
};

#endif
