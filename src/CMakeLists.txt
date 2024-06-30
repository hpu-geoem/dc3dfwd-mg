SET(dc3dfwd_SOURCES tria_finder.cc settings.cc dc_problem.cc solvers.cc dc3dfwd.cc)

ADD_EXECUTABLE(dc3dfwd ${dc3dfwd_SOURCES})
DEAL_II_SETUP_TARGET(dc3dfwd)
