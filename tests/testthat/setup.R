test_env = rlang::new_environment()

# load saved stanmodel
load("R/sysdata.rda", test_env)

# set up Stan parsing
get_parser()