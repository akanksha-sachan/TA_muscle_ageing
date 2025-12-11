library(devtools)
devtools::install_github("jishnu-lab/SLIDE")


yaml_path = "/ix/djishnu/alw399/SLIDE_py/example2/params.yaml"
input_params <- yaml::yaml.load_file(yaml_path)
SLIDE::checkDataParams(input_params)

SLIDE::optimizeSLIDE(input_params, sink_file = FALSE)

SLIDE::plotCorrelationNetworks(input_params)
