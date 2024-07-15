# set path to installed libraries
.libPaths("/ix/djishnu/Akanksha/R_libs")

# load necessary libraries
library(yaml)
library(log4r)
library(SLIDE)

# set logger
logger <- logger(threshold = "INFO", appenders = list(file_appender(file = "script.log")))

# script path
yaml_path <- "/ix/djishnu/Akanksha/SLIDE_git/runs/ERCC1_KO_MF/R_scripts/config_ercc1.yaml"
info(logger, paste("Starting script with YAML path:", yaml_path))

# Load YAML file
input_params <- tryCatch(
    {
        yaml::yaml.load_file(yaml_path)
    },
    error = function(e) {
        error(logger, paste("Failed to load YAML file:", e$message))
        stop("Failed to load YAML file")
    }
)

# Log successful loading of YAML file
info(logger, "Successfully loaded YAML file.")

# Run checkDataParams
tryCatch(
    {
        SLIDE::checkDataParams(input_params)
        info(logger, "Successfully ran checkDataParams.")
    },
    error = function(e) {
        error(logger, paste("Error in checkDataParams:", e$message))
        stop("Error in checkDataParams")
    }
)

# Run optimizeSLIDE
tryCatch(
    {
        SLIDE::optimizeSLIDE(input_params, sink_file = FALSE)
        info(logger, "Successfully ran optimizeSLIDE.")
    },
    error = function(e) {
        error(logger, paste("Error in optimizeSLIDE:", e$message))
        stop("Error in optimizeSLIDE")
    }
)

# Run plotCorrelationNetworks
tryCatch(
    {
        SLIDE::plotCorrelationNetworks(input_params)
        info(logger, "Successfully ran plotCorrelationNetworks.")
    },
    error = function(e) {
        error(logger, paste("Error in plotCorrelationNetworks:", e$message))
        stop("Error in plotCorrelationNetworks")
    }
)

# Run Cross Validation
tryCatch(
    {
        SLIDE::SLIDEcv(yaml_path, nrep = 20, k = 5)
        info(logger, "Successfully ran SLIDEcv.")
    },
    error = function(e) {
        error(logger, paste("Error in SLIDEcv:", e$message))
        stop("Error in SLIDEcv")
    }
)

# Log the end of the script
info(logger, "Script execution completed.")
