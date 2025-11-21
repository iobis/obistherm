message("\033[32m======== Preparing virtual environment ======== \033[39m\n")
if (!requireNamespace("reticulate", quietly = TRUE)) {
  install.packages("reticulate")
}
if (!reticulate::virtualenv_exists("./.venv")) {
  reticulate::virtualenv_create("./.venv")
  system('echo "RETICULATE_PYTHON=.venv/bin/python" >> .Renviron')
  system('echo ".Renviron" >> .gitignore')
}

message("\033[32m======== Installing requirements ======== \033[39m\n")
source("requirements.R")