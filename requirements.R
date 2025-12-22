##################### OBIS - species thermal information #######################
# August of 2024
# Authors: Silas C. Principe, Pieter Provoost
# Contact: s.principe@unesco.org
############################ Project Requirements ##############################

# Needed packages on CRAN
req_packages <- c(
  "terra",
  "sf",
  "h3jsr",
  "reticulate",
  "arrow",
  "parallel",
  "dplyr",
  "tidyr",
  "rerddap",
  "ncdf4",
  "fs",
  "httr",
  "storr",
  "progress"
)

# Create a function to check if is installed
is_package_installed <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

# Check which ones are not installed and install if needed:
for (i in seq_along(req_packages)) {
  if (!is_package_installed(req_packages[i])) {
    install.packages(req_packages[i])
  }
}
# Install Python packages
reticulate::py_install(envname = "./.venv", c(
  "xarray",
  "zarr",
  "copernicusmarine",
  "dask[distributed]",
  "bokeh>=3.1.0",
  "geopandas",
  "pyarrow"
), pip = TRUE)
