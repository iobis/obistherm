FROM rocker/geospatial:4.4

# Additional system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    awscli \
    python3-venv \
    libnetcdf-dev \
    libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

# Additional R packages not in rocker/geospatial
RUN install2.r --error --skipinstalled \
    storr \
    duckdb \
    arrow \
    rerddap \
    progress \
    glue \
    yaml \
    ncdf4 \
    fs \
    httr \
    h3jsr

WORKDIR /app

# Python virtual environment
RUN python3 -m venv /app/.venv && \
    /app/.venv/bin/pip install --no-cache-dir \
        xarray \
        zarr \
        copernicusmarine \
        "dask[distributed]" \
        "bokeh>=3.1.0" \
        geopandas \
        pyarrow \
        h5py \
        netCDF4

ENV RETICULATE_PYTHON=/app/.venv/bin/python

# Copy project files
COPY functions/ functions/
COPY update_temperatures.R check_files.R ./
COPY settings_docker.yml settings.yml

RUN mkdir -p temp logs

COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
