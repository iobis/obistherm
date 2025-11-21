check_user <- function(user = NULL, pwd = NULL) {
    if (Sys.getenv("COPERNICUS_USER") != "") {
        .user <- Sys.getenv("COPERNICUS_USER")
        .pwd <- Sys.getenv("COPERNICUS_PWD")
    } else {
        if (interactive()) {
            .user <- readline("Enter your user\n")
            .pwd <- readline("Enter your password\n")
        } else {
            if (any(is.null(c(user, pwd)))) {
                stop("You need to supply the copernicus user and password in settings.yml")
            } else {
                .user <- user
                .pwd <- pwd
            }
        }
    }
    .user <<- .user
    .pwd <<- .pwd
    return(invisible())
}

get_obis <- function(obis_source = NULL) {
    message("Downloading/updating OBIS dataset...")

    if (is.null(obis_source)) {
        f <- "data/obis_data"
    } else {
        f <- obis_source
    }
    fs::dir_create(f)
    # Test if aws is available
    if (!nzchar(Sys.which("aws"))) {
        stop("aws is not available on the system; install it to download OBIS data")
    }
    comm <- paste0("aws s3 sync --no-sign-request s3://obis-open-data/occurrence/ ", f)
    system(comm)

    writeLines(paste("obis_data dataset downloaded on", format(Sys.time()), "using the command", comm,
                    "\n\nCite as: Ocean Biodiversity Information System (OBIS) (25 March 2025) OBIS Occurrence Data. Intergovernmental Oceanographic Commission of UNESCO. https://doi.org/10.25607/obis.occurrence.b89117cd."),
            "data/obis-download.txt")
    return(f)
}

partition_by_year <- function(obis_source, min_year = 1982, skip = FALSE) {

    message("Partitioning dataset by year")
    require(duckdb)

    f <- dirname(obis_source) |>
        file.path("obis_partitioned")
    
    if (skip) {
        message("Skipping data preparing and using old version at ", f)
        return(f)
    }

    if (dir.exists(f)) {
        message("Cleaning folder...")
        fs::dir_delete(f)
    }

    con <- DBI::dbConnect(duckdb())

    DBI::dbSendQuery(con, paste0("
    COPY(
        SELECT _id, dataset_id, source.occurrenceID as occurrenceID, source.datasetID as datasetID,
            interpreted.aphiaid as AphiaID, interpreted.scientificName as scientificName,
            interpreted.species as species, interpreted.genus as genus, interpreted.family as family,
            interpreted.order as order, interpreted.class as class, interpreted.phylum as phylum,
            interpreted.kingdom as kingdom,
            interpreted.eventDate as eventDate, interpreted.date_start as date_start,
            interpreted.date_mid as date_mid, interpreted.date_end as date_end, interpreted.date_year as date_year,
            interpreted.decimalLongitude as decimalLongitude, interpreted.decimalLatitude as decimalLatitude,
            interpreted.coordinatePrecision as coordinatePrecision, interpreted.coordinateUncertaintyInMeters as coordinateUncertaintyInMeters,
            interpreted.minimumDepthInMeters as minimumDepthInMeters, interpreted.maximumDepthInMeters as maximumDepthInMeters,
            dropped, absence, flags
        FROM read_parquet('/Volumes/OBIS2/data/obis_data/*')
        WHERE dropped IS NOT TRUE AND date_year >= ", min_year, "
    ) TO '/Volumes/OBIS2/data/obis_partitioned'
    (FORMAT parquet, PARTITION_BY (date_year));
    "))

    DBI::dbDisconnect(con)

    return(f)
}

open_db <- function(source_folder) {
    require(duckdb)

    con <- DBI::dbConnect(duckdb())

    DBI::dbSendQuery(con, paste0(
    "
    CREATE VIEW obis AS SELECT * FROM read_parquet('", file.path(source_folder, "*/*"), "');
    "
    ))

    return(con)
}

check_depth_diff <- function(depth_original, depth_new, limit = 5) {
    diff_depth <- abs(depth_original - depth_new)
    diff_depth <- ifelse(diff_depth > limit, TRUE, FALSE)
    ifelse(is.na(diff_depth), FALSE, diff_depth)
}

start_dask <- function(browse = TRUE) {
    da <- import("dask")
    dd <- import("dask.distributed")
    client <- dd$Client()
    if (browse) {
        browseURL("http://localhost:8787/status")
    } else {
        cat("For browsing, access http://localhost:8787/status\n")
    }
    return(client)
}

decode_flag <- function(flag, collapse = TRUE) {
  
  if (!is.na(flag)) {
    flag_text <- c(
        "1"  = "date is range",
        "2"  = "GLORYS coordinate is approximated",
        "4"  = "Minimum depth differs >5m from true value",
        "8"  = "Maximum depth differs >5m from true value",
        "16" = "CoralTempSST coordinate is approximated",
        "32" = "MUR SST coordinate is approximated",
        "64" = "OSTIA SST coordinate is approximated"
    )
    
    flag_id <- as.numeric(names(flag_text))
    
    active <- flag_id[bitwAnd(flag, flag_id) > 0]
    
    if (length(active) == 0) {
        return(NA)
    }
    
    result <- flag_text[as.character(active)]
    
    if (collapse) {
        paste(result, collapse = "; ")
    } else {
        result
    }
  } else {
    return(NA)
  }

}
#decode_flag(64+32+16+8+4+2+1)

check_venv <- function() {
    use_python(".venv/bin/python", required = TRUE)
    if (!grepl("\\.venv/bin/python", py_config()$python)) stop("Check virtual environment")
    return(invisible(NULL))
}

check_ifdate <- function(files, sel_year, sel_month, start_chr = "/", target = NULL) {
    files_sb <- unlist(lapply(files, as.character), recursive = T)
    files <- files[grepl(paste0(start_chr, sel_year, sprintf("%02d", sel_month)), files_sb)]
    if (length(files) > 31) {
        stop("More days than maximum in a month")
    } else if (length(files) < 1) {
        stop("No valid files")
    } else {
        if (!is.null(target)) {
            if (length(target) == 1) {
                if (length(files) != target) stop("Number of files different than target")
            } else {
                if (length(files) < target[1]) {
                    stop("Less files than target[1]")
                }
                if (length(files) > target[2]) {
                    stop("More files than target[2]")
                }
            }
        }
        return(files)
    }
}
