######################### OBIS sea temperature dataset ########################
# January of 2024
# Author: Silas C. Principe
# Contact: s.principe@unesco.org
#
################################# Data download ################################

# Flags
# 0 = no problem
# 1 = date is range
# 2 = GLORYS coordinate is approximated
# 4 = Minimum depth closest value is more than 5 meters different than the true value
# 8 = Maximum depth closest value is more than 5 meters different than the true value
# 16 = CoralTempSST coordinate is approximated
# 32 = MUR SST coordinate is approximated
# 64 = OSTIA SST coordinate is approximated
# Flags can be summed.
# In the new version flags are later decoded and served as strings

# Load packages ----
library(terra)
library(reticulate)
library(arrow)
library(duckdb)
library(parallel)
library(dplyr)
library(tidyr)
library(storr)
# Source functions and settings
source("functions/nearby_from_nc.R")
source("functions/dates.R")
source("functions/download_data.R")
source("functions/utils.R")
source("functions/get_depths.R")
source("functions/data_load.R")
settings <- yaml::read_yaml("settings.yml", readLines.warn = FALSE)
# Load python libraries
check_venv()
cm <- import("copernicusmarine")
xr <- import("xarray")
pd <- import("pandas")

# Settings ----
check_user(settings$copernicus_user, settings$copernicus_password)
options(timeout = 999999999)
is_test <- ifelse(settings$mode == "test", TRUE, FALSE)

# Start Dask for parallel processing
start_dask(browse = ifelse(interactive(), TRUE, FALSE))

st <- storr_rds("control_storr")
outfolder <- settings$outfolder
outfolder_final <- settings$outfolder_final
fs::dir_create(outfolder)
fs::dir_create(outfolder_final)
filename <- "var=thetao"
coordnames <- c("decimalLongitude", "decimalLatitude")
mur_info <- get_mur_ds()

# Define range of dates to get information
range_year <- 1982:lubridate::year(Sys.Date())
range_month <- 1:12

# Define ranges that are available per product
glorys_range <- 1993:lubridate::year(Sys.Date())
coraltemp_range <- 1986:lubridate::year(Sys.Date())
mur_range <- 2002:lubridate::year(Sys.Date())
#https://help.marine.copernicus.eu/en/articles/4872705-what-are-the-main-differences-between-nearrealtime-and-multiyear-products
ostia_rep <- cm$open_dataset(
  dataset_id = "METOFFICE-GLO-SST-L4-REP-OBS-SST",
  username = .user,
  password = .pwd
)
ostia_rep_range <- c(ostia_rep$time$min()$values, ostia_rep$time$max()$values)
ostia_rep_range <- lubridate::year(ostia_rep_range[1]):lubridate::year(ostia_rep_range[2])
ostia_nrt <- cm$open_dataset(
  dataset_id = "METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2",
  username = .user,
  password = .pwd
)
ostia_nrt_range <- c(ostia_nrt$time$min()$values, ostia_nrt$time$max()$values)
ostia_nrt_range <- lubridate::year(ostia_nrt_range[1]):lubridate::year(ostia_nrt_range[2])
ostia_range <- seq(min(ostia_rep_range), max(ostia_nrt_range))

# Define target dataset for Copernicus, time and depths
dataset <- "cmems_mod_glo_phy_my_0.083deg_P1M-m"
product <- "glorys" # This is also used to name the files, as the "core" product
variables <- list("thetao")

ds_sample <- cm$open_dataset(
  dataset_id = dataset,
  # variables = variables,
  username = .user,
  password = .pwd,
  minimum_longitude = -10,
  maximum_longitude = 10,
  minimum_latitude = -10,
  maximum_latitude = 10
)

# Get available depths and maximum date of first dataset
depths <- ds_sample$depth$to_dataframe()
glorys1_max_date <- max(ds_sample$time$to_dataframe()[, 1])

# Open OBIS dataset ------
obis_ds <- get_obis(obis_source = settings$obis_source) |>
  partition_by_year() |>
  open_db()

if (!st$exists("log")) {
  log_df <- data.frame(
    year = rep(range_year, each = 12),
    month = rep(range_month, length(range_year)),
    status_glorys = NA,
    status_coraltemp = NA,
    status_mur = NA,
    status_ostia = NA,
    status_general = NA
  )
  st$set("log", log_df)
}

if (is_test) {
  range_year <- 2010
  range_month <- 1
}

# Get data ------
for (yr in seq_along(range_year)) {
  sel_year <- range_year[yr]
  cat("Processing year", sel_year, "\n")

  # Check if any of the year is pending
  st_yr <- lapply(st$mget(paste0(sel_year, 1:12)), \(x) !is.null(x)) |> unlist()
  if (!all(st_yr)) {
    cat(sum(!st_yr), "months to be processed.\n")
  } else {
    cat("Year already processed, skipping.\n")
    next
  }

  # Load data for a specific year and month
  obis_sel <- load_data_year(sel_year, obis_ds)
  cat(nrow(obis_sel), "total points for this year.\n")

  for (mo in seq_along(range_month)) {
    log_df <- st$get("log")

    sel_month <- range_month[mo]
    st_cod <- paste0(sel_year, sel_month)
    cat("Proccessing month", sel_month, "\n")

    obis_sel_month <- filter_data_month(obis_sel, sel_month)
    cat(nrow(obis_sel_month), "total points for this month.\n")

    if (nrow(obis_sel_month) > 0 && !st$exists(st_cod)) {
      all_vals <- data.frame(
        temp_ID = obis_sel_month$temp_ID,
        surfaceTemperature = NA,
        midTemperature = NA,
        deepTemperature = NA,
        bottomTemperature = NA,
        midDepth = NA,
        deepDepth = NA,
        minimumDepthTemperature = NA,
        maximumDepthTemperature = NA,
        minimumDepthClosestDepth = NA,
        maximumDepthClosestDepth = NA,
        coraltempSST = NA,
        murSST = NA,
        ostiaSST = NA,
        ostiaProduct = NA,
        flag = obis_sel_month$flagDate
      )

      st$set(st_cod, "started")

      obis_dataset <- obis_sel_month |>
          select(decimalLongitude, decimalLatitude, temp_ID,
            depth_min = minimumDepthInMeters, depth_max = maximumDepthInMeters
          )

      # GLORYS PRODUCT ------
      glorys_date <- as.Date(paste0(sel_year, "-", sprintf("%02d", sel_month), "-01"))
      if (sel_year %in% glorys_range && glorys_date <= glorys1_max_date) {
        cat("Downloading GLORYS\n")
        # Interim no longer available, modified conditional to proceed only if year and month are available
        # if (as.Date(paste0(sel_year, "-", sprintf("%02d", sel_month), "-01")) <= glorys1_max_date) {
        #   dataset <- "cmems_mod_glo_phy_my_0.083deg_P1M-m"
        # } else {
        #   dataset <- "cmems_mod_glo_phy_myint_0.083deg_P1M-m"
        # }

        outf_temp_glorys <- cm$get(
          dataset_id = dataset,
          username = .user,
          password = .pwd,
          filter = paste0("*mean_", sel_year, sprintf("%02d", sel_month), "*"),
          output_directory = "temp/",
          no_directories = T
        )
        outf_temp_glorys <- lapply(outf_temp_glorys$files, \(x) x$file_path) |>
          check_ifdate(sel_year, sel_month, "_", target = 1) |>
          (\(x) unlist(lapply(x, as.character), use.names = F))()

        obis_dataset <- obis_dataset |>
          mutate(
            to_remove_dmin = ifelse(is.na(depth_min), TRUE, FALSE),
            to_remove_dmax = ifelse(is.na(depth_max), TRUE, FALSE),
            depth_min = ifelse(is.na(depth_min), 0, depth_min),
            depth_max = ifelse(is.na(depth_max), 0, depth_max),
            depth_surface = 0,
            to_remove_nacord = FALSE,
            to_remove_naaprox = FALSE
          )

        # See which are NA
        cat("Processing GLORYS\n")
        coords_na <- extract_from_nc(outf_temp_glorys, "thetao", obis_dataset[, 1:2], depth = 0L)
        coords_na <- which(is.na(coords_na[, "value"]))

        if (length(coords_na) > 1) {
          new_coords <- get_nearby(outf_temp_glorys, "thetao", obis_dataset[coords_na, 1:2], mode = "25", depth = 0L)
          na_to_rm <- which(is.na(new_coords[, "value"]))
          na_approx <- which(!is.na(new_coords[, "value"]))

          obis_dataset$decimalLongitude[coords_na[na_approx]] <- new_coords$new_lon[na_approx]
          obis_dataset$decimalLatitude[coords_na[na_approx]] <- new_coords$new_lat[na_approx]
          obis_dataset$to_remove_nacord[coords_na[na_to_rm]] <- TRUE
          obis_dataset$to_remove_naaprox[coords_na[na_approx]] <- TRUE
        }

        ds_temp_glorys <- xr$open_dataset(outf_temp_glorys, chunks = "auto")
        valid_depths <- get_depths(ds_temp_glorys, obis_dataset,
         paste(sel_year, sprintf("%02d", sel_month), "01", sep = "-"))

        valid_depths <- valid_depths |>
          mutate(to_remove_ddeep = ifelse(is.na(depth_deep), TRUE, FALSE),
                 to_remove_dmid = ifelse(is.na(depth_mid), TRUE, FALSE)) |>
          mutate(depth_deep = ifelse(is.na(depth_deep), 0, depth_deep),
                 depth_mid = ifelse(is.na(depth_mid), 0, depth_mid))

        obis_dataset <- left_join(obis_dataset, valid_depths, by = "temp_ID")

        to_remove <- obis_dataset |>
          select(starts_with("to_remove"))

        obis_dataset <- obis_dataset |>
          select(!starts_with("to_remove"))

        success <- try(download_temp("glorys", ds_temp_glorys, obis_dataset, sel_month, sel_year))

        if (!inherits(success, "try-error")) {
          glorys_data <- success |>
            select(temp_ID, depth_type, value) |>
            pivot_wider(names_from = depth_type, values_from = value)

          glorys_data_depths <- success |>
            select(temp_ID, depth_type, depth) |>
            mutate(depth_type = paste0(depth_type, "_depth")) |>
            pivot_wider(names_from = depth_type, values_from = depth)

          rm(success)

          glorys_data <- left_join(glorys_data, glorys_data_depths, by = "temp_ID")
          glorys_data <- tibble::as_tibble(apply(glorys_data, 2, function(x) {
            x[is.nan(x)] <- NA
            x
          }))

          # Input data
          all_vals$surfaceTemperature[all_vals$temp_ID == glorys_data$temp_ID] <- glorys_data$depth_surface
          all_vals$midTemperature[all_vals$temp_ID == glorys_data$temp_ID] <- glorys_data$depth_mid
          all_vals$deepTemperature[all_vals$temp_ID == glorys_data$temp_ID] <- glorys_data$depth_deep
          all_vals$bottomTemperature[all_vals$temp_ID == glorys_data$temp_ID] <- glorys_data$depth_bottom

          all_vals$minimumDepthTemperature[all_vals$temp_ID == glorys_data$temp_ID] <- glorys_data$depth_min
          all_vals$maximumDepthTemperature[all_vals$temp_ID == glorys_data$temp_ID] <- glorys_data$depth_max

          all_vals$midDepth[all_vals$temp_ID == glorys_data$temp_ID] <- glorys_data$depth_mid_depth
          all_vals$deepDepth[all_vals$temp_ID == glorys_data$temp_ID] <- glorys_data$depth_deep_depth

          all_vals$minimumDepthClosestDepth[all_vals$temp_ID == glorys_data$temp_ID] <- glorys_data$depth_min_depth
          all_vals$maximumDepthClosestDepth[all_vals$temp_ID == glorys_data$temp_ID] <- glorys_data$depth_max_depth

          all_vals$minimumDepthTemperature[to_remove$to_remove_dmin] <- NA
          all_vals$minimumDepthClosestDepth[to_remove$to_remove_dmin] <- NA

          all_vals$maximumDepthTemperature[to_remove$to_remove_dmax] <- NA
          all_vals$maximumDepthClosestDepth[to_remove$to_remove_dmax] <- NA

          all_vals$midTemperature[to_remove$to_remove_dmid] <- NA
          all_vals$midDepth[to_remove$to_remove_dmid] <- NA

          all_vals$deepTemperature[to_remove$to_remove_ddeep] <- NA
          all_vals$deepDepth[to_remove$to_remove_ddeep] <- NA

          all_vals[to_remove$to_remove_nacord,c("surfaceTemperature", "midTemperature",
           "deepTemperature", "midDepth", "deepDepth", "minimumDepthTemperature", "maximumDepthTemperature",
          "minimumDepthClosestDepth", "maximumDepthClosestDepth")] <- NA

          all_vals$flag[to_remove$to_remove_naaprox] <- all_vals$flag[to_remove$to_remove_naaprox] + 2

          depth_diff_min <- check_depth_diff(obis_dataset$depth_min, all_vals$minimumDepthClosestDepth)
          depth_diff_max <- check_depth_diff(obis_dataset$depth_max, all_vals$maximumDepthClosestDepth)

          all_vals$flag[depth_diff_min] <- all_vals$flag[depth_diff_min] + 4
          all_vals$flag[depth_diff_max] <- all_vals$flag[depth_diff_max] + 8

          fs::file_delete(outf_temp_glorys)

          st$set(st_cod, c(st$get(st_cod), "glorys"))
          log_df[log_df$year == sel_year & log_df$month == sel_month, "status_glorys"] <- "concluded"
        } else {
          st$set(st_cod, c(st$get(st_cod), "glorys"))
          log_df[log_df$year == sel_year & log_df$month == sel_month, "status_glorys"] <- "failed_extract"
        }
      } else {
        log_df[log_df$year == sel_year & log_df$month == sel_month, "status_glorys"] <- "unavailable"
      }

      # CORALTEMP PRODUCT ----
      if (sel_year %in% coraltemp_range) {
        cat("Downloading CoralTemp\n")
        cttemp <- file.path(outfolder, "coraltemp.nc")

        if (nrow(obis_dataset) <= 100) {
          coraltemp_ds <- "https://coastwatch.pfeg.noaa.gov/erddap/griddap/NOAA_DHW_monthly"
          proceed <- TRUE
        } else {
          df <- try(download.file(url = paste0(
            "https://coastwatch.pfeg.noaa.gov/erddap/files/NOAA_DHW_monthly/ct5km_sst_ssta_monthly_v31_",
            sel_year, sprintf("%02d", sel_month), ".nc"),
            destfile = cttemp, method = "libcurl", mode = "wb"), silent = T)
          if (!inherits(df, "try-error")) {
            proceed <- TRUE
            coraltemp_ds <- cttemp
          } else {
            proceed <- FALSE
          }
        }

        if (proceed) {
          cat("Extracting CoralTemp\n")
          success <- try(download_temp("coraltemp", coraltemp_ds, obis_dataset, sel_month, sel_year))
        } else {
          success <- try(stop("error"), silent = T)
        }

        if (!inherits(success, "try-error")) {
          cat("Processing CoralTemp\n")

          success$value[is.nan(success$value)] <- NA

          na_to_solve <- which(is.na(success[, "value"]))

          nearby_ct <- get_nearby(coraltemp_ds, "sea_surface_temperature",
                                  obis_dataset[na_to_solve, c("decimalLongitude", "decimalLatitude")],
                                  mode = "25", depth = NULL,
                                  date = paste(sel_year, sprintf("%02d", sel_month), "15", sep = "-"))

          all_vals$coraltempSST <- success$value
          all_vals$coraltempSST[na_to_solve] <- nearby_ct$value
          all_vals$flag[na_to_solve[!is.na(nearby_ct$value)]] <-
            all_vals$flag[na_to_solve[!is.na(nearby_ct$value)]] + 16

          st$set(st_cod, c(st$get(st_cod), "coraltemp"))
          log_df[log_df$year == sel_year & log_df$month == sel_month, "status_coraltemp"] <- "concluded"
        } else {
          st$set(st_cod, c(st$get(st_cod), "coraltemp"))
          log_df[log_df$year == sel_year & log_df$month == sel_month, "status_coraltemp"] <- "failed_extract"
        }
        fs::file_delete(list.files(outfolder, full.names = T))
      } else {
        log_df[log_df$year == sel_year & log_df$month == sel_month, "status_coraltemp"] <- "unavailable"
      }


      # MUR PRODUCT ----
      if (sel_year %in% mur_range) {
        cat("Downloading MUR\n")
        murtemp <- file.path(outfolder, "murtemp.nc")

        if (nrow(obis_dataset) <= 100) {
          mur_ds <- "https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41mday"
          proceed <- TRUE
        } else {
          url_try <- c(
            paste0(
              "https://coastwatch.pfeg.noaa.gov/erddap/files/jplMURSST41mday/",
              paste0(
                sel_year, sprintf("%02d", sel_month), "01", sel_year, sprintf("%02d", sel_month),
                lubridate::days_in_month(
                  as.Date(paste0(sel_year, sprintf("%02d", sel_month), "01"),
                    format = "%Y%m%d"
                  )
                )
              ),
              "-GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc"
            )
          )
          df <- safe_download_mur(url = url_try,
            destfile = murtemp, mur_info = mur_info, method = "libcurl", mode = "wb")
          if (!inherits(df, "try-error")) {
            proceed <- TRUE
            mur_ds <- murtemp
          } else {
            proceed <- FALSE
          }
        }

        if (proceed) {
          cat("Extracting MUR\n")
          success <- try(download_temp("mur", mur_ds, obis_dataset, sel_month, sel_year))
        } else {
          success <- try(stop("error"), silent = T)
        }

        if (!inherits(success, "try-error")) {
          cat("Processing MUR\n")

          success$value[is.nan(success$value)] <- NA

          na_to_solve <- which(is.na(success[, "value"]))

          nearby_mur <- get_nearby(mur_ds, "sst",
                                  obis_dataset[na_to_solve, c("decimalLongitude", "decimalLatitude")],
                                  mode = "25", depth = NULL,
                                  date = paste(sel_year, sprintf("%02d", sel_month), "15", sep = "-"))

          all_vals$murSST <- success$value
          all_vals$murSST[na_to_solve] <- nearby_mur$value
          all_vals$flag[na_to_solve[!is.na(nearby_mur$value)]] <-
            all_vals$flag[na_to_solve[!is.na(nearby_mur$value)]] + 32

          st$set(st_cod, c(st$get(st_cod), "mur"))
          log_df[log_df$year == sel_year & log_df$month == sel_month, "status_mur"] <- "concluded"
        } else {
          st$set(st_cod, c(st$get(st_cod), "mur"))
          log_df[log_df$year == sel_year & log_df$month == sel_month, "status_mur"] <- "failed_extract"
        }
        fs::file_delete(list.files(outfolder, full.names = T))
      } else {
        log_df[log_df$year == sel_year & log_df$month == sel_month, "status_mur"] <- "unavailable"
      }


      # OSTIA PRODUCT ----
      if (sel_year %in% ostia_range) {
        cat("Downloading OSTIA\n")

        if (sel_year %in% ostia_rep_range) {
          ostia_id <- "METOFFICE-GLO-SST-L4-REP-OBS-SST"
          ostia_prod <- "REP"
        } else {
          ostia_id <- "METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2"
          ostia_prod <- "NRT"
        }
        ostia_ds <- file.path(outfolder, "ostiatemp.nc")

        df <- try(
            cm$get(
              dataset_id = ostia_id,
              username = .user,
              password = .pwd,
              filter = paste0("*", sel_year, sprintf("%02d", sel_month), "*"),
              output_directory = "temp/",
              no_directories = T
            )
          )
        if (!inherits(df, "try-error")) {
          df <- lapply(df$files, \(x) x$file_path)
          df <- try(check_ifdate(df, sel_year, sel_month, target = c(ifelse(sel_month == 2, 28, 30), 31)), silent = TRUE)
        }

        if (!inherits(df, "try-error")) {
          cat("Extracting OSTIA\n")
          df_ds <- xr$open_mfdataset(unlist(lapply(df, as.character), recursive = T))
          df_ds <- df_ds$analysed_sst
          df_ds <- df_ds$mean(dim = "time", skipna = T)
          df_ds <- df_ds$to_dataset()
          df_ds <- df_ds$expand_dims(time = pd$to_datetime(paste0(sel_year, "-", sprintf("%02d", sel_month), "-01")))
          df_ds$to_netcdf(ostia_ds)

          ostia_ds_open <- xr$open_dataset(ostia_ds)
          success <- try(download_temp("ostia", ostia_ds_open, obis_dataset, sel_month, sel_year))
        } else {
          success <- try(stop("error"), silent = T)
        }

        if (!inherits(success, "try-error")) {
          cat("Processing OSTIA\n")

          success$value[is.nan(success$value)] <- NA

          na_to_solve <- which(is.na(success[, "value"]))

          nearby_ostia <- get_nearby(ostia_ds, "analysed_sst",
                                  obis_dataset[na_to_solve, c("decimalLongitude", "decimalLatitude")],
                                  mode = "25", depth = NULL,
                                  date = paste(sel_year, sprintf("%02d", sel_month), "01", sep = "-"))
          nearby_ostia$value <- nearby_ostia$value - 273.15

          all_vals$ostiaSST <- success$value
          all_vals$ostiaSST[na_to_solve] <- nearby_ostia$value
          all_vals$flag[na_to_solve[!is.na(nearby_ostia$value)]] <-
            all_vals$flag[na_to_solve[!is.na(nearby_ostia$value)]] + 64
          all_vals$ostiaProduct[!is.na(all_vals$ostiaSST)] <- ostia_prod

          st$set(st_cod, c(st$get(st_cod), "ostia"))
          log_df[log_df$year == sel_year & log_df$month == sel_month, "status_ostia"] <- "concluded"
        } else {
          st$set(st_cod, c(st$get(st_cod), "mur"))
          log_df[log_df$year == sel_year & log_df$month == sel_month, "status_ostia"] <- "failed_extract"
        }
        fs::file_delete(list.files(outfolder, full.names = T))
      } else {
        log_df[log_df$year == sel_year & log_df$month == sel_month, "status_ostia"] <- "unavailable"
      }

      all_vals <- all_vals |> rename(obistherm_flags = flag)
      all_vals <- left_join(obis_sel_month, all_vals, by = "temp_ID") |>
        select(-temp_ID, -flagDate)

      have_data <- apply(
        all_vals |>
          select(
            surfaceTemperature, midTemperature, deepTemperature,
            bottomTemperature,
            minimumDepthTemperature, maximumDepthTemperature,
            coraltempSST, murSST, ostiaSST
          ),
        1, function(x) any(!is.na(x))
      )

      all_vals <- all_vals[have_data, ]

      if (nrow(all_vals) > 0) {
        write_parquet(all_vals, paste0(outfolder_final, "/", filename, "_year=", sel_year, "_month=", sel_month, ".parquet"))
      }

      d_files <- list.files(outfolder, full.names = T)
      d_files <- d_files[grepl(paste0(sel_year, "_", sel_month), d_files)]
      fs::file_delete(d_files)

      st$set(st_cod, c(st$get(st_cod), "done"))
      log_df[log_df$year == sel_year & log_df$month == sel_month, "status_general"] <- "concluded"
    } else {
      if (st$exists(st_cod)) {
        cat("Already done\n")
      } else {
        log_df[log_df$year == sel_year & log_df$month == sel_month, "status_general"] <- "skipped"
      }
    }
    cat("Month log:\n")
    print(log_df[log_df$year == sel_year & log_df$month == sel_month, ])
    cat("========================\n\n")
    st$set("log", log_df)
  }
}

fs::dir_create("logs")
log_df <- st$get("log")
write.csv(log_df, paste0("logs/log_", format(Sys.Date(), "%Y%m%d"), ".csv"), row.names = F)

# END