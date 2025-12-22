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
source("functions/process_functions.R")
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
safe_limit <- settings$safe_limit_points |> as.integer()
download_lib <- settings$download_library
stop_if_fail <- settings$stop_if_fail
backup_noaa <- settings$backup_noaa_files
backup_folder <- settings$backup_noaa_folder
bypass_download <- settings$bypass_if_exists

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
ctemp_info <- get_ctemp_ds()

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
  partition_by_year(skip = T) |>
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
  catg("Processing year", sel_year)

  # Check if any month of the year is pending
  st_yr <- lapply(st$mget(paste0(sel_year, 1:12)), \(x) !is.null(x)) |> unlist()
  if (!all(st_yr)) {
    catn(sum(!st_yr), "months to be processed.")
  } else {
    catn("Year already processed, skipping.")
    next
  }

  # Load data for a specific year and month
  obis_sel <- load_data_year(sel_year, obis_ds, outfolder)
  catn(nrow(obis_sel), "total points for this year.")

  for (mo in seq_along(range_month)) {
    # Get log file
    log_df <- st$get("log")

    sel_month <- range_month[mo]
    st_cod <- paste0(sel_year, sel_month)
    catn("Proccessing month", sel_month)

    if (lubridate::year(Sys.Date()) == sel_year && lubridate::month(Sys.Date()) == sel_month) {
      catn("Skipping current month")
      next
    }

    # Load data for month
    obis_sel_month_total <- filter_data_month(obis_sel, sel_month)
    catn(nrow(obis_sel_month_total), "total points for this month.")

    if (nrow(obis_sel_month_total) > 0 && !st$exists(st_cod)) {

      if (nrow(obis_sel_month_total) > safe_limit) {
        ds_parts <- split(
          seq_len(nrow(obis_sel_month_total)),
          ceiling(seq_len(nrow(obis_sel_month_total)) / safe_limit)
        )
        catn("Dataset with more rows than safety limit. Dividing in parts.")
      } else {
        ds_parts <- list(seq_len(nrow(obis_sel_month_total)))
      }
      total_parts <- length(ds_parts)
      st$set(st_cod, "started")

      for (part in seq_len(total_parts)) {

        if (total_parts > 1) catn("Processing part", part, "out of", total_parts)

        obis_sel_month <- obis_sel_month_total[ds_parts[[part]],]

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

        obis_dataset <- obis_sel_month |>
            select(decimalLongitude, decimalLatitude, temp_ID,
              depth_min = minimumDepthInMeters, depth_max = maximumDepthInMeters
            )


        # PROCESS PRODUCTS --------
        # GLORYS PRODUCT ------
        glorys_date <- as.Date(paste0(sel_year, "-", sprintf("%02d", sel_month), "-01"))
        if (sel_year %in% glorys_range && glorys_date <= glorys1_max_date) {

          if (part == 1) {
            catn("Downloading GLORYS")
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
          }

          glorys_result <- process_glorys(all_vals, obis_dataset, outf_temp_glorys, sel_month, sel_year)

          if (glorys_result$status == "succeeded") {
            all_vals <- glorys_result$result
            st$set(st_cod, c(st$get(st_cod), paste0("glorys_", part)))
          }
          rm(glorys_result)

          if (part == total_parts) {
            if (sum(grepl("glorys", st$get(st_cod))) == total_parts) {
              log_df <- log_df |> set_success(sel_year, sel_month, "status_glorys")
            } else {
              log_df <- log_df |> set_failed(sel_year, sel_month, "status_glorys")
            }
            fs::file_delete(outf_temp_glorys)
          }
        } else {
          log_df <- log_df |> set_unavailable(sel_year, sel_month, "status_glorys")
        }


        # CORALTEMP PRODUCT ----
        if (sel_year %in% coraltemp_range) {
          cttemp <- file.path(outfolder, "coraltemp.nc")

          if (part == 1) {
            catn("Downloading CoralTemp")
            coraltemp_ds <- bypass_try("coraltemp", backup_folder, sel_year, sel_month, bypass_download)
            if (is.null(coraltemp_ds)) {
              if (nrow(obis_dataset) <= 100) {
                coraltemp_ds <- "https://coastwatch.pfeg.noaa.gov/erddap/griddap/NOAA_DHW_monthly"
                proceed <- TRUE
              } else {
                df_ctemp <- safe_download_ctemp(sel_year = sel_year, sel_month = sel_month,
                  destfile = cttemp, ctemp_info = ctemp_info, method = download_lib, mode = "wb")
                if (!inherits(df_ctemp, "try-error")) {
                  proceed <- TRUE
                  coraltemp_ds <- cttemp
                } else {
                  proceed <- FALSE
                  if (stop_if_fail) stop("CoralTemp failed.")
                }
              }
            } else {
              proceed <- TRUE
            }
          }

          if (proceed) {
            coraltemp_result <- process_coraltemp(all_vals, obis_dataset, coraltemp_ds, sel_month, sel_year)
          } else {
            coraltemp_result <- list(status = "failed")
          }

          if (coraltemp_result$status == "succeeded") {
            all_vals <- coraltemp_result$result
            st$set(st_cod, c(st$get(st_cod), paste0("coraltemp_", part)))
          }
          rm(coraltemp_result)

          if (part == total_parts) {
            if (sum(grepl("coraltemp", st$get(st_cod))) == total_parts) {
              log_df <- log_df |> set_success(sel_year, sel_month, "status_coraltemp")
            } else {
              log_df <- log_df |> set_failed(sel_year, sel_month, "status_coraltemp")
            }
            if (file.exists(cttemp)) {
              backup_it(cttemp, "coraltemp", sel_year, sel_month, backup_folder, backup_noaa) |> 
                fs::file_delete()
            }
          }
        } else {
          log_df <- log_df |> set_unavailable(sel_year, sel_month, "status_coraltemp")
        }


        # MUR PRODUCT ----
        if (sel_year %in% mur_range) {
          murtemp <- file.path(outfolder, "murtemp.nc")

          if (part == 1) {
            catn("Downloading MUR")
            mur_ds <- bypass_try("mur", backup_folder, sel_year, sel_month, bypass_download)
            if (is.null(mur_ds)) {
              if (nrow(obis_dataset) <= 100) {
                mur_ds <- "https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplMURSST41mday"
                proceed <- TRUE
              } else {
                df_mur <- safe_download_mur(sel_year = sel_year, sel_month = sel_month,
                  destfile = murtemp, mur_info = mur_info, method = download_lib, mode = "wb")
                if (!inherits(df_mur, "try-error")) {
                  proceed <- TRUE
                  mur_ds <- murtemp
                } else {
                  proceed <- FALSE
                  if (stop_if_fail) stop("MUR failed.")
                }
              }
            } else {
              proceed <- TRUE
            }
          }

          if (proceed) {
            mur_result <- process_mur(all_vals, obis_dataset, mur_ds, sel_month, sel_year)
          } else {
            mur_result <- list(status = "failed")
          }

          if (mur_result$status == "succeeded") {
            all_vals <- mur_result$result
            st$set(st_cod, c(st$get(st_cod), paste0("mur_", part)))
          }
          rm(mur_result)

          if (part == total_parts) {
            if (sum(grepl("mur", st$get(st_cod))) == total_parts) {
              log_df <- log_df |> set_success(sel_year, sel_month, "status_mur")
            } else {
              log_df <- log_df |> set_failed(sel_year, sel_month, "status_mur")
            }
            if (file.exists(murtemp)) {
              backup_it(murtemp, "mur", sel_year, sel_month, backup_folder, backup_noaa) |>
                fs::file_delete()
            }
          }
        } else {
          log_df <- log_df |> set_unavailable(sel_year, sel_month, "status_mur")
        }


        # OSTIA PRODUCT ----
        if (sel_year %in% ostia_range) {

          if (sel_year %in% ostia_rep_range) {
            ostia_id <- "METOFFICE-GLO-SST-L4-REP-OBS-SST"
            ostia_prod <- "REP"
          } else {
            ostia_id <- "METOFFICE-GLO-SST-L4-NRT-OBS-SST-V2"
            ostia_prod <- "NRT"
          }
          ostia_ds <- file.path(outfolder, "ostiatemp.nc")

          if (part == 1) {
            catn("Downloading OSTIA")
            df_ostia <- try(
                cm$get(
                  dataset_id = ostia_id,
                  username = .user,
                  password = .pwd,
                  filter = paste0("*/", sel_year, sprintf("%02d", sel_month), "*"),
                  output_directory = "temp/",
                  no_directories = T
                )
              )
            if (!inherits(df_ostia, "try-error")) {
              df_ostia <- lapply(df_ostia$files, \(x) x$file_path)
              df_ostia <- try(check_ifdate(df_ostia, sel_year, sel_month,
                                           target = c(ifelse(sel_month == 2, 28, 30), 31)),
                              silent = TRUE)
              if (!inherits(df_ostia, "try-error")) {
                df_ds <- xr$open_mfdataset(unlist(lapply(df_ostia, as.character), recursive = T))
                df_ds <- df_ds$analysed_sst
                df_ds <- df_ds$mean(dim = "time", skipna = T)
                df_ds <- df_ds$to_dataset()
                df_ds <- df_ds$expand_dims(
                  time = pd$to_datetime(paste0(sel_year, "-", sprintf("%02d", sel_month), "-01"))
                )
                df_ds$to_netcdf(ostia_ds)
                fs::file_delete(unlist(lapply(df_ostia, as.character), recursive = T))
              }
            }
          }

          if (!inherits(df_ostia, "try-error")) {
            ostia_result <- process_ostia(all_vals, obis_dataset, ostia_ds, sel_month, sel_year, ostia_prod)
          } else {
            ostia_result <- list(status = "failed")
          }

          if (ostia_result$status == "succeeded") {
            all_vals <- ostia_result$result
            st$set(st_cod, c(st$get(st_cod), paste0("ostia_", part)))
          }
          rm(ostia_result)

          if (part == total_parts) {
            if (sum(grepl("ostia", st$get(st_cod))) == total_parts) {
              log_df <- log_df |> set_success(sel_year, sel_month, "status_ostia")
            } else {
              log_df <- log_df |> set_failed(sel_year, sel_month, "status_ostia")
            }
            if (file.exists(ostia_ds)) fs::file_delete(ostia_ds)
          }
        } else {
          log_df <- log_df |> set_unavailable(sel_year, sel_month, "status_ostia")
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
          if (total_parts > 1) {
            write_parquet(all_vals, 
                          file.path(outfolder_final,
                            paste0(filename, "_year=", sel_year, "_month=", sel_month, "_part=", part, ".parquet")))
          } else {
            write_parquet(all_vals, 
                          file.path(outfolder_final,
                            paste0(filename, "_year=", sel_year, "_month=", sel_month, ".parquet")))
          }
        }

        if (part == total_parts) {
          st$set(st_cod, c(st$get(st_cod), "done"))
          log_df[log_df$year == sel_year & log_df$month == sel_month, "status_general"] <- "concluded"
        }
      }
    } else {
      if (st$exists(st_cod)) {
        catn("Already done")
      } else {
        log_df[log_df$year == sel_year & log_df$month == sel_month, "status_general"] <- "skipped"
      }
    }
    catn("Month log:")
    print(log_df[log_df$year == sel_year & log_df$month == sel_month, ])
    catn("========================")
    st$set("log", log_df)
  }

  if (file.exists(file.path(outfolder, paste0("obisdata_year=", sel_year, ".parquet")))) {
    fs::file_delete(file.path(outfolder, paste0("obisdata_year=", sel_year, ".parquet")))
  }
}

fs::dir_create("logs")
log_df <- st$get("log")
write.csv(log_df, paste0("logs/log_", format(Sys.Date(), "%Y%m%d"), ".csv"), row.names = F)

# END