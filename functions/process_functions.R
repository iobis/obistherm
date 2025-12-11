#' Process GLORYS file
#'
#' @param results_df table to hold results
#' @param obis_dataset the dataset from OBIS
#' @param outf_temp_glorys GLORYS file
#' @param sel_month selected month
#' @param sel_year selected year
#'
#' @return list containing the results_df (processed) and status
#' @export
#' 
#' @details
#' Depends on functions loaded in the main code
#'
#' @examples
#' \dontrun{
#'  glorys_result <- process_glorys(results_df, obis_dataset, outf_temp_glorys, sel_month, sel_year)
#'  if (glorys_result$status == "failed") stop("Failed")
#' }
#'
process_glorys <- function(
    results_df,
    obis_dataset,
    outf_temp_glorys,
    sel_month,
    sel_year
) {
    cat("Processing GLORYS\n")

    final_result <- list(
        result = NULL,
        status = NULL
    )

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

    coords_na <- extract_from_nc(outf_temp_glorys, "thetao", obis_dataset[, 1:2], depth = 0L)
    coords_na <- which(is.na(coords_na[, "value"]))

    if (length(coords_na) > 0) {
        new_coords <- get_nearby(outf_temp_glorys, "thetao", obis_dataset[coords_na, 1:2], mode = "25", depth = 0L)
        na_to_rm <- which(is.na(new_coords[, "value"]))
        na_approx <- which(!is.na(new_coords[, "value"]))

        obis_dataset$decimalLongitude[coords_na[na_approx]] <- new_coords$new_lon[na_approx]
        obis_dataset$decimalLatitude[coords_na[na_approx]] <- new_coords$new_lat[na_approx]
        obis_dataset$to_remove_nacord[coords_na[na_to_rm]] <- TRUE
        obis_dataset$to_remove_naaprox[coords_na[na_approx]] <- TRUE
    }

    ds_temp_glorys <- xr$open_dataset(outf_temp_glorys, chunks = "auto")
    valid_depths <- get_depths(
        ds_temp_glorys, obis_dataset,
        paste(sel_year, sprintf("%02d", sel_month), "01", sep = "-")
    )

    valid_depths <- valid_depths |>
        mutate(
            to_remove_ddeep = ifelse(is.na(depth_deep), TRUE, FALSE),
            to_remove_dmid = ifelse(is.na(depth_mid), TRUE, FALSE)
        ) |>
        mutate(
            depth_deep = ifelse(is.na(depth_deep), 0, depth_deep),
            depth_mid = ifelse(is.na(depth_mid), 0, depth_mid)
        )

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

        glorys_data <- left_join(glorys_data, glorys_data_depths, by = "temp_ID")
        glorys_data <- tibble::as_tibble(apply(glorys_data, 2, function(x) {
            x[is.nan(x)] <- NA
            x
        }))

        results_df$surfaceTemperature[results_df$temp_ID == glorys_data$temp_ID] <- glorys_data$depth_surface
        results_df$midTemperature[results_df$temp_ID == glorys_data$temp_ID] <- glorys_data$depth_mid
        results_df$deepTemperature[results_df$temp_ID == glorys_data$temp_ID] <- glorys_data$depth_deep
        results_df$bottomTemperature[results_df$temp_ID == glorys_data$temp_ID] <- glorys_data$depth_bottom

        results_df$minimumDepthTemperature[results_df$temp_ID == glorys_data$temp_ID] <- glorys_data$depth_min
        results_df$maximumDepthTemperature[results_df$temp_ID == glorys_data$temp_ID] <- glorys_data$depth_max

        results_df$midDepth[results_df$temp_ID == glorys_data$temp_ID] <- glorys_data$depth_mid_depth
        results_df$deepDepth[results_df$temp_ID == glorys_data$temp_ID] <- glorys_data$depth_deep_depth

        results_df$minimumDepthClosestDepth[results_df$temp_ID == glorys_data$temp_ID] <- glorys_data$depth_min_depth
        results_df$maximumDepthClosestDepth[results_df$temp_ID == glorys_data$temp_ID] <- glorys_data$depth_max_depth

        results_df$minimumDepthTemperature[to_remove$to_remove_dmin] <- NA
        results_df$minimumDepthClosestDepth[to_remove$to_remove_dmin] <- NA

        results_df$maximumDepthTemperature[to_remove$to_remove_dmax] <- NA
        results_df$maximumDepthClosestDepth[to_remove$to_remove_dmax] <- NA

        results_df$midTemperature[to_remove$to_remove_dmid] <- NA
        results_df$midDepth[to_remove$to_remove_dmid] <- NA

        results_df$deepTemperature[to_remove$to_remove_ddeep] <- NA
        results_df$deepDepth[to_remove$to_remove_ddeep] <- NA

        results_df[to_remove$to_remove_nacord, c(
            "surfaceTemperature", "midTemperature",
            "deepTemperature", "midDepth", "deepDepth", "minimumDepthTemperature", "maximumDepthTemperature",
            "minimumDepthClosestDepth", "maximumDepthClosestDepth"
        )] <- NA

        results_df$flag[to_remove$to_remove_naaprox] <- results_df$flag[to_remove$to_remove_naaprox] + 2

        depth_diff_min <- check_depth_diff(obis_dataset$depth_min, results_df$minimumDepthClosestDepth)
        depth_diff_max <- check_depth_diff(obis_dataset$depth_max, results_df$maximumDepthClosestDepth)

        results_df$flag[depth_diff_min] <- results_df$flag[depth_diff_min] + 4
        results_df$flag[depth_diff_max] <- results_df$flag[depth_diff_max] + 8

        final_result$result <- results_df
        final_result$status <- "succeeded"
    } else {
        final_result$status <- "failed"
    }
    return(final_result)
}



#' Process CoralTemp file
#'
#' @param results_df table to hold results
#' @param obis_dataset the dataset from OBIS
#' @param coraltemp_ds CoralTemp file
#' @param sel_month selected month
#' @param sel_year selected year
#'
#' @return list containing the results_df (processed) and status
#' @export
#' 
#' @details
#' Depends on functions loaded in the main code
#'
#' @examples
#' \dontrun{
#'  coraltemp_result <- process_coraltemp(results_df, obis_dataset, coraltemp_ds, sel_month, sel_year)
#'  if (coraltemp_result$status == "failed") stop("Failed")
#' }
#'
process_coraltemp <- function(
    results_df,
    obis_dataset,
    coraltemp_ds,
    sel_month,
    sel_year
) {
    cat("Extracting CoralTemp\n")

    final_result <- list(
        result = NULL,
        status = NULL
    )

    success <- try(download_temp("coraltemp", coraltemp_ds, obis_dataset, sel_month, sel_year))

    if (!inherits(success, "try-error")) {
        cat("Processing CoralTemp\n")

        success$value[is.nan(success$value)] <- NA
        results_df$coraltempSST <- success$value

        na_to_solve <- which(is.na(success[, "value"]))

        if (length(na_to_solve) > 0) {
            nearby_ct <- get_nearby(coraltemp_ds, "sea_surface_temperature",
                                obis_dataset[na_to_solve, c("decimalLongitude", "decimalLatitude")],
                                mode = "25", depth = NULL,
                                date = paste(sel_year, sprintf("%02d", sel_month), "15", sep = "-"))

            results_df$coraltempSST[na_to_solve] <- nearby_ct$value
            results_df$flag[na_to_solve[!is.na(nearby_ct$value)]] <-
                results_df$flag[na_to_solve[!is.na(nearby_ct$value)]] + 16
        }

        final_result$result <- results_df
        final_result$status <- "succeeded"
    } else {
        final_result$status <- "failed"
    }
    return(final_result)
}


#' Process MUR file
#'
#' @param results_df table to hold results
#' @param obis_dataset the dataset from OBIS
#' @param mur_ds MUR file
#' @param sel_month selected month
#' @param sel_year selected year
#'
#' @return list containing the results_df (processed) and status
#' @export
#' 
#' @details
#' Depends on functions loaded in the main code
#'
#' @examples
#' \dontrun{
#'  mur_result <- process_mur(results_df, obis_dataset, mur_ds, sel_month, sel_year)
#'  if (mur_result$status == "failed") stop("Failed")
#' }
#'
process_mur <- function(
    results_df,
    obis_dataset,
    mur_ds,
    sel_month,
    sel_year
) {
    cat("Extracting MUR\n")

    final_result <- list(
        result = NULL,
        status = NULL
    )

    success <- try(download_temp("mur", mur_ds, obis_dataset, sel_month, sel_year))

    if (!inherits(success, "try-error")) {
        cat("Processing MUR\n")

        success$value[is.nan(success$value)] <- NA
        results_df$murSST <- success$value

        na_to_solve <- which(is.na(success[, "value"]))

        if (length(na_to_solve) > 0) {
            nearby_mur <- get_nearby(mur_ds, "sst",
                                    obis_dataset[na_to_solve, c("decimalLongitude", "decimalLatitude")],
                                    mode = "25", depth = NULL,
                                    date = paste(sel_year, sprintf("%02d", sel_month), "15", sep = "-"))

            results_df$murSST[na_to_solve] <- nearby_mur$value
            results_df$flag[na_to_solve[!is.na(nearby_mur$value)]] <-
                results_df$flag[na_to_solve[!is.na(nearby_mur$value)]] + 32
        }

        final_result$result <- results_df
        final_result$status <- "succeeded"
    } else {
        final_result$status <- "failed"
    }
    return(final_result)
}


#' Process OSTIA file
#'
#' @param results_df table to hold results
#' @param obis_dataset the dataset from OBIS
#' @param ostia_ds OSTIA file
#' @param sel_month selected month
#' @param sel_year selected year
#' @param ostia_prod name of the OSTIA product
#'
#' @return list containing the results_df (processed) and status
#' @export
#' 
#' @details
#' Depends on functions loaded in the main code
#'
#' @examples
#' \dontrun{
#'  ostia_result <- process_ostia(results_df, obis_dataset, ostia_ds_open, sel_month, sel_year)
#'  if (ostia_result$status == "failed") stop("Failed")
#' }
#'
process_ostia <- function(
    results_df,
    obis_dataset,
    ostia_ds,
    sel_month,
    sel_year,
    ostia_prod
) {
    cat("Extracting OSTIA\n")

    final_result <- list(
        result = NULL,
        status = NULL
    )

    ostia_ds_open <- xr$open_dataset(ostia_ds)
    success <- try(download_temp("ostia", ostia_ds_open, obis_dataset, sel_month, sel_year))

    if (!inherits(success, "try-error")) {
        cat("Processing OSTIA\n")

        success$value[is.nan(success$value)] <- NA
        results_df$ostiaSST <- success$value

        na_to_solve <- which(is.na(success[, "value"]))

        if (length(na_to_solve) > 0) {
            nearby_ostia <- get_nearby(ostia_ds, "analysed_sst",
                                    obis_dataset[na_to_solve, c("decimalLongitude", "decimalLatitude")],
                                    mode = "25", depth = NULL,
                                    date = paste(sel_year, sprintf("%02d", sel_month), "01", sep = "-"))
            nearby_ostia$value <- nearby_ostia$value - 273.15

            results_df$ostiaSST[na_to_solve] <- nearby_ostia$value
            results_df$flag[na_to_solve[!is.na(nearby_ostia$value)]] <-
                results_df$flag[na_to_solve[!is.na(nearby_ostia$value)]] + 64
        }

        results_df$ostiaProduct[!is.na(results_df$ostiaSST)] <- ostia_prod

        final_result$result <- results_df
        final_result$status <- "succeeded"
    } else {
        final_result$status <- "failed"
    }
    return(final_result)
}
