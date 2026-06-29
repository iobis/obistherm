######################### OBIS sea temperature dataset ########################
# January of 2024
# Author: Silas C. Principe
# Contact: s.principe@unesco.org
#
################################ Data aggregation ##############################

library(arrow)
library(dplyr)
library(reticulate)
source("functions/utils.R")
reticulate::source_python("functions/convert_geoarrow.py")
settings <- yaml::read_yaml("settings.yml", readLines.warn = FALSE)
input_folder <- settings$outfolder_final
output_folder <- settings$outfolder_aggregated
fs::dir_create(output_folder)

species_info <- build_species_info(input_folder) |>
  arrow::read_parquet()

aggregate_data <- function(input_folder, output_folder,
                           h3_resolutions = c(7L),
                           species_attrs) {
  
  cat("Aggregating as Parquet dataset (hive format)\n")

  file_schema <- schema(read_parquet(list.files(input_folder, full.names = T)[1]))
  file_schema$AphiaID <- arrow::int32()
  file_schema$minimumDepthInMeters <- double()
  file_schema$maximumDepthInMeters <- double()
  file_schema$surfaceTemperature <- double()
  file_schema$midTemperature <- double()
  file_schema$deepTemperature <- double()
  file_schema$bottomTemperature <- double()
  file_schema$midDepth <- double()
  file_schema$deepDepth <- double()
  file_schema$minimumDepthTemperature <- double()
  file_schema$maximumDepthTemperature <- double()
  file_schema$minimumDepthClosestDepth <- double()
  file_schema$maximumDepthClosestDepth <- double()
  file_schema$coraltempSST <- double()
  file_schema$murSST <- double()
  file_schema$ostiaSST <- double()
  
  ds <- open_dataset(input_folder,
                     schema = file_schema)
  
  ds |>
    select(-extractedDateStart, -extractedDateMid, -extractedDateEnd, -date_year, -dropped) |>
    rename(year = extractedDateYear, month = extractedDateMonth) |>
    mutate(year = as.integer(year), 
           month = as.integer(month),
           surfaceTemperature = round(surfaceTemperature, 2),
           midTemperature = round(midTemperature, 2),
           deepTemperature = round(deepTemperature, 2),
           bottomTemperature = round(bottomTemperature, 2),
           minimumDepthTemperature = round(minimumDepthTemperature, 2),
           maximumDepthTemperature = round(maximumDepthTemperature, 2),
           coraltempSST = round(coraltempSST, 2),
           murSST = round(murSST, 2),
           ostiaSST = round(ostiaSST, 2)) |>
    group_by(year) |>
    write_dataset(path = output_folder)

  rm(ds)
  
  cat("Converting to GeoArrow and adding H3\n")

  species_attrs <- species_attrs |>
    select(AphiaID, adultFunctionalGroup = functional_group, fg_flag = flag) |>
    filter(grepl("FG life stage: adult", fg_flag)) |>
    mutate(fg_flag = gsub("; FG life stage: adult|FG life stage: adult", "", fg_flag))

  tf <- list.files(output_folder, recursive = T, full.names = T)

  pb <- progress::progress_bar$new(total = length(tf))

  for (id in seq_along(tf)) {
    pb$tick()

    tf_content <- read_parquet(tf[id])

    depth_flag <- "Difference between minimumDepthInMeters and maximumDepthInMeters is higher than 100m"
    no_depth_flag <- "No depth information available for the record"

    tf_content2 <- tf_content |>
      mutate(obistherm_flags = decode_flag(obistherm_flags)) |>
      left_join(species_attrs, by = "AphiaID") |>
      mutate(obistherm_flags = ifelse(
        is.na(fg_flag),
        obistherm_flags,
        ifelse(
          is.na(obistherm_flags),
          fg_flag,
          paste(obistherm_flags, fg_flag, sep = "; ")
        )
      )) |>
      select(-fg_flag) |>
      mutate(
        depth_case = case_when(
          is.na(minimumDepthInMeters) & is.na(maximumDepthInMeters) ~ "all_na_case",
          is.na(minimumDepthInMeters) | is.na(maximumDepthInMeters) ~ "na_case",
          abs(maximumDepthInMeters - minimumDepthInMeters) > 100 ~ "higher_case",
          TRUE ~ "normal_case"
        )
      ) |>
      mutate(
        obistherm_flags = case_when(
          is.na(obistherm_flags) ~ NA,
          depth_case == "na_case" ~ obistherm_flags,
          depth_case == "higher_case" ~ paste(obistherm_flags, depth_flag, sep = "; "),
          depth_case == "all_na_case" ~ paste(obistherm_flags, no_depth_flag, sep = "; "),
          TRUE ~ obistherm_flags
        )
      ) |>
      select(-depth_case) |>
      mutate(AphiaID = as.integer(AphiaID))
    
    for (hr in h3_resolutions) {
      batches <- split(seq_len(nrow(tf_content)), ceiling(seq_len(nrow(tf_content)) / 10000))

      cell_values <- lapply(batches, function(bt) {
        suppressMessages(h3jsr::point_to_cell(tf_content[bt, c("decimalLongitude", "decimalLatitude")], res = hr))
      })
      cell_values <- unlist(cell_values, use.names = F)

      tf_content[[paste0("h3_", hr)]] <- cell_values
    }

    tf_content |> write_parquet(tf[id])
    rm(tf_content, batches)

    proc_res <- process_file(tf[id])

    if (!proc_res) stop("Failed processing file ", paste0("\033[47m", tf[id], "\033[49m"))
  }
  
  cat("Aggregation concluded \n")
    
  return(invisible(NULL))
  
}

# Aggregate
aggregate_data(input_folder, output_folder, species_attrs = species_info)

### END