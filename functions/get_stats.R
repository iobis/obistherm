#' Get stats of the final dataset
#'
#' @param final_folder path to final dataset
#' @param save_rds if `TRUE`, save an RDS containing results
#'
#' @return list with stats and optionally saved RDS
#' @export
#' 
#' @details 
#' Return basic stats for the total, by year and by group
#'
#' @examples
#' \dontrun{
#' get_ds_stats("results")
#' }
#'
get_ds_stats <- function(final_folder, save_rds = TRUE) {
    require(duckdb)
    require(dplyr)
    con <- dbConnect(duckdb())

    total <- dbGetQuery(con, 
        glue::glue(
            "select count(*) as total, count(surfaceTemperature) as glorys, 
                count(coraltempSST) as coraltemp, count(murSST) as mur, count(ostiaSST) as ostia
            from read_parquet('{final_folder}/*/*.parquet', hive_partitioning = true);
            "
        )
    )

    by_group <- dbGetQuery(con, 
        glue::glue(
            "select phylum, \"order\", class, family, count(*) as total
            from read_parquet('{final_folder}/*/*.parquet', hive_partitioning = true)
            group by phylum, \"order\", class, family
            "
        )
    )

    phylum <- by_group |> group_by(phylum) |> summarise(total = sum(total))
    order <- by_group |> group_by(order) |> summarise(total = sum(total))
    class <- by_group |> group_by(class) |> summarise(total = sum(total))
    family <- by_group |> group_by(family) |> summarise(total = sum(total))

    by_year <- dbGetQuery(con, 
        glue::glue(
            "select year, count(*) as total
            from read_parquet('{final_folder}/*/*.parquet', hive_partitioning = true)
            group by year
            order by year;
            "
        )
    )

    ret_obj <- list(
        total = total,
        by_group = by_group,
        phylum = phylum,
        order = order,
        class = class,
        family = family,
        by_year = by_year
    )

    if (save_rds) saveRDS(ret_obj, paste0(format(Sys.Date(), "%Y%m%d"), "_stats.rds"))

    return(ret_obj)
}

ds_stats <- get_ds_stats()