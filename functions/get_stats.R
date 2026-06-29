#' Get stats of the final dataset
#'
#' @param final_folder path to final dataset
#' @param save_rds if `TRUE`, save an RDS containing results
#' @param rds_name name to the RDS file. The date will be appended to the start.
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
get_ds_stats <- function(final_folder, save_rds = TRUE, rds_name = "_stats.rds") {
    require(duckdb)
    require(dplyr)
    set.seed(2026)
    con <- dbConnect(duckdb())

    dbSendQuery(con, glue::glue("
    create or replace view obis as
        select *
        from read_parquet('{final_folder}/*/*.parquet', hive_partitioning = true);
    "))

    total <- dbGetQuery(con, 
        glue::glue(
            "select count(*) as total, count(surfaceTemperature) as glorys, 
                count(coraltempSST) as coraltemp, count(murSST) as mur, count(ostiaSST) as ostia
            from obis;
            "
        )
    )

    by_group <- dbGetQuery(con, 
        glue::glue(
            "select phylum, \"order\", class, family, count(*) as total, sum(absence) as n_absence
            from obis
            group by phylum, \"order\", class, family
            "
        )
    )

    phylum <- by_group |> group_by(phylum) |>
        summarise(total = sum(total), total_absence = sum(n_absence))
    order <- by_group |> group_by(order) |>
        summarise(total = sum(total), total_absence = sum(n_absence))
    class <- by_group |> group_by(class) |>
        summarise(total = sum(total), total_absence = sum(n_absence))
    family <- by_group |> group_by(family) |>
        summarise(total = sum(total), total_absence = sum(n_absence))

    by_year <- dbGetQuery(con, 
        glue::glue(
            "select year, count(*) as total
            from obis
            group by year
            order by year;
            "
        )
    )

    number_absences <- dbGetQuery(con, 
        glue::glue(
            "select count(*) as total
            from obis
            where absence
            "
        )
    )

    coverage_prod <- dbGetQuery(con,
        glue::glue(
            "with counts as (
                select
                    ((surfaceTemperature is not null)::int +
                    (coraltempSST is not null)::int +
                    (murSST is not null)::int +
                    (ostiaSST is not null)::int) as n_products
                from obis
            )
            select
                count(*) filter (where n_products = 4) as all_products,
                count(*) filter (where n_products = 1) as one_product
            from counts;"
        )
    )

    coverage_depth <- dbGetQuery(con,
        glue::glue(
            "select count(*) as n_with_depth
            from obis
            where minimumDepthTemperature is not null
            or maximumDepthTemperature is not null;"
        ))

    product_comparison <- dbGetQuery(con,
        "
        select
            decimallongitude,
            decimallatitude,
            surfacetemperature,
            coraltempsst,
            mursst,
            ostiasst
        from obis
        where surfacetemperature is not null
            and coraltempsst is not null
            and mursst is not null
            and ostiasst is not null
        order by random()
        limit 100000;
        "
        )

    product_comparison$dist_coast <- obistools::lookup_xy(
        product_comparison, shoredistance = TRUE, grids = FALSE,
        areas = FALSE
    )[,1]

    product_diff <- apply(product_comparison[,4:6], 2, \(x) x - product_comparison$surfaceTemperature)

    ret_obj <- list(
        total = total,
        by_group = by_group,
        phylum = phylum,
        order = order,
        class = class,
        family = family,
        by_year = by_year,
        absences = number_absences,
        product_coverage = coverage_prod,
        depth_coverage = coverage_depth,
        product_comparison = list(
            full_data = product_comparison,
            differences = product_diff
        )
    )

    if (save_rds) saveRDS(ret_obj, paste0(format(Sys.Date(), "%Y%m%d"), rds_name))

    dbDisconnect(con)

    return(ret_obj)
}

settings <- yaml::read_yaml("settings.yml", readLines.warn = FALSE)

final_folder <- settings$outfolder_aggregated

ds_stats <- get_ds_stats(final_folder)

fm <- function(v) format(v, big.mark=",", scientific=FALSE)

ds_stats$total |> fm()
ds_stats$phylum[order(ds_stats$phylum$total, decreasing = T),] |>
    mutate(diff = total-total_absence) |>
    mutate(total = fm(total), total_absence = fm(total_absence), diff = fm(diff))
ds_stats$absences |> fm()

gp <- function(total, v2) {
    (v2 * 100)/total
}

gp(ds_stats$total$total, ds_stats$absences$total) |> round(1)

ds_stats$product_coverage |> fm()

gp(ds_stats$total$total, ds_stats$product_coverage[[1]]) |> round(1)
gp(ds_stats$total$total, ds_stats$product_coverage[[2]]) |> round(1)

fm(ds_stats$depth_coverage[[1]])
fm(ds_stats$total$total - ds_stats$depth_coverage[[1]])
gp(ds_stats$total$total, ds_stats$depth_coverage[[1]]) |> round(1)

ds_stats$product_comparison[[1]] |> head()

apply(
    ds_stats$product_comparison[[1]][,c("coraltempSST", "murSST", "ostiaSST")], 2,
    \(x) cor(x, ds_stats$product_comparison[[1]]$surfaceTemperature)
) |> round(2)

apply(
    ds_stats$product_comparison[[2]], 2,
    \(x) data.frame(avg = mean(abs(x)), sds = sd(abs(x))) |> round(2)
) 

cor(ds_stats$product_comparison[[1]][,c("surfaceTemperature", "coraltempSST", "murSST", "ostiaSST")]) |>
    min()

prod_coast <- which(ds_stats$product_comparison[[1]]$dist_coast > 0 & ds_stats$product_comparison[[1]]$dist_coast <= 100)

apply(
    ds_stats$product_comparison[[1]][prod_coast,c("coraltempSST", "murSST", "ostiaSST")], 2,
    \(x) cor(x, ds_stats$product_comparison[[1]]$surfaceTemperature[prod_coast])
) |> round(2)

apply(
    ds_stats$product_comparison[[2]][prod_coast,], 2,
    \(x) data.frame(avg = mean(abs(x)), sds = sd(abs(x))) |> round(2)
) 

cor(ds_stats$product_comparison[[1]][prod_coast,c("surfaceTemperature", "coraltempSST", "murSST", "ostiaSST")]) |>
    min()
