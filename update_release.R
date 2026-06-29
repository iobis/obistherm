######################### OBIS sea temperature dataset ########################
# June 2026
# Author: Silas C. Principe
# Contact: s.principe@unesco.org
#
########################## Update release information ##########################

release <- readLines("RELEASE.md", warn = FALSE)

gstring <- function(what, text) {
    start <- grep(paste("<!--", what, "-->"), text)
    end <- grep(paste("<!-- END", what, "-->"), text)
    return(c(start, end))
}

release <- c(
    release[1:gstring("DATE", release)[1]],
    paste("Version concluded on", paste0("**", as.character(Sys.Date()), "**")),
    release[gstring("DATE", release)[2]:length(release)]
)

find_latest <- function(dir = "data", pattern = "species_functional_attr_\\d{4}-\\d{2}-\\d{2}\\.parquet$") {
    files <- list.files(
        dir,
        pattern = pattern,
        full.names = TRUE
    )
    if (length(files) == 0) {
        return(NULL)
    }
    files[order(file.mtime(files), decreasing = TRUE)][1]
}

lf_habitat <- find_latest() |>
    (\(x) {
        gsub("\\.parquet", "", gsub("*.*_attr_", "", x))
    })()

release <- c(
    release[1:gstring("HABITAT", release)[1]],
    paste("Habitat zone information updated on", paste0("**", lf_habitat, "**")),
    release[gstring("HABITAT", release)[2]:length(release)]
)

settings <- yaml::read_yaml("settings.yml", readLines.warn = FALSE)
outfolder_final <- settings$outfolder_final

get_stats <- function(outf) {
    require(duckdb())

    con <- dbConnect(duckdb())

    outf <- file.path(outf, "*.parquet")

    dbSendQuery(con, glue::glue("create view obis as select * from read_parquet('{outf}');"))

    count <- dbGetQuery(con, "select count(*) as total, absence from obis group by absence;")
    per_group <- dbGetQuery(con, 'select count(*) as total, phylum from obis where not absence group by phylum;')
    per_sat <- dbGetQuery(con,
    "
    select
        absence,
        count(surfaceTemperature) as GLORYS,
        count(coraltempSST) as CoralTemp,
        count(murSST) as MUR,
        count(ostiaSST) as OSTIA,
    from obis
    group by absence;
    ")

    dbDisconnect(con)

    list(
        total = count,
        per_group = per_group,
        products = per_sat
    )
}

stats <- get_stats(outfolder_final)

stats$products[,1] <- ifelse(stats$products[,1], "Absence", "Presence")
colnames(stats$products)[1] <- "Record type"
for (i in 2:ncol(stats$products)) stats$products[,i] <- format(stats$products[,i], big.mark = ",", scientific = FALSE)

colnames(stats$per_group) <- c("Records", "Phylum")
stats$per_group <- stats$per_group[,c(2,1)]
stats$per_group <- stats$per_group[order(stats$per_group$Records, decreasing = TRUE),]
stats$per_group[,2] <- format(stats$per_group[,2], big.mark = ",", scientific = FALSE)

release <- c(
    release[1:gstring("DS INFO", release)[1]],
    paste(
        "## Dataset info \n### Number of records\n",
        paste("Presence records:", format(stats$total$total[!stats$total$absence], big.mark = ",", scientific = FALSE)),
        paste("Absence records:", format(stats$total$total[stats$total$absence], big.mark = ",", scientific = FALSE)),
        "\n### Records per product \n",
        paste(
            knitr::kable(stats$products),
            collapse = "\n"
        ),
        "\n### Presence records per phylum",
        paste(
            knitr::kable(stats$per_group),
            collapse = "\n"
        ),
        sep = "\n", collapse = "\n"
    ),
    release[gstring("DS INFO", release)[2]:length(release)]
)

# st <- storr::storr_rds("control_storr")
# log_ds <- st$get("log")
log_ds <- find_latest("logs", "log") |>
    read.csv()

log_ds <- log_ds[log_ds$status_general != "skipped", ]

colnames(log_ds) <- c(
    "Year", "Month", "GLORYS", "CoralTemp", "MUR", "OSTIA", "General"
)

which_failed <- apply(log_ds, 1, \(x) {
    any("failed_extract" %in% x[3:6])
})
which_failed <- which(which_failed)

if (length(which_failed) > 0) {
    f <- paste("Failed extractions:", length(which_failed))
    f <- c(f, paste("On year/month:", paste(log_ds$Year[which_failed], log_ds$Month[which_failed], sep = "/", collapse = ", ")))
} else {
    f <- "No extractions failed."
}

release <- c(
    release[1:gstring("LOG INFO", release)[1]],
    "## Log info\n",
    f, "\n",
    paste(
        knitr::kable(log_ds),
        collapse = "\n"
    ),
    release[gstring("LOG INFO", release)[2]:length(release)]
)

writeLines(release, "RELEASE.md")
