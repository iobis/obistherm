######################### OBIS sea temperature dataset ########################
# June 2026
# Author: Silas C. Principe
# Contact: s.principe@unesco.org
#
############ Update release information and README figures######################

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

# Generate README figures
source("functions/retrieve_data.R")
library(ggplot2)
out_f <- "/Volumes/OBIS2/results_final/"

ocy <- retrieve_data(scientificname = c(
    "Leptuca thayeri", "Minuca rapax"
), datasource = out_f) # Change here with your data source or NULL
# to use the S3 access point

ocy_proc <- ocy |>
    mutate(date = as.Date(paste(year, month, "01", sep = "-"))) |>
    rename(glorysSST = surfaceTemperature) |>
    tidyr::pivot_longer(cols = ends_with("SST"),
                        names_to = "SSTsource", values_to = "sst")

ggplot(ocy_proc) +
    geom_point(aes(x = date, y = sst, color = SSTsource), alpha = .5) +
    geom_smooth(aes(x = date, y = sst, group = SSTsource),
                method = "lm", color = "grey30") +
    scale_color_manual(values = c("#11b5ae", "#4046ca", "#f68512", "#de3c82")) +
    xlab(NULL) + ylab("Temperature (°C)") +
    theme_light() +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(face = "italic")
    ) +
    facet_grid(SSTsource ~ species)

ggsave("images/fiddlers.png", width = 13, height = 10)

library(arrow)
library(dplyr)
library(ggplot2)
library(sf)

ds <- open_dataset(out_f)

acanthurus <- ds |>
    filter(!absence) |>
    select(species, coraltempSST, geometry, h3_7, year, month) |>
    filter(species == "Acanthurus coeruleus") |>
    filter(!is.na(coraltempSST)) |>
    filter(year == 2014)

acanthurus <- sfarrow::read_sf_dataset(acanthurus)
st_crs(acanthurus) <- 4326

period <- data.frame(
    period = rep(c("T1", "T2", "T3", "T4"), each = 3),
    month = 1:12
)

acanthurus <- left_join(acanthurus, period)

ggplot() +
    geom_sf(data = rnaturalearth::ne_countries(returnclass = "sf"),
            color = "gray80", fill = "gray80") +
    geom_sf(data = acanthurus, aes(color = coraltempSST, shape = period), alpha = .5, size = 2) +
    scale_color_viridis_c() +
    coord_sf(xlim = c(-95, -30), ylim = c(-30, 30)) +
    ggtitle("Acanthurus choeruleus", "Year: 2014 - Product: CoralTemp (NOAA)") +
    theme_light()  +
    theme(plot.title = element_text(face = "italic"), legend.position = "bottom")

ggsave("images/acanthurus.png")

acanthurus_agg <- acanthurus |>
    sf::st_drop_geometry() |>
    mutate(h3_4 = h3jsr::get_parent(h3_7, res = 4)) |>
    group_by(h3_4) |>
    summarise(mean_sst = mean(coraltempSST))

acanthurus_agg_pol <- h3jsr::cell_to_polygon(acanthurus_agg$h3_4)
acanthurus_agg_pol <- st_as_sf(acanthurus_agg_pol)
acanthurus_agg_pol <- bind_cols(acanthurus_agg_pol, acanthurus_agg)

ggplot() +
    geom_sf(data = rnaturalearth::ne_countries(returnclass = "sf"),
            color = "gray80", fill = "gray80") +
    geom_sf(data = acanthurus_agg_pol, aes(fill = mean_sst)) +
    scale_fill_viridis_c() +
    coord_sf(xlim = c(-95, -30), ylim = c(-30, 30)) +
    ggtitle("Acanthurus choeruleus", "Year: 2014 - Product: CoralTemp (NOAA)") +
    theme_light()  +
    theme(plot.title = element_text(face = "italic"), legend.position = "bottom")

ggsave("images/acanthurus_agg.png")

reticulate::py_install(c("matplotlib", "contextily", "palettable", "lonboard", "seaborn"))
reticulate::py_run_string(glue::glue("
import geopandas as gpd
import pyarrow.dataset as ds
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import contextily as ctx
from palettable.matplotlib import Viridis_10

dataset = ds.dataset('{out_f}', format='parquet', partitioning='hive')
cols = ['geometry', 'species', 'coraltempSST', 'year']
table = dataset.to_table(
    filter=(ds.field('species') == 'Acanthurus chirurgus') &
           (ds.field('year') == 2000) &
           ~ds.field('absence'),
    columns=cols
)
gdf = gpd.GeoDataFrame.from_arrow(table)

values = gdf['coraltempSST']
norm = (values - values.min()) / (values.max() - values.min())

cmap = Viridis_10.mpl_colormap  # Sequential colormap
norm = mcolors.Normalize(vmin=values.min(), vmax=values.max())

gdf = gdf.to_crs(epsg=3857)

fig, ax = plt.subplots(figsize=(10, 8))
sc = gdf.plot(column='coraltempSST', 
              cmap=cmap, 
              norm=norm, 
              markersize=5, 
              legend=False, 
              ax=ax, 
              alpha=0.7)

ctx.add_basemap(ax, source=ctx.providers.CartoDB.Positron, zoom=5)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([]) 
cbar = fig.colorbar(sm, ax=ax, orientation='vertical')
cbar.set_label('Coral Temperature SST (°C)')

ax.set_xticks([])
ax.set_yticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.set_title('Species: Acanthurus chirurgus - Coral Temperature SST')

plt.savefig('images/species_py.png', dpi=150, bbox_inches='tight')
plt.close()
"))

reticulate::py_run_string(glue::glue('
import geopandas as gpd
import pyarrow.dataset as ds
import lonboard
from lonboard.colormap import apply_continuous_cmap
from palettable.colorbrewer.diverging import BrBG_10

dataset = ds.dataset("{out_f}", format="parquet", partitioning="hive")
cols = ["geometry", "species", "coraltempSST", "year"]
table = dataset.to_table(
    filter=(ds.field("species") == "Acanthurus chirurgus") &
           (ds.field("year") == 2000) &
           ~ds.field("absence"),
    columns=cols
)
gdf = gpd.GeoDataFrame.from_arrow(table)

values = gdf["coraltempSST"]
normalized_values = (values - values.min()) / (values.max() - values.min())

point_layer = lonboard.ScatterplotLayer.from_geopandas(gdf)
point_layer.get_radius = 10000
point_layer.radius_max_pixels = 2
point_layer.get_fill_color = apply_continuous_cmap(normalized_values, BrBG_10, alpha=0.7)

map_ = lonboard.Map(point_layer, view_state={{"longitude": -75, "latitude": 15, "zoom": 4}})
with open("temp_lonboard_map.html", "w") as f:
    f.write(map_.as_html().data)
'))

if (!requireNamespace("webshot2", quietly = TRUE)) {
    install.packages("webshot2", repos = "https://cloud.r-project.org/")
}
webshot2::webshot("temp_lonboard_map.html", "images/lonboard_py.png", delay = 5)

fs::file_delete("temp_lonboard_map.html")
