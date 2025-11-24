import pandas as pd
import geopandas as gpd
import pyarrow.parquet as pq

def process_file(path):
    """
    Convert file to GeoParquet

    Args:
        path: target Parquet file
        h3_resolutions: list of H3 resolutions

    Returns:
        file saved

    Depends:
        pandas, geopandas, pyarrow
    """

    df = pd.read_parquet(path)

    gdf = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df.decimalLongitude, df.decimalLatitude),
        crs="EPSG:4326"
    )

    gdf.to_parquet(path, index=False)

    return True