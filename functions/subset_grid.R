#' Subset grid
#'
#' @param dataset dataset to subset
#' @param grid_resolution two values vector, with resolution in degrees
#' @param crd_name coordinate names
#' 
#' @return dataset with a column workID
#' @export
#'
#' @examples
#' \dontrun{
#' ds <- subset_grid(ds)
#' }
#'
subset_grid <- function(dataset, grid_res = c(8, 8), 
                        crd_name = c("decimalLongitude", "decimalLatitude")) {
  
  grid <- terra::rast(nrow = grid_res[2], ncol = grid_res[1])
  
  grid[] <- seq_len(terra::ncell(grid))
  
  id <- terra::extract(grid, dataset[,crd_name], ID = F)
  
  dataset$workID <- id[,1]
  
  return(dataset)  
}

