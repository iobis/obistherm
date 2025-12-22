# deprecated, use nearby_from_nc.R
#' Get nearby points based on the XY coordinates
#'
#' @param coords decimal longitude and decimal latitude (in this order)
#' @param tlayer the SpatRaster
#' @param mode queen for "queen mode" or any other value for a 5x5 matrix
#'
#' @return nearby cells
#' @export
#'
#' @examples
#' \dontrun{
#' get_nearby_xy(c(-5, 6), my_rast, mode = "other")
#' }
#'
get_nearby_xy <- function(coords, tlayer, mode = "queen") {
  if (nrow(coords) > 0) {
    tcell <- cellFromXY(tlayer, as.data.frame(coords))
    if (mode == "queen") {
      adj_m <- "queen"
    } else {
      adj_m <- matrix(c(rep(1, 12), 0, rep(1, 12)), 5, 5)
    }
    adj <- adjacent(tlayer, cells = tcell, adj_m, include = T)
    
    result <- apply(adj, 1, function(x) {
      
      ext_vals <- terra::extract(tlayer, x)
      
      if (any(!is.na(ext_vals[,1]))) {
        if (sum(!is.na(ext_vals[,1])) > 1) {
          dists <- terra::nearest(vect(xyFromCell(tlayer, x[1]), crs = "EPSG:4326"),
                                  vect(xyFromCell(tlayer, x[!is.na(ext_vals[,1])]), crs = "EPSG:4326"))
          dists <- geom(dists)
          ncell <- cellFromXY(tlayer, data.frame(x = dists[,"x"], y = dists[,"y"]))
        } else {
          ncell <- x[!is.na(ext_vals[,1])]
        }
        to_return <- data.frame(xyFromCell(tlayer, ncell))
      } else {
        to_return <- data.frame(x = NA, y = NA)
      }
      
      return(to_return)
      
    })
    
    return(do.call("rbind", result)) # to change nan for NA
  } else {
    return(NULL)
  }
}