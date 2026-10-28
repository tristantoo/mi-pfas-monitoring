# Download Michigan EGLE public water supply PFAS results + hexbin polygons.
# Source: EGLE PublicWaterSupplySamplingOpenData FeatureServer (State of Michigan Open Data).
# Run: Rscript R/ingest.R   (writes data/raw/, idempotent; re-run to refresh)

library(jsonlite)

base <- "https://gisagoegle.state.mi.us/arcgis/rest/services/EGLE/PublicWaterSupplySamplingOpenData/FeatureServer"
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)

query <- function(layer, offset, geom = FALSE, fmt = "json") {
  url <- sprintf(
    "%s/%d/query?where=1%%3D1&outFields=*&returnGeometry=%s&orderByFields=OBJECTID&resultOffset=%d&resultRecordCount=2000&f=%s",
    base, layer, tolower(geom), offset, fmt
  )
  fromJSON(url)
}

expected <- fromJSON(sprintf("%s/1/query?where=1%%3D1&returnCountOnly=true&f=json", base))$count
pages <- list()
offset <- 0
repeat {
  page <- query(1, offset)
  if (!length(page$features)) break
  pages[[length(pages) + 1]] <- page$features$attributes
  offset <- offset + nrow(page$features$attributes)
  if (!isTRUE(page$exceededTransferLimit)) break
}
res <- do.call(rbind, pages)
stopifnot(nrow(res) == expected, !anyDuplicated(res$OBJECTID))

# ArcGIS dates are epoch milliseconds (UTC midnight); keep them as ISO dates.
res$SampleDate <- as.Date(as.POSIXct(res$SampleDate / 1000, origin = "1970-01-01", tz = "UTC"))
write.csv(res, "data/raw/pws_results.csv", row.names = FALSE)

hex <- readLines(sprintf(
  "%s/0/query?where=1%%3D1&outFields=HexID&returnGeometry=true&outSR=4326&f=geojson", base
), warn = FALSE)
writeLines(hex, "data/raw/hexbins.geojson")

writeLines(c(
  paste("retrieved:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("source:", base),
  paste("rows:", nrow(res)),
  paste("date_range:", paste(range(res$SampleDate, na.rm = TRUE), collapse = " to ")),
  paste("md5_results:", tools::md5sum("data/raw/pws_results.csv")),
  paste("md5_hexbins:", tools::md5sum("data/raw/hexbins.geojson"))
), "data/raw/MANIFEST.txt")
cat(readLines("data/raw/MANIFEST.txt"), sep = "\n")
