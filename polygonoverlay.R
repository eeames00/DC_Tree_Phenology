## Polygon Overlay

library(sf)
library(data.table)


shp_path <- "C:/Users/valea/OneDrive/Desktop/pheno_sensitivity/UFD_segs_w_species_v2.shp"

file.exists(shp_path)   # should be TRUE

crowns <- st_read(shp_path, quiet = TRUE)
names(crowns)
st_crs(crowns)



# 1) your polygon (GeoJSON string)
gj_txt <- '{"type":"FeatureCollection","crs":{"type":"name","properties":{"name":"EPSG:32618"}},"features":[{"type":"Feature","geometry":{"type":"Polygon","coordinates":[[[317861.78003610391,4312223.0008060876],[317862.73753167409,4308276.8546297532],[329898.60809294786,4308204.8934277985],[329958.74735368416,4312072.2047635317],[317861.78003610391,4312223.0008060876]]]}}]}'
gj_file <- tempfile(fileext = ".geojson")
writeLines(gj_txt, gj_file)
poly <- st_read(gj_file, quiet = TRUE)

# 2) crowns already read as `crowns`
# crowns <- st_read(shp_path, quiet = TRUE)

# 3) intersect (crowns are polygons, so intersects is right)
in_poly <- lengths(st_intersects(crowns, poly)) > 0
uids_in_poly <- crowns$UID[in_poly]

length(uids_in_poly)
head(uids_in_poly)

# 4) filter your NDVI tables
all_dt_poly  <- all_dt[UID %in% uids_in_poly]
tree_ts_poly <- tree_ts[UID %in% uids_in_poly]

nrow(all_dt_poly)
uniqueN(all_dt_poly$UID)


plot(st_geometry(poly), border = "red", lwd = 2)
plot(st_geometry(crowns), add = TRUE, border = NA, col = rgb(0,0,0,0.05))
plot(st_geometry(crowns[in_poly, ]), add = TRUE, border = "blue", lwd = 1)



out_dir <- "C:/Users/valea/OneDrive/Desktop/pheno_sensitivity"

fwrite(all_dt_poly,  file.path(out_dir, "compiled_crown_ndvi_polygon_subset.csv"))
saveRDS(all_dt_poly, file.path(out_dir, "compiled_crown_ndvi_polygon_subset.rds"))
saveRDS(tree_ts_poly, file.path(out_dir, "tree_timeseries_daily_polygon_subset.rds"))

