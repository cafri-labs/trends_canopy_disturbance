library(terra)
library(dplyr)
library(landscapemetrics)

#code for assessment of general disturbance trends


######## load data ########
#cafri crs is  "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
crs <- labrador.client::get_cafri_crs()
#requires annual loss brick, code for creation can be found here: https://github.com/cafri-labs/lt_tuning/blob/main/scripts/create_loss_layer.R
raster_stack <- rast("rt1_1985_2022_all_losses.tif")

#name of tuning
tuning <- "rt_1"

#adirondack forest preserve lands
afp_polygons <- vect("afp_aquisition_dates_final.gpkg")
afp_polygons <- project(afp_polygons, crs)
afp_df <- as.data.frame(afp_polygons, geom = "WKT")

raster_stack <- crop(raster_stack, afp_polygons)

#mask all water out using NYS hydrology boundaries. file can be found here: https://data.gis.ny.gov/maps/sharegisny::nys-hydrography/about
hydro <- vect("../GIS_Files/hydrography.gdb")
hydro <- project(hydro, crs)

raster_stack <- mask(raster_stack, hydro, inverse = TRUE)


######## create masked raster stack, apply magnitude threshold #########
#masking by afp areas annually
masked_rasters <- lapply(1:nlyr(raster_stack), \(i) {
  message(paste0("Working on iteration ", i))
  year_polys <- afp_df %>% dplyr::filter(Year <= (1989 +i))
  year_polys <- vect(year_polys, geom = "geometry")
  
  raster_layer <- raster_stack[[i]]
  
  raster_layer <- mask(raster_layer, year_polys)
  raster_layer[raster_layer < 50] <- NA #disturbance magnitude threshold
  return(raster_layer)
})

masked_rasters <- rast(masked_rasters) 
writeRaster(masked_rasters, paste0(tuning, "_masked_rasters.tif"), overwrite = TRUE)

################# create classes ##############
#this is for the patch metrics, all disturbance patches from a single year are 1 class
class_rasters <- lapply(1:nlyr(masked_rasters), \(i) {
  message(paste0("Working on iteration ", i))
  raster <- masked_rasters[[i]]
  raster[!is.na(raster)] <- 1989 + i
  return(raster)
})

class_rasters <- rast(class_rasters)
writeRaster(class_rasters, paste0(tuning, "_class_rasters.tif"), overwrite = TRUE)
#class_rasters <- rast("minor_disturbances_results/minor_class_stack.tif")

################## create disturbance patch raster ##############
patch_rasters <- get_patches(class_rasters, class = "all", directions = 8, to_disk = TRUE)

patch_rasters <- lapply(1:length(patch_rasters), \(i) {
  patch_rasters[[i]][[1]]
})

patch_rasters <- rast(patch_rasters)
names(patch_rasters) <- c(1990:2022)

################## calculate patch metrics ###############
options_landscapemetrics(to_disk = TRUE)
patch_metrics <- calculate_lsm(class_rasters, level = "patch", directions = 8, neighbourhood = 8 ,progress = TRUE)

#combine patches and patch metrics
patches <- lapply(1:nlyr(patch_rasters), \(i) {
  raster_layer <- patch_rasters[[i]]
  patch_polygons <- as.polygons(raster_layer, values = TRUE, dissolve = TRUE, round = FALSE)
  metrics <- patch_metrics %>% filter(layer == i)
  metrics <- tidyr::pivot_wider(metrics, 
                                id_cols = "id", 
                                names_from = "metric", 
                                values_from = "value")
  metrics$Year <- (1989 + i)
  
  #extract average magnitude
  mean_mags <- extract(masked_rasters[[i]], patch_polygons, fun = "mean", ID = TRUE)
  names(mean_mags) <- c("ID", "mean_mag")
  #combine metrics
  patch_df <- as.data.frame(patch_polygons, geom = "WKT")
  names(patch_df) <- c("id", "geometry")
  metrics <- left_join(metrics, patch_df, by = c("id" = "id"))
  metrics <- left_join(metrics, mean_mags, by = c("id" = "ID"))
  metrics <- metrics %>% filter(area >= 0.27) #patch size threshold
  vect(metrics, geom = "geometry")
})

patches <- vect(patches)
patches$unique_id <- 1:nrow(patches)

writeVector(patches, paste0(tuning, "afp_patch_vector_1990-2022.gpkg"), overwrite = TRUE)

################## calculate total annual area disturbed (area weighted) ###############

#total land area owned by state in each year
yearly_areas <- lapply(1:nlyr(raster_stack), \(i) {
  year_polys <- afp_df %>% dplyr::filter(Year <= (1989 +i))
  year_area <- sum(year_polys$Area)/10000
  
  tibble("Year" = (1989 + i) ,"total_area_ha" = year_area)
})
yearly_area <- bind_rows(yearly_areas)

patches_df <- as.data.frame(patches)
annual_dist_area <- patches_df %>% group_by(Year) %>% summarise(area_dist_ha = sum(area))
annual_dist_area$total_area_ha <-  yearly_area$total_area_ha
annual_dist_area$percent_dist <-  (annual_dist_area$area_dist_ha / annual_dist_area$total_area_ha)*100

write.csv(annual_dist_area, paste0(tuning, "_annual_disturbed_area.csv"))

#################### summarise patch metrics by year #################
yearly_summaries <- patches_df %>% group_by(Year) %>% summarise(area = mean(area), 
                                                                cai = mean(cai), 
                                                                area = mean(area),    
                                                                cai = mean(cai), 
                                                                circle = mean(cai),
                                                                contig = mean(cai),
                                                                core = mean(core),
                                                                enn = mean(enn),
                                                                frac = mean(frac),
                                                                gyrate = mean(gyrate),
                                                                ncore = mean(ncore),
                                                                para = mean(para),
                                                                perim = mean(perim),
                                                                shape = mean(shape),
                                                                mean_mag = mean(mean_mag))

yearly_summaries$percent_area_dist <- annual_dist_area$percent_dist

write.csv(yearly_summaries, paste0("_yearly_metric_summaries.csv"), row.names = FALSE)
#yearly_summaries <- read.csv("yearly_metric_summaries.csv")


########## statistical tests ##############
library(trend)

### #man kendall tests ####
mk_test_results <- lapply(2:15, \(i){
  metric <- names(yearly_summaries)[i]
  data <- pull(yearly_summaries, i)
  mk_results <- mk.test(data)
  mk_estimates <- as.data.frame(mk_results[6])
  tibble("metric" = metric,
         "man_kendall_statistic" = mk_estimates[3,1],
         "normalized_test_stat" = as.numeric(mk_results[3]),
         "p_value" = as.numeric(mk_results[2]))
})
mk_test_results <- bind_rows(mk_test_results)

mk_test_results$trend <-as.character(lapply(1:nrow(mk_test_results), \(i){
  x$trend[i] <- if(mk_test_results$man_kendall_statistic[i] > 0 & mk_test_results$p_value[i] < 0.05){
    return("increasing")
  } else if (mk_test_results$man_kendall_statistic[i] < 0 & mk_test_results$p_value[i] < 0.05){
    return('decreasing')
  } else {
    return("no trend")
  }
}))


#### sens slope ####
sens_slope_results <- lapply(2:15, \(i){
  metric <- names(yearly_summaries)[i]
  data <- pull(yearly_summaries, i)
  sens_results <- sens.slope(data)
  tibble("metric" = metric,
         "sens slope" = as.numeric(sens_results[1]))
})
sens_slope_results <- bind_rows(sens_slope_results)

mk_test_results <- left_join(mk_test_results, sens_slope_results, by = "metric")

write.csv(mk_test_results, paste0(tuning, "_man_kendall_test_results.csv"), row.names = FALSE)
