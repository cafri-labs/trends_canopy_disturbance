library(dplyr)
library(terra)
library(landscapemetrics)

#code for assessment of disturbance trends by age group
#method is the same as general trends assessment, just split by age group

######### load data ######
#cafri crs is  "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
crs <- labrador.client::get_cafri_crs()
#requires annual loss brick, code for creation can be found here: https://github.com/cafri-labs/lt_tuning/blob/main/scripts/create_loss_layer.R
raster_stack <- rast("rt1_1985_2022_all_losses.tif")

#name of tuning
tuning <- "rt_1"

#adirondack forest preserve lands
afp_polygons <- vect("afp_aquisition_dates_final.gpkg")
afp_polygons <- project(afp_polygons, crs)

raster_stack <- crop(raster_stack, afp_polygons)

#mask all water out using NYS hydrology boundaries. file can be found here: https://data.gis.ny.gov/maps/sharegisny::nys-hydrography/about
hydro <- vect("../GIS_Files/hydrography.gdb")
hydro <- project(hydro, crs)

raster_stack <- mask(raster_stack, hydro, inverse = TRUE)

afp_df <- as.data.frame(afp_polygons, geom = "WKT")

#select parcels with only precisely determined acquisition dates
selected_parcels <- afp_df %>% filter(APADATE > 0 | !is.na(DATE_VEST)) %>% filter(Year < 1990)

########### create age bins ##############
yearly_area <- selected_parcels %>% group_by(Year) %>% summarise(sum_area_ha = (sum(Area)/10000))
#bin 1
lapply(1:nrow(yearly_area), \(i){
  area <- sum(yearly_area$sum_area_ha[1:i])
  if (area > 174454){
    return(yearly_area$Year[i])
  }
})
#bin 2
x <- yearly_area %>% filter(Year > 1893)

lapply(1:nrow(x), \(i){
  area <- sum(x$sum_area_ha[1:i])
  if (area > 174454){
    return(x$Year[i])
  }
})
#bin 3
x <- yearly_area %>% filter(Year > 1919)

lapply(1:nrow(x), \(i){
  area <- sum(x$sum_area_ha[1:i])
  if (area > 174454){
    return(x$Year[i])
  }
})

########### label bins ###########
#1: <= 1893 2:1894 - 1918 3: 1919 - 1956 4: 1957 - 1989
x <- lapply(1:nrow(selected_parcels), \(i){
  if(selected_parcels$Year[i] <= 1893){
    bin <-  1
  } else if((selected_parcels$Year[i] >= 1894)&(selected_parcels$Year[i] <= 1918)){
    bin <-  2
  } else if((selected_parcels$Year[i] >= 1919)&(selected_parcels$Year[i] <= 1956)){
    bin <-  3
  } else {
    bin <-  4
  }
  tibble("bin"= bin)
})
 x <- bind_rows(x)
selected_parcels <- bind_cols(selected_parcels, x)

sp_vect <- vect(selected_parcels, geom = "geometry")
crs(sp_vect) <- crs

#writeVector(sp_vect, "selected_binned_afp_parcels.gpkg")


########### create masked rasters ###################
masked_rasters <- mask(raster_stack, sp_vect)
masked_rasters[masked_rasters < 50] <- NA #distubance magnitude threshold
writeRaster(masked_rasters, paste0(tuning, "_age_group_masked_rasters.tif"), overwrite = TRUE)
#masked_rasters <- rast("age_group_masked_rasters.tif")

############### create class raster##############
class_rasters <- lapply(1:nlyr(masked_rasters), \(i) {
  raster <- masked_rasters[[i]]
  raster[!is.na(raster)] <- 1989 + i
  return(raster)
})

class_rasters <- rast(class_rasters)
writeRaster(class_rasters, paste0(tuning, "_age_groups_class_raster.tif"), overwrite = TRUE)
#class_rasters <- rast("age_groups_class_raster_rt1.tif")

##############create disturbance patch raster###############
patch_rasters <- get_patches(class_rasters, class = "all", 8, to_disk = TRUE)

patch_rasters <- lapply(1:length(patch_rasters), \(i) {
  patch_rasters[[i]][[1]]
})

patch_rasters <- rast(patch_rasters)
names(patch_rasters) <- c(1990:2022)

##########calculate patch metrics###########
options_landscapemetrics(to_disk = TRUE)
patch_metrics <- calculate_lsm(class_rasters, level = "patch", directions = 8, neighbourhood = 8 ,progress = TRUE)

#########combine patches and patch metrics##########
patches <- lapply(1:nlyr(patch_rasters), \(i) {
  raster_layer <- patch_rasters[[i]]
  patch_polygons <- as.polygons(raster_layer, values = TRUE, dissolve = TRUE)
  #calculate metrics with landscapemetrics package
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

writeVector(patches, paste0(tuning, "_afp_precise_dates_only_patch_vector_1990-2022_minor_age_group.gpkg"), overwrite = TRUE)


##############bin patches###############
#patches <- vect("afp_precise_dates_only_patch_vector_1990-2022_rt075.gpkg")
patches_df <- as.data.frame(patches)
selected_parcels <- as.data.frame(selected_parcels)

#group_1 
group1_parcels <- selected_parcels %>% filter(Year <= 1896) 
parcels <- vect(group1_parcels, geom = "geometry")
group1 <- mask(patches, parcels)

#group 2
group2_parcels <- selected_parcels %>% filter((Year >= 1897) & (Year <= 1909))
parcels <- vect(group2_parcels, geom = "geometry")
group2 <- mask(patches, parcels)
#group 3
group3_parcels <- selected_parcels %>% filter((Year >= 1910) & (Year <= 1957))
parcels <- vect(group3_parcels, geom = "geometry")
group3 <- mask(patches, parcels)
#group 4
group4_parcels <- selected_parcels %>% filter(Year >= 1958) 
parcels <- vect(group4_parcels, geom = "geometry")
group4 <- mask(patches, parcels)

############ summarise patch metrics by year AND by age class ##############
year_sum <- function(group_vector, group_parcels){
   group_df <- as.data.frame(group_vector)
  
  #summarise patch metrics by year AND by age class
  yearly_summaries <- group_df %>% group_by(Year) %>% summarise(area = mean(area), 
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
  annual_area_disturbed <- group_df %>% group_by(Year) %>% summarise(area_dist_ha = sum(area))
  annual_area_disturbed$total_area_ha <-  sum(group_parcels$Area)/10000
  annual_area_disturbed$percent_dist <-  (annual_area_disturbed$area_dist_ha / annual_area_disturbed$total_area_ha)*100
  
  yearly_summaries$percent_area_dist <- annual_area_disturbed$percent_dist
  return(yearly_summaries)
}

group1_sum <- year_sum(group1, group1_parcels)
group2_sum <- year_sum(group2, group2_parcels)
group3_sum <- year_sum(group3, group3_parcels)
group4_sum <- year_sum(group4, group4_parcels)

write.csv(group1_sum, paste0(tuning, "_group1_metrics.csv"), row.names = FALSE)
write.csv(group2_sum, paste0(tuning, "_group2_metrics.csv"), row.names = FALSE)
write.csv(group3_sum, paste0(tuning, "_group3_metrics.csv"), row.names = FALSE)
write.csv(group4_sum, paste0(tuning, "_group4_metrics.csv"), row.names = FALSE)


############### statistical tests ############
library(trend)

#man kendall tests
mk_tester <- function(summary){
mk_test_results <- lapply(2:15, \(i){
  metric <- names(summary)[i]
  data <- pull(summary, i)
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

#sens slope
sens_slope_results <- lapply(2:15, \(i){
  metric <- names(summary)[i]
  data <- pull(summary, i)
  sens_results <- sens.slope(data)
  tibble("metric" = metric,
         "sens slope" = as.numeric(sens_results[1]))
})
sens_slope_results <- bind_rows(sens_slope_results)
mk_test_results <- left_join(mk_test_results, sens_slope_results, by = "metric")
return(mk_test_results)
}

group1_mk_results <- mk_tester(group1_sum)
group2_mk_results <- mk_tester(group2_sum)
group3_mk_results <- mk_tester(group3_sum)
group4_mk_results <- mk_tester(group4_sum)

write.csv(group1_mk_results, paste0(tuning, "_group1_mk_results.csv"), row.names = FALSE)
write.csv(group2_mk_results, paste0(tuning, "_group2_mk_results.csv"), row.names = FALSE)
write.csv(group3_mk_results, paste0(tuning, "_group3_mk_results.csv"), row.names = FALSE)
write.csv(group4_mk_results, paste0(tuning, "_group4_mk_results.csv"), row.names = FALSE)