
#### setting working directory ####
setwd("~/capital_natural/soundscape")
memory.limit(1000000)

##### loading required packages ####
library(tidyverse)
library(dplyr)
library(lubridate)
library(ggplot2)
library(RColorBrewer)
library(Hmisc)
library(corrplot)
library(rstatix)
library(fitdistrplus)
library(scales)
library(glmmTMB)
library(AICcmodavg)
library(performance)
library(MuMIn)
library(DHARMa)
library(effects)
library(broom.mixed)
library(dotwhisker)
library(cowplot)
#library(lme4)
#library(multcomp)
#library(multcompView)
#library(emmeans)
#library(MASS)
#library(car)
#library(cowplot)


# soundscape index
CN.data.raw.ed2 <- read.csv("data/AcousticIndex_102022_v3.csv", header = T)
# reordering columns
CN.data.raw.ed2 <- CN.data.raw.ed2[,c(1:7,43,8:42,44)]

#### exploratory ####
# copy
CN.data <- CN.data.raw.ed2

## checking
#head(CN.data);tail(CN.data)
#str(CN.data)
#summary(CN.data)


#### adding environmental variables ####
# location and environmental data
CN.location <- read.csv("data/pontos_soundscape.csv", header = T)


# temperature, precipitation, air humidity, wind velocity
# source: https://bdmep.inmet.gov.br/
weather <- read.csv("data/raw_environmental_data/dados_A230_H_2019-11-01_2022-05-20 - Copia.csv", sep = ";") 
weather <- weather %>% 
               mutate(Year = as.numeric(substr(Data.Medicao, 1, 4)),
                      Month = as.numeric(sub("^0+", "", substr(Data.Medicao, 6, 7))),
                      Day = as.numeric(sub("^0+", "", substr(Data.Medicao, 9, 10))),
                      Hour = as.numeric(rep(0:23, length.out=nrow(weather)))) %>% 
               dplyr::select(8:11,3:6)

# NDVI (Normalized Difference Vegetation Index)
# source: https://neo.gsfc.nasa.gov/view.php?datasetId=MOD_NDVI_M
ndvi2019 <- raster::raster("data/raw_environmental_data/MOD_NDVI_M_2019-11-01_rgb_3600x1800.TIFF")
ndvi2022 <- raster::raster("data/raw_environmental_data/MOD_NDVI_M_2022-04-01_rgb_3600x1800.TIFF")

# Global Forest Canopy Height - 2019
# source: https://glad.umd.edu/dataset/gedi/
treeheight <- raster::raster("data/raw_environmental_data/Forest_height_2019_SAM.tif")


# transforming:
# adding Date column; 
# sorting sites from CN1 to CN14;
# adding time period category variable;
# and converting columns to factors
CN.data <- CN.data %>% 
             left_join(CN.location) %>% 
             left_join(weather) %>%
             mutate(Date = parse_date_time(sub(".*_", "", CN.data$ldt.index), orders = "%Y%m%d%H%M"),
                    Local=fct_relevel(Local,paste0("CN", seq(1:14))),
                    time_period = ifelse(as.numeric(as.character(Min)) >= 39 & as.numeric(as.character(Min)) <= 50, "dawn",
                                         ifelse(as.numeric(as.character(Min)) >= 51 & as.numeric(as.character(Min)) <= 110, "day",
                                                ifelse(as.numeric(as.character(Min)) >= 111 & as.numeric(as.character(Min)) <= 122, "dusk", "night")))) %>% 
             mutate_at(vars(ID:Min, time_period), factor)

CN.data$NDVI <- ifelse(CN.data$Year==2019, raster::extract(ndvi2019, sp::SpatialPoints(CN.data[,c("long","lat")])),
                       raster::extract(ndvi2022, sp::SpatialPoints(CN.data[,c("long","lat")])))

CN.data$treeheight <- raster::extract(treeheight, sp::SpatialPoints(CN.data[,c("long","lat")]))

str(CN.data)



# reordering columns
CN.data <- CN.data[,c(1,8,2,45,46,59,3:7,60,9:43,47:51,61,62,55:58,52:54,44)]

# saving
write.csv(CN.data, "data/AcousticIndex_102022_v4.csv", row.names = F)


#### between minutes by points ####
dir.create("results")


# Audio sampling by point after removal of recordings with noise only
sampling.overview <- CN.data %>% 
                       group_by(Local, Year) %>% 
                       summarise(Start_Date = min(Date),
                                 End_Date = max(Date),
                                 Count = n()) %>% 
                       mutate(N_days = as.numeric(End_Date - Start_Date)) %>% 
                       left_join(CN.location) %>% 
                       ungroup()

# reordering columns
sampling.overview <- sampling.overview[,c(1,7:9,2:4,6,5)]
# obs:
# Locals CN11, CN12, CN13 and CN14 need check the number of days for 2022 sampling
# these points started record in april, stoped and then re-started in may
sampling.overview[sampling.overview$Local=="CN11" & sampling.overview$Year==2022,"N_days"] <- 7
sampling.overview[sampling.overview$Local=="CN12" & sampling.overview$Year==2022,"N_days"] <- 8
sampling.overview[sampling.overview$Local=="CN13" & sampling.overview$Year==2022,"N_days"] <- 8
sampling.overview[sampling.overview$Local=="CN14" & sampling.overview$Year==2022,"N_days"] <- 8


# saving
write.csv(sampling.overview, "results/sampling_by_point_by_year.csv", row.names = F)

rm(list= ls()[(ls() %in% c("CN.location","ndvi2019","ndvi2022","treeheight","weather"))])
gc()

#### study area map ####
library(raster)
library(rgdal)
library(rgeos)
library(scales)
library(sp)
library(sf)
library(cowplot)
library(ggrepel)
library(ggspatial)

std.proj <- "+proj=longlat +datum=WGS84 +units=m +no_defs"

sa <- readOGR(dsn = "C:/Users/miral/Dropbox/GIS/Politico/America do Sul/Continente", layer = "tudo")
sa <- crop(sa, extent(x=c(-83,-17), y=c(-56,15)))
proj4string(sa) <- CRS(std.proj)
sa <- spTransform(sa, crs(std.proj))
plot(sa, col="gray95")

#amz <- readOGR(dsn = "C:/Users/miral/Dropbox/GIS/Politico/Amazonia bioma", layer = "Lim_Biogeografico")
#plot(amz, add=T, col="gray85")

br <- readOGR(dsn = "C:/Users/miral/Dropbox/GIS/Politico/Divisao Politica BR", layer = "BR_Contorno")
proj4string(br) <- CRS(std.proj)
br <- spTransform(br, crs(std.proj))
plot(br, add=T, col="gray85")

pa <- readOGR(dsn = "C:/Users/miral/Dropbox/GIS/Politico/Divisao Politica BR", layer = "BRASIL")
pa <- pa[pa@data$UF=="PA",]
proj4string(pa) <- CRS(std.proj)
pa <- spTransform(pa, crs(std.proj))
plot(pa, add=T, col="gray75")

#ogrListLayers("C:/Users/miral/Dropbox/GIS/Politico/flona_carajas.kml")
carajas <- readOGR(dsn = "C:/Users/miral/Dropbox/GIS/Politico/flona_carajas.kml", layer = "sql_statement")
proj4string(carajas) <- CRS(std.proj)
carajas <- spTransform(carajas, crs(std.proj))
#plot(carajas, add=T, col="gray55")
#extent(carajas)

#ogrListLayers("C:/Users/miral/Dropbox/GIS/Politico/parna_campos_ferruginosos.kml")
ferrugem <- readOGR(dsn = "C:/Users/miral/Dropbox/GIS/Politico/parna_campos_ferruginosos.kml", layer = "PARNA_Campos_Ferruginosos")
proj4string(ferrugem) <- CRS(std.proj)
ferrugem <- spTransform(ferrugem, crs(std.proj))
#plot(ferrugem, add=T, col="gray55")
#extent(ferrugem)

uc <- rbind(carajas, ferrugem)
uc@data[,"Name"] <- c("FLONA Carajas", "PARNA Campos Ferruginosos")
plot(uc, add=T, col="gray55")
#extent(uc)

lulc.carajas <- raster("C:/Users/miral/Dropbox/GIS/Politico/carajas.tif")
lulc.carajas <- crop(lulc.carajas, extent(x=c(-51.2,-49.8), y=c(-6.8,-5.4)))

lulc.carajas[lulc.carajas[]==3] <- 1
lulc.carajas[lulc.carajas[]==4] <- 1
lulc.carajas[lulc.carajas[]==11] <- 10
lulc.carajas[lulc.carajas[]==12] <- 10
lulc.carajas[lulc.carajas[]==15] <- 14
lulc.carajas[lulc.carajas[]==39] <- 14
lulc.carajas[lulc.carajas[]==41] <- 14

lulc.carajas.df <- as.data.frame(lulc.carajas, xy = TRUE)
head(lulc.carajas.df)

breakpoints <- sort(unique(values(lulc.carajas)))
mapbiomas.legend <- c("#129912", "#BBFCAC", "#FFFFB2", "#af2a2a", "#8a2be2", "#0000ff")
labels.legend <- c("Forest Formation", "Non Forest Natural Formation", "Farming", "Urban Area", "Mining", "Water")

location <- CN.data %>% dplyr::select(Local, long, lat) %>% distinct()

g1 <- ggplot() +
  geom_raster(data = lulc.carajas.df , aes(x = x, y = y, fill = factor(classification_2021)), alpha = 0.5) + 
  scale_fill_manual(breaks = breakpoints, values = mapbiomas.legend, labels = labels.legend, name = "LULC Classes") + 
  scale_y_continuous(position = "right", sec.axis = sec_axis(~., labels = NULL)) +
  geom_sf(data = st_as_sf(uc), aes(color = Name), fill = NA, linewidth = 1.2) +
  coord_sf(crs = crs(lulc.carajas), xlim = c(-51.8, -49.7), ylim=c(-6.9, -5.399773), expand = F, label_graticule = "ES") +
  scale_color_manual(values = c("#006400", "#B8AF4F"), name = "UCs") +
  geom_point(data = location, aes(x=long, y=lat), size=3.5) +
  geom_text_repel(data = location, aes(x=long, y=lat, label=Local),
                  size = 5,
                  min.segment.length = .1,
                  force = 3,
                  max.overlaps = 99) +
  annotation_scale(location = "br", width_hint = 0.25, text_family = "serif", text_cex = 1.5, line_width = 2, style = "ticks") +
  annotation_north_arrow(location = "tr", which_north = "true", pad_x = unit(2.5, "cm"), height = unit(2.5, "cm"), width = unit(2.5, "cm"),
                         style = north_arrow_fancy_orienteering(text_family = "serif", text_size = 16)) +
  labs(x = "", y = "") +
  theme_minimal(base_family = "serif") + 
  theme(axis.text = element_text(size = 12),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18),
        legend.position = c(0, .98),
        legend.justification = c(0, 1),
        panel.background = element_blank())


g2 <- ggplot() +
  geom_sf(data = st_as_sf(sa), fill = "gray95", linewidth = 1) +
  geom_sf(data = st_as_sf(br), fill = "gray85", linewidth = .8) +
  geom_sf(data = st_as_sf(pa), fill = "gray75", linewidth = .5) +
  geom_sf(data = st_as_sf(uc), fill = "gray55", linewidth = .2) +
  theme_void()


png("results/study_area.png", width = 1440, height = 754, units = "px")

ggdraw(g1) + draw_plot(g2, x = .14, y = .4, width = 0.16, height = 0.26)

dev.off()




#### correlogram ####
png("results/correlograms.png", width = 2560, height = 2560, units = "px", bg = "transparent")
par(mfrow = c(5, 5))

# frequency bin 1 [0.3 - 4kHz], dawn [0630 - 0830]
corrplot(cor(CN.data[CN.data$time_period=="dawn", grep("fbin1", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[CN.data$time_period=="dawn", grep("fbin1", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# frequency bin 1 [0.3 - 4kHz], day [0830 - 1830]
corrplot(cor(CN.data[CN.data$time_period=="day", grep("fbin1", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[CN.data$time_period=="day", grep("fbin1", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# frequency bin 1 [0.3 - 4kHz], dusk [1830 - 2030]
corrplot(cor(CN.data[CN.data$time_period=="dusk", grep("fbin1", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[CN.data$time_period=="dusk", grep("fbin1", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# frequency bin 1 [0.3 - 4kHz], night [2030 - 0630]
corrplot(cor(CN.data[CN.data$time_period=="night", grep("fbin1", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[CN.data$time_period=="night", grep("fbin1", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# frequency bin 1 [0.3 - 4kHz], total [24h]
corrplot(cor(CN.data[, grep("fbin1", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[, grep("fbin1", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# frequency bin 2 [4 - 12kHz], dawn [0630 - 0830]
corrplot(cor(CN.data[CN.data$time_period=="dawn", grep("fbin2", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[CN.data$time_period=="dawn", grep("fbin2", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# frequency bin 2 [4 - 12kHz], day [0830 - 1830]
corrplot(cor(CN.data[CN.data$time_period=="day", grep("fbin2", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[CN.data$time_period=="day", grep("fbin2", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# frequency bin 2 [4 - 12kHz], dusk [1830 - 2030]
corrplot(cor(CN.data[CN.data$time_period=="dusk", grep("fbin2", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[CN.data$time_period=="dusk", grep("fbin2", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# frequency bin 2 [4 - 12kHz], night [2030 - 0630]
corrplot(cor(CN.data[CN.data$time_period=="night", grep("fbin2", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[CN.data$time_period=="night", grep("fbin2", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# frequency bin 2 [4 - 12kHz], total [24h]
corrplot(cor(CN.data[, grep("fbin2", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[, grep("fbin2", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# frequency bin 3 [0.3 - 12kHz], dawn [0630 - 0830]
corrplot(cor(CN.data[CN.data$time_period=="dawn", grep("fbin3", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[CN.data$time_period=="dawn", grep("fbin3", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# frequency bin 3 [0.3 - 12kHz], day [0830 - 1830]
corrplot(cor(CN.data[CN.data$time_period=="day", grep("fbin3", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         p.mat = rcorr(as.matrix(CN.data[CN.data$time_period=="day", grep("fbin3", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# frequency bin 3 [0.3 - 12kHz], dusk [1830 - 2030]
corrplot(cor(CN.data[CN.data$time_period=="dusk", grep("fbin3", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[CN.data$time_period=="dusk", grep("fbin3", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# frequency bin 3 [0.3 - 12kHz], night [2030 - 0630]
corrplot(cor(CN.data[CN.data$time_period=="night", grep("fbin3", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[CN.data$time_period=="night", grep("fbin3", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# frequency bin 3 [0.3 - 12kHz], total [24h]
corrplot(cor(CN.data[, grep("fbin3", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[, grep("fbin3", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# frequency bin 4 [12 - 22.1kHz], dawn [0630 - 0830]
corrplot(cor(CN.data[CN.data$time_period=="dawn", grep("fbin4", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[CN.data$time_period=="dawn", grep("fbin4", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# frequency bin 4 [12 - 22.1kHz], day [0830 - 1830]
corrplot(cor(CN.data[CN.data$time_period=="day", grep("fbin4", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[CN.data$time_period=="day", grep("fbin4", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# frequency bin 4 [12 - 22.1kHz], dusk [1830 - 2030]
corrplot(cor(CN.data[CN.data$time_period=="dusk", grep("fbin4", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[CN.data$time_period=="dusk", grep("fbin4", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# frequency bin 4 [12 - 22.1kHz], night [2030 - 0630]
corrplot(cor(CN.data[CN.data$time_period=="night", grep("fbin4", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[CN.data$time_period=="night", grep("fbin4", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# frequency bin 4 [12 - 22.1kHz], total [24h]
corrplot(cor(CN.data[, grep("fbin4", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[, grep("fbin4", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# total frequency [0.3 - 22.1kHz], dawn [0630 - 0830]
corrplot(cor(CN.data[CN.data$time_period=="dawn", grep("Total", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[CN.data$time_period=="dawn", grep("Total", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# total frequency [0.3 - 22.1kHz], day [0830 - 1830]
corrplot(cor(CN.data[CN.data$time_period=="day", grep("Total", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[CN.data$time_period=="day", grep("Total", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# total frequency [0.3 - 22.1kHz], dusk [1830 - 2030]
corrplot(cor(CN.data[CN.data$time_period=="dusk", grep("Total", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[CN.data$time_period=="dusk", grep("Total", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# total frequency [0.3 - 22.1kHz], night [2030 - 0630]
corrplot(cor(CN.data[CN.data$time_period=="night", grep("Total", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[CN.data$time_period=="night", grep("Total", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

dev.off()
###
# total frequency [0.3 - 22.1kHz], total [24h]
png("results/correlogram_total.png", width = 960, height = 960, units = "px", bg = "transparent")

par(mfrow = c(1, 1))
corrplot(cor(CN.data[, grep("Total", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.col = "black", tl.srt = 90, tl.pos = "ld", tl.cex = 2, cl.pos = "r", cl.cex = 2,
         # combine with significance level
         p.mat = rcorr(as.matrix(CN.data[, grep("Total", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

dev.off()
#####

# excluding variables: S2N, NPeak, ADI
# they were built from other index (e.g., S2N from BGN)
# or they were used to build other index (e.g., NPeak to AA)
# and/or were highly correlated (H and AA and ADI)
#CN.data.total <- CN.data
#
#exclude <- c("S2NTotal", "NPeakTotal", "ADITotal",
#             "S2Nfbin1", "NPeakfbin1", "ADIfbin1",
#             "S2Nfbin2", "NPeakfbin2", "ADIfbin2",
#             "S2Nfbin3", "NPeakfbin3", "ADIfbin3",
#             "S2Nfbin4", "NPeakfbin4", "ADIfbin4")
#
#CN.data <- CN.data[,!colnames(CN.data) %in% exclude]

#### difference between time period ####
# creating data frame

time.period.diff <- data.frame(index = rep(c("H", "BGN", "S2N", "AA", "ADI", "AEI"), 70),
                                frequency.bin = rep(rep(c("Total", "0.3-4", "4-12", "0.3-12", "12-22.1"), each = 6), 14),
                                local = rep(unique(CN.data$Local), each = 30),
                                effect.size = NA, n.sig = NA)

c=1
for (l in unique(CN.data$Local)) {
  
  cnx <- CN.data[CN.data$Local==l,]
  
  for (i in grep("NPeak", names(CN.data)[13:47], invert = T, value = T)) {
    
    cnx.index <- cnx[sample(1:nrow(cnx), 500, replace = T), c(i,"time_period")]
    
    if (all(is.na(cnx.index[,i])) | all(cnx.index[,i]==0)) {
      
      time.period.diff$effect.size[c] <- NA
      
      time.period.diff$n.sig[c] <- NA
    } 
    
    else {
      
      time.period.diff$effect.size[c] <- as.numeric(kruskal_effsize(cnx.index, formula(paste(i, "~ time_period"))) %>% 
                                                      dplyr::select(effsize))
      
      time.period.diff$n.sig[c] <- as.numeric(dunn_test(cnx.index, formula(paste(i, "~ time_period")), p.adjust.method = "bonferroni") %>% 
                                                summarise(n.sig = sum(p.adj.signif != "ns")))
      
    }
    
    
    c=c+1
  }
  
}


#write.csv(time.period.diff, "results/timeperioddiff.csv", row.names = F)

#checking
str(time.period.diff)

time.period.diff <- time.period.diff %>% 
                       mutate(index=factor(index, levels = c("H", "ADI", "AEI", "BGN", "S2N", "AA")),
                              frequency.bin=factor(frequency.bin, levels = c("Total", "0.3-4", "4-12", "0.3-12", "12-22.1")),
                              local=factor(local, levels = paste0("CN", seq(1:14))))
  
  
time.period.diff.hline <- time.period.diff[time.period.diff$frequency.bin=="Total", ]
  
png("results/timeperioddifflocal.png", width = 1440, height = 2560, units = "px", bg = "transparent")

time.period.diff %>%  
  ggplot(aes(frequency.bin, n.sig, colour = effect.size))+
  geom_point(size=12)+
  geom_hline(data=time.period.diff.hline, aes(yintercept=n.sig), linewidth=1, linetype='dashed') +
  scale_y_continuous(limits = c(0, 7), breaks = c(0, 2, 4, 6), expand = c(0.1,0.1)) +
  scale_color_gradient("Effect size", low = "#F4A582", high = "#831529", na.value = NA)+
  guides(size = "none",
         colour = guide_colourbar(barheight = unit(45, "cm"))) +
  ylab(expression(atop("Number of significantly different", "time period pairs")))+xlab("Frequency bins (kHz)")+
  facet_grid(local~index) +
  theme(axis.title = element_text(family = "serif", size = 44),
        axis.text.x = element_text(family = "serif", size = 36, angle = 90, hjust = 1),
        axis.text.y = element_text(family = "serif", size = 36),
        axis.line.y = element_line(linewidth = 1), axis.line.x = element_line(linewidth = 1),
        axis.ticks.y = element_line(linewidth = 1), axis.ticks.x = element_line(linewidth = 1),
        legend.title = element_text(family = "serif", size = 36),
        legend.text = element_text(family = "serif", size = 36),
        panel.background = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        panel.spacing.y = unit(2, "lines"),
        strip.text = element_text(family = "serif", size = 36))


dev.off()
#

time.period.diff.total <- data.frame(index = rep(c("H", "BGN", "S2N", "AA", "ADI", "AEI"), 5),
                                     frequency.bin = rep(c("Total", "0.3-4", "4-12", "0.3-12", "12-22.1"), each = 6),
                                     effect.size = NA, n.sig = NA)

c=1
for (i in grep("NPeak", names(CN.data)[13:47], invert = T, value = T)) {
  
  cnx.index <- CN.data %>% group_by(Local) %>% sample_n(size = 500, replace = T) %>% ungroup() %>% dplyr::select(time_period, i)
  
  if (all(is.na(cnx.index[,i])) | all(cnx.index[,i]==0)) {
    
    time.period.diff.total$effect.size[c] <- NA
    
    time.period.diff.total$n.sig[c] <- NA
  } 
  
  else {
    
    time.period.diff.total$effect.size[c] <- as.numeric(kruskal_effsize(cnx.index, formula(paste(i, "~ time_period"))) %>% 
                                                    dplyr::select(effsize))
    
    time.period.diff.total$n.sig[c] <- as.numeric(dunn_test(cnx.index, formula(paste(i, "~ time_period")), p.adjust.method = "bonferroni") %>% 
                                              summarise(n.sig = sum(p.adj.signif != "ns")))
    
  }
  c=c+1
}

#checking
str(time.period.diff.total)

time.period.diff.total <- time.period.diff.total %>% 
                          mutate(index=factor(index, levels = c("H", "ADI", "AEI", "BGN", "S2N", "AA")),
                                 frequency.bin=factor(frequency.bin, levels = c("Total", "0.3-4", "4-12", "0.3-12", "12-22.1")))


time.period.diff.total.hline <- time.period.diff.total[time.period.diff.total$frequency.bin=="Total", ]

png("results/timeperioddifftotal.png", width = 2400, height = 640, units = "px", bg = "transparent")

time.period.diff.total %>% filter(!index %in% c("ADI", "S2N")) %>% 
  ggplot(aes(frequency.bin, n.sig, colour = effect.size))+
  geom_point(size=12)+
  geom_hline(data=time.period.diff.total.hline %>% filter(!index %in% c("ADI", "S2N"))
             , aes(yintercept=n.sig), linewidth=1, linetype='dashed') +
  scale_y_continuous(limits = c(0, 7), breaks = c(0, 2, 4, 6), expand = c(0.1,0.1)) +
  scale_color_gradient("Effect size", low = "#F4A582", high = "#831529", na.value = NA)+
  guides(size = "none",
         colour = guide_colourbar(barheight = unit(12, "cm"))) +
  ylab(expression(atop("Number of significantly different", "time periods pairs")))+xlab("Frequency bins (kHz)")+
  facet_wrap(~index, nrow = 1) +
  theme(axis.title = element_text(family = "serif", size = 34),
        axis.text.x = element_text(family = "serif", size = 30, angle = 90, hjust = 1),
        axis.text.y = element_text(family = "serif", size = 30),
        axis.line.y = element_line(linewidth = 1), axis.line.x = element_line(linewidth = 1),
        axis.ticks.y = element_line(linewidth = 1), axis.ticks.x = element_line(linewidth = 1),
        legend.title = element_text(family = "serif", size = 30),
        legend.text = element_text(family = "serif", size = 30),
        panel.background = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        panel.spacing.y = unit(2, "lines"),
        strip.text = element_text(family = "serif", size = 30))


dev.off()
#

CN.data.2019 <- CN.data[CN.data$Year=="2019",]

time.period.diff.total.2019 <- data.frame(index = rep(c("H", "BGN", "S2N", "AA", "ADI", "AEI"), 5),
                                     frequency.bin = rep(c("Total", "0.3-4", "4-12", "0.3-12", "12-22.1"), each = 6),
                                     effect.size = NA, n.sig = NA)

c=1
for (i in grep("NPeak", names(CN.data)[13:47], invert = T, value = T)) {
  
  cnx.index <- CN.data.2019 %>% group_by(Local) %>% sample_n(size = 500, replace = T) %>% ungroup() %>% dplyr::select(time_period, i)
  
  if (all(is.na(cnx.index[,i])) | all(cnx.index[,i]==0)) {
    
    time.period.diff.total.2019$effect.size[c] <- NA
    
    time.period.diff.total.2019$n.sig[c] <- NA
  } 
  
  else {
    
    time.period.diff.total.2019$effect.size[c] <- as.numeric(kruskal_effsize(cnx.index, formula(paste(i, "~ time_period"))) %>% 
                                                          dplyr::select(effsize))
    
    time.period.diff.total.2019$n.sig[c] <- as.numeric(dunn_test(cnx.index, formula(paste(i, "~ time_period")), p.adjust.method = "bonferroni") %>% 
                                                    summarise(n.sig = sum(p.adj.signif != "ns")))
    
  }
  c=c+1
}

#checking
str(time.period.diff.total.2019)

time.period.diff.total.2019 <- time.period.diff.total.2019 %>% 
  mutate(index=factor(index, levels = c("H", "ADI", "AEI", "BGN", "S2N", "AA")),
         frequency.bin=factor(frequency.bin, levels = c("Total", "0.3-4", "4-12", "0.3-12", "12-22.1")))


time.period.diff.total.2019.hline <- time.period.diff.total.2019[time.period.diff.total.2019$frequency.bin=="Total", ]

png("results/timeperioddifftotal2019.png", width = 2400, height = 640, units = "px", bg = "transparent")

time.period.diff.total.2019 %>%  
  ggplot(aes(frequency.bin, n.sig, colour = effect.size))+
  geom_point(size=12)+
  geom_hline(data=time.period.diff.total.2019.hline, aes(yintercept=n.sig), linewidth=1, linetype='dashed') +
  scale_y_continuous(limits = c(0, 7), breaks = c(0, 2, 4, 6), expand = c(0.1,0.1)) +
  scale_color_gradient("Effect size", low = "#F4A582", high = "#831529", na.value = NA)+
  guides(size = "none",
         colour = guide_colourbar(barheight = unit(12, "cm"))) +
  ylab(expression(atop("Number of significantly different", "time periods pairs")))+xlab("Frequency bins (kHz)")+
  facet_wrap(~index, nrow = 1) +
  theme(axis.title = element_text(family = "serif", size = 34),
        axis.text.x = element_text(family = "serif", size = 30, angle = 90, hjust = 1),
        axis.text.y = element_text(family = "serif", size = 30),
        axis.line.y = element_line(linewidth = 1), axis.line.x = element_line(linewidth = 1),
        axis.ticks.y = element_line(linewidth = 1), axis.ticks.x = element_line(linewidth = 1),
        legend.title = element_text(family = "serif", size = 30),
        legend.text = element_text(family = "serif", size = 30),
        panel.background = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        panel.spacing.y = unit(2, "lines"),
        strip.text = element_text(family = "serif", size = 30))


dev.off()
#

CN.data.2022 <- CN.data[CN.data$Year=="2022",]

time.period.diff.total.2022 <- data.frame(index = rep(c("H", "BGN", "S2N", "AA", "ADI", "AEI"), 5),
                                          frequency.bin = rep(c("Total", "0.3-4", "4-12", "0.3-12", "12-22.1"), each = 6),
                                          effect.size = NA, n.sig = NA)

c=1
for (i in grep("NPeak", names(CN.data)[13:47], invert = T, value = T)) {
  
  cnx.index <- CN.data.2022 %>% group_by(Local) %>% sample_n(size = 39, replace = T) %>% ungroup() %>% dplyr::select(time_period, i)
  
  if (all(is.na(cnx.index[,i])) | all(cnx.index[,i]==0)) {
    
    time.period.diff.total.2022$effect.size[c] <- NA
    
    time.period.diff.total.2022$n.sig[c] <- NA
  } 
  
  else {
    
    time.period.diff.total.2022$effect.size[c] <- as.numeric(kruskal_effsize(cnx.index, formula(paste(i, "~ time_period"))) %>% 
                                                               dplyr::select(effsize))
    
    time.period.diff.total.2022$n.sig[c] <- as.numeric(dunn_test(cnx.index, formula(paste(i, "~ time_period")), p.adjust.method = "bonferroni") %>% 
                                                         summarise(n.sig = sum(p.adj.signif != "ns")))
    
  }
  c=c+1
}

#checking
str(time.period.diff.total.2022)

time.period.diff.total.2022 <- time.period.diff.total.2022 %>% 
  mutate(index=factor(index, levels = c("H", "ADI", "AEI", "BGN", "S2N", "AA")),
         frequency.bin=factor(frequency.bin, levels = c("Total", "0.3-4", "4-12", "0.3-12", "12-22.1")))


time.period.diff.total.2022.hline <- time.period.diff.total.2022[time.period.diff.total.2022$frequency.bin=="Total", ]

png("results/timeperioddifftotal2022.png", width = 2400, height = 640, units = "px", bg = "transparent")

time.period.diff.total.2022 %>%  
  ggplot(aes(frequency.bin, n.sig, colour = effect.size))+
  geom_point(size=12)+
  geom_hline(data=time.period.diff.total.2022.hline, aes(yintercept=n.sig), linewidth=1, linetype='dashed') +
  scale_y_continuous(limits = c(0, 7), breaks = c(0, 2, 4, 6), expand = c(0.1,0.1)) +
  scale_color_gradient("Effect size", low = "#F4A582", high = "#831529", na.value = NA)+
  guides(size = "none",
         colour = guide_colourbar(barheight = unit(12, "cm"))) +
  ylab(expression(atop("Number of significantly different", "time periods pairs")))+xlab("Frequency bins (kHz)")+
  facet_wrap(~index, nrow = 1) +
  theme(axis.title = element_text(family = "serif", size = 34),
        axis.text.x = element_text(family = "serif", size = 30, angle = 90, hjust = 1),
        axis.text.y = element_text(family = "serif", size = 30),
        axis.line.y = element_line(linewidth = 1), axis.line.x = element_line(linewidth = 1),
        axis.ticks.y = element_line(linewidth = 1), axis.ticks.x = element_line(linewidth = 1),
        legend.title = element_text(family = "serif", size = 30),
        legend.text = element_text(family = "serif", size = 30),
        panel.background = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        panel.spacing.y = unit(2, "lines"),
        strip.text = element_text(family = "serif", size = 30))


dev.off()
#

rm(list= ls()[(ls() %in% c("exclude","c","l","i","cnx","cnx.index"))])
gc()

#### Compute summary statistics by groups - count, min, max, mean, sd ####
#overview <- CN.data %>% 
#  group_by(Local, time_period) %>%
#  summarise(Count = n(),
#            H.Min = min(HTotal, na.rm = TRUE),
#            H.Max = max(HTotal, na.rm = TRUE),
#            H.Mean = mean(HTotal, na.rm = TRUE),
#            H.sd = sd(HTotal, na.rm = TRUE),
#            BGN.Min = min(BGNTotal, na.rm = TRUE),
#            BGN.Max = max(BGNTotal, na.rm = TRUE),
#            BGN.Mean = mean(BGNTotal, na.rm = TRUE),
#            BGN.sd = sd(BGNTotal, na.rm = TRUE),
#            AA.Min = min(AATotal, na.rm = TRUE),
#            AA.Max = max(AATotal, na.rm = TRUE),
#            AA.Mean = mean(AATotal, na.rm = TRUE),
#            AA.sd = sd(AATotal, na.rm = TRUE),
#            AEI.Min = min(AEITotal, na.rm = TRUE),
#            AEI.Max = max(AEITotal, na.rm = TRUE),
#            AEI.Mean = mean(AEITotal, na.rm = TRUE),
#            AEI.sd = sd(AEITotal, na.rm = TRUE)) %>% 
#  ungroup()
#
##distribution of values between minutes
#cols <- c("CN1" = "#eff1ed", "CN2" = "#373d20", "CN3" = "#717744", "CN4" = "#bcbd8b", "CN5" = "#766153", "CN6" = "#9cc5a1",
#          "CN7" = "#49a078", "CN8" = "#216869", "CN9" = "#f4f0bb", "CN10" = "#43291f", "CN11" = "#d36135", "CN12" = "#bc4b51",
#          "CN13" = "#427aa1", "CN14" = "#004e64")
#
#png("results/overview_byMin.png", width = 2400, height = 640, units = "px", bg = "transparent")
#
#CN.data %>% dplyr::select("Local", "Min", "HTotal", "BGNTotal", "AATotal", "AEITotal") %>% 
#  pivot_longer(cols = HTotal:AEITotal, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
#  mutate(Index = gsub("Total", "", Index)) %>% 
#  mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
#     ggplot(aes(Min, Values, group=Local,  colour = Local)) +
#         geom_point(alpha = .25, size = .95, show.legend = F) +
#         geom_smooth(linewidth=2, method = "loess", se = F) +
#         geom_vline(aes(xintercept=39), linewidth=1, linetype='dashed', color = "gray55") +
#         geom_vline(aes(xintercept=51), linewidth=1, linetype='dashed', color = "gray55") +
#         geom_vline(aes(xintercept=111), linewidth=1, linetype='dashed', color = "gray55") +
#         geom_vline(aes(xintercept=123), linewidth=1, linetype='dashed', color = "gray55") +
#         scale_x_discrete(breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
#         scale_color_manual(values = cols)+
#         facet_wrap(Index~., ncol=4, scales = "free") +
#         guides(col = guide_legend(nrow = 2, byrow = TRUE))+
#         theme_minimal() +
#         theme(axis.title = element_text(family = "serif", size = 30),
#               axis.text = element_text(family = "serif", size = 28, hjust = 1),
#               axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
#               axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2),
#               strip.text = element_text(family = "serif", size = 28),
#               legend.title = element_blank(),
#               legend.text = element_text(family = "serif", size = 22),
#               legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center")
#  
#dev.off()



#distribution of values between minutes -- considering Total frequency [0.3 - 22.1kHz]

cols <- c("CN1" = "#eff1ed", "CN2" = "#373d20", "CN3" = "#717744", "CN4" = "#bcbd8b", "CN5" = "#766153", "CN6" = "#9cc5a1",
          "CN7" = "#49a078", "CN8" = "#216869", "CN9" = "#f4f0bb", "CN10" = "#43291f", "CN11" = "#d36135", "CN12" = "#bc4b51",
          "CN13" = "#692168", "CN14" = "#692221")

p1 <- CN.data %>% dplyr::select("Local", "Min", "HTotal") %>% 
  pivot_longer(cols = HTotal, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
  mutate(Index = gsub("Total", "", Index)) %>% 
  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
  geom_point(alpha = .25, size = 1.5, show.legend = F) +
  geom_smooth(linewidth=3, method = "loess", se = F) +
  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
  scale_x_discrete("", breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
  scale_color_manual(values = cols)+
  ggtitle("H [0.3-22.1kHz]")+xlab("")+ylab("")+
  #facet_wrap(Index~., ncol=2, scales = "free") +
  #scale_y_continuous(limits = c(0.6,1.0)) +
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme_minimal() +
  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
        axis.title = element_text(family = "serif", size = 48),
        axis.text.y = element_text(family = "serif", size = 32, hjust = 1), axis.text.x = element_blank(),
        axis.line.y = element_line(size = 2), axis.line.x = element_blank(),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_blank(),
        strip.text = element_text(family = "serif", size = 48))#,
#        legend.title = element_blank(),
#        legend.text = element_text(family = "serif", size = 38),
#        legend.position = "none")
#legend.direction = "horizontal", legend.justification = "center")

p2 <- CN.data %>% dplyr::select("Local", "Min", "ADITotal") %>% 
  pivot_longer(cols = ADITotal, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
  mutate(Index = gsub("Total", "", Index)) %>% 
  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
  geom_point(alpha = .25, size = 1.5, show.legend = F) +
  geom_smooth(linewidth=3, method = "loess", se = F) +
  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
  scale_x_discrete("", breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
  scale_color_manual(values = cols)+
  ggtitle("ADI")+xlab("")+ylab("")+
  #facet_wrap(Index~., ncol=2, scales = "free") +
  #scale_y_continuous("", limits = c(0,.8)) +
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme_minimal() +
  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
        axis.title = element_text(family = "serif", size = 48),
        axis.text.y = element_text(family = "serif", size = 32, hjust = 1), axis.text.x = element_blank(),
        axis.line.y = element_line(size = 2), axis.line.x = element_blank(),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_blank(),
        strip.text = element_text(family = "serif", size = 48))#,
#        legend.title = element_blank(),
#        legend.text = element_text(family = "serif", size = 38),
#        legend.position = "none")
#legend.direction = "horizontal", legend.justification = "center")

p3 <- CN.data %>% dplyr::select("Local", "Min", "AEITotal") %>% 
  pivot_longer(cols = AEITotal, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
  mutate(Index = gsub("Total", "", Index)) %>% 
  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
  geom_point(alpha = .25, size = 1.5, show.legend = F) +
  geom_smooth(linewidth=3, method = "loess", se = F) +
  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
  scale_x_discrete("", breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
  scale_color_manual(values = cols)+
  ggtitle("AEI")+xlab("")+ylab("")+
  #facet_wrap(Index~., ncol=2, scales = "free") +
  #scale_y_continuous("", limits = c(0,.8)) +
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme_minimal() +
  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
        axis.title = element_text(family = "serif", size = 48),
        axis.text.y = element_text(family = "serif", size = 32, hjust = 1), axis.text.x = element_blank(),
        axis.line.y = element_line(size = 2), axis.line.x = element_blank(),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_blank(),
        strip.text = element_text(family = "serif", size = 48))#,
#        legend.title = element_blank(),
#        legend.text = element_text(family = "serif", size = 38),
#        legend.position = "none")
#legend.direction = "horizontal", legend.justification = "center")

p4 <- CN.data %>% dplyr::select("Local", "Min", "S2NTotal") %>% 
  pivot_longer(cols = S2NTotal, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
  mutate(Index = gsub("Total", "", Index)) %>% 
  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
  geom_point(alpha = .25, size = 1.5, show.legend = F) +
  geom_smooth(linewidth=3, method = "loess", se = F) +
  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
  scale_x_discrete(breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
  scale_color_manual(values = cols)+
  ggtitle("S2N")+xlab("")+ylab("")+
  #facet_wrap(Index~., ncol=2, scales = "free") +
  #scale_y_continuous(limits = c(-25,-10)) +
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme_minimal() +
  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
        axis.title = element_text(family = "serif", size = 48),
        axis.text.y = element_text(family = "serif", size = 32, hjust = 1), axis.text.x = element_blank(),
        axis.line.y = element_line(size = 2), axis.line.x = element_blank(),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_blank(),
        strip.text = element_text(family = "serif", size = 48))#,
#        legend.title = element_blank(),
#        legend.text = element_text(family = "serif", size = 42),
#        legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center")

p5 <- CN.data %>% dplyr::select("Local", "Min", "BGNTotal") %>% 
  pivot_longer(cols = BGNTotal, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
  mutate(Index = gsub("Total", "", Index)) %>% 
  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
  geom_point(alpha = .25, size = 1.5, show.legend = F) +
  geom_smooth(linewidth=3, method = "loess", se = F) +
  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
  scale_x_discrete(breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
  scale_color_manual(values = cols)+
  ggtitle("BGN")+xlab("")+ylab("")+
  #facet_wrap(Index~., ncol=2, scales = "free") +
  #scale_y_continuous(limits = c(-25,-10)) +
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme_minimal() +
  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
        axis.title = element_text(family = "serif", size = 48),
        axis.text.y = element_text(family = "serif", size = 32, hjust = 1), axis.text.x = element_blank(),
        axis.line.y = element_line(size = 2), axis.line.x = element_blank(),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_blank(),
        strip.text = element_text(family = "serif", size = 48))#,
#        legend.title = element_blank(),
#        legend.text = element_text(family = "serif", size = 42),
#        legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center")

p6 <- CN.data %>% dplyr::select("Local", "Min", "AATotal") %>% 
  pivot_longer(cols = AATotal, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
  mutate(Index = gsub("Total", "", Index)) %>% 
  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
  geom_point(alpha = .25, size = 1.5, show.legend = F) +
  geom_smooth(linewidth=3, method = "loess", se = F) +
  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
  scale_x_discrete(breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
  scale_color_manual(values = cols)+
  ggtitle("AA") +xlab("")+ylab("")+
  #facet_wrap(Index~., ncol=2, scales = "free") +
  #scale_y_continuous("", limits = c(0,.5)) +
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme_minimal() +
  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
        axis.title = element_text(family = "serif", size = 48),
        axis.text = element_text(family = "serif", size = 32, hjust = 1),
        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2),
        strip.text = element_text(family = "serif", size = 48))#,
#        legend.title = element_blank(),
#        legend.text = element_text(family = "serif", size = 38),
#        legend.position = "none")



#prow <- plot_grid(p1 + theme(legend.position="none"), 
#                  p2 + theme(legend.position="none"), 
#                  p3 + theme(legend.position="none"), 
#                  p4 + theme(legend.position="none"), 
#                  p5 + theme(legend.position="none"), 
#                  p6 + theme(legend.position="none"), 
#                  p7 + theme(legend.position="none"), 
#                  nrow=2, axis = "b", align = "hv")
#
#legend <- get_legend(p1 + theme(legend.title = element_blank(), legend.text = element_text(family = "serif", size = 42),
#                                  legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center"))
#
#
#png("results/overview_byMin2.png", width = 33.11024, height = 23.38583, units = "in", res = 300, pointsize = 4)
#plot_grid( prow, legend, ncol = 1, rel_heights = c(1, .2))
#dev.off()





##distribution of values between minutes -- considering frequency bin 1 [0.3 - 4kHz]
#
#pb1a <- CN.data %>% dplyr::select("Local", "Min", "Hfbin1") %>% 
#  pivot_longer(cols = Hfbin1, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
#  mutate(Index = gsub("fbin1", "", Index)) %>% 
#  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
#  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
#  geom_point(alpha = .25, size = 1.5, show.legend = F) +
#  geom_smooth(linewidth=3, method = "loess", se = F) +
#  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
#  scale_x_discrete("", breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
#  scale_color_manual(values = cols)+
#  ggtitle("H")+
#  #facet_wrap(Index~., ncol=2, scales = "free") +
#  #scale_y_continuous(limits = c(0.45,0.75)) +
#  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
#  theme_minimal() +
#  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
#        axis.title = element_text(family = "serif", size = 48),
#        axis.text.y = element_text(family = "serif", size = 32, hjust = 1), axis.text.x = element_blank(),
#        axis.line.y = element_line(size = 2), axis.line.x = element_blank(),
#        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_blank(),
#        strip.text = element_text(family = "serif", size = 48))#,
##        legend.title = element_blank(),
##        legend.text = element_text(family = "serif", size = 38),
##        legend.position = "none")
##legend.direction = "horizontal", legend.justification = "center")
#
#pb1b <- CN.data %>% dplyr::select("Local", "Min", "AEIfbin1") %>% 
#  pivot_longer(cols = AEIfbin1, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
#  mutate(Index = gsub("fbin1", "", Index)) %>% 
#  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
#  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
#  geom_point(alpha = .25, size = 1.5, show.legend = F) +
#  geom_smooth(linewidth=3, method = "loess", se = F) +
#  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
#  scale_x_discrete("", breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
#  scale_color_manual(values = cols)+
#  ggtitle("AEI")+
#  #facet_wrap(Index~., ncol=2, scales = "free") +
#  #scale_y_continuous("", limits = c(0.45,.9)) +
#  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
#  theme_minimal() +
#  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
#        axis.title = element_text(family = "serif", size = 48),
#        axis.text.y = element_text(family = "serif", size = 32, hjust = 1), axis.text.x = element_blank(),
#        axis.line.y = element_line(size = 2), axis.line.x = element_blank(),
#        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_blank(),
#        strip.text = element_text(family = "serif", size = 48))#,
##        legend.title = element_blank(),
##        legend.text = element_text(family = "serif", size = 38),
##        legend.position = "none")
##legend.direction = "horizontal", legend.justification = "center")
#
#pb1c <- CN.data %>% dplyr::select("Local", "Min", "BGNfbin1") %>% 
#  pivot_longer(cols = BGNfbin1, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
#  mutate(Index = gsub("fbin1", "", Index)) %>% 
#  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
#  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
#  geom_point(alpha = .25, size = 1.5, show.legend = F) +
#  geom_smooth(linewidth=3, method = "loess", se = F) +
#  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
#  scale_x_discrete(breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
#  scale_color_manual(values = cols)+
#  ggtitle("BGN")+
#  #facet_wrap(Index~., ncol=2, scales = "free") +
#  #scale_y_continuous(limits = c(-40,-25)) +
#  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
#  theme_minimal() +
#  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
#        axis.title = element_text(family = "serif", size = 48),
#        axis.text = element_text(family = "serif", size = 32, hjust = 1),
#        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
#        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2),
#        strip.text = element_text(family = "serif", size = 48))#,
##        legend.title = element_blank(),
##        legend.text = element_text(family = "serif", size = 42),
##        legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center")
#
#pb1d <- CN.data %>% dplyr::select("Local", "Min", "AAfbin1") %>% 
#  pivot_longer(cols = AAfbin1, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
#  mutate(Index = gsub("fbin1", "", Index)) %>% 
#  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
#  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
#  geom_point(alpha = .25, size = 1.5, show.legend = F) +
#  geom_smooth(linewidth=3, method = "loess", se = F) +
#  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
#  scale_x_discrete(breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
#  scale_color_manual(values = cols)+
#  ggtitle("AA") +
#  #facet_wrap(Index~., ncol=2, scales = "free") +
#  #scale_y_continuous("", limits = c(0.35,.6)) +
#  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
#  theme_minimal() +
#  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
#        axis.title = element_text(family = "serif", size = 48),
#        axis.text = element_text(family = "serif", size = 32, hjust = 1),
#        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
#        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2),
#        strip.text = element_text(family = "serif", size = 48))#,
##        legend.title = element_blank(),
##        legend.text = element_text(family = "serif", size = 38),
##        legend.position = "none")
#
#
#
#pb1row <- plot_grid(pb1a + theme(legend.position="none"), 
#                    pb1b + theme(legend.position="none"), 
#                    pb1c + theme(legend.position="none"), 
#                    pb1d + theme(legend.position="none"), 
#                    nrow=2, axis = "b", align = "hv")
#
#legend_b1 <- get_legend(pb1c + theme(legend.title = element_blank(), legend.text = element_text(family = "serif", size = 42),
#                                  legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center"))
#
#
#png("results/overview_byMin2_fbin1.png", width = 33.11024, height = 23.38583, units = "in", res = 300, pointsize = 4)
#plot_grid( pb1row, legend_b1, ncol = 1, rel_heights = c(1, .2))
#dev.off()





#distribution of values between minutes -- considering frequency bin 2 [4 - 12kHz]

p1fbin2 <- CN.data %>% dplyr::select("Local", "Min", "Hfbin2") %>% 
  pivot_longer(cols = Hfbin2, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
  mutate(Index = gsub("fbin2", "", Index)) %>% 
  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
  geom_point(alpha = .25, size = 1.5, show.legend = F) +
  geom_smooth(linewidth=3, method = "loess", se = F) +
  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
  scale_x_discrete("", breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
  scale_color_manual(values = cols)+
  ggtitle("H [4-12kHz]")+xlab("")+ylab("")+
  #facet_wrap(Index~., ncol=2, scales = "free") +
  #scale_y_continuous(limits = c(0.45,0.75)) +
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme_minimal() +
  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
        axis.title = element_text(family = "serif", size = 48),
        axis.text.y = element_text(family = "serif", size = 32, hjust = 1), axis.text.x = element_blank(),
        axis.line.y = element_line(size = 2), axis.line.x = element_blank(),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_blank(),
        strip.text = element_text(family = "serif", size = 48))#,
#        legend.title = element_blank(),
#        legend.text = element_text(family = "serif", size = 38),
#        legend.position = "none")
#legend.direction = "horizontal", legend.justification = "center")

p2fbin2 <- CN.data %>% dplyr::select("Local", "Min", "ADIfbin2") %>% 
  pivot_longer(cols = ADIfbin2, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
  mutate(Index = gsub("fbin2", "", Index)) %>% 
  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
  geom_point(alpha = .25, size = 1.5, show.legend = F) +
  geom_smooth(linewidth=3, method = "loess", se = F) +
  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
  scale_x_discrete("", breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
  scale_color_manual(values = cols)+
  ggtitle("ADI")+xlab("")+ylab("")+
  #facet_wrap(Index~., ncol=2, scales = "free") +
  #scale_y_continuous("", limits = c(0.45,.9)) +
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme_minimal() +
  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
        axis.title = element_text(family = "serif", size = 48),
        axis.text.y = element_text(family = "serif", size = 32, hjust = 1), axis.text.x = element_blank(),
        axis.line.y = element_line(size = 2), axis.line.x = element_blank(),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_blank(),
        strip.text = element_text(family = "serif", size = 48))#,
#        legend.title = element_blank(),
#        legend.text = element_text(family = "serif", size = 38),
#        legend.position = "none")
#legend.direction = "horizontal", legend.justification = "center")

p3fbin2 <- CN.data %>% dplyr::select("Local", "Min", "AEIfbin2") %>% 
  pivot_longer(cols = AEIfbin2, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
  mutate(Index = gsub("fbin2", "", Index)) %>% 
  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
  geom_point(alpha = .25, size = 1.5, show.legend = F) +
  geom_smooth(linewidth=3, method = "loess", se = F) +
  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
  scale_x_discrete("", breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
  scale_color_manual(values = cols)+
  ggtitle("AEI")+xlab("")+ylab("")+
  #facet_wrap(Index~., ncol=2, scales = "free") +
  #scale_y_continuous("", limits = c(0.45,.9)) +
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme_minimal() +
  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
        axis.title = element_text(family = "serif", size = 48),
        axis.text.y = element_text(family = "serif", size = 32, hjust = 1), axis.text.x = element_blank(),
        axis.line.y = element_line(size = 2), axis.line.x = element_blank(),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_blank(),
        strip.text = element_text(family = "serif", size = 48))#,
#        legend.title = element_blank(),
#        legend.text = element_text(family = "serif", size = 38),
#        legend.position = "none")
#legend.direction = "horizontal", legend.justification = "center")

p4fbin2 <- CN.data %>% dplyr::select("Local", "Min", "S2Nfbin2") %>% 
  pivot_longer(cols = S2Nfbin2, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
  mutate(Index = gsub("fbin2", "", Index)) %>% 
  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
  geom_point(alpha = .25, size = 1.5, show.legend = F) +
  geom_smooth(linewidth=3, method = "loess", se = F) +
  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
  scale_x_discrete(breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
  scale_color_manual(values = cols)+
  ggtitle("S2N")+xlab("")+ylab("")+
  #facet_wrap(Index~., ncol=2, scales = "free") +
  #scale_y_continuous(limits = c(-40,-25)) +
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme_minimal() +
  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
        axis.title = element_text(family = "serif", size = 48),
        axis.text.y = element_text(family = "serif", size = 32, hjust = 1), axis.text.x = element_blank(),
        axis.line.y = element_line(size = 2), axis.line.x = element_blank(),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_blank(),
        strip.text = element_text(family = "serif", size = 48))#,
#        legend.title = element_blank(),
#        legend.text = element_text(family = "serif", size = 42),
#        legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center")

p5fbin2 <- CN.data %>% dplyr::select("Local", "Min", "BGNfbin2") %>% 
  pivot_longer(cols = BGNfbin2, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
  mutate(Index = gsub("fbin2", "", Index)) %>% 
  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
  geom_point(alpha = .25, size = 1.5, show.legend = F) +
  geom_smooth(linewidth=3, method = "loess", se = F) +
  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
  scale_x_discrete(breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
  scale_color_manual(values = cols)+
  ggtitle("BGN")+xlab("")+ylab("")+
  #facet_wrap(Index~., ncol=2, scales = "free") +
  #scale_y_continuous(limits = c(-40,-25)) +
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme_minimal() +
  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
        axis.title = element_text(family = "serif", size = 48),
        axis.text.y = element_text(family = "serif", size = 32, hjust = 1), axis.text.x = element_blank(),
        axis.line.y = element_line(size = 2), axis.line.x = element_blank(),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_blank(),
        strip.text = element_text(family = "serif", size = 48))#,
#        legend.title = element_blank(),
#        legend.text = element_text(family = "serif", size = 42),
#        legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center")

p6fbin2 <- CN.data %>% dplyr::select("Local", "Min", "AAfbin2") %>% 
  pivot_longer(cols = AAfbin2, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
  mutate(Index = gsub("fbin2", "", Index)) %>% 
  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
  geom_point(alpha = .25, size = 1.5, show.legend = F) +
  geom_smooth(linewidth=3, method = "loess", se = F) +
  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
  scale_x_discrete(breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
  scale_color_manual(values = cols)+
  ggtitle("AA") +xlab("")+ylab("")+
  #facet_wrap(Index~., ncol=2, scales = "free") +
  #scale_y_continuous("", limits = c(0.35,.6)) +
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme_minimal() +
  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
        axis.title = element_text(family = "serif", size = 48),
        axis.text = element_text(family = "serif", size = 32, hjust = 1),
        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2),
        strip.text = element_text(family = "serif", size = 48))#,
#        legend.title = element_blank(),
#        legend.text = element_text(family = "serif", size = 38),
#        legend.position = "none")



#pfbin2row <- plot_grid(p1fbin2 + theme(legend.position="none"), 
#                       p2fbin2 + theme(legend.position="none"), 
#                       p3fbin2 + theme(legend.position="none"), 
#                       p4fbin2 + theme(legend.position="none"), 
#                       p5fbin2 + theme(legend.position="none"),
#                       p6fbin2 + theme(legend.position="none"), 
#                       p7fbin2 + theme(legend.position="none"), 
#                       nrow=2, axis = "b", align = "hv")
#
#legend_fbin2 <- get_legend(p1fbin2 + theme(legend.title = element_blank(), legend.text = element_text(family = "serif", size = 42),
#                                     legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center"))
#
#
#png("results/overview_byMin2_fbin2.png", width = 33.11024, height = 23.38583, units = "in", res = 300, pointsize = 4)
#plot_grid( pfbin2row, legend_fbin2, ncol = 1, rel_heights = c(1, .2))
#dev.off()

prow <- plot_grid(p1 + theme(legend.position="none"), 
                  p1fbin2 + theme(legend.position="none"), 
                  #p2 + theme(legend.position="none"), 
                  #p2fbin2 + theme(legend.position="none"),
                  p3 + theme(legend.position="none"), 
                  p3fbin2 + theme(legend.position="none"),
                  #p4 + theme(legend.position="none"), 
                  #p4fbin2 + theme(legend.position="none"),
                  p5 + theme(legend.position="none"), 
                  p5fbin2 + theme(legend.position="none"),
                  p6 + theme(legend.position="none"), 
                  p6fbin2 + theme(legend.position="none"),
                  ncol=2, axis = "b", align = "hv")

legend_ <- get_legend(p1 + theme(legend.title = element_blank(), legend.text = element_text(family = "serif", size = 42),
                                   legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center"))



png("results/overview_byMin3.png", width = 23.38583, height = 33.11024, units = "in", res = 300, pointsize = 4)

plot_grid( prow, legend_, ncol = 1, rel_heights = c(1, .05))

dev.off()






##distribution of values between minutes -- considering frequency bin 3 [0.3 - 12kHz]
#
#pb3a <- CN.data %>% dplyr::select("Local", "Min", "Hfbin3") %>% 
#  pivot_longer(cols = Hfbin3, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
#  mutate(Index = gsub("fbin3", "", Index)) %>% 
#  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
#  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
#  geom_point(alpha = .25, size = 1.5, show.legend = F) +
#  geom_smooth(linewidth=3, method = "loess", se = F) +
#  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
#  scale_x_discrete("", breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
#  scale_color_manual(values = cols)+
#  ggtitle("H")+
#  #facet_wrap(Index~., ncol=2, scales = "free") +
#  #scale_y_continuous(limits = c(0.45,0.75)) +
#  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
#  theme_minimal() +
#  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
#        axis.title = element_text(family = "serif", size = 48),
#        axis.text.y = element_text(family = "serif", size = 32, hjust = 1), axis.text.x = element_blank(),
#        axis.line.y = element_line(size = 2), axis.line.x = element_blank(),
#        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_blank(),
#        strip.text = element_text(family = "serif", size = 48))#,
##        legend.title = element_blank(),
##        legend.text = element_text(family = "serif", size = 38),
##        legend.position = "none")
##legend.direction = "horizontal", legend.justification = "center")
#
#pb3b <- CN.data %>% dplyr::select("Local", "Min", "AEIfbin3") %>% 
#  pivot_longer(cols = AEIfbin3, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
#  mutate(Index = gsub("fbin3", "", Index)) %>% 
#  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
#  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
#  geom_point(alpha = .25, size = 1.5, show.legend = F) +
#  geom_smooth(linewidth=3, method = "loess", se = F) +
#  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
#  scale_x_discrete("", breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
#  scale_color_manual(values = cols)+
#  ggtitle("AEI")+
#  #facet_wrap(Index~., ncol=2, scales = "free") +
#  #scale_y_continuous("", limits = c(0.45,.9)) +
#  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
#  theme_minimal() +
#  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
#        axis.title = element_text(family = "serif", size = 48),
#        axis.text.y = element_text(family = "serif", size = 32, hjust = 1), axis.text.x = element_blank(),
#        axis.line.y = element_line(size = 2), axis.line.x = element_blank(),
#        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_blank(),
#        strip.text = element_text(family = "serif", size = 48))#,
##        legend.title = element_blank(),
##        legend.text = element_text(family = "serif", size = 38),
##        legend.position = "none")
##legend.direction = "horizontal", legend.justification = "center")
#
#pb3c <- CN.data %>% dplyr::select("Local", "Min", "BGNfbin3") %>% 
#  pivot_longer(cols = BGNfbin3, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
#  mutate(Index = gsub("fbin3", "", Index)) %>% 
#  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
#  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
#  geom_point(alpha = .25, size = 1.5, show.legend = F) +
#  geom_smooth(linewidth=3, method = "loess", se = F) +
#  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
#  scale_x_discrete(breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
#  scale_color_manual(values = cols)+
#  ggtitle("BGN")+
#  #facet_wrap(Index~., ncol=2, scales = "free") +
#  #scale_y_continuous(limits = c(-40,-25)) +
#  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
#  theme_minimal() +
#  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
#        axis.title = element_text(family = "serif", size = 48),
#        axis.text = element_text(family = "serif", size = 32, hjust = 1),
#        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
#        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2),
#        strip.text = element_text(family = "serif", size = 48))#,
##        legend.title = element_blank(),
##        legend.text = element_text(family = "serif", size = 42),
##        legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center")
#
#pb3d <- CN.data %>% dplyr::select("Local", "Min", "AAfbin3") %>% 
#  pivot_longer(cols = AAfbin3, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
#  mutate(Index = gsub("fbin3", "", Index)) %>% 
#  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
#  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
#  geom_point(alpha = .25, size = 1.5, show.legend = F) +
#  geom_smooth(linewidth=3, method = "loess", se = F) +
#  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
#  scale_x_discrete(breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
#  scale_color_manual(values = cols)+
#  ggtitle("AA") +
#  #facet_wrap(Index~., ncol=2, scales = "free") +
#  #scale_y_continuous("", limits = c(0.35,.6)) +
#  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
#  theme_minimal() +
#  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
#        axis.title = element_text(family = "serif", size = 48),
#        axis.text = element_text(family = "serif", size = 32, hjust = 1),
#        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
#        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2),
#        strip.text = element_text(family = "serif", size = 48))#,
##        legend.title = element_blank(),
##        legend.text = element_text(family = "serif", size = 38),
##        legend.position = "none")
#
#
#
#pb3row <- plot_grid(pb3a + theme(legend.position="none"), 
#                    pb3b + theme(legend.position="none"), 
#                    pb3c + theme(legend.position="none"), 
#                    pb3d + theme(legend.position="none"), 
#                    nrow=2, axis = "b", align = "hv")
#
#legend_b3 <- get_legend(pb3c + theme(legend.title = element_blank(), legend.text = element_text(family = "serif", size = 42),
#                                     legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center"))
#
#
#png("results/overview_byMin2_fbin3.png", width = 33.11024, height = 23.38583, units = "in", res = 300, pointsize = 4)
#plot_grid( pb3row, legend_b3, ncol = 1, rel_heights = c(1, .2))
#dev.off()
#
#
#
#
#
##distribution of values between minutes -- considering frequency bin 4 [12 - 22.1kHz]
#
#pb4a <- CN.data %>% dplyr::select("Local", "Min", "Hfbin4") %>% 
#  pivot_longer(cols = Hfbin4, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
#  mutate(Index = gsub("fbin4", "", Index)) %>% 
#  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
#  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
#  geom_point(alpha = .25, size = 1.5, show.legend = F) +
#  geom_smooth(linewidth=3, method = "loess", se = F) +
#  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
#  scale_x_discrete("", breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
#  scale_color_manual(values = cols)+
#  ggtitle("H")+
#  #facet_wrap(Index~., ncol=2, scales = "free") +
#  #scale_y_continuous(limits = c(0.45,0.75)) +
#  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
#  theme_minimal() +
#  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
#        axis.title = element_text(family = "serif", size = 48),
#        axis.text.y = element_text(family = "serif", size = 32, hjust = 1), axis.text.x = element_blank(),
#        axis.line.y = element_line(size = 2), axis.line.x = element_blank(),
#        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_blank(),
#        strip.text = element_text(family = "serif", size = 48))#,
##        legend.title = element_blank(),
##        legend.text = element_text(family = "serif", size = 38),
##        legend.position = "none")
##legend.direction = "horizontal", legend.justification = "center")
#
#pb4b <- CN.data %>% dplyr::select("Local", "Min", "AEIfbin4") %>% 
#  pivot_longer(cols = AEIfbin4, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
#  mutate(Index = gsub("fbin4", "", Index)) %>% 
#  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
#  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
#  geom_point(alpha = .25, size = 1.5, show.legend = F) +
#  geom_smooth(linewidth=3, method = "loess", se = F) +
#  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
#  scale_x_discrete("", breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
#  scale_color_manual(values = cols)+
#  ggtitle("AEI")+
#  #facet_wrap(Index~., ncol=2, scales = "free") +
#  #scale_y_continuous("", limits = c(0.45,.9)) +
#  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
#  theme_minimal() +
#  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
#        axis.title = element_text(family = "serif", size = 48),
#        axis.text.y = element_text(family = "serif", size = 32, hjust = 1), axis.text.x = element_blank(),
#        axis.line.y = element_line(size = 2), axis.line.x = element_blank(),
#        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_blank(),
#        strip.text = element_text(family = "serif", size = 48))#,
##        legend.title = element_blank(),
##        legend.text = element_text(family = "serif", size = 38),
##        legend.position = "none")
##legend.direction = "horizontal", legend.justification = "center")
#
#pb4c <- CN.data %>% dplyr::select("Local", "Min", "BGNfbin4") %>% 
#  pivot_longer(cols = BGNfbin4, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
#  mutate(Index = gsub("fbin4", "", Index)) %>% 
#  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
#  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
#  geom_point(alpha = .25, size = 1.5, show.legend = F) +
#  geom_smooth(linewidth=3, method = "loess", se = F) +
#  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
#  scale_x_discrete(breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
#  scale_color_manual(values = cols)+
#  ggtitle("BGN")+
#  #facet_wrap(Index~., ncol=2, scales = "free") +
#  #scale_y_continuous(limits = c(-40,-25)) +
#  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
#  theme_minimal() +
#  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
#        axis.title = element_text(family = "serif", size = 48),
#        axis.text = element_text(family = "serif", size = 32, hjust = 1),
#        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
#        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2),
#        strip.text = element_text(family = "serif", size = 48))#,
##        legend.title = element_blank(),
##        legend.text = element_text(family = "serif", size = 42),
##        legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center")
#
#pb4d <- CN.data %>% dplyr::select("Local", "Min", "AAfbin4") %>% 
#  pivot_longer(cols = AAfbin4, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
#  mutate(Index = gsub("fbin4", "", Index)) %>% 
#  #mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
#  ggplot(aes(Min, Values, group=Local,  colour = Local)) +
#  geom_point(alpha = .25, size = 1.5, show.legend = F) +
#  geom_smooth(linewidth=3, method = "loess", se = F) +
#  geom_vline(aes(xintercept=39), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=51), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=111), linewidth=1.2, linetype='dashed', color = "gray55") +
#  geom_vline(aes(xintercept=123), linewidth=1.2, linetype='dashed', color = "gray55") +
#  scale_x_discrete(breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
#  scale_color_manual(values = cols)+
#  ggtitle("AA") +
#  #facet_wrap(Index~., ncol=2, scales = "free") +
#  #scale_y_continuous("", limits = c(0.35,.6)) +
#  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
#  theme_minimal() +
#  theme(title = element_text(family = "serif", size = 48), plot.title = element_text(hjust = 0.5),
#        axis.title = element_text(family = "serif", size = 48),
#        axis.text = element_text(family = "serif", size = 32, hjust = 1),
#        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
#        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2),
#        strip.text = element_text(family = "serif", size = 48))#,
##        legend.title = element_blank(),
##        legend.text = element_text(family = "serif", size = 38),
##        legend.position = "none")
#
#
#
#pb4row <- plot_grid(pb4a + theme(legend.position="none"), 
#                    pb4b + theme(legend.position="none"), 
#                    pb4c + theme(legend.position="none"), 
#                    pb4d + theme(legend.position="none"), 
#                    nrow=2, axis = "b", align = "hv")
#
#legend_b4 <- get_legend(pb4c + theme(legend.title = element_blank(), legend.text = element_text(family = "serif", size = 42),
#                                     legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center"))
#
#
#png("results/overview_byMin2_fbin4.png", width = 33.11024, height = 23.38583, units = "in", res = 300, pointsize = 4)
#plot_grid( pb4row, legend_b4, ncol = 1, rel_heights = c(1, .2))
#dev.off()
#
##



#time-series graph
#brk1 <- c(max(as.POSIXct(CN.data.total[CN.data.total$Year==2019,"Date"])),
#         min(as.POSIXct(CN.data.total[CN.data.total$Year==2022,"Date"])))
#
#brk2 <- c(max(as.POSIXct(CN.data.total[grep("202204", CN.data.total$ldt.index),"Date"])),
#          min(as.POSIXct(CN.data.total[grep("202205", CN.data.total$ldt.index),"Date"])))
#
#ggplot(CN.data.total, aes(x=as.POSIXct(Date), y=HTotal, color=Local, group=Year)) +
#  geom_line(show.legend=F) +
#  facet_wrap(Local~., ncol = 1) +
#  scale_x_datetime(date_breaks = "1 days" , date_labels = "%d") + 
#  scale_x_break(breaks = brk1, scales = .75) +
#  scale_x_break(breaks = brk2, scales = .5) +
#  #ylim(.75,1) +
#  xlab("day") + ylab("Acustic entropy")









##### modelling the effect of environmental traits on soundscape index ####
## preparing data -- regional var
#CN.data.regional %>% group_by(Local, time_period) %>% summarise(N=n()) %>% arrange(N)
CN.data.regional <- CN.data.regional[!is.na(CN.data.regional$time_period),]

## scaling environmental variables [z-score]
CN.data.regional$z_alt <- scale(CN.data.regional$alt)
CN.data.regional$z_distwater <- scale(CN.data.regional$distwater)
CN.data.regional$z_distedge <- scale(CN.data.regional$distedge)
CN.data.regional$z_distcanga <- scale(CN.data.regional$distcanga)
CN.data.regional$z_distminning <- scale(CN.data.regional$distminning)
CN.data.regional$z_treeheight <- scale(CN.data.regional$treeheight)
CN.data.regional$z_NDVI <- scale(CN.data.regional$NDVI)
CN.data.regional$z_precipitation <- scale(as.numeric(CN.data.regional$precipitation))
CN.data.regional$z_temperature <- scale(as.numeric(CN.data.regional$temperature))
CN.data.regional$z_humidity <- scale(as.numeric(CN.data.regional$humidity))
CN.data.regional$z_wind <- scale(as.numeric(CN.data.regional$wind))



#subsetting
set.seed(123)
CN.data.regional.modelfit <- CN.data.regional %>% group_by(Local, time_period) %>% sample_n(size = 48, replace = T) %>% ungroup()


#creating time step variable to account for temporal autocorrelation
CN.data.regional.modelfit <- (CN.data.regional.modelfit
                              %>% arrange(Date)
                              %>% group_by(Year)
                              %>% mutate(times = lubridate::decimal_date(Date) %% 1)
                              %>% ungroup()
                              )


CN.data.regional.modelfit <- (CN.data.regional.modelfit
                              %>% mutate(times = glmmTMB::numFactor(times))
                              )



#Response: Acoustic entropy (H) on total frequency [0.3 - 22.1kHz]
#hist(CN.data.regional.modelfit$HTotal); range(CN.data.regional.modelfit$HTotal)
#plot(family.test <- ecdf(CN.data.regional.modelfit$HTotal))
#plot(fitdist(CN.data.regional.modelfit$HTotal, "beta"))


#modH.null <- glmmTMB(HTotal ~ 1, family=beta_family(), CN.data.regional.modelfit)

options(na.action = "na.fail")
modH.full <- glmmTMB(HTotal ~ z_alt + z_distwater + z_distedge + z_distcanga + z_distminning + z_treeheight +
                       z_NDVI + z_precipitation + z_temperature + z_humidity + z_wind + (1|Local:time_period),
                     family=beta_family(), CN.data.regional.modelfit)

check_collinearity(modH.full)


modH.full.b <- glmmTMB(HTotal ~ z_alt + z_distwater + z_distedge + z_treeheight + z_NDVI + z_precipitation + 
                         z_temperature + z_humidity + z_wind + (1|Local:time_period) + ar1(times - 1|Year),
                       family=beta_family(), CN.data.regional.modelfit)

summary(modH.full.b)
#r.squaredGLMM(modH.full.b)

test.modH.full.b <- dredge(modH.full.b)
#subset(test.modH.full.b, delta < 2)
summary(get.models(test.modH.full.b, 1)[[1]])
#r.squaredGLMM(get.models(test.modH.full.b, 1)[[1]])

sink("results/supplementar/HTotal_model_summary.csv")
print(summary(get.models(test.modH.full.b, 1)[[1]]))
sink()

modH.sel.table<-as.data.frame(test.modH.full.b)[1:5,12:16]
modH.sel.table[,2:3]<- round(modH.sel.table[,2:3],2)
modH.sel.table[,4:5]<- round(modH.sel.table[,4:5],3)
names(modH.sel.table)[1] = "K"
modH.sel.table$Model<-NA
for(i in 1:nrow(modH.sel.table)) modH.sel.table$Model[i]<-as.character(formula(get.models(test.modH.full.b, i)[[1]]))[3]
modH.sel.table<-modH.sel.table[,c(6,1,2,3,4,5)]


modH.full.b.simres <- simulateResiduals(get.models(test.modH.full.b, 1)[[1]])
plot(modH.full.b.simres)
testDispersion(modH.full.b.simres)

modH.full.b.alleffects <- allEffects(get.models(test.modH.full.b, 1)[[1]], rescale.axis=F) #, residuals = T
plot(modH.full.b.alleffects)


modH.full.b.t <- broom.mixed::tidy(get.models(test.modH.full.b, 1)[[1]], effects = "fixed", conf.int = TRUE)
modH.full.b.t <- transform(modH.full.b.t, term = sprintf("%s.%s", component, term))
modH.full.b.t$group <- "H [0.3-22.1kHz]"
modH.full.b.t$col <- ifelse(modH.full.b.t$conf.low < 0 & modH.full.b.t$conf.high < 0, "#e65a49", 
                            ifelse(modH.full.b.t$conf.low > 0 & modH.full.b.t$conf.high > 0, "#55c0f0", "#dbdfde"))



#Response: Acoustic entropy (H) on frequency bin 2 [4 - 12kHz]
#hist(CN.data.regional.modelfit$Hfbin2); range(CN.data.regional.modelfit$Hfbin2)
#plot(family.test <- ecdf(CN.data.regional.modelfit$Hfbin2))
#plot(fitdist(CN.data.regional.modelfit$Hfbin2, "beta"))


#modHfbin2.null <- glmmTMB(Hfbin2 ~ 1, family=beta_family(), CN.data.regional.modelfit)

modHfbin2.full <- glmmTMB(Hfbin2 ~ z_alt + z_distwater + z_distedge + z_distcanga + z_distminning + z_treeheight +
                       z_NDVI + z_precipitation + z_temperature + z_humidity + z_wind + (1|Local:time_period),
                     family=beta_family(), CN.data.regional.modelfit)

check_collinearity(modHfbin2.full)


modHfbin2.full.b <- glmmTMB(Hfbin2 ~ z_alt + z_distwater + z_distedge + z_treeheight + z_NDVI + z_precipitation + 
                         z_temperature + z_humidity + z_wind + (1|Local:time_period) + ar1(times - 1|Year),
                       family=beta_family(), CN.data.regional.modelfit)

summary(modHfbin2.full.b)
#r.squaredGLMM(modHfbin2.full.b)

test.modHfbin2.full.b <- dredge(modHfbin2.full.b)
#subset(test.modHfbin2.full.b, delta < 2)
summary(get.models(test.modHfbin2.full.b, 1)[[1]])
#r.squaredGLMM(get.models(test.modHfbin2.full.b, 1)[[1]])

sink("results/supplementar/Hfbin2_model_summary.csv")
print(summary(get.models(test.modHfbin2.full.b, 1)[[1]]))
sink()

modHfbin2.sel.table<-as.data.frame(test.modHfbin2.full.b)[1:5,12:16]
modHfbin2.sel.table[,2:3]<- round(modHfbin2.sel.table[,2:3],2)
modHfbin2.sel.table[,4:5]<- round(modHfbin2.sel.table[,4:5],3)
names(modHfbin2.sel.table)[1] = "K"
modHfbin2.sel.table$Model<-NA
for(i in 1:nrow(modHfbin2.sel.table)) modHfbin2.sel.table$Model[i]<-as.character(formula(get.models(test.modHfbin2.full.b, i)[[1]]))[3]
modHfbin2.sel.table<-modHfbin2.sel.table[,c(6,1,2,3,4,5)]


modHfbin2.full.b.simres <- simulateResiduals(get.models(test.modHfbin2.full.b, 1)[[1]])
plot(modHfbin2.full.b.simres)
testDispersion(modHfbin2.full.b.simres)

modHfbin2.full.b.alleffects <- allEffects(get.models(test.modHfbin2.full.b, 1)[[1]], rescale.axis=F) #, residuals = T
plot(modHfbin2.full.b.alleffects)


modHfbin2.full.b.t <- broom.mixed::tidy(get.models(test.modHfbin2.full.b, 1)[[1]], effects = "fixed", conf.int = TRUE)
modHfbin2.full.b.t <- transform(modHfbin2.full.b.t, term = sprintf("%s.%s", component, term))
modHfbin2.full.b.t$group <- "H [4-12kHz]"
modHfbin2.full.b.t$col <- ifelse(modHfbin2.full.b.t$conf.low < 0 & modHfbin2.full.b.t$conf.high < 0, "#e65a49", 
                            ifelse(modHfbin2.full.b.t$conf.low > 0 & modHfbin2.full.b.t$conf.high > 0, "#55c0f0", "#dbdfde"))



#Response: Acoustic diversity index (ADI) on total frequency [0.3 - 22.1kHz]
CN.data.regional.modelfit$ADITotal.adj <- 1-log(CN.data.regional.modelfit$ADITotal)
#hist(CN.data.regional.modelfit$ADITotal.adj); range(CN.data.regional.modelfit$ADITotal.adj)
#plot(family.test <- ecdf(CN.data.regional.modelfit$ADITotal.adj))
#plot(fitdist(CN.data.regional.modelfit$ADITotal.adj, "lnorm"))


#modADI.null <- glmmTMB(ADITotal ~ 1, CN.data.regional.modelfit)

modADI.full <- glmmTMB(ADITotal.adj ~ z_alt + z_distwater + z_distedge + z_distcanga + z_distminning + z_treeheight +
                       z_NDVI + z_precipitation + z_temperature + z_humidity + z_wind + (1|Local:time_period),
                     #family=beta_family(), 
                     CN.data.regional.modelfit)

check_collinearity(modADI.full)


modADI.full.b <- glmmTMB(ADITotal.adj ~ z_alt + z_distwater + z_distedge + z_treeheight + z_NDVI + z_precipitation + 
                         z_temperature + z_humidity + z_wind + (1|Local:time_period) + ar1(times - 1|Year), ziformula = ~1,
                       #family=beta_family(), 
                       CN.data.regional.modelfit)

#fixef(modADI.full.b)
#modADI.full.c <- update(modADI.full.b, control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))

summary(modADI.full.b)
#r.squaredGLMM(modADI.full.b)

test.modADI.full.b <- dredge(modADI.full.b)
#subset(test.modADI.full.b, delta < 2)
summary(get.models(test.modADI.full.b, 1)[[1]])
#r.squaredGLMM(get.models(test.modADI.full.b, 1)[[1]])

sink("results/supplementar/ADITotal_model_summary.csv")
print(summary(get.models(test.modADI.full.b, 1)[[1]]))
sink()

modADI.sel.table<-as.data.frame(test.modADI.full.b)[1:5,12:16]
modADI.sel.table[,2:3]<- round(modADI.sel.table[,2:3],2)
modADI.sel.table[,4:5]<- round(modADI.sel.table[,4:5],3)
names(modADI.sel.table)[1] = "K"
modADI.sel.table$Model<-NA
for(i in 1:nrow(modADI.sel.table)) modADI.sel.table$Model[i]<-as.character(formula(get.models(test.modADI.full.b, i)[[1]]))[3]
modADI.sel.table<-modADI.sel.table[,c(6,1,2,3,4,5)]


modADI.full.b.simres <- simulateResiduals(get.models(test.modADI.full.b, 1)[[1]])
plot(modADI.full.b.simres)
testDispersion(modADI.full.b.simres)

modADI.full.b.alleffects <- allEffects(get.models(test.modADI.full.b, 1)[[1]], rescale.axis=F) #, residuals = T
plot(modADI.full.b.alleffects)


modADI.full.b.t <- broom.mixed::tidy(get.models(test.modADI.full.b, 1)[[1]], effects = "fixed", conf.int = TRUE)
modADI.full.b.t <- transform(modADI.full.b.t, term = sprintf("%s.%s", component, term))
modADI.full.b.t$group <- "ADI [0.3-22.1kHz]"
modADI.full.b.t$col <- ifelse(modADI.full.b.t$conf.low < 0 & modADI.full.b.t$conf.high < 0, "#e65a49", 
                            ifelse(modADI.full.b.t$conf.low > 0 & modADI.full.b.t$conf.high > 0, "#55c0f0", "#dbdfde"))



#Response: Acoustic diversity index (ADI) on total frequency [4 - 12kHz]
CN.data.regional.modelfit$ADIfbin2.adj <- 1-log(CN.data.regional.modelfit$ADIfbin2)
#hist(CN.data.regional.modelfit$ADIfbin2.adj); range(CN.data.regional.modelfit$ADIfbin2.adj)
#plot(family.test <- ecdf(CN.data.regional.modelfit$ADIfbin2.adj))
#plot(fitdist(CN.data.regional.modelfit$ADIfbin2.adj, "lnorm"))


#modADIfbin2.null <- glmmTMB(ADIfbin2.adj ~ 1, family=beta_family(), CN.data.regional.modelfit)

modADIfbin2.full <- glmmTMB(ADIfbin2.adj ~ z_alt + z_distwater + z_distedge + z_distcanga + z_distminning + z_treeheight +
                            z_NDVI + z_precipitation + z_temperature + z_humidity + z_wind + (1|Local:time_period),
                          #family=beta_family(), 
                          CN.data.regional.modelfit)

check_collinearity(modADIfbin2.full)


modADIfbin2.full.b <- glmmTMB(ADIfbin2.adj ~ z_alt + z_distwater + z_distedge + z_treeheight + z_NDVI + z_precipitation + 
                              z_temperature + z_humidity + z_wind + (1|Local:time_period) + ar1(times - 1|Year), ziformula = ~1,
                            #family=beta_family(), 
                            CN.data.regional.modelfit)

summary(modADIfbin2.full.b)
#r.squaredGLMM(modADIfbin2.full.b)

test.modADIfbin2.full.b <- dredge(modADIfbin2.full.b)
#subset(test.modADIfbin2.full.b, delta < 4)
summary(get.models(test.modADIfbin2.full.b, 1)[[1]])
#r.squaredGLMM(get.models(test.modADIfbin2.full.b, 1)[[1]])

sink("results/supplementar/ADIfbin2_model_summary.csv")
print(summary(get.models(test.modADIfbin2.full.b, 1)[[1]]))
sink()

modADIfbin2.sel.table<-as.data.frame(test.modADIfbin2.full.b)[1:5,12:16]
modADIfbin2.sel.table[,2:3]<- round(modADIfbin2.sel.table[,2:3],2)
modADIfbin2.sel.table[,4:5]<- round(modADIfbin2.sel.table[,4:5],3)
names(modADIfbin2.sel.table)[1] = "K"
modADIfbin2.sel.table$Model<-NA
for(i in 1:nrow(modADIfbin2.sel.table)) modADIfbin2.sel.table$Model[i]<-as.character(formula(get.models(test.modADIfbin2.full.b, i)[[1]]))[3]
modADIfbin2.sel.table<-modADIfbin2.sel.table[,c(6,1,2,3,4,5)]


modADIfbin2.full.b.simres <- simulateResiduals(get.models(test.modADIfbin2.full.b, 1)[[1]])
plot(modADIfbin2.full.b.simres)
testDispersion(modADIfbin2.full.b.simres)

modADIfbin2.full.b.alleffects <- allEffects(get.models(test.modADIfbin2.full.b, 1)[[1]], rescale.axis=F) #, residuals = T
plot(modADIfbin2.full.b.alleffects)


modADIfbin2.full.b.t <- broom.mixed::tidy(get.models(test.modADIfbin2.full.b, 1)[[1]], effects = "fixed", conf.int = TRUE)
modADIfbin2.full.b.t <- transform(modADIfbin2.full.b.t, term = sprintf("%s.%s", component, term))
modADIfbin2.full.b.t$group <- "ADI [4-12kHz]"
modADIfbin2.full.b.t$col <- ifelse(modADIfbin2.full.b.t$conf.low < 0 & modADIfbin2.full.b.t$conf.high < 0, "#e65a49", 
                                 ifelse(modADIfbin2.full.b.t$conf.low > 0 & modADIfbin2.full.b.t$conf.high > 0, "#55c0f0", "#dbdfde"))



#Response: Acoustic entropy (AEI) on total frequency [0.3 - 22.1kAEIz]
CN.data.regional.modelfit$AEITotal.adj <- 1-log(CN.data.regional.modelfit$AEITotal)
#hist(CN.data.regional.modelfit$AEITotal.adj); range(CN.data.regional.modelfit$AEITotal.adj)
#plot(family.test <- ecdf(CN.data.regional.modelfit$AEITotal.adj))
#plot(fitdist(CN.data.regional.modelfit$AEITotal.adj, "beta"))


#modAEI.null <- glmmTMB(AEITotal.adj ~ 1, CN.data.regional.modelfit)

modAEI.full <- glmmTMB(AEITotal.adj ~ z_alt + z_distwater + z_distedge + z_distcanga + z_distminning + z_treeheight +
                         z_NDVI + z_precipitation + z_temperature + z_humidity + z_wind + (1|Local:time_period),
                       #family=beta_family(), 
                       CN.data.regional.modelfit)

check_collinearity(modAEI.full)


modAEI.full.b <- glmmTMB(AEITotal.adj ~ z_alt + z_distwater + z_distedge + z_treeheight + z_NDVI + z_precipitation + 
                           z_temperature + z_humidity + z_wind + (1|Local:time_period) + ar1(times - 1|Year),
                         #family=beta_family(), 
                         CN.data.regional.modelfit)

summary(modAEI.full.b)
#r.squaredGLMM(modAEI.full.b)

test.modAEI.full.b <- dredge(modAEI.full.b)
#subset(test.modAEI.full.b, delta < 2)
summary(get.models(test.modAEI.full.b, 1)[[1]])
#r.squaredGLMM(get.models(test.modAEI.full.b, 1)[[1]])

sink("results/supplementar/AEITotal_model_summary.csv")
print(summary(get.models(test.modAEI.full.b, 1)[[1]]))
sink()

modAEI.sel.table<-as.data.frame(test.modAEI.full.b)[1:5,12:16]
modAEI.sel.table[,2:3]<- round(modAEI.sel.table[,2:3],2)
modAEI.sel.table[,4:5]<- round(modAEI.sel.table[,4:5],3)
names(modAEI.sel.table)[1] = "K"
modAEI.sel.table$Model<-NA
for(i in 1:nrow(modAEI.sel.table)) modAEI.sel.table$Model[i]<-as.character(formula(get.models(test.modAEI.full.b, i)[[1]]))[3]
modAEI.sel.table<-modAEI.sel.table[,c(6,1,2,3,4,5)]


modAEI.full.b.simres <- simulateResiduals(get.models(test.modAEI.full.b, 1)[[1]])
plot(modAEI.full.b.simres)
testDispersion(modAEI.full.b.simres)

modAEI.full.b.alleffects <- allEffects(get.models(test.modAEI.full.b, 1)[[1]], rescale.axis=F) #, residuals = T
plot(modAEI.full.b.alleffects)


modAEI.full.b.t <- broom.mixed::tidy(get.models(test.modAEI.full.b, 1)[[1]], effects = "fixed", conf.int = TRUE)
modAEI.full.b.t <- transform(modAEI.full.b.t, term = sprintf("%s.%s", component, term))
modAEI.full.b.t$group <- "AEI [0.3-22.1kHz]"
modAEI.full.b.t$col <- ifelse(modAEI.full.b.t$conf.low < 0 & modAEI.full.b.t$conf.high < 0, "#e65a49", 
                              ifelse(modAEI.full.b.t$conf.low > 0 & modAEI.full.b.t$conf.high > 0, "#55c0f0", "#dbdfde"))



#Response: Acoustic entropy (AEI) on frequency bin 2 [4 - 12kAEIz]
CN.data.regional.modelfit$AEIfbin2.adj <- 1-log(CN.data.regional.modelfit$AEIfbin2)
#hist(CN.data.regional.modelfit$AEIfbin2.adj); range(CN.data.regional.modelfit$AEIfbin2.adj)
#plot(family.test <- ecdf(CN.data.regional.modelfit$AEIfbin2.adj))
#plot(fitdist(CN.data.regional.modelfit$AEIfbin2.adj, "lnorm"))


#modAEIfbin2.null <- glmmTMB(AEIfbin2 ~ 1, family=beta_family(), CN.data.regional.modelfit)

modAEIfbin2.full <- glmmTMB(AEIfbin2.adj ~ z_alt + z_distwater + z_distedge + z_distcanga + z_distminning + z_treeheight +
                              z_NDVI + z_precipitation + z_temperature + z_humidity + z_wind + (1|Local:time_period),
                            #family=beta_family(), 
                            CN.data.regional.modelfit)

check_collinearity(modAEIfbin2.full)


modAEIfbin2.full.b <- glmmTMB(AEIfbin2.adj ~ z_alt + z_distwater + z_distedge + z_treeheight + z_NDVI + z_precipitation + 
                                z_temperature + z_humidity + z_wind + (1|Local:time_period) + ar1(times - 1|Year),
                              #family=beta_family(), 
                              CN.data.regional.modelfit)

summary(modAEIfbin2.full.b)
#r.squaredGLMM(modAEIfbin2.full.b)

test.modAEIfbin2.full.b <- dredge(modAEIfbin2.full.b)
#subset(test.modAEIfbin2.full.b, delta < 4)
summary(get.models(test.modAEIfbin2.full.b, 1)[[1]])
#r.squaredGLMM(get.models(test.modAEIfbin2.full.b, 1)[[1]])

sink("results/supplementar/AEIfbin2_model_summary.csv")
print(summary(get.models(test.modAEIfbin2.full.b, 1)[[1]]))
sink()

modAEIfbin2.sel.table<-as.data.frame(test.modAEIfbin2.full.b)[1:5,12:16]
modAEIfbin2.sel.table[,2:3]<- round(modAEIfbin2.sel.table[,2:3],2)
modAEIfbin2.sel.table[,4:5]<- round(modAEIfbin2.sel.table[,4:5],3)
names(modAEIfbin2.sel.table)[1] = "K"
modAEIfbin2.sel.table$Model<-NA
for(i in 1:nrow(modAEIfbin2.sel.table)) modAEIfbin2.sel.table$Model[i]<-as.character(formula(get.models(test.modAEIfbin2.full.b, i)[[1]]))[3]
modAEIfbin2.sel.table<-modAEIfbin2.sel.table[,c(6,1,2,3,4,5)]


modAEIfbin2.full.b.simres <- simulateResiduals(get.models(test.modAEIfbin2.full.b, 1)[[1]])
plot(modAEIfbin2.full.b.simres)
testDispersion(modAEIfbin2.full.b.simres)

modAEIfbin2.full.b.alleffects <- allEffects(get.models(test.modAEIfbin2.full.b, 1)[[1]], rescale.axis=F) #, residuals = T
plot(modAEIfbin2.full.b.alleffects)


modAEIfbin2.full.b.t <- broom.mixed::tidy(get.models(test.modAEIfbin2.full.b, 1)[[1]], effects = "fixed", conf.int = TRUE)
modAEIfbin2.full.b.t <- transform(modAEIfbin2.full.b.t, term = sprintf("%s.%s", component, term))
modAEIfbin2.full.b.t$group <- "AEI [4-12kHz]"
modAEIfbin2.full.b.t$col <- ifelse(modAEIfbin2.full.b.t$conf.low < 0 & modAEIfbin2.full.b.t$conf.high < 0, "#e65a49", 
                                   ifelse(modAEIfbin2.full.b.t$conf.low > 0 & modAEIfbin2.full.b.t$conf.high > 0, "#55c0f0", "#dbdfde"))



#Response: Acoustic entropy (BGN) on total frequency [0.3 - 22.1kBGNz]
CN.data.regional.modelfit$BGNTotal.adj <- (CN.data.regional.modelfit$BGNTotal)*(-1)
#hist(CN.data.regional.modelfit$BGNTotal.adj); range(CN.data.regional.modelfit$BGNTotal.adj)
#plot(family.test <- ecdf(CN.data.regional.modelfit$BGNTotal.adj))
#plot(fitdist(CN.data.regional.modelfit$BGNTotal.adj, "pois"))


#modBGN.null <- glmmTMB(BGNTotal.adj ~ 1, family=genpois(link = "log"), CN.data.modelfit)


modBGN.full <- glmmTMB(BGNTotal.adj ~ z_alt + z_distwater + z_distedge + z_distcanga + z_distminning + z_treeheight +
                         z_NDVI + z_precipitation + z_temperature + z_humidity + z_wind + (1|Local:time_period),
                       family=genpois(link = "log"), CN.data.regional.modelfit)

check_collinearity(modBGN.full)

modBGN.full.b <- glmmTMB(BGNTotal.adj ~ z_alt + z_distwater + z_distedge + z_treeheight + z_NDVI + z_precipitation + 
                           z_temperature + z_humidity + z_wind + (1|Local:time_period) + ar1(times - 1|Year),
                         family=genpois(link = "log"), CN.data.regional.modelfit)

summary(modBGN.full.b)
#r.squaredGLMM(modBGN.full.b)

test.modBGN.full.b <- dredge(modBGN.full.b)
#subset(test.modBGN.full.b, delta < 4)
summary(get.models(test.modBGN.full.b, 1)[[1]])
#r.squaredGLMM(get.models(test.modBGN.full.b, 1)[[1]])

sink("results/supplementar/BGNTotal_model_summary.csv")
print(summary(get.models(test.modBGN.full.b, 1)[[1]]))
sink()

modBGN.sel.table<-as.data.frame(test.modBGN.full.b)[1:5,12:16]
modBGN.sel.table[,2:3]<- round(modBGNsel.table[,2:3],2)
modBGN.sel.table[,4:5]<- round(modBGN.sel.table[,4:5],3)
names(modBGN.sel.table)[1] = "K"
modBGN.sel.table$Model<-NA
for(i in 1:nrow(modBGN.sel.table)) modBGN.sel.table$Model[i]<-as.character(formula(get.models(test.modBGN.full.b, i)[[1]]))[3]
modBGN.sel.table<-modBGN.sel.table[,c(6,1,2,3,4,5)]


modBGN.full.b.simres <- simulateResiduals(get.models(test.modBGN.full.b, 1)[[1]])
plot(modBGN.full.b.simres)
testDispersion(modBGN.full.b.simres)

modBGN.full.b.alleffects <- allEffects(get.models(test.modBGN.full.b, 1)[[1]], rescale.axis=F) #, residuals = T
plot(modBGN.full.b.alleffects)


modBGN.full.b.t <- broom.mixed::tidy(get.models(test.modBGN.full.b, 1)[[1]], effects = "fixed", conf.int = TRUE)
modBGN.full.b.t <- transform(modBGN.full.b.t, term = sprintf("%s.%s", component, term))
modBGN.full.b.t$group <- "BGN [0.3-22.1kHz]"
modBGN.full.b.t$col <- ifelse(modBGN.full.b.t$conf.low < 0 & modBGN.full.b.t$conf.high < 0, "#e65a49", 
                              ifelse(modBGN.full.b.t$conf.low > 0 & modBGN.full.b.t$conf.high > 0, "#55c0f0", "#dbdfde"))



#Response: Acoustic entropy (BGN) on frequency bin 2 [4 - 12kBGNz]
CN.data.regional.modelfit$BGNfbin2.adj <- (CN.data.regional.modelfit$BGNfbin2)*(-1)
#hist(CN.data.regional.modelfit$BGNfbin2.adj); range(CN.data.regional.modelfit$BGNfbin2.adj)
#plot(family.test <- ecdf(CN.data.regional.modelfit$BGNfbin2.adj))
#plot(fitdist(CN.data.regional.modelfit$BGNfbin2.adj, "pois"))


#modBGNfbin2.null <- glmmTMB(BGNfbin2 ~ 1, family=beta_family(), CN.data.regional.modelfit)

#options(na.action = "na.fail")
modBGNfbin2.full <- glmmTMB(BGNfbin2.adj ~ z_alt + z_distwater + z_distedge + z_distcanga + z_distminning + z_treeheight +
                              z_NDVI + z_precipitation + z_temperature + z_humidity + z_wind + (1|Local:time_period),
                            family=genpois(link = "log"), CN.data.regional.modelfit)

check_collinearity(modBGNfbin2.full)


modBGNfbin2.full.b <- glmmTMB(BGNfbin2.adj ~ z_alt + z_distwater + z_distedge + z_treeheight + z_NDVI + z_precipitation + 
                                z_temperature + z_humidity + z_wind + (1|Local:time_period) + ar1(times - 1|Year),
                              family=genpois(link = "log"), CN.data.regional.modelfit)

summary(modBGNfbin2.full.b)
#r.squaredGLMM(modBGNfbin2.full.b)

test.modBGNfbin2.full.b <- dredge(modBGNfbin2.full.b)
#subset(test.modBGNfbin2.full.b, delta < 4)
summary(get.models(test.modBGNfbin2.full.b, 1)[[1]])
#r.squaredGLMM(get.models(test.modBGNfbin2.full.b, 1)[[1]])

sink("results/supplementar/BGNfbin2_model_summary.csv")
print(summary(get.models(test.modBGNfbin2.full.b, 1)[[1]]))
sink()

modBGNfbin2.sel.table<-as.data.frame(test.modBGNfbin2.full.b)[1:5,12:16]
modBGNfbin2.sel.table[,2:3]<- round(modBGNfbin2.sel.table[,2:3],2)
modBGNfbin2.sel.table[,4:5]<- round(modBGNfbin2.sel.table[,4:5],3)
names(modBGNfbin2.sel.table)[1] = "K"
modBGNfbin2.sel.table$Model<-NA
for(i in 1:nrow(modBGNfbin2.sel.table)) modBGNfbin2.sel.table$Model[i]<-as.character(formula(get.models(test.modBGNfbin2.full.b, i)[[1]]))[3]
modBGNfbin2.sel.table<-modBGNfbin2.sel.table[,c(6,1,2,3,4,5)]


modBGNfbin2.full.b.simres <- simulateResiduals(get.models(test.modBGNfbin2.full.b, 1)[[1]])
plot(modBGNfbin2.full.b.simres)
testDispersion(modBGNfbin2.full.b.simres)

modBGNfbin2.full.b.alleffects <- allEffects(get.models(test.modBGNfbin2.full.b, 1)[[1]], rescale.axis=F) #, residuals = T
plot(modBGNfbin2.full.b.alleffects)


modBGNfbin2.full.b.t <- broom.mixed::tidy(get.models(test.modBGNfbin2.full.b, 1)[[1]], effects = "fixed", conf.int = TRUE)
modBGNfbin2.full.b.t <- transform(modBGNfbin2.full.b.t, term = sprintf("%s.%s", component, term))
modBGNfbin2.full.b.t$group <- "BGN [4-12kHz]"
modBGNfbin2.full.b.t$col <- ifelse(modBGNfbin2.full.b.t$conf.low < 0 & modBGNfbin2.full.b.t$conf.high < 0, "#e65a49", 
                                   ifelse(modBGNfbin2.full.b.t$conf.low > 0 & modBGNfbin2.full.b.t$conf.high > 0, "#55c0f0", "#dbdfde"))



#Response: Acoustic entropy (S2N) on total frequency [0.3 - 22.1kS2Nz]
#hist(CN.data.regional.modelfit$S2NTotal); range(CN.data.regional.modelfit$S2NTotal)
#plot(family.test <- ecdf(CN.data.regional.modelfit$S2NTotal))
#plot(fitdist(CN.data.regional.modelfit$S2NTotal, "pois"))


#modS2N.null <- glmmTMB(S2NTotal ~ 1, family=genpois(link = "log"), CN.data.modelfit)


modS2N.full <- glmmTMB(S2NTotal ~ z_alt + z_distwater + z_distedge + z_distcanga + z_distminning + z_treeheight +
                         z_NDVI + z_precipitation + z_temperature + z_humidity + z_wind + (1|Local:time_period),
                       family=genpois(link = "log"), CN.data.regional.modelfit)

check_collinearity(modS2N.full)

modS2N.full.b <- glmmTMB(S2NTotal ~ z_alt + z_distwater + z_distedge + z_treeheight + z_NDVI + z_precipitation + 
                           z_temperature + z_humidity + z_wind + (1|Local:time_period) + ar1(times - 1|Year),
                         family=genpois(link = "log"), CN.data.regional.modelfit)

summary(modS2N.full.b)
#r.squaredGLMM(modS2N.full.b)

test.modS2N.full.b <- dredge(modS2N.full.b)
#subset(test.modS2N.full.b, delta < 4)
summary(get.models(test.modS2N.full.b, 1)[[1]])
#r.squaredGLMM(get.models(test.modS2N.full.b, 1)[[1]])

sink("results/supplementar/S2NTotal_model_summary.csv")
print(summary(get.models(test.modS2N.full.b, 1)[[1]]))
sink()

modS2N.sel.table<-as.data.frame(test.modS2N.full.b)[1:5,12:16]
modS2N.sel.table[,2:3]<- round(modS2N.sel.table[,2:3],2)
modS2N.sel.table[,4:5]<- round(modS2N.sel.table[,4:5],3)
names(modS2N.sel.table)[1] = "K"
modS2N.sel.table$Model<-NA
for(i in 1:nrow(modS2N.sel.table)) modS2N.sel.table$Model[i]<-as.character(formula(get.models(test.modS2N.full.b, i)[[1]]))[3]
modS2N.sel.table<-modS2N.sel.table[,c(6,1,2,3,4,5)]


modS2N.full.b.simres <- simulateResiduals(get.models(test.modS2N.full.b, 1)[[1]])
plot(modS2N.full.b.simres)
testDispersion(modS2N.full.b.simres)

modS2N.full.b.alleffects <- allEffects(get.models(test.modS2N.full.b, 1)[[1]], rescale.axis=F) #, residuals = T
plot(modS2N.full.b.alleffects)


modS2N.full.b.t <- broom.mixed::tidy(get.models(test.modS2N.full.b, 1)[[1]], effects = "fixed", conf.int = TRUE)
modS2N.full.b.t <- transform(modS2N.full.b.t, term = sprintf("%s.%s", component, term))
modS2N.full.b.t$group <- "S2N [0.3-22.1kHz]"
modS2N.full.b.t$col <- ifelse(modS2N.full.b.t$conf.low < 0 & modS2N.full.b.t$conf.high < 0, "#e65a49", 
                              ifelse(modS2N.full.b.t$conf.low > 0 & modS2N.full.b.t$conf.high > 0, "#55c0f0", "#dbdfde"))



#Response: Acoustic entropy (S2N) on frequency bin 2 [4 - 12kS2Nz]

#modS2Nfbin2.null <- glmmTMB(S2Nfbin2 ~ 1, family=beta_family(), CN.data.regional.modelfit)

#options(na.action = "na.fail")
modS2Nfbin2.full <- glmmTMB(S2Nfbin2 ~ z_alt + z_distwater + z_distedge + z_distcanga + z_distminning + z_treeheight +
                              z_NDVI + z_precipitation + z_temperature + z_humidity + z_wind + (1|Local:time_period),
                            family=genpois(link = "log"), CN.data.regional.modelfit)

check_collinearity(modS2Nfbin2.full)


modS2Nfbin2.full.b <- glmmTMB(S2Nfbin2 ~ z_alt + z_distwater + z_distedge + z_treeheight + z_NDVI + z_precipitation + 
                                z_temperature + z_humidity + z_wind + (1|Local:time_period) + ar1(times - 1|Year),
                              family=genpois(link = "log"), CN.data.regional.modelfit)

summary(modS2Nfbin2.full.b)
#r.squaredGLMM(modS2Nfbin2.full.b)

test.modS2Nfbin2.full.b <- dredge(modS2Nfbin2.full.b)
#subset(test.modS2Nfbin2.full.b, delta < 4)
summary(get.models(test.modS2Nfbin2.full.b, 1)[[1]])
#r.squaredGLMM(get.models(test.modS2Nfbin2.full.b, 1)[[1]])

sink("results/supplementar/S2Nfbin2_model_summary.csv")
print(summary(get.models(test.modS2Nfbin2.full.b, 1)[[1]]))
sink()

modS2Nfbin2.sel.table<-as.data.frame(test.modS2Nfbin2.full.b)[1:5,12:16]
modS2Nfbin2.sel.table[,2:3]<- round(modS2Nfbin2.sel.table[,2:3],2)
modS2Nfbin2.sel.table[,4:5]<- round(modS2Nfbin2.sel.table[,4:5],3)
names(modS2Nfbin2.sel.table)[1] = "K"
modS2Nfbin2.sel.table$Model<-NA
for(i in 1:nrow(modS2Nfbin2.sel.table)) modS2Nfbin2.sel.table$Model[i]<-as.character(formula(get.models(test.modS2Nfbin2.full.b, i)[[1]]))[3]
modS2Nfbin2.sel.table<-modS2Nfbin2.sel.table[,c(6,1,2,3,4,5)]


modS2Nfbin2.full.b.simres <- simulateResiduals(get.models(test.modS2Nfbin2.full.b, 1)[[1]])
plot(modS2Nfbin2.full.b.simres)
testDispersion(modS2Nfbin2.full.b.simres)

modS2Nfbin2.full.b.alleffects <- allEffects(get.models(test.modS2Nfbin2.full.b, 1)[[1]], rescale.axis=F) #, residuals = T
plot(modS2Nfbin2.full.b.alleffects)


modS2Nfbin2.full.b.t <- broom.mixed::tidy(get.models(test.modS2Nfbin2.full.b, 1)[[1]], effects = "fixed", conf.int = TRUE)
modS2Nfbin2.full.b.t <- transform(modS2Nfbin2.full.b.t, term = sprintf("%s.%s", component, term))
modS2Nfbin2.full.b.t$group <- "S2N [4-12kHz]"
modS2Nfbin2.full.b.t$col <- ifelse(modS2Nfbin2.full.b.t$conf.low < 0 & modS2Nfbin2.full.b.t$conf.high < 0, "#e65a49", 
                                   ifelse(modS2Nfbin2.full.b.t$conf.low > 0 & modS2Nfbin2.full.b.t$conf.high > 0, "#55c0f0", "#dbdfde"))



#Response: Acoustic entropy (AA) on total frequency [0.3 - 22.1kAAz]
#hist(CN.data.regional.modelfit$AATotal); range(CN.data.regional.modelfit$AATotal)
#plot(family.test <- ecdf(CN.data.regional.modelfit$AATotal))
#plot(fitdist(CN.data.regional.modelfit$AATotal, "beta"))


#modAA.null <- glmmTMB(AATotal ~ 1, family=beta_family(), CN.data.regional.modelfit)

#options(na.action = "na.fail")
modAA.full <- glmmTMB(AATotal ~ z_alt + z_distwater + z_distedge + z_distcanga + z_distminning + z_treeheight +
                         z_NDVI + z_precipitation + z_temperature + z_humidity + z_wind + (1|Local:time_period),
                       family=beta_family(), CN.data.regional.modelfit)

check_collinearity(modAA.full)


modAA.full.b <- glmmTMB(AATotal ~ z_alt + z_distwater + z_distedge + z_treeheight + z_NDVI + z_precipitation + 
                           z_temperature + z_humidity + z_wind + (1|Local:time_period) + ar1(times - 1|Year),
                         family=beta_family(), CN.data.regional.modelfit)

summary(modAA.full.b)
#r.squaredGLMM(modAA.full.b)

test.modAA.full.b <- dredge(modAA.full.b)
#subset(test.modAA.full.b, delta < 4)
summary(get.models(test.modAA.full.b, 1)[[1]])
#r.squaredGLMM(get.models(test.modAA.full.b, 1)[[1]])

sink("results/supplementar/AATotal_model_summary.csv")
print(summary(get.models(test.modAA.full.b, 1)[[1]]))
sink()

modAA.sel.table<-as.data.frame(test.modAA.full.b)[1:5,12:16]
modAA.sel.table[,2:3]<- round(modAA.sel.table[,2:3],2)
modAA.sel.table[,4:5]<- round(modAA.sel.table[,4:5],3)
names(modAA.sel.table)[1] = "K"
modAA.sel.table$Model<-NA
for(i in 1:nrow(modAA.sel.table)) modAA.sel.table$Model[i]<-as.character(formula(get.models(test.modAA.full.b, i)[[1]]))[3]
modAA.sel.table<-modAA.sel.table[,c(6,1,2,3,4,5)]


modAA.full.b.simres <- simulateResiduals(get.models(test.modAA.full.b, 1)[[1]])
plot(modAA.full.b.simres)
testDispersion(modAA.full.b.simres)

modAA.full.b.alleffects <- allEffects(get.models(test.modAA.full.b, 1)[[1]], rescale.axis=F) #, residuals = T
plot(modAA.full.b.alleffects)


modAA.full.b.t <- broom.mixed::tidy(get.models(test.modAA.full.b, 1)[[1]], effects = "fixed", conf.int = TRUE)
modAA.full.b.t <- transform(modAA.full.b.t, term = sprintf("%s.%s", component, term))
modAA.full.b.t$group <- "AA [0.3-22.1kHz]"
modAA.full.b.t$col <- ifelse(modAA.full.b.t$conf.low < 0 & modAA.full.b.t$conf.high < 0, "#e65a49", 
                              ifelse(modAA.full.b.t$conf.low > 0 & modAA.full.b.t$conf.high > 0, "#55c0f0", "#dbdfde"))



#Response: Acoustic entropy (AA) on frequency bin 2 [4 - 12kAAz]
#hist(CN.data.regional.modelfit$AAfbin2); range(CN.data.regional.modelfit$AAfbin2)
#plot(family.test <- ecdf(CN.data.regional.modelfit$AAfbin2))
#plot(fitdist(CN.data.regional.modelfit$AAfbin2, "beta"))

CN.data.regional.modelfit$AAfbin2b <- rescale(CN.data.regional.modelfit$AAfbin2, to=c(.001,.999))

#modAAfbin2.null <- glmmTMB(AAfbin2 ~ 1, family=beta_family(), CN.data.regional.modelfit)

#options(na.action = "na.fail")
modAAfbin2.full <- glmmTMB(AAfbin2b ~ z_alt + z_distwater + z_distedge + z_distcanga + z_distminning + z_treeheight +
                              z_NDVI + z_precipitation + z_temperature + z_humidity + z_wind + (1|Local:time_period),
                            family=beta_family(), CN.data.regional.modelfit)

check_collinearity(modAAfbin2.full)


modAAfbin2.full.b <- glmmTMB(AAfbin2b ~ z_alt + z_distwater + z_distedge + z_treeheight + z_NDVI + z_precipitation + 
                                z_temperature + z_humidity + z_wind + (1|Local:time_period) + ar1(times - 1|Year),
                              family=beta_family(), CN.data.regional.modelfit)

summary(modAAfbin2.full.b)
#r.squaredGLMM(modAAfbin2.full.b)

test.modAAfbin2.full.b <- dredge(modAAfbin2.full.b)
#subset(test.modAAfbin2.full.b, delta < 4)
summary(get.models(test.modAAfbin2.full.b, 1)[[1]])
#r.squaredGLMM(get.models(test.modAAfbin2.full.b, 1)[[1]])

sink("results/supplementar/AAfbin2_model_summary.csv")
print(summary(get.models(test.modAAfbin2.full.b, 1)[[1]]))
sink()

modAAfbin2.sel.table<-as.data.frame(test.modAAfbin2.full.b)[1:5,12:16]
modAAfbin2.sel.table[,2:3]<- round(modAAfbin2.sel.table[,2:3],2)
modAAfbin2.sel.table[,4:5]<- round(modAAfbin2.sel.table[,4:5],3)
names(modAAfbin2.sel.table)[1] = "K"
modAAfbin2.sel.table$Model<-NA
for(i in 1:nrow(modAAfbin2.sel.table)) modAAfbin2.sel.table$Model[i]<-as.character(formula(get.models(test.modAAfbin2.full.b, i)[[1]]))[3]
modAAfbin2.sel.table<-modAAfbin2.sel.table[,c(6,1,2,3,4,5)]


modAAfbin2.full.b.simres <- simulateResiduals(get.models(test.modAAfbin2.full.b, 1)[[1]])
plot(modAAfbin2.full.b.simres)
testDispersion(modAAfbin2.full.b.simres)

modAAfbin2.full.b.alleffects <- allEffects(get.models(test.modAAfbin2.full.b, 1)[[1]], rescale.axis=F) #, residuals = T
plot(modAAfbin2.full.b.alleffects)


modAAfbin2.full.b.t <- broom.mixed::tidy(get.models(test.modAAfbin2.full.b, 1)[[1]], effects = "fixed", conf.int = TRUE)
modAAfbin2.full.b.t <- transform(modAAfbin2.full.b.t, term = sprintf("%s.%s", component, term))
modAAfbin2.full.b.t$group <- "AA [4-12kHz]"
modAAfbin2.full.b.t$col <- ifelse(modAAfbin2.full.b.t$conf.low < 0 & modAAfbin2.full.b.t$conf.high < 0, "#e65a49", 
                                   ifelse(modAAfbin2.full.b.t$conf.low > 0 & modAAfbin2.full.b.t$conf.high > 0, "#55c0f0", "#dbdfde"))



options(na.action = "na.omit")

sel.table <- rbind(modH.sel.table, modHfbin2.sel.table,
               modADI.sel.table, modADIfbin2.sel.table,
               modAEI.sel.table, modAEIfbin2.sel.table,
               modBGN.sel.table, modBGNfbin2.sel.table,
               modS2N.sel.table, modS2Nfbin2.sel.table,
               modAA.sel.table, modAAfbin2.sel.table)
write.csv(sel.table,"results/supplementar/selection_table.csv", row.names = F)



mod.t <- rbind(modH.full.b.t, modHfbin2.full.b.t,
               modADI.full.b.t, modADIfbin2.full.b.t,
               modAEI.full.b.t, modAEIfbin2.full.b.t,
               modBGN.full.b.t, modBGNfbin2.full.b.t,
               modS2N.full.b.t, modS2Nfbin2.full.b.t,
               modAA.full.b.t, modAAfbin2.full.b.t)

#mod.t <- mod.t %>% filter(term!="cond.(Intercept)" & 
#                          term!="cond.sd__(Intercept)" &
#                          component!="zi") %>% 
#  mutate(group=factor(group, levels = c("H [0.3-22.1kHz]", "H [4-12kHz]", 
#                                        "ADI [0.3-22.1kHz]", "ADI [4-12kHz]",
#                                        "AEI [0.3-22.1kHz]", "AEI [4-12kHz]",
#                                        "BGN [0.3-22.1kHz]", "BGN [4-12kHz]",
#                                        "S2N [0.3-22.1kHz]", "S2N [4-12kHz]",
#                                        "AA [0.3-22.1kHz]", "AA [4-12kHz]")),
#         term=factor(term, levels = c("cond.z_alt", 
#                                      "cond.z_distwater", 
#                                      "cond.z_distedge", 
#                                      "cond.z_treeheight",
#                                      "cond.z_NDVI", 
#                                      "cond.z_precipitation", 
#                                      "cond.z_temperature",
#                                      "cond.z_humidity", 
#                                      "cond.z_wind")))
#
##write.csv(mod.t, "results/effect_sizes.csv", row.names = F)


#mod.t.graph <- mod.t %>% filter(group %in% unique(grep("ADI|S2N", mod.t$group, value = T, invert = T))) %>% 
#  ggplot()+
#  geom_point(aes(x = estimate, y = term, col = col), size = 12)+
#  geom_segment(aes(x = conf.low, xend = conf.high,
#                   y = term, yend = term, col = col), size = 7)+
#  geom_vline(xintercept=0, lty=2)+
#  scale_color_manual(values = c("#55c0f0", "#dbdfde", "#e65a49"))+
#  scale_y_discrete(labels = c("cond.z_alt"="Altitude",
#                              "cond.z_distwater"="Distance from water",
#                              "cond.z_distedge"="Distance from edge",
#                              "cond.z_treeheight"="Tree height",
#                              "cond.z_NDVI"="NDVI",
#                              "cond.z_precipitation"="Precipitation",
#                              "cond.z_temperature"="Temperature",
#                              "cond.z_humidity" = "Air humidity",
#                              "cond.z_wind"="Wind speed"))+
#  facet_wrap(~group, ncol = 2)+#, labeller = labeller(group = c("H" = "H [0.3 - 22.1kHz]",
#                               #                             "Hfbin2" = "H [4 - 12kHz]",
#                               #                             "AEI" = "AEI [0.3 - 22.1kHz]",
#                               #                             "AEIfbin2" = "AEI [4 - 12kHz]",
#                               #                             "BGN" = "BGN [0.3 - 22.1kHz]",
#                               #                             "BGNb2" = "BGN [4 - 12kHz]",
#                               #                             "AA" = "AA [0.3 - 22.1kHz]",
#                               #                             "AAfbin2" = "AA [4 - 12kHz]")))+
#  xlab("")+ ylab("")+
#  theme(plot.title = element_text(family = "serif", size = 48, hjust = .5),
#        strip.text = element_text(family = "serif", size = 36, hjust = .5),
#        axis.title = element_text(family = "serif", size = 36),
#        axis.text.x = element_text(family = "serif", size = 32),
#        axis.text.y = element_text(family = "serif", size = 32),
#        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
#        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2),
#        panel.background = element_blank(), legend.position = "none")
#
#
#
#
#png("results/regionalenvironmentaleffects.png", width = 23.38583, height = 33.11024, units = "in", res = 300, pointsize = 4)
#mod.t.graph
#dev.off()






## preparing data -- local
CN.data.local.modelfit <- CN.data.regional.modelfit[!is.na(CN.data.regional.modelfit$n_trees),]


#Response: Acoustic entropy (H) on total frequency [0.3 - 22.1kHz]
#hist(CN.data.local.modelfit$HTotal); range(CN.data.local.modelfit$HTotal)
#plot(family.test <- ecdf(CN.data.local.modelfit$HTotal))
#plot(fitdist(CN.data.local.modelfit$HTotal, "beta"))


#modH.null <- glmmTMB(HTotal ~ 1, family=beta_family(), CN.data.local.modelfit)

options(na.action = "na.fail")
local.modH.full <- glmmTMB(HTotal ~ n_trees + n_vines + n_palmtrees + (1|Local:time_period),
                     family=beta_family(), CN.data.local.modelfit)

check_collinearity(local.modH.full)


local.modH.full.b <- glmmTMB(HTotal ~ n_trees + n_vines + n_palmtrees + (1|Local:time_period) + ar1(times - 1|Year),
                       family=beta_family(), CN.data.local.modelfit)

summary(local.modH.full.b)
#r.squaredGLMM(modH.full.b)

test.local.modH.full.b <- dredge(modH.full.b)
#subset(test.modH.full.b, delta < 2)
summary(get.models(test.local.modH.full.b, 1)[[1]])
#r.squaredGLMM(get.models(test.modH.full.b, 1)[[1]])

sink("results/supplementar/HTotal_local_model_summary.csv")
print(summary(get.models(test.local.modH.full.b, 1)[[1]]))
sink()

local.modH.sel.table<-as.data.frame(test.local.modH.full.b)[1:5,6:10]
local.modH.sel.table[,2:3]<- round(local.modH.sel.table[,2:3],2)
local.modH.sel.table[,4:5]<- round(local.modH.sel.table[,4:5],3)
names(local.modH.sel.table)[1] = "K"
local.modH.sel.table$Model<-NA
for(i in 1:nrow(local.modH.sel.table)) local.modH.sel.table$Model[i]<-as.character(formula(get.models(test.local.modH.full.b, i)[[1]]))[3]
local.modH.sel.table<-local.modH.sel.table[,c(6,1,2,3,4,5)]


local.modH.full.b.simres <- simulateResiduals(get.models(test.local.modH.full.b, 1)[[1]])
plot(local.modH.full.b.simres)
testDispersion(local.modH.full.b.simres)

local.modH.full.b.alleffects <- allEffects(get.models(test.local.modH.full.b, 1)[[1]], rescale.axis=F) #, residuals = T
plot(local.modH.full.b.alleffects)


local.modH.full.b.t <- broom.mixed::tidy(get.models(test.local.modH.full.b, 1)[[1]], effects = "fixed", conf.int = TRUE)
local.modH.full.b.t <- transform(local.modH.full.b.t, term = sprintf("%s.%s", component, term))
local.modH.full.b.t$group <- "H [0.3-22.1kHz]"
local.modH.full.b.t$col <- ifelse(local.modH.full.b.t$conf.low < 0 & local.modH.full.b.t$conf.high < 0, "#e65a49", 
                            ifelse(local.modH.full.b.t$conf.low > 0 & local.modH.full.b.t$conf.high > 0, "#55c0f0", "#dbdfde"))



#Response: Acoustic entropy (H) on frequency bin 2 [4 - 12kHz]
#hist(CN.data.local.modelfit$Hfbin2); range(CN.data.local.modelfit$Hfbin2)
#plot(family.test <- ecdf(CN.data.local.modelfit$Hfbin2))
#plot(fitdist(CN.data.local.modelfit$Hfbin2, "beta"))


#modHfbin2.null <- glmmTMB(Hfbin2 ~ 1, family=beta_family(), CN.data.local.modelfit)

local.modHfbin2.full <- glmmTMB(Hfbin2 ~ n_trees + n_vines + n_palmtrees + (1|Local:time_period),
                          family=beta_family(), CN.data.local.modelfit)

check_collinearity(local.modHfbin2.full)


local.modHfbin2.full.b <- glmmTMB(Hfbin2 ~ n_trees + n_vines + n_palmtrees + (1|Local:time_period) + ar1(times - 1|Year),
                            family=beta_family(), CN.data.local.modelfit)

summary(local.modHfbin2.full.b)
#r.squaredGLMM(modHfbin2.full.b)

test.local.modHfbin2.full.b <- dredge(local.modHfbin2.full.b)
#subset(test.modHfbin2.full.b, delta < 2)
summary(get.models(test.local.modHfbin2.full.b, 1)[[1]])
#r.squaredGLMM(get.models(test.modHfbin2.full.b, 1)[[1]])

sink("results/supplementar/Hfbin2_local_model_summary.csv")
print(summary(get.models(test.local.modHfbin2.full.b, 1)[[1]]))
sink()

local.modHfbin2.sel.table<-as.data.frame(test.local.modHfbin2.full.b)[1:5,6:10]
local.modHfbin2.sel.table[,2:3]<- round(local.modHfbin2.sel.table[,2:3],2)
local.modHfbin2.sel.table[,4:5]<- round(local.modHfbin2.sel.table[,4:5],3)
names(local.modHfbin2.sel.table)[1] = "K"
local.modHfbin2.sel.table$Model<-NA
for(i in 1:nrow(local.modHfbin2.sel.table)) local.modHfbin2.sel.table$Model[i]<-as.character(formula(get.models(test.local.modHfbin2.full.b, i)[[1]]))[3]
local.modHfbin2.sel.table<-local.modHfbin2.sel.table[,c(6,1,2,3,4,5)]


local.modHfbin2.full.b.simres <- simulateResiduals(get.models(test.local.modHfbin2.full.b, 1)[[1]])
plot(local.modHfbin2.full.b.simres)
testDispersion(local.modHfbin2.full.b.simres)

local.modHfbin2.full.b.alleffects <- allEffects(get.models(test.local.modHfbin2.full.b, 1)[[1]], rescale.axis=F) #, residuals = T
plot(local.modHfbin2.full.b.alleffects)


local.modHfbin2.full.b.t <- broom.mixed::tidy(get.models(test.local.modHfbin2.full.b, 1)[[1]], effects = "fixed", conf.int = TRUE)
local.modHfbin2.full.b.t <- transform(local.modHfbin2.full.b.t, term = sprintf("%s.%s", component, term))
local.modHfbin2.full.b.t$group <- "H [4-12kHz]"
local.modHfbin2.full.b.t$col <- ifelse(local.modHfbin2.full.b.t$conf.low < 0 & local.modHfbin2.full.b.t$conf.high < 0, "#e65a49", 
                                 ifelse(local.modHfbin2.full.b.t$conf.low > 0 & local.modHfbin2.full.b.t$conf.high > 0, "#55c0f0", "#dbdfde"))



##Response: Acoustic diversity index (ADI) on total frequency [0.3 - 22.1kHz]
#CN.data.local.modelfit$ADITotal.adj <- 1-log(CN.data.local.modelfit$ADITotal)
##hist(CN.data.local.modelfit$ADITotal.adj); range(CN.data.local.modelfit$ADITotal.adj)
##plot(family.test <- ecdf(CN.data.local.modelfit$ADITotal.adj))
##plot(fitdist(CN.data.local.modelfit$ADITotal.adj, "lnorm"))
#
#
##modADI.null <- glmmTMB(ADITotal ~ 1, CN.data.local.modelfit)
#
#local.modADI.full <- glmmTMB(ADITotal.adj ~ n_trees + n_vines + n_palmtrees + (1|Local:time_period),
#                       #family=beta_family(), 
#                       CN.data.local.modelfit)
#
#check_collinearity(local.modADI.full)
#
#
#local.modADI.full.b <- glmmTMB(ADITotal.adj ~ n_trees + n_vines + n_palmtrees + (1|Local:time_period) + ar1(times - 1|Year), ziformula = ~1,
#                         #family=beta_family(), 
#                         CN.data.local.modelfit)
#
##fixef(modADI.full.b)
#local.modADI.full.c <- update(local.modADI.full.b, control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
#
#summary(local.modADI.full.b)
##r.squaredGLMM(modADI.full.b)
#
#local.modADI.full.c <- update(local.modADI.full.b, control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
##Warning message: Model convergence problem; extreme or very small eigenvalues detected.



#Response: Acoustic entropy (AEI) on total frequency [0.3 - 22.1kAEIz]
CN.data.local.modelfit$AEITotal.adj <- 1-log(CN.data.local.modelfit$AEITotal)
#hist(CN.data.local.modelfit$AEITotal.adj); range(CN.data.local.modelfit$AEITotal.adj)
#plot(family.test <- ecdf(CN.data.local.modelfit$AEITotal.adj))
#plot(fitdist(CN.data.local.modelfit$AEITotal.adj, "norm"))


#modAEI.null <- glmmTMB(AEITotal.adj ~ 1, CN.data.local.modelfit)

local.modAEI.full <- glmmTMB(AEITotal.adj ~ n_trees + n_vines + n_palmtrees + (1|Local:time_period),
                       #family=beta_family(), 
                       CN.data.local.modelfit)

check_collinearity(local.modAEI.full)


local.modAEI.full.b <- glmmTMB(AEITotal.adj ~ n_trees + n_vines + n_palmtrees + (1|Local:time_period) + ar1(times - 1|Year),
                         #family=beta_family(), 
                         CN.data.local.modelfit)

summary(local.modAEI.full.b)
#r.squaredGLMM(modAEI.full.b)

test.local.modAEI.full.b <- dredge(local.modAEI.full.b)
#subset(test.modAEI.full.b, delta < 2)
summary(get.models(test.local.modAEI.full.b, 1)[[1]])
#r.squaredGLMM(get.models(test.modAEI.full.b, 1)[[1]])

sink("results/supplementar/AEITotal_local_model_summary.csv")
print(summary(get.models(test.local.modAEI.full.b, 1)[[1]]))
sink()

local.modAEI.sel.table<-as.data.frame(test.local.modAEI.full.b)[1:5,6:10]
local.modAEI.sel.table[,2:3]<- round(local.modAEI.sel.table[,2:3],2)
local.modAEI.sel.table[,4:5]<- round(local.modAEI.sel.table[,4:5],3)
names(local.modAEI.sel.table)[1] = "K"
local.modAEI.sel.table$Model<-NA
for(i in 1:nrow(local.modAEI.sel.table)) local.modAEI.sel.table$Model[i]<-as.character(formula(get.models(test.local.modAEI.full.b, i)[[1]]))[3]
local.modAEI.sel.table<-local.modAEI.sel.table[,c(6,1,2,3,4,5)]


local.modAEI.full.b.simres <- simulateResiduals(get.models(test.local.modAEI.full.b, 1)[[1]])
plot(local.modAEI.full.b.simres)
testDispersion(local.modAEI.full.b.simres)

local.modAEI.full.b.alleffects <- allEffects(get.models(test.local.modAEI.full.b, 1)[[1]], rescale.axis=F) #, residuals = T
plot(local.modAEI.full.b.alleffects)


local.modAEI.full.b.t <- broom.mixed::tidy(get.models(test.local.modAEI.full.b, 1)[[1]], effects = "fixed", conf.int = TRUE)
local.modAEI.full.b.t <- transform(local.modAEI.full.b.t, term = sprintf("%s.%s", component, term))
local.modAEI.full.b.t$group <- "AEI [0.3-22.1kHz]"
local.modAEI.full.b.t$col <- ifelse(local.modAEI.full.b.t$conf.low < 0 & local.modAEI.full.b.t$conf.high < 0, "#e65a49", 
                              ifelse(local.modAEI.full.b.t$conf.low > 0 & local.modAEI.full.b.t$conf.high > 0, "#55c0f0", "#dbdfde"))



#Response: Acoustic entropy (AEI) on frequency bin 2 [4 - 12kAEIz]
CN.data.local.modelfit$AEIfbin2.adj <- 1-log(CN.data.local.modelfit$AEIfbin2)
#hist(CN.data.local.modelfit$AEIfbin2.adj); range(CN.data.local.modelfit$AEIfbin2.adj)
#plot(family.test <- ecdf(CN.data.local.modelfit$AEIfbin2.adj))
#plot(fitdist(CN.data.local.modelfit$AEIfbin2.adj, "lnorm"))


#modAEIfbin2.null <- glmmTMB(AEIfbin2 ~ 1, family=beta_family(), CN.data.local.modelfit)

local.modAEIfbin2.full <- glmmTMB(AEIfbin2.adj ~ n_trees + n_vines + n_palmtrees + (1|Local:time_period),
                            #family=beta_family(), 
                            CN.data.local.modelfit)

check_collinearity(local.modAEIfbin2.full)


local.modAEIfbin2.full.b <- glmmTMB(AEIfbin2.adj ~ n_trees + n_vines + n_palmtrees + (1|Local:time_period) + ar1(times - 1|Year),
                              #family=beta_family(), 
                              CN.data.local.modelfit)

summary(local.modAEIfbin2.full.b)
#r.squaredGLMM(modAEIfbin2.full.b)

test.local.modAEIfbin2.full.b <- dredge(local.modAEIfbin2.full.b)
#subset(test.modAEIfbin2.full.b, delta < 4)
summary(get.models(test.local.modAEIfbin2.full.b, 1)[[1]])
#r.squaredGLMM(get.models(test.modAEIfbin2.full.b, 1)[[1]])

sink("results/supplementar/AEIfbin2_local_model_summary.csv")
print(summary(get.models(test.local.modAEIfbin2.full.b, 1)[[1]]))
sink()

local.modAEIfbin2.sel.table<-as.data.frame(test.local.modAEIfbin2.full.b)[1:5,6:10]
local.modAEIfbin2.sel.table[,2:3]<- round(local.modAEIfbin2.sel.table[,2:3],2)
local.modAEIfbin2.sel.table[,4:5]<- round(local.modAEIfbin2.sel.table[,4:5],3)
names(local.modAEIfbin2.sel.table)[1] = "K"
local.modAEIfbin2.sel.table$Model<-NA
for(i in 1:nrow(local.modAEIfbin2.sel.table)) local.modAEIfbin2.sel.table$Model[i]<-as.character(formula(get.models(test.local.modAEIfbin2.full.b, i)[[1]]))[3]
local.modAEIfbin2.sel.table<-local.modAEIfbin2.sel.table[,c(6,1,2,3,4,5)]


local.modAEIfbin2.full.b.simres <- simulateResiduals(get.models(test.local.modAEIfbin2.full.b, 1)[[1]])
plot(local.modAEIfbin2.full.b.simres)
testDispersion(local.modAEIfbin2.full.b.simres)

local.modAEIfbin2.full.b.alleffects <- allEffects(get.models(test.local.modAEIfbin2.full.b, 1)[[1]], rescale.axis=F) #, residuals = T
plot(local.modAEIfbin2.full.b.alleffects)


local.modAEIfbin2.full.b.t <- broom.mixed::tidy(get.models(test.local.modAEIfbin2.full.b, 1)[[1]], effects = "fixed", conf.int = TRUE)
local.modAEIfbin2.full.b.t <- transform(local.modAEIfbin2.full.b.t, term = sprintf("%s.%s", component, term))
local.modAEIfbin2.full.b.t$group <- "AEI [4-12kHz]"
local.modAEIfbin2.full.b.t$col <- ifelse(local.modAEIfbin2.full.b.t$conf.low < 0 & local.modAEIfbin2.full.b.t$conf.high < 0, "#e65a49", 
                                   ifelse(local.modAEIfbin2.full.b.t$conf.low > 0 & local.modAEIfbin2.full.b.t$conf.high > 0, "#55c0f0", "#dbdfde"))



#Response: Acoustic entropy (BGN) on total frequency [0.3 - 22.1kBGNz]
CN.data.local.modelfit$BGNTotal.adj <- (CN.data.local.modelfit$BGNTotal)*(-1)
#hist(CN.data.local.modelfit$BGNTotal.adj); range(CN.data.local.modelfit$BGNTotal.adj)
#plot(family.test <- ecdf(CN.data.local.modelfit$BGNTotal.adj))
#plot(fitdist(CN.data.local.modelfit$BGNTotal.adj, "pois"))


#modBGN.null <- glmmTMB(BGNTotal.adj ~ 1, family=genpois(link = "log"), CN.data.modelfit)


local.modBGN.full <- glmmTMB(BGNTotal.adj ~ n_trees + n_vines + n_palmtrees + (1|Local:time_period),
                       family=genpois(link = "log"), CN.data.local.modelfit)

check_collinearity(local.modBGN.full)

local.modBGN.full.b <- glmmTMB(BGNTotal.adj ~ n_trees + n_vines + n_palmtrees + (1|Local:time_period) + ar1(times - 1|Year),
                         family=genpois(link = "log"), CN.data.local.modelfit)

summary(local.modBGN.full.b)
#r.squaredGLMM(modBGN.full.b)

test.local.modBGN.full.b <- dredge(local.modBGN.full.b)
#subset(test.modBGN.full.b, delta < 4)
summary(get.models(test.local.modBGN.full.b, 1)[[1]])
#r.squaredGLMM(get.models(test.modBGN.full.b, 1)[[1]])

sink("results/supplementar/BGNTotal_local_model_summary.csv")
print(summary(get.models(test.local.modBGN.full.b, 1)[[1]]))
sink()

local.modBGN.sel.table<-as.data.frame(test.local.modBGN.full.b)[1:5,6:10]
local.modBGN.sel.table[,2:3]<- round(local.modBGN.sel.table[,2:3],2)
local.modBGN.sel.table[,4:5]<- round(local.modBGN.sel.table[,4:5],3)
names(local.modBGN.sel.table)[1] = "K"
local.modBGN.sel.table$Model<-NA
for(i in 1:nrow(local.modBGN.sel.table)) local.modBGN.sel.table$Model[i]<-as.character(formula(get.models(test.local.modBGN.full.b, i)[[1]]))[3]
local.modBGN.sel.table<-local.modBGN.sel.table[,c(6,1,2,3,4,5)]


local.modBGN.full.b.simres <- simulateResiduals(get.models(test.local.modBGN.full.b, 1)[[1]])
plot(local.modBGN.full.b.simres)
testDispersion(local.modBGN.full.b.simres)

local.modBGN.full.b.alleffects <- allEffects(get.models(test.local.modBGN.full.b, 1)[[1]], rescale.axis=F) #, residuals = T
plot(local.modBGN.full.b.alleffects)


local.modBGN.full.b.t <- broom.mixed::tidy(get.models(test.local.modBGN.full.b, 1)[[1]], effects = "fixed", conf.int = TRUE)
local.modBGN.full.b.t <- transform(local.modBGN.full.b.t, term = sprintf("%s.%s", component, term))
local.modBGN.full.b.t$group <- "BGN [0.3-22.1kHz]"
local.modBGN.full.b.t$col <- ifelse(local.modBGN.full.b.t$conf.low < 0 & local.modBGN.full.b.t$conf.high < 0, "#e65a49", 
                              ifelse(local.modBGN.full.b.t$conf.low > 0 & local.modBGN.full.b.t$conf.high > 0, "#55c0f0", "#dbdfde"))



#Response: Acoustic entropy (BGN) on frequency bin 2 [4 - 12kBGNz]
CN.data.local.modelfit$BGNfbin2.adj <- (CN.data.local.modelfit$BGNfbin2)*(-1)
#hist(CN.data.local.modelfit$BGNfbin2.adj); range(CN.data.local.modelfit$BGNfbin2.adj)
#plot(family.test <- ecdf(CN.data.local.modelfit$BGNfbin2.adj))
#plot(fitdist(CN.data.local.modelfit$BGNfbin2.adj, "pois"))


#modBGNfbin2.null <- glmmTMB(BGNfbin2 ~ 1, family=beta_family(), CN.data.local.modelfit)

#options(na.action = "na.fail")
local.modBGNfbin2.full <- glmmTMB(BGNfbin2.adj ~ n_trees + n_vines + n_palmtrees + (1|Local:time_period),
                            family=genpois(link = "log"), CN.data.local.modelfit)

check_collinearity(local.modBGNfbin2.full)


local.modBGNfbin2.full.b <- glmmTMB(BGNfbin2.adj ~ n_trees + n_vines + n_palmtrees + (1|Local:time_period) + ar1(times - 1|Year),
                              family=genpois(link = "log"), CN.data.local.modelfit)

summary(local.modBGNfbin2.full.b)
#r.squaredGLMM(modBGNfbin2.full.b)

test.local.modBGNfbin2.full.b <- dredge(local.modBGNfbin2.full.b)
#subset(test.modBGNfbin2.full.b, delta < 4)
summary(get.models(test.local.modBGNfbin2.full.b, 1)[[1]])
#r.squaredGLMM(get.models(test.modBGNfbin2.full.b, 1)[[1]])

sink("results/supplementar/BGNfbin2_local_model_summary.csv")
print(summary(get.models(test.local.modBGNfbin2.full.b, 1)[[1]]))
sink()

local.modBGNfbin2.sel.table<-as.data.frame(test.local.modBGNfbin2.full.b)[1:5,6:10]
local.modBGNfbin2.sel.table[,2:3]<- round(local.modBGNfbin2.sel.table[,2:3],2)
local.modBGNfbin2.sel.table[,4:5]<- round(local.modBGNfbin2.sel.table[,4:5],3)
names(local.modBGNfbin2.sel.table)[1] = "K"
local.modBGNfbin2.sel.table$Model<-NA
for(i in 1:nrow(local.modBGNfbin2.sel.table)) local.modBGNfbin2.sel.table$Model[i]<-as.character(formula(get.models(test.local.modBGNfbin2.full.b, i)[[1]]))[3]
local.modBGNfbin2.sel.table<-local.modBGNfbin2.sel.table[,c(6,1,2,3,4,5)]


local.modBGNfbin2.full.b.simres <- simulateResiduals(get.models(test.local.modBGNfbin2.full.b, 1)[[1]])
plot(local.modBGNfbin2.full.b.simres)
testDispersion(local.modBGNfbin2.full.b.simres)

local.modBGNfbin2.full.b.alleffects <- allEffects(get.models(test.local.modBGNfbin2.full.b, 1)[[1]], rescale.axis=F) #, residuals = T
plot(local.modBGNfbin2.full.b.alleffects)


local.modBGNfbin2.full.b.t <- broom.mixed::tidy(get.models(test.local.modBGNfbin2.full.b, 1)[[1]], effects = "fixed", conf.int = TRUE)
local.modBGNfbin2.full.b.t <- transform(local.modBGNfbin2.full.b.t, term = sprintf("%s.%s", component, term))
local.modBGNfbin2.full.b.t$group <- "BGN [4-12kHz]"
local.modBGNfbin2.full.b.t$col <- ifelse(local.modBGNfbin2.full.b.t$conf.low < 0 & local.modBGNfbin2.full.b.t$conf.high < 0, "#e65a49", 
                                   ifelse(local.modBGNfbin2.full.b.t$conf.low > 0 & local.modBGNfbin2.full.b.t$conf.high > 0, "#55c0f0", "#dbdfde"))



#Response: Acoustic entropy (S2N) on total frequency [0.3 - 22.1kS2Nz]
#hist(CN.data.local.modelfit$S2NTotal); range(CN.data.local.modelfit$S2NTotal)
#plot(family.test <- ecdf(CN.data.local.modelfit$S2NTotal))
#plot(fitdist(CN.data.local.modelfit$S2NTotal, "pois"))


#modS2N.null <- glmmTMB(S2NTotal ~ 1, family=genpois(link = "log"), CN.data.modelfit)


local.modS2N.full <- glmmTMB(S2NTotal ~ n_trees + n_vines + n_palmtrees + (1|Local:time_period),
                       family=genpois(link = "log"), CN.data.local.modelfit)

check_collinearity(local.modS2N.full)

local.modS2N.full.b <- glmmTMB(S2NTotal ~ n_trees + n_vines + n_palmtrees + (1|Local:time_period) + ar1(times - 1|Year),
                         family=genpois(link = "log"), CN.data.local.modelfit)

summary(local.modS2N.full.b)
#r.squaredGLMM(modS2N.full.b)

test.local.modS2N.full.b <- dredge(local.modS2N.full.b)
#subset(test.modS2N.full.b, delta < 4)
summary(get.models(test.local.modS2N.full.b, 1)[[1]])
#r.squaredGLMM(get.models(test.modS2N.full.b, 1)[[1]])

sink("results/supplementar/S2NTotal_local_model_summary.csv")
print(summary(get.models(test.local.modS2N.full.b, 1)[[1]]))
sink()

local.modS2N.sel.table<-as.data.frame(test.local.modS2N.full.b)[1:5,6:10]
local.modS2N.sel.table[,2:3]<- round(local.modS2N.sel.table[,2:3],2)
local.modS2N.sel.table[,4:5]<- round(local.modS2N.sel.table[,4:5],3)
names(local.modS2N.sel.table)[1] = "K"
local.modS2N.sel.table$Model<-NA
for(i in 1:nrow(local.modS2N.sel.table)) local.modS2N.sel.table$Model[i]<-as.character(formula(get.models(test.local.modS2N.full.b, i)[[1]]))[3]
local.modS2N.sel.table<-local.modS2N.sel.table[,c(6,1,2,3,4,5)]


local.modS2N.full.b.simres <- simulateResiduals(get.models(test.local.modS2N.full.b, 1)[[1]])
plot(local.modS2N.full.b.simres)
testDispersion(local.modS2N.full.b.simres)

local.modS2N.full.b.alleffects <- allEffects(get.models(test.local.modS2N.full.b, 1)[[1]], rescale.axis=F) #, residuals = T
plot(local.modS2N.full.b.alleffects)


local.modS2N.full.b.t <- broom.mixed::tidy(get.models(test.local.modS2N.full.b, 1)[[1]], effects = "fixed", conf.int = TRUE)
local.modS2N.full.b.t <- transform(local.modS2N.full.b.t, term = sprintf("%s.%s", component, term))
local.modS2N.full.b.t$group <- "S2N [0.3-22.1kHz]"
local.modS2N.full.b.t$col <- ifelse(local.modS2N.full.b.t$conf.low < 0 & local.modS2N.full.b.t$conf.high < 0, "#e65a49", 
                              ifelse(local.modS2N.full.b.t$conf.low > 0 & local.modS2N.full.b.t$conf.high > 0, "#55c0f0", "#dbdfde"))



#Response: Acoustic entropy (S2N) on frequency bin 2 [4 - 12kS2Nz]

#modS2Nfbin2.null <- glmmTMB(S2Nfbin2 ~ 1, family=beta_family(), CN.data.local.modelfit)

#options(na.action = "na.fail")
local.modS2Nfbin2.full <- glmmTMB(S2Nfbin2 ~ n_trees + n_vines + n_palmtrees + (1|Local:time_period),
                            family=genpois(link = "log"), CN.data.local.modelfit)

check_collinearity(local.modS2Nfbin2.full)


local.modS2Nfbin2.full.b <- glmmTMB(S2Nfbin2 ~ n_trees + n_vines + n_palmtrees + (1|Local:time_period) + ar1(times - 1|Year),
                              family=genpois(link = "log"), CN.data.local.modelfit)

summary(local.modS2Nfbin2.full.b)
#r.squaredGLMM(modS2Nfbin2.full.b)

test.local.modS2Nfbin2.full.b <- dredge(local.modS2Nfbin2.full.b)
#subset(test.modS2Nfbin2.full.b, delta < 4)
summary(get.models(test.local.modS2Nfbin2.full.b, 1)[[1]])
#r.squaredGLMM(get.models(test.modS2Nfbin2.full.b, 1)[[1]])

sink("results/supplementar/S2Nfbin2_local_model_summary.csv")
print(summary(get.models(test.local.modS2Nfbin2.full.b, 1)[[1]]))
sink()

local.modS2Nfbin2.sel.table<-as.data.frame(test.local.modS2Nfbin2.full.b)[1:5,6:10]
local.modS2Nfbin2.sel.table[,2:3]<- round(local.modS2Nfbin2.sel.table[,2:3],2)
local.modS2Nfbin2.sel.table[,4:5]<- round(local.modS2Nfbin2.sel.table[,4:5],3)
names(local.modS2Nfbin2.sel.table)[1] = "K"
local.modS2Nfbin2.sel.table$Model<-NA
for(i in 1:nrow(local.modS2Nfbin2.sel.table)) local.modS2Nfbin2.sel.table$Model[i]<-as.character(formula(get.models(test.local.modS2Nfbin2.full.b, i)[[1]]))[3]
local.modS2Nfbin2.sel.table<-local.modS2Nfbin2.sel.table[,c(6,1,2,3,4,5)]


local.modS2Nfbin2.full.b.simres <- simulateResiduals(get.models(test.local.modS2Nfbin2.full.b, 1)[[1]])
plot(local.modS2Nfbin2.full.b.simres)
testDispersion(local.modS2Nfbin2.full.b.simres)

local.modS2Nfbin2.full.b.alleffects <- allEffects(get.models(test.local.modS2Nfbin2.full.b, 1)[[1]], rescale.axis=F) #, residuals = T
plot(local.modS2Nfbin2.full.b.alleffects)


local.modS2Nfbin2.full.b.t <- broom.mixed::tidy(get.models(test.local.modS2Nfbin2.full.b, 1)[[1]], effects = "fixed", conf.int = TRUE)
local.modS2Nfbin2.full.b.t <- transform(local.modS2Nfbin2.full.b.t, term = sprintf("%s.%s", component, term))
local.modS2Nfbin2.full.b.t$group <- "S2N [4-12kHz]"
local.modS2Nfbin2.full.b.t$col <- ifelse(local.modS2Nfbin2.full.b.t$conf.low < 0 & local.modS2Nfbin2.full.b.t$conf.high < 0, "#e65a49", 
                                   ifelse(local.modS2Nfbin2.full.b.t$conf.low > 0 & local.modS2Nfbin2.full.b.t$conf.high > 0, "#55c0f0", "#dbdfde"))



#Response: Acoustic entropy (AA) on total frequency [0.3 - 22.1kAAz]
#hist(CN.data.local.modelfit$AATotal); range(CN.data.local.modelfit$AATotal)
#plot(family.test <- ecdf(CN.data.local.modelfit$AATotal))
#plot(fitdist(CN.data.local.modelfit$AATotal, "beta"))


#modAA.null <- glmmTMB(AATotal ~ 1, family=beta_family(), CN.data.local.modelfit)

#options(na.action = "na.fail")
local.modAA.full <- glmmTMB(AATotal ~ n_trees + n_vines + n_palmtrees + (1|Local:time_period),
                      family=beta_family(), CN.data.local.modelfit)

check_collinearity(local.modAA.full)


local.modAA.full.b <- glmmTMB(AATotal ~ n_trees + n_vines + n_palmtrees + (1|Local:time_period) + ar1(times - 1|Year),
                        family=beta_family(), CN.data.local.modelfit)

summary(local.modAA.full.b)
#r.squaredGLMM(modAA.full.b)

test.local.modAA.full.b <- dredge(local.modAA.full.b)
#subset(test.modAA.full.b, delta < 4)
summary(get.models(test.local.modAA.full.b, 1)[[1]])
#r.squaredGLMM(get.models(test.modAA.full.b, 1)[[1]])

sink("results/supplementar/AATotal_local_model_summary.csv")
print(summary(get.models(test.local.modAA.full.b, 1)[[1]]))
sink()

local.modAA.sel.table<-as.data.frame(test.local.modAA.full.b)[1:5,6:10]
local.modAA.sel.table[,2:3]<- round(local.modAA.sel.table[,2:3],2)
local.modAA.sel.table[,4:5]<- round(local.modAA.sel.table[,4:5],3)
names(local.modAA.sel.table)[1] = "K"
local.modAA.sel.table$Model<-NA
for(i in 1:nrow(local.modAA.sel.table)) local.modAA.sel.table$Model[i]<-as.character(formula(get.models(test.local.modAA.full.b, i)[[1]]))[3]
local.modAA.sel.table<-local.modAA.sel.table[,c(6,1,2,3,4,5)]


local.modAA.full.b.simres <- simulateResiduals(get.models(test.local.modAA.full.b, 1)[[1]])
plot(local.modAA.full.b.simres)
testDispersion(local.modAA.full.b.simres)

local.modAA.full.b.alleffects <- allEffects(get.models(test.local.modAA.full.b, 1)[[1]], rescale.axis=F) #, residuals = T
plot(local.modAA.full.b.alleffects)


local.modAA.full.b.t <- broom.mixed::tidy(get.models(test.local.modAA.full.b, 1)[[1]], effects = "fixed", conf.int = TRUE)
local.modAA.full.b.t <- transform(local.modAA.full.b.t, term = sprintf("%s.%s", component, term))
local.modAA.full.b.t$group <- "AA [0.3-22.1kHz]"
local.modAA.full.b.t$col <- ifelse(local.modAA.full.b.t$conf.low < 0 & local.modAA.full.b.t$conf.high < 0, "#e65a49", 
                             ifelse(local.modAA.full.b.t$conf.low > 0 & local.modAA.full.b.t$conf.high > 0, "#55c0f0", "#dbdfde"))



#Response: Acoustic entropy (AA) on frequency bin 2 [4 - 12kAAz]
#hist(CN.data.local.modelfit$AAfbin2); range(CN.data.local.modelfit$AAfbin2)
#plot(family.test <- ecdf(CN.data.local.modelfit$AAfbin2))
#plot(fitdist(CN.data.local.modelfit$AAfbin2, "beta"))

CN.data.local.modelfit$AAfbin2b <- rescale(CN.data.local.modelfit$AAfbin2, to=c(.001,.999))

#modAAfbin2.null <- glmmTMB(AAfbin2 ~ 1, family=beta_family(), CN.data.local.modelfit)

#options(na.action = "na.fail")
local.modAAfbin2.full <- glmmTMB(AAfbin2b ~ n_trees + n_vines + n_palmtrees + (1|Local:time_period),
                           family=beta_family(), CN.data.local.modelfit)

check_collinearity(local.modAAfbin2.full)


local.modAAfbin2.full.b <- glmmTMB(AAfbin2b ~ n_trees + n_vines + n_palmtrees + (1|Local:time_period) + ar1(times - 1|Year),
                             family=beta_family(), CN.data.local.modelfit)

summary(local.modAAfbin2.full.b)
#r.squaredGLMM(modAAfbin2.full.b)

test.local.modAAfbin2.full.b <- dredge(local.modAAfbin2.full.b)
#subset(test.modAAfbin2.full.b, delta < 4)
summary(get.models(test.local.modAAfbin2.full.b, 1)[[1]])
#r.squaredGLMM(get.models(test.modAAfbin2.full.b, 1)[[1]])

sink("results/supplementar/AAfbin2_local_model_summary.csv")
print(summary(get.models(test.local.modAAfbin2.full.b, 1)[[1]]))
sink()

local.modAAfbin2.sel.table<-as.data.frame(test.local.modAAfbin2.full.b)[1:5,6:10]
local.modAAfbin2.sel.table[,2:3]<- round(local.modAAfbin2.sel.table[,2:3],2)
local.modAAfbin2.sel.table[,4:5]<- round(local.modAAfbin2.sel.table[,4:5],3)
names(local.modAAfbin2.sel.table)[1] = "K"
local.modAAfbin2.sel.table$Model<-NA
for(i in 1:nrow(local.modAAfbin2.sel.table)) local.modAAfbin2.sel.table$Model[i]<-as.character(formula(get.models(test.local.modAAfbin2.full.b, i)[[1]]))[3]
local.modAAfbin2.sel.table<-local.modAAfbin2.sel.table[,c(6,1,2,3,4,5)]


local.modAAfbin2.full.b.simres <- simulateResiduals(get.models(test.local.modAAfbin2.full.b, 1)[[1]])
plot(local.modAAfbin2.full.b.simres)
testDispersion(local.modAAfbin2.full.b.simres)

local.modAAfbin2.full.b.alleffects <- allEffects(get.models(test.local.modAAfbin2.full.b, 1)[[1]], rescale.axis=F) #, residuals = T
plot(local.modAAfbin2.full.b.alleffects)


local.modAAfbin2.full.b.t <- broom.mixed::tidy(get.models(test.local.modAAfbin2.full.b, 1)[[1]], effects = "fixed", conf.int = TRUE)
local.modAAfbin2.full.b.t <- transform(local.modAAfbin2.full.b.t, term = sprintf("%s.%s", component, term))
local.modAAfbin2.full.b.t$group <- "AA [4-12kHz]"
local.modAAfbin2.full.b.t$col <- ifelse(local.modAAfbin2.full.b.t$conf.low < 0 & local.modAAfbin2.full.b.t$conf.high < 0, "#e65a49", 
                                  ifelse(local.modAAfbin2.full.b.t$conf.low > 0 & local.modAAfbin2.full.b.t$conf.high > 0, "#55c0f0", "#dbdfde"))



options(na.action = "na.omit")

sel.table <- rbind(local.modH.sel.table, local.modHfbin2.sel.table,
                   #local.modADI.sel.table, local.modADIfbin2.sel.table,
                   local.modAEI.sel.table, local.modAEIfbin2.sel.table,
                   local.modBGN.sel.table, local.modBGNfbin2.sel.table,
                   local.modS2N.sel.table, local.modS2Nfbin2.sel.table,
                   local.modAA.sel.table, local.modAAfbin2.sel.table)
write.csv(sel.table,"results/supplementar/local_model_selection_table.csv", row.names = F)



mod.t <- rbind(mod.t,
               local.modH.full.b.t, local.modHfbin2.full.b.t,
               #local.modADI.full.b.t, local.modADIfbin2.full.b.t,
               local.modAEI.full.b.t, local.modAEIfbin2.full.b.t,
               local.modBGN.full.b.t, local.modBGNfbin2.full.b.t,
               local.modS2N.full.b.t, local.modS2Nfbin2.full.b.t,
               local.modAA.full.b.t, local.modAAfbin2.full.b.t)

mod.t <- mod.t %>% filter(term!="cond.(Intercept)" & 
                            term!="cond.sd__(Intercept)" &
                            component!="zi") %>% 
  mutate(group=factor(group, levels = c("H [0.3-22.1kHz]", "H [4-12kHz]", 
                                        "ADI [0.3-22.1kHz]", "ADI [4-12kHz]",
                                        "AEI [0.3-22.1kHz]", "AEI [4-12kHz]",
                                        "BGN [0.3-22.1kHz]", "BGN [4-12kHz]",
                                        "S2N [0.3-22.1kHz]", "S2N [4-12kHz]",
                                        "AA [0.3-22.1kHz]", "AA [4-12kHz]")),
         term=factor(term, levels = c("cond.n_trees",
                                      "cond.n_vines",
                                      "cond.n_palmtrees",
                                      "cond.z_alt", 
                                      "cond.z_distwater", 
                                      "cond.z_distedge", 
                                      "cond.z_treeheight",
                                      "cond.z_NDVI", 
                                      "cond.z_precipitation", 
                                      "cond.z_temperature",
                                      "cond.z_humidity", 
                                      "cond.z_wind")))

#write.csv(mod.t, "results/effect_sizes.csv", row.names = F)


mod.t.graph <- mod.t %>% filter(group %in% unique(grep("ADI|S2N", mod.t$group, value = T, invert = T))) %>% 
  ggplot()+
  geom_point(aes(x = estimate, y = term, col = col), size = 12)+
  geom_segment(aes(x = conf.low, xend = conf.high,
                   y = term, yend = term, col = col), size = 7)+
  geom_vline(xintercept=0, lty=2)+
  scale_color_manual(values = c("#55c0f0", "#dbdfde", "#e65a49"))+
  scale_y_discrete(labels = c("cond.n_trees"="Number of trees",
                              "cond.n_vines"="Number of vines",
                              "cond.n_palmtrees"="Number of palmtrees",
                              "cond.z_alt"="Altitude",
                              "cond.z_distwater"="Distance from water",
                              "cond.z_distedge"="Distance from edge",
                              "cond.z_treeheight"="Tree height",
                              "cond.z_NDVI"="NDVI",
                              "cond.z_precipitation"="Precipitation",
                              "cond.z_temperature"="Temperature",
                              "cond.z_humidity" = "Air humidity",
                              "cond.z_wind"="Wind speed"))+
  facet_wrap(~group, ncol = 2)+#, labeller = labeller(group = c("H" = "H [0.3 - 22.1kHz]",
  #                             "Hfbin2" = "H [4 - 12kHz]",
  #                             "AEI" = "AEI [0.3 - 22.1kHz]",
  #                             "AEIfbin2" = "AEI [4 - 12kHz]",
  #                             "BGN" = "BGN [0.3 - 22.1kHz]",
  #                             "BGNb2" = "BGN [4 - 12kHz]",
  #                             "AA" = "AA [0.3 - 22.1kHz]",
  #                             "AAfbin2" = "AA [4 - 12kHz]")))+
  xlab("")+ ylab("")+
  theme(plot.title = element_text(family = "serif", size = 48, hjust = .5),
        strip.text = element_text(family = "serif", size = 36, hjust = .5),
        axis.title = element_text(family = "serif", size = 36),
        axis.text.x = element_text(family = "serif", size = 32),
        axis.text.y = element_text(family = "serif", size = 32),
        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2),
        panel.background = element_blank(), legend.position = "none")




png("results/environmentaleffects.png", width = 23.38583, height = 33.11024, units = "in", res = 300, pointsize = 4)
mod.t.graph
dev.off()

#=================================|end










































































#==============================| previous modeling approach

#subsetting
CN.data.regional.modelfit <- CN.data.regional %>% group_by(Local, time_period) %>% sample_n(size = 48, replace = T) %>% ungroup()

#hist(CN.data.regional.modelfit$AATotal); range(CN.data.regional.modelfit$AATotal)
#plot(family.test <- ecdf(CN.data.regional.modelfit$AATotal))
#plot(fitdist(CN.data.regional.modelfit$AATotal, "beta"))
res1.AA <- glmmTMB(AATotal ~ distwater + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
summary(res1.AA)   

res2.AA <- glmmTMB(AATotal ~ distedge + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
summary(res2.AA)    

#res3.AA <- glmmTMB(AATotal ~ distcanga + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
#summary(res3.AA)     

#res4.AA <- glmmTMB(AATotal ~ distminning + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
#summary(res4.AA)      

#res5.AA <- glmmTMB(AATotal ~ treeheight + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
#summary(res5.AA)       

#res6.AA <- glmmTMB(AATotal ~ NDVI + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
#summary(res6.AA)        

#res7.AA <- glmmTMB(AATotal ~ precipitation + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
#summary(res7.AA)         

res8.AA <- glmmTMB(AATotal ~ temperature + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
summary(res8.AA)          

res9.AA <- glmmTMB(AATotal ~ humidity + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
summary(res9.AA)           

#res10.AA <- glmmTMB(AATotal ~ wind + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
#summary(res10.AA)            

res11.AA <- glmmTMB(AATotal ~ distwater + humidity + distedge + temperature + (1|Local:time_period), 
                    family=beta_family(), CN.data.regional.modelfit)
summary(res11.AA) 

res11.AA.simres <- simulateResiduals(res11.AA)
plot(res11.AA.simres)
res11.AA.alleffects <- allEffects(res11.AA, residuals = T)
plot(res11.AA.alleffects)

res11.AA.t <- broom.mixed::tidy(res11.AA, conf.int = TRUE)
res11.AA.t <- transform(res11.AA.t, term = sprintf("%s.%s", component, term))
res11.AA.t$group <- "AA"
res11.AA.t$col <- ifelse(res11.AA.t$conf.low < 0 & res11.AA.t$conf.high < 0, "#e65a49", 
                           ifelse(res11.AA.t$conf.low > 0 & res11.AA.t$conf.high > 0, "#55c0f0", "#dbdfde"))

res11.AA.graph <- ggplot(res11.AA.t[1:5, ])+
  geom_point(aes(x = estimate, y = term, col = col), size = 5)+
  geom_segment(aes(x = conf.low, xend = conf.high,
                   y = term, yend = term, col = col), size = 2)+
  geom_vline(xintercept=0, lty=2)+
  scale_color_manual(values = c("#55c0f0", "#dbdfde", "#e65a49"))+
  scale_y_discrete(labels = c("cond.(Intercept)"="Intercept", "cond.distwater"="Distance from water",
                              "cond.humidity" = "Air humidity", "cond.distedge"="Distance from edge",
                              "cond.temperature"="Temperature"))+
  labs(title = "AA") + xlab("") + ylab("") +
  theme(plot.title = element_text(family = "serif", size = 26, hjust = .5),
        axis.title = element_text(family = "serif", size = 22),
        axis.text.x = element_text(family = "serif", size = 20),
        axis.text.y = element_text(family = "serif", size = 20),
        axis.line.y = element_line(size = 1), axis.line.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1), axis.ticks.x = element_line(size = 1),
        panel.background = element_blank(), legend.position = "none")


##

#hist(CN.data.regional.modelfit$AEITotal); range(CN.data.regional.modelfit$AEITotal)
#plot(family.test <- ecdf(CN.data.regional.modelfit$AEITotal))
#plot(fitdist(CN.data.regional.modelfit$AEITotal, "beta"))


res1.AEI <- glmmTMB(AEITotal ~ distwater + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
summary(res1.AEI)   

#res2.AEI <- glmmTMB(AEITotal ~ distedge + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
#summary(res2.AEI)    

#res3.AEI <- glmmTMB(AEITotal ~ distcanga + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
#summary(res3.AEI)     

#res4.AEI <- glmmTMB(AEITotal ~ distminning + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
#summary(res4.AEI)      

#res5.AEI <- glmmTMB(AEITotal ~ treeheight + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
#summary(res5.AEI)       

res6.AEI <- glmmTMB(AEITotal ~ NDVI + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
summary(res6.AEI)        

#res7.AEI <- glmmTMB(AEITotal ~ precipitation + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
#summary(res7.AEI)         

res8.AEI <- glmmTMB(AEITotal ~ temperature + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
summary(res8.AEI)          

res9.AEI <- glmmTMB(AEITotal ~ humidity + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
summary(res9.AEI)           

#res10.AEI <- glmmTMB(AEITotal ~ wind + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
#summary(res10.AEI)            

res11.AEI <- glmmTMB(AEITotal ~ distwater + NDVI + temperature * humidity + (1|Local:time_period), 
                    family=beta_family(), CN.data.regional.modelfit)
summary(res11.AEI) 

res11.AEI.simres <- simulateResiduals(res11.AEI)
plot(res11.AEI.simres)
plot(allEffects(res11.AEI, residuals = T))

res11.AEI.t <- broom.mixed::tidy(res11.AEI, conf.int = TRUE)
res11.AEI.t <- transform(res11.AEI.t, term = sprintf("%s.%s", component, term))
res11.AEI.t$group <- "AEI"
res11.AEI.t$col <- ifelse(res11.AEI.t$conf.low < 0 & res11.AEI.t$conf.high < 0, "#e65a49", 
                         ifelse(res11.AEI.t$conf.low > 0 & res11.AEI.t$conf.high > 0, "#55c0f0", "#dbdfde"))

res11.BGN.graph <- ggplot(res11.AEI.t[1:6, ])+
  geom_point(aes(x = estimate, y = term, col = col), size = 5)+
  geom_segment(aes(x = conf.low, xend = conf.high,
                   y = term, yend = term, col = col), size = 2)+
  geom_vline(xintercept=0, lty=2)+
  scale_color_manual(values = c("#55c0f0", "#dbdfde", "#e65a49"))+
  scale_y_discrete(labels = c("cond.(Intercept)"="Intercept", "cond.distwater"="Distance from water",
                              "cond.NDVI"="NDVI", "cond.humidity" = "Air humidity",
                              "cond.temperature"="Temperature",
                              "cond.temperature:humidity"="Temperature * Air humidity"))+
  labs(title = "AEI") + xlab("") + ylab("") +
  theme(plot.title = element_text(family = "serif", size = 26, hjust = .5),
        axis.title = element_text(family = "serif", size = 22),
        axis.text.x = element_text(family = "serif", size = 20),
        axis.text.y = element_text(family = "serif", size = 20),
        axis.line.y = element_line(size = 1), axis.line.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1), axis.ticks.x = element_line(size = 1),
        panel.background = element_blank(), legend.position = "none")


##

#hist((CN.data.regional.modelfit$BGNTotal*(-1))); range((CN.data.regional.modelfit$BGNTotal*(-1)))
#plot(fitdist((CN.data.regional.modelfit$BGNTotal*(-1)), "pois"))  

res1.BGN <- glmmTMB((BGNTotal*(-1)) ~ distwater + (1|Local:time_period), family=poisson(link = "identity"), CN.data.regional.modelfit)
summary(res1.BGN)   

#res2.BGN <- glmmTMB((BGNTotal*(-1)) ~ distedge + (1|Local:time_period), family=poisson(link = "identity"), CN.data.regional.modelfit)
#summary(res2.BGN)    

res3.BGN <- glmmTMB((BGNTotal*(-1)) ~ distcanga + (1|Local:time_period), family=poisson(link = "identity"), CN.data.regional.modelfit)
summary(res3.BGN)     

#res4.BGN <- glmmTMB((BGNTotal*(-1)) ~ distminning + (1|Local:time_period), family=poisson(link = "identity"), CN.data.regional.modelfit)
#summary(res4.BGN)      

#res5.BGN <- glmmTMB((BGNTotal*(-1)) ~ treeheight + (1|Local:time_period), family=poisson(link = "identity"), CN.data.regional.modelfit)
#summary(res5.BGN)       

res6.BGN <- glmmTMB((BGNTotal*(-1)) ~ NDVI + (1|Local:time_period), family=poisson(link = "identity"), CN.data.regional.modelfit)
summary(res6.BGN)        

#res7.BGN <- glmmTMB((BGNTotal*(-1)) ~ precipitation + (1|Local:time_period), family=poisson(link = "identity"), CN.data.regional.modelfit)
#summary(res7.BGN)         

#res8.BGN <- glmmTMB((BGNTotal*(-1)) ~ temperature + (1|Local:time_period), family=poisson(link = "identity"), CN.data.regional.modelfit)
#summary(res8.BGN)          

#res9.BGN <- glmmTMB((BGNTotal*(-1)) ~ humidity + (1|Local:time_period), family=poisson(link = "identity"), CN.data.regional.modelfit)
#summary(res9.BGN)           

#res10.BGN <- glmmTMB((BGNTotal*(-1)) ~ wind + (1|Local:time_period), family=poisson(link = "identity"), CN.data.regional.modelfit)
#summary(res10.BGN)            

res11.BGN <- glmmTMB((BGNTotal*(-1)) ~ distwater + distcanga + NDVI + (1|Local:time_period), 
                     family=poisson(link = "identity"), CN.data.regional.modelfit)

summary(res11.BGN) 

res11.BGN.simres <- simulateResiduals(res11.BGN)
plot(res11.BGN.simres)
plot(allEffects(res11.BGN, residuals = T))

res11.BGN.t <- broom.mixed::tidy(res11.BGN, conf.int = TRUE)
res11.BGN.t <- transform(res11.BGN.t, term = sprintf("%s.%s", component, term))
res11.BGN.t$group <- "BGN"
res11.BGN.t$col <- ifelse(res11.BGN.t$conf.low < 0 & res11.BGN.t$conf.high < 0, "#e65a49", 
                          ifelse(res11.BGN.t$conf.low > 0 & res11.BGN.t$conf.high > 0, "#55c0f0", "#dbdfde"))

res11.BGN.graph <- ggplot(res11.BGN.t[1:4, ])+
  geom_point(aes(x = estimate, y = term, col = col), size = 5)+
  geom_segment(aes(x = conf.low, xend = conf.high,
                   y = term, yend = term, col = col), size = 2)+
  geom_vline(xintercept=0, lty=2)+
  scale_color_manual(values = c("#55c0f0", "#dbdfde", "#e65a49"))+
  scale_y_discrete(labels = c("cond.(Intercept)"="Intercept", "cond.distwater"="Distance from water",
                                               "cond.NDVI"="NDVI", "cond.distcanga" = "Distance from canga"))+
  labs(title = "BGN") + xlab("") + ylab("") +
  theme(plot.title = element_text(family = "serif", size = 26, hjust = .5),
        axis.title = element_text(family = "serif", size = 22),
        axis.text.x = element_text(family = "serif", size = 20),
        axis.text.y = element_text(family = "serif", size = 20),
        axis.line.y = element_line(size = 1), axis.line.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1), axis.ticks.x = element_line(size = 1),
        panel.background = element_blank(), legend.position = "none")

##

#hist(CN.data.regional.modelfit$HTotal); range(CN.data.regional.modelfit$HTotal)
#plot(fitdist(CN.data.regional.modelfit$HTotal, "beta"))
res1.H <- glmmTMB(HTotal ~ distwater + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
summary(res1.H)   

res2.H <- glmmTMB(HTotal ~ distedge + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
summary(res2.H)    

#res3.H <- glmmTMB(HTotal ~ distcanga + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
#summary(res3.H)     

#res4.H <- glmmTMB(HTotal ~ distminning + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
#summary(res4.H)      

res5.H <- glmmTMB(HTotal ~ treeheight + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
summary(res5.H)       

#res6.H <- glmmTMB(HTotal ~ NDVI + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
#summary(res6.H)        

#res7.H <- glmmTMB(HTotal ~ precipitation + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
#summary(res7.H)         

res8.H <- glmmTMB(HTotal ~ temperature + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
summary(res8.H)          

res9.H <- glmmTMB(HTotal ~ humidity + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
summary(res9.H)           

#res10.H <- glmmTMB(HTotal ~ wind + (1|Local:time_period), family=beta_family(), CN.data.regional.modelfit)
#summary(res10.H)            

res11.H <- glmmTMB(HTotal ~ distwater + humidity + distedge + temperature + treeheight + (1|Local:time_period), 
                   family=beta_family(), CN.data.regional.modelfit)
summary(res11.H) 

res11.H.simres <- simulateResiduals(res11.H)
plot(res11.H.simres)
res11.H.alleffects <- allEffects(res11.H, residuals = T)
plot(res11.H.alleffects)

res11.H.t <- broom.mixed::tidy(res11.H, conf.int = TRUE)
res11.H.t <- transform(res11.H.t, term = sprintf("%s.%s", component, term))
res11.H.t$group <- "H"
res11.H.t$col <- ifelse(res11.H.t$conf.low < 0 & res11.H.t$conf.high < 0, "#e65a49", 
                          ifelse(res11.H.t$conf.low > 0 & res11.H.t$conf.high > 0, "#55c0f0", "#dbdfde"))

res11.H.graph <- ggplot(res11.H.t[1:6, ])+
  geom_point(aes(x = estimate, y = term, col = col), size = 5)+
  geom_segment(aes(x = conf.low, xend = conf.high,
                   y = term, yend = term, col = col), size = 2)+
  geom_vline(xintercept=0, lty=2)+
  scale_color_manual(values = c("#55c0f0", "#dbdfde", "#e65a49"))+
  scale_y_discrete(labels = c("cond.(Intercept)"="Intercept", "cond.distwater"="Distance from water",
                              "cond.humidity" = "Air humidity", "cond.distedge"="Distance from edge",
                              "cond.temperature"="Temperature",
                              "cond.treeheight"="Tree height"))+
  labs(title = "H") + xlab("") + ylab("") +
  theme(plot.title = element_text(family = "serif", size = 26, hjust = .5),
        axis.title = element_text(family = "serif", size = 22),
        axis.text.x = element_text(family = "serif", size = 20),
        axis.text.y = element_text(family = "serif", size = 20),
        axis.line.y = element_line(size = 1), axis.line.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1), axis.ticks.x = element_line(size = 1),
        panel.background = element_blank(), legend.position = "none")


##

png("results/regionalenvironmentaleffects.png", width = 1280, height = 720, units = "px", bg = "transparent")
plot_grid(res11.H.graph, res11.AEI.graph, res11.BGN.graph, res11.AA.graph, nrow=2, align = "hv")
dev.off()

rm(list= ls()[(ls() %in% grep("res", ls(),value = T))])
gc()
#####



##### modelling the effect of environmental traits on soundscape index ####
## preparing data -- local
CN.data.local <- CN.data[!is.na(CN.data$n_trees),c(1:16,44:46)]

#subsetting
CN.data.local.modelfit <- CN.data.local %>% group_by(Local, time_period) %>% sample_n(size = 45, replace = T) %>% ungroup()

#hist(CN.data.local.modelfit$AATotal); range(CN.data.local.modelfit$AATotal)
#plot(fitdist(CN.data.local.modelfit$AATotal, "beta"))
res1.AA <- glmmTMB(AATotal ~ n_trees + (1|Local:time_period), family=beta_family(), CN.data.local.modelfit)
summary(res1.AA)   

#res2.AA <- glmmTMB(AATotal ~ n_vines + (1|Local:time_period), family=beta_family(), CN.data.local.modelfit)
#summary(res2.AA)    

res3.AA <- glmmTMB(AATotal ~ n_palmtrees + (1|Local:time_period), family=beta_family(), CN.data.local.modelfit)
summary(res3.AA)     

res4.AA <- glmmTMB(AATotal ~ n_trees + n_palmtrees + (1|Local:time_period), family=beta_family(), CN.data.local.modelfit)
summary(res4.AA)      


res4.AA.simres <- simulateResiduals(res4.AA)
plot(res4.AA.simres)
res4.AA.alleffects <- allEffects(res4.AA, residuals = T)
plot(res4.AA.alleffects)

res4.AA.t <- broom.mixed::tidy(res4.AA, conf.int = TRUE)
res4.AA.t <- transform(res4.AA.t, term = sprintf("%s.%s", component, term))
res4.AA.t$group <- "AA"
res4.AA.t$col <- ifelse(res4.AA.t$conf.low < 0 & res4.AA.t$conf.high < 0, "#e65a49", 
                        ifelse(res4.AA.t$conf.low > 0 & res4.AA.t$conf.high > 0, "#55c0f0", "#dbdfde"))

res4.AA.graph <- ggplot(res4.AA.t[1:3, ])+
  geom_point(aes(x = estimate, y = term, col = col), size = 5)+
  geom_segment(aes(x = conf.low, xend = conf.high,
                   y = term, yend = term, col = col), size = 2)+
  geom_vline(xintercept=0, lty=2)+
  scale_color_manual(values = c("#55c0f0", "#dbdfde", "#e65a49"))+
  scale_y_discrete(labels = c("cond.(Intercept)"="Intercept", "cond.n_trees"="Number of trees",
                              "cond.n_palmtrees"="Number of palmtrees"))+
  labs(title = "AA") + xlab("") + ylab("") +
  theme(plot.title = element_text(family = "serif", size = 26, hjust = .5),
        axis.title = element_text(family = "serif", size = 22),
        axis.text.x = element_text(family = "serif", size = 20),
        axis.text.y = element_text(family = "serif", size = 20),
        axis.line.y = element_line(size = 1), axis.line.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1), axis.ticks.x = element_line(size = 1),
        panel.background = element_blank())


##

#hist(CN.data.local.modelfit$AEITotal); range(CN.data.local.modelfit$AEITotal)
#plot(fitdist(CN.data.local.modelfit$AEITotal, "beta"))

res1.AEI <- glmmTMB(AEITotal ~ n_trees + (1|Local:time_period), family=beta_family(), CN.data.local.modelfit)
summary(res1.AEI)   

res2.AEI <- glmmTMB(AEITotal ~ n_vines + (1|Local:time_period), family=beta_family(), CN.data.local.modelfit)
summary(res2.AEI)    

res3.AEI <- glmmTMB(AEITotal ~ n_palmtrees + (1|Local:time_period), family=beta_family(), CN.data.local.modelfit)
summary(res3.AEI)     

res4.AEI <- glmmTMB(AEITotal ~ n_trees + n_palmtrees + (1|Local:time_period), family=beta_family(), CN.data.local.modelfit)
summary(res4.AEI)      


res4.AEI.simres <- simulateResiduals(res4.AEI)
plot(res4.AEI.simres)
res4.AEI.alleffects <- allEffects(res4.AEI, residuals = T)
plot(res4.AEI.alleffects)

res4.AEI.t <- broom.mixed::tidy(res4.AEI, conf.int = TRUE)
res4.AEI.t <- transform(res4.AEI.t, term = sprintf("%s.%s", component, term))
res4.AEI.t$group <- "AEI"
res4.AEI.t$col <- ifelse(res4.AEI.t$conf.low < 0 & res4.AEI.t$conf.high < 0, "#e65a49", 
                        ifelse(res4.AEI.t$conf.low > 0 & res4.AEI.t$conf.high > 0, "#55c0f0", "#dbdfde"))

res4.AEI.graph <- ggplot(res4.AEI.t[1:3, ])+
  geom_point(aes(x = estimate, y = term, col = col), size = 5)+
  geom_segment(aes(x = conf.low, xend = conf.high,
                   y = term, yend = term, col = col), size = 2)+
  geom_vline(xintercept=0, lty=2)+
  scale_color_manual(values = c("#55c0f0", "#dbdfde", "#e65a49"))+
  scale_y_discrete(labels = c("cond.(Intercept)"="Intercept", "cond.n_trees"="Number of trees",
                              "cond.n_palmtrees"="Number of palmtrees"))+
  labs(title = "AEI") + xlab("") + ylab("") +
  theme(plot.title = element_text(family = "serif", size = 26, hjust = .5),
        axis.title = element_text(family = "serif", size = 22),
        axis.text.x = element_text(family = "serif", size = 20),
        axis.text.y = element_text(family = "serif", size = 20),
        axis.line.y = element_line(size = 1), axis.line.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1), axis.ticks.x = element_line(size = 1),
        panel.background = element_blank())


##

#hist((CN.data.local.modelfit$BGNTotal*(-1))); range((CN.data.local.modelfit$BGNTotal*(-1)))
#plot(fitdist((CN.data.local.modelfit$BGNTotal*(-1)), "pois"))  

res1.BGN <- glmmTMB((BGNTotal*(-1)) ~ n_trees + (1|Local:time_period), family=poisson(link = "identity"), CN.data.local.modelfit)
summary(res1.BGN)   

res2.BGN <- glmmTMB((BGNTotal*(-1)) ~ n_vines + (1|Local:time_period), family=poisson(link = "identity"), CN.data.local.modelfit)
summary(res2.BGN)    

res3.BGN <- glmmTMB((BGNTotal*(-1)) ~ n_palmtrees + (1|Local:time_period), family=poisson(link = "identity"), CN.data.local.modelfit)
summary(res3.BGN)     

res4.BGN <- glmmTMB((BGNTotal*(-1)) ~ n_trees + n_palmtrees + (1|Local:time_period), family=poisson(link = "identity"), CN.data.local.modelfit)
summary(res4.BGN)      

res4.BGN.simres <- simulateResiduals(res4.BGN)
plot(res4.BGN.simres)
plot(allEffects(res4.BGN, residuals = T))

res4.BGN.t <- broom.mixed::tidy(res4.BGN, conf.int = TRUE)
res4.BGN.t <- transform(res4.BGN.t, term = sprintf("%s.%s", component, term))
res4.BGN.t$group <- "BGN"
res4.BGN.t$col <- ifelse(res4.BGN.t$conf.low < 0 & res4.BGN.t$conf.high < 0, "#e65a49", 
                         ifelse(res4.BGN.t$conf.low > 0 & res4.BGN.t$conf.high > 0, "#55c0f0", "#dbdfde"))

res4.BGN.graph <- ggplot(res4.BGN.t[1:3, ])+
  geom_point(aes(x = estimate, y = term, col = col), size = 5)+
  geom_segment(aes(x = conf.low, xend = conf.high,
                   y = term, yend = term, col = col), size = 2)+
  geom_vline(xintercept=0, lty=2)+
  scale_color_manual(values = c("#55c0f0", "#dbdfde", "#e65a49"))+
  scale_y_discrete(labels = c("cond.(Intercept)"="Intercept", "cond.n_trees"="Number of trees",
                              "cond.n_palmtrees"="Number of palmtrees"))+
  labs(title = "BGN") + xlab("Coefficient") + ylab("") +
  theme(plot.title = element_text(family = "serif", size = 26, hjust = .5),
        axis.title = element_text(family = "serif", size = 22),
        axis.text.x = element_text(family = "serif", size = 20),
        axis.text.y = element_text(family = "serif", size = 20),
        axis.line.y = element_line(size = 1), axis.line.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1), axis.ticks.x = element_line(size = 1),
        panel.background = element_blank())

##

#hist(CN.data.local.modelfit$HTotal); range(CN.data.local.modelfit$HTotal)
#plot(fitdist(CN.data.local.modelfit$HTotal, "beta"))
res1.H <- glmmTMB(HTotal ~ n_trees + (1|Local:time_period), family=beta_family(), CN.data.local.modelfit)
summary(res1.H)   

res2.H <- glmmTMB(HTotal ~ n_vines + (1|Local:time_period), family=beta_family(), CN.data.local.modelfit)
summary(res2.H)    

res3.H <- glmmTMB(HTotal ~ n_palmtrees + (1|Local:time_period), family=beta_family(), CN.data.local.modelfit)
summary(res3.H)     

res4.H <- glmmTMB(HTotal ~ n_trees + n_palmtrees + (1|Local:time_period), family=beta_family(), CN.data.local.modelfit)
summary(res4.H)      

res4.H.simres <- simulateResiduals(res4.H)
plot(res4.H.simres)
res4.H.alleffects <- allEffects(res4.H, residuals = T)
plot(res4.H.alleffects)

res4.H.t <- broom.mixed::tidy(res4.H, conf.int = TRUE)
res4.H.t <- transform(res4.H.t, term = sprintf("%s.%s", component, term))
res4.H.t$group <- "H"
res4.H.t$col <- ifelse(res4.H.t$conf.low < 0 & res4.H.t$conf.high < 0, "#e65a49", 
                         ifelse(res4.H.t$conf.low > 0 & res4.H.t$conf.high > 0, "#55c0f0", "#dbdfde"))

res4.H.graph <- ggplot(res4.H.t[1:3, ])+
  geom_point(aes(x = estimate, y = term, col = col), size = 5)+
  geom_segment(aes(x = conf.low, xend = conf.high,
                   y = term, yend = term, col = col), size = 2)+
  geom_vline(xintercept=0, lty=2)+
  scale_color_manual(values = c("#55c0f0", "#dbdfde", "#e65a49"))+
  scale_y_discrete(labels = c("cond.(Intercept)"="Intercept", "cond.n_trees"="Number of trees",
                              "cond.n_palmtrees"="Number of palmtrees"))+
  labs(title = "H") + xlab("Coefficient") + ylab("") +
  theme(plot.title = element_text(family = "serif", size = 26, hjust = .5),
        axis.title = element_text(family = "serif", size = 22),
        axis.text.x = element_text(family = "serif", size = 20),
        axis.text.y = element_text(family = "serif", size = 20),
        axis.line.y = element_line(size = 1), axis.line.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1), axis.ticks.x = element_line(size = 1),
        panel.background = element_blank())


##

png("results/localenvironmentaleffects.png", width = 1280, height = 720, units = "px", bg = "transparent")
plot_grid(res4.AA.graph, res4.AEI.graph, res4.BGN.graph, res4.H.graph, nrow=2)
dev.off()

rm(list= ls()[(ls() %in% grep("res", ls(),value = T))])
gc()
#####






