
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

sa <- readOGR(dsn = "C:/Users/miral/Dropbox/GIS/Politico/America do Sul/Continente", layer = "tudo")
sa <- crop(sa, extent(x=c(-83,-17), y=c(-56,15)))
plot(sa, col="gray95")

#amz <- readOGR(dsn = "C:/Users/miral/Dropbox/GIS/Politico/Amazonia bioma", layer = "Lim_Biogeografico")
#plot(amz, add=T, col="gray85")

br <- readOGR(dsn = "C:/Users/miral/Dropbox/GIS/Politico/Divisao Politica BR", layer = "BR_Contorno")
plot(br, add=T, col="gray85")

pa <- readOGR(dsn = "C:/Users/miral/Dropbox/GIS/Politico/Divisao Politica BR", layer = "BRASIL")
pa <- pa[pa@data$UF=="PA",]
plot(pa, add=T, col="gray75")

#ogrListLayers("C:/Users/miral/Dropbox/GIS/Politico/flona_carajas.kml")
carajas <- readOGR(dsn = "C:/Users/miral/Dropbox/GIS/Politico/flona_carajas.kml", layer = "sql_statement")
#plot(carajas, add=T, col="gray55")
#extent(carajas)

#ogrListLayers("C:/Users/miral/Dropbox/GIS/Politico/parna_campos_ferruginosos.kml")
ferrugem <- readOGR(dsn = "C:/Users/miral/Dropbox/GIS/Politico/parna_campos_ferruginosos.kml", layer = "PARNA_Campos_Ferruginosos")
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

# frequency bin 1, dawn
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

# frequency bin 1, day
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

# frequency bin 1, dusk
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

# frequency bin 1, night
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

# frequency bin 1, total
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

# frequency bin 2, dawn
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

# frequency bin 2, day
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

# frequency bin 2, dusk
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

# frequency bin 2, night
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

# frequency bin 2, total
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

# frequency bin 3, dawn
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

# frequency bin 3, day
corrplot(cor(CN.data[CN.data$time_period=="day", grep("fbin3", names(CN.data))], method = "spearman", use = "pairwise.complete.obs"),
         method = "color", order = "alphabet", col = brewer.pal(n = 8, name = "RdBu"), type = "lower",
         mar = c(0, 1, 0, 0), 
         number.font = 3, number.cex = 3,
         addCoef.col = "black",
         tl.pos = "n", cl.pos = "n",
         p.mat = rcorr(as.matrix(CN.data[CN.data$time_period=="day", grep("fbin3", names(CN.data))]), type = "spearman")$P, 
         sig.level = 0.05, insig = "blank",
         diag = F)

# frequency bin 3, dusk
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

# frequency bin 3, night
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

# frequency bin 3, total
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

# frequency bin 4, dawn
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

# frequency bin 4, day
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

# frequency bin 4, dusk
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

# frequency bin 4, night
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

# frequency bin 4, total
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

# total frequency, dawn
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

# total frequency, day
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

# total frequency, dusk
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

# total frequency, night
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
# total frequency, total
png("results/correlogram_total.png", width = 960, height = 960, units = "px", bg = "transparent")

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

# excluding variables: H, S2N, NPeak, AEI
# they were built from other index (e.g., S2N from BGN)
# or they were used to build other index (e.g., NPeak to AA)
# and/or were highly correlated (H and AA and ADI)
CN.data.total <- CN.data

exclude <- c("S2NTotal", "NPeakTotal", "ADITotal",
             "S2Nfbin1", "NPeakfbin1", "ADIfbin1",
             "S2Nfbin2", "NPeakfbin2", "ADIfbin2",
             "S2Nfbin3", "NPeakfbin3", "ADIfbin3",
             "S2Nfbin4", "NPeakfbin4", "ADIfbin4")

CN.data <- CN.data[,!colnames(CN.data) %in% exclude]

#### difference between time period ####
# creating data frame

time.period.diff <- data.frame(index = rep(c("H", "BGN", "AA", "AEI"), 70),
                                frequency.bin = rep(rep(c("Total", "0.3-4", "4-12", "0.3-12", "12-22.1"), each = 4), 14),
                                local = rep(unique(CN.data$Local), each = 20),
                                effect.size = NA, n.sig = NA)

c=1
for (l in unique(CN.data$Local)) {
  
  cnx <- CN.data[CN.data$Local==l,]
  
  for (i in names(CN.data)[13:32]) {
    
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
                       mutate(index=factor(index, levels = c("H", "AEI", "BGN", "AA")),
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

time.period.diff.total <- data.frame(index = rep(c("H", "BGN", "AA", "AEI"), 5),
                                     frequency.bin = rep(c("Total", "0.3-4", "4-12", "0.3-12", "12-22.1"), each = 4),
                                     effect.size = NA, n.sig = NA)

c=1
for (i in names(CN.data)[13:32]) {
  
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
                          mutate(index=factor(index, levels = c("H", "AEI", "BGN", "AA")),
                                 frequency.bin=factor(frequency.bin, levels = c("Total", "0.3-4", "4-12", "0.3-12", "12-22.1")))


time.period.diff.total.hline <- time.period.diff.total[time.period.diff.total$frequency.bin=="Total", ]

png("results/timeperioddifftotal.png", width = 2400, height = 640, units = "px", bg = "transparent")

time.period.diff.total %>%  
  ggplot(aes(frequency.bin, n.sig, colour = effect.size))+
  geom_point(size=12)+
  geom_hline(data=time.period.diff.total.hline, aes(yintercept=n.sig), linewidth=1, linetype='dashed') +
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

time.period.diff.total.2019 <- data.frame(index = rep(c("H", "BGN", "AA", "AEI"), 5),
                                     frequency.bin = rep(c("Total", "0.3-4", "4-12", "0.3-12", "12-22.1"), each = 4),
                                     effect.size = NA, n.sig = NA)

c=1
for (i in names(CN.data.2019)[13:32]) {
  
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
  mutate(index=factor(index, levels = c("H", "AEI", "BGN", "AA")),
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

time.period.diff.total.2022 <- data.frame(index = rep(c("H", "BGN", "AA", "AEI"), 5),
                                          frequency.bin = rep(c("Total", "0.3-4", "4-12", "0.3-12", "12-22.1"), each = 4),
                                          effect.size = NA, n.sig = NA)

c=1
for (i in names(CN.data.2022)[13:32]) {
  
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
  mutate(index=factor(index, levels = c("H", "AEI", "BGN", "AA")),
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

#distribution of values between minutes
cols <- c("CN1" = "#eff1ed", "CN2" = "#373d20", "CN3" = "#717744", "CN4" = "#bcbd8b", "CN5" = "#766153", "CN6" = "#9cc5a1",
          "CN7" = "#49a078", "CN8" = "#216869", "CN9" = "#f4f0bb", "CN10" = "#43291f", "CN11" = "#d36135", "CN12" = "#bc4b51",
          "CN13" = "#427aa1", "CN14" = "#004e64")

png("results/overview_byMin.png", width = 2400, height = 640, units = "px", bg = "transparent")

CN.data %>% dplyr::select("Local", "Min", "HTotal", "BGNTotal", "AATotal", "AEITotal") %>% 
  pivot_longer(cols = HTotal:AEITotal, names_to = "Index", values_to = "Values", values_ptypes = numeric()) %>% 
  mutate(Index = gsub("Total", "", Index)) %>% 
  mutate(Index = factor(Index, levels = c("H", "AEI", "BGN", "AA"))) %>% 
     ggplot(aes(Min, Values, group=Local,  colour = Local)) +
         geom_point(alpha = .25, size = .95, show.legend = F) +
         geom_smooth(linewidth=2, method = "loess", se = F) +
         geom_vline(aes(xintercept=39), linewidth=1, linetype='dashed', color = "gray55") +
         geom_vline(aes(xintercept=51), linewidth=1, linetype='dashed', color = "gray55") +
         geom_vline(aes(xintercept=111), linewidth=1, linetype='dashed', color = "gray55") +
         geom_vline(aes(xintercept=123), linewidth=1, linetype='dashed', color = "gray55") +
         scale_x_discrete(breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+  
         scale_color_manual(values = cols)+
         facet_wrap(Index~., ncol=4, scales = "free") +
         guides(col = guide_legend(nrow = 2, byrow = TRUE))+
         theme_minimal() +
         theme(axis.title = element_text(family = "serif", size = 30),
               axis.text = element_text(family = "serif", size = 28, hjust = 1),
               axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
               axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2),
               strip.text = element_text(family = "serif", size = 28),
               legend.title = element_blank(),
               legend.text = element_text(family = "serif", size = 22),
               legend.position = "bottom", legend.direction = "horizontal", legend.justification = "center")
  
dev.off()
#



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
## preparing data -- regional
CN.data.regional <- CN.data[!is.na(CN.data$precipitation),c(1:16,33:43)]

CN.data.regional$alt <- scale(CN.data.regional$alt, center = F)
CN.data.regional$distwater <- scale(CN.data.regional$distwater, center = F)
CN.data.regional$distedge <- scale(CN.data.regional$distedge, center = F)
CN.data.regional$distcanga <- scale(CN.data.regional$distcanga, center = F)
CN.data.regional$distminning <- scale(CN.data.regional$distminning, center = F)
CN.data.regional$treeheight <- scale(CN.data.regional$treeheight, center = F)
CN.data.regional$NDVI <- scale(CN.data.regional$NDVI, center = F)
CN.data.regional$precipitation <- scale(as.numeric(CN.data.regional$precipitation), center = F)
CN.data.regional$temperature <- scale(as.numeric(CN.data.regional$temperature), center = F)
CN.data.regional$humidity <- scale(as.numeric(CN.data.regional$humidity), center = F)
CN.data.regional$wind <- scale(as.numeric(CN.data.regional$wind), center = F)


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






