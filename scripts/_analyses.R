
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
# excluding NPeak variable
# it was used to calculate acustic activity [AA*]
# total correlation
CN.data <- CN.data[,grep("NPeak", colnames(CN.data), invert = T)]

# excluding S2N variable
# it was calculated base on background noise [BGN*]
# total correlation
CN.data <- CN.data[,grep("S2N", colnames(CN.data), invert = T)]

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
               select(Year:Hour,PRECIPITACAO.TOTAL..HORARIO.mm.:VENTO..VELOCIDADE.HORARIA.m.s.)

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
CN.data <- CN.data[,c(1,8,2,35,36,46,3:7,47,9:33,37:41,49,48,42:45,34)]

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


#### correlogram ####
png("results/correlograms.png", width = 1920, height = 1920, units = "px", bg = "transparent")
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



#### difference between time period ####
# creating data frame

time.period.diff <- data.frame(index = rep(c("H", "BGN", "AA", "ADI", "AEI"), 70),
                                frequency.bin = rep(rep(c("Total", "0.3-4", "4-12", "0.3-12", "12-22.1"), each = 5), 14),
                                local = rep(unique(CN.data$Local), each = 25),
                                effect.size = NA, n.sig = NA)

c=1
for (l in unique(CN.data$Local)) {
  
  cnx <- CN.data[CN.data$Local==l,]
  
  for (i in names(CN.data)[13:37]) {
    
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
                       mutate(index=factor(index, levels = c("AA", "ADI", "AEI", "BGN", "H")),
                              frequency.bin=factor(frequency.bin, levels = c("Total", "0.3-4", "4-12", "0.3-12", "12-22.1")),
                              local=factor(local, levels = paste0("CN", seq(1:14))))
  
  
time.period.diff.hline <- time.period.diff[time.period.diff$frequency.bin=="Total", ]
  
png("results/timeperioddifflocal.png", width = 1440, height = 2560, units = "px", bg = "transparent")

time.period.diff %>%  
  ggplot(aes(frequency.bin, n.sig, colour = effect.size))+
  geom_point(size=12)+
  geom_hline(data=time.period.diff.hline, aes(yintercept=n.sig), size=1, linetype='dashed') +
  scale_y_continuous(limits = c(0, 7), breaks = c(0, 2, 4, 6), expand = c(0.1,0.1)) +
  scale_color_gradient("Effect size", low = "#F4A582", high = "#831529", na.value = NA)+
  guides(size = "none",
         colour = guide_colourbar(barheight = unit(45, "cm"))) +
  ylab(expression(atop("Number of significantly different", "time period pairs")))+xlab("Frequency bins (kHz)")+
  facet_grid(local~index) +
  theme(axis.title = element_text(family = "serif", size = 44),
        axis.text.x = element_text(family = "serif", size = 36, angle = 90, hjust = 1),
        axis.text.y = element_text(family = "serif", size = 36),
        axis.line.y = element_line(size = 1), axis.line.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1), axis.ticks.x = element_line(size = 1),
        legend.title = element_text(family = "serif", size = 36),
        legend.text = element_text(family = "serif", size = 36),
        panel.background = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        panel.spacing.y = unit(2, "lines"),
        strip.text = element_text(family = "serif", size = 36))


dev.off()
#

time.period.diff.total <- data.frame(index = rep(c("H", "BGN", "AA", "ADI", "AEI"), 5),
                                     frequency.bin = rep(c("Total", "0.3-4", "4-12", "0.3-12", "12-22.1"), each = 5),
                                     effect.size = NA, n.sig = NA)

c=1
for (i in names(CN.data)[13:37]) {
  
  cnx.index <- CN.data %>% group_by(Local) %>% sample_n(size = 500, replace = T) %>% ungroup() %>% select(time_period, i)
  
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
                          mutate(index=factor(index, levels = c("AA", "ADI", "AEI", "BGN", "H")),
                                 frequency.bin=factor(frequency.bin, levels = c("Total", "0.3-4", "4-12", "0.3-12", "12-22.1")))


time.period.diff.total.hline <- time.period.diff.total[time.period.diff.total$frequency.bin=="Total", ]

png("results/timeperioddifftotal.png", width = 2400, height = 640, units = "px", bg = "transparent")

time.period.diff.total %>%  
  ggplot(aes(frequency.bin, n.sig, colour = effect.size))+
  geom_point(size=12)+
  geom_hline(data=time.period.diff.total.hline, aes(yintercept=n.sig), size=1, linetype='dashed') +
  scale_y_continuous(limits = c(0, 7), breaks = c(0, 2, 4, 6), expand = c(0.1,0.1)) +
  scale_color_gradient("Effect size", low = "#F4A582", high = "#831529", na.value = NA)+
  guides(size = "none",
         colour = guide_colourbar(barheight = unit(12, "cm"))) +
  ylab(expression(atop("Number of significantly different", "time periods pairs")))+xlab("Frequency bins (kHz)")+
  facet_wrap(~index, nrow = 1) +
  theme(axis.title = element_text(family = "serif", size = 34),
        axis.text.x = element_text(family = "serif", size = 30, angle = 90, hjust = 1),
        axis.text.y = element_text(family = "serif", size = 30),
        axis.line.y = element_line(size = 1), axis.line.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1), axis.ticks.x = element_line(size = 1),
        legend.title = element_text(family = "serif", size = 30),
        legend.text = element_text(family = "serif", size = 30),
        panel.background = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        panel.spacing.y = unit(2, "lines"),
        strip.text = element_text(family = "serif", size = 30))


dev.off()
#
#####


##### modelling the effect of environmental traits on soundscape index ####
#preparing data
CN.data2 <- CN.data[!is.na(CN.data$precipitation),]
CN.data2$alt <- scale(CN.data2$alt, center = F)
CN.data2$distwater <- scale(CN.data2$distwater, center = F)
CN.data2$distedge <- scale(CN.data2$distedge, center = F)
CN.data2$distcanga <- scale(CN.data2$distcanga, center = F)
CN.data2$distminning <- scale(CN.data2$distminning, center = F)
CN.data2$treeheight <- scale(CN.data2$treeheight, center = F)
CN.data2$NDVI <- scale(CN.data2$NDVI, center = F)
CN.data2$precipitation <- scale(as.numeric(CN.data2$precipitation), center = F)
CN.data2$temperature <- scale(as.numeric(CN.data2$temperature), center = F)
CN.data2$humidity <- scale(as.numeric(CN.data2$humidity), center = F)
CN.data2$wind <- scale(as.numeric(CN.data2$wind), center = F)
CN.data2$month_year <- as.factor(paste0(CN.data2$Month, CN.data2$Year))

#hist((CN.data2$BGNTotal*(-1))); range((CN.data2$BGNTotal*(-1)))
#plot(fitdist((CN.data2$BGNTotal*(-1)), "pois"))  

res1.BGN <- glmmTMB((BGNTotal*(-1)) ~ distwater
                  + (1|Local:time_period), 
                  family=poisson(link = "identity"), CN.data2)

summary(res1.BGN)   

res2.BGN <- glmmTMB((BGNTotal*(-1)) ~ distedge
                    + (1|Local:time_period), 
                    family=poisson(link = "identity"), CN.data2)

summary(res2.BGN)    

res3.BGN <- glmmTMB((BGNTotal*(-1)) ~ distcanga
                    + (1|Local:time_period), 
                    family=poisson(link = "identity"), CN.data2)

summary(res3.BGN)     

res4.BGN <- glmmTMB((BGNTotal*(-1)) ~ distminning
                    + (1|Local:time_period), 
                    family=poisson(link = "identity"), CN.data2)

summary(res4.BGN)      

res5.BGN <- glmmTMB((BGNTotal*(-1)) ~ treeheight
                    + (1|Local:time_period), 
                    family=poisson(link = "identity"), CN.data2)

summary(res5.BGN)       

res6.BGN <- glmmTMB((BGNTotal*(-1)) ~ NDVI
                    + (1|Local:time_period), 
                    family=poisson(link = "identity"), CN.data2)

summary(res6.BGN)        

res7.BGN <- glmmTMB((BGNTotal*(-1)) ~ precipitation
                    + (1|Local:time_period), 
                    family=poisson(link = "identity"), CN.data2)

summary(res7.BGN)         

res8.BGN <- glmmTMB((BGNTotal*(-1)) ~ temperature
                    + (1|Local:time_period), 
                    family=poisson(link = "identity"), CN.data2)

summary(res8.BGN)          

res9.BGN <- glmmTMB((BGNTotal*(-1)) ~ humidity
                    + (1|Local:time_period), 
                    family=poisson(link = "identity"), CN.data2)

summary(res9.BGN)           

res10.BGN <- glmmTMB((BGNTotal*(-1)) ~ wind
                    + (1|Local:time_period), 
                    family=poisson(link = "identity"), CN.data2)

summary(res10.BGN)            

res11.BGN <- glmmTMB((BGNTotal*(-1)) ~ distwater * treeheight * NDVI
                     + precipitation + temperature + humidity + wind
                     + distcanga + distminning
                     + (1|Local:time_period), 
                     family=poisson(link = "identity"), CN.data2)

summary(res11.BGN) 

res11.BGN.simres <- simulateResiduals(res11.BGN)
plot(res11.BGN.simres)
plot(allEffects(res11.BGN))

res11.BGN.t <- broom.mixed::tidy(res11.BGN, conf.int = TRUE)
res11.BGN.t <- transform(res11.BGN.t, term = sprintf("%s.%s", component, term))

dwplot(res11.BGN.t, by_2sd=F,
       dot_args = list(size = 5),
       whisker_args = list(size = 2))+
  geom_vline(xintercept=0, lty=2)+
  xlab("Coefficient") + ylab("") +
  theme(axis.title = element_text(family = "serif", size = 22),
        axis.text.x = element_text(family = "serif", size = 20),
        axis.text.y = element_text(family = "serif", size = 20),
        axis.line.y = element_line(size = 1), axis.line.x = element_line(size = 1),
        axis.ticks.y = element_line(size = 1), axis.ticks.x = element_line(size = 1))


##

























































































#Compute summary statistics by groups - count, min, max, mean, sd:
overview <- CN.data %>% 
              group_by(Local, Year) %>%
              summarise(Count = n(),
                        H.Min = min(H, na.rm = TRUE),
                        H.Max = max(H, na.rm = TRUE),
                        H.Mean = mean(H, na.rm = TRUE),
                        H.sd = sd(H, na.rm = TRUE),
                        BGN.Min = min(BGN, na.rm = TRUE),
                        BGN.Max = max(BGN, na.rm = TRUE),
                        BGN.Mean = mean(BGN, na.rm = TRUE),
                        BGN.sd = sd(BGN, na.rm = TRUE),
                        S2N.Min = min(S2N, na.rm = TRUE),
                        S2N.Max = max(S2N, na.rm = TRUE),
                        S2N.Mean = mean(S2N, na.rm = TRUE),
                        S2N.sd = sd(S2N, na.rm = TRUE),
                        NPeak.Min = min(NPeak, na.rm = TRUE),
                        NPeak.Max = max(NPeak, na.rm = TRUE),
                        NPeak.Mean = mean(NPeak, na.rm = TRUE),
                        NPeak.sd = sd(NPeak, na.rm = TRUE),
                        AA.Min = min(AA, na.rm = TRUE),
                        AA.Max = max(AA, na.rm = TRUE),
                        AA.Mean = mean(AA, na.rm = TRUE),
                        AA.sd = sd(AA, na.rm = TRUE),
                        ADI.Min = min(ADI, na.rm = TRUE),
                        ADI.Max = max(ADI, na.rm = TRUE),
                        ADI.Mean = mean(ADI, na.rm = TRUE),
                        ADI.sd = sd(ADI, na.rm = TRUE),
                        AEI.Min = min(AEI, na.rm = TRUE),
                        AEI.Max = max(AEI, na.rm = TRUE),
                        AEI.Mean = mean(AEI, na.rm = TRUE),
                        AEI.sd = sd(AEI, na.rm = TRUE)) %>% 
              ungroup() #%>% 
#pivot_longer(Count:AEI.sd, names_to = "Index summary", values_to = "value")

write.csv(overview, "exploratory/descriptive.csv", row.names = F)


#Visualize your data
ggplot(CN.data, aes(Min, H))+
  geom_boxplot()+
  scale_x_discrete(breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", ""))+
  facet_wrap(~Local, ncol = 2)



#Compute one-way ANOVA test
values.aov<-data.frame(Point=rep(paste0("CN", seq(1:14)), each=7),
                       Index=rep(colnames(CN.data[,c(9:15)]), 14),
                       Df=NA, "Sum sq"=NA, "Mean sq"=NA, "F value"=NA, "Pr(>F)"=NA)

j=1
for(p in unique(CN.data$Point)){
  point_X <- CN.data[which(CN.data$Point==p),]
  
  for(i in 9:15){
    res.aov <- aov(point_X[,i] ~ as.factor(Min), data = point_X)
    values.aov[j,3:7]<-summary(res.aov)[[1]][1,]
    j=j+1
  }
}

write.csv(values.aov, "exploratory/ANOVA_Min.csv", row.names = F)

# Check ANOVA assumptions 
#plot(res.aov, 1) #Homogeneity of variances
#plot(res.aov, 2) #Normality
#shapiro.test(residuals(res.aov ))

# Multiple pairwise-comparison between the means of groups
#res.tukey<-TukeyHSD(res.aov)
#res.tukey<-data.frame(res.tukey$`as.factor(Min)`)
#res.tukey<-res.tukey[which(res.tukey$p.adj<0.05),]


#### between points / minutes & days random effects ####
#plot(fitdist(CN.data$H, "beta"))
res.H <- glmmTMB(H ~ Point + (1|Min) + (1|Day), CN.data, family=beta_family())
summary(res.H)

marginal.H <- emmeans(res.H, ~Point, type = "response")
#pairs(marginal.H, adjust="tukey")
marginal.H <- cld(marginal.H, alpha=0.05, Letters=letters, adjust="tukey")
marginal.H <- marginal.H[order(factor(marginal.H$Point, levels = paste0("CN", seq(1:14)))),]
marginal.H$Point <- ordered(marginal.H$Point, levels = paste0("CN", seq(1:14)))


g1H <- ggplot(CN.data, aes(ordered(CN.data$Point, levels = paste0("CN", seq(1:14))), H))+
  geom_boxplot(outlier.shape = NA, size=2.5)+ 
  annotate("text", x=1:14, y=max(CN.data$H)+0.008, label=gsub(" ","", marginal.H$.group), size=12)+
  coord_cartesian(ylim=c(.75, 1))+ylab("H")+xlab("")+
  theme(axis.title = element_text(family = "serif", size = 48),
        axis.text.x = element_text(family = "serif", size = 48, angle = 45, hjust = 1),
        axis.text.y = element_text(family = "serif", size = 48),
        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2))

g2H <- ggplot(CN.data, aes(factor(Min), H, 
                           group = ordered(CN.data$Point, levels = paste0("CN", seq(1:14))),
                           colour = ordered(CN.data$Point, levels = paste0("CN", seq(1:14)))))+
  stat_summary(fun.y = mean, geom = "line", size=2.5)+ #, show.legend = F
  scale_x_discrete(breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+     
  scale_color_manual(values = pals::glasbey(n=14))+
  coord_cartesian(ylim=c(.75, 1))+ylab("H")+xlab("")+
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme(axis.title = element_text(family = "serif", size = 48),
        axis.text = element_text(family = "serif", size = 48, hjust = 1),
        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2),
        legend.title = element_blank(),
        legend.text = element_text(family = "serif", size = 36),
        legend.position = c(.5,.15), legend.direction = "horizontal", legend.justification = "center")


plot_grid(g1H, g2H, ncol=1, labels = c("A", "B"), label_size=60, label_fontfamily = "serif")



#plot(fitdist(CN.data$BGN*(-1), "pois"))
res.BGN <- glmmTMB(BGN*(-1) ~ Point + (1|Min) + (1|Day), CN.data, family=poisson(link = "identity"))
summary(res.BGN)

marginal.BGN <- emmeans(res.BGN, ~Point, type = "response")
marginal.BGN <- cld(marginal.BGN, alpha=0.05, Letters=letters, adjust="tukey")
marginal.BGN <- marginal.BGN[order(factor(marginal.BGN$Point, levels = paste0("CN", seq(1:14)))),]
marginal.BGN$Point <- ordered(marginal.BGN$Point, levels = paste0("CN", seq(1:14)))


g1BGN <- ggplot(CN.data, aes(ordered(CN.data$Point, levels = paste0("CN", seq(1:14))), BGN))+
  geom_boxplot(outlier.shape = NA, size=2.5)+ 
  annotate("text", x=1:14, y=max(CN.data$BGN)+2, label=gsub(" ","", marginal.BGN$.group), size=12)+
  coord_cartesian(ylim=c(-30, -3))+ylab("BGN")+xlab("")+
  theme(axis.title = element_text(family = "serif", size = 48),
        axis.text.x = element_text(family = "serif", size = 48, angle = 45, hjust = 1),
        axis.text.y = element_text(family = "serif", size = 48),
        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2))

g2BGN <- ggplot(CN.data, aes(factor(Min), BGN, 
                             group = ordered(CN.data$Point, levels = paste0("CN", seq(1:14))),
                             colour = ordered(CN.data$Point, levels = paste0("CN", seq(1:14)))))+
  stat_summary(fun.y = mean, geom = "line", size=2.5)+ #, show.legend = F
  scale_x_discrete(breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+     
  scale_color_manual(values = pals::glasbey(n=14))+
  coord_cartesian(ylim=c(-30, -3))+ylab("BGN")+xlab("")+
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme(axis.title = element_text(family = "serif", size = 48),
        axis.text = element_text(family = "serif", size = 48, hjust = 1),
        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2),
        legend.title = element_blank(),
        legend.text = element_text(family = "serif", size = 36),
        legend.position = c(.5,.15), legend.direction = "horizontal", legend.justification = "center")


plot_grid(g1BGN, g2BGN, ncol=1, labels = c("A", "B"), label_size=60, label_fontfamily = "serif")



#plot(fitdist(CN.data$S2N, "pois"))
res.S2N <- glmmTMB(S2N ~ Point + (1|Min) + (1|Day), CN.data, family=poisson(link = "identity"))
summary(res.S2N)

marginal.S2N <- emmeans(res.S2N, ~Point, type = "response")
marginal.S2N <- cld(marginal.S2N, alpha=0.05, Letters=letters, adjust="tukey")
marginal.S2N <- marginal.S2N[order(factor(marginal.S2N$Point, levels = paste0("CN", seq(1:14)))),]
marginal.S2N$Point <- ordered(marginal.S2N$Point, levels = paste0("CN", seq(1:14)))


g1S2N <- ggplot(CN.data, aes(ordered(CN.data$Point, levels = paste0("CN", seq(1:14))), S2N))+
  geom_boxplot(outlier.shape = NA, size=2.5)+ 
  annotate("text", x=1:14, y=max(CN.data$S2N)+2, label=gsub(" ","", marginal.S2N$.group), size=12)+
  coord_cartesian(ylim=c(0, 30))+ylab("S2N")+xlab("")+
  theme(axis.title = element_text(family = "serif", size = 48),
        axis.text.x = element_text(family = "serif", size = 48, angle = 45, hjust = 1),
        axis.text.y = element_text(family = "serif", size = 48),
        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2))

g2S2N <- ggplot(CN.data, aes(factor(Min), S2N, 
                             group = ordered(CN.data$Point, levels = paste0("CN", seq(1:14))),
                             colour = ordered(CN.data$Point, levels = paste0("CN", seq(1:14)))))+
  stat_summary(fun.y = mean, geom = "line", size=2.5)+ #, show.legend = F
  scale_x_discrete(breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+     
  scale_color_manual(values = pals::glasbey(n=14))+
  coord_cartesian(ylim=c(0, 30))+ylab("S2N")+xlab("")+
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme(axis.title = element_text(family = "serif", size = 48),
        axis.text = element_text(family = "serif", size = 48, hjust = 1),
        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2),
        legend.title = element_blank(),
        legend.text = element_text(family = "serif", size = 36),
        legend.position = c(.5,.15), legend.direction = "horizontal", legend.justification = "center")


plot_grid(g1S2N, g2S2N, ncol=1, labels = c("A", "B"), label_size=60, label_fontfamily = "serif")



#plot(fitdist(CN.data$NPeak, "pois"))
res.NPeak <- glmmTMB(NPeak ~ Point + (1|Min) + (1|Day), CN.data, family=poisson(link = "log"))
summary(res.NPeak)

marginal.NPeak <- emmeans(res.NPeak, ~Point, type = "response")
marginal.NPeak <- cld(marginal.NPeak, alpha=0.05, Letters=letters, adjust="tukey")
marginal.NPeak <- marginal.NPeak[order(factor(marginal.NPeak$Point, levels = paste0("CN", seq(1:14)))),]
marginal.NPeak$Point <- ordered(marginal.NPeak$Point, levels = paste0("CN", seq(1:14)))


g1NPeak <- ggplot(CN.data, aes(ordered(CN.data$Point, levels = paste0("CN", seq(1:14))), NPeak))+
  geom_boxplot(outlier.shape = NA, size=2.5)+ 
  annotate("text", x=1:14, y=450000, label=gsub(" ","", marginal.NPeak$.group), size=12)+
  scale_y_continuous(labels = scales::comma)+         
  coord_cartesian(ylim=c(0, 500000))+ylab("N. Peak")+xlab("")+
  theme(axis.title = element_text(family = "serif", size = 48),
        axis.text.x = element_text(family = "serif", size = 48, angle = 45, hjust = 1),
        axis.text.y = element_text(family = "serif", size = 48),
        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2))

g2NPeak <- ggplot(CN.data, aes(factor(Min), NPeak, 
                               group = ordered(CN.data$Point, levels = paste0("CN", seq(1:14))),
                               colour = ordered(CN.data$Point, levels = paste0("CN", seq(1:14)))))+
  stat_summary(fun.y = mean, geom = "line", size=2.5)+ #, show.legend = F
  scale_x_discrete(breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+     
  scale_color_manual(values = pals::glasbey(n=14))+
  scale_y_continuous(labels = scales::comma)+
  coord_cartesian(ylim=c(0, 500000))+ylab("N. Peak")+xlab("")+
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme(axis.title = element_text(family = "serif", size = 48),
        axis.text = element_text(family = "serif", size = 48, hjust = 1),
        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2),
        legend.title = element_blank(),
        legend.text = element_text(family = "serif", size = 36),
        legend.position = c(.5,.88), legend.direction = "horizontal", legend.justification = "center")


plot_grid(g1NPeak, g2NPeak, ncol=1, labels = c("A", "B"), label_size=60, label_fontfamily = "serif")



n<-length(CN.data$AA)
CN.data$AA <- (CN.data$AA *(n -1) +0.5)/n
#plot(fitdist(CN.data$AA, "beta"))
res.AA <- glmmTMB(AA ~ Point + (1|Min) + (1|Day), CN.data, family=beta_family())
summary(res.AA)

marginal.AA <- emmeans(res.AA, ~Point, type = "response")
marginal.AA <- cld(marginal.AA, alpha=0.05, Letters=letters, adjust="tukey")
marginal.AA <- marginal.AA[order(factor(marginal.AA$Point, levels = paste0("CN", seq(1:14)))),]
marginal.AA$Point <- ordered(marginal.AA$Point, levels = paste0("CN", seq(1:14)))


g1AA <- ggplot(CN.data, aes(ordered(CN.data$Point, levels = paste0("CN", seq(1:14))), AA))+
  geom_boxplot(outlier.shape = NA, size=2.5)+ 
  annotate("text", x=1:14, y=.45, label=gsub(" ","", marginal.AA$.group), size=12)+
  coord_cartesian(ylim=c(0, .5))+ylab("Acoustic Activity")+xlab("")+
  theme(axis.title = element_text(family = "serif", size = 48),
        axis.text.x = element_text(family = "serif", size = 48, angle = 45, hjust = 1),
        axis.text.y = element_text(family = "serif", size = 48),
        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2))

g2AA <- ggplot(CN.data, aes(factor(Min), AA, 
                            group = ordered(CN.data$Point, levels = paste0("CN", seq(1:14))),
                            colour = ordered(CN.data$Point, levels = paste0("CN", seq(1:14)))))+
  stat_summary(fun.y = mean, geom = "line", size=2.5)+ #, show.legend = F
  scale_x_discrete(breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+     
  scale_color_manual(values = pals::glasbey(n=14))+
  coord_cartesian(ylim=c(0, .5))+ylab("Acoustic Activity")+xlab("")+
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme(axis.title = element_text(family = "serif", size = 48),
        axis.text = element_text(family = "serif", size = 48, hjust = 1),
        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2),
        legend.title = element_blank(),
        legend.text = element_text(family = "serif", size = 36),
        legend.position = c(.5,.88), legend.direction = "horizontal", legend.justification = "center")


plot_grid(g1AA, g2AA, ncol=1, labels = c("A", "B"), label_size=60, label_fontfamily = "serif")



#plot(fitdist(CN.data$ADI, "norm"))
res.ADI <- glmmTMB(ADI ~ Point + (1|Min) + (1|Day), CN.data, family=gaussian())
summary(res.ADI)

marginal.ADI <- emmeans(res.ADI, ~Point, type = "response")
marginal.ADI <- cld(marginal.ADI, alpha=0.05, Letters=letters, adjust="tukey")
marginal.ADI <- marginal.ADI[order(factor(marginal.ADI$Point, levels = paste0("CN", seq(1:14)))),]
marginal.ADI$Point <- ordered(marginal.ADI$Point, levels = paste0("CN", seq(1:14)))


g1ADI <- ggplot(CN.data, aes(ordered(CN.data$Point, levels = paste0("CN", seq(1:14))), ADI))+
  geom_boxplot(outlier.shape = NA, size=2.5)+ 
  annotate("text", x=1:14, y=max(CN.data$ADI)+0.1, label=gsub(" ","", marginal.ADI$.group), size=12)+
  coord_cartesian(ylim=c(0, 3))+ylab("Acoustic Diversity")+xlab("")+
  theme(axis.title = element_text(family = "serif", size = 48),
        axis.text.x = element_text(family = "serif", size = 48, angle = 45, hjust = 1),
        axis.text.y = element_text(family = "serif", size = 48),
        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2))

g2ADI <- ggplot(CN.data, aes(factor(Min), ADI, 
                             group = ordered(CN.data$Point, levels = paste0("CN", seq(1:14))),
                             colour = ordered(CN.data$Point, levels = paste0("CN", seq(1:14)))))+
  stat_summary(fun.y = mean, geom = "line", size=2.5)+ #, show.legend = F
  scale_x_discrete(breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+     
  scale_color_manual(values = pals::glasbey(n=14))+
  coord_cartesian(ylim=c(0, 3))+ylab("Acoustic Diversity")+xlab("")+
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme(axis.title = element_text(family = "serif", size = 48),
        axis.text = element_text(family = "serif", size = 48, hjust = 1),
        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2),
        legend.title = element_blank(),
        legend.text = element_text(family = "serif", size = 36),
        legend.position = c(.5,.88), legend.direction = "horizontal", legend.justification = "center")


plot_grid(g1ADI, g2ADI, ncol=1, labels = c("A", "B"), label_size=60, label_fontfamily = "serif")



#plot(fitdist(CN.data$AEI[complete.cases(CN.data$AEI)], "beta"))
res.AEI <- glmmTMB(AEI ~ Point + (1|Min) + (1|Day), CN.data[complete.cases(CN.data$AEI),], family=beta_family(link = "logit"))
summary(res.AEI)

marginal.AEI <- emmeans(res.AEI, ~Point, type = "response")
marginal.AEI <- cld(marginal.AEI, alpha=0.05, Letters=letters, adjust="tukey")
marginal.AEI <- marginal.AEI[order(factor(marginal.AEI$Point, levels = paste0("CN", seq(1:14)))),]
marginal.AEI$Point <- ordered(marginal.AEI$Point, levels = paste0("CN", seq(1:14)))


g1AEI <- ggplot(CN.data, aes(ordered(CN.data$Point, levels = paste0("CN", seq(1:14))), AEI))+
  geom_boxplot(outlier.shape = NA, size=2.5)+ 
  annotate("text", x=1:14, y=max(CN.data$AEI, na.rm = T)+0.1, label=gsub(" ","", marginal.AEI$.group), size=12)+
  coord_cartesian(ylim=c(0, 1.1))+ylab("Acoustic Evenness")+xlab("")+
  theme(axis.title = element_text(family = "serif", size = 48),
        axis.text.x = element_text(family = "serif", size = 48, angle = 45, hjust = 1),
        axis.text.y = element_text(family = "serif", size = 48),
        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2))

g2AEI <- ggplot(CN.data, aes(factor(Min), AEI, 
                             group = ordered(CN.data$Point, levels = paste0("CN", seq(1:14))),
                             colour = ordered(CN.data$Point, levels = paste0("CN", seq(1:14)))))+
  stat_summary(fun.y = mean, geom = "line", size=2.5)+ #, show.legend = F
  scale_x_discrete(breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", "23:00"))+     
  scale_color_manual(values = pals::glasbey(n=14))+
  coord_cartesian(ylim=c(0, 1.1))+ylab("Acoustic Evenness")+xlab("")+
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme(axis.title = element_text(family = "serif", size = 48),
        axis.text = element_text(family = "serif", size = 48, hjust = 1),
        axis.line.y = element_line(size = 2), axis.line.x = element_line(size = 2),
        axis.ticks.y = element_line(size = 2), axis.ticks.x = element_line(size = 2),
        legend.title = element_blank(),
        legend.text = element_text(family = "serif", size = 36),
        legend.position = c(.5,.15), legend.direction = "horizontal", legend.justification = "center")


plot_grid(g1AEI, g2AEI, ncol=1, labels = c("A", "B"), label_size=60, label_fontfamily = "serif")



#



































































fp.prop1C<-data.frame(ID=rep("SDE1", 6144), Time=rep(seq(1:24), each=256), Freq = rep(unique(fp.f$Freq), 24), z=NA)
#fp.prop2C<-data.frame(ID=rep("SDE2", 6144), Time=rep(seq(1:24), each=256), Freq = rep(unique(fp.f$Freq), 24), z=NA)


cont<-24
j=1
for (m in 1:24) {
  for (z in unique(fp.f$Freq)) {
    fp.prop1C$z[which(fp.prop1C$ID==grep("SDE1", fp.prop1C$ID, value = T) & 
                        fp.prop1C$Time==m & fp.prop1C$Freq==z)] <- sum(fp.f$Peak[which(fp.f$ID==grep("SED1", fp.f$ID, value = T) & 
                                                                                         fp.f$Min==j & 
                                                                                         fp.f$Freq==z)], 
                                                                       fp.f$Peak[which(fp.f$ID==grep("SED1", fp.f$ID, value = T) & 
                                                                                         fp.f$Min==j+1 & 
                                                                                         fp.f$Freq==z)],
                                                                       fp.f$Peak[which(fp.f$ID==grep("SED1", fp.f$ID, value = T) & 
                                                                                         fp.f$Min==j+2 & 
                                                                                         fp.f$Freq==z)], 
                                                                       fp.f$Peak[which(fp.f$ID==grep("SED1", fp.f$ID, value = T) & 
                                                                                         fp.f$Min==j+3 & 
                                                                                         fp.f$Freq==z)],
                                                                       fp.f$Peak[which(fp.f$ID==grep("SED1", fp.f$ID, value = T) & 
                                                                                         fp.f$Min==j+4 & 
                                                                                         fp.f$Freq==z)], 
                                                                       fp.f$Peak[which(fp.f$ID==grep("SED1", fp.f$ID, value = T) & 
                                                                                         fp.f$Min==j+5 & 
                                                                                         fp.f$Freq==z)], na.rm = T)/(length(fp.f$Peak[which(fp.f$ID==grep("SED1", fp.f$ID, value = T) & 
                                                                                                                                              fp.f$Min==j & 
                                                                                                                                              fp.f$Freq==z)]) + 
                                                                                                                       length(fp.f$Peak[which(fp.f$ID==grep("SED1", fp.f$ID, value = T) & 
                                                                                                                                                fp.f$Min==j+1 & 
                                                                                                                                                fp.f$Freq==z)]) +
                                                                                                                       length(fp.f$Peak[which(fp.f$ID==grep("SED1", fp.f$ID, value = T) & 
                                                                                                                                                fp.f$Min==j+2 & 
                                                                                                                                                fp.f$Freq==z)]) +
                                                                                                                       length(fp.f$Peak[which(fp.f$ID==grep("SED1", fp.f$ID, value = T) & 
                                                                                                                                                fp.f$Min==j+3 & 
                                                                                                                                                fp.f$Freq==z)]) +
                                                                                                                       length(fp.f$Peak[which(fp.f$ID==grep("SED1", fp.f$ID, value = T) & 
                                                                                                                                                fp.f$Min==j+4 & 
                                                                                                                                                fp.f$Freq==z)]) +
                                                                                                                       length(fp.f$Peak[which(fp.f$ID==grep("SED1", fp.f$ID, value = T) & 
                                                                                                                                                fp.f$Min==j+5 & 
                                                                                                                                                fp.f$Freq==z)]))
    
  }
  j=j+6
  cont=cont-1
  cat('\n>faltam', cont, 'horas <\n')
  
}

fp.prop<-rbind(fp.prop1C, fp.prop2C)

write.csv(fp.prop, "piloto_201908/prop_frequency_peak_byH.csv", row.names = F)


library(tidyverse)
library(ggplot2)
library(rayshader)

asu1 = ggplot(fp.prop[c(1:6144),], aes(Time, Freq, fill = z)) +
  geom_raster(aes(group=Freq), interpolate = TRUE, show.legend = F) +
  scale_x_continuous(expand=c(0,0),breaks=seq(0,24,2)) +
  scale_y_continuous(expand=c(0,0),breaks=seq(0,10,2.5),limits=c(0,11)) +
  scale_fill_gradientn(colours=c("blue", "green", "yellow", "red")) +
  labs(tag = "SDE1", x="Hour", y= "Frequency (kHz)") + 
  theme(text = element_text(family = "serif", size = 42),
        axis.title.y = element_text(angle=270),
        plot.tag.position = c(0.9, 0.9), 
        plot.tag = element_text(color = "white", face = "bold"))

asu1



plot_gg(asu1, multicore=TRUE, height = 11, width = 15, scale = 75, zoom=0.9, phi=30, theta=-30)
#render_snapshot()
#render_camera(zoom="resolu??o",phi=inclina??o,theta=angulo)


plot_gg(asu1, multicore=TRUE, height = 13, width = 17, scale = 300, zoom=0.57, phi=40, theta=-43, windowsize=c(1900,1024))

pdf("piloto_201908/asu1.pdf", height = 13, width = 17, bg = "transparent")  
render_snapshot()
dev.off()

#














i=1
for (i in 1:9) {
  multiple_sounds(directory = l[i], resultfile = paste("piloto_201908/adi_", n[i], ".csv", sep=""), 
                  soundindex = "acoustic_diversity", max_freq = 22500, db_threshold = -50, 
                  no_cores = 4)
  
  multiple_sounds(directory = l[i], resultfile = paste("piloto_201908/aei_", n[i], ".csv", sep=""), 
                  soundindex = "acoustic_evenness", max_freq = 22500, db_threshold = -50, 
                  no_cores = 4)
  
  multiple_sounds(directory = l[i], resultfile = paste("piloto_201908/bio_", n[i], ".csv", sep=""), 
                  soundindex = "bioacoustic_index", max_freq = 22500, fft_w = 512, 
                  no_cores = 4)
}










































############################
#                          #
#  background noise (BGN)  #
#       power (POW)        #
#                          #
############################



cont <- length(l)
#i=1
for (i in 1:length(l)) {
  
  ss = spectro(normalize(downsample(readWave(l[i]), 44100), unit = "1"), wl = 512, plot = F)
  mm =
    melt(ss$amp, value.name = "Amplitude") %>%
    dplyr::select(FrequencyIndex = Var1, TimeIndex = Var2, Amplitude)
  ff =
    melt(ss$freq, value.name = "Frequency") %>%
    dplyr::mutate(FrequencyIndex = row_number(), Frequency = Frequency)
  tt =
    melt(ss$time, value.name = "Time") %>%
    dplyr::mutate(TimeIndex = row_number())
  sp =
    mm %>%
    dplyr::left_join(ff, by = "FrequencyIndex") %>%
    dplyr::left_join(tt, by = "TimeIndex") %>%
    dplyr::select(Time, Frequency, Amplitude)
  
  sp$Amplitude<-round(sp$Amplitude, 0)
  
  index1 <- data.frame(file=n[i], Frequency=unique(sp$Frequency), BGN=NA, POW=NA)
  
  for (z in unique(sp$Frequency)) {
    
    index1$BGN[which(index1$Frequency==z)] <- getmode(sp$Amplitude[which(sp$Frequency==z)])
    index1$POW[which(index1$Frequency==z)] <- max(sp$Amplitude[which(sp$Frequency==z)]) - index1$BGN[which(index1$Frequency==z)]
    
  }  
  
  write.table(index1, "piloto_201908/SDE2C_BGN_POW_total.csv", append = T, row.names = F, col.names = F)
  cont<-cont-1
  cat('\n>finalizado', n[i], 'faltam', cont, '<\n')  
}



SDE2C_index <- read.csv("piloto_201908/SDE2C_BGN_POW_total.csv", sep = " ", header = F)
colnames(SDE2C_index)<-c("File", "Frequency", "BGN", "POW")
SDE2C_index$Dia <- rep(c(1:3), each=36864)
SDE2C_index$Min <- rep(c(1:144), each=256)
SDE2C_index<-SDE2C_index[,c(1,5,6,2,3,4)]
write.csv(SDE2C_index, "piloto_201908/SDE2C_BGN_POW_total.csv", row.names = F)


SDE2C_index_vf <- SDE2C_index[c(1:36864),-2]
SDE2C_index_vf[,c(4,5)]<- NA
SDE2C_index_vf$File<-grep("SED2C", unlist(strsplit(as.character(unique(SDE2C_index_vf$File)), "_")), value = T)


for (ii in unique(SDE2C_index$Min)) {
  min_da_vez <- SDE2C_index[which(SDE2C_index$Min==ii),]
  
  for (zz in unique(min_da_vez$Frequency)) {
    SDE2C_index_vf$BGN[which(SDE2C_index_vf$Min==ii & SDE2C_index_vf$Frequency==zz)] <- round(mean(min_da_vez$BGN[which(min_da_vez$Frequency==zz)]),0)
    SDE2C_index_vf$POW[which(SDE2C_index_vf$Min==ii & SDE2C_index_vf$Frequency==zz)] <- round(mean(min_da_vez$POW[which(min_da_vez$Frequency==zz)]),0)
  }
}

write.csv(SDE1C_index_vf, "piloto_201908/SDE2C_BGN_POW_mean.csv", row.names = F)

rm(list= ls()[(ls() %in% c("l", "n","i", "ii", "z", "zz", "cont", "ss", "mm", "ff", "tt", "sp", "min_da_vez"))])   






index.f<-rbind(SDE1C_index_vf, SDE2C_index_vf)

summary(index.f)
hb<-hist(index.f$BGN, breaks = 20)
hp<-hist(index.f$POW)

for (y in 1:length(index.f$File)) {
  if (index.f$BGN[y] > -29) { index.f$aa1[y]<-1 } else if (index.f$POW[y] > 10) { index.f$aa1[y]<-1 } else { index.f$aa1[y]<-0 }
  if (index.f$BGN[y] > -28) { index.f$aa2[y]<-1 } else if (index.f$POW[y] > 10) { index.f$aa2[y]<-1 } else { index.f$aa2[y]<-0 }
  if (index.f$BGN[y] > -26) { index.f$aa3[y]<-1 } else if (index.f$POW[y] > 10) { index.f$aa3[y]<-1 } else { index.f$aa3[y]<-0 }
  if (index.f$BGN[y] > -24) { index.f$aa4[y]<-1 } else if (index.f$POW[y] > 10) { index.f$aa4[y]<-1 } else { index.f$aa4[y]<-0 }
  if (index.f$BGN[y] > -20) { index.f$aa5[y]<-1 } else if (index.f$POW[y] > 10) { index.f$aa5[y]<-1 } else { index.f$aa5[y]<-0 }
  if (index.f$BGN[y] > -29) { index.f$aa6[y]<-1 } else if (index.f$POW[y] > 12) { index.f$aa6[y]<-1 } else { index.f$aa6[y]<-0 }
  if (index.f$BGN[y] > -28) { index.f$aa7[y]<-1 } else if (index.f$POW[y] > 12) { index.f$aa7[y]<-1 } else { index.f$aa7[y]<-0 }
  if (index.f$BGN[y] > -26) { index.f$aa8[y]<-1 } else if (index.f$POW[y] > 12) { index.f$aa8[y]<-1 } else { index.f$aa8[y]<-0 }
  if (index.f$BGN[y] > -24) { index.f$aa9[y]<-1 } else if (index.f$POW[y] > 12) { index.f$aa9[y]<-1 } else { index.f$aa9[y]<-0 }
  if (index.f$BGN[y] > -20) { index.f$aa10[y]<-1 } else if (index.f$POW[y] > 12) { index.f$aa10[y]<-1 } else { index.f$aa10[y]<-0 }
  if (index.f$BGN[y] > -29) { index.f$aa11[y]<-1 } else if (index.f$POW[y] > 14) { index.f$aa11[y]<-1 } else { index.f$aa11[y]<-0 }
  if (index.f$BGN[y] > -28) { index.f$aa12[y]<-1 } else if (index.f$POW[y] > 14) { index.f$aa12[y]<-1 } else { index.f$aa12[y]<-0 }
  if (index.f$BGN[y] > -26) { index.f$aa13[y]<-1 } else if (index.f$POW[y] > 14) { index.f$aa13[y]<-1 } else { index.f$aa13[y]<-0 }
  if (index.f$BGN[y] > -24) { index.f$aa14[y]<-1 } else if (index.f$POW[y] > 14) { index.f$aa14[y]<-1 } else { index.f$aa14[y]<-0 }
  if (index.f$BGN[y] > -20) { index.f$aa15[y]<-1 } else if (index.f$POW[y] > 14) { index.f$aa15[y]<-1 } else { index.f$aa15[y]<-0 }
  if (index.f$BGN[y] > -29) { index.f$aa16[y]<-1 } else if (index.f$POW[y] > 16) { index.f$aa16[y]<-1 } else { index.f$aa16[y]<-0 }
  if (index.f$BGN[y] > -28) { index.f$aa17[y]<-1 } else if (index.f$POW[y] > 16) { index.f$aa17[y]<-1 } else { index.f$aa17[y]<-0 }
  if (index.f$BGN[y] > -26) { index.f$aa18[y]<-1 } else if (index.f$POW[y] > 16) { index.f$aa18[y]<-1 } else { index.f$aa18[y]<-0 }
  if (index.f$BGN[y] > -24) { index.f$aa19[y]<-1 } else if (index.f$POW[y] > 16) { index.f$aa19[y]<-1 } else { index.f$aa19[y]<-0 }
  if (index.f$BGN[y] > -20) { index.f$aa20[y]<-1 } else if (index.f$POW[y] > 16) { index.f$aa20[y]<-1 } else { index.f$aa20[y]<-0 }
  if (index.f$BGN[y] > -29) { index.f$aa21[y]<-1 } else if (index.f$POW[y] > 18) { index.f$aa21[y]<-1 } else { index.f$aa21[y]<-0 }
  if (index.f$BGN[y] > -28) { index.f$aa22[y]<-1 } else if (index.f$POW[y] > 18) { index.f$aa22[y]<-1 } else { index.f$aa22[y]<-0 }
  if (index.f$BGN[y] > -26) { index.f$aa23[y]<-1 } else if (index.f$POW[y] > 18) { index.f$aa23[y]<-1 } else { index.f$aa23[y]<-0 }
  if (index.f$BGN[y] > -24) { index.f$aa24[y]<-1 } else if (index.f$POW[y] > 18) { index.f$aa24[y]<-1 } else { index.f$aa24[y]<-0 }
  if (index.f$BGN[y] > -20) { index.f$aa25[y]<-1 } else if (index.f$POW[y] > 18) { index.f$aa25[y]<-1 } else { index.f$aa25[y]<-0 }
  if (index.f$BGN[y] > -29) { index.f$aa26[y]<-1 } else if (index.f$POW[y] > 20) { index.f$aa26[y]<-1 } else { index.f$aa26[y]<-0 }
  if (index.f$BGN[y] > -28) { index.f$aa27[y]<-1 } else if (index.f$POW[y] > 20) { index.f$aa27[y]<-1 } else { index.f$aa27[y]<-0 }
  if (index.f$BGN[y] > -26) { index.f$aa28[y]<-1 } else if (index.f$POW[y] > 20) { index.f$aa28[y]<-1 } else { index.f$aa28[y]<-0 }
  if (index.f$BGN[y] > -24) { index.f$aa29[y]<-1 } else if (index.f$POW[y] > 20) { index.f$aa29[y]<-1 } else { index.f$aa29[y]<-0 }
  if (index.f$BGN[y] > -20) { index.f$aa30[y]<-1 } else if (index.f$POW[y] > 20) { index.f$aa30[y]<-1 } else { index.f$aa30[y]<-0 }
  if (index.f$BGN[y] > -29) { index.f$aa31[y]<-1 } else if (index.f$POW[y] > 22) { index.f$aa31[y]<-1 } else { index.f$aa31[y]<-0 }
  if (index.f$BGN[y] > -28) { index.f$aa32[y]<-1 } else if (index.f$POW[y] > 22) { index.f$aa32[y]<-1 } else { index.f$aa32[y]<-0 }
  if (index.f$BGN[y] > -26) { index.f$aa33[y]<-1 } else if (index.f$POW[y] > 22) { index.f$aa33[y]<-1 } else { index.f$aa33[y]<-0 }
  if (index.f$BGN[y] > -24) { index.f$aa34[y]<-1 } else if (index.f$POW[y] > 22) { index.f$aa34[y]<-1 } else { index.f$aa34[y]<-0 }
  if (index.f$BGN[y] > -20) { index.f$aa35[y]<-1 } else if (index.f$POW[y] > 22) { index.f$aa35[y]<-1 } else { index.f$aa35[y]<-0 }
  if (index.f$BGN[y] > -29) { index.f$aa36[y]<-1 } else if (index.f$POW[y] > 24) { index.f$aa36[y]<-1 } else { index.f$aa36[y]<-0 }
  if (index.f$BGN[y] > -28) { index.f$aa37[y]<-1 } else if (index.f$POW[y] > 24) { index.f$aa37[y]<-1 } else { index.f$aa37[y]<-0 }
  if (index.f$BGN[y] > -26) { index.f$aa38[y]<-1 } else if (index.f$POW[y] > 24) { index.f$aa38[y]<-1 } else { index.f$aa38[y]<-0 }
  if (index.f$BGN[y] > -24) { index.f$aa39[y]<-1 } else if (index.f$POW[y] > 24) { index.f$aa39[y]<-1 } else { index.f$aa39[y]<-0 }
  if (index.f$BGN[y] > -20) { index.f$aa40[y]<-1 } else if (index.f$POW[y] > 24) { index.f$aa40[y]<-1 } else { index.f$aa40[y]<-0 }
  if (index.f$BGN[y] > -29) { index.f$aa41[y]<-1 } else if (index.f$POW[y] > 26) { index.f$aa41[y]<-1 } else { index.f$aa41[y]<-0 }
  if (index.f$BGN[y] > -28) { index.f$aa42[y]<-1 } else if (index.f$POW[y] > 26) { index.f$aa42[y]<-1 } else { index.f$aa42[y]<-0 }
  if (index.f$BGN[y] > -26) { index.f$aa43[y]<-1 } else if (index.f$POW[y] > 26) { index.f$aa43[y]<-1 } else { index.f$aa43[y]<-0 }
  if (index.f$BGN[y] > -24) { index.f$aa44[y]<-1 } else if (index.f$POW[y] > 26) { index.f$aa44[y]<-1 } else { index.f$aa44[y]<-0 }
  if (index.f$BGN[y] > -20) { index.f$aa45[y]<-1 } else if (index.f$POW[y] > 26) { index.f$aa45[y]<-1 } else { index.f$aa45[y]<-0 }
  if (index.f$BGN[y] > -29) { index.f$aa46[y]<-1 } else if (index.f$POW[y] > 28) { index.f$aa46[y]<-1 } else { index.f$aa46[y]<-0 }
  if (index.f$BGN[y] > -28) { index.f$aa47[y]<-1 } else if (index.f$POW[y] > 28) { index.f$aa47[y]<-1 } else { index.f$aa47[y]<-0 }
  if (index.f$BGN[y] > -26) { index.f$aa48[y]<-1 } else if (index.f$POW[y] > 28) { index.f$aa48[y]<-1 } else { index.f$aa48[y]<-0 }
  if (index.f$BGN[y] > -24) { index.f$aa49[y]<-1 } else if (index.f$POW[y] > 28) { index.f$aa49[y]<-1 } else { index.f$aa49[y]<-0 }
  if (index.f$BGN[y] > -20) { index.f$aa50[y]<-1 } else if (index.f$POW[y] > 28) { index.f$aa50[y]<-1 } else { index.f$aa50[y]<-0 }  
  if (index.f$BGN[y] > -29) { index.f$aa51[y]<-1 } else if (index.f$POW[y] > 30) { index.f$aa51[y]<-1 } else { index.f$aa51[y]<-0 }
  if (index.f$BGN[y] > -28) { index.f$aa52[y]<-1 } else if (index.f$POW[y] > 30) { index.f$aa52[y]<-1 } else { index.f$aa52[y]<-0 }
  if (index.f$BGN[y] > -26) { index.f$aa53[y]<-1 } else if (index.f$POW[y] > 30) { index.f$aa53[y]<-1 } else { index.f$aa53[y]<-0 }
  if (index.f$BGN[y] > -24) { index.f$aa54[y]<-1 } else if (index.f$POW[y] > 30) { index.f$aa54[y]<-1 } else { index.f$aa54[y]<-0 }
  if (index.f$BGN[y] > -20) { index.f$aa55[y]<-1 } else if (index.f$POW[y] > 30) { index.f$aa55[y]<-1 } else { index.f$aa55[y]<-0 }
}

write.csv(index.f, "piloto_201908/acoustically_active_min.csv", row.names = F)

#aa1 == BGN > -29 (90%) ; POW > 10 *159
#aa2 == BGN > -28 (92%) ; POW > 10 *175
#aa3 == BGN > -26 (95%) ; POW > 10 *233
#aa4 == BGN > -24 (97%) ; POW > 10 *282
#aa5 == BGN > -20 (99%) ; POW > 10 *394

#aa6 == BGN > -29 (90%) ; POW > 12 *1382
#aa7 == BGN > -28 (92%) ; POW > 12 *1550
#aa8 == BGN > -26 (95%) ; POW > 12 *1868
#aa9 == BGN > -24 (97%) ; POW > 12 *2189
#aa10 == BGN > -20 (99%) ; POW > 12 *2800

#aa11 == BGN > -29 (90%) ; POW > 14 *5409
#aa12 == BGN > -28 (92%) ; POW > 14 *5862
#aa13 == BGN > -26 (95%) ; POW > 14 *6676
#aa14 == BGN > -24 (97%) ; POW > 14 *7397
#aa15 == BGN > -20 (99%) ; POW > 14 *8580

#aa16 == BGN > -29 (90%) ; POW > 16 *12904
#aa17 == BGN > -28 (92%) ; POW > 16 *13759
#aa18 == BGN > -26 (95%) ; POW > 16 *15139
#aa19 == BGN > -24 (97%) ; POW > 16 *16170
#aa20 == BGN > -20 (99%) ; POW > 16 *17651

#aa21 == BGN > -29 (90%) ; POW > 18 *22676
#aa22 == BGN > -28 (92%) ; POW > 18 *23908
#aa23 == BGN > -26 (95%) ; POW > 18 *25670
#aa24 == BGN > -24 (97%) ; POW > 18 *26903
#aa25 == BGN > -20 (99%) ; POW > 18 *28510

#aa26 == BGN > -29 (90%) ; POW > 20 *33135
#aa27 == BGN > -28 (92%) ; POW > 20 *34605
#aa28 == BGN > -26 (95%) ; POW > 20 *36602
#aa29 == BGN > -24 (97%) ; POW > 20 *37944
#aa30 == BGN > -20 (99%) ; POW > 20 *39613

#aa31 == BGN > -29 (90%) ; POW > 22 *43142
#aa32 == BGN > -28 (92%) ; POW > 22 *44706
#aa33 == BGN > -26 (95%) ; POW > 22 *46819
#aa34 == BGN > -24 (97%) ; POW > 22 *48210
#aa35 == BGN > -20 (99%) ; POW > 22 *49887

#aa36 == BGN > -29 (90%) ; POW > 24 *51275
#aa37 == BGN > -28 (92%) ; POW > 24 *52885
#aa38 == BGN > -26 (95%) ; POW > 24 *55033
#aa39 == BGN > -24 (97%) ; POW > 24 *56434
#aa40 == BGN > -20 (99%) ; POW > 24 *58111

#aa41 == BGN > -29 (90%) ; POW > 26 *57276
#aa42 == BGN > -28 (92%) ; POW > 26 *58899
#aa43 == BGN > -26 (95%) ; POW > 26 *61055
#aa44 == BGN > -24 (97%) ; POW > 26 *62456
#aa45 == BGN > -20 (99%) ; POW > 26 *64133

#aa46 == BGN > -29 (90%) ; POW > 28 *61043
#aa47 == BGN > -28 (92%) ; POW > 28 *62672
#aa48 == BGN > -26 (95%) ; POW > 28 *64828
#aa49 == BGN > -24 (97%) ; POW > 28 *66229
#aa50 == BGN > -20 (99%) ; POW > 28 *67906

#aa51 == BGN > -29 (90%) ; POW > 30 *63313
#aa52 == BGN > -28 (92%) ; POW > 30 *64942
#aa53 == BGN > -26 (95%) ; POW > 30 *67098
#aa54 == BGN > -24 (97%) ; POW > 30 *68499
#aa55 == BGN > -20 (99%) ; POW > 30 *70176


rm(list= ls()[(ls() %in% c("y", "hb", "hp"))])  

Sm <- data.frame(File=rep(c("SDE1","SDE2"), each=144), Min=rep(seq(1:144),2))
Sm.cols.text <- "Sm"; Sm.cols.number <- seq(1:55); Sm.cols <- paste(Sm.cols.text, Sm.cols.number, sep = "")
Sm[Sm.cols]<-NA


for (a in 3:57) {
  for (m in 1:144) {
    Sm[which(Sm$File==grep("SDE1", Sm$File, value = T) & Sm$Min==m),a] <- round((sum(index.f[which(index.f$File==grep("SED1", index.f$File, value = T) & index.f$Min==m), a+3])/256), 2)
    Sm[which(Sm$File==grep("SDE2", Sm$File, value = T) & Sm$Min==m),a] <- round((sum(index.f[which(index.f$File==grep("SED2", index.f$File, value = T) & index.f$Min==m), a+3])/256), 2)
  }  
}

write.csv(index.f, "piloto_201908/soundscape_saturation_unchose.csv", row.names = F)

rm(list= ls()[(ls() %in% c("a", "m", "Sm.cols.text", "Sm.cols.number", "Sm.cols"))])   


