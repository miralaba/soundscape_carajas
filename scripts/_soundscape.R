
#### setting working directory ####
setwd("~/capital_natural/soundscape")
memory.limit(1000000)

##### loading required packages ####
library(tuneR)
library(seewave)
library(soundecology)
library(parallel)
library(reshape2)
library(tidyverse)
library(dplyr)
library(lubridate)
library(scales)
library(multcomp)
library(multcompView)
library(emmeans)
library(MASS)
library(fitdistrplus)
library(lme4)
library(glmmTMB)
library(car)
library(MuMIn)
library(ggplot2)
library(cowplot)
library(TTR)


##### import data ####
list.sound.files <- list.files(path = getwd(), pattern = ".WAV$", full.names = T, recursive = T)
# checking
head(list.sound.files); tail(list.sound.files)

##### editing data ####
# converting list to dataframe
list.sound.files.df <- data.frame(list.sound.files)

list.sound.files.df <- separate(list.sound.files.df, 
                                list.sound.files, 
                                as.character(1:9), 
                                "/|[.]|-", 
                                extra = "merge")

# some editions in name of files -- column 8
list.sound.files.df$`10` <- ifelse(list.sound.files.df$`7`==2022, gsub("_", "", list.sound.files.df$`8`), list.sound.files.df$`8`)
list.sound.files.df$`10` <- as.numeric(substr(list.sound.files.df$`10`, 1, 12))
list.sound.files.df$`10` <- ifelse(list.sound.files.df$`7`==2019, list.sound.files.df$`10`-301, list.sound.files.df$`10`)
list.sound.files.df$`10` <- str_replace_all(list.sound.files.df$`10`, 
                                            c("9700" = "2100", "9710" = "2110", "9720" = "2120", "9730" = "2130", "9740" = "2140", "9750" = "2150",
                                              "9800" = "2200", "9810" = "2210", "9820" = "2220", "9830" = "2230", "9840" = "2240", "9850" = "2250",
                                              "9900" = "2300", "9910" = "2310", "9920" = "2320", "9930" = "2330", "9940" = "2340", "9950" = "2350"))
list.sound.files.df$`10` <- paste0(list.sound.files.df$`5`, "_", list.sound.files.df$`10`)
list.sound.files.df$ID <- 1:nrow(list.sound.files.df)

# keeping variables -- sample point and date
CN.data.raw <- list.sound.files.df[,c(11,10,5,7)]
names(CN.data.raw) <- c("ID", "ldt.index", "Local", "Year")

# adding variables -- month
CN.data.raw$Month <- as.numeric(substr(sub(".*_", "", CN.data.raw$ldt.index),5,6))

# adding variables -- day
CN.data.raw$Day <- as.numeric(substr(sub(".*_", "", CN.data.raw$ldt.index),7,8))

# adding variables -- hours
CN.data.raw$Hour <- as.numeric(substr(sub(".*_", "", CN.data.raw$ldt.index),9,10))

# adding variables -- minutes
minute.index <- c("*0000$", "*0010$", "*0020$", "*0030$", "*0040$", "*0050$",
                  "*0100$", "*0110$", "*0120$", "*0130$", "*0140$", "*0150$",
                  "*0200$", "*0210$", "*0220$", "*0230$", "*0240$", "*0250$",
                  "*0300$", "*0310$", "*0320$", "*0330$", "*0340$", "*0350$",
                  "*0400$", "*0410$", "*0420$", "*0430$", "*0440$", "*0450$",
                  "*0500$", "*0510$", "*0520$", "*0530$", "*0540$", "*0550$",
                  "*0600$", "*0610$", "*0620$", "*0630$", "*0640$", "*0650$",
                  "*0700$", "*0710$", "*0720$", "*0730$", "*0740$", "*0750$",
                  "*0800$", "*0810$", "*0820$", "*0830$", "*0840$", "*0850$",
                  "*0900$", "*0910$", "*0920$", "*0930$", "*0940$", "*0950$",
                  "*1000$", "*1010$", "*1020$", "*1030$", "*1040$", "*1050$",
                  "*1100$", "*1110$", "*1120$", "*1130$", "*1140$", "*1150$",
                  "*1200$", "*1210$", "*1220$", "*1230$", "*1240$", "*1250$",
                  "*1300$", "*1310$", "*1320$", "*1330$", "*1340$", "*1350$",
                  "*1400$", "*1410$", "*1420$", "*1430$", "*1440$", "*1450$",
                  "*1500$", "*1510$", "*1520$", "*1530$", "*1540$", "*1550$",
                  "*1600$", "*1610$", "*1620$", "*1630$", "*1640$", "*1650$",
                  "*1700$", "*1710$", "*1720$", "*1730$", "*1740$", "*1750$",
                  "*1800$", "*1810$", "*1820$", "*1830$", "*1840$", "*1850$",
                  "*1900$", "*1910$", "*1920$", "*1930$", "*1940$", "*1950$",
                  "*2000$", "*2010$", "*2020$", "*2030$", "*2040$", "*2050$",
                  "*2100$", "*2110$", "*2120$", "*2130$", "*2140$", "*2150$",
                  "*2200$", "*2210$", "*2220$", "*2230$", "*2240$", "*2250$",
                  "*2300$", "*2310$", "*2320$", "*2330$", "*2340$", "*2350$")

for (t in 1:length(minute.index)) {
  min_da_vez<-grep(minute.index[t], CN.data.raw$ldt.index, value = T)
  CN.data.raw$Min[CN.data.raw$ldt.index %in% min_da_vez]<-t-1
}

# adding variables -- path to sound files
CN.data.raw$path2file <- list.sound.files

# adding variables -- soundscape index total frequency: 0.3 - 22.05kHz
CN.data.raw$HTotal=NA; CN.data.raw$BGNTotal=NA; CN.data.raw$S2NTotal=NA; CN.data.raw$NPeakTotal=NA; CN.data.raw$AATotal=NA; CN.data.raw$ADITotal=NA; CN.data.raw$AEITotal=NA

# adding variables -- soundscape index frequency bin 1: 0.3 - 4 
CN.data.raw$Hfbin1=NA; CN.data.raw$BGNfbin1=NA; CN.data.raw$S2Nfbin1=NA; CN.data.raw$NPeakfbin1=NA; CN.data.raw$AAfbin1=NA; CN.data.raw$ADIfbin1=NA; CN.data.raw$AEIfbin1=NA

# adding variables -- soundscape index frequency bin 2: 4 - 12 
CN.data.raw$Hfbin2=NA; CN.data.raw$BGNfbin2=NA; CN.data.raw$S2Nfbin2=NA; CN.data.raw$NPeakfbin2=NA; CN.data.raw$AAfbin2=NA; CN.data.raw$ADIfbin2=NA; CN.data.raw$AEIfbin2=NA

# adding variables -- soundscape index frequency bin 3: 0.3 - 12 
CN.data.raw$Hfbin3=NA; CN.data.raw$BGNfbin3=NA; CN.data.raw$S2Nfbin3=NA; CN.data.raw$NPeakfbin3=NA; CN.data.raw$AAfbin3=NA; CN.data.raw$ADIfbin3=NA; CN.data.raw$AEIfbin3=NA

# adding variables -- soundscape index frequency bin 4: 12 - 22.05 
CN.data.raw$Hfbin4=NA; CN.data.raw$BGNfbin4=NA; CN.data.raw$S2Nfbin4=NA; CN.data.raw$NPeakfbin4=NA; CN.data.raw$AAfbin4=NA; CN.data.raw$ADIfbin4=NA; CN.data.raw$AEIfbin4=NA

# reordering columns
CN.data.raw <- CN.data.raw[,c(1,3:8,10:44,2,9)]

# saving
dir.create("data")
write.csv(CN.data.raw, "data/AcousticIndex_112022_v1.csv", na = "NA", row.names = F)

rm(list= ls()[!(ls() %in% c("CN.data.raw"))])
gc()
## checking
#ex1 <- readWave(l[1])
#spectro(ex1)
#listen(ex1)


#function to get modal values
#getmode <- function(v) {
#  uniqv <- unique(v)
#  uniqv[which.max(tabulate(match(v, uniqv)))]
#}

##### getting soundscape indexes for total frequency ####

for (i in 1:nrow(CN.data.raw)) {
  
  if(file.info(CN.data.raw$path2file[i])$size<5242880) next
  
  S <- normalize(downsample(readWave(CN.data.raw$path2file[i]), 44100), unit = "1")
  Stotal <- ffilter(S, from = 300, to = 22050, output="Wave")
  Sbin1 <- ffilter(S, from = 300, to = 4000, output="Wave")
  Sbin2 <- ffilter(S, from = 4000, to = 12000, output="Wave")
  Sbin3 <- ffilter(S, from = 300, to = 12000, output="Wave")
  Sbin4 <- ffilter(S, from = 12000, to = 22050, output="Wave")
  
  # total entropy of a time wave
  CN.data.raw$HTotal[i] <- H(Stotal, wl = 512)
  CN.data.raw$Hfbin1[i] <- H(Sbin1, wl = 512)
  CN.data.raw$Hfbin2[i] <- H(Sbin2, wl = 512)
  CN.data.raw$Hfbin3[i] <- H(Sbin3, wl = 512)
  CN.data.raw$Hfbin4[i] <- H(Sbin4, wl = 512)
  
  # acoustic diversity index
  CN.data.raw$ADITotal[i] <- acoustic_diversity(Stotal)$adi_left
  CN.data.raw$ADIfbin1[i] <- acoustic_diversity(Sbin1)$adi_left
  CN.data.raw$ADIfbin2[i] <- acoustic_diversity(Sbin2)$adi_left
  CN.data.raw$ADIfbin3[i] <- acoustic_diversity(Sbin3)$adi_left
  CN.data.raw$ADIfbin4[i] <- acoustic_diversity(Sbin4)$adi_left
  
  
  # acoustic evenness index
  CN.data.raw$AEITotal[i] <- acoustic_evenness(Stotal)$aei_left
  CN.data.raw$AEIfbin1[i] <- acoustic_evenness(Sbin1)$aei_left
  CN.data.raw$AEIfbin2[i] <- acoustic_evenness(Sbin2)$aei_left
  CN.data.raw$AEIfbin3[i] <- acoustic_evenness(Sbin3)$aei_left
  CN.data.raw$AEIfbin4[i] <- acoustic_evenness(Sbin4)$aei_left
  
  # spectrogram (freqency X time X amplitude) -- Total
  spectrogram.data.total <- spectro(Stotal, wl = 512, plot = F)
  
  # isolating components: amplitude
  amp.total <- melt(spectrogram.data.total$amp, value.name = "Amplitude") %>%
                dplyr::select(FrequencyIndex = Var1, TimeIndex = Var2, Amplitude)
  # isolating components: frequency
  freq.total <- melt(spectrogram.data.total$freq, value.name = "Frequency") %>%
                mutate(FrequencyIndex = row_number(), Frequency = Frequency)
  # isolating components: time
  tm.total <- melt(spectrogram.data.total$time, value.name = "Time") %>%
                mutate(TimeIndex = row_number())
  
  # consolidating spectrogram data into dataframe
  spectrogram.df.total <- amp.total %>%
                            left_join(freq.total, by = "FrequencyIndex") %>%
                            left_join(tm.total, by = "TimeIndex") %>%
                            dplyr::select(Time, Frequency, Amplitude)
  
  # frequeny precision: 3 decimals
  spectrogram.df.total$Frequency <- round(spectrogram.df.total$Frequency, 3)
  
  # amplitude rescale to dB (-50, -3) and precision to no decimals
  spectrogram.df.total$Amplitude.rescaled <- rescale(spectrogram.df.total$Amplitude, to = c(-50,-3))
  spectrogram.df.total$Amplitude.rescaled <- round(spectrogram.df.total$Amplitude.rescaled, 0)
  
  # distribution of amplitude values
  Amp.intesities.total <- hist(spectrogram.df.total$Amplitude.rescaled, breaks = 999, plot=F)
  # maximum simple moving average value -- background noise
  CN.data.raw$BGNTotal[i] <- round(Amp.intesities.total$breaks[which.max(SMA(Amp.intesities.total$counts, n=5))], 0)
  
  # difference between maximum amplitude and background noise -- sound to noise ratio
  CN.data.raw$S2NTotal[i] <- max(spectrogram.df.total$Amplitude.rescaled) - CN.data.raw$BGNTotal[i]
  
  # N times amplitude is 3dB greater than background noise
  spectrogram.df.total$NPeak <- NA
  spectrogram.df.total$NPeak <- ifelse(spectrogram.df.total$Amplitude.rescaled>CN.data.raw$BGNTotal[i]+3, 1,0)
  CN.data.raw$NPeakTotal[i] <- sum(spectrogram.df.total$NPeak, na.rm = T)
  
  # proportion of cells in spectrogram with amplitude 3dB greater than background noise -- acoustic activity
  CN.data.raw$AATotal[i] <- round(CN.data.raw$NPeakTotal[i]/nrow(spectrogram.df.total),3)
  
  # spectrogram (freqency X time X amplitude) -- fbin1
  spectrogram.data.fbin1 <- spectro(Sbin1, wl = 512, plot = F)
  
  # isolating components: amplitude
  amp.fbin1 <- melt(spectrogram.data.fbin1$amp, value.name = "Amplitude") %>%
    dplyr::select(FrequencyIndex = Var1, TimeIndex = Var2, Amplitude)
  # isolating components: frequency
  freq.fbin1 <- melt(spectrogram.data.fbin1$freq, value.name = "Frequency") %>%
    mutate(FrequencyIndex = row_number(), Frequency = Frequency)
  # isolating components: time
  tm.fbin1 <- melt(spectrogram.data.fbin1$time, value.name = "Time") %>%
    mutate(TimeIndex = row_number())
  
  # consolidating spectrogram data into dataframe
  spectrogram.df.fbin1 <- amp.fbin1 %>%
    left_join(freq.fbin1, by = "FrequencyIndex") %>%
    left_join(tm.fbin1, by = "TimeIndex") %>%
    dplyr::select(Time, Frequency, Amplitude)
  
  # frequeny precision: 3 decimals
  spectrogram.df.fbin1$Frequency <- round(spectrogram.df.fbin1$Frequency, 3)
  
  # amplitude rescale to dB (-50, -3) and precision to no decimals
  spectrogram.df.fbin1$Amplitude.rescaled <- rescale(spectrogram.df.fbin1$Amplitude, to = c(-50,-3))
  spectrogram.df.fbin1$Amplitude.rescaled <- round(spectrogram.df.fbin1$Amplitude.rescaled, 0)
  
  # distribution of amplitude values
  Amp.intesities.fbin1 <- hist(spectrogram.df.fbin1$Amplitude.rescaled, breaks = 999, plot=F)
  # maximum simple moving average value -- background noise
  CN.data.raw$BGNfbin1[i] <- round(Amp.intesities.fbin1$breaks[which.max(SMA(Amp.intesities.fbin1$counts, n=5))], 0)
  
  # difference between maximum amplitude and background noise -- sound to noise ratio
  CN.data.raw$S2Nfbin1[i] <- max(spectrogram.df.fbin1$Amplitude.rescaled) - CN.data.raw$BGNfbin1[i]
  
  # N times amplitude is 3dB greater than background noise
  spectrogram.df.fbin1$NPeak <- NA
  spectrogram.df.fbin1$NPeak <- ifelse(spectrogram.df.fbin1$Amplitude.rescaled>CN.data.raw$BGNfbin1[i]+3, 1,0)
  CN.data.raw$NPeakfbin1[i] <- sum(spectrogram.df.fbin1$NPeak, na.rm = T)
  
  # proportion of cells in spectrogram with amplitude 3dB greater than background noise -- acoustic activity
  CN.data.raw$AAfbin1[i] <- round(CN.data.raw$NPeakfbin1[i]/nrow(spectrogram.df.fbin1),3)
  
  # spectrogram (freqency X time X amplitude) -- fbin2
  spectrogram.data.fbin2 <- spectro(Sbin2, wl = 512, plot = F)
  
  # isolating components: amplitude
  amp.fbin2 <- melt(spectrogram.data.fbin2$amp, value.name = "Amplitude") %>%
    dplyr::select(FrequencyIndex = Var1, TimeIndex = Var2, Amplitude)
  # isolating components: frequency
  freq.fbin2 <- melt(spectrogram.data.fbin2$freq, value.name = "Frequency") %>%
    mutate(FrequencyIndex = row_number(), Frequency = Frequency)
  # isolating components: time
  tm.fbin2 <- melt(spectrogram.data.fbin2$time, value.name = "Time") %>%
    mutate(TimeIndex = row_number())
  
  # consolidating spectrogram data into dataframe
  spectrogram.df.fbin2 <- amp.fbin2 %>%
    left_join(freq.fbin2, by = "FrequencyIndex") %>%
    left_join(tm.fbin2, by = "TimeIndex") %>%
    dplyr::select(Time, Frequency, Amplitude)
  
  # frequeny precision: 3 decimals
  spectrogram.df.fbin2$Frequency <- round(spectrogram.df.fbin2$Frequency, 3)
  
  # amplitude rescale to dB (-50, -3) and precision to no decimals
  spectrogram.df.fbin2$Amplitude.rescaled <- rescale(spectrogram.df.fbin2$Amplitude, to = c(-50,-3))
  spectrogram.df.fbin2$Amplitude.rescaled <- round(spectrogram.df.fbin2$Amplitude.rescaled, 0)
  
  # distribution of amplitude values
  Amp.intesities.fbin2 <- hist(spectrogram.df.fbin2$Amplitude.rescaled, breaks = 999, plot=F)
  # maximum simple moving average value -- background noise
  CN.data.raw$BGNfbin2[i] <- round(Amp.intesities.fbin2$breaks[which.max(SMA(Amp.intesities.fbin2$counts, n=5))], 0)
  
  # difference between maximum amplitude and background noise -- sound to noise ratio
  CN.data.raw$S2Nfbin2[i] <- max(spectrogram.df.fbin2$Amplitude.rescaled) - CN.data.raw$BGNfbin2[i]
  
  # N times amplitude is 3dB greater than background noise
  spectrogram.df.fbin2$NPeak <- NA
  spectrogram.df.fbin2$NPeak <- ifelse(spectrogram.df.fbin2$Amplitude.rescaled>CN.data.raw$BGNfbin2[i]+3, 1,0)
  CN.data.raw$NPeakfbin2[i] <- sum(spectrogram.df.fbin2$NPeak, na.rm = T)
  
  # proportion of cells in spectrogram with amplitude 3dB greater than background noise -- acoustic activity
  CN.data.raw$AAfbin2[i] <- round(CN.data.raw$NPeakfbin2[i]/nrow(spectrogram.df.fbin2),3)
  
  # spectrogram (freqency X time X amplitude) -- fbin3
  spectrogram.data.fbin3 <- spectro(Sbin3, wl = 512, plot = F)
  
  # isolating components: amplitude
  amp.fbin3 <- melt(spectrogram.data.fbin3$amp, value.name = "Amplitude") %>%
    dplyr::select(FrequencyIndex = Var1, TimeIndex = Var2, Amplitude)
  # isolating components: frequency
  freq.fbin3 <- melt(spectrogram.data.fbin3$freq, value.name = "Frequency") %>%
    mutate(FrequencyIndex = row_number(), Frequency = Frequency)
  # isolating components: time
  tm.fbin3 <- melt(spectrogram.data.fbin3$time, value.name = "Time") %>%
    mutate(TimeIndex = row_number())
  
  # consolidating spectrogram data into dataframe
  spectrogram.df.fbin3 <- amp.fbin3 %>%
    left_join(freq.fbin3, by = "FrequencyIndex") %>%
    left_join(tm.fbin3, by = "TimeIndex") %>%
    dplyr::select(Time, Frequency, Amplitude)
  
  # frequeny precision: 3 decimals
  spectrogram.df.fbin3$Frequency <- round(spectrogram.df.fbin3$Frequency, 3)
  
  # amplitude rescale to dB (-50, -3) and precision to no decimals
  spectrogram.df.fbin3$Amplitude.rescaled <- rescale(spectrogram.df.fbin3$Amplitude, to = c(-50,-3))
  spectrogram.df.fbin3$Amplitude.rescaled <- round(spectrogram.df.fbin3$Amplitude.rescaled, 0)
  
  # distribution of amplitude values
  Amp.intesities.fbin3 <- hist(spectrogram.df.fbin3$Amplitude.rescaled, breaks = 999, plot=F)
  # maximum simple moving average value -- background noise
  CN.data.raw$BGNfbin3[i] <- round(Amp.intesities.fbin3$breaks[which.max(SMA(Amp.intesities.fbin3$counts, n=5))], 0)
  
  # difference between maximum amplitude and background noise -- sound to noise ratio
  CN.data.raw$S2Nfbin3[i] <- max(spectrogram.df.fbin3$Amplitude.rescaled) - CN.data.raw$BGNfbin3[i]
  
  # N times amplitude is 3dB greater than background noise
  spectrogram.df.fbin3$NPeak <- NA
  spectrogram.df.fbin3$NPeak <- ifelse(spectrogram.df.fbin3$Amplitude.rescaled>CN.data.raw$BGNfbin3[i]+3, 1,0)
  CN.data.raw$NPeakfbin3[i] <- sum(spectrogram.df.fbin3$NPeak, na.rm = T)
  
  # proportion of cells in spectrogram with amplitude 3dB greater than background noise -- acoustic activity
  CN.data.raw$AAfbin3[i] <- round(CN.data.raw$NPeakfbin3[i]/nrow(spectrogram.df.fbin3),3)
  
  # spectrogram (freqency X time X amplitude) -- fbin4
  spectrogram.data.fbin4 <- spectro(Sbin4, wl = 512, plot = F)
  
  # isolating components: amplitude
  amp.fbin4 <- melt(spectrogram.data.fbin4$amp, value.name = "Amplitude") %>%
    dplyr::select(FrequencyIndex = Var1, TimeIndex = Var2, Amplitude)
  # isolating components: frequency
  freq.fbin4 <- melt(spectrogram.data.fbin4$freq, value.name = "Frequency") %>%
    mutate(FrequencyIndex = row_number(), Frequency = Frequency)
  # isolating components: time
  tm.fbin4 <- melt(spectrogram.data.fbin4$time, value.name = "Time") %>%
    mutate(TimeIndex = row_number())
  
  # consolidating spectrogram data into dataframe
  spectrogram.df.fbin4 <- amp.fbin4 %>%
    left_join(freq.fbin4, by = "FrequencyIndex") %>%
    left_join(tm.fbin4, by = "TimeIndex") %>%
    dplyr::select(Time, Frequency, Amplitude)
  
  # frequeny precision: 3 decimals
  spectrogram.df.fbin4$Frequency <- round(spectrogram.df.fbin4$Frequency, 3)
  
  # amplitude rescale to dB (-50, -3) and precision to no decimals
  spectrogram.df.fbin4$Amplitude.rescaled <- rescale(spectrogram.df.fbin4$Amplitude, to = c(-50,-3))
  spectrogram.df.fbin4$Amplitude.rescaled <- round(spectrogram.df.fbin4$Amplitude.rescaled, 0)
  
  # distribution of amplitude values
  Amp.intesities.fbin4 <- hist(spectrogram.df.fbin4$Amplitude.rescaled, breaks = 999, plot=F)
  # maximum simple moving average value -- background noise
  CN.data.raw$BGNfbin4[i] <- round(Amp.intesities.fbin4$breaks[which.max(SMA(Amp.intesities.fbin4$counts, n=5))], 0)
  
  # difference between maximum amplitude and background noise -- sound to noise ratio
  CN.data.raw$S2Nfbin4[i] <- max(spectrogram.df.fbin4$Amplitude.rescaled) - CN.data.raw$BGNfbin4[i]
  
  # N times amplitude is 3dB greater than background noise
  spectrogram.df.fbin4$NPeak <- NA
  spectrogram.df.fbin4$NPeak <- ifelse(spectrogram.df.fbin4$Amplitude.rescaled>CN.data.raw$BGNfbin4[i]+3, 1,0)
  CN.data.raw$NPeakfbin4[i] <- sum(spectrogram.df.fbin4$NPeak, na.rm = T)
  
  # proportion of cells in spectrogram with amplitude 3dB greater than background noise -- acoustic activity
  CN.data.raw$AAfbin4[i] <- round(CN.data.raw$NPeakfbin4[i]/nrow(spectrogram.df.fbin4),3)
  
  write.csv(CN.data.raw, "data/AcousticIndex_102022_v1.csv", na = "NA", row.names = F)
  cat('\n>finalizado', CN.data.raw$ldt.index[i], 'faltam', nrow(CN.data.raw)-i, '<\n')
  
  rm(list= ls()[!(ls() %in% c("CN.data.raw", "i"))])
  gc()
  
}

#


#### Cleaning based on data ####

CN.data.raw <- read.csv("data/AcousticIndex_102022_v1.csv", header = T)
## checking
#head(CN.data.raw);tail(CN.data.raw)
#str(CN.data.raw)
#summary(CN.data.raw)

# excluding data based on size file [<5MB -- see line 130]
CN.data.raw.ed <- CN.data.raw[!is.na(CN.data.raw$AEITotal),]

# saving
write.csv(CN.data.raw.ed, "data/AcousticIndex_102022_v2.csv", na = "NA", row.names = F)
#

# excluding audios based on background noise [<= -40 or >= -10 dB]
CN.data.raw.ed2 <- CN.data.raw.ed

hist(CN.data.raw.ed2$BGNTotal)
range(CN.data.raw.ed2$BGNTotal)
nrow(CN.data.raw.ed2[CN.data.raw.ed2$BGNTotal <= -40 | CN.data.raw.ed2$BGNTotal >= -10, ])
CN.data.raw.ed2[CN.data.raw.ed2$BGNTotal <= -40  | CN.data.raw.ed2$BGNTotal >= -10, "ldt.index"]

exclude <- CN.data.raw.ed2[CN.data.raw.ed2$BGNTotal <= -40 | CN.data.raw.ed2$BGNTotal >= -10, "path2file"]
## checking
#(s <- exclude[1])
#play(readWave(s))
#
CN.data.raw.ed2 <- CN.data.raw.ed2[!CN.data.raw.ed2$path2file %in% exclude,]

# excluding audios based on acoustic activity [<= 0.001]
hist(CN.data.raw.ed2$AATotal)
range(CN.data.raw.ed2$AATotal)
nrow(CN.data.raw.ed2[CN.data.raw.ed2$AATotal <= .001, ])
CN.data.raw.ed2[CN.data.raw.ed2$AATotal <= .001, "ldt.index"]

exclude <- CN.data.raw.ed2[CN.data.raw.ed2$AATotal <= .001, "path2file"]
## checking
#(s <- exclude[1])
#play(readWave(s))
#
CN.data.raw.ed2 <- CN.data.raw.ed2[!CN.data.raw.ed2$path2file %in% exclude,]

write.csv(CN.data.raw.ed2, "data/AcousticIndex_102022_v3.csv", na = "NA", row.names = F)

#rm(list= ls()[!(ls() %in% c("CN.data.raw.ed2"))])
#gc()
#


