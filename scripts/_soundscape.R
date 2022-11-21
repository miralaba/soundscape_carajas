
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
  
  write.csv(CN.data.raw, "AcousticIndex_102022_v1.csv", na = "NA", row.names = F)
  cat('\n>finalizado', CN.data.raw$ldt.index[i], 'faltam', nrow(CN.data.raw)-i, '<\n')
  
  rm(list= ls()[!(ls() %in% c("CN.data.raw", "i"))])
  gc()
  
}

#


#### Cleaning based on data ####

CN.data.raw <- read.csv("AcousticIndex_102022_v1.csv", header = T)
## checking
#head(CN.data.raw);tail(CN.data.raw)
#str(CN.data.raw)
#summary(CN.data.raw)

# excluding data based on size file [<5MB -- see line 130]
CN.data <- CN.data.raw[!is.na(CN.data.raw$AEI),]

# saving
write.csv(CN.data, "AcousticIndex_102022_v2.csv", na = "NA", row.names = F)
#

# excluding audios based on background noise [<= -40 or >= -10 dB]
CN.data.raw <- CN.data

hist(CN.data.raw$BGN)
range(CN.data.raw$BGN)
nrow(CN.data.raw[CN.data.raw$BGN <= -40 | CN.data.raw$BGN >= -10, ])
CN.data.raw[CN.data.raw$BGN <= -40  | CN.data.raw$BGN >= -10, "ldt.index"]

exclude <- CN.data.raw[CN.data.raw$BGN <= -40 | CN.data.raw$BGN >= -10, "path2file"]
## checking
#(s <- exclude[1])
#play(readWave(s))
#
CN.data.raw <- CN.data.raw[!CN.data.raw$path2file %in% exclude,]

# excluding audios based on acoustic activity [<= 0.001]
hist(CN.data.raw$AA)
range(CN.data.raw$AA)
nrow(CN.data.raw[CN.data.raw$AA <= .001, ])
CN.data.raw[CN.data.raw$AA <= .001, "ldt.index"]

exclude <- CN.data.raw[CN.data.raw$AA <= .001, "path2file"]
## checking
#(s <- exclude[1])
#play(readWave(s))
#
CN.data.raw <- CN.data.raw[!CN.data.raw$path2file %in% exclude,]

write.csv(CN.data.raw, "AcousticIndex_102022_v3.csv", na = "NA", row.names = F)

rm(list= ls()[!(ls() %in% c("CN.data.raw"))])
gc()
#


# reordering columns
#CN.data.raw <- read.csv("AcousticIndex_102022_v3.csv", header = T)
CN.data.raw <- CN.data.raw[,-c(13:14,20:23)]
CN.data.raw <- CN.data.raw[,c(1:7,15,8:12,16:17,13:14)]
names(CN.data.raw)[14:15] <- c("ADI", "AEI")

#### exploratory ####

CN.data <- CN.data.raw

CN.data <- CN.data %>% 
  mutate_at(vars(ID:Min), factor) %>%
  mutate(Date = parse_date_time(sub(".*_", "", CN.data$ldt.index), orders = "%Y%m%d%H%M")) %>% 
  mutate(Local=fct_relevel(Local,paste0("CN", seq(1:14))))

str(CN.data)

#### between minutes by points ####
dir.create("exploratory")

#Visualize your data
ggplot(CN.data, aes(Min, H))+
  geom_boxplot()+
  scale_x_discrete(breaks=c("0","36","72", "108", "143"), labels=c("00:00", "06:00", "12:00", "18:00", ""))+
  facet_wrap(~Local, ncol = 2)


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


