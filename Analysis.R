#### **TERRESTRIAL** ----

### Acute Fear Summer 2019 Proj
### N.R. Sommer & R.W. Buchkowski

require(tidyverse)
require(lubridate)

#### Read-in Data ####
IDtable = read_csv("Data/IDtable_13Sep2019.csv") %>%
  filter(!(Individual =="K" & Species =="TRRA")) %>% # remove lost isopod
  filter(!(Individual =="B_remove" & Species =="Lynx")) %>% # remove lost lynx
  rename(Start_0 = Start,End_0 = End, WWt = `Wet Weight (g)`) %>%
  gather(-Species:-Treatment, -WWt:-Year, key=Type, value=time) %>%
  separate(Type, into=c("Type", "Rep"), sep="_") %>%
  spread(key=Type, value=time)

respdata = read_csv("Data/19July2019_crickets.csv") %>% 
  select(Time:Flow) %>%
  mutate(Date = 19200) %>% # Date = yday
  bind_rows(
    read_csv("Data/20July2019_crickets.csv") %>% 
      select(Time:Flow) %>%
      mutate(Date = 19201)
  ) %>%
  bind_rows(
    read_csv("Data/21July2019_crickets.csv") %>% 
      select(Time:Flow) %>%
      mutate(Date = 19202)
  ) %>%
  bind_rows(
    read_csv("Data/22July2019_isopod_fear.csv") %>% 
      select(Time:Flow) %>%
      mutate(Date = 19203)
  ) %>%
  bind_rows(
    read_csv("Data/isopod_fear_24July2019.csv") %>% 
      select(Time:Flow) %>%
      mutate(Date = 19205)
  ) %>%
  bind_rows(
    read_csv("Data/26July2019_spiders.csv") %>% 
      select(Time:Flow) %>%
      mutate(Date = 19207)
  ) %>%
  bind_rows(
    read_csv("Data/29July2019.csv") %>% 
      select(Time:Flow) %>%
      rename(`Probe Temp` = Temperature) %>%
      mutate(Date = 19210)
  ) %>%
  bind_rows(
    read_csv("Data/30July2019.csv") %>% 
      select(Time:Flow) %>%
      mutate(Date = 19211)
  ) %>%
  bind_rows(
    read_csv("Data/31July2019.csv") %>% 
      select(Time:Flow) %>%
      mutate(Date = 19212)
  ) %>%
  bind_rows(
    read_csv("Data/1Aug2019.csv") %>% 
      select(Time:Flow) %>%
      mutate(Date = 19213)
  ) %>%
  bind_rows(
    read_csv("Data/20Aug2019.csv") %>% 
      select(Time:Flow) %>%
      mutate(Date = 19232)
  ) %>%
  bind_rows(
    read_csv("Data/21Aug2019.csv") %>% 
      select(Time:Flow) %>%
      mutate(Date = 19233)
  ) %>%
  bind_rows(
    read_csv("Data/23Aug2019.csv") %>% 
      select(Time:Flow) %>%
      mutate(Date = 19235)
  ) %>%
  bind_rows(
    read_csv("Data/24Aug2019.csv") %>% 
      select(Time:Flow) %>%
      mutate(Date = 19236)
  ) %>%
  bind_rows(
    read_csv("Data/28Aug2019.csv") %>% 
      select(Time:Flow) %>%
      mutate(Date = 19240)
  ) %>%
  bind_rows(
    read_csv("Data/29Aug2019.csv") %>% 
      select(Time:Flow) %>%
      mutate(Date = 19241)
  ) %>%
  bind_rows(
    read_csv("Data/3Sep2019.csv") %>% 
      select(Time:Flow) %>%
      mutate(Date = 19246)
  ) %>%
  bind_rows(
    read_csv("Data/4Sep2019.csv") %>% 
      select(Time:Flow) %>%
      mutate(Date = 19247)
  ) %>%
  bind_rows(
    read_csv("Data/5Sep2019.csv") %>% 
      select(Time:Flow) %>%
      mutate(Date = 19248)
  ) %>%
  bind_rows(
    read_csv("Data/7Sep2019.csv") %>% 
      select(Time:Flow) %>%
      mutate(Date = 19250)
  ) %>%
  bind_rows(
    read_csv("Data/8Sep2019.csv") %>% 
      select(Time:Flow) %>%
      mutate(Date = 19251)
  ) %>%
  bind_rows(
    read_csv("Data/9Sep2019.csv") %>% 
      select(Time:Flow) %>%
      mutate(Date = 19252)
  ) %>%
  bind_rows(
    read_csv("Data/10Sep2019.csv") %>% 
      select(Time:Flow) %>%
      mutate(Date = 19253)
  ) %>%
  bind_rows(
    read_csv("Data/11Sep2019.csv") %>% 
      select(Time:Flow) %>%
      mutate(Date = 19254)
  ) %>%
  bind_rows(
    read_csv("Data/12Sep2019.csv") %>% 
      select(Time:Flow) %>%
      mutate(Date = 19255))

colnames(respdata) = c("Time", "QO2", "FlowV",
                       "CO2", "Temp", "Flow",
                       "Julian", "Flow2")

respdata = respdata %>% select(-QO2, -FlowV, -Flow2)

NN = dim(IDtable)[1]

#### Calculate respiration rates ####

runeach = function(XX = 1, IDTABLE = IDtable, RESPDATA = respdata){ 
  select = RESPDATA$Time > as.numeric(IDTABLE[XX,"Start"]) & 
    RESPDATA$Time < as.numeric(IDTABLE[XX,"End"]) &
    RESPDATA$Julian == as.numeric(IDTABLE[XX, "Julian"])
  
  selected = RESPDATA[select,]
  
  PRESSURE = IDTABLE[XX,"Pressure"]*25.4
  
  m1 = lm(CO2 ~ Time, data=selected)
  slope = coefficients(m1)[2]
  r2 = summary(m1)$adj.r.squared
  TEMP=mean(selected$Temp)
  Weight = as.numeric(IDTABLE[XX,"WWt"])
  WEIGHT = ifelse(is.na(Weight), 1, Weight)
  WEIGHT0 = ifelse(is.na(Weight), 0, Weight)
  
  if(IDTABLE$Julian[XX] < 19236) { 
    vol = 60.79535394
  } else {
    if(IDTABLE$Julian[XX] > 19236) {
      vol = 53.7489539
    } else {
      if(IDTABLE$Julian[XX] == 19236 && IDTABLE$Species[XX] == "TRRA") {
        vol = 60.79535394
      } else {
        if(IDTABLE$Julian[XX] == 19236 && IDTABLE$Species[XX] == "Lynx") {
          vol = 53.7489539 }
      }
    }
  }

  volume = ((vol-WEIGHT0)*((TEMP+273)/273)*760/PRESSURE)/1000
  resprate = as.numeric(slope*volume*60/(WEIGHT))
  c(r2 = r2, resprate = resprate)
}

output = cbind(IDtable, t(sapply(seq(1, NN, 1), FUN = runeach)))

# diagnose errors by line, if any
for(kk in 1:NN){ 
  aaaa = runeach(kk)
  print(kk)
}

# corrected CO2 rate by blanks
blank = output %>% filter(Individual=="Blank") %>%
  group_by(Julian) %>%
  summarize(blank = mean(resprate))

fdata = output %>% filter((Treatment %in% c("B", "C", "N"))) %>%
  left_join(blank) %>%
  mutate(resp_rate = resprate - blank)

# get mean of corrected CO2 (avg of 3 obs) 
sumfdata = fdata %>% group_by(Species, Individual, Treatment) %>%
  summarize(resp_rate = mean(resp_rate))

fdata %>% group_by(Species, Individual, Treatment) %>%
  summarize(resp_ratesd = sd(resp_rate)) %>%
  left_join(sumfdata) %>%
  mutate(cv = resp_ratesd/resp_rate)

### Cleaning ####

# Check blanks
blank %>% ggplot(aes(x=Julian, y = blank)) + geom_point() # good

# Check consistency of observation length 
fdata %>% mutate(End-Start == 200) %>% View() 
# Observation lenght incorrect only for ONAS T, treatment N. End is 41360, start is 41150. 

# Check R2 for resp rates
fdata %>% filter(Species == "Cricket") %>% 
  ggplot(aes(x=resp_rate, y = r2, color = Individual)) + geom_point() # good

fdata %>% filter(Species == "Lynx") %>% 
  ggplot(aes(x=resp_rate, y = r2, color = Individual)) + geom_point() # bad

fdata %>% filter(Species == "MEFE") %>% 
  ggplot(aes(x=resp_rate, y = r2, color = Individual)) + geom_point() # good

fdata %>% filter(Species == "ONAS") %>% 
  ggplot(aes(x=resp_rate, y = r2, color = Individual)) + geom_point() # outlier?

fdata %>% filter(Species == "Phiddipus") %>% 
  ggplot(aes(x=resp_rate, y = r2, color = Individual)) + geom_point() # outlier?

fdata %>% filter(Species == "TRRA") %>% 
  ggplot(aes(x=resp_rate, y = r2, color = Individual)) + geom_point() # good

### Mixed-effects models ####
# NRS to-do: use lme4 and bootstrap confidence values 
# bootMER(), perhaps predictInterval() 

require(lme4)
require(boot)

# Filter data by species for mixed effects models 
cricketdf = sumfdata %>% filter(Species=="Cricket")
MEFEdf = sumfdata %>% filter(Species== "MEFE")
PHIDdf = sumfdata %>% filter(Species=="Phiddipus")
lynxdf = sumfdata %>% filter(Species=="Lynx")
ONASdf = sumfdata %>% filter(Species=="ONAS")
TRRAdf = sumfdata %>% filter(Species=="TRRA")

cricketdf$Treatment <- factor(cricketdf$Treatment, levels = c("N", "C", "B"))
MEFEdf$Treatment <- factor(MEFEdf$Treatment, levels = c("N", "C", "B"))
PHIDdf$Treatment <- factor(PHIDdf$Treatment, levels = c("N","C","B"))
lynxdf$Treatment  <- factor(lynxdf$Treatment, levels = c("N", "C", "B"))
ONASdf$Treatment <- factor(ONASdf$Treatment, levels = c("N", "C", "B"))
TRRAdf$Treatment <- factor(TRRAdf$Treatment, levels = c("N", "C", "B"))

# LME by species

cricketm1 = lmer(resp_rate~Treatment + (1|Individual), data=cricketdf)
plot(cricketm1)
summary(cricketm1)
# messing around with bootMer
cricketm1_boot <- bootMer(x=cricketm1,FUN=fixef,nsim=200) 
boot.ci(cricketm1_boot,type="perc",index=1)

MEFEm1 = lmer(resp_rate~Treatment + (1|Individual), data = MEFEdf)
plot(MEFEm1)
summary(MEFEm1) 

PHIDm1 = lmer(resp_rate~Treatment + (1|Individual), data = PHIDdf)
plot(PHIDm1)
summary(PHIDm1)

lynxm1 = lmer(resp_rate~Treatment + (1|Individual), data = lynxdf)
plot(lynxm1)
summary(lynxm1)

ONASm1 = lmer(resp_rate~Treatment + (1|Individual), data = ONASdf)
plot(ONASm1)
summary(ONASm1)

TRRAm1 = lmer(resp_rate~Treatment + (1|Individual), data=TRRAdf)
plot(TRRAm1)
summary(TRRAm1)

# Does removing low r2 and negative respiration change results? First pass says no, results consistent

PHIDdf2 = fdata %>% filter(Species == "Phiddipus") %>% 
  filter(!resp_rate < -1 & !r2 < .25)
PHIDdf2$Treatment <- factor(PHIDdf$Treatment, levels = c("N","C","B"))
PHIDm2 = lmer(resp_rate~Treatment + (1|Individual), data = PHIDdf2)
plot(PHIDm2)
summary(PHIDm2)

ONASdf2 = fdata %>% filter(Species == "ONAS") %>%
  filter(!resp_rate < -1 & !r2 <.25)
ONASdf2$Treatment <- factor(ONASdf$Treatment, levels = c("N","C","B"))
ONASm2 = lmer(resp_rate~Treatment + (1|Individual), data = ONASdf2)
plot(ONASm2)
summary(ONASm2)

lynxdf2 = fdata %>% filter(Species == "Lynx") %>%
  filter(!resp_rate < -1 & !r2 <.25)
lynxdf2$Treatment <- factor(lynxdf$Treatment, levels = c("N","C","B"))
lynxm2 = lmer(resp_rate~Treatment + (1|Individual), data = lynxdf2)
plot(lynxm2)
summary(lynxm2)

#### Figures ####
# to-do: for final figure, overlay mean data points for each treatment
minID = fdata %>% group_by(Species, Individual) %>%
  summarize(Start = min(Start)) %>%
  rename(Start0 = Start)

fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
  group_by(Species, Individual, Treatment) %>%
  summarize(Start = min(Start)) %>%
  left_join(sumfdata) %>%
  left_join(minID) %>%
  mutate(Start = Start -Start0) %>%
  ggplot(aes(x=Start, y=resp_rate, shape=Treatment)) + 
  geom_point(size=1) + geom_line(aes(group=Individual)) + theme_classic() +
  facet_grid(.~Species)

fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
  group_by(Species, Individual, Treatment) %>%
  ggplot(aes(x=Treatment, y=resp_rate, group=Individual)) +
  geom_line(data = fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
              group_by(Species, 
                       Treatment, 
                       Individual) %>%
                    summarise(AvRR = mean(resp_rate)),
             aes(x = Treatment, y = AvRR, group=Individual)) + 
  geom_point(data = fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
               group_by(Species, 
                        Treatment, 
                        Individual) %>%
               summarise(AvRR = mean(resp_rate)),
             aes(x = Treatment, y = AvRR, group=Individual)) + 
  theme_classic() +
  labs(y = "uL CO2 g-1 min-1") +
  facet_grid(.~Species)

 ggsave("Fear_results_21Feb2020.png", plot = last_plot())
 

### **AQUATIC** ----
 ### CLEANING ----
 ### Main code written by A. Arietta. Updated by Y. Alshwairikh
 ### Last update 2020 Mar 16 by Andis
 
 library(marelac)
 library(tidyverse)
 library(lubridate)
 library(ggalt)
 library(zoo)
 
 setwd("C://Users//Andis//Desktop//RworkingProjects//RASY_PredationResp//JulXX")
 
 RespDataPrep <- function(ExperimentData){
   
   X <- read_tsv(ExperimentData, skip = 38, col_names = FALSE, cols_only(X1 = "?", X2 = "c", X3 = "n", X4 = "n", X6 = "n", X7 = "n", X9 = "n", X10 = "n", X12 = "n", X13 = "n", X15 = "n", X16 = "n", X18 = "n", X19 = "n", X21 = "n", X22 = "n", X24 = "n", X25 = "n", X27 = "n", X28 = "n")) # Read in the FILENAME_raw.txt file and extract experiment data
   colnames(X) <- c("DateTime", "Phase", "Baro", "Sal", "CH1.temp", "CH1.O2persat", "CH2.temp", "CH2.O2persat", "CH3.temp", "CH3.O2persat", "CH4.temp", "CH4.O2persat", "CH5.temp", "CH5.O2persat", "CH6.temp", "CH6.O2persat", "CH7.temp", "CH7.O2persat", "CH8.temp", "CH8.O2persat") # Rename the columns
   
   AutoResp.convert.DateTime <- function(i){
     Xdates <- (strsplit(as.character(i), "/"))
     Xdates <- matrix(unlist(Xdates), ncol=4, byrow=TRUE)
     Xdates <- as.data.frame(Xdates)
     Date <- paste(Xdates[,3], Xdates[,1], Xdates[,2], sep = "-")
     Time <- as.character(Xdates$V4)
     datetime <- paste(Date, Time, sep = "T")
     DT <- as_datetime(datetime)
     DT} # This function converts the original datetime column from the text file into lubridate datetime format
   
   AutoResp.combine.channels <- function(data, colTemp, colO2, ch){
     data %>% 
       select(DateTime:Sal, "Temp" = colTemp, "O.2persat_raw" = colO2) %>% 
       mutate(CH = ch) %>% 
       mutate(DateTime = AutoResp.convert.DateTime(DateTime)) %>% 
       mutate(sec = ((hour(DateTime)*3600) + (minute(DateTime)*60) + (second(DateTime)))) %>% 
       mutate(O.2persat = rollmean(O.2persat_raw, k = 5, align = "right", fill = NA))
   } # This function reformats the temp and O2 columns from the wide-format raw text file so that they can be combined into a tidy-format dataframe in the next step.
   # It also converts the DateTime into a variable "sec" which is the second from midnight
   # It also computes a rolling average of raw O2 values to reduce noise using a window value of 5 seconds and aligning the new data to the right retains only averages of measurement data to the end of the measurement period.
   
   Xfull <- AutoResp.combine.channels(X, "CH1.temp", "CH1.O2persat", 1) %>%
     bind_rows(AutoResp.combine.channels(X, "CH2.temp", "CH2.O2persat", 2)) %>%
     bind_rows(AutoResp.combine.channels(X, "CH3.temp", "CH3.O2persat", 3)) %>%
     bind_rows(AutoResp.combine.channels(X, "CH4.temp", "CH4.O2persat", 4)) %>%
     bind_rows(AutoResp.combine.channels(X, "CH5.temp", "CH5.O2persat", 5)) %>%
     bind_rows(AutoResp.combine.channels(X, "CH6.temp", "CH6.O2persat", 6)) %>%
     bind_rows(AutoResp.combine.channels(X, "CH7.temp", "CH7.O2persat", 7)) %>%
     bind_rows(AutoResp.combine.channels(X, "CH8.temp", "CH8.O2persat", 8))# This combines all of the data from each chamer into a tidy-format table.
   
   Xfull$Resp.DOY <- yday(Xfull$DateTime) # Creates a column of DOY date in order to combine with larval data
   
   Xfull$O2.mg_per_L <- marelac::gas_satconc(S = Xfull$Sal, t = Xfull$Temp, P = Xfull$Baro/1000, species = 'O2') * 1e-6 * marelac::molweight('O2') * 1e3 * Xfull$O.2persat/100 # Converts O2 readings from percent saturation to mg/L using solubility conversions from the 'marelac' package
   
   Xinfo <- read_tsv(ExperimentData, skip = 25, n_max = 8, col_names = FALSE, cols_only(X1 = "c", X2 = "n", X3 = "n", X4 = "n")) # Reads in the pertinent portion of the info data output form AutoResp
   colnames(Xinfo) <- c("CH", "CHvol.mL", "Tubevol.mL", "Mass.g") # Rename the columns
   Xinfo$CH <- seq(1:length(Xinfo$CH)) # Adds a numerical variable for the channel numer
   
   X <- left_join(Xfull, Xinfo, by = "CH") # Joins the experiment data and info
   
   X <- filter(X, CH != 5) # Remove CH 5 that was not used in the experiment
 } # Preps the raw data
 RespDataSummary <- function(X, Type){
   tab <- c()
   for(i in unique(X$CH)){
     for(j in unique(filter(X, grepl("M", Phase))$Phase)){
       dat <- filter(X, CH == i, Phase == j) %>%
         filter((sec - sec[1]) > 200) # This bit excludes the first 200 seconds of the measurement phase.
       lm1 <- lm(O2.mg_per_L ~ sec, data = dat)
       sum1 <- summary(lm1)
       lo1 <- loess(O2.mg_per_L ~ sec, data = dat, span = 0.5)
       sum2 <- summary(lo1)
       slope <- sum1$coefficients[2,1]
       R2 <- sum1$r.squared
       df <- data.frame("CH" = i, "Phase" = j, "Slope" = slope, "R2" = R2, "RMSE.lm" = sqrt(mean((sum1$residuals)^2)), "RMSE.lo" = sqrt(mean((sum2$residuals)^2)), "Mass.g" = mean(dat$Mass.g), "Temp.C" = mean(dat$Temp), sec = mean(dat$sec), CHvol.mL = mean(dat$CHvol.mL), Tubevol.mL = mean(dat$Tubevol.mL), Resp.DOY = mean(dat$Resp.DOY))
       tab <- bind_rows(tab, df)
     }
   }
   tab$Type <- Type # Adds a variable for the run type (i.e. "SMR", "Treat", "Pre" test, "Post" test)
   
   tab$O2.mg_h <- (((tab$Slope * -1) * (((tab$CHvol.mL+tab$Tubevol.mL) - tab$Mass.g) / 1000)) * 3600) # Converts O2 reduction into MO2 by scaling to the water volume of the chamber.
   
   tab
 } # Creates summary of measurements
 RespPrepPlot1 <- function(RespData){
   ggplot(RespData, aes(x = sec, y = O2.mg_per_L, color = as.factor(CH))) +
     geom_line() +
     theme_minimal()
 } # Visualize the entirety of the run for quality control
 RespPrepPlot2 <- function(RespData){
   RespData %>% filter(grepl('M', Phase)) %>%
     ggplot(aes(x = sec, y = O2.mg_per_L, color = as.factor(CH))) +
     geom_line() +
     geom_smooth(se = FALSE, col = 1) +
     facet_grid(vars(CH), vars(Phase), scales = "free") +
     theme_minimal()
 } # Visualize the measurement slopes
 Resp.InputQC <- function(RespSum){
   RespSum %>% group_by(CH) %>% summarise(Mass = mean(Mass.g), Baro = mean(Baro), Temp = mean(Temp), Resp.DOY = mean(Resp.DOY), Sal = mean(Sal), CHvol = mean(CHvol.mL), Tubevol = mean(Tubevol.mL))} # Summary of input variables to check for errors.
 RespSumPlot1 <- function(RespSum){
   RespSum %>%
     ggplot(aes(x = Slope, y = RMSE.lm - RMSE.lo, col = Phase))+
     geom_point(size = 4) +
     geom_encircle() +
     theme_minimal()
 }
 
 
 # Define the run of interest:
 Pre <- RespDataPrep(  "201907XX_PredatorExp_blankPre_raw.txt")
 Acc <- RespDataPrep(  "201907XX_PredatorExp_acc_raw.txt")
 SMR <- RespDataPrep(  "201907XX_PredatorExp_smr_None_raw.txt")
 Treat <- RespDataPrep("201907XX_PredatorExp_smr_Treat_raw.txt")
 Post <- RespDataPrep( "201907XX_PredatorExp_blankPost_raw.txt")
 
 RespPrepPlot1(Pre) # Quality control for raw data
 RespPrepPlot2(Pre) # Quality control for all slopes
 Resp.InputQC(Pre)
 
 RespPrepPlot1(Post) # Quality control for raw data
 RespPrepPlot2(Post) # Quality control for all slopes
 Resp.InputQC(Post)
 
 RespPrepPlot1(Acc) # Quality control for raw data
 RespPrepPlot2(Acc) # Quality control for all slopes
 Resp.InputQC(Acc)
 
 RespPrepPlot1(SMR) # Quality control for raw data
 RespPrepPlot2(SMR) # Quality control for all slopes
 Resp.InputQC(SMR)
 
 RespPrepPlot1(Treat) # Quality control for raw data
 RespPrepPlot2(Treat) # Quality control for all slopes
 Resp.InputQC(Treat)
 
 # Create a summary of the measurement values
 SumPre <- RespDataSummary(Pre, "Pre")
 SumAcc <- RespDataSummary(Acc, "Acc")
 SumSMR <- RespDataSummary(SMR, "SMR")
 SumTreat <- RespDataSummary(Treat, "Treat")
 SumPost <- RespDataSummary(Post, "Post")
 
 RespSumPlot1(SumPre) # Quality control for slope values
 RespSumPlot1(SumAcc) # Quality control for slope values
 RespSumPlot1(SumSMR) # Quality control for slope values
 RespSumPlot1(SumTreat) # Quality control for slope values
 RespSumPlot1(SumPost) # Quality control for slope values
 
 # Define background respiration ====
 ## Background respiration accumulates over time as the animals are in the baths. We measure the background respiration rate in the clean chambers before introducing larvae and again at the end of the experiment. We assume a linear increase in background respiration over time and calculate the change in background respiration over the time span of the experiment. This allows us to estimate the amount of 02 consumption at each timepoint due to background respiration and can discount this from the larval MO2 estimates.
 
 SumPre <- SumPre %>% filter(Phase != "M1") %>% group_by(Type, CH) %>% summarise(avMO2 = mean(O2.mg_h), minMO2 = min(O2.mg_h), maxMO2 = max(O2.mg_h), avR2 = mean(R2), minR2= min(R2), maxR2 = max(R2), Temp.C = mean(Temp.C), sec = mean(sec))
 
 SumPost <- SumPost %>% filter(Phase != "M1") %>% group_by(Type, CH) %>% summarise(avMO2 = mean(O2.mg_h), minMO2 = min(O2.mg_h), maxMO2 = max(O2.mg_h), avR2 = mean(R2), minR2= min(R2), maxR2 = max(R2), Temp.C = mean(Temp.C), sec = mean(sec))
 
 (SumPost$avMO2 - SumPre$avMO2) / (SumPost$sec - SumPre$sec)
 
 BG.corrrection <- left_join(
   (SumPre %>% group_by() %>% select(CH, avMO2, sec) %>% rename(PreMO2 = avMO2, PreSec = sec)), 
   (SumPost %>% group_by() %>% select(CH, avMO2, sec) %>% rename(PostMO2 = avMO2, PostSec = sec)), by = "CH") %>%
   mutate(MO2.Cor = (PostMO2 - PreMO2)/(PostSec - PreSec)) %>%
   select(CH, PreSec, PostSec, MO2.Cor)
 
 (BG.corrrection$PostSec - BG.corrrection$PreSec) * BG.corrrection$MO2.Cor # In order to correct for the background respiration, I need to subtract the time in sec of the phase to get the number of seconds since the beginning of the experiment and then multiply by the correction factor.
 
 # Define experimental respiration rates =====
 
 # . No Treatment analysis: =====
 
 SumSMR <- left_join(SumSMR, BG.corrrection, by = "CH") %>% mutate(corO2.mg_h = O2.mg_h - ((sec - PreSec) * MO2.Cor), corO2.mg_h.g = corO2.mg_h / Mass.g)
 
 # Check for influence of temperature over time:
 SumSMR %>% filter(Phase != "M1") %>%
   ggplot(aes(x = Temp.C, corO2.mg_h, color = as.factor(CH))) +
   geom_point()+
   geom_smooth(method = 'loess', se = FALSE)
 summary(lm(corO2.mg_h ~ Temp.C * as.factor(CH), data = SumSMR))
 
 # Check for low quality measurements 
 SumSMR %>% filter(Phase != "M1") %>%
   ggplot(aes(x = CH, y = corO2.mg_h.g)) +
   geom_jitter(width = 0.1, aes(color = ifelse(RMSE.lm - RMSE.lo >= 0.01, 2, 3))) +
   scale_color_identity()
 
 FinalSMR <- SumSMR  %>% select(Resp.DOY, Type, CH, Phase, corO2.mg_h, corO2.mg_h.g, R2, Mass.g, Temp.C, sec, CHvol.mL, Tubevol.mL)
 
 # . Treatment SMR analysis: =====
 
 SumTreat <- left_join(SumTreat, BG.corrrection, by = "CH") %>% mutate(corO2.mg_h = O2.mg_h - ((sec - PreSec) * MO2.Cor), corO2.mg_h.g = corO2.mg_h / Mass.g)
 
 # Check for influence of temperature over time:
 SumTreat %>% filter(Phase != "M1") %>%
   ggplot(aes(x = Temp.C, corO2.mg_h, color = as.factor(CH))) +
   geom_point()+
   geom_smooth(method = 'loess', se = FALSE)
 summary(lm(corO2.mg_h ~ Temp.C * CH, data = SumTreat))
 
 # Check for low quality measurements
 SumTreat %>% filter(Phase != "M1") %>%
   ggplot(aes(x = CH, y = corO2.mg_h.g)) +
   geom_jitter(width = 0.2, aes(color = ifelse(RMSE.lm - RMSE.lo >= 0.01, 2, 3))) +
   scale_color_identity()
 
 FinalTreat <- SumTreat  %>% select(Resp.DOY, Type, CH, Phase, corO2.mg_h, corO2.mg_h.g, R2, Mass.g, Temp.C, sec, CHvol.mL, Tubevol.mL)
 
 # . Combined analysis: ====
 
 Final <- bind_rows(SumSMR, SumTreat)
 
 ggplot(Final, aes(x = CH, y = corO2.mg_h.g, col = Type)) +
   geom_jitter(width = 0.2) +
   theme_minimal()
 
 summary(lm(corO2.mg_h.g ~ Type * as.factor(CH), data = Final))
 
 write.csv(Final, "JulXX.csv")
 #Repeat for the other experiment days
 
#### COMPILE ----
 ### Main code written by A. Arietta. Updated by Y. Alshwairikh
 ###Last update 02.20.2020
 
 library(lme4)
 library(lmerTest)
 library(tidyverse)
 theme_set(theme_minimal())
 
 setwd("C://Users//Andis//Desktop//RworkingProjects//RASY_PredationResp//data")
 
 Jul11 <- read.csv("raw/Jul11/Jul11.csv")
 Jul12 <- read.csv("raw/Jul12/Jul12.csv")
 Jul15 <- read.csv("raw/Jul15/Jul15.csv")
 Jul16 <- read.csv("raw/Jul16/Jul16.csv")
 Jul17 <- read.csv("raw/Jul17/Jul17.csv")
 Jul18 <- read.csv("raw/Jul18/Jul18.csv")
 Jul19 <- read.csv("raw/Jul19/Jul19.csv")
 Jul22 <- read.csv("raw/Jul22/Jul22.csv")
 Jul23 <- read.csv("raw/Jul23/Jul23.csv")
 
 X <- rbind(Jul11, Jul12, Jul15, Jul16, Jul17, Jul18, Jul19, Jul22, Jul23)
 #bind the clean data into one master file
 
 ## Analyse and correct for temperature ====
 
 X %>% group_by(Type, Resp.DOY, CH) %>% summarise(Temp.C = mean(Temp.C)) %>% group_by(Type) %>% summarise(Temp.C = mean(Temp.C))
 # There is about a 0.2 degree difference between SMR and Treatment.
 
 X %>% filter(Type == "SMR" | Type == "Treat") %>%
   ggplot(aes( x = Temp.C, y = corO2.mg_h.g, col = Resp.DOY, group = Resp.DOY)) +
   geom_point() +
   geom_smooth(method = 'lm') +
   geom_smooth(se = F) +
   facet_grid(cols = vars(Type), rows = vars(Resp.DOY), scales = "free")
 # Visualize the temp change over the course of the run. It doesn't seem like fluctuations in temperature have a huge impact on the metabolic rates.
 
 summary(lmer(corO2.mg_h.g ~ Temp.C + (1 | Resp.DOY/Type), data = X)) # There is no evidence that the range of temperatures we considered impacted metabolic rates.
 
 summary(lmer(corO2.mg_h.g ~ Type + (1 | Resp.DOY), data = X)) # There is also no evidence that the range of temperatures we considered impacted metabolic rates within a trial between treatments.
 
 X %>% group_by(Type, Resp.DOY) %>% summarise(maxT = max(Temp.C), minT = min(Temp.C), difT = maxT - minT) %>% select(-Resp.DOY) %>% group_by() %>% summarise(mean(difT), sd(difT))
 
 X %>% group_by(Type, Resp.DOY) %>% summarise(maxT = max(Temp.C), minT = min(Temp.C), difT = maxT - minT) %>% select(-Resp.DOY) %>% group_by(Type) %>% summarise(mean(difT), sd(difT)) # The temps varaied by 0.35 degrees (sd = 0.08) for SMR and 0.25 degrees (sd = 0.06) for Treatments
 
 X %>% group_by(Type, Resp.DOY) %>% summarise(meanT = mean(Temp.C)) %>%
   spread(Type, meanT) %>%
   mutate(difT = Treat - SMR) %>%
   summarise(mean(difT), sd(difT)) # Between SMR and Treatments the temps were always warmer. On average, the Treatment trials were 0.18 degrees warmer (sd = 0.12)
 
 X %>% group_by(Type, Resp.DOY) %>%
   ggplot(aes(x = as.factor(Resp.DOY), y = Temp.C, fill = Type)) +
   geom_boxplot(aes(group = interaction(Type, Resp.DOY)), alpha = 0.8, size = 1.1) +
   scale_fill_manual(values = c("dodger blue", "orange red")) +
   theme_minimal()
 
 hist(-X$MO2.Cor)
 max(-X$MO2.Cor)
 min(-X$MO2.Cor)
 mean(-X$MO2.Cor)
 
 hist(X$Temp.C)
 
 ggplot(X, aes(x = Temp.C, fill = Type)) +
   geom_density(col = NA, alpha = 0.8) +
   scale_fill_manual(values = c("dodger blue", "orange red")) +
   theme_minimal()
 
 ## . Analyze metabolic rate change by temperature. ====
 
 ggplot(X, aes(x = Temp.C, corO2.mg_h.g)) +
   geom_point() +
   geom_smooth(method = 'lm')
 
 summary(lmer(corO2.mg_h.g ~ Temp.C + Mass.g + (1 | CH) + (1 | Type), data = X)) # There is no strong effect of temperature, so we can just account for this in the model.
 
 X <- X %>% mutate(Ex = ifelse(CH == 1 | CH == 2 | CH == 3 | CH == 4, "DF", "CUE")) %>% 
   rename(MO2 = corO2.mg_h.g)
 
 write.csv(X, "processed/masterData_PredResp.csv")
 #master data with assignment of treatment condition (Ex= DF or CUE) and MO2 adjusted for temperature
 
 ####
 ####
 ####
 
 # OLD CODE
 
 ####
 ####
 ####
 
 # We can look within our data to estimate the effect of temperature on oxygen consumption rates by seeing how temperature changes MO2 within each treatment. We also are respective to individual differences using random intercepts for each channel and treatment.
 # Although temperature does not have a significant effect, we can still use the estimated value of 0.00639 units of MO2 increase for every 1 degress temp increase to correct between the trials.
 
 ## . . Correct the values for temperature. ====
 X <- X %>% left_join(
   X %>% group_by(Type, Resp.DOY) %>% summarise(meanT = mean(Temp.C)) %>%
     spread(Type, meanT) %>%
     mutate(difT = Treat - SMR) %>%
     mutate(T.cor.factor = -difT * 0.00639), by = "Resp.DOY") %>%
   mutate(MO2 = ifelse(Type == "Treat", corO2.mg_h.g + T.cor.factor, corO2.mg_h.g))
 
 X %>% group_by(Type, Resp.DOY) %>%
   ggplot(aes(x = as.factor(Resp.DOY), y = MO2, fill = Type)) +
   geom_boxplot(aes(group = interaction(Type, Resp.DOY)), alpha = 0.8, size = 1.1) +
   scale_fill_manual(values = c("dodger blue", "orange red")) +
   theme_minimal() # These are all of the corrected MO2 values, but we still need to remove outliers and average across individuals.
 
 
 ## . Analyze metabolic rate change by temperature. ====
 ## Originally, I had thought we might use the values from a different experiment to estimate a temperature correction, but I'm not sure if mixing data from a different experiment with different animals is appropriate.
 
 ExpMetRate <- read.csv("raw/2019_ExpMetabRates20190801.csv", header = TRUE)
 ExpMetRate %>% filter(Mass > 0) %>%
   ggplot(aes(x = Temp.C, y = O2.mg_h.g)) + geom_point() + geom_smooth(method = 'lm')
 
 summary(lm(O2.mg_h.g ~ Temp.C, data = filter(ExpMetRate, Mass > 0)))
 summary(lm(O2.mg_h.g ~ Temp.C + GS + CH + Mass, data = filter(ExpMetRate, Mass > 0)))
 # We know from our other (no predators) study that respiration should increase by 0.0539 mgO2/h*g per degree increase in Temp.
 # So, our null hypothesis should be that the Treatment measurments will be (0.0539 * difference in mean temperature) units higher than the SMR.
 
 ## . . Correct the values for temperature. ====
 X <- X %>% left_join(
   X %>% group_by(Type, Resp.DOY) %>% summarise(meanT = mean(Temp.C)) %>%
     spread(Type, meanT) %>%
     mutate(difT = Treat - SMR) %>%
     mutate(T.cor.factor = -difT * 0.0539), by = "Resp.DOY") %>%
   mutate(MO2 = ifelse(Type == "Treat", corO2.mg_h.g + T.cor.factor, corO2.mg_h.g))
 
 X %>% group_by(Type, Resp.DOY) %>%
   ggplot(aes(x = as.factor(Resp.DOY), y = MO2, fill = Type)) +
   geom_boxplot(aes(group = interaction(Type, Resp.DOY)), alpha = 0.8, size = 1.1) +
   scale_fill_manual(values = c("dodger blue", "orange red")) +
   theme_minimal() # These are all of the corrected MO2 values, but we still need to remove outliers and average across individuals.
 
 X <- X %>% mutate(Ex = ifelse(CH == 1 | CH == 2 | CH == 3 | CH == 4, "DF", "CUE"))
 
 write.csv(X, "processed/masterData_PredResp.csv")
 #master data with assignment of treatment condition (Ex= DF or CUE) and MO2 adjusted for temperature
 
 
### ANALYSIS ----
 library(lmer)
 library(parameters)
 library(tidyverse)
 library(MuMIn)
 
 ## Read in the comiled data ====
 
 X <- read.csv("processed/masterData_PredResp.csv") %>%
   mutate(RMSE.dif = RMSE.lm - RMSE.lo) %>% # Calculate a metric of linearity from the difference in loess and linear model fits
   mutate(Tx = ifelse(Type == "SMR", as.character(Type), as.character(Ex))) %>% # Make a new factor that combines the treatment levels and SMR baseline
   mutate(Tx = fct_relevel(Tx, "SMR", "CUE", "DF")) # Reorder the factors so that they make sense in the figures.
 
 ## Refine the dataset ====
 # The main goal here is to make sure that we are only using the true minimum resting MO2 and excluding measurement phases when there are either problems with the sensor readings or the animal has increase respiration.
 
 X %>% ggplot(aes(x = RMSE.dif)) +
   geom_histogram() # These values look pretty good. Most of the loess fits match up to the linear fits really well.
 
 RMSE.dif.threshold <- as.numeric(quantile((X$RMSE.dif), .95)[1]) # 95% of the observations fall within 0.009 difference between the linear and loess regression fits, so we will exclude any observations with fits above this value
 
 X %>% ggplot(aes(x = MO2, y = RMSE.dif)) +
   geom_jitter(aes(col = ifelse(RMSE.dif >= RMSE.dif.threshold, 2, 1))) +
   scale_color_identity()
 
 X %>% filter(RMSE.dif < RMSE.dif.threshold & MO2 > 0) %>%
   ggplot(aes(x = MO2)) +
   geom_histogram(binwidth = 0.001) # This looks like pretty reasonable distribution
 
 ## . Remove outliers ====
 X <- X %>% filter(RMSE.dif < RMSE.dif.threshold & MO2 > 0) # We remove all of the high outliers with RSME differeces above 0.009
 
 ## . Restrict to the lower half of observations for each individual within a trial =====
 # Given this, we'll average the lowest 50% of the data for each individual.
 X2 <- X %>% 
   left_join(X %>% group_by(Resp.DOY, Type, CH) %>% summarise(med.MO2 = median(MO2)), by = c("Resp.DOY", "Type", "CH")) %>%
   filter(MO2 <= med.MO2) # Now, we take the lowest 50% of the values for each channel in each trial by treatment
 
 X2 %>% group_by(Type, Resp.DOY) %>%
   ggplot(aes(x = as.factor(Resp.DOY), y = MO2, fill = Type)) +
   geom_boxplot(aes(group = interaction(Type, Resp.DOY)), alpha = 0.8, size = 1.1) +
   scale_fill_manual(values = c("dodger blue", "orange red")) +
   theme_minimal() # Plot of baseline versus either treatment.
 
 X2 %>% group_by(Tx, Resp.DOY) %>%
   ggplot(aes(x = as.factor(Resp.DOY), y = MO2, fill = Tx)) +
   geom_boxplot(aes(group = interaction(Tx, Resp.DOY)), alpha = 0.8, size = 1.1) +
   scale_fill_manual(values = c("dodger blue", "orange red", "dark green")) +
   theme_minimal()
 
 ## Fit the model ====
 
 mod1 <- lmer(MO2 ~ Tx + Temp.C + (1 | Resp.DOY/CH), data = X2)
 mod1.sum <- summary(mod1)
 mod1.sum # Temperature has no strong directional effect but we will leave it in to account for it. There is no strong difference between the baseline and olfactory cue, but the olfactory plus visual cue is lower.
 
 ## . Confidence intervals ====
 # This will take a while to run
 
 mod1.ci <- model_parameters(mod1, ci = 0.95, bootstrap = TRUE, iterations = 1000)
 mod1.ci
 mod1.ci$Coefficient # The values are small so we need to output them explicitly
 mod1.ci$CI_low
 mod1.ci$CI_high
 
 ## . R2 values ====
 r.squaredGLMM(mod1)
 
 ## Plot the figure ====
 
 X2.2 <- X2 %>% 
   mutate(fit.re = predict(mod1, re.form = NULL)) %>% 
   group_by(Tx, Resp.DOY) %>% 
   summarise(MO2 = mean(fit.re)) # This first predicts the value for each trial based on the random effect (essentially an individual-level average). These values are the same for all phases and individuals so we just need to average them down to a single value to get our 9 trial estimates for plotting.
 
 X2.3 <- X2 %>% 
   mutate(fit = predict(mod1, re.form = ~0)) %>% 
   group_by(Tx) %>% 
   summarise(MO2 = mean(fit) + mod1.sum$sigma) # This predicts the value excluding the random effects. These values are the same for all phases and individuals and trials, so again, we just need to average to collapse these into single values. The predictions actually EXCLUDE the random effects rather than average them, so we add the residual sigma from the model to shift the values up to the average for plotting.
 
 ggplot(X2.2, aes(x = Tx, y = MO2)) +
   geom_line(aes(group = Resp.DOY), size = 1, col = "grey40") +
   geom_line(data = X2.3, aes(x = Tx, y = MO2, group = 1), col = "orange red", size = 2) # Need to change the labels in illustrator
 
 ####
 ####
 ####
 # OLD CODE
 ####
 ####
 ####
 
 ## This is just an alternate way to plot the figure calculating means by hand (not the ideal method)
 avSMR <- X2 %>% select(Resp.DOY, Type, Ex, CH, Phase, Mass.g, Temp.C, MO2) %>%
   filter(Type == "SMR") %>% group_by(Resp.DOY, CH) %>%
   summarise(MO2 = mean(MO2)) %>% group_by(Resp.DOY) %>%
   summarise(SMR = mean(MO2))
 
 X2.1 <- X2 %>% select(Resp.DOY, Type, Ex, CH, Phase, Mass.g, Temp.C, MO2) %>%
   group_by(Resp.DOY, Type, Ex, CH) %>%
   summarise(MO2 = mean(MO2)) %>%
   group_by() %>%
   spread(Type, MO2) %>%
   mutate(Treat.dif = Treat - SMR) %>%
   group_by(Resp.DOY, Ex) %>%
   summarise(Treat.dif = mean(Treat.dif)) %>%
   group_by() %>%
   spread(Ex, Treat.dif) %>%
   left_join(avSMR, by = "Resp.DOY") %>%
   mutate(CUE = SMR + CUE, DF = SMR + DF) %>%
   gather(SMR, CUE, DF, key = "Tx", value = MO2) %>%
   mutate(Tx = fct_relevel(Tx, "SMR", "CUE", "DF"))
 
 ggplot(X2.1, aes(x = Tx, y = MO2)) +
   geom_line(aes(group = Resp.DOY), size = 1, col = "grey40") +
   geom_line(data = X2.3, aes(x = Tx, y = MO2, group = 1), col = "orange red", size = 2)
 
 # Difference in mean of trucated MO2. ====
 X3 <- X2 %>% group_by(Resp.DOY, Ex, CH) %>% select(Type, MeanMO2.g) %>% spread(Type, MeanMO2.g)
 
 X5 <- X2 %>%
   mutate(avg.temp = mean(Temp.C))
 
 summary(lm((Treat - SMR) ~ Ex, data = X3))
 
 
 ggplot(X2, aes(x = MeanMO2.g, fill = Type)) +
   geom_density(adjust = 1, col = NA, alpha = .3) +
   geom_vline(xintercept = 0, lty = 2) +
   theme_minimal() # Exposure to predator had no impact on respiration rates.
 summary(lm(MeanMO2.g ~ Type, data = X2))
 #This is linear model of average O2 vs. type (basically is there difference
 #in meanO2 between SMR and treatment phase?) #no there's not
 
 summary(lmer(MeanMO2.g ~ Type + Temp.C + (1|Resp.DOY/CH), data = X2))
 #Here we use a MM with temp as fixed effect, and day/ch as random effect
 #We expect every CH within day will have its own intercept
 #again, no difference in meanO2 in SMR vs. treatment
 
 ggplot(X2, aes(x = Type, y = MeanMO2.g, col = Ex)) +
   geom_point() +
   stat_summary(aes(group = Ex), geom = "line", fun.y = mean, size = 2) +
   theme_minimal() # There is no huge difference in the type of exposure.
 summary(lm(MeanMO2.g ~ Type * Ex, data = X2))
 summary(lmer(MeanMO2.g ~ Type * Ex + Temp.C + (1|Resp.DOY/CH), data = X2)) 
 # However there is an interesting difference in that adding the visual 
 #stimulus decreases respiration relative to the chemical cue alone (p = 0.2).
 
 #Start YA
 hist(X2$MeanMO2.g)
 hist(X$MO2.mg_h.g)
 
 X4 <- X2 %>% 
   filter(Type == "Treat")
 
 summary(lm(MeanMO2.g ~ Ex + Temp.C, data = X4))
 summary(lmer(MeanMO2.g ~ Ex + Temp.C + (1|Resp.DOY/CH), data = X4))
 
 summary(lmer(MeanMO2.g ~ Ex + Temp.C + (1|indiv), data = X4 %>% mutate(indiv = paste0(Resp.DOY,"ch",CH))))
 
 X4 <- X4 %>% mutate(indiv = paste0(Resp.DOY,"ch",CH))
 
 X4 %>% filter(Resp.DOY == 193)
 
 X4 %>% group_by(Resp.DOY, Type, Ex) %>% summarise_all(list(mean)) %>%
   lm( ~ Ex + Temp.C, data = .,) %>%
   summary(.)
 
 X4 %>% ggplot(aes(x = Ex, y = dif))
 #End YA
 
 ggplot(X3, aes(x = (Treat - SMR), fill = Ex)) +
   geom_density(adjust = 1, col = NA, alpha = .3) +
   geom_vline(xintercept = 0, lty = 2) +
   theme_minimal()
 
 p1 <- ggplot(X3, aes(x = Ex, y = (Treat - SMR))) +
   geom_boxplot() +
   geom_jitter(width = 0.1) +
   theme_minimal()
 p1 + labs(y = "mean MO2 Treat - mean MO2 SMR") + 
   scale_x_discrete(labels=c("Cue only", "Predator + Cue"))
 
 # Difference in median MO2. ====
 
 X4 <- X2 %>% group_by(Resp.DOY, Ex, CH) %>% select(Type, MedMO2.g) %>% spread(Type, MedMO2.g)
 summary(lm((SMR-Treat) ~ Ex, data = X4))
 
 ggplot(X4, aes(x = Ex, y = (Treat - SMR))) +
   geom_boxplot() +
   geom_jitter(width = 0.1) +
   theme_minimal()
 
 ggplot(X2, aes(x = Type, y = MedMO2.g, col = Ex)) +
   geom_point() +
   stat_summary(aes(group = Ex), geom = "line", fun.y = mean, size = 2) +
   theme_minimal()
 
 ggplot(X4, aes(x = (Treat - SMR), fill = Ex)) +
   geom_density(adjust = 1, col = NA, alpha = .3) +
   geom_vline(xintercept = 0, lty = 2) +
   theme_minimal()
 
 # Difference in mean of lowest 50% of truncated MO2 measurements
 
 X %>% group_by(Resp.DOY, Type, Ex, CH) %>% tally() %>% ggplot(aes(x = n)) + geom_histogram()
 
 X %>% filter(R2 >= 0.95) %>% group_by(Resp.DOY, Type, Ex, CH) %>% tally() %>% ggplot(aes(x = n)) + geom_histogram()
 
 X.QC.1 <- X %>% filter(R2 >= 0.95) %>% filter(rank(MO2.mg_h.g) != min(rank(MO2.mg_h.g))) %>% filter(rank(MO2.mg_h.g) != max(rank(MO2.mg_h.g))) %>% filter(rank(MO2.mg_h.g) <= (median(rank(X))+0.5)) %>% group_by(Resp.DOY, Type, Ex, CH) %>% tally() %>% arrange(n)
 
 X.QC.1 %>% ggplot(aes(x = n)) + geom_histogram()
 
 X.QC.2 <- X %>% filter(R2 >= 0.95) %>% # Excludes measurements with R2 below 0.95
   group_by(Resp.DOY, Type, Ex, CH) %>% filter(rank(MO2.mg_h.g) != min(rank(MO2.mg_h.g))) %>% filter(rank(MO2.mg_h.g) != max(rank(MO2.mg_h.g))) %>% # Truncates the data by removing the highest and lowest values.
   filter(rank(MO2.mg_h.g) <= (median(rank(X))+0.5)) %>% # Selects the lowest 50% of MO2 measurments for each individual in each treatment.
   group_by(Resp.DOY, Type, Ex, CH) %>% tally() %>% arrange(n)
 
 X.QC.2 %>% ggplot(aes(x = n)) + geom_histogram()
 
 X5 <- X %>% filter(R2 >= 0.95) %>% # Excludes measurements with R2 below 0.95
   group_by(Resp.DOY, Type, Ex, CH) %>% filter(rank(MO2.mg_h.g) != min(rank(MO2.mg_h.g))) %>% filter(rank(MO2.mg_h.g) != max(rank(MO2.mg_h.g))) %>% # Truncates the data by removing the highest and lowest values.
   filter(rank(MO2.mg_h.g) <= (median(rank(X))+0.5)) %>% # Selects the lowest 50% of MO2 measurments for each individual in each treatment.
   summarise(MeanMO2.g = mean(MO2.mg_h.g), Temp.C = mean(Temp.C), meanR2 = mean(R2), Mass.g = mean(Mass.g)) # Excludes the top and bottom values in the dataset then selects the lowest 50% of values.
 
 # View(left_join(X5, X.QC.2, by = c("Resp.DOY", "Type", "Ex", "CH")) %>% group_by(Ex) %>% arrange(desc(meanR2)))
 
 X5.1 <- X5 %>% group_by(Resp.DOY, Ex, CH) %>% select(Type, MeanMO2.g) %>% spread(Type, MeanMO2.g) %>% na.omit() %>% mutate(difMO2 = (Treat+0.0113) - SMR)
 
 X5.Temp <- X5 %>% group_by(Resp.DOY, Ex, CH) %>% select(Type, Temp.C) %>% spread(Type, Temp.C) %>% na.omit() %>% mutate(difTemp = Treat - SMR)
 
 X5.1 <- left_join(select(X5.1, Resp.DOY, Ex, CH, difMO2), select(X5.Temp, Resp.DOY, Ex, CH, difTemp), by = c("Resp.DOY", "Ex", "CH")) 
 
 # View(arrange(X5.1, desc(difMO2)))
 
 ggplot(X5.1, aes(x = difTemp, y = difMO2)) +
   geom_point() +
   geom_smooth()+
   theme_minimal()
 
 ggplot(X5.1, aes(x = Ex, y = difMO2)) +
   geom_boxplot() +
   geom_jitter(width = 0.1) +
   theme_minimal()
 
 summary(lm(difMO2 ~ Ex, data = X5.1))
 summary(lmer(difMO2 ~ Ex + difTemp + (1|CH), data = X5.1))
 
 
 ggplot(X5, aes(x = Type, y = MeanMO2.g, col = Ex)) +
   geom_point() +
   stat_summary(aes(group = Ex), geom = "line", fun.y = mean, size = 2) +
   theme_minimal()
 
 ggplot(X5.1, aes(x = difMO2, fill = Ex)) +
   geom_density(adjust = 1, col = NA, alpha = .3) +
   geom_vline(xintercept = 0, lty = 2) +
   theme_minimal()
 
 summary(lmer(difMO2 ~ Ex + (1|CH), data = X5.1))
 summary(lmer(difMO2 ~ Ex + difTemp + (1|CH), data = X5.1))
 # About a 1.5% reduction in MO2 with jus the Cue, and a 10% reduction with Cue + Dragon Fly
 X5 %>% group_by(Type) %>% summarise(grandmeanMO2.g = mean(MeanMO2.g, na.rm = TRUE))
 
 X5 %>% filter(Ex == "DF")
 
 X5.1 %>% group_by(Ex) %>% summarise(sd(difMO2))
 sd(X5.1$difMO2)