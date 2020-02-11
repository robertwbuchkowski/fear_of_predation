### Acute Fear Summer 2019 Proj
### N.R. Sommer & R.W. Buchkowski

# Andis has a min 10 obs per individual and removes observations in the lower quartile (ish) of R^2
# He also does this because they assume they're taking a baseline "resting" respiration rate.
# We should probably remove individuals that are obvious outliers, but correcting for R^2 would reduce the total number of observations down to one or two per individual. I would rather retain that variation than reduce the number of observations per individual.... 

# Andis suggested going back and looking at the raw
# CO2 ~ Time*Treatment + (Treatment | Individual)

require(tidyverse)
require(lubridate)

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

runeach_ADJ = function(XX = 1, IDTABLE = IDtable, RESPDATA = respdata) { 
  select = RESPDATA$Time > as.numeric(IDTABLE[XX,"Start"]) & 
    RESPDATA$Time < as.numeric(IDTABLE[XX,"End"]) &
    RESPDATA$Julian == as.numeric(IDTABLE[XX, "Julian"])
  
  selected = RESPDATA[select,]
  
  PRESSURE = IDTABLE[XX,"Pressure"]*25.4
  
  m1 = lm(CO2 ~ Time, data=selected)
  slope = coefficients(m1)[2]
  r2 = summary(m1)$adj.r.squared
  
  IDTABLE$Species = as.character(IDTABLE[XX, "Species"])
  if(IDTABLE$Species[XX] == "Blank") {
    slope_blank = coefficients(lm(CO2 ~ Time, data=selected))[2]
  }
  
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
  resprate = as.numeric((slope-slope_blank)*volume*60/(WEIGHT))
  c(r2 = r2, resprate = resprate)
}


output = cbind(IDtable, t(sapply(seq(1, NN, 1), FUN = runeach)))

for(kk in 1:NN){ # diagnose errors by line
  aaaa = runeach(kk)
  print(kk)
}

blank = output %>% filter(!(Treatment %in% c("B", "C", "N"))) %>%
  group_by(Julian) %>%
  summarize(blank = mean(resprate))

fdata = output %>% filter((Treatment %in% c("B", "C", "N"))) %>%
  left_join(blank) %>%
  mutate(resp_rate = resprate - blank)

sumfdata = fdata %>% group_by(Species, Individual, Treatment) %>%
  summarize(resp_rate = mean(resp_rate))

fdata %>% group_by(Species, Individual, Treatment) %>%
  summarize(resp_ratesd = sd(resp_rate)) %>%
  left_join(sumfdata) %>%
  mutate(cv = resp_ratesd/resp_rate)

dev.new()
dev.off()

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

library(lmerTest)
cricketm1 = lmer(resp_rate~Treatment + (1|Individual), data=cricketdf)
summary(cricketm1)
MEFEm1 = lmer(resp_rate~Treatment + (1|Individual), data = MEFEdf)
summary(MEFEm1) 
boxplot(MEFEdf$resp_rate ~ MEFEdf$Treatment)
PHIDm1 = lmer(resp_rate~Treatment + (1|Individual), data = PHIDdf)
summary(PHIDm1)
boxplot(PHIDdf$resp_rate ~ PHIDdf$Treatment)
lynxm1 = lmer(resp_rate~Treatment + (1|Individual), data = lynxdf)
summary(lynxm1)
ONASm1 = lmer(resp_rate~Treatment + (1|Individual), data = ONASdf)
summary(ONASm1)
TRRAm1 = lmer(resp_rate~Treatment + (1|Individual), data=TRRAdf)
summary(TRRAm1)

minID = fdata %>% group_by(Species, Individual) %>%
  summarize(Start = min(Start)) %>%
  rename(Start0 = Start)

fdata %>% group_by(Species, Individual, Treatment) %>%
  summarize(Start = min(Start)) %>%
  left_join(sumfdata) %>%
  left_join(minID) %>%
  mutate(Start = Start -Start0) %>%
  ggplot(aes(x=Start, y=resp_rate, shape=Treatment)) + 
  geom_point(size=1) + geom_line(aes(group=Individual)) + theme_classic() +
  facet_grid(.~Species)

fdata %>% group_by(Species, Individual, Treatment) %>%
  ggplot(aes(x=Treatment, y=resp_rate, group=Individual)) +
  geom_line(data = fdata %>% group_by(Species, 
                                      Treatment, 
                                      Individual) %>%
                    summarise(AvRR = mean(resp_rate)),
             aes(x = Treatment, y = AvRR, group=Individual)) + 
  geom_point(data = fdata %>% group_by(Species, 
                                       Treatment, 
                                       Individual) %>%
               summarise(AvRR = mean(resp_rate)),
             aes(x = Treatment, y = AvRR, group=Individual)) + 
  theme_classic() +
  labs(y = "uL CO2 g-1 min-1") +
  facet_grid(.~Species)

ggsave("Fear_results_Jan2020.png", plot = last_plot())

