### Acute Fear Summer 2019 Proj
### N.R. Sommer & R.W. Buchkowski

require(tidyverse)
require(lubridate)

#### READ-IN DATA ####
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

#### CALCULATE RESPIRATION RATES ####

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

### CLEANING ####

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

### MIXED EFFECTS MODELS ####
# to-do: use lme4 and bootstrap confidence values 
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

#### FIGURES ####
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