# CODE FOR MS: Prey functional types and physiological responses to acute predation risk

#### **TERRESTRIAL** ----
# Main code by R.W. Buchkowski & N.R. Sommer.
# Last update 10 Oct 2020 by N.R. Sommer

require(tidyverse)
require(lubridate)

#### READ-IN DATA ####

## Treatment code:
# N = no cue
# C = chemical cue
# B = both chemical & visual

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

### Function to calculate respiration rates ----
runeach = function(XX = 1, IDTABLE = IDtable, RESPDATA = respdata){ 
  select = RESPDATA$Time > as.numeric(IDTABLE[XX,"Start"]) & 
    RESPDATA$Time < as.numeric(IDTABLE[XX,"End"]) &
    RESPDATA$Julian == as.numeric(IDTABLE[XX, "Julian"])
  
  selected = RESPDATA[select,]
  
  PRESSURE = IDTABLE[XX,"Pressure"]*25.4
  
  m1 = lm(CO2 ~ Time, data=selected)
  slope = coefficients(m1)[2]
  r2 = summary(m1)$adj.r.squared
  TEMP = mean(selected$Temp)
  Weight = as.numeric(IDTABLE[XX,"WWt"])
  WEIGHT = ifelse(is.na(Weight), 1, Weight)
  WEIGHT0 = ifelse(is.na(Weight), 0, Weight)
  
  if(IDTABLE$Julian[XX] < 19236) { 
    vol = 60.79535394 # Volume by date, based on the tubing length of respirometer
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

### Calculate respiration rates ####

output = cbind(IDtable, t(sapply(seq(1, NN, 1), FUN = runeach)))

### CLEANING ----

# Diagnose errors by line, if any
## This takes time
for(kk in 1:NN){ 
  aaaa = runeach(kk)
  print(kk)
}

# Calculate and check blanks
blank = output %>% filter(Individual=="Blank") %>%
  group_by(Julian) %>%
  summarize(blank = mean(resprate))

blank %>% ggplot(aes(x=Julian, y = blank)) + geom_point() # looks good

# Join output and blanks
fdata = output %>% filter((Treatment %in% c("B", "C", "N"))) %>%
  left_join(blank)

# Check consistency of observation length 
fdata %>% mutate(End-Start == 200) %>% View() 

# Observation length incorrect for ONAS T, treatment N. End is 41360, start is 41150.
fdata$End[713] = 41350 # End time modified to 41350

### Respiration rates corrected by blanks ----

fdata = fdata %>% mutate(resp_rate = resprate - blank)

### Calculates mean of corrected CO2 (avg of 3 obs) ----
sumfdata = fdata %>% group_by(Species, Individual, Treatment) %>%
  summarize(resp_rate = mean(resp_rate))

fdata %>% group_by(Species, Individual, Treatment) %>%
  summarize(resp_ratesd = sd(resp_rate)) %>%
  left_join(sumfdata) %>%
  mutate(cv = resp_ratesd/resp_rate)

# Check R2 for respirations rates to get a sense of how linear CO2 accumulation was over measurement period
fdata %>% filter(Species == "Cricket") %>% 
  ggplot(aes(x=resp_rate, y = r2, color = Individual)) + geom_point() # good

fdata %>% filter(Species == "Lynx") %>% 
  ggplot(aes(x=resp_rate, y = r2, color = Individual)) + geom_point() # bad

fdata %>% filter(Species == "MEFE") %>% 
  ggplot(aes(x=resp_rate, y = r2, color = Individual)) + geom_point() # good

fdata %>% filter(Species == "ONAS") %>% 
  ggplot(aes(x=resp_rate, y = r2, color = Individual)) + geom_point() # outliers?

fdata %>% filter(Species == "Phiddipus") %>% 
  ggplot(aes(x=resp_rate, y = r2, color = Individual)) + geom_point() # outliers?

fdata %>% filter(Species == "TRRA") %>% 
  ggplot(aes(x=resp_rate, y = r2, color = Individual)) + geom_point() # one outlier

# Even with low R2 on CO2 accumulation, we have no reason to believe these were not real measurements.

### Calculate the average temperature variance and sd during trials ---

MM = dim(IDtable)[1]

tempparse = function(XX = 1, 
                     IDTABLE = IDtable, 
                     RESPDATA = respdata %>% select(Time, Temp, Julian)){ 
  select = RESPDATA$Time > as.numeric(IDTABLE[XX,"Start"]) & 
      RESPDATA$Time < as.numeric(IDTABLE[XX,"End"]) &
      RESPDATA$Julian == as.numeric(IDTABLE[XX, "Julian"])
  
  selected = RESPDATA[select,]
  
  tempvar = var(selected$Temp)
  tempsd = sd(selected$Temp)
    
  c(tempvar = tempvar, tempsd = tempsd)
}

tempparsed = cbind(IDtable, t(sapply(seq(1, MM, 1), FUN = tempparse)))
tempparsed = tempparsed %>% filter(!WWt == "NA") # remove blanks
mean(tempparsed$tempvar)
mean(tempparsed$tempsd)

### ANALYSIS ####

require(lme4)
require(parameters)
require(MuMIn)

## Renaming treatment codes for analyses & graphing
# N --> Baseline
# C --> Olfactory
# B --> Visual+Olfactory

sumfdata <- sumfdata %>%
  mutate(Treatment = replace(Treatment, Treatment == "N", "Baseline")) %>%
  mutate(Treatment = replace(Treatment, Treatment =="C", "Olfactory")) %>%
  mutate(Treatment = replace(Treatment, Treatment == "B", "Visual+Olfactory"))

fdata <- fdata %>%
  mutate(Treatment = replace(Treatment, Treatment == "N", "Baseline")) %>%
  mutate(Treatment = replace(Treatment, Treatment =="C", "Olfactory")) %>%
  mutate(Treatment = replace(Treatment, Treatment == "B", "Visual+Olfactory"))

# Filter data by species for mixed effects models 
cricketdf = sumfdata %>% filter(Species=="Cricket") 
MEFEdf = sumfdata %>% filter(Species== "MEFE")
PHIDdf = sumfdata %>% filter(Species=="Phiddipus")
lynxdf = sumfdata %>% filter(Species=="Lynx")
ONASdf = sumfdata %>% filter(Species=="ONAS")
TRRAdf = sumfdata %>% filter(Species=="TRRA")

cricketdf$Treatment <- factor(cricketdf$Treatment, 
                              levels = c("Baseline", "Olfactory", "Visual+Olfactory"))
MEFEdf$Treatment <- factor(MEFEdf$Treatment, 
                           levels = c("Baseline", "Olfactory", "Visual+Olfactory"))
PHIDdf$Treatment <- factor(PHIDdf$Treatment, 
                           levels = c("Baseline", "Olfactory", "Visual+Olfactory"))
lynxdf$Treatment  <- factor(lynxdf$Treatment, 
                            levels = c("Baseline", "Olfactory", "Visual+Olfactory"))
ONASdf$Treatment <- factor(ONASdf$Treatment, 
                           levels = c("Baseline", "Olfactory", "Visual+Olfactory"))
TRRAdf$Treatment <- factor(TRRAdf$Treatment, 
                           levels = c("Baseline", "Olfactory", "Visual+Olfactory"))

### Fit models by species ----

# CRICKETS #
cricketm1 = lmer(resp_rate~Treatment + (1|Individual), data=cricketdf)
plot(cricketm1)
summary(cricketm1)
r.squaredGLMM(cricketm1)

cricketm1.ci <- model_parameters(cricketm1, ci = 0.95, bootstrap = TRUE, iterations = 1000)
cricketm1.ci

# GRASSHOPPERS #
MEFEm1 = lmer(resp_rate~Treatment + (1|Individual), data = MEFEdf)
plot(MEFEm1)
summary(MEFEm1)
r.squaredGLMM(MEFEm1)

MEFEm1.ci <- model_parameters(MEFEm1, ci = 0.95, bootstrap = TRUE, iterations = 1000)
MEFEm1.ci

# PHIDDIPUS #
PHIDm1 = lmer(resp_rate~Treatment + (1|Individual), data = PHIDdf)
plot(PHIDm1)
summary(PHIDm1)
r.squaredGLMM(PHIDm1)

PHIDm1.ci <- model_parameters(PHIDm1, ci = 0.95, bootstrap = TRUE, iterations = 1000)
PHIDm1.ci

# LYNX # 
lynxm1 = lmer(resp_rate~Treatment + (1|Individual), data = lynxdf)
plot(lynxm1)
summary(lynxm1)
r.squaredGLMM(lynxm1)

lynxm1.ci <- model_parameters(lynxm1, ci = 0.95, bootstrap = TRUE, iterations = 1000)
lynxm1.ci

# ONAS # 
ONASm1 = lmer(resp_rate~Treatment + (1|Individual), data = ONASdf)
plot(ONASm1)
summary(ONASm1)
r.squaredGLMM(ONASm1)

ONASm1.ci <- model_parameters(ONASm1, ci = 0.95, bootstrap = TRUE, iterations = 1000)
ONASm1.ci

# TRRA # 
TRRAm1 = lmer(resp_rate~Treatment + (1|Individual), data=TRRAdf)
plot(TRRAm1)
summary(TRRAm1)
r.squaredGLMM(TRRAm1)

TRRAm1.ci <- model_parameters(TRRAm1, ci = 0.95, bootstrap = TRUE, iterations = 1000)
TRRAm1.ci

### Cleaning II: Does removing low r2 and negative respiration change results? ---- 
  # From plots above, look at PHID, ONAS, lynx

PHIDdf2 = fdata %>% filter(Species == "Phiddipus") %>% 
  filter(!resp_rate < 0 & !r2 < .25)
PHIDm2 = lmer(resp_rate~Treatment + (1|Individual), data = PHIDdf2)
plot(PHIDm2)
summary(PHIDm2)
PHIDm2.ci <- model_parameters(PHIDm2, ci = 0.95, bootstrap = TRUE, iterations = 1000)
r.squaredGLMM(PHIDm2)
# Comparision - no real differences. Remove these points.
PHIDm2.ci
PHIDm1.ci

ONASdf2 = fdata %>% filter(Species == "ONAS") %>%
  filter(!resp_rate < 0 & !r2 <.25)
ONASm2 = lmer(resp_rate~Treatment + (1|Individual), data = ONASdf2)
plot(ONASm2)
summary(ONASm2)
ONASm2.ci <- model_parameters(ONASm2, ci = 0.95, bootstrap = TRUE, iterations = 1000)
r.squaredGLMM(ONASm2)
# Comparison - no real differences. Remove these points
ONASm2.ci
ONASm1.ci

lynxdf2 = fdata %>% filter(Species == "Lynx") %>%
  filter(!resp_rate < 0 & !r2 <.25)
lynxm2 = lmer(resp_rate~Treatment + (1|Individual), data = lynxdf2)
plot(lynxm2)
summary(lynxm2)
lynxm2.ci <- model_parameters(lynxm2, ci = 0.95, bootstrap = TRUE, iterations = 1000)
r.squaredGLMM(lynxm2)
# Comparison - no real differences. Remove these points
lynxm2.ci
lynxm1.ci

#### Figures ####

detach("package:parameters", unload=TRUE)
detach("package:MuMIn", unload=TRUE)
detach("package:tidyverse", unload=TRUE)
detach("package:lubridate", unload=TRUE)

require("plyr")
require("lattice")
require("ggplot2")
require("dplyr")
require("Rmisc")
require("devtools")
require("gghalves")


# MEFE & Crickets

fdata %>% filter(!resp_rate < 0 & !r2 <.25) %>% 
  mutate(Species = replace(Species, Species == "Cricket", "GRPE")) %>%
  filter(Species %in% c("GRPE", "MEFE")) %>%
  group_by(Species, Individual, Treatment) %>%
  ggplot(aes(x=Treatment, y=resp_rate)) +
  geom_line(data = fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
              mutate(Species = replace(Species, Species == "Cricket", "GRPE")) %>%
              filter(Species %in% c("GRPE", "MEFE")) %>%
              group_by(Species,
                       Treatment, 
                       Individual) %>%
                    summarise(AvRR = mean(resp_rate)),
             aes(x = Treatment, y = AvRR, group=Individual), 
                 color="light gray", 
                 size = 0.5,
                 position = position_dodge(width=0.2)) + 
  geom_point(data = fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
              mutate(Species = replace(Species, Species == "Cricket", "GRPE")) %>%
              filter(Species %in% c("GRPE", "MEFE")) %>%
              group_by(Species,
                       Treatment, 
                       Individual) %>%
              summarise(AvRR = mean(resp_rate)) %>%
               filter(Treatment == "Olfactory"),
            aes(y = AvRR, group = Individual), 
            color="dodgerblue", 
            position = position_dodge(width=0.2)) + 
  geom_point(data = fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
               mutate(Species = replace(Species, Species == "Cricket", "GRPE")) %>%
               filter(Species %in% c("GRPE", "MEFE")) %>%
               group_by(Species,
                        Treatment, 
                        Individual) %>%
               summarise(AvRR = mean(resp_rate)) %>%
               filter(Treatment == "Baseline"),
             aes(y = AvRR, group = Individual), 
             color="darkgreen", 
             position = position_dodge(width=0.2)) +
  geom_point(data = fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
               mutate(Species = replace(Species, Species == "Cricket", "GRPE")) %>%
               filter(Species %in% c("GRPE", "MEFE")) %>%
               group_by(Species,
                        Treatment, 
                        Individual) %>%
               summarise(AvRR = mean(resp_rate)) %>%
               filter(Treatment == "Visual+Olfactory"),
             aes(y = AvRR, group = Individual), 
             color= "darkorange", 
             position = position_dodge(width=0.2)) +
  
  # GRPE VISUAL + OLFACTORY
  geom_half_boxplot(data = fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
                      mutate(Species = replace(Species, Species == "Cricket", "GRPE")) %>%
                      filter(Treatment == "Visual+Olfactory") %>%
                      filter(Species == "GRPE") %>%
                      group_by(Species,
                               Individual,
                               Treatment) %>%
                      summarise(AvRR = mean(resp_rate)),
                   aes(x = Treatment, y = AvRR),
                   position = position_nudge(x = -2.25), 
                   side = "l", outlier.shape = NA, 
                   center = TRUE, errorbar.draw = FALSE, width = .2,
                   fill = 'darkorange', alpha = 0.6) +
  
  # GRPE OLFACTORY
  geom_half_boxplot(data = fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
                      mutate(Species = replace(Species, Species == "Cricket", "GRPE")) %>%
                      filter(Treatment == "Olfactory") %>%
                      filter(Species == "GRPE") %>%
                      group_by(Species,
                               Individual,
                               Treatment) %>%
                      summarise(AvRR = mean(resp_rate)),
                    aes(x = Treatment, y = AvRR),
                    position = position_nudge(x = -1.4), 
                    side = "l", outlier.shape = NA, 
                    center = TRUE, errorbar.draw = FALSE, width = .2,
                    fill = 'dodgerblue', alpha = 0.6) +
  
  # GRPE BASELINE
  geom_half_boxplot(data = fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
                      mutate(Species = replace(Species, Species == "Cricket", "GRPE")) %>%
                      filter(Treatment == "Baseline") %>%
                      filter(Species == "GRPE") %>%
                      group_by(Species,
                               Individual,
                               Treatment) %>%
                      summarise(AvRR = mean(resp_rate)),
                    aes(x = Treatment, y = AvRR),
                    position = position_nudge(x = -.55), 
                    side = "l", outlier.shape = NA, 
                    center = TRUE, errorbar.draw = FALSE, width = .2,
                    fill = 'darkgreen', alpha = 0.6) +
  
  # MEFE VISUAL + OLFACTORY
  geom_half_boxplot(data = fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
                      filter(Treatment == "Visual+Olfactory") %>%
                      filter(Species == "MEFE") %>%
                      group_by(Species,
                               Individual,
                               Treatment) %>%
                      summarise(AvRR = mean(resp_rate)),
                    aes(x = Treatment, y = AvRR),
                    position = position_nudge(x = -2.25), 
                    side = "l", outlier.shape = NA, 
                    center = TRUE, errorbar.draw = FALSE, width = .2,
                    fill = 'darkorange', alpha = 0.6) +
  
  # MEFE OLFACTORY
  geom_half_boxplot(data = fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
                      filter(Treatment == "Olfactory") %>%
                      filter(Species == "MEFE") %>%
                      group_by(Species,
                               Individual,
                               Treatment) %>%
                      summarise(AvRR = mean(resp_rate)),
                    aes(x = Treatment, y = AvRR),
                    position = position_nudge(x = -1.4), 
                    side = "l", outlier.shape = NA, 
                    center = TRUE, errorbar.draw = FALSE, width = .2,
                    fill = 'dodgerblue', alpha = 0.6) +
  
  # MEFE BASELINE
  geom_half_boxplot(data = fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
                      mutate(Species = replace(Species, Species == "Cricket", "GRPE")) %>%
                      filter(Treatment == "Baseline") %>%
                      filter(Species == "MEFE") %>%
                      group_by(Species,
                               Individual,
                               Treatment) %>%
                      summarise(AvRR = mean(resp_rate)),
                    aes(x = Treatment, y = AvRR),
                    position = position_nudge(x = -.55), 
                    side = "l", outlier.shape = NA, 
                    center = TRUE, errorbar.draw = FALSE, width = .2,
                    fill = 'darkgreen', alpha = 0.6) +
  
  # OVERLAY MEANS
  geom_line(data = fdata %>% filter(!resp_rate < 0 & !r2 <.25) %>%
              mutate(Species = replace(Species, Species == "Cricket", "GRPE")) %>%
              filter(Species %in% c("GRPE", "MEFE")) %>%
              group_by(Species, 
                       Treatment) %>%
              summarise(AvRR = mean(resp_rate)), 
            aes(x=Treatment, y = AvRR, group = Species), color = "black", size = .8) +
  
  geom_point(data = fdata %>% filter(!resp_rate < 0 & !r2 <.25) %>%
               mutate(Species = replace(Species, Species == "Cricket", "GRPE")) %>%
               filter(Species %in% c("GRPE", "MEFE")) %>%
               group_by(Species, 
                        Treatment) %>%
               summarise(AvRR = mean(resp_rate)), 
             aes(x=Treatment, y = AvRR, group = Species), color = "black", size = 2) +
  
  theme_classic() +
  coord_cartesian(xlim = c(3.5,.75)) +
  labs(y = "Respiration rate (uL CO2 g-1 min-1)", x = "Predator Cue") +
  facet_grid(.~Species)

# ONAS & TRRA
fdata %>% filter(!resp_rate < 0 & !r2 <.25) %>% 
  filter(Species %in% c("ONAS", "TRRA")) %>% 
  group_by(Species, Individual, Treatment) %>%
  ggplot(aes(x=Treatment, y=resp_ratel)) +
  geom_line(data = fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
              filter(Species %in% c("ONAS", "TRRA")) %>%
              group_by(Species,
                       Treatment, 
                       Individual) %>%
              summarise(AvRR = mean(resp_rate)),
            aes(x = Treatment, y = AvRR, group=Individual), 
            color="light gray", 
            size = 0.5,
            position = position_dodge(width=0.2)) + 
  
  geom_point(data = fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
               filter(Species %in% c("ONAS", "TRRA")) %>%
               group_by(Species,
                        Treatment, 
                        Individual) %>%
               summarise(AvRR = mean(resp_rate)) %>%
               filter(Treatment == "Olfactory"),
             aes(y = AvRR, group = Individual), 
             color="dodgerblue", 
             position = position_dodge(width=0.2)) + 
  geom_point(data = fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
               filter(Species %in% c("ONAS", "TRRA")) %>%
               group_by(Species,
                        Treatment, 
                        Individual) %>%
               summarise(AvRR = mean(resp_rate)) %>%
               filter(Treatment == "Baseline"),
             aes(y = AvRR, group = Individual), 
             color="darkgreen", 
             position = position_dodge(width=0.2)) +
  geom_point(data = fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
               filter(Species %in% c("ONAS", "TRRA")) %>%
               group_by(Species,
                        Treatment, 
                        Individual) %>%
               summarise(AvRR = mean(resp_rate)) %>%
               filter(Treatment == "Visual+Olfactory"),
             aes(y = AvRR, group = Individual), 
             color= "darkorange", 
             position = position_dodge(width=0.2)) +
  
  # ONAS VISUAL + OLFACTORY
  geom_half_boxplot(data = fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
                      filter(Treatment == "Visual+Olfactory") %>%
                      filter(Species == "ONAS") %>%
                      group_by(Species,
                               Individual,
                               Treatment) %>%
                      summarise(AvRR = mean(resp_rate)),
                    aes(x = Treatment, y = AvRR),
                    position = position_nudge(x = -2.25), 
                    side = "l", outlier.shape = NA, 
                    center = TRUE, errorbar.draw = FALSE, width = .2,
                    fill = 'darkorange', alpha = 0.6) +
  
  # ONAS OLFACTORY
 
  geom_half_boxplot(data = fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
                      filter(Treatment == "Olfactory") %>%
                      filter(Species == "ONAS") %>%
                      group_by(Species,
                               Individual,
                               Treatment) %>%
                      summarise(AvRR = mean(resp_rate)),
                    aes(x = Treatment, y = AvRR),
                    position = position_nudge(x = -1.4), 
                    side = "l", outlier.shape = NA, 
                    center = TRUE, errorbar.draw = FALSE, width = .2,
                    fill = 'dodgerblue', alpha = 0.6) +
  
  # ONAS BASELINE

  geom_half_boxplot(data = fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
                      filter(Treatment == "Baseline") %>%
                      filter(Species == "ONAS") %>%
                      group_by(Species,
                               Individual,
                               Treatment) %>%
                      summarise(AvRR = mean(resp_rate)),
                    aes(x = Treatment, y = AvRR),
                    position = position_nudge(x = -.55), 
                    side = "l", outlier.shape = NA, 
                    center = TRUE, errorbar.draw = FALSE, width = .2,
                    fill = 'darkgreen', alpha = 0.6) +
  
  # TRRA VISUAL + OLFACTORY

  geom_half_boxplot(data = fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
                      filter(Treatment == "Visual+Olfactory") %>%
                      filter(Species == "TRRA") %>%
                      group_by(Species,
                               Individual,
                               Treatment) %>%
                      summarise(AvRR = mean(resp_rate)),
                    aes(x = Treatment, y = AvRR),
                    position = position_nudge(x = -2.25), 
                    side = "l", outlier.shape = NA, 
                    center = TRUE, errorbar.draw = FALSE, width = .2,
                    fill = 'darkorange', alpha = 0.6) +
  
  # TRRA OLFACTORY

  geom_half_boxplot(data = fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
                      filter(Treatment == "Olfactory") %>%
                      filter(Species == "TRRA") %>%
                      group_by(Species,
                               Individual,
                               Treatment) %>%
                      summarise(AvRR = mean(resp_rate)),
                    aes(x = Treatment, y = AvRR),
                    position = position_nudge(x = -1.4), 
                    side = "l", outlier.shape = NA, 
                    center = TRUE, errorbar.draw = FALSE, width = .2,
                    fill = 'dodgerblue', alpha = 0.6) +
  
  # TRRA BASELINE

  geom_half_boxplot(data = fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
                      mutate(Species = replace(Species, Species == "Cricket", "GRPE")) %>%
                      filter(Treatment == "Baseline") %>%
                      filter(Species == "TRRA") %>%
                      group_by(Species,
                               Individual,
                               Treatment) %>%
                      summarise(AvRR = mean(resp_rate)),
                    aes(x = Treatment, y = AvRR),
                    position = position_nudge(x = -.55), 
                    side = "l", outlier.shape = NA, 
                    center = TRUE, errorbar.draw = FALSE, width = .2,
                    fill = 'darkgreen', alpha = 0.6) +
  
  # OVERLAY MEANS
  geom_point(data = fdata %>% filter(!resp_rate < 0 & !r2 <.25) %>%
               filter(Species %in% c("ONAS", "TRRA")) %>%
               group_by(Species, 
                        Treatment) %>%
               summarise(AvRR = mean(resp_rate)), 
             aes(x=Treatment, y = AvRR, group = Species), color = "black", size = 2) +
  geom_line(data = fdata %>% filter(!resp_rate < 0 & !r2 <.25) %>%
              filter(Species %in% c("ONAS", "TRRA")) %>%
              group_by(Species, 
                       Treatment) %>%
              summarise(AvRR = mean(resp_rate)), 
            aes(x=Treatment, y = AvRR, group = Species), color = "black", size = .8) +
  
  theme_classic() +
  coord_cartesian(xlim = c(3.5,.75)) +
  labs(y = "Respiration rate (uL CO2 g-1 min-1)", x = "Predator Cue") +
  facet_grid(.~Species)

# LYNX & PHID
fdata %>% filter(!resp_rate < -.5 & !r2 <.25) %>% 
  mutate(Species = replace(Species, Species == "Lynx", "OXSA")) %>%
  mutate(Species = replace(Species, Species == "Phiddipus", "PHspp")) %>%
  filter(Species %in% c("OXSA", "PHspp")) %>%
  group_by(Species, 
           Individual, 
           Treatment) %>%
  ggplot(aes(x=Treatment, y=resp_rate)) +
  
  geom_line(data = fdata %>% filter(!resp_rate < -.5 & !r2 <.25) %>% 
              mutate(Species = replace(Species, Species == "Lynx", "OXSA")) %>%
              mutate(Species = replace(Species, Species == "Phiddipus", "PHspp")) %>%
              filter(Species %in% c("OXSA", "PHspp")) %>%
              group_by(Species,
                       Treatment, 
                       Individual) %>%
              summarise(AvRR = mean(resp_rate)),
            aes(x = Treatment, y = AvRR, group=Individual), 
            color="light gray", 
            size = 0.5,
            position = position_dodge(width=0.2)) + 
  
  geom_point(data = fdata %>% filter(!resp_rate < -.5 & !r2 <.25) %>% 
               mutate(Species = replace(Species, Species == "Lynx", "OXSA")) %>%
               mutate(Species = replace(Species, Species == "Phiddipus", "PHspp")) %>%
               filter(Species %in% c("OXSA", "PHspp")) %>%
               group_by(Species,
                        Treatment, 
                        Individual) %>%
               summarise(AvRR = mean(resp_rate)) %>%
               filter(Treatment == "Olfactory"),
             aes(y = AvRR, group = Individual), 
             color="dodgerblue", 
             position = position_dodge(width=0.2)) + 
  geom_point(data = fdata %>% filter(!resp_rate < -.5 & !r2 <.25) %>% 
               mutate(Species = replace(Species, Species == "Lynx", "OXSA")) %>%
               mutate(Species = replace(Species, Species == "Phiddipus", "PHspp")) %>%
               filter(Species %in% c("OXSA", "PHspp")) %>%
               group_by(Species,
                        Treatment, 
                        Individual) %>%
               summarise(AvRR = mean(resp_rate)) %>%
               filter(Treatment == "Baseline"),
             aes(y = AvRR, group = Individual), 
             color="darkgreen", 
             position = position_dodge(width=0.2)) +
  geom_point(data = fdata %>% filter(!resp_rate < -.5 & !r2 <.25) %>% 
               mutate(Species = replace(Species, Species == "Lynx", "OXSA")) %>%
               mutate(Species = replace(Species, Species == "Phiddipus", "PHspp")) %>%
               filter(Species %in% c("OXSA", "PHspp")) %>%
               group_by(Species,
                        Treatment, 
                        Individual) %>%
               summarise(AvRR = mean(resp_rate)) %>%
               filter(Treatment == "Visual+Olfactory"),
             aes(y = AvRR, group = Individual), 
             color= "darkorange", 
             position = position_dodge(width=0.2)) +
  
  # OXSA VISUAL + OLFACTORY

  geom_half_boxplot(data = fdata %>% filter(!resp_rate < -1 & !r2 <.25) %>% 
                      mutate(Species = replace(Species, Species == "Lynx", "OXSA")) %>%
                      filter(Treatment == "Visual+Olfactory") %>%
                      filter(Species == "OXSA") %>%
                      group_by(Species,
                               Individual,
                               Treatment) %>%
                      summarise(AvRR = mean(resp_rate)),
                    aes(x = Treatment, y = AvRR),
                    position = position_nudge(x = -2.25), 
                    side = "l", outlier.shape = NA, 
                    center = TRUE, errorbar.draw = FALSE, width = .2,
                    fill = 'darkorange', alpha = 0.6) +
  
  # OXSA OLFACTORY

  geom_half_boxplot(data = fdata %>% filter(!resp_rate < -.5 & !r2 <.25) %>% 
                      mutate(Species = replace(Species, Species == "Lynx", "OXSA")) %>%
                      filter(Treatment == "Olfactory") %>%
                      filter(Species == "OXSA") %>%
                      group_by(Species,
                               Individual,
                               Treatment) %>%
                      summarise(AvRR = mean(resp_rate)),
                    aes(x = Treatment, y = AvRR),
                    position = position_nudge(x = -1.4), 
                    side = "l", outlier.shape = NA, 
                    center = TRUE, errorbar.draw = FALSE, width = .2,
                    fill = 'dodgerblue', alpha = 0.6) +
  
  # OXSA BASELINE

  geom_half_boxplot(data = fdata %>% filter(!resp_rate < -.5 & !r2 <.25) %>% 
                      mutate(Species = replace(Species, Species == "Lynx", "OXSA")) %>%
                      filter(Treatment == "Olfactory") %>%
                      filter(Species == "OXSA") %>%
                      group_by(Species,
                               Individual,
                               Treatment) %>%
                      summarise(AvRR = mean(resp_rate)),
                    aes(x = Treatment, y = AvRR),
                    position = position_nudge(x = -1.55), 
                    side = "l", outlier.shape = NA, 
                    center = TRUE, errorbar.draw = FALSE, width = .2,
                    fill = 'darkgreen', alpha = 0.6) +
  
  # PHsp VISUAL + OLFACTORY

  geom_half_boxplot(data = fdata %>% filter(!resp_rate < 0 & !r2 <.25) %>% 
                      mutate(Species = replace(Species, Species == "Phiddipus", "PHspp")) %>%
                      filter(Treatment == "Visual+Olfactory") %>%
                      filter(Species == "PHspp") %>%
                      group_by(Species,
                               Individual,
                               Treatment) %>%
                      summarise(AvRR = mean(resp_rate)),
                    aes(x = Treatment, y = AvRR),
                    position = position_nudge(x = -2.25), 
                    side = "l", outlier.shape = NA, 
                    center = TRUE, errorbar.draw = FALSE, width = .2,
                    fill = 'darkorange', alpha = 0.6) +
  
  # PHsp OLFACTORY

  geom_half_boxplot(data = fdata %>% filter(!resp_rate < 0 & !r2 <.25) %>% 
                      mutate(Species = replace(Species, Species == "Phiddipus", "PHspp")) %>%
                      filter(Treatment == "Olfactory") %>%
                      filter(Species == "PHspp") %>%
                      group_by(Species,
                               Individual,
                               Treatment) %>%
                      summarise(AvRR = mean(resp_rate)),
                    aes(x = Treatment, y = AvRR),
                    position = position_nudge(x = -1.4), 
                    side = "l", outlier.shape = NA, 
                    center = TRUE, errorbar.draw = FALSE, width = .2,
                    fill = 'dodgerblue', alpha = 0.6) +
  
  # PHsp BASELINE

  geom_half_boxplot(data = fdata %>% filter(!resp_rate < 0 & !r2 <.25) %>% 
                      mutate(Species = replace(Species, Species == "Phiddipus", "PHspp")) %>%
                      filter(Treatment == "Baseline") %>%
                      filter(Species == "PHspp") %>%
                      group_by(Species,
                               Individual,
                               Treatment) %>%
                      summarise(AvRR = mean(resp_rate)),
                    aes(x = Treatment, y = AvRR),
                    position = position_nudge(x = -.55), 
                    side = "l", outlier.shape = NA, 
                    center = TRUE, errorbar.draw = FALSE, width = .2,
                    fill = 'darkgreen', alpha = 0.6) +
  
  # OVERLAY MEANS
  geom_point(data = fdata %>% filter(!resp_rate < 0 & !r2 <.25) %>%
               mutate(Species = replace(Species, Species == "Lynx", "OXSA")) %>%
               mutate(Species = replace(Species, Species == "Phiddipus", "PHspp")) %>%
               filter(Species %in% c("OXSA", "PHspp")) %>%
               group_by(Species, 
                        Treatment) %>%
               summarise(AvRR = mean(resp_rate)), 
             aes(x=Treatment, y = AvRR, group = Species), color = "black", size = 2) +
  geom_line(data = fdata %>% filter(!resp_rate < 0 & !r2 <.25) %>%
              mutate(Species = replace(Species, Species == "Lynx", "OXSA")) %>%
              mutate(Species = replace(Species, Species == "Phiddipus", "PHspp")) %>%
              filter(Species %in% c("OXSA", "PHspp")) %>%
              group_by(Species, 
                       Treatment) %>%
              summarise(AvRR = mean(resp_rate)), 
            aes(x=Treatment, y = AvRR, group = Species), color = "black", size = 1.1) +
 
  theme_classic() +
  coord_cartesian(xlim = c(3.5,.75)) +
  labs(y = "Respiration rate (uL CO2 g-1 min-1)", x = "Predator Cue") +
  facet_grid(.~Species)

### **AQUATIC** ----
 ### CLEANING FUNCTIONS ----
 ### Main code written by A. Arietta.
 ### Last update 2020 Sep 1 by N.R. Sommer
 
require(marelac)
require(tidyverse)
require(readr)
require(lubridate)
require(ggalt)
require(zoo)
 
RespDataPrep <- function(ExperimentData){

  # Read in the FILENAME_raw.txt file and extract experiment data
  X <- read_tsv(ExperimentData, skip = 38, 
                 col_names = FALSE, 
                 cols_only(X1 = "?", X2 = "c", X3 = "n", X4 = "n", X6 = "n", 
                           X7 = "n", X9 = "n", X10 = "n", X12 = "n", X13 = "n", 
                           X15 = "n", X16 = "n", X18 = "n", X19 = "n", X21 = "n", 
                           X22 = "n", X24 = "n", X25 = "n", X27 = "n", X28 = "n")) 
  # Rename the columns
  colnames(X) <- c("DateTime", "Phase", "Baro", "Sal", "CH1.temp", "CH1.O2persat", 
                   "CH2.temp", "CH2.O2persat", "CH3.temp", "CH3.O2persat", "CH4.temp", 
                   "CH4.O2persat", "CH5.temp", "CH5.O2persat", "CH6.temp", "CH6.O2persat", 
                   "CH7.temp", "CH7.O2persat", "CH8.temp", "CH8.O2persat") 
  
  # Function to convert  original datetime column fromtext file into lubridate datetime 
   AutoResp.convert.DateTime <- function(i){
     Xdates <- (strsplit(as.character(i), "/"))
     Xdates <- matrix(unlist(Xdates), ncol=4, byrow=TRUE)
     Xdates <- as.data.frame(Xdates)
     Date <- paste(Xdates[,3], Xdates[,1], Xdates[,2], sep = "-")
     Time <- as.character(Xdates$V4)
     datetime <- paste(Date, Time, sep = "T")
     DT <- as_datetime(datetime)
     DT} 
   # This function reformats the temp and O2 columns from the wide-format raw text file so that they can be combined into a tidy-format dataframe in the next step.
   # It also converts the DateTime into a variable "sec" which is the second from midnight
   # It also computes a rolling average of raw O2 values to reduce noise using a window value of 5 seconds and aligning the new data to the right retains only averages of measurement data to the end of the measurement period.
   AutoResp.combine.channels <- function(data, colTemp, colO2, ch){
     data %>% 
       select(DateTime:Sal, "Temp" = colTemp, "O.2persat_raw" = colO2) %>% 
       mutate(CH = ch) %>% 
       mutate(DateTime = AutoResp.convert.DateTime(DateTime)) %>% 
       mutate(sec = ((hour(DateTime)*3600) + 
                       (minute(DateTime)*60) + 
                       (second(DateTime)))) %>% 
       mutate(O.2persat = rollmean(O.2persat_raw, k = 5, align = "right", fill = NA))
   } 

   # This combines all of the data from each chamer into a tidy-format table.
   Xfull <- AutoResp.combine.channels(X, "CH1.temp", "CH1.O2persat", 1) %>%
     bind_rows(AutoResp.combine.channels(X, "CH2.temp", "CH2.O2persat", 2)) %>%
     bind_rows(AutoResp.combine.channels(X, "CH3.temp", "CH3.O2persat", 3)) %>%
     bind_rows(AutoResp.combine.channels(X, "CH4.temp", "CH4.O2persat", 4)) %>%
     bind_rows(AutoResp.combine.channels(X, "CH5.temp", "CH5.O2persat", 5)) %>%
     bind_rows(AutoResp.combine.channels(X, "CH6.temp", "CH6.O2persat", 6)) %>%
     bind_rows(AutoResp.combine.channels(X, "CH7.temp", "CH7.O2persat", 7)) %>%
     bind_rows(AutoResp.combine.channels(X, "CH8.temp", "CH8.O2persat", 8)) 
   
   Xfull$Resp.DOY <- yday(Xfull$DateTime) # Creates a column of DOY date in order to combine with larval data
   
   Xfull$O2.mg_per_L <- marelac::gas_satconc(S = Xfull$Sal, 
                                             t = Xfull$Temp, 
                                             P = Xfull$Baro/1000, 
                                             species = 'O2') * 1e-6 * marelac::molweight('O2') * 1e3 * Xfull$O.2persat/100 
   # Converts O2 readings from percent saturation to mg/L using solubility conversions from the 'marelac' package
   
   Xinfo <- read_tsv(ExperimentData, 
                     skip = 25, 
                     n_max = 8, 
                     col_names = FALSE, 
                     cols_only(X1 = "c", X2 = "n", X3 = "n", X4 = "n")) 
   # Reads in the pertinent portion of the info data output form AutoResp
   
   colnames(Xinfo) <- c("CH", "CHvol.mL", "Tubevol.mL", "Mass.g") # Rename the columns
   Xinfo$CH <- seq(1:length(Xinfo$CH)) # Adds a numerical variable for the channel numer
   
   X <- left_join(Xfull, Xinfo, by = "CH") # Joins the experiment data and info
   
   X <- filter(X, CH != 5) # Remove chamber 5 that was not used in the experiment
 } # Preps the raw data

RespDataSummary <- function(X, Type){
   tab <- c()
   for(i in unique(X$CH)){
     for(j in unique(filter(X, grepl("M", Phase))$Phase)){
       dat <- filter(X, CH == i, Phase == j) %>%
         filter((sec - sec[1]) > 200) # Excludes the first 200 seconds of the measurement.
       lm1 <- lm(O2.mg_per_L ~ sec, data = dat)
       sum1 <- summary(lm1)
       lo1 <- loess(O2.mg_per_L ~ sec, data = dat, span = 0.5)
       sum2 <- summary(lo1)
       slope <- sum1$coefficients[2,1]
       R2 <- sum1$r.squared
       df <- data.frame("CH" = i, "Phase" = j, "Slope" = slope, "R2" = R2, 
                        "RMSE.lm" = sqrt(mean((sum1$residuals)^2)), 
                        "RMSE.lo" = sqrt(mean((sum2$residuals)^2)), 
                        "Mass.g" = mean(dat$Mass.g), 
                        "Temp.C" = mean(dat$Temp), 
                        sec = mean(dat$sec), 
                        CHvol.mL = mean(dat$CHvol.mL), 
                        Tubevol.mL = mean(dat$Tubevol.mL), 
                        Resp.DOY = mean(dat$Resp.DOY))
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
   RespSum %>% group_by(CH) %>% summarise(Mass = mean(Mass.g), 
                                          Baro = mean(Baro), 
                                          Temp = mean(Temp), 
                                          Resp.DOY = mean(Resp.DOY), 
                                          Sal = mean(Sal), 
                                          CHvol = mean(CHvol.mL), 
                                          Tubevol = mean(Tubevol.mL))} # Summary of input variables to check for errors.
RespSumPlot1 <- function(RespSum){
   RespSum %>%
     ggplot(aes(x = Slope, y = RMSE.lm - RMSE.lo, col = Phase)) +
     geom_point(size = 4) +
     geom_encircle() +
     theme_minimal()
 }
 
 
# CLEANING OF AQUATIC DATA BY DATE ----
  ### 11 Jul ----
# Define the run of interest:
Pre <- RespDataPrep("Data/20190711_PredatorExp_blankPre_raw.txt")
Acc <- RespDataPrep("Data/20190711_PredatorExp_acc_raw.txt")
SMR <- RespDataPrep("Data/20190711_PredatorExp_smr_None_raw.txt")
Treat <- RespDataPrep("Data/20190711_PredatorExp_smr_Treat_raw.txt")
Post <- RespDataPrep("Data/20190711_PredatorExp_blankPost_raw.txt")

RespPrepPlot1(Pre) # Quality control for raw data
RespPrepPlot2(Pre) # Quality control for all slopes
Resp.InputQC(Pre)

Pre <- filter(Pre, Phase != "M4") # Removing extra, short measurment phase

RespPrepPlot1(Post) # Quality control for raw data
RespPrepPlot2(Post) # Quality control for all slopes
Resp.InputQC(Post)

Post <- filter(Post, Phase != "M4") # Removing extra, short measurment phase

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

SumSMR <- anti_join(SumSMR, data.frame("CH" = 8, "Phase" = "M5")) # Remove bad measurements

SumTreat <- RespDataSummary(Treat, "Treat")

SumTreat <- anti_join(SumTreat, data.frame("CH" = 8, "Phase" = "M8")) # Remove bad measurements

SumPost <- RespDataSummary(Post, "Post")

RespSumPlot1(SumPre) # Quality control for slope values
RespSumPlot1(SumAcc) # Quality control for slope values
RespSumPlot1(SumSMR) # Quality control for slope values
RespSumPlot1(SumTreat) # Quality control for slope values
RespSumPlot1(SumPost) # Quality control for slope values

# Define background respiration ===
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

# Define experimental respiration rates ===

# . No Treatment analysis: ===

SumSMR <- left_join(SumSMR, BG.corrrection, by = "CH") %>% mutate(corO2.mg_h = O2.mg_h - ((sec - PreSec) * MO2.Cor), corO2.mg_h.g = corO2.mg_h / Mass.g)

# Check for influence of temperature over time:
SumSMR %>% filter(Phase != "M1") %>%
  ggplot(aes(x = Temp.C, corO2.mg_h, color = as.factor(CH))) +
  geom_point()+
  geom_smooth(method = 'loess', se = FALSE)
summary(lm(corO2.mg_h ~ Temp.C * as.factor(CH), data = SumSMR))

# Check for low quality measurements under R2 = 0.95
SumSMR %>% filter(Phase != "M1") %>%
  ggplot(aes(x = CH, y = corO2.mg_h.g)) +
  geom_jitter(width = 0.1, aes(color = ifelse(RMSE.lm - RMSE.lo >= 0.01, 2, 3))) +
  scale_color_identity()

FinalSMR <- SumSMR  %>% select(Resp.DOY, Type, CH, Phase, corO2.mg_h, corO2.mg_h.g, R2, Mass.g, Temp.C, sec, CHvol.mL, Tubevol.mL)

# . Treatment SMR analysis: ===

SumTreat <- left_join(SumTreat, BG.corrrection, by = "CH") %>% mutate(corO2.mg_h = O2.mg_h - ((sec - PreSec) * MO2.Cor), corO2.mg_h.g = corO2.mg_h / Mass.g)

# Check for influence of temperature over time:
SumTreat %>% filter(Phase != "M1") %>%
  ggplot(aes(x = Temp.C, corO2.mg_h, color = as.factor(CH))) +
  geom_point()+
  geom_smooth(method = 'loess', se = FALSE)
summary(lm(corO2.mg_h ~ Temp.C * CH, data = SumTreat))

# Check for low quality measurements under R2 = 0.95
SumTreat %>% filter(Phase != "M1") %>%
  ggplot(aes(x = CH, y = corO2.mg_h.g)) +
  geom_jitter(width = 0.2, aes(color = ifelse(RMSE.lm - RMSE.lo >= 0.01, 2, 3))) +
  scale_color_identity()

FinalTreat <- SumTreat  %>% select(Resp.DOY, Type, CH, Phase, corO2.mg_h, corO2.mg_h.g, R2, Mass.g, Temp.C, sec, CHvol.mL, Tubevol.mL)

# . Combined analysis: ==

Final <- bind_rows(SumSMR, SumTreat)

ggplot(Final, aes(x = CH, y = corO2.mg_h.g, col = Type)) +
  geom_jitter(width = 0.2) +
  theme_minimal()

summary(lm(corO2.mg_h.g ~ Type * as.factor(CH), data = Final))

write.csv(Final, "Jul11.csv")

## REPEAT ##
 
### 12 Jul ----
# Define the run of interest:
Pre <- RespDataPrep("Data/20190712_PredatorExp_blankPre_raw.txt")
Acc <- RespDataPrep("Data/20190712_PredatorExp_acc_raw.txt")
SMR <- RespDataPrep("Data/20190712_PredatorExp_smr_None_raw.txt")
Treat <- RespDataPrep("Data/20190712_PredatorExp_smr_Treat_raw.txt")
Post <- RespDataPrep("Data/20190712_PredatorExp_blankPost_raw.txt")

RespPrepPlot1(Pre) # Quality control for raw data
RespPrepPlot2(Pre) # Quality control for all slopes
Resp.InputQC(Pre)

Pre <- filter(Pre, Phase != "M1") # Removing wonky first measurement phase

RespPrepPlot1(Post) # Quality control for raw data
RespPrepPlot2(Post) # Quality control for all slopes
Resp.InputQC(Post)

Post <- filter(Post, Phase != "M1") # Removing wonky first measurement phase

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

SumTreat <- anti_join(SumTreat, data.frame("CH" = rep(6:8, each = 3), "Phase" = rep(c("M1", "M2", "M3"), 3))) # Remove Phases 1,2,3 for CH 6,7,8 during the treatment runs.

SumPost <- RespDataSummary(Post, "Post")

RespSumPlot1(SumPre) # Quality control for slope values
RespSumPlot1(SumAcc) # Quality control for slope values
RespSumPlot1(SumSMR) # Quality control for slope values
RespSumPlot1(SumTreat) # Quality control for slope values
RespSumPlot1(SumPost) # Quality control for slope values

# Define background respiration ===

SumPre <- SumPre %>% filter(Phase != "M1") %>% group_by(Type, CH) %>% summarise(avMO2 = mean(O2.mg_h), minMO2 = min(O2.mg_h), maxMO2 = max(O2.mg_h), avR2 = mean(R2), minR2= min(R2), maxR2 = max(R2), Temp.C = mean(Temp.C), sec = mean(sec))

SumPost <- SumPost %>% filter(Phase != "M1") %>% group_by(Type, CH) %>% summarise(avMO2 = mean(O2.mg_h), minMO2 = min(O2.mg_h), maxMO2 = max(O2.mg_h), avR2 = mean(R2), minR2= min(R2), maxR2 = max(R2), Temp.C = mean(Temp.C), sec = mean(sec))

(SumPost$avMO2 - SumPre$avMO2) / (SumPost$sec - SumPre$sec)

BG.corrrection <- left_join(
  (SumPre %>% group_by() %>% select(CH, avMO2, sec) %>% rename(PreMO2 = avMO2, PreSec = sec)), 
  (SumPost %>% group_by() %>% select(CH, avMO2, sec) %>% rename(PostMO2 = avMO2, PostSec = sec)), by = "CH") %>%
  mutate(MO2.Cor = (PostMO2 - PreMO2)/(PostSec - PreSec)) %>%
  select(CH, PreSec, PostSec, MO2.Cor)

(BG.corrrection$PostSec - BG.corrrection$PreSec) * BG.corrrection$MO2.Cor # In order to correct for the background respiration, I need to subtract the time in sec of the phase to get the number of seconds since the beginning of the experiment and then multiply by the correction factor.

# Define experimental respiration rates ===

# . No Treatment analysis: ===

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

# . Treatment SMR analysis: ===

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

# . Combined analysis: ==

Final <- bind_rows(SumSMR, SumTreat)

ggplot(Final, aes(x = CH, y = corO2.mg_h.g, col = Type)) +
  geom_jitter(width = 0.2) +
  theme_minimal()

summary(lm(corO2.mg_h.g ~ Type * as.factor(CH), data = Final))

write.csv(Final, "Jul12.csv")
 
 ### 15 Jul ----
Pre <- RespDataPrep("Data/20190715_PredatorExp_blankPre_raw.txt")
Acc <- RespDataPrep("Data/20190715_PredatorExp_acc_raw.txt")
SMR <- RespDataPrep("Data/20190715_PredatorExp_smr_None_raw.txt")
Treat <- RespDataPrep("Data/20190715_PredatorExp_smr_Treat_raw.txt")
Post <- RespDataPrep("Data/20190715_PredatorExp_blankPost_raw.txt")

RespPrepPlot1(Pre) # Quality control for raw data
RespPrepPlot2(Pre) # Quality control for all slopes
Resp.InputQC(Pre)

Pre <- filter(Pre, Phase != "M5") # Removing extra, short measurment phase

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

# Define background respiration ==

SumPre <- SumPre %>% filter(Phase != "M1") %>% group_by(Type, CH) %>% summarise(avMO2 = mean(O2.mg_h), minMO2 = min(O2.mg_h), maxMO2 = max(O2.mg_h), avR2 = mean(R2), minR2= min(R2), maxR2 = max(R2), Temp.C = mean(Temp.C), sec = mean(sec))

SumPost <- SumPost %>% filter(Phase != "M1") %>% group_by(Type, CH) %>% summarise(avMO2 = mean(O2.mg_h), minMO2 = min(O2.mg_h), maxMO2 = max(O2.mg_h), avR2 = mean(R2), minR2= min(R2), maxR2 = max(R2), Temp.C = mean(Temp.C), sec = mean(sec))

(SumPost$avMO2 - SumPre$avMO2) / (SumPost$sec - SumPre$sec)

BG.corrrection <- left_join(
  (SumPre %>% group_by() %>% select(CH, avMO2, sec) %>% rename(PreMO2 = avMO2, PreSec = sec)), 
  (SumPost %>% group_by() %>% select(CH, avMO2, sec) %>% rename(PostMO2 = avMO2, PostSec = sec)), by = "CH") %>%
  mutate(MO2.Cor = (PostMO2 - PreMO2)/(PostSec - PreSec)) %>%
  select(CH, PreSec, PostSec, MO2.Cor)

(BG.corrrection$PostSec - BG.corrrection$PreSec) * BG.corrrection$MO2.Cor # In order to correct for the background respiration, I need to subtract the time in sec of the phase to get the number of seconds since the beginning of the experiment and then multiply by the correction factor.

# Define experimental respiration rates ===

# . No Treatment analysis: ===

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

# . Treatment SMR analysis: ===

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

# . Combined analysis: ==

Final <- bind_rows(SumSMR, SumTreat)

ggplot(Final, aes(x = CH, y = corO2.mg_h.g, col = Type)) +
  geom_jitter(width = 0.2) +
  theme_minimal()

summary(lm(corO2.mg_h.g ~ Type * as.factor(CH), data = Final))

write.csv(Final, "Jul15.csv")
 
### 16 Jul ----
Pre <- RespDataPrep("Data/20190716_PredatorExp_blankPre_raw.txt")
Acc <- RespDataPrep("Data/20190716_PredatorExp_acc_raw.txt")
SMR <- RespDataPrep("Data/20190716_PredatorExp_smr_None_raw.txt")
Treat <- RespDataPrep("Data/20190716_PredatorExp_smr_Treat_raw.txt")
Post <- RespDataPrep("Data/20190716_PredatorExp_blankPost_raw.txt")

RespPrepPlot1(Pre) # Quality control for raw data
RespPrepPlot2(Pre) # Quality control for all slopes
Resp.InputQC(Pre)

Pre <- filter(Pre, Phase != "M4") # Removing extra, short measurment phase

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

# Define background respiration ==

SumPre <- SumPre %>% filter(Phase != "M1") %>% group_by(Type, CH) %>% summarise(avMO2 = mean(O2.mg_h), minMO2 = min(O2.mg_h), maxMO2 = max(O2.mg_h), avR2 = mean(R2), minR2= min(R2), maxR2 = max(R2), Temp.C = mean(Temp.C), sec = mean(sec))

SumPost <- SumPost %>% filter(Phase != "M1") %>% group_by(Type, CH) %>% summarise(avMO2 = mean(O2.mg_h), minMO2 = min(O2.mg_h), maxMO2 = max(O2.mg_h), avR2 = mean(R2), minR2= min(R2), maxR2 = max(R2), Temp.C = mean(Temp.C), sec = mean(sec))

(SumPost$avMO2 - SumPre$avMO2) / (SumPost$sec - SumPre$sec)

BG.corrrection <- left_join(
  (SumPre %>% group_by() %>% select(CH, avMO2, sec) %>% rename(PreMO2 = avMO2, PreSec = sec)), 
  (SumPost %>% group_by() %>% select(CH, avMO2, sec) %>% rename(PostMO2 = avMO2, PostSec = sec)), by = "CH") %>%
  mutate(MO2.Cor = (PostMO2 - PreMO2)/(PostSec - PreSec)) %>%
  select(CH, PreSec, PostSec, MO2.Cor)

(BG.corrrection$PostSec - BG.corrrection$PreSec) * BG.corrrection$MO2.Cor # In order to correct for the background respiration, I need to subtract the time in sec of the phase to get the number of seconds since the beginning of the experiment and then multiply by the correction factor.

# Define experimental respiration rates ===

# . No Treatment analysis: ===

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

# . Treatment SMR analysis: ===

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

# . Combined analysis: ==

Final <- bind_rows(SumSMR, SumTreat)

ggplot(Final, aes(x = CH, y = corO2.mg_h.g, col = Type)) +
  geom_jitter(width = 0.2) +
  theme_minimal()

summary(lm(corO2.mg_h.g ~ Type * as.factor(CH), data = Final))

write.csv(Final, "Jul16.csv")
 
 
 ### 17 Jul ----
# Define the run of interest:
Pre <- RespDataPrep("Data/20190717_PredatorExp_blankPre_raw.txt")
Acc <- RespDataPrep("Data/20190717_PredatorExp_acc_raw.txt")
SMR <- RespDataPrep("Data/20190717_PredatorExp_smr_None_raw.txt")
Treat <- RespDataPrep("Data/20190717_PredatorExp_smr_Treat_raw.txt")
Post <- RespDataPrep( "Data/20190717_PredatorExp_blankPost_raw.txt")

RespPrepPlot1(Pre) # Quality control for raw data
RespPrepPlot2(Pre) # Quality control for all slopes
Resp.InputQC(Pre)

Pre <- filter(Pre, Phase != "M4") # Remove erroneous, extra measurement period

RespPrepPlot1(Post) # Quality control for raw data
RespPrepPlot2(Post) # Quality control for all slopes
Resp.InputQC(Post)

Post <- filter(Post, Phase != "M4") # Remove erroneous, extra measurement period

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

# Define background respiration ==

SumPre <- SumPre %>% filter(Phase != "M1") %>% group_by(Type, CH) %>% summarise(avMO2 = mean(O2.mg_h), minMO2 = min(O2.mg_h), maxMO2 = max(O2.mg_h), avR2 = mean(R2), minR2= min(R2), maxR2 = max(R2), Temp.C = mean(Temp.C), sec = mean(sec))

SumPost <- SumPost %>% filter(Phase != "M1") %>% group_by(Type, CH) %>% summarise(avMO2 = mean(O2.mg_h), minMO2 = min(O2.mg_h), maxMO2 = max(O2.mg_h), avR2 = mean(R2), minR2= min(R2), maxR2 = max(R2), Temp.C = mean(Temp.C), sec = mean(sec))

(SumPost$avMO2 - SumPre$avMO2) / (SumPost$sec - SumPre$sec)

BG.corrrection <- left_join(
  (SumPre %>% group_by() %>% select(CH, avMO2, sec) %>% rename(PreMO2 = avMO2, PreSec = sec)), 
  (SumPost %>% group_by() %>% select(CH, avMO2, sec) %>% rename(PostMO2 = avMO2, PostSec = sec)), by = "CH") %>%
  mutate(MO2.Cor = (PostMO2 - PreMO2)/(PostSec - PreSec)) %>%
  select(CH, PreSec, PostSec, MO2.Cor)

(BG.corrrection$PostSec - BG.corrrection$PreSec) * BG.corrrection$MO2.Cor # In order to correct for the background respiration, I need to subtract the time in sec of the phase to get the number of seconds since the beginning of the experiment and then multiply by the correction factor.

# Define experimental respiration rates ===

# . No Treatment analysis: ===

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

# . Treatment SMR analysis: ===

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

# . Combined analysis: ===

Final <- bind_rows(SumSMR, SumTreat)

ggplot(Final, aes(x = CH, y = corO2.mg_h.g, col = Type)) +
  geom_jitter(width = 0.2) +
  theme_minimal()

summary(lm(corO2.mg_h.g ~ Type * as.factor(CH), data = Final))

write.csv(Final, "Jul17.csv")
 
 ### 18 Jul ----

# Define the run of interest:
Pre <- RespDataPrep(  "Data/20190718_PredatorExp_blankPre_raw.txt")
Acc <- RespDataPrep(  "Data/20190718_PredatorExp_acc_raw.txt")
SMR <- RespDataPrep(  "Data/20190718_PredatorExp_smr_None_raw.txt")
Treat <- RespDataPrep("Data/20190718_PredatorExp_smr_Treat_raw.txt")
Post <- RespDataPrep( "Data/20190718_PredatorExp_blankPost_raw.txt")

RespPrepPlot1(Pre) # Quality control for raw data
RespPrepPlot2(Pre) # Quality control for all slopes
Resp.InputQC(Pre)

Pre <- filter(Pre, Phase != "M5") %>% filter(Phase != "M6") # Removing extra, short measurment phase

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

Treat <- filter(Treat, Phase != "M9") # Remove extra measurement phase

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

# Define background respiration ==

SumPre <- SumPre %>% filter(Phase != "M1") %>% group_by(Type, CH) %>% summarise(avMO2 = mean(O2.mg_h), minMO2 = min(O2.mg_h), maxMO2 = max(O2.mg_h), avR2 = mean(R2), minR2= min(R2), maxR2 = max(R2), Temp.C = mean(Temp.C), sec = mean(sec))

SumPost <- SumPost %>% filter(Phase != "M1") %>% group_by(Type, CH) %>% summarise(avMO2 = mean(O2.mg_h), minMO2 = min(O2.mg_h), maxMO2 = max(O2.mg_h), avR2 = mean(R2), minR2= min(R2), maxR2 = max(R2), Temp.C = mean(Temp.C), sec = mean(sec))

(SumPost$avMO2 - SumPre$avMO2) / (SumPost$sec - SumPre$sec)

BG.corrrection <- left_join(
  (SumPre %>% group_by() %>% select(CH, avMO2, sec) %>% rename(PreMO2 = avMO2, PreSec = sec)), 
  (SumPost %>% group_by() %>% select(CH, avMO2, sec) %>% rename(PostMO2 = avMO2, PostSec = sec)), by = "CH") %>%
  mutate(MO2.Cor = (PostMO2 - PreMO2)/(PostSec - PreSec)) %>%
  select(CH, PreSec, PostSec, MO2.Cor)

(BG.corrrection$PostSec - BG.corrrection$PreSec) * BG.corrrection$MO2.Cor # In order to correct for the background respiration, I need to subtract the time in sec of the phase to get the number of seconds since the beginning of the experiment and then multiply by the correction factor.

# Define experimental respiration rates ===

# . No Treatment analysis: ===

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

# . Treatment SMR analysis: ===

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

# . Combined analysis: ===

Final <- bind_rows(SumSMR, SumTreat)

ggplot(Final, aes(x = CH, y = corO2.mg_h.g, col = Type)) +
  geom_jitter(width = 0.2) +
  theme_minimal()

summary(lm(corO2.mg_h.g ~ Type * as.factor(CH), data = Final))

write.csv(Final, "Jul18.csv")
 
 
 ### 19 Jul ----
# Define the run of interest:
Pre <- RespDataPrep(  "Data/20190719_PredatorExp_blankPre_raw.txt")
Acc <- RespDataPrep(  "Data/20190719_PredatorExp_acc_raw.txt")
SMR <- RespDataPrep(  "Data/20190719_PredatorExp_smr_None_raw.txt")
Treat <- RespDataPrep("Data/20190719_PredatorExp_smr_Treat_raw.txt")
Post <- RespDataPrep( "Data/20190719_PredatorExp_blankPost_raw.txt")

RespPrepPlot1(Pre) # Quality control for raw data
RespPrepPlot2(Pre) # Quality control for all slopes
Resp.InputQC(Pre)

Pre <- filter(Pre, Phase != "M4") # Removing extra, short measurment phase

RespPrepPlot1(Post) # Quality control for raw data
RespPrepPlot2(Post) # Quality control for all slopes
Resp.InputQC(Post)

Post <- filter(Post, Phase != "M4") # Removing extra, short measurment phase

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

SumSMR <- anti_join(SumSMR, data.frame("CH" = 8, "Phase" = "M5")) # Remove bad measurements

SumTreat <- RespDataSummary(Treat, "Treat")

SumTreat <- anti_join(SumTreat, data.frame("CH" = 8, "Phase" = "M8")) # Remove bad measurements

SumPost <- RespDataSummary(Post, "Post")

RespSumPlot1(SumPre) # Quality control for slope values
RespSumPlot1(SumAcc) # Quality control for slope values
RespSumPlot1(SumSMR) # Quality control for slope values
RespSumPlot1(SumTreat) # Quality control for slope values
RespSumPlot1(SumPost) # Quality control for slope values

# Define background respiration ==

SumPre <- SumPre %>% filter(Phase != "M1") %>% group_by(Type, CH) %>% summarise(avMO2 = mean(O2.mg_h), minMO2 = min(O2.mg_h), maxMO2 = max(O2.mg_h), avR2 = mean(R2), minR2= min(R2), maxR2 = max(R2), Temp.C = mean(Temp.C), sec = mean(sec))

SumPost <- SumPost %>% filter(Phase != "M1") %>% group_by(Type, CH) %>% summarise(avMO2 = mean(O2.mg_h), minMO2 = min(O2.mg_h), maxMO2 = max(O2.mg_h), avR2 = mean(R2), minR2= min(R2), maxR2 = max(R2), Temp.C = mean(Temp.C), sec = mean(sec))

(SumPost$avMO2 - SumPre$avMO2) / (SumPost$sec - SumPre$sec)

BG.corrrection <- left_join(
  (SumPre %>% group_by() %>% select(CH, avMO2, sec) %>% rename(PreMO2 = avMO2, PreSec = sec)), 
  (SumPost %>% group_by() %>% select(CH, avMO2, sec) %>% rename(PostMO2 = avMO2, PostSec = sec)), by = "CH") %>%
  mutate(MO2.Cor = (PostMO2 - PreMO2)/(PostSec - PreSec)) %>%
  select(CH, PreSec, PostSec, MO2.Cor)

(BG.corrrection$PostSec - BG.corrrection$PreSec) * BG.corrrection$MO2.Cor # In order to correct for the background respiration, I need to subtract the time in sec of the phase to get the number of seconds since the beginning of the experiment and then multiply by the correction factor.

# Define experimental respiration rates ===

# . No Treatment analysis: ===

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

# . Treatment SMR analysis: ===

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

# . Combined analysis: ===

Final <- bind_rows(SumSMR, SumTreat)

ggplot(Final, aes(x = CH, y = corO2.mg_h.g, col = Type)) +
  geom_jitter(width = 0.2) +
  theme_minimal()

summary(lm(corO2.mg_h.g ~ Type * as.factor(CH), data = Final))

write.csv(Final, "Jul19.csv")
 
### 22 Jul ----
# Define the run of interest:
Pre <- RespDataPrep(  "Data/20190722_PredatorExp_blankPre_raw.txt")
Acc <- RespDataPrep(  "Data/20190722_PredatorExp_acc_raw.txt")
SMR <- RespDataPrep(  "Data/20190722_PredatorExp_smr_None_raw.txt")
Treat <- RespDataPrep("Data/20190722_PredatorExp_smr_Treat_raw.txt")
Post <- RespDataPrep( "Data/20190722_PredatorExp_blankPost_raw.txt")

RespPrepPlot1(Pre) # Quality control for raw data
RespPrepPlot2(Pre) # Quality control for all slopes
Resp.InputQC(Pre)

Pre <- filter(Pre, Phase != "M4")

RespPrepPlot1(Post) # Quality control for raw data
RespPrepPlot2(Post) # Quality control for all slopes
Resp.InputQC(Post)

Post <- filter(Post, Phase != "M4")

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

# Define background respiration ===

SumPre <- SumPre %>% filter(Phase != "M1") %>% group_by(Type, CH) %>% summarise(avMO2 = mean(O2.mg_h), minMO2 = min(O2.mg_h), maxMO2 = max(O2.mg_h), avR2 = mean(R2), minR2= min(R2), maxR2 = max(R2), Temp.C = mean(Temp.C), sec = mean(sec))

SumPost <- SumPost %>% filter(Phase != "M1") %>% group_by(Type, CH) %>% summarise(avMO2 = mean(O2.mg_h), minMO2 = min(O2.mg_h), maxMO2 = max(O2.mg_h), avR2 = mean(R2), minR2= min(R2), maxR2 = max(R2), Temp.C = mean(Temp.C), sec = mean(sec))

(SumPost$avMO2 - SumPre$avMO2) / (SumPost$sec - SumPre$sec)

BG.corrrection <- left_join(
  (SumPre %>% group_by() %>% select(CH, avMO2, sec) %>% rename(PreMO2 = avMO2, PreSec = sec)), 
  (SumPost %>% group_by() %>% select(CH, avMO2, sec) %>% rename(PostMO2 = avMO2, PostSec = sec)), by = "CH") %>%
  mutate(MO2.Cor = (PostMO2 - PreMO2)/(PostSec - PreSec)) %>%
  select(CH, PreSec, PostSec, MO2.Cor)

(BG.corrrection$PostSec - BG.corrrection$PreSec) * BG.corrrection$MO2.Cor # In order to correct for the background respiration, I need to subtract the time in sec of the phase to get the number of seconds since the beginning of the experiment and then multiply by the correction factor.

# Define experimental respiration rates ===

# . No Treatment analysis: ===

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

# . Treatment SMR analysis: ===

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

# . Combined analysis: ==

Final <- bind_rows(SumSMR, SumTreat)

ggplot(Final, aes(x = CH, y = corO2.mg_h.g, col = Type)) +
  geom_jitter(width = 0.2) +
  theme_minimal()

summary(lm(corO2.mg_h.g ~ Type * as.factor(CH), data = Final))

write.csv(Final, "Jul22.csv")
 
### 23 Jul ----
# Define the run of interest:
Pre <- RespDataPrep(  "Data/20190723_PredatorExp_blankPre_raw.txt")
Acc <- RespDataPrep(  "Data/20190723_PredatorExp_acc_raw.txt")
SMR <- RespDataPrep(  "Data/20190723_PredatorExp_smr_None_raw.txt")
Treat <- RespDataPrep("Data/20190723_PredatorExp_smr_Treat_raw.txt")
Post <- RespDataPrep( "Data/20190723_PredatorExp_blankPost_raw.txt")

RespPrepPlot1(Pre) # Quality control for raw data
RespPrepPlot2(Pre) # Quality control for all slopes
Resp.InputQC(Pre)

Pre <- filter(Pre, Phase != "M4")

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

Treat <- filter(Treat, Phase != "M9")

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

# Define background respiration ==

SumPre <- SumPre %>% filter(Phase != "M1") %>% group_by(Type, CH) %>% summarise(avMO2 = mean(O2.mg_h), minMO2 = min(O2.mg_h), maxMO2 = max(O2.mg_h), avR2 = mean(R2), minR2= min(R2), maxR2 = max(R2), Temp.C = mean(Temp.C), sec = mean(sec))

SumPost <- SumPost %>% filter(Phase != "M1") %>% group_by(Type, CH) %>% summarise(avMO2 = mean(O2.mg_h), minMO2 = min(O2.mg_h), maxMO2 = max(O2.mg_h), avR2 = mean(R2), minR2= min(R2), maxR2 = max(R2), Temp.C = mean(Temp.C), sec = mean(sec))

(SumPost$avMO2 - SumPre$avMO2) / (SumPost$sec - SumPre$sec)

BG.corrrection <- left_join(
  (SumPre %>% group_by() %>% select(CH, avMO2, sec) %>% rename(PreMO2 = avMO2, PreSec = sec)), 
  (SumPost %>% group_by() %>% select(CH, avMO2, sec) %>% rename(PostMO2 = avMO2, PostSec = sec)), by = "CH") %>%
  mutate(MO2.Cor = (PostMO2 - PreMO2)/(PostSec - PreSec)) %>%
  select(CH, PreSec, PostSec, MO2.Cor)

(BG.corrrection$PostSec - BG.corrrection$PreSec) * BG.corrrection$MO2.Cor # In order to correct for the background respiration, I need to subtract the time in sec of the phase to get the number of seconds since the beginning of the experiment and then multiply by the correction factor.

# Define experimental respiration rates ===

# . No Treatment analysis: ===

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

# . Treatment SMR analysis: ===

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

# . Combined analysis: ==

Final <- bind_rows(SumSMR, SumTreat)

ggplot(Final, aes(x = CH, y = corO2.mg_h.g, col = Type)) +
  geom_jitter(width = 0.2) +
  theme_minimal()

summary(lm(corO2.mg_h.g ~ Type * as.factor(CH), data = Final))

write.csv(Final, "Jul23.csv")
 
 #### COMPILE ----
 ### Main code written by A. Arietta.
 ### Last update 2020 1 Sep by N.R. Sommer
 
require(lme4)
require(lmerTest)
require(tidyverse)
theme_set(theme_minimal())
 
 Jul11 <- read.csv("Jul11.csv")
 Jul12 <- read.csv("Jul12.csv")
 Jul15 <- read.csv("Jul15.csv")
 Jul16 <- read.csv("Jul16.csv")
 Jul17 <- read.csv("Jul17.csv")
 Jul18 <- read.csv("Jul18.csv")
 Jul19 <- read.csv("Jul19.csv")
 Jul22 <- read.csv("Jul22.csv")
 Jul23 <- read.csv("Jul23.csv")
 
 X <- rbind(Jul11, Jul12, Jul15, Jul16, Jul17, Jul18, Jul19, Jul22, Jul23)
 # Bind the clean data into one master file
 
 ## Analyse and correct for temperature ====
 
 X %>% group_by(Type, Resp.DOY, CH) %>% 
   summarise(Temp.C = mean(Temp.C)) %>% 
   group_by(Type) %>% 
   summarise(Temp.C = mean(Temp.C))
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
 
 ## Analyze metabolic rate change by temperature. ====
 
 ggplot(X, aes(x = Temp.C, corO2.mg_h.g)) +
   geom_point() +
   geom_smooth(method = 'lm')
 
 summary(lmer(corO2.mg_h.g ~ Temp.C + Mass.g + (1 | CH) + (1 | Type), data = X))
 
 X <- X %>% mutate(Ex = ifelse(CH == 1 | CH == 2 | CH == 3 | CH == 4, "DF", "CUE")) %>% 
   rename(MO2 = corO2.mg_h.g)

write.csv(X, "masterData_PredResp.csv")
# Master data with assignment of treatment condition (Ex= DF or CUE) and MO2 adjusted for temperature
 
### ANALYSIS ----
require(lmer)
require(parameters)
require(tidyverse)
require(MuMIn)
 
 ## Read in the compiled data ====
 
X <- read.csv("masterData_PredResp.csv") %>%
   mutate(RMSE.dif = RMSE.lm - RMSE.lo) %>% # Calculate a metric of linearity from the difference in loess and linear model fits
   mutate(Tx = ifelse(Type == "SMR", as.character(Type), as.character(Ex))) # Make a new factor that combines the treatment levels and SMR baseline


X <- X %>% 
  mutate(Tx = replace(Tx, Tx == "CUE", "Olfactory")) %>% # Rename treatments for consistency
  mutate(Tx = replace(Tx, Tx == "DF", "Olfactory+Visual")) %>%
  mutate(Tx = fct_relevel(Tx, "SMR", "Olfactory", "Olfactory+Visual")) # Reorder the factors so that they make sense in the figures.
 
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
 
## Remove outliers ====
X <- X %>% filter(RMSE.dif < RMSE.dif.threshold & MO2 > 0) # We remove all of the high outliers with RSME differeces above 0.009
 
 ## Restrict to the lower half of observations for each individual within a trial =====
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
 mod1.sum # Temperature has no strong directional effect but we will leave it in to account for it.
 
 ## Confidence intervals ====
 # This will take a while to run
 
 mod1.ci <- model_parameters(mod1, ci = 0.95, bootstrap = TRUE, iterations = 1000)
 mod1.ci
 mod1.ci$Coefficient # The values are small so we need to output them explicitly
 mod1.ci$CI_low
 mod1.ci$CI_high
 
 ## R2 values ====
 r.squaredGLMM(mod1)
 
 ## FIGURE ====
 
 X2.2 <- X2 %>% 
   mutate(fit.re = predict(mod1, re.form = NULL)) %>% 
   group_by(Tx, Resp.DOY) %>% 
   summarise(MO2 = mean(fit.re)) # This first predicts the value for each trial based on the random effect (essentially an individual-level average). These values are the same for all phases and individuals so we just need to average them down to a single value to get our 9 trial estimates for plotting.
 
 X2.3 <- X2 %>% 
   mutate(fit = predict(mod1, re.form = ~0)) %>% 
   group_by(Tx) %>% 
   summarise(MO2 = mean(fit) + mod1.sum$sigma) # This predicts the value excluding the random effects. These values are the same for all phases and individuals and trials, so again, we just need to average to collapse these into single values. The predictions actually EXCLUDE the random effects rather than average them, so we add the residual sigma from the model to shift the values up to the average for plotting.
 
 X2.2 <- X2.2 %>% mutate(Tx = recode(Tx, SMR = "Baseline"))
 X2.3 <- X2.3 %>% mutate(Tx = recode(Tx, SMR = "Baseline"))
 
 ggplot(data = X2.2,
        aes(x = Tx, y = MO2)) +
    
   # BASELINE
    geom_point(data = X2.2 %>% 
                filter(Tx == "Baseline"),
              aes(), size = 2, col = "darkgreen") +
   geom_half_boxplot(data = X2.2 %>% 
                     filter(Tx == "Baseline"),
                     aes(y=MO2), 
                     position = position_nudge(x = -.6),
                     outlier.shape = NA, 
                     center = TRUE, errorbar.draw = FALSE, width = .2,
                     fill = 'darkgreen', alpha = 0.6) +

   # OLFACTORY
   geom_point(data = X2.2 %>% filter(Tx == "Olfactory"),
              aes(), size = 2, col = "dodgerblue") +
   geom_half_boxplot(data = X2.2 %>% 
                       filter(Tx == "Olfactory"),
                     aes(y=MO2), 
                     position = position_nudge(x = -1.42),
                     outlier.shape = NA, 
                     center = TRUE, errorbar.draw = FALSE, width = .2,
                     fill = 'dodgerblue', alpha = 0.6) +
   
   # OLFACTORY+VISUAL
   geom_point(data = X2.2 %>% filter(Tx == "Olfactory+Visual"),
              aes(), size = 2, col = "darkorange") +
   geom_half_boxplot(data = X2.2 %>% 
                       filter(Tx == "Olfactory+Visual"),
                     aes(y=MO2), 
                     position = position_nudge(x = -2.25),
                     outlier.shape = NA, 
                     center = TRUE, errorbar.draw = FALSE, width = .2,
                     fill = 'darkorange', alpha = 0.6) + 

   # MEAN OVERLAY
   geom_line(data = X2.3 %>% mutate(Tx = recode(Tx, SMR = "Baseline")),
             aes(x = Tx, y = MO2, group = 1), 
             col = "black", 
             size = 1.1) +
   geom_point(data = X2.3, 
               aes(x = Tx, y = MO2, group = 1), 
               col = "black", size = 2) +
    
   coord_cartesian(xlim = c(3,.75)) +
    theme_classic() +
    labs(y = "Metabolic rate (O2 * water volume mL /g)", x = "Predator Cue")
 