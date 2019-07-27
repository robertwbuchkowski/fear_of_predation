# Analysis of Cricket

require(tidyverse)

## Analysis of Isopod data from June 2019 ------

IDtable = read_csv("Data/IDtable_June2019.csv") %>% 
  filter(Treatment %in% c("N", "V", "S", "BLANK", "BLANK_SP")) %>%
  rename(Start_0 = Start,End_0 = End) %>%
  gather(-Individual, -Treatment, -WWt, -Cotton, key=Type, value=time) %>%
  separate(Type, into=c("Type", "Rep"), sep="_") %>%
  select(-Cotton) %>%
  spread(key=Type, value=time)

respdata = read_csv("Data/june2019_isopods.csv") %>% select(Time:Flow)

colnames(respdata) = c("Time", "QO2", "FlowV",
                       "CO2", "Temp", "Flow")

respdata = respdata %>% select(-QO2, -FlowV)

vol = 60.79535394
PRESSURE = 25.4*29.6

NN = dim(IDtable)[1]

XX = 1

runeach = function(XX = 1, IDTABLE = IDtable, RESPDATA = respdata){
  select = RESPDATA$Time > as.numeric(IDTABLE[XX,"Start"]) & 
    RESPDATA$Time < as.numeric(IDTABLE[XX,"End"])
  
  selected = RESPDATA[select,]
  
  m1 = lm(CO2~ Time, data=selected)
  slope = coefficients(m1)[2]
  r2 = summary(m1)$adj.r.squared
  TEMP=mean(selected$Temp)
  Weight = as.numeric(IDTABLE[XX,"WWt"])
  WEIGHT = ifelse(Weight==0, 1,Weight)
  volume = ((vol-Weight)*((TEMP+273)/273)*760/PRESSURE)/1000
  resprate = as.numeric(slope*volume*60/(WEIGHT))
  
  # plot(CO2~Time, data=selected, type="l", main = paste0(IDtable[XX,c(1,2,4)]))
  # points(predict(m1)~selected$Time, type="l", col="red")
  
  c(r2 = r2, resprate = resprate)
  
}

output = cbind(IDtable, t(sapply(seq(1, NN, 1), FUN = runeach)))

blank = output %>% filter(!(Treatment %in% c("N", "S", "V"))) %>%
  group_by(Treatment) %>%
  summarize(blank = mean(resprate)) %>%
  mutate(Treatment = ifelse(Treatment=="BLANK", "N", "S"))

blank = rbind(blank, blank[blank$Treatment=="N",])
blank[3,1] = "V"

fdata = output %>% filter((Treatment %in% c("N", "S", "V"))) %>%
  left_join(blank) %>%
  mutate(resp_rate = resprate -blank)

sumfdata = fdata %>% group_by(Individual, Treatment) %>%
  summarize(resp_rate = mean(resp_rate))

fdata %>% group_by(Individual, Treatment) %>%
  summarize(resp_ratesd = sd(resp_rate)) %>%
  left_join(sumfdata) %>%
  mutate(cv = resp_ratesd/resp_rate)

  
fdata %>% ggplot(aes(x=Treatment, y=resp_rate, color=Individual)) +
  geom_point(data = sumfdata, size=6, shape = 18) +
  geom_line(data = sumfdata, aes(group=Individual)) + 
  geom_jitter(aes(shape=Rep)) + theme_classic()


minID = fdata %>% group_by(Individual) %>%
  summarize(Start = min(Start)) %>%
  rename(Start0 = Start)

fdata %>% group_by(Individual, Treatment) %>%
  summarize(Start = min(Start)) %>%
  left_join(sumfdata) %>%
  left_join(minID) %>%
  mutate(Start = Start -Start0) %>%
  ggplot(aes(x=Start, y=resp_rate, color=Treatment, shape=Individual)) + 
  geom_point(size=4) + theme_classic()


## Analsis of Cricket and Isopod Data from July 2019 -----


IDtable = read_csv("Data/IDtable_22July2019.csv") %>%
  filter(!(Individual =="K" & Species =="TRRA")) %>% # get rid of lost isopod
  rename(Start_0 = Start,End_0 = End, WWt = `Wet Weight (g)`) %>%
  gather(-Species:-Treatment, -WWt:-Year, key=Type, value=time) %>%
  separate(Type, into=c("Type", "Rep"), sep="_") %>%
  spread(key=Type, value=time)

respdata = read_csv("Data/19July2019_crickets.csv") %>% 
  select(Time:Flow) %>%
  mutate(Day = 19) %>%
  bind_rows(
    read_csv("Data/20July2019_crickets.csv") %>% 
      select(Time:Flow) %>%
      mutate(Day = 20)
  ) %>%
  bind_rows(
    read_csv("Data/21July2019_crickets.csv") %>% 
      select(Time:Flow) %>%
      mutate(Day = 21)
  ) %>%
  bind_rows(
    read_csv("Data/22July2019_isopod_fear.csv") %>% 
      select(Time:Flow) %>%
      mutate(Day = 22)
  ) %>%
  bind_rows(
    read_csv("Data/isopod_fear_24July2019.csv") %>% 
      select(Time:Flow) %>%
      mutate(Day = 24)
  ) %>%
  bind_rows(
    read_csv("Data/26July2019_spiders.csv") %>% 
      select(Time:Flow) %>%
      mutate(Day = 26)
  )



colnames(respdata) = c("Time", "QO2", "FlowV",
                       "CO2", "Temp", "Flow",
                       "Day", "Flow2")

respdata = respdata %>% select(-QO2, -FlowV, -Flow2)

vol = 60.79535394
# PRESSURE = 25.4*29.6

NN = dim(IDtable)[1]

runeach = function(XX = 1, IDTABLE = IDtable, RESPDATA = respdata){
  select = RESPDATA$Time > as.numeric(IDTABLE[XX,"Start"]) & 
    RESPDATA$Time < as.numeric(IDTABLE[XX,"End"]) &
    RESPDATA$Day == as.numeric(IDTABLE[XX, "Day"])
  
  selected = RESPDATA[select,]
  
  PRESSURE = IDTABLE[XX,"Pressure"]*25.4
  
  m1 = lm(CO2~ Time, data=selected)
  slope = coefficients(m1)[2]
  r2 = summary(m1)$adj.r.squared
  TEMP=mean(selected$Temp)
  Weight = as.numeric(IDTABLE[XX,"WWt"])
  WEIGHT = ifelse(is.na(Weight), 1,Weight)
  WEIGHT0= ifelse(is.na(Weight), 0,Weight)
  volume = ((vol-WEIGHT0)*((TEMP+273)/273)*760/PRESSURE)/1000
  resprate = as.numeric(slope*volume*60/(WEIGHT))
  
  # plot(CO2~Time, data=selected, type="l", main = paste0(IDtable[XX,c(1,2,4)]))
  # points(predict(m1)~selected$Time, type="l", col="red")
  
  c(r2 = r2, resprate = resprate)
  
}

output = cbind(IDtable, t(sapply(seq(1, NN, 1), FUN = runeach)))

blank = output %>% filter(!(Treatment %in% c("B", "C", "N"))) %>%
  group_by(Day) %>%
  summarize(blank = mean(resprate))

fdata = output %>% filter((Treatment %in% c("B", "C", "N"))) %>%
  left_join(blank) %>%
  mutate(resp_rate = resprate -blank)

sumfdata = fdata %>% group_by(Species, Individual, Treatment) %>%
  summarize(resp_rate = mean(resp_rate))

fdata %>% group_by(Species, Individual, Treatment) %>%
  summarize(resp_ratesd = sd(resp_rate)) %>%
  left_join(sumfdata) %>%
  mutate(cv = resp_ratesd/resp_rate)

png("Plot/Fear_results_26July2019.png", width=10, height=5, units="in", res=300)
fdata %>% ungroup() %>% ggplot(aes(x=Treatment, y=resp_rate, color=Individual)) +
  geom_point(data = sumfdata, size=6, shape = 23) +
  geom_line(data = sumfdata, aes(group=Individual)) +
  geom_jitter(aes(shape=Rep), alpha=0.8, width = 0.2, height= 0.2) + 
  theme_classic() +
  facet_wrap(.~Species, scales="free") +
  scale_x_discrete(breaks=c("N", "C", "B"),
                   limits=c("N", "C", "B"),
                   labels=c("Nothing", "Cue", "Visual + Cue"))
dev.off()

library(nlme)
cricketdf = sumfdata %>% filter(Species=="Cricket")
TRRAdf = sumfdata %>% filter(Species=="TRRA")
PHIDdf = sumfdata %>% filter(Species=="Phiddipus")

PHIDdf$Treatment <- factor(PHIDdf$Treatment, levels=c("N","C","B"))

cricketm1 = lme(resp_rate~Treatment, random=~1|Individual, data=cricketdf); summary(cricketm1)
TRRAm1 = lme(resp_rate~Treatment, random=~1|Individual, data=TRRAdf); summary(TRRAm1)
PHIDm1 = lme(resp_rate~Treatment, random=~1|Individual, data=PHIDdf); summary(PHIDm1)

minID = fdata %>% group_by(Species, Individual) %>%
  summarize(Start = min(Start)) %>%
  rename(Start0 = Start)

fdata %>% group_by(Species, Individual, Treatment) %>%
  summarize(Start = min(Start)) %>%
  left_join(sumfdata) %>%
  left_join(minID) %>%
  mutate(Start = Start -Start0) %>%
  ggplot(aes(x=Start, y=resp_rate, color=Individual, shape=Treatment)) + 
  geom_point(size=4) + geom_line(aes(group=Individual)) + theme_classic() +
  facet_grid(.~Species)
