# Analysis of Cricket

require(tidyverse)

IDtable = read_csv("Data/IDtable.csv") %>% 
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
