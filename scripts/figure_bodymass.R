
library(tidyverse)
library(here)
library(readxl)
library(patchwork)
library(ggdark)
bodymass<-read_xlsx(here("data","NES_Male_foraging.2026v3.xlsx"))


# Summarize data ----------------------------------------------------------


masschange<-bodymass |>
  filter(!is.na(ARRMASS)) |>
  select(ID, DEPDATE, `ARR DATE`, SEASON,DEPMASS,ARRMASS,MASSGAIN,DAYSATSEA,ADGAIN,LEANGAIN,`Mj/day`) |>
  pivot_longer(-c(ID, DEPDATE,`ARR DATE`,MASSGAIN,SEASON)) 

bodymass<-bodymass |>
  group_by(ID)|>
  arrange(DEPDATE)|>
  mutate(DeployType=c("Departure","Arrival"),
         PercGain=MASSGAIN/DEPMASS)

deployTypeSum<-bodymass |>
  group_by(DeployType)|>
  summarise(across(where(is.numeric), mean, na.rm=T))

SeasonSum<-bodymass |>
  group_by(SEASON)|>
  summarise(across(where(is.numeric), mean, na.rm=T))

# Plots -------------------------------------------------------------------

GMass<-ggplot(data=masschange |>filter(name=="DEPMASS"| name=="ARRMASS"), aes(x=reorder(name, value), y=value))+
  geom_line(aes(group=ID,color=MASSGAIN),linewidth=0.5, alpha=0.65)+
  paletteer::scale_color_paletteer_c("grDevices::Teal", direction=-1, name="Mass gain (kg)")+
  ggdist::stat_pointinterval()+
ggthemes::theme_few()+  
  scale_x_discrete(labels=c("Departure","Arrival"),expand=c(0,0))+
 coord_cartesian(xlim=c(0.95,2.05))+
  theme(panel.grid.major = element_line(linewidth=0.3, colour="grey", linetype = "dashed"),
        axis.text.y = element_text(size = 16),
        axis.text.x=element_text(size=16),
        axis.title = element_text(size = 18),
        axis.ticks.x = element_blank(),
        legend.position =c(0.85,0.145),
        legend.background = element_blank())+
  xlab(NULL)+
  ylab("Body mass (kg)")

G2<-  ggplot(data=masschange |>filter(name=="Mj/day"), aes(x=1, y=value))+
  ggdist::stat_pointinterval()+
  #geom_boxplot()+
  labs(y="Energy gain per day at sea (MJ)", x=NULL)+
  ggthemes::theme_few()+  
  theme(panel.grid.major = element_line(linewidth=0.3, colour="grey", linetype = "dashed"),
        axis.text.y = element_text(size = 16),
        axis.text.x=element_blank(),
        axis.title = element_text(size = 18),
        axis.ticks.x = element_blank(),
        legend.position =c(0.85,0.145))
  
GMass+G2+plot_layout(widths=c(1,0.2))

GMass_epoc<-ggplot(data=masschange |>filter(name=="DEPMASS"| name=="ARRMASS"), aes(x=reorder(name, value), y=value))+
  geom_line(aes(group=ID,color=MASSGAIN),linewidth=2, alpha=0.65)+
  paletteer::scale_color_paletteer_c("grDevices::Teal", direction=-1, name="Mass gain (kg)")+
  ggdist::stat_pointinterval(interval_size_range = c(2, 2))+
  theme_bw()+ 
  scale_x_discrete(labels=c("Departure","Arrival"),expand=c(0,0))+
  coord_cartesian(xlim=c(0.95,2.05))+
  theme(panel.grid.major = element_line(linewidth=0.3, colour="grey", linetype = "dashed"),
        axis.text.y = element_text(size = 22),
        axis.text.x=element_text(size=22),
        axis.title = element_text(size = 22),
        axis.ticks.x = element_blank(),
        legend.position =c(0.85,0.15),
        legend.background = element_blank(), legend.text=element_text(size=18), legend.title=element_text(size=22),
        plot.margin = margin(t = 0, r = 0.25, b = 0.25 ,l = 0.25, unit = "in"))+
  xlab(NULL)+
  ylab("Body mass (kg)")

ggsave(here("figures","EPOC_Mass Figure.png"), width=8, height=6.5)
