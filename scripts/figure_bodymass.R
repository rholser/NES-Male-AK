
library(tidyverse)
library(here)
library(readxl)
library(patchwork)
library(ggdark)

bodymass<-read_xlsx(here("data-raw","NES_Male_foraging.2026v4.xlsx")) |>
  janitor::clean_names()


# Summarize data ----------------------------------------------------------

masschange<-bodymass |>
  filter(!is.na(arrmass)) |>
  select(id, depdate, arr_date, season,depmass,arrmass,massgain,daysatsea,adgain,leangain,mj_day) |>
  pivot_longer(-c(id, depdate,arr_date,massgain,season)) 

bodymass<-bodymass |>
  group_by(id)|>
  arrange(depdate)|>
  mutate(DeployType=c("Departure","Arrival"),
         PercGain=massgain/depmass) |>
  ungroup()

deployTypeSum<-bodymass |>
  group_by(DeployType)|>
  summarise(across(where(is.numeric), mean, na.rm=T))

SeasonSum<-bodymass |>
  group_by(season)|>
  summarise(across(where(is.numeric), mean, na.rm=T))

# Output summary ----------------------------------------------------------
output<-bodymass |>
  select(depmass, arrmass, nrg_mj, mj_day) |>
  set_names(c("depmass_kg","arrmass_kg","MJ_total","MJ_day"))

write_csv(output, here("output", "Summary of mass and energy changes.csv"))


# Plots -------------------------------------------------------------------

GMass<-ggplot(data=masschange |>filter(name=="depmass"| name=="arrmass"), aes(x=reorder(name, value), y=value))+
  geom_line(aes(group=id,color=massgain),linewidth=0.75, alpha=0.65)+
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

G2<-  ggplot(data=masschange |>filter(name=="mj_day"), aes(x=1, y=value))+
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
  
G<-GMass+G2+plot_layout(widths=c(1,0.2))

ggsave(here("figures","Mass Figure.png"), G,width=8, height=6.5)


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
