
# Estimate energy intake and biomass consumption 

library(here)
library(tidyverse)
library(readxl)

# Read in parameter values ------------------------------------------------
parms<-read_xlsx(here("data","Bioenergetic parameters.xlsx"))

for (i in 1:nrow(parms)){
  assign(parms$Parameter[i],as.numeric(parms$Mean[i]))
  assign(paste(parms$Parameter[i],"SD",sep="_"),as.numeric(parms$SD[i]))
  assign(paste(parms$Parameter[i],"Low",sep="_"),as.numeric(parms$Low[i]))
  assign(paste(parms$Parameter[i],"Upr",sep="_"),as.numeric(parms$Upp[i]))
}

# Create data frame --------------------------------------------------------
numseals<-10000 # Just a number higher than the population to draw from

set.seed(38383)
males<-tibble(sealid=1:numseals, mass=NA, pbtrip=NA,pmtrip=NA,bmr=NA,fmrmult=NA,fmr=NA,
              me=ME, hif=NA, egainpb=NA, egainpm=NA, netdailyepb=NA, netdailyepm=NA,
              medailyepb=NA, medailyepm=NA, grossdailyepb=NA, grossdailyepm=NA,
              grossepb=NA, grossepm=NA) |>
  mutate(mass=truncnorm::rtruncnorm(numseals, mean=bodyMass, sd=bodyMass_SD, a=bodyMass_Low, b=bodyMass_Upr),
         pbtrip=truncnorm::rtruncnorm(numseals,mean=tripDur,sd=tripDur_SD, a=tripDur_Low, b=tripDur_Upr),
         pmtrip=truncnorm::rtruncnorm(numseals,mean=tripDur,sd=tripDur_SD, a=tripDur_Low, b=tripDur_Upr),
         bmr=0.293*mass^0.75, fmr=runif(numseals, min=FMR_Low, max=FMR_Upr)*mass^0.75,
         hif=runif(numseals, min=HIF_Low, max=HIF_Upr), 
         egainpb=truncnorm::rtruncnorm(numseals, a=eGain_Low, b=eGain_Upr,mean=eGain, sd=eGain_SD),
         egainpm=truncnorm::rtruncnorm(numseals, a=eGain_Low, b=eGain_Upr,mean=eGain, sd=eGain_SD),
         netdailyepb=egainpb+fmr,  netdailyepm=egainpm+fmr,
         medailyepb=netdailyepb/(1-hif),medailyepm=netdailyepm/(1-hif),
         grossdailyepb=medailyepb/me,grossdailyepm=medailyepm/me,
         grossepb=grossdailyepb*pbtrip, grossepm=grossdailyepm*pmtrip, grosse=grossepb+grossepm)
         

# Create population -------------------------------------------------------

# right now replicates random draws based on two different population sizes using 100 reps

set.seed(5767)
malespop<-males |>
  nest()|>
  uncount(2) |>
  mutate(popsize=c(2763,3316),population=map2(data,popsize,~.x|>infer::rep_slice_sample(.y,reps=100, replace = TRUE)))|>
  unnest(population) |>
  select(-data)
  
# Summary statistics used to populate Table 2
malespopSum <- malespop |>
  group_by(popsize,replicate)|>
  summarise(popGrossE=sum(grosse),dailyfmr=mean(fmr),dailyfmrsd=sd(fmr),dailyegain=mean((egainpb+egainpm)/2),
            dailyegainsd=sd((egainpb+egainpm)/2),dailynet=mean((netdailyepb+netdailyepm)/2),
            dailynetsd=sd((netdailyepb+netdailyepm)/2),grossdaily=mean((grossdailyepb+grossdailyepm)/2),
            dailygrosssd=sd((grossdailyepb+grossdailyepm)/2), propfmrgrossdaily=dailyfmr/grossdaily,
            propegaingrossdaily=dailyegain/grossdaily,
            yearlyfmr=mean(fmr*pbtrip+fmr*pmtrip),yearlyfmrsd=sd(fmr*pbtrip+fmr*pmtrip),
            yearlyegain=mean(egainpb*pbtrip+egainpm*pmtrip),yearlyegainsd=sd(egainpb*pbtrip+egainpm*pmtrip),
            yearlygross=mean(grossdailyepb*pbtrip+grossdailyepm*pmtrip),yearlygrosssd=sd(grossdailyepb*pbtrip+grossdailyepm*pmtrip))|> 
  ungroup()

# Population summary
malespopSum|>group_by(popsize)|> summarise(MPop=mean(popGrossE)/1000, SDPPop=sd(popGrossE)/1000)

# Individual smmary
summary(malespopSum$grossdaily)
summary(malespopSum$dailygrosssd)
summary(malespopSum$propfmrgrossdaily)
summary(malespopSum$propegaingrossdaily)
summary(malespopSum$yearlygross)
summary(malespopSum$yearlygrosssd)


# Save output -------------------------------------------------------------

saveRDS(malespop, here("output","Bioenergetic model output.rds"))

# Jaw accelerometer  seals-----------------------------------------------------

jaw<-data.frame(seal=c("J914", "J916", "G841"),mass=c(1241.5,1429.4, 1226.9), egain=c(9408,7945,3928), trip=c(127,116,127),
                feedingevents=c(0.56, 0.77, 0.55), focalforagedays=c(53,56,52), foragedays=)|>
  mutate(fmr=trip*(FMR*mass^0.75), me=(fmr+egain)*(1-HIF), gross=me/ME, focalnrg=gross*feedingevents/focalforagedays)
