#first run 'SIH_varyphylo' and 'rum_simphylo'

#Create dataframe
require(dplyr) 
Data_storage_total<-summarise(group_by(Data_storage, Dispersal, Scale), Mean_SR=mean(SR,na.rm=T), SD_SR=sd(SR,na.rm=T), Mean_Biomass=mean(Biomass,na.rm=T), SD_Biomass=sd(Biomass,na.rm=T), Mean_PD=mean(PD,na.rm=T), SD_PD=sd(PD,na.rm=T), Mean_MPD_abund=mean(MPD_abund,na.rm=T), SD_MPD_abund=sd(MPD_abund,na.rm=T), Mean_MPD_pa=mean(MPD_pa,na.rm=T), SD_MPD_pa=sd(MPD_pa,na.rm=T), Mean_MNTD_abund=mean(MNTD_abund,na.rm=T), SD_MNTD_abund=sd(MNTD_abund,na.rm=T),Mean_beta_MPDabund=mean(beta_MPDabund,na.rm=T), SD_beta_MPDabund=sd(beta_MPDabund,na.rm=T), 
Mean_beta_MNTDabund=mean(beta_MNTDabund,na.rm=T), SD_beta_MNTDabund=sd(beta_MNTDabund,na.rm=T))

#replace 'SR' in the function below with PD/MPDpa/MPDabund/MNTDabund
require(ggplot2) #need to define x and y within aes in ggplot
ggplot(Data_storage_total,aes(x=Dispersal,y=Mean_SR,color=Scale,group=Scale,fill=Scale,alpha=0.1))+
  geom_line(size=2)+ #plots data as lines
  geom_ribbon(aes(ymin=Mean_SR-SD_SR,ymax=Mean_SR+SD_SR),width=0.1)+
  scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines
  
  
  #Plot Biomass
ggplot(Data_storage_total,aes(x=Dispersal,y=Mean_Biomass))+ #defines the variables that you are plotting
  geom_line(size=2)+ #plots data as lines
  facet_grid(Scale~.,scale="free")+ #separates out local and regional into two panels
  geom_errorbar(aes(ymin=Mean_Biomass-SD_Biomass,ymax=Mean_Biomass+SD_Biomass),width=0.1)+
  #facet_grid(ReplicateNume~.,scale='free')
  #geom_errorbar(data = Function_rates_sd)+
  scale_x_log10(breaks=DispV)+ #sets x axis to log10 scale
  theme_bw(base_size = 18)+ #gets rid of grey background
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #removes grid lines
  
#example of how to plot data colour-coded and connected by dispersal
ggplot(Data_storage_total,aes(x=Mean_SR,y=Mean_Biomass,color=factor(Dispersal),group=Scale))+
geom_point()+ #plots data as points
geom_path()+
facet_grid(Scale~.,scale="free")+
theme_bw(base_size = 18)+ #gets rid of grey background
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#AND MANY MORE WHOO