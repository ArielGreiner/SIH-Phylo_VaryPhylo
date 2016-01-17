#first run 'SIH_varyphylo' and 'rum_simphylo'

#Create dataframe
require(dplyr) 
Data_storage_total<-summarise(group_by(Data_storage, Dispersal, Scale), Mean_SR=mean(SR,na.rm=T), SD_SR=sd(SR,na.rm=T), Mean_Biomass=mean(Biomass,na.rm=T), SD_Biomass=sd(Biomass,na.rm=T), Mean_PD=mean(PD,na.rm=T), SD_PD=sd(PD,na.rm=T), Mean_MPD_abund=mean(MPD_abund,na.rm=T), SD_MPD_abund=sd(MPD_abund,na.rm=T), Mean_MPD_pa=mean(MPD_pa,na.rm=T), SD_MPD_pa=sd(MPD_pa,na.rm=T), Mean_MNTD_abund=mean(MNTD_abund,na.rm=T), SD_MNTD_abund=sd(MNTD_abund,na.rm=T),Mean_MNTD_pa=mean(MNTD_pa,na.rm=T), SD_MNTD_pa=sd(MNTD_pa,na.rm=T), Mean_beta_MPDabund=mean(beta_MPDabund,na.rm=T), SD_beta_MPDabund=sd(beta_MPDabund,na.rm=T), 
Mean_beta_MNTDabund=mean(beta_MNTDabund,na.rm=T), SD_beta_MNTDabund=sd(beta_MNTDabund,na.rm=T), Mean_sesMPD_abund_z = mean(sesMPD_abund_z,na.rm=T), SD_sesMPD_abund_z = sd(sesMPD_abund_z,na.rm=T),
Mean_sesMPD_abund_p = mean(sesMPD_abund_p,na.rm=T), SD_sesMPD_abund_p = sd(sesMPD_abund_p,na.rm=T), Mean_sesMNTD_abund_z = mean(sesMNTD_abund_z,na.rm=T), SD_sesMNTD_abund_z = sd(sesMNTD_abund_z,na.rm=T), Mean_sesMNTD_abund_p = mean(sesMNTD_abund_p,na.rm=T), SD_sesMNTD_abund_p = sd(sesMNTD_abund_p,na.rm=T), sum_phylogeven_mpd = sum(phylogeven_mpd), sum_phylogeven_mntd = sum(phylogeven_mntd), sum_phylogcluster_mpd = sum(phylogcluster_mpd), sum_phylogcluster_mntd = sum(phylogcluster_mntd), sum_pe_pc_mpd = sum(phylogeven_mpd + phylogcluster_mpd), sum_pe_pc_mntd = sum(phylogeven_mntd + phylogcluster_mntd))



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