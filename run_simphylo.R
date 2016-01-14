nspecies<-7 #the number of species
npatches<-10 #the number of patches
nreplicates<-2 #number of replicates
nfunctions<-1

DispV<-c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1) #the dispersal rates 

Data_storage<-data.frame(SR=NA,Biomass=NA,Biomass_CV=NA,PD=NA,MPD_abund=NA,MPD_pa=NA,MNTD_abund=NA,beta_MPDabund =NA, beta_MNTDabund=NA, sesMPD_abund = NA, sesMNTDabund = NA,Dispersal=rep(DispV,each=nreplicates),ReplicateNum=factor(1:nreplicates),Scale=rep(c("Local","Regional"),each=length(DispV)*nreplicates)) #building the data frame

MTraits<-t(matrix(1,nspecies,nfunctions))

for(j in 1:nreplicates){
#runs the SIH model at all dispersal rates in DispV and saves the abundances and productivity in a list
set.seed(j)
eff_sd <- 0.005
eff_mean <- 0.2
SIH_data<-sapply(DispV,SIH_function,species=nspecies,patches=npatches,eff_vary=T,eff_sd=eff_sd, eff_mean = eff_mean)


#removed the SES measures b/c on 1.12.2016 they weren't working 
for(i in 1:length(DispV)){
 #BIOMASS
  LFunc_rate1<-apply(SIH_data[["Abund",i]],1,function(M){MTraits%*%M}) #calculates the amount of each function in each patch based on the abundances from the SIH model
  LFunc_rate<-array(t(LFunc_rate1),dim=c(nrow(SIH_data[["Abund",1]]),nfunctions,npatches)) #arranges the functions into an array that is easier to use
  Data_storage$Biomass[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(apply(LFunc_rate,3,colMeans))
  Data_storage$Biomass[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-colMeans(apply(LFunc_rate,2,rowSums)) #calculates and saves the mean regoinal rate for each function
  
	#calculate species richness at the local and regional scale	
  Data_storage$SR[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"] <-mean(rowMeans(apply(SIH_data[["Abund",i]]>0,3,rowSums))) 
  Data_storage$SR[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mean(rowSums(apply(SIH_data[["Abund",i]],2,rowMeans)>0)) #regional SR of all species in each time step
  #At the local scale...	
  com_data<-t(SIH_data[["Abund",i]][400,,])
  colnames(com_data)<-1:nspecies
  Data_storage$PD[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(pd(com_data,SIH_data[["phylo",i]])$PD) 
  Data_storage$MPD_abund[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(mpd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = T)) 
  Data_storage$MPD_pa[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(mpd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = F))
  Data_storage$MNTD_abund[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Local"]<-mean(mntd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = T))

  #beta measures
  colnames(com_data) <- paste('species', 1:nspecies)
     interspec_mat <- cophenetic(SIH_data[["phylo",i]])
     colnames(interspec_mat) <- paste('species', 1:7)
     rownames(interspec_mat) <- paste('species', 1:7)
  Data_storage$beta_MPDabund[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"] <- mean(comdist(com_data,interspec_mat,abundance.weighted=TRUE))
Data_storage$beta_MNTDabund[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"] <- mean(comdistnt(com_data,interspec_mat,abundance.weighted=TRUE,exclude.conspecifics=TRUE))
#^gives weird warning messages...

  #At the regional scale
  com_data<-matrix(colSums(t(SIH_data[["Abund",i]][400,,])),1,nspecies)
  colnames(com_data)<-1:nspecies
  Data_storage$PD[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-pd(com_data,SIH_data[["phylo",i]])$PD
  Data_storage$MPD_abund[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mpd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = T)
  Data_storage$MPD_pa[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mpd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = F)
  Data_storage$MNTD_abund[Data_storage$Dispersal==DispV[i] & Data_storage$ReplicateNum==j & Data_storage$Scale == "Regional"]<-mntd(com_data,cophenetic(SIH_data[["phylo",i]]),abundance.weighted = T) 
  
#storing the phylogenies and their mpd and mntd values, ignoring the real community
fakecom <- matrix(Data = rep(1,nspecies), ncol = nspecies, nrow = 1)
colnames(fakecom)<-1:nspecies
mpdphylo <- matrix(NA, nrow = length(DispV), ncol = nreplicates)
mntdphylo <- matrix(NA, nrow = length(DispV), ncol = nreplicates)
#mpdphylo <- rep(NA,nreplicates)
#mntdphylo <- rep(NA,nreplicates)
phlgs <- array(data = NA, dim=c(nspecies, nspecies, nreplicates, length(DispV)))
placeholder = cophenetic(rcoal(nspecies))
colnames(placeholder)<-1:nspecies
rownames(placeholder)<-1:nspecies
mpdphylo[i,j] = mpd(fakecom,placeholder, abundance.weighted = F)
mntdphylo[i,j] = mntd(fakecom,placeholder, abundance.weighted = F)
phlgs[,,i,j] <- placeholder
   
  }
  }
