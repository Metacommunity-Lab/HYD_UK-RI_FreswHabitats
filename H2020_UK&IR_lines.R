#################################################
######### UK metacommunity assembly lines #######
# Freshwaer environment was classified in:
# "ephemeral": <10% of water events along 30 years of sattelite images
# "temporal" : 10%< water events<90%
# "permanent":  water events>90%
################################################
# Performance of species was depenent on freshwater environment
# 10 x 10Kms grid for UK was constructed form satellite images
# Images pixels has a resolution of 10x10mts
# Each cell in the grid report the amount of ephemeral, temporal, and permanent water cover
# ephemeral, temporal, and permanent water from each were considered as different communities in the simulation
# Each environment has a different filter for species performance (estimated in "H2020_UK_Fitness_set_up.R")
################################################
# Number of individuals in each cell is propotional to the are covered b

################################################
M.UK<-M.UK.cent.D50.4  # matrix with cells centroids and cells centrality for a dispersal ability D50=4kms
id.0<-which(M.UK.cent.D50.4$`EF(1).T(2).P(3)`==0)
M.UK<-M.UK.cent.D50.4[-id.0,]
M.dist.UK<-as.matrix(dist(M.UK[,2:3], method = "euclidean"))

Spp.pool<-rep(1,200) # Uniform species abundance distribution is assumed for the pool
mod.temp<-rep(1,nrow(M.UK))
mod.temp[1:10]<-2
Uk.metacomm<-iteration.model(replicas=14, Meta.pool=Spp.pool, m.pool=0.01, Js=M.UK[,6], id.module=M.UK[,7], filter.env=Filter,
                          M.dist=M.dist.UK, D50=4, m.max=1,
                          id.fixed=NULL, D50.fixed=0, m.max.fixed=0, comm.fixed=Spp.pool,
                          Lottery=T, it=100, prop.dead.by.it=0.05, id.obs=1:nrow(M.UK))

UK.out<-resume.out(Uk.metacomm)
names(UK.out)
UK.out.median<-UK.out[[1]]
N<-nrow(M.dist.UK)
S<-UK.out.median[11:(10+N)]
tempo<-cbind(M.UK.cent.D50.4[id.0,],0)
colnames(tempo)[ncol(tempo)]<-"S"
M.UK.S<-cbind(M.UK,S)
M.UK.S<-rbind(M.UK.S,tempo)
which(M.UK.S$grado.in.D50.4==min(M.UK.S$grado.in.D50.4))
M.UK.S<-M.UK.S[-2191,]
which(M.UK.S$grado.in.D50.4==min(M.UK.S$grado.in.D50.4))
M.UK.S<-M.UK.S[-5585,]
which(M.UK.S$grado.in.D50.4==min(M.UK.S$grado.in.D50.4))
M.UK.S<-M.UK.S[-4462,]
4462
5585
id.1<-which(M.UK[,7]==1)
id.2<-which(M.UK[,7]==2)
id.3<-which(M.UK[,7]==3)


#####################################################################################################
id.1<-which(M.UK[,7]==1)
id.2<-which(M.UK[,7]==2)
id.3<-which(M.UK[,7]==3)
id_without_Ephemeral<-c(id.2,id.3)
id_without_temporal<-c(id.1,id.3)
id_without_permanent<-c(id.1,id.2)
id_permanet<-c(id.3)
id_temporal<-c(id.2)
id_ephemeral<-c(id.1)

M.UK.S<-cbind(M.UK,S)
Uk.metacomm_without_Ephi<-iteration.model(replicas=14, Meta.pool=Spp.pool, m.pool=0.01, Js=M.UK[id_without_Ephemeral,6], 
                             id.module=M.UK[id_without_Ephemeral,7], filter.env=Filter[,id_without_Ephemeral],
                             M.dist=M.dist.UK[id_without_Ephemeral,id_without_Ephemeral], D50=4, m.max=1,
                             id.fixed=NULL, D50.fixed=0, m.max.fixed=0, comm.fixed=Spp.pool,
                             Lottery=T, it=100, prop.dead.by.it=0.05, id.obs=1:length(id_without_Ephemeral))

Uk.metacomm_without_Ephi.out<-resume.out(Uk.metacomm_without_Ephi)
S.without_Ephemeral<-Uk.metacomm_without_Ephi.out[[1]][10:(9+length(id_without_Ephemeral))]
M.UK.S<-cbind(M.UK.S,0)
M.UK.S[id_without_Ephemeral,ncol(M.UK.S)]<-S.without_Ephemeral
colnames(M.UK.S)[ncol(M.UK.S)]<-"S.without_Ephemeral"

Uk.metacomm_without_temp<-iteration.model(replicas=14, Meta.pool=Spp.pool, m.pool=0.01, Js=M.UK[id_without_temporal,6], 
                                          id.module=M.UK[id_without_temporal,7], filter.env=Filter[,id_without_temporal],
                                          M.dist=M.dist.UK[id_without_temporal,id_without_temporal], D50=4, m.max=1,
                                          id.fixed=NULL, D50.fixed=0, m.max.fixed=0, comm.fixed=Spp.pool,
                                          Lottery=T, it=100, prop.dead.by.it=0.05, id.obs=1:length(id_without_temporal))
Uk.metacomm_without_temp.out<-resume.out(Uk.metacomm_without_temp)
S.without_temp<-Uk.metacomm_without_temp.out[[1]][10:(9+length(id_without_temporal))]
M.UK.S<-cbind(M.UK.S,0)
M.UK.S[id_without_temporal,ncol(M.UK.S)]<-S.without_temp
colnames(M.UK.S)[ncol(M.UK.S)]<-"S.without_temp"

Uk.metacomm_without_permanent<-iteration.model(replicas=14, Meta.pool=Spp.pool, m.pool=0.01, Js=M.UK[id_without_permanent,6], 
                                          id.module=M.UK[id_without_permanent,7], filter.env=Filter[,id_without_permanent],
                                          M.dist=M.dist.UK[id_without_permanent,id_without_permanent], D50=4, m.max=1,
                                          id.fixed=NULL, D50.fixed=0, m.max.fixed=0, comm.fixed=Spp.pool,
                                          Lottery=T, it=100, prop.dead.by.it=0.05, id.obs=1:length(id_without_permanent))
Uk.metacomm_without_permanent.out<-resume.out(Uk.metacomm_without_permanent)
S.without_perm<-Uk.metacomm_without_permanent.out[[1]][10:(9+length(id_without_permanent))]
M.UK.S<-cbind(M.UK.S,0)
M.UK.S[id_without_permanent,ncol(M.UK.S)]<-S.without_perm
colnames(M.UK.S)[ncol(M.UK.S)]<-"S.without_perm"

Uk.metacomm_without_only_ephemeral<-iteration.model(replicas=14, Meta.pool=Spp.pool, m.pool=0.01, Js=M.UK[id_ephemeral,6], 
                                               id.module=M.UK[id_ephemeral,7], filter.env=Filter[,id_ephemeral],
                                               M.dist=M.dist.UK[id_ephemeral,id_ephemeral], D50=4, m.max=1,
                                               id.fixed=NULL, D50.fixed=0, m.max.fixed=0, comm.fixed=Spp.pool,
                                               Lottery=T, it=100, prop.dead.by.it=0.05, id.obs=1:length(id_ephemeral))
Uk.metacomm_without_only_ephemeral.out<-resume.out(Uk.metacomm_without_only_ephemeral)
S_ephy<-Uk.metacomm_without_only_ephemeral.out[[1]][10:(9+length(id_ephemeral))]
M.UK.S<-cbind(M.UK.S,0)
M.UK.S[id_ephemeral,ncol(M.UK.S)]<-S_ephy
colnames(M.UK.S)[ncol(M.UK.S)]<-"S_only_ephemeral"

Uk.metacomm_without_only_temporal<-iteration.model(replicas=14, Meta.pool=Spp.pool, m.pool=0.01, Js=M.UK[id_temporal,6], 
                                                    id.module=M.UK[id_temporal,7], filter.env=Filter[,id_temporal],
                                                    M.dist=M.dist.UK[id_temporal,id_temporal], D50=4, m.max=1,
                                                    id.fixed=NULL, D50.fixed=0, m.max.fixed=0, comm.fixed=Spp.pool,
                                                    Lottery=T, it=100, prop.dead.by.it=0.05, id.obs=1:length(id_temporal))
Uk.metacomm_without_only_temporal.out<-resume.out(Uk.metacomm_without_only_temporal)
S_temporal<-Uk.metacomm_without_only_temporal.out[[1]][10:(9+length(id_temporal))]
M.UK.S<-cbind(M.UK.S,0)
M.UK.S[id_temporal,ncol(M.UK.S)]<-S_temporal
colnames(M.UK.S)[ncol(M.UK.S)]<-"S_only_temporal"

Uk.metacomm_without_only_permanent<-iteration.model(replicas=14, Meta.pool=Spp.pool, m.pool=0.01, Js=M.UK[id_permanet,6], 
                                                   id.module=M.UK[id_permanet,7], filter.env=Filter[,id_permanet],
                                                   M.dist=M.dist.UK[id_permanet,id_permanet], D50=4, m.max=1,
                                                   id.fixed=NULL, D50.fixed=0, m.max.fixed=0, comm.fixed=Spp.pool,
                                                   Lottery=T, it=100, prop.dead.by.it=0.05, id.obs=1:length(id_permanet))

Uk.metacomm_without_only_permanent.out<-resume.out(Uk.metacomm_without_only_permanent)
S_permanent<-Uk.metacomm_without_only_permanent.out[[1]][10:(9+length(id_permanet))]
M.UK.S<-cbind(M.UK.S,0)
M.UK.S[id_permanet,ncol(M.UK.S)]<-S_permanent
colnames(M.UK.S)[ncol(M.UK.S)]<-"S_only_permanent"
