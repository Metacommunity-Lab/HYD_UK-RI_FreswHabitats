#setwd("~/Documents/Analisis_Datos/A_H2020_Networks/EUROPA_EN_10X10")
#load("~/Documents/Analisis_Datos/A_H2020_Networks/EUROPA_EN_10X10/Europa_J.RData")
#save.image("~/Documents/Analisis_Datos/A_H2020_Networks/EUROPA_EN_10X10/Paisaje_UK.RData")
save.image("~/Desktop/H2020_UK/UK_space.RData")

#setwd("~/Desktop/H2020_UK")

head(M.TOTAL.corregida)
#M.TOTAL.corregida[-127670,]->M.Europa
####################################
#################################################
## ACA ARRANCA ESTIMACION DE JJs
####### AGREGA JJ con J=a*A^b1 con a = a=Jmax/(Amax^b2)
area.max.europa<-max(M.Europa$FRESHWATER) ### areas en porcentaje
area.max.europa
J.max<-400
b.ef<-0.5
J.efimero<-ceiling((J.max/(area.max.europa)^b.ef)*M.Europa$FRESHWATER.efimeros^b.ef)
J.temporal<-ceiling((J.max/(area.max.europa)^b.ef)*M.Europa$FRESHWATER.temporales^b.ef)
J.permanente<-ceiling((J.max/(area.max.europa)^b.ef)*M.Europa$FRESHWATER.permanentes^b.ef)
J.freshwater<-ceiling((J.max/(area.max.europa)^b.ef)*M.Europa$FRESHWATER^b.ef)

M.Europa.J<-data.frame(M.Europa,J.freshwater,J.efimero,J.temporal,J.permanente)
head(M.Europa.J) ### esta MATRIZ ES CON los jj por sistema y total
dim(M.Europa.J)

#save.image("~/Documents/Analisis_Datos/A_H2020_Networks/EUROPA_EN_10X10/Europa_J.RData")

unique(M.Europa.J$ISO3_CODE)
a1<-M.Europa.J[which(M.Europa.J$ISO3_CODE=="GBR"),]
a2<-M.Europa.J[which(M.Europa.J$ISO3_CODE=="IRL"),]
dim(a1)[1]+dim(a2)[1]
M.UK<-rbind(a2,a1)
head(M.UK)
dim(M.UK)

### plot
par(mar=c(2.5,2.5,1.5,1.5),bty="l",cex=1.8)
plot(M.Europa.J$CENTROID_Y~M.Europa.J$CENTROID_X,cex=0.1,col="light gray",pch=15)
points(M.UK$CENTROID_Y~M.UK$CENTROID_X,cex=0.1,col="orange",pch=15)
unique(M.UK$ISO3_CODE)

############### MATRIZ DISTANCIA DE LATU.LONG A DISTANCIA EN METROS
library(geosphere) # Charging package to caluclate distances
library(sp)
library(stars)
coord.xy<-cbind(M.UK$CENTROID_X,M.UK$CENTROID_Y)
colnames(coord.xy)<-c("CENTROID_X "," CENTROID_Y")
FreshW_coord <- st_as_sf(as.data.frame(coord.xy),coords = c("CENTROID_X "," CENTROID_Y")) # Transforming to sf format
FreshW_coord_spatial <- as(FreshW_coord, "Spatial") # Convert to "old" format to caluclate to distances
aa_en_metros<- distm(FreshW_coord_spatial,FreshW_coord_spatial, fun = distGeo) #### Calculate the distance matrix (in METRES!!!) from a lat/long values 

dim(M.UK); dim(aa_en_metros)
colnames(aa_en_metros)<-M.UK$PageName
rownames(aa_en_metros)<-M.UK$PageName
aa_en_km<-aa_en_metros/1000
aa_en_km[1:10,1:10]

### aa_en_kms la matriz de distancia con buffer de Suecia
# correcion en MATRIZ DISTANCIA para mij =1 de centr-centr contiguos, y mij=0.5 centr-centr en diagonal
aa_en_km.min_10<-ifelse(aa_en_km<10,10,aa_en_km)
## corregida para que de centroide-centroide CONTIGUO mij = 1
Distancia.UK<-(aa_en_km.min_10-10) ## matriz corregida; distancia 0 consigo mismo y con vecino
Distancia.UK[1:10,1:10] ## matriz de distancia corregida en kilometros

################################################################################################################################################
############################################################################################
####### D50.4; ##### opcion donde D50 es 4 KM (a celda diagonal de la matriz)
############################################################################################
################################################################################################################################################
rm(D50)
D50<-4
b.mig=(-log(.5)/(D50)) ## -log(y)/D50 #### D50= diagonal cetr-cent + distancia en KM
par(mar=c(4,4,.5,.5),cex=1.4,bty="l")
plot(exp(-b.mig*(0:200))~I(0:200),ylim=c(0,1),type="l",lwd=4,col="dark green",ylab="Dispersal probability",xlab="Distance (Km)")
abline(h=0.5)
exp(-b.mig*4)
b.mig_D50.4<-b.mig
rm(b.mig)

##########################################
## exp(-b^Dij)*Ji
dim(Distancia.UK); dim(M.UK)
exp(-b.mig_D50.4*Distancia.UK)*M.UK$J.efimero->mm.efimero.D50.4 ## 1
exp(-b.mig_D50.4*Distancia.UK)*M.UK$J.temporal->mm.temporal.D50.4 ## 2
exp(-b.mig_D50.4*Distancia.UK)*M.UK$J.permanente->mm.permanente.D50.4 ## 3
exp(-b.mig_D50.4*Distancia.UK)*M.UK$J.freshwater->mm.fw.D50.4 ## 
##########################################
id.fw<-which(M.UK$FRESHWATER!=0) ### 4384 son agua
id.efimero<-which(M.UK$FRESHWATER.efimeros!=0) ### 2220 fw efimeros
id.temporal<-which(M.UK$FRESHWATER.temporales!=0) ### 4527 fw temporales
id.permanentes<-which(M.UK$FRESHWATER.permanentes!=0) ### 4733 fw permanentes
######################################################################

#########
####
detach("package:sna", unload = TRUE)
library(igraph) ### PARA el betwenness y closeness estimados en igraph los weights son considerados distancias, por eso ponemos 1/dist como weight
##  weighted betweenness centrality per freshwater
##### bc.efimero.D50.4
graph.adjacency(1/(mm.efimero.D50.4[id.efimero,id.efimero]),mode="directed",weighted = TRUE)->g0.efimero.D50.4
betweenness(g0.efimero.D50.4,directed = TRUE)->bc.efimero.D50.4 ### SUS WEIGHTS SON DISTANCIAS; ACA el peso deberia ser 1/mij y esta incluido en la construccion del grafo
rm(g0.efimero.D50.4)

graph.adjacency(1/(mm.temporal.D50.4[id.temporal,id.temporal]),mode="directed",weighted = TRUE)->g0.temporal.D50.4
betweenness(g0.temporal.D50.4,directed = TRUE)->bc.temporal.D50.4
rm(g0.temporal.D50.4)

save.image("~/Documents/Analisis_Datos/A_H2020_Networks/EUROPA_EN_10X10/Paisaje_UK.RData")

graph.adjacency(1/(mm.permanente.D50.4[id.permanentes,id.permanentes]),mode="directed",weighted = TRUE)->g0.permanente.D50.4
betweenness(g0.permanente.D50.4,directed = TRUE)->bc.permanente.D50.4
rm(g0.permanente.D50.4)

graph.adjacency(1/(mm.fw.D50.4[id.fw,id.fw]),mode="directed",weighted = TRUE)->g0.fw.D50.4
betweenness(g0.fw.D50.4,directed = TRUE)->bc.fw.D50.4
rm(g0.fw.D50.4)

save.image("~/Documents/Analisis_Datos/A_H2020_Networks/EUROPA_EN_10X10/Paisaje_UK.RData")

######################################################################
##### GRADO
detach("package:igraph", unload = TRUE)
library("sna") #### PARA EL degree estimada en sna, los weights son considerados pesos
degree(mm.efimero.D50.4[id.efimero,id.efimero],gmode="digraph",cmode="indegree")->grado.in.efimero.D50.4
degree(mm.efimero.D50.4[id.efimero,id.efimero],gmode="digraph",cmode="outdegree")->grado.out.efimero.D50.4
degree(mm.efimero.D50.4[id.efimero,id.efimero],gmode="digraph",cmode="freeman")->grado.all.efimero.D50.4
colnames(M.UK)
head(cbind(M.UK[id.efimero,c(1:9,c(10,11),c(14,15))],grado.in.efimero.D50.4,grado.out.efimero.D50.4,grado.all.efimero.D50.4))

degree(mm.temporal.D50.4[id.temporal,id.temporal],gmode="digraph",cmode="indegree")->grado.in.temporal.D50.4
degree(mm.temporal.D50.4[id.temporal,id.temporal],gmode="digraph",cmode="outdegree")->grado.out.temporal.D50.4
degree(mm.temporal.D50.4[id.temporal,id.temporal],gmode="digraph",cmode="freeman")->grado.all.temporal.D50.4
head(cbind(M.UK[id.temporal,c(1:9,c(10,12),c(14,16))],grado.in.temporal.D50.4,grado.out.temporal.D50.4,grado.all.temporal.D50.4))

degree(mm.permanente.D50.4[id.permanentes,id.permanentes],gmode="digraph",cmode="indegree")->grado.in.permanente.D50.4
degree(mm.permanente.D50.4[id.permanentes,id.permanentes],gmode="digraph",cmode="outdegree")->grado.out.permanente.D50.4
degree(mm.permanente.D50.4[id.permanentes,id.permanentes],gmode="digraph",cmode="freeman")->grado.all.permanente.D50.4
head(cbind(M.UK[id.permanentes,c(1:9,c(10,13),c(14,17))],grado.in.permanente.D50.4,grado.out.permanente.D50.4,grado.all.permanente.D50.4))

degree(mm.fw.D50.4[id.fw,id.fw],gmode="digraph",cmode="indegree")->grado.in.fw.D50.4
degree(mm.fw.D50.4[id.fw,id.fw],gmode="digraph",cmode="outdegree")->grado.out.fw.D50.4
degree(mm.fw.D50.4[id.fw,id.fw],gmode="digraph",cmode="freeman")->grado.all.fw.D50.4
head(cbind(M.UK[id.fw,c(1:9,c(10),c(14))],grado.in.fw.D50.4,grado.out.fw.D50.4,grado.all.fw.D50.4))

save.image("~/Documents/Analisis_Datos/A_H2020_Networks/EUROPA_EN_10X10/Paisaje_UK.RData")

##################################################
####################################################################################################
colnames(M.UK)
M.Efimeros.D50.4<-cbind(M.UK[id.efimero,c(1:9,c(11),c(15))],rep(1,length(id.efimero)),grado.in.efimero.D50.4,grado.out.efimero.D50.4,grado.all.efimero.D50.4,bc.efimero.D50.4)
M.Temporal.D50.4<-cbind(M.UK[id.temporal,c(1:9,c(12),c(16))],rep(2,length(id.temporal)),grado.in.temporal.D50.4,grado.out.temporal.D50.4,grado.all.temporal.D50.4,bc.temporal.D50.4)
M.Permanentes.D50.4<-cbind(M.UK[id.permanentes,c(1:9,c(13),c(17))],rep(3,length(id.permanentes)),grado.in.permanente.D50.4,grado.out.permanente.D50.4,grado.all.permanente.D50.4,bc.permanente.D50.4)
M.FW.D50.4.UK<-(cbind(M.UK[id.fw,c(1:9,c(10),c(14))],rep(100,length(id.fw)),grado.in.fw.D50.4,grado.out.fw.D50.4,grado.all.fw.D50.4,bc.fw.D50.4))

colnames(M.Efimeros.D50.4)[10:16]
colnames(M.Efimeros.D50.4)[10:16]<-c("FW.Area","J.uk","EF(1).T(2).P(3)","grado.in.D50.4",
                                     "grado.out.D50.4","grado.all.D50.4",
                                     "bc.D50.4")
colnames(M.Efimeros.D50.4)[10:16]->colnames(M.Temporal.D50.4)[10:16]
colnames(M.Efimeros.D50.4)[10:16]->colnames(M.Permanentes.D50.4)[10:16]
colnames(M.Efimeros.D50.4)[10:16]->colnames(M.FW.D50.4.UK)[10:16]

M.LAND<-cbind(M.UK[-id.fw,c(1:9)],0,0,0,0,0,0,0)
colnames(M.LAND)[10:16]<-colnames(M.Efimeros.D50.4)[10:16]
head(M.LAND)

M.final.UK<-rbind(M.Efimeros.D50.4,M.Temporal.D50.4,M.Permanentes.D50.4,M.LAND)
head(M.final.UK)
dim(M.final.UK) ## 11721
unique(M.final.UK$ISO3_CODE)
unique(M.final.UK$`EF(1).T(2).P(3)`)

M.final.FW.UK<-rbind(M.FW.D50.4.UK,M.LAND)
unique(M.final.FW.UK$ISO3_CODE)
unique(M.final.FW.UK$`EF(1).T(2).P(3)`)


###########################################################################################
            ####### ####### ####### PLOTS ####### ####### #######
###########################################################################################
## M.final.UK esta es la matriz completa dierenciando cuerpos de agua y tierra
## M.final.FW.UK esta es la matriz de agua (toda junta) y tiera 
M.final.UK->MM
M.final.FW.UK->MM.FW


#### In_Degree_MassEffect #### D50=4Km
#mtext("In-degree (Mass effect) \n", outer = TRUE, cex = 1.2)
par(mar=c(2.5,2.5,1.5,1.5),bty="l",cex=1.8)
par(mfrow=c(3,4),oma = c(0, 0, 0, 0))
plot(MM$CENTROID_Y~MM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Ephemeral (4Km)")
ii <- cut(MM$grado.in.D50.4[which(MM$`EF(1).T(2).P(3)`==1)], 
          breaks = seq(min(MM$grado.in.D50.4[which(MM$`EF(1).T(2).P(3)`==1)]), 
                       max(MM$grado.in.D50.4[which(MM$`EF(1).T(2).P(3)`==1)]), len = 600), include.lowest = TRUE)
points(MM$CENTROID_Y[which(MM$`EF(1).T(2).P(3)`==1)]~MM$CENTROID_X[which(MM$`EF(1).T(2).P(3)`==1)],cex=0.4,
       col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
text(-8,60,"IN DEGREE\nMass effect", )
rm(ii)

plot(MM$CENTROID_Y~MM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Temporal (4Km)")
ii <- cut(MM$grado.in.D50.4[which(MM$`EF(1).T(2).P(3)`==2)], 
          breaks = seq(min(MM$grado.in.D50.4[which(MM$`EF(1).T(2).P(3)`==2)]), 
                       max(MM$grado.in.D50.4[which(MM$`EF(1).T(2).P(3)`==2)]), len = 600), include.lowest = TRUE)
points(MM$CENTROID_Y[which(MM$`EF(1).T(2).P(3)`==2)]~MM$CENTROID_X[which(MM$`EF(1).T(2).P(3)`==2)],cex=0.4,
       col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
rm(ii)

plot(MM$CENTROID_Y~MM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Permanent (4Km)")
ii <- cut(MM$grado.in.D50.4[which(MM$`EF(1).T(2).P(3)`==3)], 
          breaks = seq(min(MM$grado.in.D50.4[which(MM$`EF(1).T(2).P(3)`==3)]), 
                       max(MM$grado.in.D50.4[which(MM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
points(MM$CENTROID_Y[which(MM$`EF(1).T(2).P(3)`==3)]~MM$CENTROID_X[which(MM$`EF(1).T(2).P(3)`==3)],cex=0.4,
       col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)

rm(ii)

plot(MM.FW$CENTROID_Y~MM.FW$CENTROID_X,col="gray",cex=0.4,pch=15,main="All-freshwater (4Km)")
ii <- cut(MM.FW$grado.in.D50.4, 
          breaks = seq(min(MM.FW$grado.in.D50.4), 
                       max(MM.FW$grado.in.D50.4), len = 600), include.lowest = TRUE)
points(MM.FW$CENTROID_Y~MM.FW$CENTROID_X,cex=0.4,
       col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)


#### Out_Degree_SourcePatchEffect ####
rm(ii)
plot(MM$CENTROID_Y~MM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Ephemeral (4Km)")
ii <- cut(MM$grado.out.D50.4[which(MM$`EF(1).T(2).P(3)`==1)], 
          breaks = seq(min(MM$grado.out.D50.4[which(MM$`EF(1).T(2).P(3)`==1)]), 
                       max(MM$grado.out.D50.4[which(MM$`EF(1).T(2).P(3)`==1)]), len = 600), include.lowest = TRUE)
points(MM$CENTROID_Y[which(MM$`EF(1).T(2).P(3)`==1)]~MM$CENTROID_X[which(MM$`EF(1).T(2).P(3)`==1)],cex=0.4,
       col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
rm(ii)
text(-8,60,"OUT DEGREE\nsource patches \nin mass effect")


plot(MM$CENTROID_Y~MM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Temporal (4Km)")
ii <- cut(MM$grado.out.D50.4[which(MM$`EF(1).T(2).P(3)`==2)], 
          breaks = seq(min(MM$grado.out.D50.4[which(MM$`EF(1).T(2).P(3)`==2)]), 
                       max(MM$grado.out.D50.4[which(MM$`EF(1).T(2).P(3)`==2)]), len = 600), include.lowest = TRUE)
points(MM$CENTROID_Y[which(MM$`EF(1).T(2).P(3)`==2)]~MM$CENTROID_X[which(MM$`EF(1).T(2).P(3)`==2)],cex=0.4,
       col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
rm(ii)

plot(MM$CENTROID_Y~MM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Permanent (4Km)")
ii <- cut(MM$grado.out.D50.4[which(MM$`EF(1).T(2).P(3)`==3)], 
          breaks = seq(min(MM$grado.out.D50.4[which(MM$`EF(1).T(2).P(3)`==3)]), 
                       max(MM$grado.out.D50.4[which(MM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
points(MM$CENTROID_Y[which(MM$`EF(1).T(2).P(3)`==3)]~MM$CENTROID_X[which(MM$`EF(1).T(2).P(3)`==3)],cex=0.4,
       col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)

rm(ii)

plot(MM.FW$CENTROID_Y~MM.FW$CENTROID_X,col="gray",cex=0.4,pch=15,main="All-freshwater (4Km)")
ii <- cut(MM.FW$grado.out.D50.4, 
          breaks = seq(min(MM.FW$grado.out.D50.4), 
                       max(MM.FW$grado.out.D50.4), len = 600), include.lowest = TRUE)
points(MM.FW$CENTROID_Y~MM.FW$CENTROID_X,cex=0.4,
       col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
#mtext("Out-degree (Source effect) \n", outer = TRUE, cex = 1.2)

##### bc
rm(ii)
plot(MM$CENTROID_Y~MM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Ephemeral (4Km)")
ii <- cut(MM$bc.D50.4[which(MM$`EF(1).T(2).P(3)`==1)], 
          breaks = seq(min(MM$bc.D50.4[which(MM$`EF(1).T(2).P(3)`==1)]), 
                       max(MM$bc.D50.4[which(MM$`EF(1).T(2).P(3)`==1)]), len = 600), include.lowest = TRUE)
points(MM$CENTROID_Y[which(MM$`EF(1).T(2).P(3)`==1)]~MM$CENTROID_X[which(MM$`EF(1).T(2).P(3)`==1)],cex=0.4,
       col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
rm(ii)
text(-10,60,"Betweenness\n steping stone \nlandscape dispersal", pos = 4)
plot(MM$CENTROID_Y~MM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Temporal (4Km)")
ii <- cut(MM$bc.D50.4[which(MM$`EF(1).T(2).P(3)`==2)], 
          breaks = seq(min(MM$bc.D50.4[which(MM$`EF(1).T(2).P(3)`==2)]), 
                       max(MM$bc.D50.4[which(MM$`EF(1).T(2).P(3)`==2)]), len = 600), include.lowest = TRUE)
points(MM$CENTROID_Y[which(MM$`EF(1).T(2).P(3)`==2)]~MM$CENTROID_X[which(MM$`EF(1).T(2).P(3)`==2)],cex=0.4,
       col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)
rm(ii)

plot(MM$CENTROID_Y~MM$CENTROID_X,col="gray",cex=0.4,pch=15,main="Permanent (4Km)")
ii <- cut(MM$bc.D50.4[which(MM$`EF(1).T(2).P(3)`==3)], 
          breaks = seq(min(MM$bc.D50.4[which(MM$`EF(1).T(2).P(3)`==3)]), 
                       max(MM$bc.D50.4[which(MM$`EF(1).T(2).P(3)`==3)]), len = 600), include.lowest = TRUE)
points(MM$CENTROID_Y[which(MM$`EF(1).T(2).P(3)`==3)]~MM$CENTROID_X[which(MM$`EF(1).T(2).P(3)`==3)],cex=0.4,
       col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)

rm(ii)

plot(MM.FW$CENTROID_Y~MM.FW$CENTROID_X,col="gray",cex=0.4,pch=15,main="All-freshwater (4Km)")
ii <- cut(MM.FW$bc.D50.4, 
          breaks = seq(min(MM.FW$bc.D50.4), 
                       max(MM.FW$bc.D50.4), len = 600), include.lowest = TRUE)
points(MM.FW$CENTROID_Y~MM.FW$CENTROID_X,cex=0.4,
       col = colorRampPalette(c("yellow", "orange", "red"))(599)[ii],pch=15)

