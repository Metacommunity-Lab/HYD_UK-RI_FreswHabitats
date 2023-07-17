# in a logistc f(x)/1+ f(x) when f(x)= a + bx p(y)=0.5 is observed in x) -a/b
# Larger b implies steeper transitions
# when a<0 and b<0 derivate is positive
# when a>0 and b<0 derivate is negative

###############
# ESTIMATES FILTER BY ENVIRONMENT
# These lines generates vector of species perfomances on different enviornments
# It is assumed that there is a species pool with S species. Some species do better in one environemnt (e.g. ephimeral),...
# and other do better in other environments (e.g. temporal or permanent)

S<-200       # Number of species
b1=0.05; a1<- -S*0.65*b1
Species<-1:S
Perf.1<- exp(a1+b1*Species)/(1+exp(a1+b1*Species)) # Performance of species pool 1


b2= -0.05; a2<- -S*0.35*b2
Perf.2<- exp(a2+b2*Species)/(1+exp(a2+b2*Species))# performance of species pool 2

Perf.matrix<-cbind(Perf.1, Perf.2)
colnames(Perf.matrix)<-c(20,1)

y<-c(rep(0,19),1,1,rep(0,19),rep(1,58),rep(0,4),rep(1,58),rep(0,19),1,1,rep(0,19))
y[9]<-1;y[187]<-1 #;y[99:101]<-0;
summary(glm(y~Species+I(Species^2), family = binomial))
p<-coefficients(glm(y~Species+I(Species^2), family = binomial))
#b3.1= 55; b3.2= -0.28; a3<- -0.002#-S*0.85*b2
#Perf.2<- exp(a3+b3.1*Species+b3.2*Species^2)/(1+exp(a3+b3.1*Species+b3.2*Species^2))
Perf.3<- exp(p[1]+p[2]*Species+p[3]*Species^2)/(1+exp(p[1]+p[2]*Species+p[3]*Species^2))# performance species pool 3

par(mfrow=c(1,1), bty="l", mar=c(4,4,.5,.5))
plot(Perf.1~Species, ylab = "Performance", xlab = "Species Id", lwd=3, col="navy", type="l", ylim=c(0,1))
points(Perf.3~Species, lwd=3, col="olivedrab"  , type="l", lty=1)
points(Perf.2~Species, lwd=3, col="tomato2", type="l", lty=1)
text(10,.9,"Ephimeral", col="tomato2")
text(100,.9,"Temporal", col="olivedrab")
text(190,.9,"Permanent", col="navy")

Perf.matrix<-cbind(Perf.1, Perf.2, Perf.3)       # build a matrix with columns "performance of species pool" in each one of the different enviornments
colnames(Perf.matrix)<-c(1:3)

# NOTE: species performance in different environemnet could be estimated by fitting logistc model to species occurrence among environments
########################################################################
#ESTIMATES FILTER(SUITABILITY) MATRIX
#M.UK.cent.D50.4 is a matrix of UK freshwater cover
id.0<-which(M.UK.cent.D50.4$`EF(1).T(2).P(3)`==0)
M.UK<-M.UK.cent.D50.4[-id.0,]

Filter<-NULL
for(i in 1:nrow(M.UK)){
  Filter<-cbind(Filter, Perf.matrix[,which(colnames(Perf.matrix)==M.UK[i,7])])
}
UK<-M.UK.cent.D50.4 
Filter.UK<-NULL
for(i in 1:nrow(UK)){
  Filter.UK<-cbind(Filter.UK, Perf.matrix[,which(colnames(Perf.matrix)==UK[i,7])])
}
dim(Filter.UK)
dim(UK)

