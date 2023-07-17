library(ggplot2)
library(gridExtra)
library(maps)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)
library(viridis)


# Nice colors? CUNILLERA_palette is what you need
source("C:/Users/David CM/Dropbox/DAVID DOC/LLAM al DIA/CUNILLERA_palette.R")

# We set the workind directory in the folder for Fifures (to avoud any conflict with real models results)
setwd("C:/Users/David CM/Dropbox/H2020/Paisaje_UK/Figuras UK_DCM")

# These lines below are copied and edited from the script "H2020_UK_lines.R" where the main plots are built in pdf
#format. 
# Previous steps to built the map from UK+Ireland

world <- ne_countries(scale = "medium", returnclass = "sf")
worldmap = map_data('world')

# Load the main Dataset
# We load "UK_space.Rdata". The original file is found in the folder UK_10x10. It is the output of the null model that 
#was used in the paper and it is the basic information to built the main plots in the paper (Figure 2,3,4) 
load("UK_space.Rdata")

# Creation of some filters to select specific "groups" of data (also used in the plotting script)
id.1<-which(M.UK[,7]==1) #Ephemerals
id_ephemeral<-c(id.1)
id.2<-which(M.UK[,7]==2) #Temporary
id_temporal<-c(id.2)
id.3<-which(M.UK[,7]==3) #Permanent
id_permanet<-c(id.3)

id_without_Ephemeral<-c(id.2,id.3)
id_without_temporal<-c(id.1,id.3)
id_without_permanent<-c(id.1,id.2)

# M.UK contains ALL the dataset, we need to only select cells with "water" on them. 
M.UK$`EF(1).T(2).P(3)` # This column contains the IDS of each "cell" that has each of the three types of habitats

# We filter the M.UK matrix removing the cells that do not belong to any of the habitats and that terefore do not have an "S"
#assigned (no water = no species)
M.UK.S<-cbind(M.UK[-which(M.UK$`EF(1).T(2).P(3)`==0),],S)

### Figure 2 Ephemeral WATERS################################################################################  

#Ephemeral waters with temporal and permanent water dispersal
# Small histogram located at the top ____
small_histo<-ggplot(data.frame(M.UK.S[id_ephemeral,]))+
             aes(S)+
             geom_histogram(col="grey70", fill="grey70", alpha=1,binwidth = 1)+
             labs(x="Local Richness", y="counts")+
             #geom_density(aes(y =..density..* (nrow(data.frame(M.UK.S[id_ephemeral,])))),alpha = 0.1, fill = "navy")+
             theme_classic()

# PLot Ephemeral with the others
A <- ggplot() + geom_polygon(data = worldmap, aes(x = long, y = lat, group = group),fill = 'gray90', color = 'black')+ 
  coord_fixed(ratio = 1.5, xlim = c(-10.5,2), ylim = c(50, 59))+
  theme_void()+
  geom_point(data = data.frame(M.UK.S[id_ephemeral,]), 
             aes(x = CENTROID_X, y =CENTROID_Y, color = (S_ephy)),shape=15, size=1.25,alpha=0.6) +
  labs(title="A)", subtitle = "Ephemeral without the others")+
  scale_color_viridis(name="S. unique Eph.",limits=c(0,35))+
  #scale_color_CUNILLERA(palette = "wildfire",discrete = F,name="S. unique Eph.")+ 
  annotation_custom(ggplotGrob(small_histo),  xmin = -1.5, xmax = 4.2,  ymin = 56.1, ymax = 59.2)


#Ephemeral waters without the others ____
# Small histogram located at the top
small_histo<-ggplot(data.frame(M.UK.S[id_ephemeral,]))+
             aes(S_ephy)+
             geom_histogram(col="grey70", fill="grey70", alpha=1, binwidth = 1)+
             labs(x="Local Richness", y="counts")+
            #geom_density(aes(y =..density..* (nrow(data.frame(M.UK.S[id_ephemeral,])))),alpha = 0.1, fill = "navy")+
             theme_classic()

# PLot Ephemeral alone
B <-ggplot()+
     geom_polygon(data = worldmap, aes(x = long, y = lat, group = group),fill = 'gray90', color = 'black') + 
     coord_fixed(ratio = 1.5, xlim = c(-10.5,2), ylim = c(50, 59)) +theme_void()+
     geom_point(data = data.frame(M.UK.S[id_ephemeral,]),aes(x = CENTROID_X, y =CENTROID_Y,
                                                             color = (S)),shape=15, size=1.25,alpha=0.6)+ 
     labs(title="B)", subtitle = "Ephemeral with Temporal and Permanent")+
     scale_color_viridis(name="S. all Eph.", limits=c(0,35))+
     #scale_color_CUNILLERA(palette = "wildfire",discrete = F, name="S. all Eph.")+  
     #scale_color_gradient(low="yellow",high="red", guide_colourbar(title ="S"))+
     annotation_custom(ggplotGrob(small_histo),  xmin = -1.5, xmax = 4.2,  ymin = 56.1, ymax = 59.2)


# Difference in richness ____
# Data calculation
M.UK.S[,13]<-0
M.UK.S[id_ephemeral,13] <- log10((M.UK.S$S[id_ephemeral]/S_ephy))
colnames(M.UK.S)[13]<-"Diff.Ephe.log.ratio"

# Small histogram located at the top
small_histo<-ggplot(data.frame(M.UK.S[id_ephemeral,]))+aes(Diff.Ephe.log.ratio)+
               geom_histogram(col="grey70", fill="grey70", alpha=1, bins = 15)+
               labs(x="Richness ratio", y="counts")+
               #geom_density(aes(y =..density..*120),alpha = 0.1, fill = "navy")+
               theme_classic()



# PLot Logratio of the difference S/Seph
C <-ggplot()+
    geom_polygon(data = worldmap, aes(x = long, y = lat, group = group),fill = 'gray90', color = 'black') + 
    coord_fixed(ratio = 1.5, xlim = c(-10.5,2), ylim = c(50, 59)) +theme_void()+
    geom_point(data = data.frame(M.UK.S[id_ephemeral,]), 
    aes(x = CENTROID_X, y =CENTROID_Y, color = (Diff.Ephe.log.ratio)),shape=15, size=1.25,alpha=0.6) + 
    labs(title="C)", subtitle = "Richness log ratio")+
    scale_color_gradient2(midpoint = 0,
                          low=viridis(1,direction = -1),
                          mid="white",
                          high=viridis(1,direction = 1),
                          guide_colourbar(title ="SLR"))+
    annotation_custom(ggplotGrob(small_histo), xmin = -1.5, xmax = 4.2,  ymin = 56.1, ymax = 59.2)


png(filename ="Figure2_Ephemeral.png", 
    width = 900*4.5, height = 600*4.5, 
    units = "px",res = 300)
grid.arrange(A,B,C, ncol=3, nrow=1)
dev.off()
  

### Figure 3 Temporal WATERS################################################################################  

#Temporal waters with temporal and permanent water dispersal
# Small histogram located at the top ____
small_histo<-ggplot(data.frame(M.UK.S[id_temporal,]))+
  aes(S)+
  geom_histogram(col="grey70", fill="grey70", alpha=1,binwidth = 1)+
  labs(x="Local Richness", y="counts")+
  #geom_density(aes(y =..density..* (nrow(data.frame(M.UK.S[id_temporal,])))),alpha = 0.1, fill = "navy")+
  theme_classic()

# PLot Ephemeral with the others
A <- ggplot() + geom_polygon(data = worldmap, aes(x = long, y = lat, group = group),fill = 'gray90', color = 'black')+ 
  coord_fixed(ratio = 1.5, xlim = c(-10.5,2), ylim = c(50, 59))+
  theme_void()+
  geom_point(data = data.frame(M.UK.S[id_temporal,]), 
             aes(x = CENTROID_X, y =CENTROID_Y, color = (S_temporal)),shape=15, size=1.25,alpha=0.6) +
  labs(title="A)", subtitle = "Temporal without the others")+
  scale_color_viridis(name="S. unique Temp.",limits=c(0,55))+
  #scale_color_CUNILLERA(palette = "wildfire",discrete = F,name="S. unique Temp.")+ 
  annotation_custom(ggplotGrob(small_histo),  xmin = -1.5, xmax = 4.2,  ymin = 56.1, ymax = 59.2)

  #Ephemeral waters without the others ____
  # Small histogram located at the top
small_histo<-ggplot(data.frame(M.UK.S[id_temporal,]))+
  aes(S_temporal)+
  geom_histogram(col="grey70", fill="grey70", alpha=1, binwidth = 1)+
  labs(x="Local Richness", y="counts")+
  #geom_density(aes(y =..density..* (nrow(data.frame(M.UK.S[id_temporal,])))),alpha = 0.1, fill = "navy")+
  theme_classic()

# PLot Ephemeral alone
B <-ggplot()+
  geom_polygon(data = worldmap, aes(x = long, y = lat, group = group),fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.5, xlim = c(-10.5,2), ylim = c(50, 59)) +theme_void()+
  geom_point(data = data.frame(M.UK.S[id_temporal,]),aes(x = CENTROID_X, y =CENTROID_Y,
                                                          color = (S)),shape=15, size=1.25,alpha=0.6)+ 
  labs(title="B)", subtitle = "Temporal with Ephemeral and Permanent")+
  scale_color_viridis(name="S. all Temp.",limits=c(0,55))+
  #scale_color_CUNILLERA(palette = "wildfire",discrete = F, name="S. all Temp.")+  
  #scale_color_gradient(low="yellow",high="red", guide_colourbar(title ="S"))+
  annotation_custom(ggplotGrob(small_histo),  xmin = -1.5, xmax = 4.2,  ymin = 56.1, ymax = 59.2)

# Difference in richness ____
# Data calculation
M.UK.S[,13]<-0
M.UK.S[id_temporal,13] <- log10((M.UK.S$S[id_temporal]/S_temporal))
colnames(M.UK.S)[13]<-"Diff.Temp.log.ratio"

# Small histogram located at the top
small_histo<-ggplot(data.frame(M.UK.S[id_temporal,]))+aes(Diff.Temp.log.ratio)+
  geom_histogram(col="grey70", fill="grey70", alpha=1, bins = 15)+
  labs(x="Richness ratio", y="counts")+
  #geom_density(aes(y =..density..*200),alpha = 0.1, fill = "navy")+
  theme_classic()

# PLot Logratio of the difference S/Seph
C <-ggplot()+
  geom_polygon(data = worldmap, aes(x = long, y = lat, group = group),fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.5, xlim = c(-10.5,2), ylim = c(50, 59)) +theme_void()+
  geom_point(data = data.frame(M.UK.S[id_temporal,]), 
             aes(x = CENTROID_X, y =CENTROID_Y, color = (Diff.Temp.log.ratio)),shape=15, size=1.25,alpha=0.6) + 
  labs(title="C)", subtitle = "Richness log ratio")+
  scale_color_gradient2(midpoint = 0,
                        low=viridis(1,direction = -1),
                        mid="white",
                        high=viridis(1,direction = 1),
                        guide_colourbar(title ="SLR"))+
  annotation_custom(ggplotGrob(small_histo), xmin = -1.5, xmax = 4.2,  ymin = 56.1, ymax = 59.2)



png(filename ="Figure3_Temporal.png", 
    width = 900*4.5, height = 600*4.5, 
    units = "px",res = 300)
grid.arrange(A,B,C, ncol=3, nrow=1)
dev.off()


### Figure 4 Permanent WATERS################################################################################  

#Ephemeral waters with temporal and permanent water dispersal
# Small histogram located at the top ____
small_histo<-ggplot(data.frame(M.UK.S[id_permanet,]))+
  aes(S)+
  geom_histogram(col="grey70", fill="grey70", alpha=1,binwidth = 1)+
  labs(x="Local Richness", y="counts")+
  #geom_density(aes(y =..density..* (nrow(data.frame(M.UK.S[id_permanet,])))),alpha = 0.1, fill = "navy")+
  theme_classic()

# PLot Ephemeral with the others
A <-ggplot() + geom_polygon(data = worldmap, aes(x = long, y = lat, group = group),fill = 'gray90', color = 'black')+ 
  coord_fixed(ratio = 1.5, xlim = c(-10.5,2), ylim = c(50, 59))+
  theme_void()+
  geom_point(data = data.frame(M.UK.S[id_permanet,]), 
             aes(x = CENTROID_X, y =CENTROID_Y, color = (S_permanent)),shape=15, size=1.25,alpha=0.6) +
  labs(title="A)", subtitle = "Permanent without the others")+
  scale_color_viridis(name="S. unique Perm.",limits=c(0,80))+
  #scale_color_CUNILLERA(palette = "wildfire",discrete = F,name="S. unique Perm.")+ 
  annotation_custom(ggplotGrob(small_histo),  xmin = -1.5, xmax = 4.2,  ymin = 56.1, ymax = 59.2) 

  #Ephemeral waters without the others ____
  # Small histogram located at the top
  small_histo<-ggplot(data.frame(M.UK.S[id_permanet,]))+
    aes(S_permanent)+
    geom_histogram(col="grey70", fill="grey70", alpha=1,binwidth = 1)+
    labs(x="Local Richness", y="counts")+
    #geom_density(aes(y =..density..* (nrow(data.frame(M.UK.S[id_permanet,])))),alpha = 0.1, fill = "navy")+
    theme_classic()


# PLot Ephemeral alone
B <-ggplot()+
  geom_polygon(data = worldmap, aes(x = long, y = lat, group = group),fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.5, xlim = c(-10.5,2), ylim = c(50, 59)) +theme_void()+
  geom_point(data = data.frame(M.UK.S[id_permanet,]),aes(x = CENTROID_X, y =CENTROID_Y,
                                                          color = (S)),shape=15, size=1.25,alpha=0.6)+ 
  labs(title="B)", subtitle = "Permanent with Ephemeral and Temporal")+
  scale_color_viridis(name="S. all Perm.",limits=c(0,80))+
  #scale_color_CUNILLERA(palette = "wildfire",discrete = F, name="S. all Perm.")+  
  #scale_color_gradient(low="yellow",high="red", guide_colourbar(title ="S"))+
  annotation_custom(ggplotGrob(small_histo),  xmin = -1.5, xmax = 4.2,  ymin = 56.1, ymax = 59.2)


# Difference in richness ____
# Data calculation
M.UK.S[,13]<-0
M.UK.S[id_permanet,13] <- log10((M.UK.S$S[id_permanet]/S_permanent))
colnames(M.UK.S)[13]<-"Diff.Perm.log.ratio"

# Small histogram located at the top
small_histo<-ggplot(data.frame(M.UK.S[id_permanet,]))+aes(Diff.Perm.log.ratio)+
  geom_histogram(col="grey70", fill="grey70", alpha=1,bins = 15)+
  labs(x="Richness ratio", y="counts")+
  #geom_density(aes(y =..density..*90),alpha = 0.1, fill = "navy")+
  theme_classic()

# PLot Logratio of the difference S/Seph
C <-ggplot()+
  geom_polygon(data = worldmap, aes(x = long, y = lat, group = group),fill = 'gray90', color = 'black') + 
  coord_fixed(ratio = 1.5, xlim = c(-10.5,2), ylim = c(50, 59)) +theme_void()+
  geom_point(data = data.frame(M.UK.S[id_permanet,]), 
             aes(x = CENTROID_X, y =CENTROID_Y, color = (Diff.Perm.log.ratio)),shape=15, size=1.25,alpha=0.6) + 
  labs(title="C)", subtitle = "Richness log ratio")+
  scale_color_gradient2(midpoint = 0,
                        low=viridis(1,direction = -1),
                        mid="white",
                        high=viridis(1,direction = 1),
                        guide_colourbar(title ="SLR"))+
  annotation_custom(ggplotGrob(small_histo), xmin = -1.5, xmax = 4.2,  ymin = 56.1, ymax = 59.2)



png(filename ="Figure4_Permanent.png", 
    width = 900*4.5, height = 600*4.5, 
    units = "px",res = 300)
grid.arrange(A,B,C, ncol=3, nrow=1)
dev.off()


# Scatterplot S/Sephy ____
D_eph <- ggplot(data.frame(M.UK.S[id_ephemeral,]))+ aes(x=S_ephy, y=S)+
  geom_point(aes(fill=S),shape=22,size=3, col="grey20", alpha=0.2)+
  geom_smooth(method="lm", color="black", size=1.5)+
  geom_abline(intercept=0, slope=1, color="red", size=1, linetype=2)+
  scale_fill_viridis()+
  #scale_fill_CUNILLERA(palette = "wildfire",discrete = F)+
  labs(title="A)", subtitle = "Richness relationship (S. all Eph./S. unique Eph.)",
       y="S. all Eph.",x="S. unique Eph.")+
  geom_text(aes(y=40, x=10), label=paste("Slope=", 
                                        as.character(round(lm(M.UK.S[id_ephemeral,]$S~S_ephy)$coefficients[2],2)),
                                        sep=" "))+
  scale_y_continuous(limits = c(0,40))+
  scale_x_continuous(limits = c(0,40))+
  theme_classic()+
  theme(plot.margin=margin(t = 2, r = 2.5, b = 0.5, l = 2.5, unit = "cm"),
        legend.position = "none")

# Scatterplot S/Sephy ____
D_temp <- ggplot(data.frame(M.UK.S[id_temporal,]))+ aes(x=S_temporal, y=S)+
  geom_point(aes(fill=S),shape=22,size=3, col="grey20", alpha=0.2)+
  geom_smooth(method="lm", color="black", size=1.5)+
  geom_abline(intercept=0, slope=1, color="red", size=1, linetype=2)+
  scale_fill_viridis()+
  #scale_fill_CUNILLERA(palette = "wildfire",discrete = F)+
  labs(title="B)", subtitle = "Richness relationship (S. all Temp./S. unique Temp.)",
       y="S. all Temp.",x="S. unique Temp.")+
  geom_text(aes(y=60, x=10), label=paste("Slope=", 
                                        as.character(round(lm(M.UK.S[id_temporal,]$S~S_temporal)$coefficients[2],2)),
                                        sep=" "))+
  scale_y_continuous(limits = c(0,60))+
  scale_x_continuous(limits = c(0,60))+
  theme_classic()+
  theme(plot.margin=margin(t = 2, r = 2.5, b = 0.5, l = 2.5, unit = "cm"),
        legend.position = "none")

# Scatterplot S/Sperm ____
D_perm <- ggplot(data.frame(M.UK.S[id_permanet,]))+ aes(x=S_permanent, y=S)+
  geom_point(aes(fill=S),shape=22,size=3, col="grey20", alpha=0.2)+
  geom_smooth(method="lm", color="black", size=1.5)+
  geom_abline(intercept=0, slope=1, color="red", size=1, linetype=2)+
  scale_fill_viridis()+
  #scale_fill_CUNILLERA(palette = "wildfire",discrete = F)+
  labs(title="C)", subtitle = "Richness relationship (S. all Perm./S. unique Perm.)",
       y="S. all Perm.",x="S. unique Perm.")+
  geom_text(aes(y=100, x=10), label=paste("Slope=", 
                                        as.character(round(lm(M.UK.S[id_permanet,]$S~S_permanent)$coefficients[2],2)),
                                        sep=" "))+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  theme_classic()+
  theme(plot.margin=margin(t = 2, r = 2.5, b = 0.5, l = 2.5, unit = "cm"),
        legend.position = "none")

png(filename ="Figure5_Rel.png", 
    width = 1000*7, height = 300*4.5, 
    units = "px",res = 300)
grid.arrange(D_eph,D_temp,D_perm, ncol=3, nrow=1)
dev.off()

table_result <- rbind(
as.matrix(round(summary(lm(M.UK.S[id_ephemeral,]$S~S_ephy))[[4]],3)),
as.matrix(round(summary(lm(M.UK.S[id_temporal,]$S~S_temporal))[[4]],3)),
as.matrix(round(summary(lm(M.UK.S[id_permanet,]$S~S_permanent))[[4]],3))
)

rownames(table_result) <- NULL   
write.table(table_result,"table_result.txt",sep = ",")



# Comparison between lat and long 

Plot_geogr_relationships <- gridExtra::grid.arrange(
# Ephemeral
gridExtra::arrangeGrob(
data.frame(lat=data.frame(M.UK.S[id_ephemeral,])[,5],long=data.frame(M.UK.S[id_ephemeral,])[,4],
           S=data.frame(M.UK.S[id_ephemeral,])[,12],Type="ephem_All") %>% 
bind_rows(data.frame(lat=data.frame(M.UK.S[id_ephemeral,])[,5],long=data.frame(M.UK.S[id_ephemeral,])[,4],
          S=S_ephy,Type="ephem_Only"))%>% 
  ggplot()+
  geom_jitter(aes(x=lat, y=S, fill=Type), shape=21, alpha=0.1, colour="black")+
  geom_density2d(aes(x=lat, y=S, colour=Type),alpha=.5, size=1.5, linetype=2)+ 
  geom_smooth(aes(x=lat, y=S, colour=Type), method = "lm")+
  scale_fill_viridis(discrete = T)+
  scale_colour_viridis(discrete = T, label=c("S. all Eph.", "S. unique Eph."))+
  guides(fill="none")+labs(x="Latitude", y="Species richness")+
  geom_text(aes(x=59,y=33), colour=viridis(1,direction = -1),
            label=paste("Slope=",
            round(as.numeric(lm(S_ephy~
            data.frame(M.UK.S[id_ephemeral,])[,5])[[1]][2]),2),
            sep=" "))+
  geom_text(aes(x=59,y=35), colour=viridis(1,direction = 1),
            label=paste("Slope=",
            round(as.numeric(lm(data.frame(M.UK.S[id_ephemeral,])[,12]~
            data.frame(M.UK.S[id_ephemeral,])[,5])[[1]][2]),2),
            sep=" "))+
  theme_classic(),

data.frame(lat=data.frame(M.UK.S[id_ephemeral,])[,5],long=data.frame(M.UK.S[id_ephemeral,])[,4],
           S=data.frame(M.UK.S[id_ephemeral,])[,12],Type="ephem_All") %>% 
  bind_rows(data.frame(lat=data.frame(M.UK.S[id_ephemeral,])[,5],long=data.frame(M.UK.S[id_ephemeral,])[,4],
                       S=S_ephy,Type="ephem_Only"))%>% 
  ggplot()+
  geom_jitter(aes(x=long, y=S, fill=Type), shape=21, alpha=0.1, colour="black")+
  geom_density2d(aes(x=long, y=S, colour=Type),alpha=.3, size=1.5, linetype=2)+ 
  geom_smooth(aes(x=long, y=S, colour=Type), method = "lm", size=2, linetype=1)+
  scale_fill_viridis(discrete = T)+
  scale_colour_viridis(discrete = T, label=c("S. all Eph.", "S. unique Eph."))+
  guides(fill="none")+labs(x="Longitude", y="Species richness")+
  geom_text(aes(x=-1,y=33), colour=viridis(1,direction = -1),
            label=paste("Slope=",
            round(as.numeric(lm(S_ephy~
            data.frame(M.UK.S[id_ephemeral,])[,4])[[1]][2]),2),
            sep=" "))+
  geom_text(aes(x=-1,y=35), colour=viridis(1,direction = 1),
            label=paste("Slope=",
            round(as.numeric(lm(data.frame(M.UK.S[id_ephemeral,])[,12]~
            data.frame(M.UK.S[id_ephemeral,])[,4])[[1]][2]),2),
            sep=" "))+
  theme_classic(),
ncol=2, top="Ephemeral"),

# Temporal 
gridExtra::arrangeGrob(
  data.frame(lat=data.frame(M.UK.S[id_temporal,])[,5],long=data.frame(M.UK.S[id_temporal,])[,4],
             S=data.frame(M.UK.S[id_temporal,])[,12],Type="tempo_All") %>% 
    bind_rows(data.frame(lat=data.frame(M.UK.S[id_temporal,])[,5],long=data.frame(M.UK.S[id_temporal,])[,4],
                         S=S_temporal,Type="tempo_Only"))%>% 
    ggplot()+
    geom_jitter(aes(x=lat, y=S, fill=Type), shape=21, alpha=0.1, colour="black")+
    geom_density2d(aes(x=lat, y=S, colour=Type),alpha=.5, size=1.5, linetype=2)+ 
    geom_smooth(aes(x=lat, y=S, colour=Type), method = "lm")+
    scale_fill_viridis(discrete = T)+
    scale_colour_viridis(discrete = T, label=c("S. all Temp.", "S. unique Temp."))+
    guides(fill="none")+labs(x="Latitude", y="Species richness")+
    geom_text(aes(x=59,y=50), colour=viridis(1,direction = -1),
              label=paste("Slope=",
              round(as.numeric(lm(S_temporal~
              data.frame(M.UK.S[id_temporal,])[,5])[[1]][2]),2),
              sep=" "))+
    geom_text(aes(x=59,y=53), colour=viridis(1,direction = 1),
              label=paste("Slope=",
              round(as.numeric(lm(data.frame(M.UK.S[id_temporal,])[,12]~
              data.frame(M.UK.S[id_temporal,])[,5])[[1]][2]),2),
              sep=" "))+
    theme_classic(),
  
  data.frame(lat=data.frame(M.UK.S[id_temporal,])[,5],long=data.frame(M.UK.S[id_temporal,])[,4],
             S=data.frame(M.UK.S[id_temporal,])[,12],Type="tempo_All") %>% 
    bind_rows(data.frame(lat=data.frame(M.UK.S[id_temporal,])[,5],long=data.frame(M.UK.S[id_temporal,])[,4],
                         S=S_temporal,Type="tempo_Only"))%>% 
    ggplot()+
    geom_jitter(aes(x=long, y=S, fill=Type), shape=21, alpha=0.1, colour="black")+
    geom_density2d(aes(x=long, y=S, colour=Type),alpha=.3, size=1.5, linetype=2)+ 
    geom_smooth(aes(x=long, y=S, colour=Type), method = "lm", size=2, linetype=1)+
    scale_fill_viridis(discrete = T)+
    scale_colour_viridis(discrete = T, label=c("S. all Temp.", "S. unique Temp."))+
    guides(fill="none")+labs(x="Longitude", y="Species richness")+
    geom_text(aes(x=0,y=50), colour=viridis(1,direction = -1),
              label=paste("Slope=",
              round(as.numeric(lm(S_temporal~
              data.frame(M.UK.S[id_temporal,])[,4])[[1]][2]),2),
              sep=" "))+
    geom_text(aes(x=0,y=53), colour=viridis(1,direction = 1),
              label=paste("Slope=",
              round(as.numeric(lm(data.frame(M.UK.S[id_temporal,])[,12]~
              data.frame(M.UK.S[id_temporal,])[,4])[[1]][2]),2),
              sep=" "))+
    theme_classic(),
  ncol=2, top="Temporal"),

# Permanent
gridExtra::arrangeGrob(
  data.frame(lat=data.frame(M.UK.S[id_permanet,])[,5],long=data.frame(M.UK.S[id_permanet,])[,4],
             S=data.frame(M.UK.S[id_permanet,])[,12],Type="perma_All") %>% 
    bind_rows(data.frame(lat=data.frame(M.UK.S[id_permanet,])[,5],long=data.frame(M.UK.S[id_permanet,])[,4],
                         S=S_permanent,Type="perma_Only"))%>% 
    ggplot()+
    geom_jitter(aes(x=lat, y=S, fill=Type), shape=21, alpha=0.1, colour="black")+
    geom_density2d(aes(x=lat, y=S, colour=Type),alpha=.5, size=1.5, linetype=2)+ 
    geom_smooth(aes(x=lat, y=S, colour=Type), method = "lm")+
    scale_fill_viridis(discrete = T)+
    scale_colour_viridis(discrete = T, label=c("S. all Eph.", "S. unique Eph."))+
    guides(fill="none")+labs(x="Latitude", y="Species richness")+
    geom_text(aes(x=59,y=70), colour=viridis(1,direction = -1),
              label=paste("Slope=",
              round(as.numeric(lm(S_permanent~
              data.frame(M.UK.S[id_permanet,])[,5])[[1]][2]),2),
              sep=" "))+
    geom_text(aes(x=59,y=75), colour=viridis(1,direction = 1),
              label=paste("Slope=",
              round(as.numeric(lm(data.frame(M.UK.S[id_permanet,])[,12]~
              data.frame(M.UK.S[id_permanet,])[,5])[[1]][2]),2),
              sep=" "))+
    theme_classic(),
  
  data.frame(lat=data.frame(M.UK.S[id_permanet,])[,5],long=data.frame(M.UK.S[id_permanet,])[,4],
             S=data.frame(M.UK.S[id_permanet,])[,12],Type="perma_All") %>% 
    bind_rows(data.frame(lat=data.frame(M.UK.S[id_permanet,])[,5],long=data.frame(M.UK.S[id_permanet,])[,4],
                         S=S_permanent,Type="perma_Only"))%>% 
    ggplot()+
    geom_jitter(aes(x=long, y=S, fill=Type), shape=21, alpha=0.1, colour="black")+
    geom_density2d(aes(x=long, y=S, colour=Type),alpha=.3, size=1.5, linetype=2)+ 
    geom_smooth(aes(x=long, y=S, colour=Type), method = "lm", size=2, linetype=1)+
    scale_fill_viridis(discrete = T)+
    scale_colour_viridis(discrete = T, label=c("S. all Perm.", "S. unique Perm."))+
    guides(fill="none")+labs(x="Longitude", y="Species richness")+
    geom_text(aes(x=0,y=70), colour=viridis(1,direction = -1),
    label=paste("Slope=",
    round(as.numeric(lm(S_permanent~
                        data.frame(M.UK.S[id_permanet,])[,4])[[1]][2]),2),
    sep=" "))+
    geom_text(aes(x=0,y=75), colour=viridis(1,direction = 1),
    label=paste("Slope=",
    round(as.numeric(lm(data.frame(M.UK.S[id_permanet,])[,12]~
                        data.frame(M.UK.S[id_permanet,])[,4])[[1]][2]),2),
    sep=" "))+
    theme_classic(),
  ncol=2, top="Permanent"), 
nrow=3)

png(filename ="FigureSup_Rel.png", 
    width = 600*7, height = 600*7, 
    units = "px",res = 300)
grid.arrange(Plot_geogr_relationships)
dev.off()


table_result <- rbind(
#Long
as.matrix(round(summary(lm(S_ephy~data.frame(M.UK.S[id_ephemeral,])[,4]))[[4]],3)),
as.matrix(round(summary(lm(data.frame(M.UK.S[id_ephemeral,])[,12]~data.frame(M.UK.S[id_ephemeral,])[,4]))[[4]],3)),
#Lat
as.matrix(round(summary(lm(S_ephy~data.frame(M.UK.S[id_ephemeral,])[,5]))[[4]],3)),
as.matrix(round(summary(lm(data.frame(M.UK.S[id_ephemeral,])[,12]~data.frame(M.UK.S[id_ephemeral,])[,5]))[[4]],3)),

#Long
as.matrix(round(summary(lm(S_temporal~data.frame(M.UK.S[id_temporal,])[,4]))[[4]],3)),
as.matrix(round(summary(lm(data.frame(M.UK.S[id_temporal,])[,12]~data.frame(M.UK.S[id_temporal,])[,4]))[[4]],3)),
#Lat
as.matrix(round(summary(lm(S_temporal~data.frame(M.UK.S[id_temporal,])[,5]))[[4]],3)),
as.matrix(round(summary(lm(data.frame(M.UK.S[id_temporal,])[,12]~data.frame(M.UK.S[id_temporal,])[,5]))[[4]],3)),

#Long
as.matrix(round(summary(lm(S_permanent~data.frame(M.UK.S[id_permanet,])[,4]))[[4]],3)),
as.matrix(round(summary(lm(data.frame(M.UK.S[id_permanet,])[,12]~data.frame(M.UK.S[id_permanet,])[,4]))[[4]],3)),
#Lat
as.matrix(round(summary(lm(S_permanent~data.frame(M.UK.S[id_permanet,])[,5]))[[4]],3)),
as.matrix(round(summary(lm(data.frame(M.UK.S[id_permanet,])[,12]~data.frame(M.UK.S[id_permanet,])[,5]))[[4]],3))
)
rownames(table_result) <- NULL   
write.table(table_result,"table_result.txt",sep = ",")

# Figure 1 Fitness curves ####
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
Perf.3<- exp(p[1]+p[2]*Species+p[3]*Species^2)/(1+exp(p[1]+p[2]*Species+p[3]*Species^2))# performance species pool 3

#Plot

png(filename ="Figure1_Performance.png", 
    width = 700*3, height = 500*3, 
    units = "px",res = 300)
data.frame(Species,Perf.1,Perf.2,Perf.3)%>%
  ggplot()+
  geom_line(aes(x=Species, y=Perf.1), size=3,colour=cividis(3)[1],alpha=1)+
  geom_line(aes(x=Species, y=Perf.2), size=3,colour=cividis(3)[2],alpha=1)+
  geom_line(aes(x=Species, y=Perf.3), size=3,colour=cividis(3)[3],alpha=1)+
  scale_y_continuous(expand = c(0.02,0.02),breaks = c(0,0.25,0.5,0.75,1),labels = c(0,0.25,0.5,0.75,1))+
  scale_x_continuous(expand = c(0.01,0.01),breaks = c(1,25,50,75,100,125,150,175,200),labels = c(1,25,50,75,100,125,150,175,200))+
  labs(title = "Species habitat performance", y="Affinity")+
  annotate(x=25,y=1, label="Ephemeral",geom="label", fill=cividis(3)[2], colour="white",alpha=1, size=5)+
  annotate(x=100,y=0.9,label="Temporal", geom="label", fill=cividis(3)[3], colour="white",alpha=1, size=5)+
  annotate(x=175,y=1,label="Permanent",geom="label", fill=cividis(3)[1], colour="white",alpha=1, size=5)+
  theme_classic()
dev.off()
