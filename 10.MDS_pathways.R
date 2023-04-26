#Scripts agafats dels que vaig fer per al paper de PeerJ del 2018
#SCRIPTS per fer el MDS i els simper de les prediccions funcionals
library(vegan)
library (MASS)
library(scales)
library(ggplot2)
#remotes::install_github("AckerDWM/gg3D")
#library(spdep)

#scales es per fer els colors transparents amb la funcio alpha de forma m?s facil
#spdep es per rotar els punts del COI

setwd('E:/styela/Molecular/Microbioma/microbioma/')

#### data preparation ####

mx1 <-read.csv("E:/styela/Molecular/Microbioma/microbioma/Coses Xavier/tax4fun2/Elena_19410_prediction/Functional_prediction.csv",stringsAsFactors = F)

mx1

# #si ens volem quedar nomes amb els 100 primers ESVs que simper ha detectat entre juvenil i adult
# juvadult<-read.csv("taula_simper_KO_juv_adult.csv",stringsAsFactors = F)
# mx1<-mx1[juvadult$row_name,]



#eliminem les columnes no numeriques
mx1<-mx1[,2:105]



#transpose the data, in order to have the samples as observations (rows) and the sequences as variables (columns)
samples <- t(mx1)



##quantitative
#hem de passar a relative frequencies
# mx2<-as.matrix(mx1)
# mx3<-prop.table(mx2,2)
# sample<-t(mx3)

#read codes
codes<-read.table("metadata_filtered.txt", sep="\t", stringsAsFactors = F)



###MDS


#segons si es vol fer sqrt o sqrt(sqrt) triar un dels seguents:
# sampless<-sqrt(samples)
# sampless<-sqrt(sqrt(samples))

#metanmds <- metaMDS (sampless, trymax=500, distance="jaccard",binary=T,autotransform=F)
#si volem tenir m?s dimensions
#metanmds <- metaMDS (sampless,k=5,trymax=500, distance="bray",autotransform=F)
#reserva<-metanmds

metanmds <- metaMDS(samples,k=3,trymax=1000, distance="bray",autotransform=F)
#ordiplot(metanmds,dis="sites")
metanmds$stress



#el fem rotar mig radian perque quedi mes maco
#NOTA: amb sqrt no cal rotar, ho he fet nomes per sqrt(sqrt))
# metanmds$points<-Rotation(metanmds$points,0.5)
# ordiplot(metanmds,dis="sites")



#per guardar el plot
#jpeg(file="MDS_path_19410.jpg",res=600,quality=100,width=4000,height=4000)
#pdf(file="MDS_ellipse_jaccard_95percent_stress0.0803.pdf")
#pdf(file="Fig 4.pdf")

#plot(metanmds$points,type="n", xlab="", ylab="", xlim=c(-1,0.5),ylim=c(-0.5,0.5), main=paste("12S","\n","stress =",round(metanmds$stress,digits=3)))

#plot(metanmds$points,type="n", xlab="", ylab="",xaxt="n",yaxt="n", main=paste("COI ","\n","stress =",round(metanmds$stress,digits=3)))

#per plotejar diverses dimensions
#metanmds$points<-reserva$points[,4:5] #posar les columnes de les dimensions que el vulguin plotejar


#plot(metanmds$points,type="n", xlab="", ylab="",xaxt="n",yaxt="n",main=paste("Pathways ","\n","stress =",round(metanmds$stress,digits=3)))

MDS<-as.data.frame(metanmds$points)
MDS['POP'] <- codes$POP
MDS['STAGE'] <- codes$STAGE
MDS['TISSUE'] <- codes$TISSUE


pdf("MDS_functions_1vs2.pdf", width = 7, height = 6)
ggplot(data = MDS, aes(x = MDS1, y = MDS2))+
  geom_point(data = MDS, size = 4, aes(color=STAGE, shape=POP, fill=TISSUE, alpha=1), stroke = 1.5)+
  geom_hline(yintercept = 0, lty = "dotted", color="grey", cex=1) +    geom_vline(xintercept = 0, lty = "dotted", color="grey", cex=1) +
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_fill_manual(values= c("#58D68D", "#EC7063","#F4D03F","#85C1E9"), labels=c("GILL","DIGESTIVE","TUNIC","WATER"))+
  scale_shape_manual(values=c(25,23,21)) +
  scale_color_manual(values=c("purple4","#DB7093","steelblue")) +
  theme_classic()+
  labs(title = paste("Stress  ", 0.0644, sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))   
dev.off()

pdf("MDS_functions_3vs2.pdf", width = 7, height = 6)
ggplot(data = MDS, aes(x = MDS3, y = MDS2))+
  geom_point(data = MDS, size = 4, aes(color=STAGE, shape=POP, fill=TISSUE, alpha=1), stroke = 1.5)+
  geom_hline(yintercept = 0, lty = "dotted", color="grey", cex=1) +    geom_vline(xintercept = 0, lty = "dotted", color="grey", cex=1) +
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_fill_manual(values= c("#58D68D", "#EC7063","#F4D03F","#85C1E9"), labels=c("GILL","DIGESTIVE","TUNIC","WATER"))+
  scale_shape_manual(values=c(25,23,21)) +
  scale_color_manual(values=c("purple4","#DB7093","steelblue")) +
  theme_classic()+
  labs(title = paste("Stress  ", 0.0644, sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))   
dev.off()


mx2 <-read.csv("E:/styela/Molecular/Microbioma/microbioma/Coses Xavier/tax4fun2/Elena_19410_prediction/Pathway_prediction.csv",stringsAsFactors = F)

mx2

#eliminem les columnes no numeriques
mx2<-mx2[,2:105]

#transpose the data, in order to have the samples as observations (rows) and the sequences as variables (columns)
samplesp <- t(mx2)


metanmdsp <- metaMDS(samplesp,k=3,trymax=1000, distance="bray",autotransform=F)
metanmdsp$stress


MDSp<-as.data.frame(metanmdsp$points)
MDSp['POP'] <- codes$POP
MDSp['STAGE'] <- codes$STAGE
MDSp['TISSUE'] <- codes$TISSUE


pdf("MDS_pathways_1vs2.pdf", width = 7, height = 6)
ggplot(data = MDSp, aes(x = MDS1, y = MDS2))+
  geom_point(data = MDSp, size = 4, aes(color=STAGE, shape=POP, fill=TISSUE, alpha=1), stroke = 1.5)+
  geom_hline(yintercept = 0, lty = "dotted", color="grey", cex=1) +    geom_vline(xintercept = 0, lty = "dotted", color="grey", cex=1) +
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_fill_manual(values= c("#58D68D", "#EC7063","#F4D03F","#85C1E9"), labels=c("GILL","DIGESTIVE","TUNIC","WATER"))+
  scale_shape_manual(values=c(25,23,21)) +
  scale_color_manual(values=c("purple4","#DB7093","steelblue")) +
  theme_classic()+
  labs(title = paste("Stress  ", 0.0366, sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))   
dev.off()

pdf("MDS_pathways_3vs2.pdf", width = 7, height = 6)
ggplot(data = MDSp, aes(x = MDS3, y = MDS2))+
  geom_point(data = MDSp, size = 4, aes(color=STAGE, shape=POP, fill=TISSUE, alpha=1), stroke = 1.5)+
  geom_hline(yintercept = 0, lty = "dotted", color="grey", cex=1) +    geom_vline(xintercept = 0, lty = "dotted", color="grey", cex=1) +
  theme( legend.position="none", 
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         panel.background = element_rect(fill = "transparent", colour = NA))+
  scale_fill_manual(values= c("#58D68D", "#EC7063","#F4D03F","#85C1E9"), labels=c("GILL","DIGESTIVE","TUNIC","WATER"))+
  scale_shape_manual(values=c(25,23,21)) +
  scale_color_manual(values=c("purple4","#DB7093","steelblue")) +
  theme_classic()+
  labs(title = paste("Stress  ", 0.0366, sep = ""))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(colour = "black", size = 12, face = "bold")) +
  theme(axis.title = element_text(size = 16, colour = "black", family = "Helvetica", face = "bold"))   
dev.off()