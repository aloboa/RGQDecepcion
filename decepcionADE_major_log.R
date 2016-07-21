#' Discriminant Analysis of Geochemical composition vs. Eruption Phase (IIa)
#' ========================================================================
#' Major compounds Decepcion Island (Adelina Geyer 2014)
#' -------------------------------
#' 
#' * Agustin.Lobo@ictja.csic.es
#' * Data by ageyertraver@gmail.com
#' * 20160721
#' 
#' **Goal: do samples defined by their geochemical composition cluster according to field-diagnosed eruption phase?**
#' 
knitr::opts_chunk$set(fig.path='figure/decepcionADE_major-')
require(MASS)
require(ggplot2)
require(reshape2)
require(plyr)
require(scales)
require(robCompositions)
require(subselect)

rwd <- "/media/alobo/LACIE500/Adelina/GQDecepcion/RGQDecepcion"
dirdata1 <- "/media/alobo/LACIE500/Adelina/GQDecepcion/GQDecepcionData"

#' ## 1. Data handling
#' ### 1.1 Data reading
decepade.ori <- read.csv(file.path(dirdata1,"resultados_analiticos_Decepcion.csv"),header=TRUE,stringsAsFactors = FALSE)
head(decepade.ori,1)
names(decepade.ori)
names(decepade.ori)[1:4] <- c("Sample","Longitude","Latitude","Phase")
decepade.ori$Longitude <- -decepade.ori$Longitude
decepade.ori$Latitude  <- -decepade.ori$Latitude
names(decepade.ori)[16] <- "LOI"
names(decepade.ori)[c(8,20)] <- "FeO*"
table(decepade.ori$Phase)
decepade.ori$Phase <- factor(decepade.ori$Phase,levels=c("pre","syn","post","dique","Mush","pipe"))
#' Check percentages were calculated without LOI
a.suma <- apply(decepade.ori[,5:14],1,sum)
a.suma
a.percent <- decepade.ori[,5:14]*100/t(a.suma)
a.percent[1:3,]
decepade.ori[1:3,17:26]
dim(decepade.ori)
cbind(1:ncol(decepade.ori),names(decepade.ori))
decepade <- decepade.ori[,c(1:14,28:74)]
decepade[,5:14] <- a.percent
options(width=180)
head(decepade,1)
summary(decepade)
options(width=80)

#' ### 1.2 Data transformation: Centered Log-Ratio
#+ fig.width=10, fig.height=8
barplot(apply(decepade[,5:14],1,sum))
decepade.cenLR <- cenLR(decepade[,5:14])
decepadecen <- decepade[,1:14]
decepadecen[,5:14] <- decepade.cenLR$x.clr
head(decepadecen,1)
decepadecen.melt <- melt(decepadecen,id=1:4,variable.name="Compound",value.name="LRC_Concentration")
head(decepadecen.melt,10)
ggplot(data=decepadecen.melt) +
    geom_histogram(aes(x=LRC_Concentration)) +
    facet_wrap(~Compound,scales="free")
pairs(decepadecen[,-(1:6)])

#' ## 2. LDA of Major Compounds
decepadecen.lda <- lda(decepadecen[,5:14], grouping=decepadecen$Phase)
decepadecen.lda <- lda(prcomp(decepadecen[,5:14])$x, grouping=decepadecen$Phase)
decepadecen.lda <- lda(prcomp(decepadecen[,5:14])$x[,-10], grouping=decepadecen$Phase)
decepadecen.ldaclas <- predict(decepadecen.lda)$class
decepadecen.ldc <- predict(decepadecen.lda)$x
table(decepadecen$Phase)
adecolorines <- as.character(mapvalues(decepadecen.ldaclas,from=c("pre","syn","post","dique","Mush","pipe"),
                                    to=c("green","red","blue","brown","pink","black")))
#+ out.width='50%', fig.show='hold'
#bmp("adeLDAmajor1.bmp",width=969,height=800)
plot(decepadecen.ldc, type="n", ,main="LDA Plot")
text(decepadecen.ldc,col=adecolorines, labels =decepadecen$Phase, cex=0.75)
legend("topright",title="Predicted Phases",title.col="black",
       legend=c("pre","syn","post","dique","Mush","pipe"),
       text.col=c("green","red","blue","brown","pink","black")
       ,cex=0.75)
#dev.off()
#bmp("adeLDAmajor2.bmp",width=969,height=800)
plot(decepadecen.ldc,col=adecolorines,pch=" ",main="LDA Plot")
text(decepadecen.ldc,col=adecolorines, labels =decepadecen$Sample, cex=0.75)
legend("topright",title="Predicted Phases",title.col="black",
       legend=c("pre","syn","post","dique","Mush","pipe"), 
       pch=19,col=c("green","red","blue","brown","pink","black"),bty="o",cex=0.75)
#dev.off()
