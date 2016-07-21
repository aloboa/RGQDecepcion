#' Discriminant Analysis of Geochemical composition vs. Eruption Phase (I)
#' ========================================================================
#' Decepcion Island (Smellie 2001)
#' -------------------------------
#' * Agustin.Lobo@ictja.csic.es
#' * Data by ageyertraver@gmail.com
#' * 20160720
#' 
#' **Goal: do samples defined by their geochemical composition cluster according to field-diagnosed eruption phase?**
#' 
knitr::opts_chunk$set(fig.path='figure/smellie-')
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
smellie <- read.csv(file.path(dirdata1,"smellie2001.csv"),header=TRUE,stringsAsFactors = FALSE)
head(smellie,1)
smellie$Longitude <- -smellie$Longitude
smellie$Latitude <- -smellie$Latitude

smellie.phases <- read.csv(file.path(dirdata1,"smellie_phases.csv"),header=TRUE,stringsAsFactors = FALSE)
head(smellie.phases,1)
#' Make sure the ordering is the same:
diag(table(smellie$Sample,smellie.phases$sample))
smellie$Phase <- smellie.phases$phase
smellie$Phase <- factor(smellie$Phase,levels=c("pre","syn","post"))
head(smellie,1)
dim(smellie)
smellie <- smellie[,c(1,17,2:16)]
head(smellie,1)
smellie[1:10,1:6]
summary(smellie)
#' [Fe0] is 0 in all samples, we discard it:
smellie <- smellie[,-9]
smellie.melt <- melt(smellie,id=1:6,variable.name="Compound",value.name="Concentration")

#' ### 1.2 Data transformation: Centered Log-Ratio
#+ fig.width=10, fig.height=8
barplot(apply(smellie[,6:15],1,sum))
smellie.cenLR <- cenLR(smellie[,6:15])
smelliecen <- smellie
smelliecen[,6:15] <- smellie.cenLR$x.clr
head(smelliecen,1)
smelliecen.melt <- melt(smelliecen,id=1:6,variable.name="Compound",value.name="LRC_Concentration")
head(smelliecen.melt,10)
ggplot(data=smelliecen.melt) +
    geom_histogram(aes(x=LRC_Concentration)) +
    facet_wrap(~Compound,scales="free")
pairs(smelliecen[,-(1:6)])

#' ## 2. LDA
smelliecen.lda <- lda(smelliecen[,6:16], grouping=smelliecen$Phase)
smelliecen.lda <- lda(prcomp(smelliecen[,6:16])$x, grouping=smelliecen$Phase)
smelliecen.lda <- lda(prcomp(smelliecen[,6:16])$x[,-11], grouping=smelliecen$Phase)
smelliecen.ldaclas <- predict(smelliecen.lda)$class
smelliecen.ldc <- predict(smelliecen.lda)$x
colorines <- as.character(mapvalues(smelliecen.ldaclas,from=c("pre","syn","post"),
                       to=c("green","red","blue")))
#+ fig.width=10, fig.height=8, fig.scap="Colors: observed phases; labels: predicted phases."
plot(smelliecen.ldc, type="n", ,main="LDA Plot")
text(smelliecen.ldc,col=colorines, labels =smelliecen$Phase, cex=0.75)
legend("topright",title="Predicted Phases",title.col="black",legend=c("pre","syn","post"),text.col=c("green","red","blue"),bty="o",cex=0.75)
#' Labels refer to observed eruptive phases, while colors refer to LDA-predicted eruptive phases. Note 3 discrepancies:
#' 
#' * 1 observed "syn" predicted as "pre"
#' * 1 observed "pre" predicted as "post"
#' * 1 observed "syn" predicted as "post".  
#+ fig.width=10, fig.height=8, fig.scap="Colors: observed phases; labels: samples."
plot(smelliecen.ldc,col=colorines,pch=" ",main="LDA Plot")
text(smelliecen.ldc,col=colorines, labels =smelliecen$Sample, cex=0.75)
legend("topright",title="Predicted Phases",title.col="black",legend=c("pre","syn","post"), pch=19,col=c("green","red","blue"),bty="o",cex=0.75)
#' According to this plot, the samples correspond to:
#' 
#' * 751.2 observed "syn" predicted as "pre"
#' * 793.9 observed "pre" predicted as "post"
#' * 751.3 observed "syn" predicted as "post"
#' 
#' ### 2.2 Centroids
#+ fig.width=12, fig.height=8, fig.scap="Colors: observed phases; labels: predicted phases."
ggplot(data=smelliecen.melt) +
    geom_boxplot(aes(x=Compound,y=LRC_Concentration)) +
    facet_wrap(~Phase)
ggplot(data=smellie.melt) +
    geom_boxplot(aes(x=Compound,y=Concentration)) +
    facet_wrap(~Phase)

#' ## 3. MANOVA
#' Check multi-variate significance
summary(manova(data.matrix(smelliecen[,-(1:5)]) ~ smelliecen$Phase),test="Wilks")
#' This seems to be a consequence of the co-llinearity again. So we use the same solution:
summary(manova(prcomp(smelliecen[,6:16])$x[,-11]~ smelliecen$Phase),test="Wilks")
#' **Phase groups significantly differ in terms of Geochemical composition**
#' 
#' ## 4. Selection of variables
#+ eval=FALSE
smelliecenHmat <- ldaHmat(smelliecen[,-(1:5)], grouping=smelliecen$Phase)
smelliecen.eleaps <- eleaps(smelliecenHmat$mat,kmin=2,kmax=10,H=smelliecenHmat$H,r=smelliecenHmat$r,crit="tau2",timelimit=15)
smelliecen.eleaps
#' Plot of quality of subsets
#+ fig.width=12, fig.height=6, fig.show='hold'
par(mfrow=c(1,2))
plot(2:10,smelliecen.eleaps$bestvalues,xlab="Nb. of variables", ylab="Tau2 value")
plot(3:10,diff(smelliecen.eleaps$bestvalues)*100/smelliecen.eleaps$bestvalues[-9],
     xlab="Nb. of variables", ylab="Delta of Tau2")
#' Table of subsets
options(width=180)
subsets <- smelliecen.eleaps$bestsets
subsets[subsets==0] <- NA
subsets2 <- names(smelliecen)[subsets+5]
dim(subsets2) <- dim(subsets)
colnames(subsets2) <- colnames(subsets)
subsets2
options(width=80)

#+ out.width='50%', fig.show='hold'
ggplot(data=data.frame(smelliecen.ldc,Phase=smelliecen$Phase)) +
    geom_text(aes(x=LD1,y=LD2,label=Phase,col=Phase),size=3) +
    ggtitle("LDA with all variables")
for(i in c(5,7,10)){
    selvar <- smelliecen.eleaps$bestsets[i-1,]
    selvar <- selvar[selvar!=0]+5
    delme <- paste(names(smelliecen)[selvar], collapse=" ")
    smelliecensub.lda <- lda(smelliecen[,selvar],grouping=smelliecen$Phase)
    smelliecensub.ldc <- data.matrix(predict(smelliecensub.lda)$x)
    smelliecensub.ldc <- data.frame(smelliecen[,1:5],smelliecensub.ldc)
    print(ggplot(data=smelliecensub.ldc) +
              geom_text(aes(x=LD1,y=LD2,label=Phase,col=Phase),size=3) +
              ggtitle(paste("LDA with", delme)))
}
