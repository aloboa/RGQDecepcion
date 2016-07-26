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

#' ## 2. LDA of Major Compounds (all Phases considered)
decepadecen.lda <- lda(decepadecen[,5:14], grouping=decepadecen$Phase)
decepadecen.pc <- prcomp(decepadecen[,5:14])$x
decepadecen.lda <- lda(decepadecen.pc[,-10], grouping=decepadecen$Phase)
decepadecen.lda <- lda(decepadecen.pc[,-10], grouping=decepadecen$Phase)
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

#' ## 2. LDA of Major Compounds (Phases pre and pipe discarded at model fit )
#' We exclude from the model:
#' 
#' * All "pre", as they do not cluster at all (they cannot be considered a group in GQ terms)
#' * DI-48: conflicting diagnostic ("pre" for some authors, "post" for some others)
#' * DI-66: early phase? weird as "post" in any case
#' * All samples classified with a prob < 0.95
#' 
#' All samples will be reprojected to the LDA space to assess where thy go
decepadecen$sel <- 1
decepadecen$sel[decepadecen$Phase=="pre" | decepadecen$Phase== "pipe"] <- 0
decepadecen$sel[decepadecen$Sample=="DI-48"] <- 0
decepadecen$sel[decepadecen$Sample=="DI-66"] <- 0
decepadecensel <- decepadecen[decepadecen$sel==1,]
decepadecensel$Phase <- droplevels(decepadecensel$Phase)
decepadecensel.lda <- lda(decepadecen.pc[decepadecen$sel==1,-10], grouping=decepadecensel$Phase)
decepadecensel.ldaclas <- predict(decepadecensel.lda)$class
decepadecensel.ldc <- predict(decepadecensel.lda)$x
decepadecensel.ldp <- predict(decepadecensel.lda)$posterior
decepadecensel.ldp <- apply(decepadecensel.ldp,1,max)
#' Now we discard samples classified with prob<0.95
decepadecensel$sel[decepadecensel.ldp<0.95] <-0
bad <- decepadecensel$Sample[decepadecensel.ldp<0.95]
exclude <- decepadecen$Sample%in%bad
decepadecen$sel[exclude] <-0 #we need this to keep a record of all excluded samples
decepadecensel <- decepadecen[decepadecen$sel==1,]
decepadecensel$Phase <- droplevels(decepadecensel$Phase)
#' And repeat the lda:
decepadecensel.lda <- lda(decepadecen.pc[decepadecen$sel==1,-10], grouping=decepadecensel$Phase)
decepadecensel.ldaclas <- predict(decepadecensel.lda)$class
decepadecensel.ldc <- predict(decepadecensel.lda)$x
decepadecensel.ldp <- predict(decepadecensel.lda)$posterior
decepadecensel.ldp <- apply(decepadecensel.ldp,1,max)
table(decepadecensel$Phase)
adecolorines2 <- as.character(mapvalues(decepadecensel.ldaclas,from=c("syn","post","dique","Mush"),
                                       to=c("red","blue","brown","cyan")))
#+ out.width='50%', fig.show='hold'
#bmp("adeLDAmajor1.bmp",width=969,height=800)
plot(decepadecensel.ldc, type="n", ,main="LDA Plot")
text(decepadecensel.ldc,col=adecolorines2, labels =decepadecensel$Phase, cex=0.7)
legend("bottomleft",title="Predicted Phases",title.col="black",
       legend=c("syn","post","dique","Mush"),
       text.col=c("red","blue","brown","pink")
       ,cex=0.75)
#dev.off()
plot(decepadecensel.ldc, type="n", ,main="LDA Plot (probabilities)")
text(decepadecensel.ldc,col=adecolorines2, labels = round(decepadecensel.ldp,2), cex=0.7)
legend("bottomleft",title="Predicted Phases",title.col="black",
       legend=c("syn","post","dique","Mush"),
       text.col=c("red","blue","brown","pink")
       ,cex=0.75)
#dev.off()
#bmp("adeLDAmajor2.bmp",width=969,height=800)
plot(decepadecensel.ldc,col=adecolorines2,pch=" ",main="LDA Plot (labels)")
text(decepadecensel.ldc,col=adecolorines2, labels =decepadecensel$Sample, cex=0.5)
legend("bottomleft",title="Predicted Phases",title.col="black",
       legend=c("syn","post","dique","Mush"), 
       pch=19,col=c("red","blue","brown","pink"),bty="o",cex=0.75)
#dev.off()

#' ### 3.2. Composition of Phases
#+ fig.width=12, fig.height=10
decepadesel <- decepade[decepadecen$sel==1,c(1:4,5:14)] 
decepadesel$class <- decepadecensel.ldaclas
decepadesel.melt <- melt(decepadesel,id=c(1:4,15),variable.name="Compound",value.name="Concentration")
head(decepadesel.melt,3)
ggplot(data=decepadesel.melt) +
    geom_boxplot(aes(x=Compound,y=Concentration)) +
    facet_wrap(~Phase) +
    ggtitle("Composition (major compounds) of modeled Phases")

#' ### 3.3 Selection of variables
#TBD

#' ### 3.4 Predict with all samples
decepadecensel.ldaclasall <- predict(decepadecensel.lda,decepadecen.pc[,-10])$class
decepadecensel.ldcall <- predict(decepadecensel.lda,decepadecen.pc[,-10])$x
decepadecensel.ldpall <- predict(decepadecensel.lda,decepadecen.pc[,-10])$posterior
decepadecensel.ldpall <- apply(decepadecensel.ldpall,1,max)
adecolorines3 <- as.character(mapvalues(decepadecensel.ldaclasall,from=c("syn","post","dique","Mush"),
                                        to=c("red","blue","brown","cyan")))

#+ fig.width=10, fig.height=10
#bmp("adeLDAmajor1.bmp",width=969,height=800)
plot(decepadecensel.ldcall, type="n", ,main="LDA Plot",xlim=c(-6,5))
text(decepadecensel.ldcall[decepadecen$sel==1,],col=adecolorines3[decepadecen$sel==1], labels =decepadecen$Phase[decepadecen$sel==1], cex=1)
text(decepadecensel.ldcall[decepadecen$sel==0,], col=adecolorines3[decepadecen$sel==0],labels =decepadecen$Phase[decepadecen$sel==0], cex=0.75)
text(decepadecensel.ldcall[decepadecen$sel==0,], labels =round(decepadecensel.ldpall[decepadecen$sel==0],2), 
     pos=4, offset=0.5,cex=0.75)
text(decepadecensel.ldcall[decepadecen$sel==0,], labels =decepadecen$Sample[decepadecen$sel==0], 
     pos=1, offset=0.3,cex=0.75)
#'
#' * Larger Phase labels indicate samples used to fit the LDA
#' * Smaller Phase labels indicate samples not used in the fit
#'  * Color indicates predicted Phase
#'  * Text indicates field diagnosed Phase
#'  * Number (in black to the right) indicates probability of the prediction
#'  * Sample ID is displayed below the Phase label
