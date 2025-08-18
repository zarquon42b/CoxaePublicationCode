require(Morpho)
require(Rvcg)
require(mesheR)
require(rgl)
require(RvtkStatismo)

## load necessary custom functions and reference data
source("./functions.r")
load("./Data/reference.RDATA")
predictorrow <- readRDS("./Data/predictorrow.RDS")


## load landmarks of extant taxa
coxae_data_hum <- readRDS("./Data/coxae_data_homo.RDS")
nhum <- dim(coxae_data_hum)[3]
coxae_data_pan <- readRDS("./Data/coxae_data_pan.RDS")
npan <- dim(coxae_data_pan)[3]
coxae_data_gor <- readRDS("./Data/coxae_data_gorilla.RDS")
ngor <- dim(coxae_data_gor)[3]


### load vertices of response region
pubVerts <- readRDS("./Data/pubVerts.RDS")

### load fossil data
fossilArray <- readRDS("./Data/fossilArray.RDS")
testingMesh <- list.files("./TestingData/",full.names=T)
fossilMesh <- list()
for(i in 1:5)     
     fossilMesh[[i]] <- vcgImport(testingMesh[grep(dimnames(fossilArray)[[3]][i],testingMesh)])

### Predict using Human SSM
## load SSM
load("./Data/humMod.RDATA")


### reduce SSM to 10 PCs
humMod10 <- statismoReducedVariance(humModRed,npc=10)
## Predict os pubis for all fossils using only predictor landmarks
predictedFossilsHum <- list()
for (i in 1:dim(fossilArray)[3]) {
    f.lm <- fossilArray[,,i]
    f.m <- fossilMesh[[i]]
    predictedFossilsHum[[i]] <- predictFossil(f.lm,humMod10,humMod.lm,humProcDef,use.lm = predictorrow,useweights = F,target.m = f.m,evalverts=pubVerts)
}
names(predictedFossilsHum) <- dimnames(fossilArray)[[3]]

## Visualize Predictions - in this example the first prediction
visualizePredict(predictedFossilsHum[[1]])


### Predict using Human/Pan SSM with a shifting weighted mean
load("./Data/humpanMod.RDATA")
### reduce SSM to 10 PCs
humpanMod10 <- statismoReducedVariance(humpanModRed,npc=10)
## run Procrustes on predictor landmarks of training data
humpanProcDef <- procSym(bindArr(coxae_data_hum,coxae_data_pan,along = 3)[predictorrow,,])
meanlistHP  <- list(hum=arrMean3(humpanProcDef$rotated[,,1:nhum]),pan=arrMean3(humpanProcDef$rotated[,,(nhum+1):(nhum+npan)]))


predictedFossilsHumPan <- list()
for (i in 1:dim(fossilArray)[3]) {
    f.lm <- fossilArray[,,i]
    f.m <- fossilMesh[[i]]
    predictedFossilsHumPan[[i]] <- predictFossil(f.lm,humpanMod10,humpanMod.lm,humpanProcDef,means=meanlistHP,meanMesh=humpanMList,use.lm = predictorrow,useweights = T,target.m = f.m,evalverts=pubVerts)
}

names(predictedFossilsHumPan) <- dimnames(fossilArray)[[3]]
## Visualize Predictions - in this example the first prediction
visualizePredict(predictedFossilsHumPan[[1]])

### Predict using Human/Pan/Gorilla SSM with a shifting weighted mean

load("./Data/humpangorMod.RDATA")
humpangorMod10 <- statismoReducedVariance(humpangorModRed,npc=10)
humpangorProcDef <- procSym(bindArr(coxae_data_hum,coxae_data_pan,coxae_data_gor,along = 3)[predictorrow,,])
meanlistHPG  <- list(hum=arrMean3(humpangorProcDef$rotated[,,1:nhum]),pan=arrMean3(humpangorProcDef$rotated[,,(nhum+1):(nhum+npan)]),gor=arrMean3(humpangorProcDef$rotated[,,(nhum+npan+1):(nhum+npan+ngor)]))


predictedFossilsHumPanGor <- list()
for (i in 1:dim(fossilArray)[3]) {
    f.lm <- fossilArray[,,i]
    f.m <- fossilMesh[[i]]
    predictedFossilsHumPanGor[[i]] <- predictFossil(f.lm,humpangorMod10,humpangorMod.lm,humpangorProcDef,means=meanlistHPG,meanMesh=humpangorMList,use.lm = predictorrow,useweights = T,target.m = f.m,evalverts=pubVerts,getMod=F)
}

names(predictedFossilsHumPanGor) <- dimnames(fossilArray)[[3]]
visualizePredict(predictedFossilsHumPanGor[[1]])


### compute Errors
### load fossil surface meshes, registered to the SSMs
fossilMatched <- list()
testingMatched <- list.files("./matchedFossils/",full.names=T)
for(i in 1:5) {    
    fossilMatched[[i]] <- vcgImport(testingMatched[grep(dimnames(fossilArray)[[3]][i],testingMatched)])
}

pubref <- rmVertex(reference.m,pubVerts,keep = T)
set.seed(42)
subsample <- fastKmeans(pubref,k=1000)$selected
subind <- vcgKDtree(reference.m,vert2points(pubref)[subsample,],k=1)$index


## Human SSM
predictedFossilsHumError <- list()
for (i in 1:length(predictedFossilsHum)) {
    predictedFossilsHumError[[i]] <- computeErrorNew(predictedFossilsHum[[i]]$fossil2tar,fossilMatched[[i]],subsample=subind,iterations=10,tol=0,silent=T)
}

## Human/Pan SSM
predictedFossilsHumPanError <- list()
for (i in 1:length(predictedFossilsHumPan)) {
    predictedFossilsHumPanError[[i]] <- computeErrorNew(predictedFossilsHumPan[[i]]$fossil2tar,fossilMatched[[i]],subsample=subind,iterations=10,tol=0,silent=T)
}

## Human/Pan/Gorilla SSM
predictedFossilsHumPanGorError <- list()
for (i in 1:length(predictedFossilsHumPanGor)) {
    predictedFossilsHumPanGorError[[i]] <- computeErrorNew(predictedFossilsHumPanGor[[i]]$fossil2tar,fossilMatched[[i]],subsample=subind,iterations=10,tol=0,silent=T)
}

## extract mean errors
## NOTE as random subsampling is involved, the error values can slightly vary (also compared to the numbers in the article)

humModError <- sapply(predictedFossilsHumError,function(x) x <- x$meanerror)
humpanModError <- sapply(predictedFossilsHumPanError,function(x) x <- x$meanerror)
humpangorModError <- sapply(predictedFossilsHumPanGorError,function(x) x <- x$meanerror)

newErrorTable <- data.frame(hum=humModError,humpan=humpanModError, humpangor=humpangorModError)
row.names(newErrorTable) <- dimnames(fossilArray)[[3]]
boxplot(newErrorTable)


### Visualize Metaspecies by sampling from hum/pan/gor SSMs
threescore <- GetPCScores(humpangorMod)[,1:3]
myscores <- NULL
for (i in 1:100) myscores <- rbind(myscores,ComputeCoefficients(humpangorMod10, DrawSample(humpangorMod10)))

allscores <- rbind(threescore,myscores[,1:3])
groups <- c(rep("hum",nhum),rep("pan",npan),rep("gor",ngor))
groups <- c(groups,rep("Sample",100))
car::sp(allscores[,1],allscores[,2],group=groups,smooth=F,pch=c(15:18),regLine=F,cex=2,xlab="PC1",ylab="PC2",col= c( "magenta" , "cyan" ,"orange","green"))


