#' ## Necessary objects in global workspace
#' 
#' fossil.lm:    landmarks on fossil remain
#' mymod:        SSM
#' mymod.lm:     landmarks on SSM mean
#' proc:         Procrustes Analysis based on landmarks contained in fossil.lm
#' use.lm:       indices of subset to use in fossil.lm
#' means:        list of means from different groups contained in SSM and proc
#' meanMesh:     list of vertices of meshes belonging to means in the same order
#' useweights:   shift SSM mean to weighted average from meanMesh based on inverse Procrustes distance#'
#' 
predictFossil <- function(fossil.lm,  mymod, mymod.lm, proc, use.lm=NULL, meanMesh=NULL, means=NULL, useweights=T,target.m=NULL,subs=500,evalverts=NULL,getMod=F,uprange=1) {
    weightsL=NULL
    meanmesh=NULL
    align2tar=NULL
    error <- NULL
    outmod=NULL
    fossilPred2=NULL
    fossilPred2Back=NULL
    if (useweights) {
        fossil2proc <- align2procSym(proc,fossil.lm[use.lm,])
        dist2mean <- sapply(means,kendalldist,fossil2proc)
        weightsL <- 1/dist2mean
        weightsL <- weightsL/sum(weightsL)
        meanmesh <- lapply(1:length(means),function(x) x <- weightsL[x]*meanMesh[[x]])
        meanmesh <- apply(list2array(meanmesh),1:2,sum)
        ##meanmesh <- meanMeshes(meanmesh)
        SetMeanVector(mymod) <- as.vector(t(meanmesh))
    }
    
    fossilPredict <- PredictSample(mymod,lmModel=mymod.lm[use.lm,],lmDataset=fossil.lm[use.lm,])
    if (!is.null(target.m))
        align2tar <- icp(fossilPredict,target.m,subsample=2000,iterations=100,rho=pi/2,uprange = uprange)
    if (!is.null(evalverts)) {
        fossilPub <- rmVertex(align2tar,evalverts,keep = T)
        fossiPubSample <- vcgSample(fossilPub,SampleNum = subs)
        foss2tarErr <- vcgClostKD(fossilPub,target.m)
        error <- list()
        error$mean <- mean(abs(foss2tarErr$quality))
        error$quant <- quantile(abs(foss2tarErr$quality),probs = c(.05,.5,.95))
        error$raw <- foss2tarErr$quality
    }
    
    if (getMod) {
        f2modTrans <- computeTransform(mymod.lm[use.lm,],fossil.lm[use.lm,])
        f2mod.lm <- applyTransform(fossil.lm,f2modTrans)
       
        postmod <- statismoConstrainModel(mymod,sample =f2mod.lm[use.lm,],pt=mymod.lm[use.lm,],ptValueNoise = 2,computeScores = F)
        if (!is.null(target.m)) {
            postMean <- DrawMean(postmod)
            postMean.lm <- transferPoints(mymod.lm,DrawMean(mymod),postMean)
            f2mod <- icp(target.m,postMean,lm1=fossil.lm, lm2=postMean.lm,iterations=50,subsample = 500,rhotol = pi/2,uprange=uprange,getTransform = T)
        }       
        outmod=list(model=postmod,f2mod=f2mod)
    }
 
    out <- list(predict=fossilPredict,mod=mymod,wt=weightsL,mm=meanmesh,fossil2tar=align2tar,target.m=target.m,error=error,postMod=outmod,predict2=fossilPred2Back)
    class(out) <- "FossilPredict"
    return(out)
    
}

visualizePredict <- function(x,col="green",colpred="orange",window=FALSE,usealigned=T,use2=FALSE) {
    if (window)
        open3d()
    else {
        if (rgl.cur())
            clear3d()
    }
    if (usealigned)
        wire3d(x$fossil2tar,col=colpred,specular=1)
    else
        wire3d(x$predict2,col=colpred,specular=1)

    wire3d(x$target.m,col=col,specular=1)

}

exportMod <- function(x,file) {
    statismoSaveModel(x$postMod$model,file)
    vcgPlyWrite(x$postMod$f2mod$mesh,file)
    
}

exportPrediction <- function(x,filename = dataname,modname=NULL,colpred="orange") {
    dataname <- deparse(substitute(x))
    filename <- path.expand(as.character(filename))
    x$target.m <- colorMesh(x$target.m,col="white")
    x$fossil2tar <- colorMesh(x$fossil2tar,col=colpred)
    vcgPlyWrite(x$target.m,filename = paste0(filename,"_orig"))
    vcgPlyWrite(x$fossil2tar,filename = paste0(filename,"_prediction","_",modname))

}

exportPredictionFromList <- function(x,i,folder=".",modname=NULL,colpred="orange") {
    filename=paste0(folder,"/",names(x)[i])
    print(filename)
    exportPrediction(x[[i]],filename=filename,modname=modname,colpred=colpred)
  #  png(paste0(filename,"_",modname,".png"),width=500,height=500)
 #   boxplot(abs(x[[i]]$error$raw),main=names(x))
  #  dev.off()
    
}

### compute error my allowing to slide 1000 subsampled vertices on the registered target shape.
computeErrorNew <- function(x,target,subsample=NULL, slide=T,...) {
    if (is.null(subsample)) {
        subsample <- 1:nverts(x)
    }
    x.coo <- vert2points(x)[subsample,]
    target.coo <- vert2points(target)[subsample,]
    if (slide) {
        x.coo.slide <- relaxLM(x.coo,SMvector = 1:nrow(x.coo),surp=1:nrow(x.coo),mesh=x,reference=target.coo,...)
    }
    error <- sqrt(rowSums
    ((target.coo-x.coo.slide)^2))
    meanerror <- mean(error)
    return(list(x.slide=x.coo.slide,target=target.coo,error=error,meanerror=meanerror))
    
    
    
}
