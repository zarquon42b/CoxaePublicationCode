plotColramp <- function(x) {
    tol <- x$params$tol
    colramp <- x$colramp
    diffo <- ((colramp[[2]][2] - colramp[[2]][1])/2)
    image(colramp[[1]], colramp[[2]][-1] - diffo, t(colramp[[3]][1, 
                                                                 -1]) - diffo, col = colramp[[4]], useRaster = TRUE, ylab = "Distance in mm", 
          xlab = "", xaxt = "n")
    if (!is.null(tol)) {
        if (sum(abs(tol)) != 0) {
            image(colramp[[1]], c(tol[1], tol[2]), matrix(c(tol[1], 
                                                            tol[2]), 1, 1), col = "green", useRaster = TRUE, 
                  add = TRUE)
        }
    }
}


reconstructPooled <- function(mylm,mymesh,proc) {
    atts <- attributes(proc)
    scale <- atts$scale
    
    lm2proc <- align2procSym(mylm,proc)
    
    distH <- kendalldist(lm2proc,humMean)
    
}

getSSMerror <- function() {
    origarrAll <<- meshlist2array(matched)##meshslider$dataslide[,,]
    origarr <<- origarrAll[pubVerts,,]
    
    perVertexError <<- sapply(1:n, function(x) x <<- sqrt(rowSums((vert2points(predSSM[[x]])-vert2points(matched[[x]]))[pubVerts,])^2))
    perVertexErrorMean <<- rowMeans(perVertexError)
    perVertexErrorPlainSQ <<- sapply(1:n, function(x) x <<- (rowSums(((vert2points(predSSM[[x]])-vert2points(matched[[x]]))[pubVerts,]))^2))
    RMSEPlain <<- sqrt(mean(perVertexErrorPlainSQ))
    PVEPlain <<- mean(perVertexError)
    PVEPlain95 <<- quantile(colMeans(perVertexError),probs = .95)
    distvec <<- rep(1000,nverts(reference))
    distvec[pubVerts] <<- perVertexErrorMean
    mDPlain <<- meshDist(reference,distvec=distvec,tol=1,to=max(distvec[pubVerts])+2,plot = F)
    
    ## ERROR TPS
    perVertexErrorTPS <<- sapply(1:n, function(x) x <<- sqrt(rowSums((vert2points(predSSMTPS[[x]])-vert2points(aligned[[x]]))[pubVerts,])^2))
    perVertexErrorMeanTPS <<- rowMeans(perVertexErrorTPS)
    perVertexErrorTPSSQ <<- sapply(1:n, function(x) x <<- rowSums(((vert2points(predSSMTPS[[x]])-vert2points(aligned[[x]]))[pubVerts,])^2))
    RMSETPS <<- sqrt(mean(perVertexErrorTPSSQ))
    PVETPS <<- mean(perVertexErrorTPS)
    PVETPS95 <<- quantile(colMeans(perVertexErrorTPS),probs = .95)
    distvecTPS <<- rep(1000,nverts(reference))
    distvecTPS[pubVerts] <<- perVertexErrorMeanTPS
    mDTPS <<- meshDist(reference,distvec=distvecTPS,tol=1,to=max(distvec[pubVerts])+2,plot = F)
    
    ##perVertexError <<- sapply(1:n, function(x) x <<- sqrt(sum((vert2points(predSSM[[x]])-vert2points(alignedIlium[[x]]))[pubVerts,])^2))
    
### ERROR ESTIMATION: align pubis separately
    predSSMarr <<- meshlist2array(predSSM)[pubVerts,,]
    predalign <<- procSym(bindArr(predSSMarr,origarrAll[pubVerts,,],along=3),CSinit=F,scale = F)
    perVertexErrorAlign <<- sapply(1:n, function(x) x <<- sqrt(rowSums(predalign$rotated[,,x]-predalign$rotated[,,x+n])^2))
    perVertexErrorMeanAlign <<- rowMeans(perVertexErrorAlign)
    perVertexErrorAlignSQ <<- sapply(1:n, function(x) x <<- (rowSums(predalign$rotated[,,x]-predalign$rotated[,,x+n])^2))
    RMSEAlign <<- sqrt(mean(perVertexErrorAlignSQ))
    PVEAlign <<- mean(perVertexErrorAlign)
    PVEAlign95 <<- quantile(colMeans(perVertexErrorAlign),probs = .95)
    distvecAlign <<- rep(1000,nverts(reference))
    distvecAlign[pubVerts] <<- perVertexErrorMeanAlign
    mDalign <<- meshDist(reference,distvec=distvecAlign,tol=1,to=max(distvec[pubVerts])+2,plot = T)
    ##plot(hclust(dist(predalign$PCscores),method="ward.D"))
    
### ERROR ESTIMATION: align pubis separately TPS
    predSSMarrTPS <<- meshlist2array(predSSMTPS)[pubVerts,,]
    dimnames(predSSMarrTPS)[[3]] <<- paste0(dimnames(predSSMarrTPS)[[3]],"_pred")
    predalignTPS <<- procSym(bindArr(predSSMarrTPS,origarr,along=3),CSinit=F,scale = F)
    perVertexErrorAlignTPS <<- sapply(1:n, function(x) x <<- sqrt(rowSums(predalignTPS$rotated[,,x]-predalignTPS$rotated[,,x+n])^2))
    perVertexErrorAlignTPSSQ <<- sapply(1:n, function(x) x <<- (rowSums(predalignTPS$rotated[,,x]-predalignTPS$rotated[,,x+n])^2))
    RMSEAlignTPS <<- sqrt(mean(perVertexErrorAlignTPSSQ))
    PVEAlignTPS <<- mean(perVertexErrorAlignTPS)
    PVEAlignTPS95 <<- quantile(colMeans(perVertexErrorAlignTPS),probs = .95)
    perVertexErrorMeanAlignTPS <<- rowMeans(perVertexErrorAlignTPS)
    distvecAlignTPS <<- rep(1000,nverts(reference))
    distvecAlignTPS[pubVerts] <<- perVertexErrorMeanAlignTPS
    mDalignTPS <<- meshDist(reference,distvec=distvecAlignTPS,tol=1,to=max(distvec[pubVerts])+2,plot = F)
    ##plot(hclust(dist(predalignTPS$PCscores),method="ward.D"))
    
### ERROR ESTIMATION: realign entire structure 
    predSSMarrAll <<- meshlist2array(predSSM)[,,]
    predalignAll <<- procSym(bindArr(predSSMarrAll,origarrAll,along=3),CSinit=F,scale = F)
    perVertexErrorAlignAll <<- sapply(1:n, function(x) x <<- sqrt(rowSums(predalignAll$rotated[pubVerts,,x]-predalignAll$rotated[pubVerts,,x+n])^2))
    perVertexErrorMeanAlignAll <<- rowMeans(perVertexErrorAlignAll)
    ## RMSE
    perVertexErrorAllSQ <<- sapply(1:n, function(x) x <<- (rowSums(predalignAll$rotated[pubVerts,,x]-predalignAll$rotated[pubVerts,,x+n])^2))
    RMSEAll <<- sqrt(mean(perVertexErrorAllSQ))
    PVEAll <<- mean(perVertexErrorAlignAll)
    PVEAll95 <<- quantile(colMeans(perVertexErrorAlignAll),probs = .95)
    
    ## Distance vector
    distvecAlignAll <- rep(1000,nverts(reference))
    distvecAlignAll[pubVerts] <- perVertexErrorMeanAlignAll
    mDAll <<- meshDist(reference,distvec=distvecAlignAll,tol=1,to=max(distvec[pubVerts])+2,plot = T)
    
### ERROR ESTIMATION: realign entire structure TPS
    predSSMarrAllTPS <<- meshlist2array(predSSMTPS)[,,]
    predalignAllTPS <<- procSym(bindArr(predSSMarrAllTPS,origarrAll,along=3),CSinit=F,scale = F)
    perVertexErrorAlignAllTPS <<- sapply(1:n, function(x) x <<- sqrt(rowSums(predalignAllTPS$rotated[pubVerts,,x]-predalignAllTPS$rotated[pubVerts,,x+n])^2))
    perVertexErrorMeanAlignAllTPS <<- rowMeans(perVertexErrorAlignAllTPS)
    ## RMSE
    perVertexErrorAllTPSSQ <<- sapply(1:n, function(x) x <<- (rowSums(predalignAllTPS$rotated[pubVerts,,x]-predalignAllTPS$rotated[pubVerts,,x+n])^2))
    RMSEAllTPS <<- sqrt(mean(perVertexErrorAllTPSSQ))
    PVEAllTPS <<- mean(perVertexErrorAlignAllTPS)
    PVEAllTPS95 <<- quantile(colMeans(perVertexErrorAlignAllTPS),probs = .95)
    ## Distance vector
    distvecAlignAllTPS <<- rep(1000,nverts(reference))
    distvecAlignAllTPS[pubVerts]<<- perVertexErrorMeanAlignAllTPS
    mDAllTPS <<- meshDist(reference,distvec=distvecAlignAllTPS,tol=1,to=max(distvec[pubVerts])+2,plot = F)
    
    SSMError <<- c(PVEPlain,PVETPS,PVEAlign,PVEAlignTPS,PVEAll,PVEAllTPS)
    SSMError95 <<- c(PVEPlain95,PVETPS95,PVEAlign95,PVEAlignTPS95,PVEAll95,PVEAllTPS95)
}


getErrorMargins <- function(x,target,verts=predVertsInfo$respVerts) {
    out <- vcgClostKD(x,target,sign = FALSE)$quality
    out <- out[verts]
    return(list(mean=mean(out),q95=quantile(out,probs = .95)))
    
}


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
predictFossil <- function(fossil.lm,  mymod, mymod.lm, proc, use.lm=NULL, meanMesh=NULL, means=NULL, useweights=T,target.m=NULL,subs=500,evalverts=NULL,getMod=F,uprange=1,realign2mod=FALSE) {
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

    if(realign2mod) {
        mesh2mod <- icp(fossilPredict,DrawMean(mymod),lm1=fossil.lm[use.lm,], lm2=mymod.lm[use.lm,],iterations=50,subsample = 500,rhotol = pi/2,getTransform = T,uprange=uprange)
        fossilPred2 <- PredictSample(mymod,lmModel=mymod.lm[use.lm,],lmDataset=mesh2mod$landmarks,align=FALSE)
      ##  mm <- statismoConstrainModelSafe(mymod,pt=mymod.lm[use.lm,],sample=mesh2mod$landmarks,computeScores=T,ptValueNoise=2,sdmax = 7)
        fossilPred2Back <- icp(applyTransform(fossilPred2,mesh2mod$transform,T),target.m,subsample=2000,iterations=100,rho=pi/2,uprange = uprange)

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
