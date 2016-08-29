## load required packages
require(limma);require(edgeR); require(EMMREML);require(RColorBrewer);require(ggplot2)
require(gridExtra); require(grid); require(doParallel)

#Load in data.
load("~/Dropbox/SGE_data/MANUSCRIPT_FEB2016/Science submission/revisions/github_code/macaque_status_cell.RData")
# There are 5 R objects
# read_counts:  raw read counts for all 440 samples after removal of genes with a median RPKM < 2
dim(read_counts)
# read_count_per_cell: list of raw read counts for each cell type as its own dataframe in a list. Genes with a median RPKM < 2 within a cell type have been removed from the data frame for that cell type:
names(read_count_per_cell)
lapply(read_count_per_cell,dim)
# info: metadata for all 440 samples
dim(info)
# kin: kinship matrix for all 45 individuals
dim(kin)
# kin_per_cell: kinship matrices for each cell type. This is necessary because there is a different number of samples in each cell type because not all samples passed QC.
names(kin_per_cell)
lapply(kin_per_cell,dim)


## make the cell specific metadata and order it ("info_cell")
info_cell=lapply(unique(info$cell),function(x){return(info[as.character(rownames(info)[info$cell==x]),])})
names(info_cell)=unique(info$cell)
info_cell=lapply(info_cell,function(x){tmp=x;rownames(tmp)=tmp$sample; tmp=tmp[order(rownames(tmp)),];return(tmp)})

## make the Z matrix for EMMA for each cell type:
Z_matrix=lapply(names(info_cell),function(x){mat=matrix(nrow=nrow(info_cell[[x]]),ncol=length(unique(info_cell[[x]][,"ID"]))); colnames(mat)=unique(info_cell[[x]][,"ID"]); rownames(mat)=info_cell[[x]][,"ID"]
for (r in 1:nrow(mat)) {
  for (c in 1:ncol(mat)) {
    if (rownames(mat)[r]==colnames(mat)[c]) {mat[r,c]=1}else {mat[r,c]=0}} }
return(mat)})
names(Z_matrix)=c("CD4","CD8","CD14","CD16","CD20")

#################################################################################################
###                         PCA of cell specific expression
#################################################################################################

## voom normalize and remove group effects using limma
design <- model.matrix(~0+info$group)
dge <- DGEList(counts=read_counts)
dge <- calcNormFactors(dge)
v <- voom(dge,design,plot=FALSE); rm(dge)
fit <-lmFit(v,design); rm(design)
residuals=residuals.MArrayLM(object=fit, v);rm(v);rm(fit)
pca<-prcomp(cor(residuals))
rm(residuals)
qplot(pca$x[,1],pca$x[,2],col=as.factor(info[rownames(pca$x),"cell"]))+xlab(paste("PC1","\nPVE =",round(summary(pca)$imp[2,1],2)))+ylab(paste("PC2","\nPVE =",round(summary(pca)$imp[2,2],2)))+scale_colour_brewer("Cell",palette="Set1")

#################################################################################################
###                         model effects of rank (elo) in each cell
#################################################################################################

## Voom normalize and remove group ('batch') effects using limma: 
residuals=lapply(names(info_cell),function(x){
  design <- model.matrix(~0+info_cell[[x]]$group)
  dge <- DGEList(counts=read_count_per_cell[[x]])
  dge <- calcNormFactors(dge)
  v <- voom(dge,design,plot=FALSE)
  fit <-lmFit(v,design)
  return(residuals.MArrayLM(object=fit, v))
});names(residuals)=names(info_cell)

## Model the effects of rank (elo) and age using EMMA 
EMMA_elo=lapply(c("CD4","CD8","CD14","CD16","CD20"),function(x){
tmp=t(apply(residuals[[x]],1,function(y){
  design=model.matrix(~scale(info_cell[[x]]$elo)+info_cell[[x]]$age)
  K=as.matrix(kin_per_cell[[x]])
  Z=as.matrix(Z_matrix[[x]])
  emma=emmreml(y=y,X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
  p=emma$pvalbeta ## store p-values
  varb=emma$varbetahat ## store variance of BLUPs
  b=emma$betahat ## store BLUPs
  return(c(b,varb,p[,"none"])) ## return results for each gene 
})); colnames(tmp)=c('β_intercept','β_elo','β_age','var_intercept','var_elo','var_age','pval_intercept','pval_elo','pval_age'); return(tmp)}); names(EMMA_elo)=c("CD4","CD8","CD14","CD16","CD20")

## plot histograms of rank effects
par(mfrow=c(2,3))
lapply(names(EMMA_elo),function(x){
ps=EMMA_elo[[x]][,"pval_elo"]
hist(ps,col=alpha("grey",0.5),breaks=100,freq=F,main=x,xlim=c(0,0.5))
})


## Permute the effect of rank on multiple cores using the "parallel" package
## NOTE 1: This can take a long time, so it is advisable to use 
clus <- makeCluster(32)
library(parallel); library(doParallel)
registerDoParallel(cores=32)  

EMMA_permute_rank_pvalues=lapply(names(residuals),function(pp){
  nperm=1000
  for (sim in 1:nperm){
    cell=pp
    clusterExport(clus,c("residuals","info_cell","Z_matrix","kin_per_cell"))
    elo=sample(scale(info_cell[[cell]]$elo)) # permute ranks. NOTE 2: To permute the effects of age or weight, adjust this line to permute age or weight instead of rank (elo)
    tmp=parApply(clus,residuals[[cell]],1,function(y){
      library(EMMREML)
      design=model.matrix(~elo+info_cell[[cell]]$age) ## if permuting age or weight, adjust the design matrix accordingly.
      K=as.matrix(kin_per_cell[[cell]])
      Z=as.matrix(Z_matrix[[cell]])
      emma=emmreml(y=y,X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
      return(emma$pvalbeta[,"none"])}) ## store p-values: column 1 = intercept_pvalue, column 2 = permuted_rank_pvalue, column 3 = age_pvalue
    print(sim)
    if (sim==1) final=t(tmp)
    else final=rbind(final,t(tmp))
  }
  colnames(final)=c("intercept_pvalue","permuted_rank_pvalue","age_pvalue")
  return(final)
}
)
names(EMMA_permute_rank_pvalues)=names(residuals)

## calculate permuted FDRs


#################################################################################################
###                       rank effect in each phase (nested models)
#################################################################################################

EMMA_nested_elo=lapply(c("CD4","CD8","CD14","CD16","CD20"),function(x){
  tmp=t(apply(residuals[[x]],1,function(y){
    design=model.matrix(~info_cell[[x]]$phase:scale(info_cell[[x]]$elo)+info_cell[[x]]$age) 
    K=as.matrix(kin_per_cell[[x]])
    Z=as.matrix(Z_matrix[[x]])
    emma=emmreml(y=y,X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
    p=emma$pvalbeta
    varb=emma$varbetahat
    b=emma$betahat
    return(c(b,varb,p[,"none"]))
  }))
  colnames(tmp)=c('β_intercept','β_age','β_elo_SGE1','β_elo_SGE2','var_intercept','var_age','var_elo_SGE1','var_elo_SGE2','pval_intercept','pval_age','pval_elo_SGE1','pval_elo_SGE2')
  return(tmp)}); names(EMMA_nested_elo)=c("CD4","CD8","CD14","CD16","CD20")

## Plot histograms of p-value distributions
par(mfcol=c(2,5))
lapply(names(EMMA_nested_elo),function(x){
  ps=EMMA_nested_elo[[x]][,"pval_elo_SGE1"]
  hist(ps,col=alpha("grey",0.5),breaks=100,freq=F,main=paste0(x,"\n Phase 1"),xlim=c(0,0.5))
  ps=EMMA_nested_elo[[x]][,"pval_elo_SGE2"]
  hist(ps,col=alpha("grey",0.5),breaks=100,freq=F,main=paste0(x,"\n Phase 2"),xlim=c(0,0.5))
})

#################################################################################################
###                  mediation effects of received aggression and grooming
#################################################################################################


## First, run the models with either group mean centered aggression or group mean centered grooming 
## grooming model
EMMA_groom=lapply(c("CD4","CD8","CD14","CD16","CD20"),function(x){
tmp=t(apply(residuals[[x]],1,function(y){
  design=model.matrix(~scale(info_cell[[x]]$elo)+info_cell[[x]]$age+info_cell[[x]]$mean_centered_grooming)
  K=as.matrix(kin_per_cell[[x]])
  Z=as.matrix(Z_matrix[[x]])
  emma=emmreml(y=y,X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
  p=emma$pvalbeta
  varb=emma$varbetahat
  b=emma$betahat
  return(c(b,varb,p[,"none"]))
})); colnames(tmp)=c('β_intercept','β_elo','β_age','β_groom','var_intercept','var_elo','var_age','var_groom','pval_intercept','pval_elo','pval_age','pval_groom'); return(tmp)}); names(EMMA_groom)=c("CD4","CD8","CD14","CD16","CD20")


EMMA_aggression=lapply(c("CD4","CD8","CD14","CD16","CD20"),function(x){
  tmp=t(apply(residuals[[x]],1,function(y){
    design=model.matrix(~scale(info_cell[[x]]$elo)+info_cell[[x]]$age+info_cell[[x]]$mean_centered_aggression)
    K=as.matrix(kin_per_cell[[x]])
    Z=as.matrix(Z_matrix[[x]])
    emma=emmreml(y=y,X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
    p=emma$pvalbeta
    varb=emma$varbetahat
    b=emma$betahat
    return(c(b,varb,p[,"none"]))
  })); colnames(tmp)=c('β_intercept','β_elo','β_age','β_aggression','var_intercept','var_elo','var_age','var_aggression','pval_intercept','pval_elo','pval_age','pval_aggression'); return(tmp)}); names(EMMA_aggression)=c("CD4","CD8","CD14","CD16","CD20")

## Next, bootstrap the data to calculate confidence intervals of the indirect effect
clus <- makeCluster(16)
registerDoParallel(cores=16)  

niter=1000  ## NOTE: This will take a long time. We found that it was best to run these in groups of 10 and to parallelize across mutliple machines with multiple nodes each. Make sure that you set a different seed, using set.seed(), for each new job. Otherwise you will end up with the same data from each run...
system.time(for (t in 1:niter){
  if (t==1) {bootstrap_results1=lapply(names(info_cell), function(cell){
    C=info_cell[[cell]] 
    n=sample(rownames(C),replace=T) ## sample with replacement
    C=C[n,]  ## covariate matrix
    E=residuals[[cell]][,n] ## expression matrix
    Zmat=matrix(nrow=nrow(C),ncol=length(unique(C$ID)))
    colnames(Zmat)=unique(C$ID)
    rownames(Zmat)=C$ID
    for (r in 1:nrow(Zmat)) {
      for (c in 1:ncol(Zmat)) {
        if (rownames(Zmat)[r]==colnames(Zmat)[c]) {Zmat[r,c]=1}else {Zmat[r,c]=0}} }
    K=kin[unique(C$ID),unique(C$ID)]
    ## run models
    clusterExport(clus,varlist=c("C","E","Zmat","K"),envir=environment())
    ## model 1: age + rank
    m1=t(parApply(clus,E,1,function(y){library(EMMREML);emma=emmreml(y=y,X=model.matrix(~C$age+scale(C$elo)),Z=as.matrix(Zmat),K=as.matrix(K),varbetahat=T); return(c(emma$betahat[3],emma$varbetahat[3]))}))
    ## model 2: age + rank + grooming
    m2=t(parApply(clus,E,1,function(y){library(EMMREML);emma=emmreml(y=y,X=model.matrix(~C$age+scale(C$elo)+C$mean_centered_grooming),Z=as.matrix(Zmat),K=as.matrix(K),varbetahat=T);return(c(emma$betahat[3:4],emma$varbetahat[3:4]))}))
    ## model 3: age + rank + aggression
    m3=t(parApply(clus,E,1,function(y){library(EMMREML);emma=emmreml(y=y,X=model.matrix(~C$age+scale(C$elo)+C$mean_centered_aggression),Z=as.matrix(Zmat),K=as.matrix(K),varbetahat=T);return(c(emma$betahat[3:4],emma$varbetahat[3:4]))}))
    ## model a1: grooming as a function of age and rank
    a1=emmreml(y=C$mean_centered_grooming,X=model.matrix(~scale(C$elo)+C$age),Z=as.matrix(Zmat),K=as.matrix(K))$betahat[2]
    var_a1=emmreml(y=C$mean_centered_grooming,X=model.matrix(~scale(C$elo)+C$age),Z=as.matrix(Zmat),K=as.matrix(K),varbetahat=T)$varbetahat[2]
    ## model a2: aggression as a function of age and rank
    a2=emmreml(y=C$mean_centered_aggression,X=model.matrix(~scale(C$elo)+C$age),Z=as.matrix(Zmat),K=as.matrix(K))$betahat[2]
    var_a2=emmreml(y=C$mean_centered_aggression,X=model.matrix(~scale(C$elo)+C$age),Z=as.matrix(Zmat),K=as.matrix(K),varbetahat=T)$varbetahat[2]
    ## combine the data into a dataframe that will be returned from the lapply function
    data=cbind(m1,m3,m4,rep(a1,nrow(m1)),rep(a2,nrow(m1)),rep(var_a1,nrow(m1)),rep(var_a2,nrow(m1)))
    colnames(data)=c("c","var_c","c_prime_groom_SM","b1_groom_SM","var_c_prime_groom_SM","var_b1_groom_SM","c_prime_agg_SM","b1_agg_SM","var_c_prime_agg_SM","var_b1_agg_SM","a1","a2","var_a1","var_a2")
    return(data)
  })
  names(bootstrap_results1)=names(info_cell)
  bootstrap_results=bootstrap_results1}
  if (t!=1) {
    bootstrap_results1=lapply(names(info_cell), function(cell){
      C=info_cell[[cell]]
      n=sample(rownames(C),replace=T) ## sample with replacement
      C=C[n,]  ## covariate matrix
      E=residuals[[cell]][,n] ## expression matrix
      Zmat=matrix(nrow=nrow(C),ncol=length(unique(C$ID)))
      colnames(Zmat)=unique(C$ID)
      rownames(Zmat)=C$ID
      for (r in 1:nrow(Zmat)) {
        for (c in 1:ncol(Zmat)) {
          if (rownames(Zmat)[r]==colnames(Zmat)[c]) {Zmat[r,c]=1}else {Zmat[r,c]=0}} }
      K=kin[unique(C$ID),unique(C$ID)]
      ## run models
      clusterExport(clus,varlist=c("C","E","Zmat","K"),envir=environment())
      ## model 1: age + rank
      m1=t(parApply(clus,E,1,function(y){library(EMMREML);emma=emmreml(y=y,X=model.matrix(~C$age+scale(C$elo)),Z=as.matrix(Zmat),K=as.matrix(K),varbetahat=T); return(c(emma$betahat[3],emma$varbetahat[3]))}))
      ## model 2: age + rank + grooming
      m2=t(parApply(clus,E,1,function(y){library(EMMREML);emma=emmreml(y=y,X=model.matrix(~C$age+scale(C$elo)+C$mean_centered_grooming),Z=as.matrix(Zmat),K=as.matrix(K),varbetahat=T);return(c(emma$betahat[3:4],emma$varbetahat[3:4]))}))
      ## model 3: age + rank + aggression
      m3=t(parApply(clus,E,1,function(y){library(EMMREML);emma=emmreml(y=y,X=model.matrix(~C$age+scale(C$elo)+C$mean_centered_aggression),Z=as.matrix(Zmat),K=as.matrix(K),varbetahat=T);return(c(emma$betahat[3:4],emma$varbetahat[3:4]))}))
      ## model a1: grooming as a function of age and rank
      a1=emmreml(y=C$mean_centered_grooming,X=model.matrix(~scale(C$elo)+C$age),Z=as.matrix(Zmat),K=as.matrix(K))$betahat[2]
      var_a1=emmreml(y=C$mean_centered_grooming,X=model.matrix(~scale(C$elo)+C$age),Z=as.matrix(Zmat),K=as.matrix(K),varbetahat=T)$varbetahat[2]
      ## model a2: aggression as a function of age and rank
      a2=emmreml(y=C$mean_centered_aggression,X=model.matrix(~scale(C$elo)+C$age),Z=as.matrix(Zmat),K=as.matrix(K))$betahat[2]
      var_a2=emmreml(y=C$mean_centered_aggression,X=model.matrix(~scale(C$elo)+C$age),Z=as.matrix(Zmat),K=as.matrix(K),varbetahat=T)$varbetahat[2]
      ## combine the data into a dataframe that will be returned from the lapply function
      data=cbind(m1,m3,m4,rep(a1,nrow(m1)),rep(a2,nrow(m1)),rep(var_a1,nrow(m1)),rep(var_a2,nrow(m1)))
      colnames(data)=c("c","var_c","c_prime_groom_SM","b1_groom_SM","var_c_prime_groom_SM","var_b1_groom_SM","c_prime_agg_SM","b1_agg_SM","var_c_prime_agg_SM","var_b1_agg_SM","a1","a2","var_a1","var_a2")
      return(data)
    })
    names(bootstrap_results1)=names(info_cell)
    bootstrap_results2=lapply(names(bootstrap_results),function(cell){cbind(bootstrap_results[[cell]],bootstrap_results1[[cell]])}); names(bootstrap_results2)=names(bootstrap_results);bootstrap_results=bootstrap_results2}
})

### Write the data to a text file for posterity (and to combine if you parallelized across multiple machines)
lapply(names(bootstrap_results), function(x){write.table(as.matrix(bootstrap_results[[x]]),col.names=T,quote=F,row.names=T,file=paste0(x,"_boot_mean_center.SEEDNUMBER.txt"))})  

## calculate the median and 95% CIs of the indirect effects (IDE) of aggression and grooming
## aggression
aggression_IDE=lapply(bootstrap_results,function(x){t(apply(x,1,function(z){
  y=as.numeric(z)
  a=y[seq(1,niter,by=14)]/sqrt(y[seq(2,niter,by=14)]) ## standardized effect of rank in model 1 (age + rank)
  b=y[seq(7,niter,by=14)]/sqrt(y[seq(9,niter,by=14)]) ## standardized effect of rank in model 3 (age + rank + aggression)
  aggression_IDE=a-b
  q1=quantile(aggression_IDE,c(0.025,0.5,0.975))
  return(q1)
}))})

## grooming
grooming_IDE=lapply(bootstrap_results,function(x){t(apply(x,1,function(z){
  y=as.numeric(z)
  a=y[seq(1,niter,by=14)]/sqrt(y[seq(2,niter,by=14)]) ## standardized effect of rank in model 1 (age + rank)
  b=y[seq(3,niter,by=14)]/sqrt(y[seq(5,niter,by=14)]) ## standardized effect of rank in model 2 (age + rank + grooming)
  grooming_IDE=a-b
  q1=quantile(grooming_IDE,c(0.025,0.5,0.975))## 
  return(q1)
}))})

## Count the number of genes that with a bootstrapped 95% CI that does not include 0 (so are significant at p<0.05)
lapply(aggression_IDE,function(x){sum(x[,3]*x[,1]>0)})
lapply(grooming_IDE,function(x){sum(x[,3]*x[,1]>0)})

## NOTE: be sure to filter these to only include those genes that are signficantly affected by rank (as determined from the permutation-based FDR above)
