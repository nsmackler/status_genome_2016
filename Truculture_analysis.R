####################
## Load dependencies
####################
library(ggplot2)
library(limma)
library(edgeR)
library(EMMREML)
library(cobs)

###############################################################################
## Set number of random permutations to get empiric null p-values distributions
###############################################################################
iters=2

#############################################################################################################
## Load input files: reads, metadata and kinship matrixes in three data.frames (reads, cols, K, respectively)
#############################################################################################################

load(url("https://github.com/nsmackler/status_genome_2016/blob/master/Truculture_data.RData?raw=true"))

#################################################################################
## Set up metadata:
## 1)Get tissue composition principal components
## 2)Relevel condition factor so as to put the baseline at negative controls (NC)
## 3)Standradize Elo
## 4)Mean-center age.
## 5)Declare flag variable, instrumental for random reshuffling
#################################################################################

#1)Get the two first principal components of tissue composition data to include them as covariates.
tc_vector=c("TC1","TC2","TC3","TC4","TC5","TC6","TC7","TC8","TC9","TC10","TC11")
#TC1: Granulocytes (PMN)
#TC2: Classical monocytes: (CD14+ CD16-) %
#TC3: CD14+ Activated monocytes: (CD14+ CD16+) %
#TC4: CD14- Activated monocytes: (CD14- CD16+) %
#TC5: CD4+ Helper T cells: (CD3+ CD4+) %
#TC6: CD8+ Cytotoxic T cells: (CD3+ CD8+) %
#TC7: CD4+ CD8+ Double positive T cells: (CD3+ CD4+ CD8+) %
#TC8: CD8- B cells (CD3-, CD20+ CD8-) %
#TC9: CD8+ B cells (CD3-, CD20+ CD8+) %
#TC10: Natural killer T-lymphocytes: (CD3+ CD16+) %
#TC11: Natural killer cells: (CD3- CD16+) %

#Get PCA of these
exp <- cols[,tc_vector]
exp <- data.frame(scale(exp))
pca <-prcomp(exp,scale=T)
cols<-data.frame(cols,PC1=pca$x[,1],PC2=pca$x[,2])

##2)Relevel condition factor so as to put the baseline in control condition (NC)
cols$Condition <- relevel(cols$Condition, ref = "NC")

##3)Standardize Elo scores for animals' hierarchy.
cols$Elo=scale(cols$Elo)

##4)Mean center age.
cols$Age=cols$Age-mean(cols$Age)

##5)Declare flag variable (1 at the first entry per animal, 0 in the rest, instrumental for Elo reshuffling)
cols$flag=0
for(i in 1:length(levels(cols$Animal)))
{
	cols$flag[which(cols$Animal == levels(cols$Animal)[i])[1]]=1
}

##################################
## Remove group batch-like effects
##################################

#Get dge object
dge <- DGEList(counts=reads)
dge <- calcNormFactors(dge)

#Remove group effects using limma and extract residuals.
design = model.matrix(~0+Group,data=cols)
v <- voom(dge,design,plot=FALSE)
fit <-lmFit(v,design)
fit <- eBayes(fit)
exp<-residuals.MArrayLM(object=fit, v)

#order genes alphabetically to ensure a criterium for ordering rows all the time.
exp=exp[order(rownames(exp)),]

###########################################################################################
## Run linear model: nested (for Elo effects in LPS+ and LPS-)
###########################################################################################
title="Elo_and_LPS_effects"

#Declare full, nested design for fixed effects.
design = model.matrix(~Condition+Condition:Elo+Condition:Age+Condition:PC1+Condition:PC2,data=cols)

#Declare object res_full to store in it the model results: beta coefficients, standard deviations and p values
res_full=exp[,1:(3*ncol(design))]
colnames(res_full)[1:ncol(design)]=paste0("beta_",colnames(design))
colnames(res_full)[(ncol(design)+1):(2*ncol(design))]=paste0("sdev_",colnames(design))
colnames(res_full)[((2*ncol(design))+1):(3*ncol(design))]=paste0("p_value_",colnames(design))

#Declare object random_effects to store in it the individual-wise u random effects
random_effects=exp[,1:length(levels(cols$Animal))]
colnames(random_effects)=levels(cols$Animal)

#Define matrix Z describing the sample-to-individual mapping
Z=matrix(rep(0,nrow(cols)*ncol(K)),nrow=nrow(cols),ncol=ncol(K))
rownames(Z)=rownames(cols)
colnames(Z)=colnames(K)
for(i in 1:ncol(Z))
{
	set=which(cols$Animal == colnames(Z)[i])
	Z[set,i]=1
}

#Fit a model for each gene using emmreml
for(i in 1:nrow(exp))
{
	emma=emmreml(y=exp[i,],X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
    random_effects[i,]=t(emma$uhat)
	res_full[i,]=t(c(emma$betahat,emma$varbetahat,emma$pvalbeta[,"none"]))
}

####################################################################################
## Regress out age, tissue composition and individual random effects to get PCA plot
####################################################################################

reduced_expression_matrix <- exp -
res_full[,"beta_ConditionNC:Age"]%*%t(design[,"ConditionNC:Age"]) -
res_full[,"beta_ConditionLPS:Age"]%*%t(design[,"ConditionLPS:Age"]) -
res_full[,"beta_ConditionNC:PC1"]%*%t(design[,"ConditionNC:PC1"]) -
res_full[,"beta_ConditionLPS:PC1"]%*%t(design[,"ConditionLPS:PC1"]) -
res_full[,"beta_ConditionNC:PC2"]%*%t(design[,"ConditionNC:PC2"]) -
res_full[,"beta_ConditionLPS:PC2"]%*%t(design[,"ConditionLPS:PC2"])

design = model.matrix(~0+Animal,data=cols)

for(i in 1:ncol(design))
{
    reduced_expression_matrix = reduced_expression_matrix -random_effects[,i]%*%t(design[,i])
}

##################
## plot PC1 vs PC2
##################

## Retrieve PC's from pearson correlation matrix on the reduced expression matrix
corr_mat <- cor(reduced_expression_matrix,method="pearson")
pca <-prcomp(corr_mat,scale=T)

## Build a data frame with PC1,PC2,condition and ELo.
PCA_table <- data.frame(PC1=pca$x[,1],PC2=pca$x[,2])
metadata=cols[,c("Elo","Condition")]
PCA_table=merge(PCA_table,metadata,by=0)

#Flip PC1 (conventional, just to put LPS at the right side of the plot)
PCA_table$PC1=-PCA_table$PC1

#Plot
p=ggplot(PCA_table,aes(x=PC1, y=PC2)) +
geom_point(data=PCA_table[which(PCA_table$Condition %in% "NC"),], aes(x=PC1, y=PC2, color=Elo), size=3) +
geom_point(data=PCA_table[which(PCA_table$Condition %in% "LPS"),], aes(x=PC1, y=PC2, fill=Elo), color="white",shape=21, size=4) +
scale_color_gradient(name="Elo (NC)",low="deepskyblue", high="blue4",breaks=c(0,900,1800)) +
scale_fill_gradient(name="Elo (LPS)",low="palegreen", high="darkgreen",breaks=c(0,900,1800))+
ylab("PC2 (22.5% variance)")+
xlab("PC1 (63.8% variance)")+
guides(colour = guide_colorbar(barwidth = 2, barheight = 8,title.position="top",title.vjust=1)) +
guides(fill = guide_colorbar(barwidth = 2, barheight = 8,title.position="top",title.vjust=1,text.hjust=0)) +
theme(axis.text.y   = element_text(size=12,color="black"),
axis.text.x   = element_text(size=12,color="black"),
axis.title.x  = element_text(size=12),
axis.title.y  = element_text(size=12),
panel.background = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.line = element_line(colour = "black"),
panel.border = element_rect(colour = "black", fill=NA, size=1),
legend.title=element_text(size=12),
legend.text=element_text(size=12))

#Save plot
pdf(paste0("results/",title,"/PC1_PC2.pdf"),width = 5.5, height = 4)
print(p)
dev.off()

#Save PCA table
write.table(PCA_table,paste0("results/",title,"/PCA_table.txt"))

################################################################################################################################
## Run iters iterations of the model after permutting Elo ratings to retrieve an empiric distribution of p-values for rank effects
################################################################################################################################

cols_random<-cols

for(iter in 1:iters)
{
	print(iter)
	
    #Permute Elo among the set of samples flagged as 1 (one sample per individual)
	cols_random$Elo[which(cols_random$flag==1)]=sample(cols_random$Elo[which(cols_random$flag==1)])
    
    #Complete the permuted Elos of the other samples (this way the mapping individual-to-rank is conserved in the permutations)
	for(i in 1:length(levels(cols_random$Animal)))
	{
		set=which(cols_random$Animal==levels(cols_random$Animal)[i] & cols_random$flag==0)
		set_ref=which(cols_random$Animal==levels(cols_random$Animal)[i] & cols_random$flag==1)
		if(length(set)==1)
		cols_random$Elo[set]=cols_random$Elo[set_ref]
	}	
	
    #Declare null (i.e. based on permutations) nested design for fixed effects.
	design = model.matrix(~Condition+Condition:Elo+Condition:Age+Condition:PC1+Condition:PC2,data=cols_random)
    
    #Declare object res_null to store in it the permutations' p-values:
	res_null=exp[,1:(ncol(design))]
	colnames(res_null)[1:ncol(design)]=paste0("p_value_",colnames(design))
	
    #Fit a model for each gene using emmreml
	for(i in 1:nrow(exp))
	{
		emma=emmreml(y=exp[i,],X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
		res_null[i,]=t(c(emma$pvalbeta[,"none"]))
	}
    
    #we register p-values of the associations to Elo at NC and LPS alone.
	if(iter==1)
	{
		shuffled_elos_pvals_NC <-data.frame(x=res_null[,"p_value_ConditionNC:Elo"])
		shuffled_elos_pvals_LPS <-data.frame(x=res_null[,"p_value_ConditionLPS:Elo"])
		
		rownames(shuffled_elos_pvals_NC)=rownames(res_null)
		rownames(shuffled_elos_pvals_LPS)=rownames(res_null)
	} else {
		shuffled_elos_pvals_NC <- cbind(shuffled_elos_pvals_NC,x=res_null[,"p_value_ConditionNC:Elo"])
		shuffled_elos_pvals_LPS <- cbind(shuffled_elos_pvals_LPS,x=res_null[,"p_value_ConditionLPS:Elo"])
	}
}

###################################################################################################
## Run iters iterations of the model after permuting Condition to retrieve an empirical distribution
## of p-values for LPS stimulation effect (at average rank, age and tissue composition)
###################################################################################################

cols_random<-cols

for(iter in 1:iters)
{
    print(iter)
    
    #Permute Elo among the set of samples flagged as 1 (one sample per individual)
    cols_random$Elo[which(cols_random$flag==1)]=sample(cols_random$Elo[which(cols_random$flag==1)])
    
    #Complete the permuted Elos of the other samples (this way the mapping individual-to-rank is conserved in the permutations.
    for(i in 1:length(levels(cols_random$Animal)))
    {
        set_ref=which(cols_random$Animal==levels(cols_random$Animal)[i])
        if(length(set_ref)>1)
        {
            set_rand=sample(set_ref)
            cols_random$Condition[set_ref]=cols$Condition[set_rand]
        }else{
            if(runif(1)<0.5){
                cols_random$Condition[set_ref]="NC"}else{cols_random$Condition[set_ref]="LPS"
                }
        }
    }

    #Declare null (i.e. based on permutations) nested design for fixed effects.
    design = model.matrix(~Condition+Condition:Elo+Condition:Age+Condition:PC1+Condition:PC2,data=cols_random)
    
    #Declare object res_null to store in it the permutations'  p-values:
    res_null=exp[,1:(ncol(design))]
    colnames(res_null)[1:ncol(design)]=paste0("p_value_",colnames(design))
    
    #Fit a model for each gene using emmreml
    for(i in 1:nrow(exp))
    {
        emma=emmreml(y=exp[i,],X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
        res_null[i,]=t(c(emma$pvalbeta[,"none"]))
    }
    
    #we register p-values of the associations to Elo at NC and LPS alone.
    if(iter==1)
    {
        shuffled_elos_pvals_LPS_effect <-data.frame(x=res_null[,"p_value_ConditionLPS"])
        rownames(shuffled_elos_pvals_LPS_effect)=rownames(res_null)
    } else {
        shuffled_elos_pvals_LPS_effect <- cbind(shuffled_elos_pvals_LPS_effect,x=res_null[,"p_value_ConditionLPS"])
    }
}

###################################################
## Register tables with p-values after permutations
###################################################

# Create subfolder to store the tables (and other outputs)
comando=paste0("mkdir -p results/",title)
system(comando)

# Write permutations' p-values tables.
write.table(shuffled_elos_pvals_NC,paste0(getwd(), "/results/",title,"/permuted_pvalues_Elo_at_NC.txt"))
write.table(shuffled_elos_pvals_LPS,paste0(getwd(), "/results/",title,"/permuted_pvalues_Elo_at_LPS.txt"))
write.table(shuffled_elos_pvals_LPS_effect,paste0(getwd(), "/results/",title,"/permuted_pvalues_LPS_effect.txt"))

############################################################
## Correct for Multiple testing and filter results to report
############################################################

perm.fdr=function(input_df,perm_df,Pvals_col_name,name){
    
    pvals_index=which(colnames(input_df)==Pvals_col_name)
    ro<-input_df[order(input_df[,pvals_index]),]
    p_obs <- data.frame(pvalue=ro[,pvals_index])
    p_vector<-matrix(as.matrix(perm_df),ncol=1)
    p_vector=data.frame(p_vector[order(p_vector)])
    
    F<-p_obs[,1]
    F_o<-p_obs[,1]
    pi_hat<-p_obs[,1]
    
    j=1
    observed<-length(p_obs[,1])
    randoms<-length(p_vector[,1])
    
    for(i in 1:observed)
    {
        repeat
        {
            if((p_vector[j,1]<p_obs[i,1])&j<randoms){j<-j+1}else{break}
        }
        F[i]=i/observed
        F_o[i]=(j-1)/randoms
        if(F_o[i]<1){pi_hat[i]=(1-F[i])/(1-F_o[i])}else{pi_hat[i]=1}
    }
    tabla <-data.frame(pi_hat,pval=p_obs[,1])
    
    tabla[1,]=c(1,0)
    last_percentile_average=mean(tabla$pi_hat[as.integer(min((length(tabla[,1])*0.99),(nrow(tabla)-1)):length(tabla[,1]))])
    tabla[nrow(tabla),]=c(last_percentile_average,1)
    constraint_matrix=as.matrix(data.frame(c(0,2),c(0,1),c(1,0)))
    f_hat<-suppressWarnings(cobs(tabla$pval,tabla$pi_hat,constraint="convex",pointwise=constraint_matrix,maxiter=1000,print.warn=FALSE,print.mesg=FALSE))
    
    f_hat_serie=f_hat$fitted
    pi_o=f_hat_serie[length(f_hat_serie)]
    pi_o=min(pi_o,1)
    pi_o=max(pi_o,0)
    
    Fdr_ST_perm=pi_o*F_o/F
    
    for(i in 1:length(p_obs[,1]))
    {
        Fdr_ST_perm[i]=pi_o*F_o[i]/F[i]
        if(i>1)
        {
            for(j in 1:(i-1))
            {
                if(Fdr_ST_perm[i-j]>Fdr_ST_perm[i]){Fdr_ST_perm[i-j]=Fdr_ST_perm[i]}else{break}
            }
        }
        if(Fdr_ST_perm[i]>1)  Fdr_ST_perm[i]=1
    }
    
    fdrs_df <-data.frame(ro,q_ST_perm=Fdr_ST_perm)
    rownames(fdrs_df)=rownames(ro)
    colnames(fdrs_df)[ncol(fdrs_df)]=paste0("fdr_",name)
    
    return(fdrs_df)
}
res_full=perm.fdr(data.frame(res_full),shuffled_elos_pvals_NC,"p_value_ConditionNC.Elo","Elo_at_NC")
res_full=perm.fdr(data.frame(res_full),shuffled_elos_pvals_LPS,"p_value_ConditionLPS.Elo","Elo_at_LPS")
res_full=perm.fdr(data.frame(res_full),shuffled_elos_pvals_LPS_effect,"p_value_ConditionLPS","LPS_effect")


res_full=res_full[order(rownames(res_full)),c(
"beta_ConditionLPS","beta_ConditionNC.Elo","beta_ConditionLPS.Elo",
"sdev_ConditionLPS","sdev_ConditionNC.Elo","sdev_ConditionLPS.Elo",
"p_value_ConditionLPS","p_value_ConditionNC.Elo","p_value_ConditionLPS.Elo",
"fdr_LPS_effect","fdr_Elo_at_NC","fdr_Elo_at_LPS")]

colnames(res_full)=c(
"beta_LPS_effect","beta_Elo_at_NC","beta_Elo_at_LPS",
"sdev_LPS_effect","sdev_Elo_at_NC","sdev_Elo_at_LPS",
"p_value_LPS_effect","p_value_Elo_at_NC","p_value_Elo_at_LPS",
"fdr_LPS_effect","fdr_Elo_at_NC","fdr_Elo_at_LPS")


###############################
## Write full model results
###############################

#Fixed effects coefficients, standard errors, pvalues and fdrs for Elo effects within condition and LPS stimulation effects:
write.table(res_full,paste0("results/",title,"/results_Elo_and_LPS_effects.txt"))

########################################################
## Run linear model modeling explicit interactions terms
########################################################
title="Interactions_model"

#Re-declare full design for fixed effects. We now model interactions explicitly:
design = model.matrix(~Condition+Condition*Elo+Condition*Age+Condition*PC1+Condition*PC2,data=cols)

#Re-declare object res_full to store in it the model results: beta coefficients, standard deviations and p values
res_full=exp[,1:(3*ncol(design))]
colnames(res_full)[1:ncol(design)]=paste0("beta_",colnames(design))
colnames(res_full)[(ncol(design)+1):(2*ncol(design))]=paste0("sdev_",colnames(design))
colnames(res_full)[((2*ncol(design))+1):(3*ncol(design))]=paste0("p_value_",colnames(design))

#Fit a model for each gene using emmreml
for(i in 1:nrow(exp))
{
    emma=emmreml(y=exp[i,],X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
    res_full[i,]=t(c(emma$betahat,emma$varbetahat,emma$pvalbeta[,"none"]))
}


#################################################################################################################################
## Run iters iterations of the model after permutting Elo scores to retrieve an empiric distribution of p-values for interactions
#################################################################################################################################

cols_random<-cols

for(iter in 1:iters)
{
    print(iter)
    
    #Permute Elo among the set of samples flagged as 1 (one sample per individual)
    cols_random$Elo[which(cols_random$flag==1)]=sample(cols_random$Elo[which(cols_random$flag==1)])
    
    #Complete the poermuted Elos of the other samples (this way the mapping individual-to-rank is conserved in the permutations.
    for(i in 1:length(levels(cols_random$Animal)))
    {
        set=which(cols_random$Animal==levels(cols_random$Animal)[i] & cols_random$flag==0)
        set_ref=which(cols_random$Animal==levels(cols_random$Animal)[i] & cols_random$flag==1)
        if(length(set)==1)
        cols_random$Elo[set]=cols_random$Elo[set_ref]
    }
    
    #Declare null (i.e. based on permutations) nested design for fix effects.
    design = model.matrix(~Condition+Condition*Elo+Condition*Age+Condition*PC1+Condition*PC2,data=cols_random)
    
    #Declare object res_null to store in it the permutations'  p-values:
    res_null=exp[,1:(ncol(design))]
    colnames(res_null)[1:ncol(design)]=paste0("p_value_",colnames(design))
    
    #Fit a model for each gene using emmreml
    for(i in 1:nrow(exp))
    {
        emma=emmreml(y=exp[i,],X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
        res_null[i,]=t(c(emma$pvalbeta[,"none"]))
    }
    
    #we register p-values of the interaction between LPS stimulation and Elo:
    if(iter==1)
    {
        shuffled_elos_pvals_interaction <-data.frame(x=res_null[,"p_value_ConditionLPS:Elo"])
        rownames(shuffled_elos_pvals_interaction)=rownames(res_null)
    } else {
        shuffled_elos_pvals_interaction <- cbind(shuffled_elos_pvals_interaction,x=res_null[,"p_value_ConditionLPS:Elo"])
    }
}

###################################################
## Register tables with p-values after permutations
###################################################

# Create subfolder to store the tables (and other outputs)
comando=paste0("mkdir -p results/",title)
system(comando)

# Write permutations' p-values table.
write.table(shuffled_elos_pvals_interaction,paste0(getwd(), "/results/",title,"/permuted_pvalues_interaction.txt"))

###############################
## Correct for Multiple testing
###############################

res_full=perm.fdr(data.frame(res_full),shuffled_elos_pvals_interaction,"p_value_ConditionLPS.Elo","interaction")
res_full=res_full[order(rownames(res_full)),c("beta_ConditionLPS.Elo","sdev_ConditionLPS.Elo","p_value_ConditionLPS.Elo","fdr_interaction")]
colnames(res_full)=c("beta_interaction_LPSxElo","sdev_interaction_LPSxElo","p_value_interaction_LPSxElo","fdr_interaction_LPSxElo")

###############################
## Write full model results
###############################

#Fixed effects coefficients, standard errors, pvalues and fdrs for rank effects:
write.table(res_full,paste0("results/",title,"/results_interactions.txt"))








