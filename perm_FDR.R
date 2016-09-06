############################################################
## Function for calculating permutation-based false discovery rates based on a 
## generalization of the false discovery rate method of Storey and Tibshirani (2003)
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
