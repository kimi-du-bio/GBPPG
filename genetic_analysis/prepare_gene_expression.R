library(edgeR)

dat =read.csv("gene.csv",head=T)
ID = read.table('qc.fam')

idx = match(ID$V1,dat[,1])
dat = dat[idx,]
na.num = apply(dat[,-1],2,function(x){
  sum(x>0.5)
})
sum(na.num > 99) ## 16037
idx = na.num > 99
dat2 = dat[,-1]
dat3 = dat2[,idx]
ID = as.character(dat[,1])
exp = data.frame(taxa=ID,dat3)
write.table(exp,"pheno_exp_16037_ori.txt",col.names = T,row.names = F,quote=F,sep="\t")
exp2 = exp
exp3 = exp2[,-1]

quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)

  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }

  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}
exp4 = quantile_normalisation(exp3)
data = data.frame(FID=exp$taxa,IID=exp$taxa,father=rep(0,295),mother=rep(0,295),sex=rep(0,295),exp4)
write.table(data,"pheno_exp_16037.norm.txt",col.names = T,row.names = F,quote = F,sep="\t")
