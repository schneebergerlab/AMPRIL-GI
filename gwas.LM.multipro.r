# nohup R CMD BATCH gwas_LM.r
args <- commandArgs(trailingOnly = TRUE)


##load genotype and kinship data
rdata <- args[1]
phenoFile <- args[2]
#outdir <- args[3]
outp <- args[3]
ncpu <- as.integer(args[4])

#load("../../../../../data/1001/GWAS/gwas.data.RData")
load(rdata)

###phenotype data
#infile<-'phenotype.GWAS.txt'

phenoData<-read.table(phenoFile,header=T,stringsAsFactors=F)
colnames(phenoData)<-c("Genotype","Phenotype")
#print(head(phenoData))

DATA <- merge(phenoData, genodata, by.x = "Genotype", by.y = "Genotype",sort=T)

ngen<-nrow(DATA)
r_ind<-match(DATA$Genotype,rownames(RD_all))
RD<-RD_all[r_ind,r_ind]

## caculate the coefficient to change the mixed linear model to simple linear model
kinship<-(1-RD)
library(BGLR)
traitY <- DATA$Phenotype
ETA <- list(list(K=kinship,model="RKHS"))
fm <- BGLR(y=traitY,
           ETA=ETA,
           nIter=10000,
           burnIn=1000,
           saveAt=paste(outp,"BGLR_GBLUP_",sep="."),
           verbose=FALSE)

var.gen <- fm$ETA[[1]]$varU
var.err <- fm$varE
I<-diag(1,nrow(RD),ncol(RD))
decomposition<-eigen(kinship)
D<-diag(decomposition$values)

TU<-t(decomposition$vectors)
inv_V<-solve(sqrt(D*var.gen + I*var.err))
coef<-inv_V %*% TU
print(var.gen)
print(var.err)


### Initialize the result
outfile <- paste(outp,"gwas.result.txt",sep=".")
start<-data.frame('Number','Marker','MAF','Estimate','P_value','adjust_r','pos','chr')
write.table(start,file=outfile,quote=F,row.names=F,col.names=F,sep='\t')

### prepare the data for GWAS
onlyDATA <-data.frame(Y=coef %*% traitY,
                      intercept=coef %*% rep(1,ngen))


### function of GWAS
marker_pos<-grep('SNP',colnames(DATA))
GWAS_LM<-function(m){
  marker<-coef %*% DATA[,m]
  snp<-colnames(DATA)[m]
  chr<-allgeno[snp,]$CHR
  pos<-allgeno[snp,]$POS
  maf<-mean(DATA[,m])/2
  subDATA<-data.frame(onlyDATA,marker)
  lm_marker<-lm(Y ~ -1 + intercept + marker, data= subDATA)
  lm_sum<-summary(lm_marker)
  adj_r<-summary(lm_marker)$adj.r.squared
  #res_marker<-c(maker_name,P_alle,P_value,Estimate)
  Pvalue = lm_sum$coefficients['marker',4]
  Estimate = lm_sum$coefficients['marker',1]
  gwas<-c(m,snp,maf,Estimate,Pvalue,adj_r,pos,chr)
  #write.table(gwas,file=sprintf('%s/gwas_%s',outfile,m),quote=F,row.names=F,col.names=F,sep='\t',append=T)
  result<-as.character(gwas)
  return(result)
  cat('Marker:',m,'\n')
}

## run GWAS 
library(parallel)
PARALLELE<-TRUE

N_CPU_CORE <- ncpu 
## N_CPU_CORE the number of CPU used for GWAS
if(PARALLELE){
  res<-mclapply(marker_pos,mc.cores=N_CPU_CORE,FUN=GWAS_LM,mc.preschedule = T)
}else{
  res<-lapply(marker_pos,FUN=GWAS_LM)
}
reuslt_all<-do.call(rbind,res)
colnames(reuslt_all)<-c('Number','Marker','MAF','Estimate','P_value','adjust_r','pos','chr')
write.table(reuslt_all,file=outfile,quote=F,row.names=F,col.names=F,sep='\t', append = TRUE)

#rm(list=ls(all=TRUE))