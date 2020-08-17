#setwd("/netscratch/dep_coupland/grp_schneeberger/projects/AMPRILdenovo/population/candiGeneInPop/AT1G71920-AT5G10330/GWAS")
args <- commandArgs(trailingOnly = TRUE)
inFile <-args[1]
outFile <- args[2]
####################################
#gwas result
gwas<-read.table(inFile,fill=T,stringsAsFactors = F,header = T)
gwas<-gwas[which(gwas$MAF>0.05 & gwas$MAF<0.95),]
dim(gwas)
rownames(gwas)<-gwas$Marker



manhattan <- function(data, colors=c("gray10", "gray50"),ymax="max", limitchromosomes=chr_level,suggestiveline=-log10(5e-06), genomewideline=-log10(5e-8),annotate=NULL, ...) 
{
  d=data
  #if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
  
  #if (any(1:length(limitchromosomes))) d=d[d$CHR %in% limitchromosomes, ]
  d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
  d$logp = -log10(d$P)
  d$pos=NA
  ticks=NULL
  lastbase=0
 
  
  if (ymax=="max") ymax<-ceiling(max(d$logp))
  if (ymax<8) ymax<-8
  
  cat('ymax is:','\t',ymax,'\n')
  
  chr<-unique(as.character(d$CHR))
  numchroms=length(chr)
  colors <- rainbow(numchroms)
  #print(paste(numchroms,chr,colors,sep=' '))
  if (numchroms==1) {
    d$pos=d$BP
    ticks=floor(length(d$pos))/2+1
  } else {
    for (i in 1:numchroms) {
      if (i==1) {
        chr_i<-chr[i]
        d[d$CHR==chr_i, ]$pos<-d[d$CHR==chr_i, ]$BP
      } else {
        chr_i<-chr[i]
        lastbase=lastbase+tail(subset(d,CHR==chr[i-1])$BP, 1)
        d[d$CHR==chr_i, ]$pos<-d[d$CHR==chr_i, ]$BP+lastbase
        print(lastbase)
      }
      ichr<-d[d$CHR==chr_i, ]
      ticks=c(ticks,(ichr$pos[1]+ichr$pos[nrow(ichr)])/2 )
    }
    print(ticks)
  }
  print(tail(d))
  if (numchroms==1) {
    colors<-rainbow(8)
    with(d, plot(pos, logp, ylim=c(0,ymax), xlim=c(d$BP[1],d$BP[nrow(d)]),ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",chr,"position"), col=colors[2],...))
  }else {
    with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab="Chromosome", cex.axis=1.5,cex.lab=1.5,xaxt="n",bty="l", type="n", ...))
    axis(1, at=ticks, lab=chr, cex.axis=2,...)
    icol=1
    for (i in chr) {
      #icol=as.numeric(gsub('chr','',gsub('H','',i)))
      with(d[d$CHR==i, ],points(pos, logp, col=colors[icol], ...))
      print(paste(icol,i,colors[icol],sep=' '))
      icol=icol+1
      #print(head(d[d$CHR==i, ]))
    }
    
  }
  
  if (!is.null(annotate)) {
    d.annotate=d[which(d$SNP %in% annotate), ]
    with(d.annotate, points(pos, logp, col="green3", ...)) 
  }
  
  #if (suggestiveline) abline(h=suggestiveline, col="blue",lty=2)
  #if (genomewideline) abline(h=genomewideline, col="red",lty=6)
}



###########################
#construct a data for plot
pdf(paste(outFile,".manhattan.pdf",sep=""))
par(mar=c(5,5,2,1))
gwas_data<-data.frame(factor(gsub('chr','',gwas$chr),ordered = T),gwas$pos,gwas$P_value)
colnames(gwas_data)<-c('CHR','BP','P')
### manhattan plot
manhattan(gwas_data,pch=20)
abline(h=-log10(0.05/nrow(gwas)),col='black',lty=3,lwd=2) ##threshold
dev.off()

png(paste(outFile,".manhattan.png",sep=""))
par(mar=c(5,5,2,1))
gwas_data<-data.frame(factor(gsub('chr','',gwas$chr),ordered = T),gwas$pos,gwas$P_value)
colnames(gwas_data)<-c('CHR','BP','P')
### manhattan plot
manhattan(gwas_data,pch=20)
abline(h=-log10(0.05/nrow(gwas)),col='black',lty=3,lwd=2) ##threshold
dev.off()

tiff(paste(outFile,".manhattan.tiff",sep=""), units="in", width=5, height=5, res=300)
# insert ggplot code
par(mar=c(5,5,2,1))
gwas_data<-data.frame(factor(gsub('chr','',gwas$chr),ordered = T),gwas$pos,gwas$P_value)
colnames(gwas_data)<-c('CHR','BP','P')
### manhattan plot
manhattan(gwas_data,pch=20)
abline(h=-log10(0.05/nrow(gwas)),col='black',lty=3,lwd=2) ##threshold
dev.off()

res<-gwas
library(qqman)
png(paste(outFile,".qqplot.png",sep=""))
qq(res$P_value, main="Q-Q plot of GWAS p-value",  pch=18, col = "blue4", cex=1.5, las=1)
dev.off()

tiff(paste(outFile,".qqplot.tiff",sep=""),units="in",width=5, heigh=5,res=300)
qq(res$P_value, main="Q-Q plot of GWAS p-value",  pch=18, col = "blue4", cex=1.5, las=1)
dev.off()

pdf(paste(outFile,".qqplot.pdf",sep=""))
qq(res$P_value, main="Q-Q plot of GWAS p-value",  pch=18, col = "blue4", cex=1.5, las=1)
dev.off()
