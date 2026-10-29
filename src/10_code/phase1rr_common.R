suppressPackageStartupMessages(library(edgeR));suppressPackageStartupMessages(library(data.table));suppressPackageStartupMessages(library(jsonlite))
setDTthreads(1)
R<-Sys.getenv('PHASE1RR_ROOT');O<-file.path(R,'04_screen_reconstruction/phase1rr');Q<-file.path(R,'06_qc/phase1rr');T<-file.path(R,'07_tables/phase1rr')
write_tsv<-function(x,p){if(grepl('\\.gz$',p)){con<-gzfile(p,'wt');write.table(x,con,sep='\t',row.names=FALSE,quote=FALSE,na='NA');close(con)}else fwrite(x,p,sep='\t',na='NA')}
make_design<-function(labels){stopifnot(length(labels)>=5,length(unique(labels))==2);model.matrix(~factor(labels,levels=c(0,1)))}
fit_dataset<-function(cts,keep,gene,labels=c(0,0,0,1,1,1),full=FALSE){
 stopifnot(all(is.finite(cts)),all(cts>=0),all(cts==floor(cts)),length(keep)==nrow(cts),ncol(cts)==length(labels))
 des<-make_design(labels);d<-DGEList(counts=cts,group=labels);d<-d[keep,,keep.lib.sizes=TRUE];d<-calcNormFactors(d,method='TMM');d<-estimateDisp(d,des,robust=TRUE);fit<-glmQLFit(d,des,robust=TRUE,legacy=FALSE);test<-glmQLFTest(fit,coef=2);tab<-test$table;std<-sign(tab$logFC)*sqrt(pmax(0,tab$F));std[tab$logFC==0 & tab$F==0]<-0
 guide<-data.table(mean_baseline_raw=rowMeans(cts[keep,which(labels==0),drop=FALSE]),guide=rownames(cts)[keep],gene=gene[keep],log2FC=tab$logFC,SE_QL_equivalent=ifelse(tab$F>0,abs(tab$logFC)/sqrt(tab$F),NA_real_),QL_F=tab$F,P=tab$PValue,standardized_effect=std,logCPM=tab$logCPM,NB_tagwise=d$tagwise.dispersion,QL_var_post=fit$var.post)
 gs<-guide[,.(PRIMARY_GENE_SCORE=median(standardized_effect),median_log2FC=median(log2FC),trimmed_mean_log2FC=mean(log2FC,trim=.2),guide_sign_concordance=max(mean(log2FC>0),mean(log2FC<0)),n_guides=.N,mean_baseline_raw=mean(mean_baseline_raw)),by=gene]
 gs<-gs[n_guides>=3 & !gene %chin% c('LacZ','EGFP','luciferase')];setorder(gs,gene)
 if(any(!is.finite(gs$PRIMARY_GENE_SCORE)))stop('Nonfinite gene score')
 if(full)return(list(guide=guide,gene=gs,d=d,fit=fit,design=des))
 list(gene=gs,common_dispersion=d$common.dispersion)
}
read_contrast<-function(cid){spec<-fromJSON(file.path(R,'01_analysis_freeze/CONTRAST_SPECIFICATIONS.json'));s<-spec[spec$contrast==cid,];x<-fread(file.path(R,s$raw));cts<-as.matrix(x[,-(1:2)]);rownames(cts)<-x$sgRNA;stopifnot(identical(colnames(cts),unlist(s$samples)),ncol(cts)==6);keep<-rowSums(cts)>=60 & rowSums(cts>=10)>=3 & rowSums(cts[,1:3])>=30;list(spec=s,counts=cts,gene=x$Gene,keep=keep)}
gen_counts<-function(st,mu=st$mu0,copula=FALSE){
 if(!copula)return(matrix(rnbinom(length(mu),mu=as.vector(mu),size=rep(1/st$phi,ncol(mu))),nrow(mu),dimnames=dimnames(mu)))
 z<-matrix(rnorm(length(mu)),nrow(mu));u<-matrix(rnorm(length(st$gene_levels)*ncol(mu)),length(st$gene_levels));z<-sqrt(st$rho)*u[st$gene_id,,drop=FALSE]+sqrt(1-st$rho)*z;p<-pnorm(z);matrix(qnbinom(p,mu=as.vector(mu),size=rep(1/st$phi,ncol(mu))),nrow(mu),dimnames=dimnames(mu))
}
noise_metrics<-function(cts,st){
 # Center within each condition so gene/guide condition effects do not count as culture noise.
 y<-sweep(cts,2,st$raw_depth/mean(st$raw_depth),'/');m<-cbind(matrix(rowMeans(y[,1:3]),nrow(y),3),matrix(rowMeans(y[,4:6]),nrow(y),3));v<-(rowSums((y-m)^2))/4
 rs<-(y-m)/sqrt(pmax(m+st$phi*m^2,1));rs<-rs[st$keep,,drop=FALSE];eligible<-st$gene[st$keep] %in% st$eligible_genes;rs<-rs[eligible,,drop=FALSE];gg<-st$gene[st$keep][eligible];sums<-rowsum(rs,gg,reorder=TRUE);sq<-rowsum(rs^2,gg,reorder=TRUE);n<-as.numeric(table(factor(gg,levels=rownames(sums))));gm<-sums/n
 coherent<-median(rowMeans(gm^2));pair<-sum((sums^2-sq)/(n*(n-1)))/sum(sq/n)
 list(var_quantile=quantile(v[st$keep],c(.25,.5,.75,.9,.99)),var_bins=tapply(v[st$keep],st$ab_bin,median),coherent=coherent,pair=pair,zero=colMeans(cts==0))
}
