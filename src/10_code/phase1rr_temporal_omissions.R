args<-commandArgs(trailingOnly=TRUE);Sys.setenv(PHASE1RR_ROOT=normalizePath(args[1]));source(file.path(Sys.getenv('PHASE1RR_ROOT'),'10_code/phase1rr_common.R'))
for(cid in c('mTKO_B16_T13','mTKO_MC38_T12')){fn<-file.path(T,paste0(cid,'.OBSERVED_GENE_SCORES.tsv'));if(file.exists(fn))next;xx<-read_contrast(cid);fit<-fit_dataset(xx$counts,xx$keep,xx$gene,full=TRUE);write_tsv(fit$gene,fn);write_tsv(fit$guide,file.path(O,paste0(cid,'.GUIDE_QL.tsv.gz')));cat('TEMPORAL',cid,nrow(fit$gene),'\n');flush.console()}
for(cid in c('mTKO_4T1_T25','mTKO_B16_T17','mTKO_MC38_T15')){xx<-read_contrast(cid)
 for(k in 1:6){fn<-file.path(T,paste0(cid,'.OMIT_',colnames(xx$counts)[k],'.GENE_SCORES.tsv'));if(file.exists(fn))next
  labs<-c(0,0,0,1,1,1)[-k];fit<-fit_dataset(xx$counts[,-k],xx$keep,xx$gene,labels=labs,full=TRUE);write_tsv(fit$gene,fn);write_tsv(data.table(contrast=cid,omitted_sample=colnames(xx$counts)[k],sample=colnames(xx$counts)[-k],condition=labs,fit$d$samples),file.path(Q,paste0(cid,'.OMIT_',colnames(xx$counts)[k],'.TMM_DESIGN.tsv')))
 };cat('OMISSIONS',cid,'complete\n');flush.console()
}
