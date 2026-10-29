args<-commandArgs(TRUE); root<-normalizePath(args[1]);out<-normalizePath(args[2]); cid<-'mTKO_B16_T17';src<-file.path(root,'04_screen_reconstruction/phase1rr')
for(stage in c('TRAIN','HELDOUT','SPIKE')){
 z<-readRDS(file.path(src,paste0(cid,'.',stage,'.rds')))
 for(n in c('score','effect')){con<-file(file.path(out,paste0(stage,'_',n,'.bin')),'wb');writeBin(as.double(z[[n]]),con,size=8);close(con)}
 write.table(data.frame(gene=z$gene),file.path(out,paste0(stage,'_genes.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
 write.table(data.frame(simulation=seq_along(z$seeds),seed=z$seeds),file.path(out,paste0(stage,'_seeds.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
}
comp<-file.path(root,'09_reports/bmc_genomics_v2/comparator')
for(stage in c('HELDOUT','SPIKE')){
 B<-if(stage=='HELDOUT')100 else 60
 for(n in c('P_A','P_B')){
  mat<-sapply(1:B,function(b)readRDS(file.path(comp,paste0(stage,'_',b,'.rds')))[[n]])
  con<-file(file.path(out,paste0(stage,'_',n,'.bin')),'wb');writeBin(as.double(mat),con,size=8);close(con)
 }
}
cat('EXPORT_COMPLETE\n')
