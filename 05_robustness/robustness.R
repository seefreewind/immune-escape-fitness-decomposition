args<-commandArgs(trailingOnly=TRUE);Sys.setenv(PHASE1RR_ROOT=normalizePath(args[1]));source(file.path(Sys.getenv('PHASE1RR_ROOT'),'10_code/phase1rr_common.R'));cid<-args[2]
x<-read_contrast(cid);full<-fit_dataset(x$counts,x$keep,x$gene,full=TRUE);g<-full$gene;g[,sample_direction_fraction:=0];g[,sample_half_effect_fraction:=0]
for(j in 1:6){v<-fit_dataset(x$counts[,-j],x$keep,x$gene,labels=c(0,0,0,1,1,1)[-j])$gene;v<-v[match(g$gene,gene)];g[,sample_direction_fraction:=sample_direction_fraction+(median_log2FC*v$median_log2FC>0)/6];g[,sample_half_effect_fraction:=sample_half_effect_fraction+(abs(v$median_log2FC)>=abs(median_log2FC)/2)/6]}
gg<-full$guide[, {effects<-log2FC; f<-median(effects); del<-sapply(seq_along(effects),function(k)median(effects[-k]));list(guide_direction_fraction=mean(f*del>0),guide_half_effect_fraction=mean(abs(del)>=abs(f)/2))},by=gene]
g<-merge(g,gg,by='gene');g[,sample_robust:=sample_direction_fraction>=1-1e-12 & sample_half_effect_fraction>=1-1e-12];g[,guide_robust:=guide_direction_fraction==1 & guide_half_effect_fraction>=.8]
write_tsv(g,file.path(T,paste0(cid,'.RELEASE_ROBUSTNESS.tsv')))
