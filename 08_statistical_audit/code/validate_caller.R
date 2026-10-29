args<-commandArgs(TRUE);root<-normalizePath(args[1]);out<-normalizePath(args[2]);Sys.setenv(PHASE1RR_ROOT=root);source(file.path(root,'10_code/phase1rr_common.R'));cid<-'mTKO_B16_T17'
st<-readRDS(file.path(O,paste0(cid,'.GENERATOR_STATE.rds')));tr<-readRDS(file.path(O,paste0(cid,'.TRAIN.rds')))
sp<-fread(file.path(Q,paste0(cid,'.PRESPECIFIED_SPIKE_GENES.tsv')));f<-setNames(2^sp$true_log2FC,sp$gene)[st$gene];f[is.na(f)]<-1;mu<-st$mu0;mu[,4:6]<-mu[,4:6]*f;mu<-sweep(mu,2,colSums(mu)/st$raw_depth,'/')
audit<-list()
for(stage in c('HELDOUT','SPIKE')){
 ref<-readRDS(file.path(O,paste0(cid,'.',stage,'.rds')));z<-readRDS(file.path(root,'09_reports/bmc_genomics_v2/comparator',paste0(stage,'_1.rds')))
 set.seed(ref$seeds[1]);cts<-gen_counts(st,if(stage=='SPIKE')mu else st$mu0,copula=st$method!='INDEPENDENT_NB');fit<-fit_dataset(cts,st$keep,st$gene,full=TRUE)
 stopifnot(identical(fit$gene$gene,tr$gene))
 pb<-fit$guide[,.(P=min(1,.N*min(P))),by=gene];pb<-pb$P[match(tr$gene,pb$gene)]
 delta<-max(abs(pb-z$P_B));ds<-max(abs(fit$gene$PRIMARY_GENE_SCORE-ref$score[,1]));stopifnot(delta<1e-10,ds<1e-10)
 # Verify BH against an independent direct ordered implementation.
 ord<-order(pb);qs<-pmin(1,rev(cummin(rev(pb[ord]*length(pb)/seq_along(pb)))));q<-numeric(length(pb));q[ord]<-qs
 audit[[stage]]<-data.table(stage=stage,simulation=1,same_genes=TRUE,score_max_error=ds,Bonferroni_minP_max_error=delta,BH_max_error=max(abs(q-p.adjust(pb,'BH'))),count_rows=nrow(cts),guide_rows=nrow(fit$guide),genes=length(pb),seed=ref$seeds[1])
}
fwrite(rbindlist(audit),file.path(out,'CALLER_IMPLEMENTATION_REGENERATION.tsv'),sep='\t');cat('TWO_EXACT_REGENERATIONS_PASS\n')
