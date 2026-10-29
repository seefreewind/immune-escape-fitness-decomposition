args<-commandArgs(trailingOnly=TRUE);Sys.setenv(PHASE1RR_ROOT=normalizePath(args[1]));source(file.path(Sys.getenv('PHASE1RR_ROOT'),'10_code/phase1rr_common.R'));cid<-args[2];set.seed(20261019+match(cid,c('mTKO_4T1_T25','mTKO_B16_T17','mTKO_MC38_T15'))*100000)
cat('START',cid,format(Sys.time()),'\n');flush.console();xx<-read_contrast(cid);cts<-xx$counts;start<-proc.time()[3];res<-fit_dataset(cts,xx$keep,xx$gene,full=TRUE);cat('Observed fit seconds',proc.time()[3]-start,'\n');flush.console()
write_tsv(res$guide,file.path(O,paste0(cid,'.GUIDE_QL.tsv.gz')));write_tsv(res$gene,file.path(T,paste0(cid,'.OBSERVED_GENE_SCORES.tsv')));write_tsv(data.table(sample=colnames(cts),res$d$samples),file.path(Q,paste0(cid,'.TMM.tsv')))
# Full library generator including low-abundance guides, with dispersion fallback for non-inference rows.
phi<-rep(res$d$common.dispersion,nrow(cts));phi[xx$keep]<-res$d$tagwise.dispersion;phi<-pmax(phi,1e-8);off<-log(colSums(cts)*res$d$samples$norm.factors)
bf<-glmFit(cts[,1:3],matrix(1,3,1),dispersion=phi,offset=matrix(off[1:3],nrow(cts),3,byrow=TRUE),prior.count=0)
w<-rowMeans(sweep(bf$fitted.values,2,exp(off[1:3]),'/'));w<-w/sum(w);mu0<-outer(w,colSums(cts));dimnames(mu0)<-dimnames(cts)
pf<-glmFit(cts,res$design,dispersion=phi,offset=matrix(off,nrow(cts),6,byrow=TRUE),prior.count=0);mup<-sweep(pf$fitted.values,2,colSums(pf$fitted.values),'/');mup<-sweep(mup,2,colSums(cts),'*')
ab<-rowMeans(cts[xx$keep,1:3]);abbin<-pmin(10,ceiling(rank(ab,ties.method='first')/length(ab)*10));st<-list(contrast=cid,counts=cts,keep=xx$keep,gene=xx$gene,raw_depth=colSums(cts),phi=phi,mu0=mu0,mup=mup,observed=res$gene,common_dispersion=res$d$common.dispersion,ab_bin=abbin,eligible_genes=res$gene$gene,gene_levels=unique(xx$gene),gene_id=match(xx$gene,unique(xx$gene)),rho=0,method='INDEPENDENT_NB',seed=.Random.seed)
obs<-noise_metrics(cts,st);rr<-list();nullsum<-matrix(0,nrow(cts),6);preds<-list();ndep<-matrix(NA_real_,100,6);abq<-matrix(NA_real_,100,4);refdisp<-list()
for(b in 1:100){null<-gen_counts(st);pred<-gen_counts(st,st$mup);nullsum<-nullsum+null;ndep[b,]<-colSums(null);abq[b,]<-quantile(rowMeans(null[,1:3]),c(.1,.5,.9,.99));preds[[b]]<-noise_metrics(pred,st)
 if(b<=3){z<-fit_dataset(pred,st$keep,st$gene,full=TRUE);refdisp[[b]]<-list(common=z$d$common.dispersion,bins=tapply(z$d$tagwise.dispersion,st$ab_bin,median))}
}
avg<-function(name)Reduce('+',lapply(preds,'[[',name))/length(preds);rows<-list();add<-function(metric,observed,simulated,criterion,pass,scope='MANDATORY'){rows[[length(rows)+1]]<<-data.table(contrast=cid,generator=st$method,metric=metric,observed=observed,simulated=simulated,criterion=criterion,PASS=pass,scope=scope)}
simmean<-rowMeans(nullsum[,1:3]/100);obsmean<-rowMeans(cts[,1:3]);rho<-cor(obsmean,simmean,method='spearman');add('baseline_mean_spearman',1,rho,'>=.95',rho>=.95)
qobs<-quantile(obsmean,c(.1,.5,.9,.99));qsim<-colMeans(abq);for(k in 1:4)add(paste0('baseline_abundance_q',c(10,50,90,99)[k]),qobs[k],qsim[k],'relative difference<=.30',abs(qsim[k]-qobs[k])/max(qobs[k],1)<=.30)
for(k in 1:6){ratio<-ndep[,k]/st$raw_depth[k];add(paste0('mean_depth_',colnames(cts)[k]),1,mean(ratio),'[.95,1.05]',mean(ratio)>=.95 & mean(ratio)<=1.05);f<-mean(ratio>=.85 & ratio<=1.15);add(paste0('depth_coverage_',colnames(cts)[k]),1,f,'>=.95 within +/-15%',f>=.95);add(paste0('predictive_zero_',colnames(cts)[k]),obs$zero[k],avg('zero')[k],'absolute difference<=.03',abs(obs$zero[k]-avg('zero')[k])<=.03);add(paste0('direct_null_zero_',colnames(cts)[k]),colMeans(cts==0)[k],pnbinom(0,mu=mu0[,k],size=1/phi)|>mean(),'diagnostic, means differ',NA,'DIAGNOSTIC')}
vq<-avg('var_quantile');for(k in 1:5)add(paste0('within_variance_q',c(25,50,75,90,99)[k]),obs$var_quantile[k],vq[k],'sim/observed [.5,2]',vq[k]/obs$var_quantile[k]>=.5 & vq[k]/obs$var_quantile[k]<=2)
r<-cor(obs$var_bins,avg('var_bins'));add('variance_abundance_profile_correlation',1,r,'>=.90',r>=.90)
cd<-mean(sapply(refdisp,'[[','common'));add('refit_common_dispersion',st$common_dispersion,cd,'sim/observed [.5,2]',cd/st$common_dispersion>=.5 & cd/st$common_dispersion<=2)
db<-tapply(phi[st$keep],st$ab_bin,median);sb<-Reduce('+',lapply(refdisp,'[[','bins'))/3;r<-cor(db,sb);add('dispersion_abundance_profile_correlation',1,r,'>=.80',r>=.80)
add('gene_residual_coherence',obs$coherent,avg('coherent'),'observed/sim [.67,1.5]',obs$coherent/avg('coherent')>=.67 & obs$coherent/avg('coherent')<=1.5)
add('within_gene_residual_pair_correlation',obs$pair,avg('pair'),'absolute difference<=.10',abs(obs$pair-avg('pair'))<=.10)
add('guide_identity_number_mapping',nrow(cts),nrow(mu0),'identical',identical(dim(cts),dim(mu0)) & length(st$gene)==nrow(cts))
audit<-rbindlist(rows);marg<-audit[scope=='MANDATORY' & !metric %chin% c('gene_residual_coherence','within_gene_residual_pair_correlation')];st$fidelity_pass<-all(audit[scope=='MANDATORY',PASS]);st$marginal_fidelity_pass<-all(marg$PASS);st$coherence_observed<-obs$coherent;st$coherence_predictive<-avg('coherent');st$pair_observed<-obs$pair;st$pair_predictive<-avg('pair');st$generator_diagnostics<-audit
# Fallback only if needed; estimate pooled residual rank correlation without using any gene labels as positives.
if(!st$fidelity_pass && st$marginal_fidelity_pass){st$fallback_needed<-TRUE}else st$fallback_needed<-FALSE
saveRDS(st,file.path(O,paste0(cid,'.GENERATOR_STATE.rds')));write_tsv(audit,file.path(Q,paste0(cid,'.NULL_FIDELITY.tsv')))
write_tsv(data.table(sim=1:100,ndep),file.path(Q,paste0(cid,'.NULL_SIMULATED_DEPTHS.tsv')));cat('DONE FIDELITY',cid,'PASS',st$fidelity_pass,'marginal',st$marginal_fidelity_pass,'fallback',st$fallback_needed,'\n');print(audit[scope=='MANDATORY' & !PASS]);flush.console()
