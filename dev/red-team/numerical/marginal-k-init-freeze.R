## AFTER-FIX verification of MARGINAL-K-INIT-001 (init logLik marginal-aware).
## Same 12 seeds as the before run; marginal_k only (sampled_k untouched by the fix).
## Confirm WHOLE-CHAIN unfreeze: per-chain sd for tl/sigma/p all > 0 and means unstick
## from init (tl=0.1*nEdge~3.0, p=0.5, sigma=0.5). Preview rank spread.
setwd("C:/Users/pjjg18/GitHub/worktrees/mkp/marginal-k")
suppressPackageStartupMessages({ pkgload::load_all(getwd(), quiet = TRUE); library(ape); library(TreeTools) })
EXPSTEPS_FIXED <- 1.4; TREE_SHAPE <- 2; K_MAX_PRIOR <- 30L; N_TIP <- 16L; N_CHAR <- 100L
.simTree <- function(nTip){tr<-ape::rtree(nTip,tip.label=paste0("t",seq_len(nTip)));tl<-stats::rgamma(1,TREE_SHAPE,TREE_SHAPE/EXPSTEPS_FIXED);g<-stats::rgamma(nrow(tr$edge),1);tr$edge.length<-tl*(g/sum(g));tr}
.simJCchar <- function(tree,kTrue){nTip<-length(tree$tip.label);states<-integer(2L*nTip-1L);states[nTip+1L]<-sample.int(kTrue,1L)-1L;edges<-tree$edge;el<-tree$edge.length;for(e in seq_len(nrow(edges))){pa<-edges[e,1L];ch<-edges[e,2L];t<-el[e];pSame<-1/kTrue+(1-1/kTrue)*exp(-kTrue*t/(kTrue-1));if(runif(1L)<pSame)states[ch]<-states[pa] else states[ch]<-sample(setdiff(seq.int(0L,kTrue-1L),states[pa]),1L)};states[seq_len(nTip)]}
.canon <- function(v){uv<-sort(unique(v));out<-match(v,uv)-1L;attr(out,"kObs")<-length(uv);out}
mk_one <- function(seed){set.seed(seed);p<-rbeta(1,1,1);rls<-rgamma(1,1,1);tr<-.simTree(N_TIP);tl<-sum(tr$edge.length);u<-rgeom(N_CHAR,p);kT<-pmin(2L+u,K_MAX_PRIOR);m<-matrix(NA_integer_,N_TIP,N_CHAR,dimnames=list(tr$tip.label,NULL));ko<-integer(N_CHAR);for(j in seq_len(N_CHAR)){repeat{cv<-.canon(.simJCchar(tr,kT[j]));if(attr(cv,"kObs")>=2L)break};m[,j]<-cv;ko[j]<-attr(cv,"kObs")};list(tr=tr,mat=m,p_true=p,rls_true=rls,tl_true=tl,kTrue=kT,nedge=nrow(tr$edge))}
run_chain <- function(d){mkd<-MkPrimeData(TreeTools::MatrixToPhyDat(d$mat));st<-d$tr;st$edge.length<-rep_len(0.1,nrow(d$tr$edge));model<-suppressMessages(MkPrimeModel(coding="variable",nCat=1L,kPrimePrior="geometric",likelihoodMode="marginal_k",priorVariant="unconditional",kprimeHyperA=1,kprimeHyperB=1,expSteps=EXPSTEPS_FIXED));mcmc<-MkPrimeMCMC(nIter=12000L,thin=60L,minWarmup=4000L,maxWarmup=4000L,autoTune=FALSE,nRuns=1L,nChains=1L);suppressMessages(suppressWarnings(RunMkPrime(mkd,st,model=model,mcmc=mcmc,fixTopology=TRUE,overwrite=TRUE)))$samples}
seeds <- 20260529L + 0:11
rows <- list()
for (s in seeds) {
  d <- mk_one(s); sm <- tryCatch(run_chain(d), error=function(e) NULL)
  if (is.null(sm)) { cat(sprintf("seed %d ERROR\n", s)); next }
  L<-nrow(sm); tl<-sm[,"tree_length"]; rl<-sm[,"rate_log_sd"]; p<-sm[,"p"]
  rows[[length(rows)+1]] <- data.frame(seed=s, L=L, initTL=0.1*d$nedge,
    tl_sd=sd(tl), tl_mean=mean(tl), tl_true=d$tl_true, tl_rank=sum(tl<d$tl_true),
    rl_sd=sd(rl), rl_mean=mean(rl), rl_true=d$rls_true, rl_rank=sum(rl<d$rls_true),
    p_sd=sd(p),  p_mean=mean(p),   p_true=d$p_true,    p_rank=sum(p<d$p_true))
  cat(sprintf("seed %d L=%d | tl sd=%.3f mean=%.2f(tru %.2f) rk=%d | sig sd=%.3f rk=%d | p sd=%.3f mean=%.2f(tru %.2f) rk=%d\n",
      s,L,sd(tl),mean(tl),d$tl_true,sum(tl<d$tl_true),sd(rl),sum(rl<d$rls_true),sd(p),mean(p),d$p_true,sum(p<d$p_true)))
}
df <- do.call(rbind, rows); saveRDS(df, "C:/Users/pjjg18/AppData/Local/Temp/claude/cache002-after.rds")
L <- df$L[1]; frozen <- df$tl_sd<1e-3 & df$rl_sd<1e-3 & df$p_sd<1e-3
extf <- function(r) mean(r==0L | r==L)
cat(sprintf("\n=== AFTER FIX (marginal_k, n=%d, 12k SBC config) ===\n", nrow(df)))
cat(sprintf("frozen (all-3 sd~0): %d/%d = %.0f%%   (BEFORE fix: 7/12 = 58%%)\n", sum(frozen), nrow(df), 100*mean(frozen)))
cat(sprintf("mean per-chain sd: tl=%.3f  sigma=%.3f (prior sd 1.0)  p=%.3f\n", mean(df$tl_sd), mean(df$rl_sd), mean(df$p_sd)))
cat(sprintf("rank-at-extreme frac: tl=%.0f%%  sigma=%.0f%%  p=%.0f%%   (BEFORE: ~57/57/59%%)\n", 100*extf(df$tl_rank),100*extf(df$rl_rank),100*extf(df$p_rank)))
cat(sprintf("tl_mean/tl_true median=%.2f (1.0=calibrated; init/truth would be >>1)\n", median(df$tl_mean/df$tl_true)))
cat("ranks tl :", paste(df$tl_rank,collapse=" "), "\n")
cat("ranks sig:", paste(df$rl_rank,collapse=" "), "\n")
cat("ranks p  :", paste(df$p_rank,collapse=" "), "\n")
