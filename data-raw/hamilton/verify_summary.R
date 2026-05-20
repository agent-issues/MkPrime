.libPaths(c("/nobackup/pjjg18/mkp-study/lib", .libPaths()))
for (f in c("/nobackup/pjjg18/mkp-study/summary/mk_t01_r01.rds",
            "/nobackup/pjjg18/mkp-study/summary/mkp_geo_t01_r01.rds")) {
  x <- readRDS(f)
  cat("== ", basename(f), " ==\n", sep = "")
  cat("  arm=", x$arm, " n_samples=", x$n_samples,
      " n_trees_total=", x$n_trees_total,
      " n_trees_kept=", x$n_trees_kept, "\n", sep = "")
  cat("  scalar_means:\n"); print(round(x$scalar_means, 3))
  cat("  length(br_means)=", length(x$br_means),
      " length(kp_means)=", length(x$kp_means),
      " p_mean=", x$p_mean, "\n", sep = "")
  cat("  CID: length=", length(x$cid),
      " mean=", round(mean(x$cid), 3),
      " range=[", round(min(x$cid),3), ",", round(max(x$cid),3), "]\n",
      sep = "")
}
