# Calculate the shallow-deep difference in d13c for LGM and Holocene
# and the difference between the two time periods, with 95% confidence intervals

# load samples filtered to LGM and Holocene ----------------------------------

load('smp_corr_prox.rds')

# uncomment to restrict analysis to one ocean [At,Pa]
# smp_corr_prox <- smp_corr_prox[smp_corr_prox$ocean == "At", ]

# find differences, t-test ------------------------------------------------

# get values for Hol
shallow.vals.hol <- smp_corr_prox[smp_corr_prox$mean.calc_depth <= 50 & smp_corr_prox$period == "Holocene", "mean.d13c_corr"]
deep.vals.hol <- smp_corr_prox[smp_corr_prox$mean.calc_depth >= 500 & smp_corr_prox$period == "Holocene", "mean.d13c_corr"]

# and for LGM
shallow.vals.lgm <- smp_corr_prox[smp_corr_prox$mean.calc_depth <= 50 & smp_corr_prox$period == "LGM", "mean.d13c_corr"]
deep.vals.lgm <- smp_corr_prox[smp_corr_prox$mean.calc_depth >= 500 & smp_corr_prox$period == "LGM", "mean.d13c_corr"]

# Calculate sample sizes for paired data (the size of the smaller dataset)
sample.size.hol <- min(length(shallow.vals.hol), length(deep.vals.hol))
sample.size.lgm <- min(length(shallow.vals.lgm), length(deep.vals.lgm))

my.mean.diff.hol <- vector(mode = "double")
my.se.diff.hol   <- vector(mode = "double")
my.mean.diff.lgm <- vector(mode = "double")
my.se.diff.lgm   <- vector(mode = "double")

# run the sampling lots of times
for (i in 1:1000) {
  # sample the bigger dataset to the size of the smaller for paired data
  shallow.vals.hol_sampled <- shallow.vals.hol[sample(1:length(shallow.vals.hol), sample.size.hol)]
  deep.vals.hol_sampled <- deep.vals.hol[sample(1:length(deep.vals.hol), sample.size.hol)]
  
  # calc diff for Hol
  differences.hol <- shallow.vals.hol_sampled - deep.vals.hol_sampled
  
  # sample the bigger dataset to the size of the smaller for paired data
  shallow.vals.lgm_sampled <- shallow.vals.lgm[sample(1:length(shallow.vals.lgm), sample.size.lgm)]
  deep.vals.lgm_sampled <- deep.vals.lgm[sample(1:length(deep.vals.lgm), sample.size.lgm)]
  
  # calc diff for LGM
  differences.lgm <- shallow.vals.lgm_sampled - deep.vals.lgm_sampled
  
  # add this run's data to running total
  my.mean.diff.hol <- append(my.mean.diff.hol, mean(differences.hol))
  my.se.diff.hol   <- append(my.se.diff.hol, sd(differences.hol)/sqrt(length(differences.hol)))
  my.mean.diff.lgm <- append(my.mean.diff.lgm, mean(differences.lgm))
  my.se.diff.lgm   <- append(my.se.diff.lgm, sd(differences.lgm)/sqrt(length(differences.lgm)))
}

# calc means of each sample run
p.val <- quantile(p.vals, c(.05,.95))
mean.diff.hol <- mean(my.mean.diff.hol )
se.diff.hol   <- mean(my.se.diff.hol   )
mean.diff.lgm <- mean(my.mean.diff.lgm )
se.diff.lgm   <- mean(my.se.diff.lgm   )

print(paste("Hol shallow-deep difference: ", round(mean.diff.hol,1),"±",round(se.diff.hol*1.96,1)," ‰"))
print(paste("LGM shallow-deep difference: ", round(mean.diff.lgm,1),"±",round(se.diff.lgm*1.96,1)," ‰"))
print(paste("Hol-LGM difference: ", round(mean.diff.hol-mean.diff.lgm,1),"±",round(sqrt((se.diff.hol*1.96)^2+(se.diff.lgm*1.96)^2),1)," ‰"))
