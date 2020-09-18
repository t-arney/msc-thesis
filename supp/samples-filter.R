
# make a dataframe which includes the filtered list of samples from either the 
# LGM (19 - 23 ka) or the Holocene (0-6 kyr), not including samples from 
# The polar regions or the Mediterranean. Correct d13c for planktics, and for
# LGM sea level for benthics. Group very closely separated sites, average per site

load('../data/qpid.rds')

samples <- samples[!is.na(samples$d13c) ,c("site_name","d13c","age","taxon_id")]
sites   <- sites[ ,c("site_name","lat","lon","water_depth","ocean_basin")]
taxa    <- taxa[ ,c("taxon_id","grouping","genus","species","habitat")]

samples.lgm <- samples[samples$age >= 19 & samples$age <= 23, ] # peterson.etal2014 operational definition
samples.hol <- samples[samples$age > 0.25 & samples$age <= 6, ]
samples.lgm$period <- "LGM"
samples.hol$period <- "Holocene"

smp <- rbind(samples.hol, samples.lgm)

# merge in the site location info and taxon info
smp <- merge(x = smp, y = sites, all.x = T, all.y = F)
smp <- merge(x = smp, y = taxa, all.x = T, all.y = F)

# create a new column for the major (parent) ocean (i.e. strip modifiers)
smp <- cbind(smp, ocean = factor(gsub("[NEWTS](?!o)", "", smp$ocean_basin, perl = T)))

# remove unreliable polar regions and Mediterranean
smp <- smp[!smp$ocean %in% c("Ar","Med", "So"), ]

# merge in assumed calcification depth of planktics from literature
pl_depths <- read.csv('pl_depths.csv')
smp <- merge(x = smp, y = pl_depths, all.x = T, all.y = F)

# save version before corrections and aggregation
lgm_hol_samples <- smp
save(lgm_hol_samples, file = "lgm_hol_samples.rds")




# corrections (d13c Pl, & LGM SL Bn) --------------------------------------


# select G ruber (not pink), apply correction value from spero.etal2003
ruber_corr <- smp[smp$taxon_id == "G_ruber_w" | smp$taxon_id == "G_ruber", ]
ruber_corr$d13c_corr <- ruber_corr$d13c + 0.94 
ruber_corr$grouping <- "G. ruber"

# select reliable epibenthics
cibs <- smp[smp$genus %in% c("Cibicides", "Cibicidoides", "Planulina"), ]
cibs$grouping <- "Cib. spp."

# correct for LGM shallower sea level (125 +/- 5 m BPSL @murray-wallace.woodroffe2014)
cibs$calc_depth <- cibs$water_depth
cibs[cibs$period == "LGM", "calc_depth"]         <- cibs[cibs$period == "LGM", "water_depth"] - 125 
cibs[cibs$period == "LGM", "calc_depth_shallow"] <- cibs[cibs$period == "LGM", "calc_depth"]  - 5 
cibs[cibs$period == "LGM", "calc_depth_deep"]    <- cibs[cibs$period == "LGM", "calc_depth"]  + 5 
cibs$d13c_corr <- cibs$d13c # correction factor of zero


smp_corr <- rbind(ruber_corr, cibs)


# remove deep benthics
smp_corr <- smp_corr[smp_corr$calc_depth <= 1000, ]



# Group by proximity ------------------------------------------------------

# assign each site to a 10x10deg grid cell
smp_corr$latlon <- paste(floor(smp_corr$lat/10)*10, floor(smp_corr$lon/10)*10, sep = "_")


# standard error on the mean function for aggregation
std.err <- function (x) {
  x <- if (is.vector(x) || is.factor(x)) x else as.double(x)
  x <- na.omit(x)
  sd(x)/sqrt(length(x))
}

# take the mean for each site, planktics and benthics
smp_corr_agg <- collapse::collap(X = smp_corr,
     by = ~ site_name + lat + lon + ocean + ocean_basin + period + grouping,
     custom = list(mean = c("d13c_corr", "calc_depth", "calc_depth_shallow", "calc_depth_deep"),
     std.err = "d13c_corr"))

# take the mean for each latlon "cell", planktics and benthics
smp_corr_prox <- collapse::collap(X = smp_corr,
     by = ~ latlon + ocean + ocean_basin + period + grouping,
     custom = list(mean = c("d13c_corr", "calc_depth", "calc_depth_shallow", "calc_depth_deep"),
     std.err = "d13c_corr"))



# Save to file ------------------------------------------------------------


# LGM and Holocene G. ruber (not pink, corrected d13c) and 
# Cib. spp. (corrected for LGM SL), average per site
save(smp_corr_agg, file = "smp_corr_agg.rds")

# LGM and Holocene G. ruber (not pink, corrected d13c) and 
# Cib. spp. (corrected for LGM SL), average per 10deg lat/lon cell
save(smp_corr_prox, file = "smp_corr_prox.rds")
