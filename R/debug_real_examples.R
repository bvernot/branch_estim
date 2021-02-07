


library(here)
library(cowplot)
source(here('R/estim_branchpoints_fns.R'))
source(here('R/read_data_fns.R'))
source(here('R/read_sims_fns.R'))




####################################################3
##
## now load some DATA
##
####################################################3


sims.dat.archaics <- generate_sims_dat(here('data/sims_og_newrun.txt'), agCols = c('v_gt', 'c_gt', 'a_gt', 'd_gt'))
sims.dat.archaics.no_v <- generate_sims_dat(here('data/sims_og_newrun.txt'), agCols = c('c_gt', 'a_gt', 'd_gt'))
sims.dat.only_v <- generate_sims_dat(here('data/sims_og_newrun.txt'), agCols = c('v_gt'))
# sims.dat.archaics_mh <- generate_sims_dat(here('data/sims_og_newrun.txt'), agCols = c('v_gt', 'c_gt', 'a_gt', 'd_gt', 'mh'))
sims.dat.neand <- generate_sims_dat(here('data/sims_og_newrun.txt'), agCols = c('v_gt', 'c_gt', 'a_gt'))

# Mez1_R5661 Mez2_A9180

dt.sed.analysis.archaics <- read_and_process_genos(here('data/all_simple_gts.mez.tsv.gz'), f_mh.col = 'f_mh.yri',
                                                   agCols = c('v_gt', 'c_gt', 'a_gt', 'd_gt'),
                                                   site.cats = 'all')
dt.sed.analysis.neand.all <- read_and_process_genos(here('data/all_simple_gts.mez.tsv.gz'), f_mh.col = 'f_mh.yri',
                                                agCols = c('v_gt', 'c_gt', 'a_gt'),
                                                site.cats = 'all')
dt.sed.analysis.neand.not_poly_arc <- read_and_process_genos(here('data/all_simple_gts.mez.tsv.gz'), f_mh.col = 'f_mh.yri',
                                                             agCols = c('v_gt', 'c_gt', 'a_gt'),
                                                             site.cats = c('poly_neand', 'mh_seg_arc_fixed0', 'sed_qc_hominin'))
dt.sed.analysis.archaics.swede <- read_and_process_genos(here('data/all_simple_gts.anon_swede_30k.tsv.gz'), f_mh.col = 'f_mh.yri',
                                                             agCols = c('v_gt', 'c_gt', 'a_gt', 'd_gt'),
                                                             site.cats = 'all')

dt.sed.analysis.archaics.lateN <- read_and_process_genos(here('data/all_simple_gts.spy_vind_goy_les.tsv.gz'), 
                                                         f_mh.col = 'f_mh.yri',
                                                         agCols = c('v_gt', 'c_gt', 'a_gt', 'd_gt'),
                                                         site.cats = 'all')
dt.sed.analysis.archaics.no_v.lateN <- read_and_process_genos(here('data/all_simple_gts.spy_vind_goy_les.tsv.gz'), 
                                                         f_mh.col = 'f_mh.yri',
                                                         agCols = c('c_gt', 'a_gt', 'd_gt'),
                                                         site.cats = 'all')
dt.sed.analysis.only_v.lateN <- read_and_process_genos(here('data/all_simple_gts.spy_vind_goy_les.tsv.gz'),
                                                              f_mh.col = 'f_mh.yri',
                                                              agCols = c('v_gt'),
                                                              site.cats = 'all')


dt.sed.analysis.neand.not_poly_arc[lib == 'Mez1_R5661', .N]
dt.sed.analysis.neand.all[lib == 'Mez1_R5661', .N]
dt.sed.analysis.archaics[lib == 'Mez1_R5661', .N]


####################################################3
##
## now run some ANALYSES
##
####################################################3


ret.mez1.archaics <- run_simple_analysis(dt.sed.analysis.archaics[lib == 'Mez1_R5661'],
                                         sims.dat.archaics, max.iter = 30, nbootstraps = 10)
ret.mez1.neand <- run_simple_analysis(dt.sed.analysis.neand.all[lib == 'Mez1_R5661'], 
                                      sims.dat.neand, max.iter = 30, nbootstraps = 10)
ret.mez1.neand.not_poly_arc <- run_simple_analysis(dt.sed.analysis.neand.not_poly_arc[lib == 'Mez1_R5661'], 
                                      sims.dat.neand, max.iter = 30, nbootstraps = 10)

ret.mez2.archaics <- run_simple_analysis(dt.sed.analysis.archaics[lib == 'Mez2_A9180'],
                                         sims.dat.archaics, max.iter = 30, nbootstraps = 10)

dt.sed.analysis.archaics.lateN[, .N, lib]
#                 lib     N
# 1:      Goyet_A9229 58373
# 2: Les_Cottes_A9230 28765
# 3:           A15867  3341
# 4:         SpyR5556 24506
# 5: Vindija_G1_A9348 28786

ret.goyet.archaics <- run_simple_analysis(dt.sed.analysis.archaics.lateN[lib == 'Goyet_A9229'],
                                     sims.dat.archaics, max.iter = 30, nbootstraps = 10)
ret.lescot.archaics <- run_simple_analysis(dt.sed.analysis.archaics.lateN[lib == 'Les_Cottes_A9230'],
                                     sims.dat.archaics, max.iter = 30, nbootstraps = 10)
ret.spy.archaics <- run_simple_analysis(dt.sed.analysis.archaics.lateN[lib == 'SpyR5556'],
                                     sims.dat.archaics, max.iter = 30, nbootstraps = 10)
ret.vind.archaics <- run_simple_analysis(dt.sed.analysis.archaics.lateN[lib == 'Vindija_G1_A9348'],
                                     sims.dat.archaics, max.iter = 30, nbootstraps = 10)
ret.vind.archaics.no_v <- run_simple_analysis(dt.sed.analysis.archaics.no_v.lateN[lib == 'Vindija_G1_A9348'],
                                              sims.dat.archaics.no_v, max.iter = 30, nbootstraps = 10)
ret.vind.only_v <- run_simple_analysis(dt.sed.analysis.only_v.lateN[lib == 'Vindija_G1_A9348'],
                                       sims.dat.only_v, max.iter = 30, nbootstraps = 10)
ret.vind.only_v.hets <- run_simple_analysis(dt.sed.analysis.only_v.lateN[lib == 'Vindija_G1_A9348' & v_gt == 1],
                                       sims.dat.only_v, max.iter = 30, nbootstraps = 10, set.faunal_prop = 0, set.mh_contam = 0)


pct_swede <- .1
dt.sed.analysis.tmp <- dt.sed.analysis.archaics[lib == 'Mez1_R5661']
dt.sed.analysis.tmp[, is.archaic := T]
dt.sed.analysis.tmp.swede <- dt.sed.analysis.archaics.swede[sample(.N, dt.sed.analysis.tmp[, .N*pct_swede])]
# dt.sed.analysis.tmp.swede.rg <- 
dt.sed.analysis.tmp.swede[deam53 == T]$rg <- dt.sed.analysis.tmp[deam53 == T, unique(rg)]
dt.sed.analysis.tmp.swede[deam53 == F]$rg <- dt.sed.analysis.tmp[deam53 == F, unique(rg)]
dt.sed.analysis.tmp.swede[, is.archaic := F]
dt.sed.analysis.tmp <- rbind(dt.sed.analysis.tmp, dt.sed.analysis.tmp.swede)
dt.sed.analysis.tmp[, .(mh_contam = sum(!is.archaic)/.N, .N), rg]

ret.swede.10pct <- run_simple_analysis(dt.sed.analysis.tmp, sims.dat.archaics, max.iter = 30, nbootstraps = 10)
ret.swede.10pct$em.theta
ret.swede.10pct$real.theta <- dt.sed.analysis.tmp[, .(mh_contam = sum(!is.archaic)/.N, .N), rg]


pct_swede <- .1
dt.sed.analysis.tmp <- dt.sed.analysis.archaics[lib == 'Mez1_R5661']
dt.sed.analysis.tmp[, rg := paste0(rg, '___', 1:2)]
dt.sed.analysis.tmp[, is.archaic := T]
dt.sed.analysis.tmp[, .N, rg]
dt.sed.analysis.tmp.swede <- dt.sed.analysis.archaics.swede[sample(.N, dt.sed.analysis.tmp[, .N*pct_swede])]
# dt.sed.analysis.tmp.swede.rg <- 
n.sam <- dt.sed.analysis.tmp.swede[deam53 == T, .N]
dt.sed.analysis.tmp.swede[deam53 == T]$rg <- dt.sed.analysis.tmp[deam53 == T, sample(unique(rg), n.sam, prob = c(.2,.8), replace = T)]
n.sam <- dt.sed.analysis.tmp.swede[deam53 == F, .N]
dt.sed.analysis.tmp.swede[deam53 == F]$rg <- dt.sed.analysis.tmp[deam53 == F, sample(unique(rg), n.sam, prob = c(.2,.8), replace = T)]
dt.sed.analysis.tmp.swede[, is.archaic := F]
dt.sed.analysis.tmp <- rbind(dt.sed.analysis.tmp, dt.sed.analysis.tmp.swede)
dt.sed.analysis.tmp[, .(mh_contam = sum(!is.archaic)/.N, .N), rg]

ret.swede.10pct_split <- run_simple_analysis(dt.sed.analysis.tmp, sims.dat.archaics, max.iter = 30, nbootstraps = 10)
ret.swede.10pct_split$em.theta
ret.swede.10pct_split$real.theta <- dt.sed.analysis.tmp[, .(mh_contam = sum(!is.archaic)/.N, .N), rg]

# merge(ret.swede.10pct_split$em$c$dt.theta, ret.swede.10pct_split$real.theta, suffixes = c('.em', '.real'), by='rg')
ggplot(merge(ret.swede.10pct_split$em$c$dt.theta, ret.swede.10pct_split$real.theta, suffixes = c('.em', '.real'), by='rg'),
       aes(x=mh_contam.real, y=mh_contam.em)) + geom_point() +
  geom_abline(slope=1)
ggplot(merge(ret.swede.10pct$em$c$dt.theta, ret.swede.10pct$real.theta, suffixes = c('.em', '.real'), by='rg'),
       aes(x=mh_contam.real, y=mh_contam.em)) + geom_point() +
  geom_abline(slope=1)



contam.split.iters <- list()
for (c.iter in sprintf("contam_iter_%d", 1:100)) {
  pct_swede <- .1
  dt.sed.analysis.tmp <- dt.sed.analysis.archaics[lib == 'Mez1_R5661']
  dt.sed.analysis.tmp[, rg := paste0(rg, '___', 1:2)]
  dt.sed.analysis.tmp[, is.archaic := T]
  dt.sed.analysis.tmp[, .N, rg]
  dt.sed.analysis.tmp.swede <- dt.sed.analysis.archaics.swede[sample(.N, dt.sed.analysis.tmp[, .N*pct_swede])]
  # dt.sed.analysis.tmp.swede.rg <- 
  n.sam <- dt.sed.analysis.tmp.swede[deam53 == T, .N]
  dt.sed.analysis.tmp.swede[deam53 == T]$rg <- dt.sed.analysis.tmp[deam53 == T, sample(unique(rg), n.sam, prob = c(.2,.8), replace = T)]
  n.sam <- dt.sed.analysis.tmp.swede[deam53 == F, .N]
  dt.sed.analysis.tmp.swede[deam53 == F]$rg <- dt.sed.analysis.tmp[deam53 == F, sample(unique(rg), n.sam, prob = c(.2,.8), replace = T)]
  dt.sed.analysis.tmp.swede[, is.archaic := F]
  dt.sed.analysis.tmp <- rbind(dt.sed.analysis.tmp, dt.sed.analysis.tmp.swede)
  dt.sed.analysis.tmp[, .(mh_contam = sum(!is.archaic)/.N, .N), rg]
  
  ret.swede.10pct_split_iter <- run_simple_analysis(dt.sed.analysis.tmp, sims.dat.archaics, max.iter = 30, nbootstraps = 10)
  ret.swede.10pct_split_iter$em.theta
  ret.swede.10pct_split_iter$real.theta <- dt.sed.analysis.tmp[, .(mh_contam = sum(!is.archaic)/.N, .N), rg]
  contam.split.iters[[c.iter]] <- ret.swede.10pct_split_iter
}

# names(contam.split.iters)
# c.iter <- names(contam.split.iters)[1]
# contam.split.iters[[c.iter]]
contam.split.iters.theta = foreach (c.iter = names(contam.split.iters), .combine=rbind) %do% {
  ret <- contam.split.iters[[c.iter]]
  ret.max <- ret$em$max
  for (b in c('v', 'c', 'anc_1')) {
    if (ret.max$man.max.ll < ret$em[[b]]$man.max.ll) {
      ret.max <- ret$em[[b]]
    }
  }
  cat(ret.max$branchtime)
  dt = merge(ret.max$dt.theta, contam.split.iters[[c.iter]]$real.theta, by='rg', suffixes = c('.em', '.real'))
  dt[rg == 'Mez1_R5661_rg_1_FALSE___1', mh_contam.real := mh_contam.real + 0.02933404] ## this just roughly corrects for mh contam in the real mez1 library
  dt[rg == 'Mez1_R5661_rg_1_FALSE___2', mh_contam.real := mh_contam.real + 0.03207317] ## this just roughly corrects for mh contam in the real mez1 library
  dt[, c.iter := c.iter]
  dt[, branch := ret.max$branch]
  dt[, branchtime := ret.max$branchtime]
  dt[, ll := ret.max$man.max.ll]
}
ggplot(contam.split.iters.theta, aes(x=mh_contam.real, y=mh_contam.em, color=rg)) + geom_point(alpha=.3) + geom_abline(slope=1) +
  ggtitle('Mez1 adding MH contamination [swede] in four read groups') + 
  xlab('Modern human contamination [real]') +
  ylab('Estimated mh contamination') + theme_light() + theme(legend.position = 'none')
ggsave('figures/mez1_contam_from_swede.max.pdf', width=7, height=4)
ggplot(contam.split.iters.theta, aes(y=faunal_prop, x=mh_contam.real, color=rg)) + geom_boxplot() +
  ggtitle('Mez1 adding MH contamination [swede] in four read groups')
# ggsave('figures/mez1_contam_from_swede.faunal.pdf')

ret.mez1.archaics$p1 + ggtitle('Mez1 branch estimates with contam') + theme_classic() + xlab('branch time') + ylab('log-likelihood') + scale_color_discrete('Branch') +
  geom_point(data=contam.split.iters.theta, aes(x=branchtime, color=branch), alpha=.1, y=max(ret.mez1.archaics$grid$ll)) +
  geom_point(y=max(ret.mez1.archaics$grid$ll), x=ret.mez1.archaics$em$c$branchtime, color='red')
ggsave('figures/mez1__ll_plot.contam.pdf', width=7, height=4)



############
## plot low cov neand results
ret.goyet.archaics$em$c$man.max.ll
ret.goyet.archaics$em$anc_1$man.max.ll
ret.goyet.archaics$em$v$man.max.ll
ret.goyet.archaics$em$max$man.max.ll

get_max_deets <- function(ret) {
  ret.max <- ret$em$max
  for (b in c('v', 'c', 'anc_1')) {
    if (ret.max$man.max.ll < ret$em[[b]]$man.max.ll) ret.max <- ret$em[[b]]
  }
  ret.max.deets <- ret.max$dt.theta
  ret.max.deets[, branch := ret.max$branch]
  ret.max.deets[, branchtime := ret.max$branchtime]
  ret.max.deets
}
a <- get_max_deets(ret.goyet.archaics); a; a

low.cov.deets <- rbind(get_max_deets(ret.goyet.archaics),
                       get_max_deets(ret.lescot.archaics),
                       get_max_deets(ret.spy.archaics),
                       get_max_deets(ret.mez1.archaics),
                       get_max_deets(ret.mez2.archaics),
                       get_max_deets(ret.vind.archaics))
low.cov.deets[, lib := tstrsplit(rg, '_')[1]]

ggplot(low.cov.deets[rg %like% 'TRUE'], aes(x=branch, y=branchtime, color=lib)) + geom_point() +
  geom_text_repel(aes(label=lib)) +
  geom_point(data=sims.dat.archaics$branch.bounds[branch %in% c('c', 'v')], aes(y=t.low), pch='-', color='black', size=10) +
  geom_point(data=sims.dat.archaics$branch.bounds[branch %in% c('c', 'v')], aes(y=t.high), pch='-', color='black', size=10) +
  xlab('branch') + ylab('branch time of sample') + ggtitle('branchtime estimates low-coverage Neandertals')
ggsave('figures/low_cov_archaics_branchtime.pdf')

low.cov.deets.flat <- merge(low.cov.deets[rg %like% 'TRUE'], low.cov.deets[rg %like% 'FALSE'], by='lib', suffixes = c('.deam', '.non_deam'))
ggplot(low.cov.deets, aes(x=lib, color=(rg %like% 'TRUE'))) + geom_point(aes(y=mh_contam)) +
  coord_flip() + scale_color_discrete('deam') + theme_classic() +
  ggtitle('contamination estimates in low-coverage Neandertals') +
  geom_linerange(data=low.cov.deets.flat, aes(ymin=mh_contam.deam, ymax=mh_contam.non_deam), color='grey')
ggsave('figures/low_cov_archaics_mhcontam.pdf', width=7, height=4)

low.cov.deets.grid <- rbind(data.table(ret.goyet.archaics$grid.bootstraps, lib='goyet'),
                            data.table(ret.lescot.archaics$grid.bootstraps, lib='les_cot'),
                            data.table(ret.spy.archaics$grid.bootstraps, lib='spy'),
                            data.table(ret.mez1.archaics$grid.bootstraps, lib='mez1'),
                            data.table(ret.mez2.archaics$grid.bootstraps, lib='mez2'),
                            data.table(ret.vind.archaics$grid.bootstraps, lib='vind'))

ggplot(low.cov.deets.grid[is.max == T], aes(x=my.branch, y=my.t, color=lib)) + geom_point() +
  geom_point(data=sims.dat.archaics$branch.bounds[branch %in% c('c', 'v')], aes(y=t.low, x=branch), pch='-', color='black', size=10) +
  geom_point(data=sims.dat.archaics$branch.bounds[branch %in% c('c', 'v')], aes(y=t.high, x=branch), pch='-', color='black', size=10) +
  facet_wrap(~lib) + xlab('branch') + ylab('branch time of sample') + ggtitle('branchtime estimates over 10 bootstraps')
ggsave('figures/low_cov_archaics_mhcontam.pdf')


vind.deets.grid <- rbind(data.table(ret.vind.archaics$grid.bootstraps, lib='vind'),
                            data.table(ret.vind.only_v$grid.bootstraps, lib='vind.only_v'))

ggplot(vind.deets.grid[is.max == T], aes(x=my.branch, y=my.t, color=lib)) + geom_point() +
  geom_point(data=sims.dat.archaics$branch.bounds[branch %in% c('c', 'v')], aes(y=t.low, x=branch), pch='-', color='black', size=10) +
  geom_point(data=sims.dat.archaics$branch.bounds[branch %in% c('c', 'v')], aes(y=t.high, x=branch), pch='-', color='black', size=10) +
  facet_wrap(~lib)


low.cov.deets.grid[, bs.iter := rep(1:(.N/60), each=60)]

ggplot(low.cov.deets.grid[lib == 'goyet' & step.x == 1], aes(x=my.t, y=ll, color=my.branch, lty=lib, 
                                                             group=interaction(bs.iter, lib, my.branch))) + geom_line()
ggplot(low.cov.deets.grid, 
       aes(x=my.t, y=ll, color=my.branch, lty=lib, 
           group=interaction(bs.iter, lib, my.branch))) + geom_line()

ggplot(ret.mez1.archaics$grid, aes(x=my.t, y=ll, color=my.branch)) + geom_line()
ggsave('figures/mez1_archaics_grid_ll.pdf')




### plot ll surfaces
ret.spy.archaics$p1 + ggtitle('Spy') + theme_classic() + xlab('branch time') + ylab('log-likelihood') + scale_color_discrete('Branch')
ggsave('figures/spy_ll_plot.pdf', width=7, height=4)
ret.mez1.archaics$p1 + ggtitle('Mez1') + theme_classic() + xlab('branch time') + ylab('log-likelihood') + scale_color_discrete('Branch')
ggsave('figures/mez1_ll_plot.pdf', width=7, height=4)


### plot ll surfaces in a grid
plot_grid(ret.goyet.archaics$p1 + ggtitle('Goyet') + theme_classic() + xlab('branch time') + ylab('log-likelihood') + scale_color_discrete('Branch'),
          ret.spy.archaics$p1 + ggtitle('Spy') + theme_classic() + xlab('branch time') + ylab('log-likelihood') + scale_color_discrete('Branch'),
          ret.mez1.archaics$p1 + ggtitle('Mez1') + theme_classic() + xlab('branch time') + ylab('log-likelihood') + scale_color_discrete('Branch'),
          ret.mez2.archaics$p1 + ggtitle('Mez2') + theme_classic() + xlab('branch time') + ylab('log-likelihood') + scale_color_discrete('Branch'),
          ncol = 2, align='vh')
ggsave('figures/low_cov__ll_plot.pdf', width=7, height=4)

### plot ll surfaces in a grid w/ bootstraps
plot_grid(ret.goyet.archaics$p1.boots + ggtitle('Goyet') + theme_classic() + xlab('branch time') + ylab('log-likelihood') + scale_color_discrete('Branch'),
          ret.spy.archaics$p1.boots + ggtitle('Spy') + theme_classic() + xlab('branch time') + ylab('log-likelihood') + scale_color_discrete('Branch'),
          ret.mez1.archaics$p1.boots + ggtitle('Mez1') + theme_classic() + xlab('branch time') + ylab('log-likelihood') + scale_color_discrete('Branch'),
          ret.mez2.archaics$p1.boots + ggtitle('Mez2') + theme_classic() + xlab('branch time') + ylab('log-likelihood') + scale_color_discrete('Branch'),
          ncol = 2, align='vh')
ggsave('figures/low_cov__ll_plot.boots.pdf', width=7, height=4)

ret.mez1.archaics$p1.boots + ggtitle('Mez1') + theme_classic() + xlab('branch time') + ylab('log-likelihood') + scale_color_discrete('Branch')
ggsave('figures/mez1__ll_plot.boots.pdf', width=7, height=4)







####################3333
### plot sediment results!?

all_files = list.files(here('analysis/'), pattern='simple_analysis.*.full.RDS')
layer_files = 
  list.files(here('analysis/'), pattern='simple_analysis.*all.full.RDS')

layer.list <- list()
for(this.f in layer_files) {
  tag = gsub('simple_analysis.', replacement = '', this.f)
  tag = gsub('.libs.all.full.RDS', replacement = '', tag)
  cat(tag, '\n')
  layer.list[[tag]] <- readRDS(here(sprintf('analysis/%s', this.f)))
}

libs.list <- list()
plot.list <- list()
for(this.f in all_files) {
  tag = gsub('simple_analysis.', replacement = '', this.f)
  tag = gsub('.full.RDS', replacement = '', tag)
  layer = gsub('.libs.*', replacement = '', tag)
  tag = gsub('.*libs.', replacement = '', tag)
  cat(layer, tag, '\n')
  libs.list[[tag]] <- readRDS(here(sprintf('analysis/%s', this.f)))
  libs.list[[tag]]$layer <- layer
  if (tag == layer) plot.list[[tag]] <- libs.list[[tag]]$p1.boots + ggtitle(tag) + theme_classic() + theme(legend.position = "none") + xlab('branchtime')
}

names(layer.list)
plot_grid(layer.list[['geII_2']]$p1.boots + xlim(.44,1.5) + ylim(-3600,NA) + ggtitle('Pit 2: layer 2') + theme_classic() + theme(legend.position = "none") + xlab('branchtime'),
          layer.list[['geI_2']]$p1.boots + xlim(.44,1.5) + ylim(-8500,NA) + ggtitle('Pit 1: layer 2') + theme_classic() + theme(legend.position = "none") + xlab('branchtime'),
          layer.list[['geI_3']]$p1.boots + xlim(.44,1.5) + ylim(-6000,NA) + ggtitle('Pit 1: layer 3') + theme_classic() + theme(legend.position = "none") + xlab('branchtime'),
          layer.list[['geI_3-4']]$p1.boots + xlim(.44,1.5) + ylim(-390,NA) + ggtitle('Pit 1: layer 3-4') + theme_classic() + theme(legend.position = "none") + xlab('branchtime'),
          layer.list[['sp4424']]$p1.boots + xlim(.44,1.5) + ylim(-975,NA) + ggtitle('Pit 1: layer 4 sp4424') + theme_classic() + theme(legend.position = "none") + xlab('branchtime'),
          layer.list[['geI_4']]$p1.boots + xlim(.44,1.5) + ylim(-720,NA) + ggtitle('Pit 1: layer 4') + theme_classic() + theme(legend.position = "none") + xlab('branchtime'),
          layer.list[['geI_5']]$p1.boots + xlim(.44,1.5) + ylim(-1735,NA) + ggtitle('Pit 1: layer 5') + theme_classic() + theme(legend.position = "none") + xlab('branchtime'),
          ncol = 3, align='vh')
ggsave(here('figures/estatuas_layers_top_libs.pdf'), width=9, height=6)


plot.list <- list()
for (tag in c("D5266", "D5276", 'A15858', 'A15884', 'A15895', 'A15922')) {
  plot.list[[tag]] <- libs.list[[tag]]$p1.boots + 
    ggtitle(tag) + 
    theme_classic() + theme(legend.position = "none") + xlab('branchtime')
}
plot_grid(plotlist = plot.list)
ggsave(here('figures/den_chag_layers_top_libs.pdf'), width=7, height=4)



##########################333
#### FAUNAL CONTAM


contam.split.iters.faunal <- list()
for (c.iter in sprintf("contam_iter_%d", 41:50)) {
  pct_swede <- .1
  dt.sed.analysis.tmp <- dt.sed.analysis.archaics[lib == 'Mez1_R5661']
  dt.sed.analysis.tmp[, rg := paste0(rg, '___', 1:2)]
  dt.sed.analysis.tmp[, is.archaic := T]
  dt.sed.analysis.tmp[, .N, rg]
  dt.sed.analysis.tmp.swede <- dt.sed.analysis.archaics.swede[sample(.N, dt.sed.analysis.tmp[, .N*pct_swede])]
  # dt.sed.analysis.tmp.swede.rg <- 
  n.sam <- dt.sed.analysis.tmp.swede[deam53 == T, .N]
  dt.sed.analysis.tmp.swede[deam53 == T]$rg <- dt.sed.analysis.tmp[deam53 == T, sample(unique(rg), n.sam, prob = c(.2,.8), replace = T)]
  n.sam <- dt.sed.analysis.tmp.swede[deam53 == F, .N]
  dt.sed.analysis.tmp.swede[deam53 == F]$rg <- dt.sed.analysis.tmp[deam53 == F, sample(unique(rg), n.sam, prob = c(.2,.8), replace = T)]
  dt.sed.analysis.tmp.swede[, is.archaic := F]
  dt.sed.analysis.tmp <- rbind(dt.sed.analysis.tmp, dt.sed.analysis.tmp.swede)
  dt.sed.analysis.tmp[, .(mh_contam = sum(!is.archaic)/.N, .N), rg]
  # pct_faunal = runif(1,0,.3)
  pct_faunal = runif(1,.3,.8)
  dt.sed.analysis.tmp[sample(.N, .N*pct_faunal), sed_gt := 0]
  
  ret.swede.10pct_split_iter <- run_simple_analysis(dt.sed.analysis.tmp, sims.dat.archaics, max.iter = 30, nbootstraps = 3)
  ret.swede.10pct_split_iter$em.theta
  ret.swede.10pct_split_iter$real.theta.mh <- dt.sed.analysis.tmp[, .(mh_contam = sum(!is.archaic)/.N, .N), rg]
  ret.swede.10pct_split_iter$real.theta.fp <- dt.sed.analysis.tmp[, .(faunal_prop = pct_faunal, .N), rg]
  contam.split.iters.faunal[[c.iter]] <- ret.swede.10pct_split_iter
}


contam.split.iters.faunal.theta = foreach (c.iter = names(contam.split.iters.faunal), .combine=rbind) %do% {
  ret <- contam.split.iters.faunal[[c.iter]]
  ret.max <- ret$em$max
  ret.max <- ret$em$max
     for (b in c('v', 'c', 'anc_1')) {
         if (ret.max$man.max.ll < ret$em[[b]]$man.max.ll) {
             cat(c.iter, b, ret$em[[b]]$man.max.ll, '\n')
             # ret.max <- ret$em[[b]]
             }
       }
  cat(ret.max$branchtime)
  dt = merge(ret.max$dt.theta, ret$real.theta.fp, by='rg', suffixes = c('.em', '.real'))
  dt = merge(dt, ret$real.theta.mh, by='rg', suffixes = c('.em', '.real'))
  dt[rg == 'Mez1_R5661_rg_1_FALSE___1', mh_contam.real := mh_contam.real + 0.02933404] ## this just roughly corrects for mh contam in the real mez1 library
  dt[rg == 'Mez1_R5661_rg_1_FALSE___2', mh_contam.real := mh_contam.real + 0.03207317] ## this just roughly corrects for mh contam in the real mez1 library
  dt[, c.iter := c.iter]
  dt[, branch := ret.max$branch]
  dt[, branchtime := ret.max$branchtime]
  dt[, ll := ret.max$man.max.ll]
}
ggplot(contam.split.iters.faunal.theta, aes(x=faunal_prop.real, y=faunal_prop.em, color=rg)) + geom_point(alpha=.3) + geom_abline(slope=1) +
  ggtitle('Mez1 adding MH contamination [swede] and Faunal DNA in four read groups') + 
  xlab('Faunal proportion [real]') +
  ylab('Estimated faunal proportion') + theme_light() + theme(legend.position = 'none')
ggsave('figures/mez1_contam_from_swede.fp.max.pdf', width=7, height=4)
ggplot(contam.split.iters.faunal.theta, aes(x=mh_contam.real, y=mh_contam.em, color=faunal_prop.real)) + geom_point(alpha=.5) + geom_abline(slope=1) +
  ggtitle('Mez1 adding MH contamination [swede] and Faunal DNA in four read groups') + 
  xlab('Modern human contamination [real]') +
  ylab('Estimated mh contamination') + theme_light()
ggsave('figures/mez1_contam_from_swede.mh_with_fp.max.pdf', width=7, height=4)
ggplot(contam.split.iters.faunal.theta, aes(x=faunal_prop.real, y=branchtime, color=faunal_prop.real)) + geom_point(alpha=.3) + 
  # geom_abline(slope=1) +
  ggtitle('Mez1 adding MH contamination [swede] in four read groups') + 
  xlab('Faunal proportion [real]') +
  ylab('Estimated branchtime') + theme_light() + theme(legend.position = 'none')
ggsave('figures/mez1_contam_from_swede.fp.branchtime_vs_fp.max.pdf', width=7, height=4)
# ggplot(contam.split.iters.faunal.theta, aes(y=faunal_prop.em, x=mh_contam.real, color=rg)) + geom_boxplot() +
#   ggtitle('Mez1 adding MH contamination [swede] in four read groups')
# ggsave('figures/mez1_contam_from_swede.faunal.pdf')

ret.mez1.archaics$p1 + ggtitle('Mez1 branch estimates with contam') + theme_classic() + xlab('branch time') + ylab('log-likelihood') + scale_color_discrete('Branch') +
  geom_point(data=contam.split.iters.theta, aes(x=branchtime, color=branch), alpha=.1, y=max(ret.mez1.archaics$grid$ll)) +
  geom_point(y=max(ret.mez1.archaics$grid$ll), x=ret.mez1.archaics$em$c$branchtime, color='red')
ggsave('figures/mez1__ll_plot.contam.pdf', width=7, height=4)




library(doParallel)
registerDoParallel(cores=8)
getDoParWorkers()


contam.split.ll <- data.table()
for (pct_swede in c(.1,.3,.5,.8)) {
  cat(pct_swede, '\n')
  dt.sed.analysis.tmp <- dt.sed.analysis.archaics[lib == 'Mez1_R5661']
  dt.sed.analysis.tmp[, rg := paste0(rg, '___', 1:2)]
  dt.sed.analysis.tmp[, is.archaic := T]
  dt.sed.analysis.tmp[, .N, rg]
  dt.sed.analysis.tmp.swede <- dt.sed.analysis.archaics.swede[sample(.N, dt.sed.analysis.tmp[, .N*pct_swede])]
  # dt.sed.analysis.tmp.swede.rg <- 
  n.sam <- dt.sed.analysis.tmp.swede[deam53 == T, .N]
  dt.sed.analysis.tmp.swede[deam53 == T]$rg <- dt.sed.analysis.tmp[deam53 == T, sample(unique(rg), n.sam, prob = c(.2,.8), replace = T)]
  n.sam <- dt.sed.analysis.tmp.swede[deam53 == F, .N]
  dt.sed.analysis.tmp.swede[deam53 == F]$rg <- dt.sed.analysis.tmp[deam53 == F, sample(unique(rg), n.sam, prob = c(.2,.8), replace = T)]
  dt.sed.analysis.tmp.swede[, is.archaic := F]
  dt.sed.analysis.tmp <- rbind(dt.sed.analysis.tmp, dt.sed.analysis.tmp.swede)
  dt.sed.analysis.tmp[, .(mh_contam = sum(!is.archaic)/.N, .N), rg]
  
  tmp.ll <- grid_t_em_theta(dt.sed.analysis.tmp,
                                       sims.dat.archaics,
                                       # my.branches = c('v'), bins.t = 5,
                                       my.branches = c('v','anc_1','c'), bins.t = 5,
                                       max.iter = 2, ll.converge = .01, nsteps = 2, ll.thresh = 10)
  tmp.ll$real.mh_contam <- dt.sed.analysis.tmp[, sum(!is.archaic)/.N]
  tmp.ll$real.nsites <- dt.sed.analysis.tmp[, .N]
  contam.split.ll <- rbind(contam.split.ll, tmp.ll)
  cat('real:', dt.sed.analysis.tmp[, sum(!is.archaic)/.N], '\n')
}

dt.sed.analysis.tmp[, rg := 'hey']
dt.sed.analysis.tmp <- dt.sed.analysis.tmp[sample(.N-1)]
tmp.ll <- grid_t_em_theta(dt.sed.analysis.tmp,
                          sims.dat.archaics,
                          # my.branches = c('v'), bins.t = 5,
                          my.branches = c('v','anc_1','c'), bins.t = 5,
                          max.iter = 2, ll.converge = .01, nsteps = 2, ll.thresh = 10)
tmp.ll$real.mh_contam <- dt.sed.analysis.tmp[, sum(!is.archaic)/.N]
tmp.ll$real.nsites <- dt.sed.analysis.tmp[, .N]
contam.split.ll <- rbind(contam.split.ll, tmp.ll, fill=T)

contam.split.ll[, rel.ll := man.max.ll - max.ll]
## adding more human contam smooths out the ll
## but also adding more rg's tightens it up (i.e., more granularity in the data)
ggplot(contam.split.ll[rg != 'hey'], aes(x=branchtime, y=rel.ll, group=interaction(branch,real.mh_contam), color = paste(real.mh_contam)))  + geom_hline(aes(yintercept = max(rel.ll)-10)) + 
  geom_smooth(se = F)








####
# # I also tried running w/o the sites that are fixed N and poly archaic, didn't make any difference
# dt.sed.analysis.mez1.neand.not_poly_arc.em <- sed_EM_allbranch(dt.sed.analysis.neand.not_poly_arc[lib == 'Mez1_R5661'], 
#                                                       sims.dat.neand, 
#                                                       branches = c('c','v','anc_1'),
#                                                       err_rate = 0.001,
#                                                       max.iter = 30, ll.converge = 1e-6,
#                                                       set.faunal_prop = 'estim',
#                                                       set.mh_contam = 'estim',
#                                                       p_h_method = 'simple')
# 
# dt.sed.analysis.mez1.neand.not_poly_arc.em.theta <- merge(dt.sed.analysis.neand.not_poly_arc[lib == 'Mez1_R5661', .N, rg], dt.sed.analysis.mez1.neand.not_poly_arc.em$max$dt.theta)
# dt.sed.analysis.mez1.neand.not_poly_arc.em.theta[, sum(mh_contam * N / sum(N))]
# 
# dt.sed.analysis.mez1.neand.not_poly_arc.gridll = sed_grid_search(dt.sed.analysis.neand.not_poly_arc[lib == 'Mez1_R5661'], 
#                                                         sims.dat.neand, my.branch = c('v','c','anc_1', 'anc_2', 'anc_3', 'a', 'd'),
#                                                         err_rate = 0.001, p_h_method = 'simple', nsteps = 2,
#                                                         range.mh_contam = dt.sed.analysis.mez1.neand.not_poly_arc.em.theta[, sum(mh_contam * N / sum(N))],
#                                                         range.faunal_prop = dt.sed.analysis.mez1.neand.not_poly_arc.em.theta[, sum(faunal_prop * N / sum(N))])
# 
# p2 <- ggplot(dt.sed.analysis.mez1.neand.not_poly_arc.gridll, aes(x=my.t, y=ll, color=my.branch)) + geom_line() + 
#   geom_point(data=dt.sed.analysis.mez1.neand.not_poly_arc.gridll[is.max == T]) + ggtitle('grid search full likelihood')
# p2
