


library(here)

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


############3
## debug strange results with vindija only_v

# ret.vind.only_v <- run_simple_analysis(dt.sed.analysis.only_v.lateN[lib == 'Vindija_G1_A9348'],
#                                        sims.dat.only_v, max.iter = 30, nbootstraps = 10)

ret.vind.only_v$em$v$ll.trace
ret.vind.only_v$em$anc_1$ll.trace
ret.vind.only_v$em$v$man.ll.trace
ret.vind.only_v$em$anc_1$man.ll.trace

# > ret.vind.only_v$em$v$ll.trace
# [1] -5335.910 -3897.158 -3840.239 -3833.052 -3831.874 -3831.659 -3831.606 -3831.601
# > ret.vind.only_v$em$anc_1$ll.trace
# [1] -5317.363 -3862.077 -3825.394 -3823.057 -3822.819 -3822.811 -3822.809 -3822.809 ## with reg ll, this gives anc_1 as the best branch
# > ret.vind.only_v$em$v$man.ll.trace
# [1] -4684.888 -3528.673 -3508.780 -3508.135 -3508.113 -3508.112 -3508.112 -3508.112 ## with man.ll, it gives v - which is correct
# > ret.vind.only_v$em$anc_1$man.ll.trace
# [1] -4719.193 -3535.112 -3520.137 -3519.948 -3519.945 -3519.945 -3519.945 -3519.945

ret.spy.archaics$em$v$ll.trace
ret.spy.archaics$em$anc_1$ll.trace
ret.spy.archaics$em$v$man.ll.trace
ret.spy.archaics$em$anc_1$man.ll.trace



dt.sed.analysis.em.ret.vind.only_v <- sed_EM_allbranch(dt.sed.analysis.only_v.lateN[lib == 'Vindija_G1_A9348'],
                                                       sims.dat.only_v, 
                                                       branches = 'anc_1',
                                                       err_rate = 0.001,
                                                       max.iter = 30, ll.converge = 1e-6,
                                                       set.faunal_prop = 'estim',
                                                       set.mh_contam = 'estim',
                                                       p_h_method = 'simple')





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
