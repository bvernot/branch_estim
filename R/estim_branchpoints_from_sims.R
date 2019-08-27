options(width=200)
# install.packages('R.utils')
# install.packages('argparse')
library(argparse)
library(here)

# options(error=recover)
# options(error=NULL)


# create parser object
parser <- ArgumentParser()

parser$add_argument("-v", "--verbose", action="store_true", default=F,
                    help="Print extra output")
parser$add_argument("-nc", "--ncores", type='integer', default=1,
                    help="Number of cores to use. [not currently used?]")
parser$add_argument("-sims", "--sims.dat", required = T,
                    help="")
parser$add_argument("-gts", "--simple-gts", required=T,
                    help="all_simple_gts.tsv.gz")
parser$add_argument("-libs", "--libs", required=F,
                    help="One or more libraries to group together. The EM still treats them as separate read groups.  To treat these as a single read group, use --merge-libs.")
parser$add_argument("-merge-libs", "--merge-libs", action="store_true", default=F,
                    help="Merge all requested libraries into a single read group (or entire file, if --libs not given).")
parser$add_argument("-tag", "--tag", required=F, default='none',
                    help="One (or more, in the future) tags for this analysis.")
parser$add_argument("-niter", "--num-iters", type='integer', default=100,
                    help="Number of EM iterations")
parser$add_argument("-n-qc1", "--n-qc1", type='integer', default=1000,
                    help="Artificially add N QC sites that are DERIVED in all hominins. These are used for calculating faunal proportions, and have to be artificially added to simulated data.")
parser$add_argument("-n-qc0", "--n-qc0", type='integer', default=1000,
                    help="Artificially add N QC sites that are ANCESTRAL in all hominins. These are not present/useful in real data, so this should mostly be used for debugging.")
parser$add_argument("-add-contam", "--add-contam", type='double', default=0, nargs='+',
                    help="Artificially add contamination in these proportions")
parser$add_argument("-add-faunal", "--add-faunal", type='double', default=0, nargs='+',
                    help='Artificially add faunal "contamination" in these proportions')
parser$add_argument("-rg-props", "--rg-props", type='double', default=1, nargs='+',
                    help="Randomly split the simulations into read groups with these proportions")
parser$add_argument("-sites", "--site-cat", required=F, default = 'all',
                    help='Site categories to use [not currently implemented]')
parser$add_argument("-method", "--sim-method", required=F, default = 'simple',
                    help='Site categories to use [not currently implemented]')
parser$add_argument("-prefix", "--prefix", required=T,
                    help="Prefix for output files.")

if (interactive()) {
  # args <- parser$parse_args(strsplit('-gts ~/Downloads/all_simple_gts.tsv.gz --sims ~/Documents/soil_dna_capture/sims.dat.RDS -libs A17273 --prefix what -nc 2', split = ' ')[[1]])
  # args <- parser$parse_args(strsplit('-gts "~/GoogleDrive/branch_point_estimates/all_simple_gts.tsv.gz" --sims "~/GoogleDrive/branch_point_estimates/sims.dat.RDS" -libs A20281 --prefix what -nc 2 -sites all -tag hey', split = ' ')[[1]])
  args <- parser$parse_args(strsplit('-gts ~/Downloads/all_simple_gts.small.tsv.gz --sims ~/Downloads/sims.dat.RDS -libs A20281 --prefix what -nc 2 -sites all -tag hey', split = ' ')[[1]])
  args <- parser$parse_args(strsplit('-gts ~/Downloads/hey.gt.txt.gz --sims ~/Downloads/sims.dat.RDS --prefix what -nc 2 -sites all -tag hey', split = ' ')[[1]])
  args <- parser$parse_args(strsplit('-gts ~/Downloads/test_sims_v_0.7601_ALL.gt.txt.gz --sims ~/Downloads/sims.dat.RDS --prefix what -nc 2 -sites all -tag hey', split = ' ')[[1]])
  args <- parser$parse_args(strsplit('-gts ~/GoogleDrive/branch_point_esimates/data/test_sims_v_0.7601_ALL.gt.txt.gz --sims ~/GoogleDrive/branch_point_esimates/data/sims.dat.RDS --prefix what -nc 2 -sites all -tag hey', split = ' ')[[1]])
  # args <- parser$parse_args(strsplit('--splits ~/Documents/index_cross_contam/data/ludovic/splitstats_ludovic_orlando_001.myformat3.txt -nhits 100 --prefix splitstats_ludovic_orlando_001 -nc 2 --sources 150', split = ' ')[[1]])
} else {
  args <- parser$parse_args()
}


if (length(args$add_contam) == 1) args$add_contam <- rep(args$add_contam, length(args$rg_props))
if (length(args$add_faunal) == 1) args$add_faunal <- rep(args$add_faunal, length(args$rg_props))

if ( length(args$rg_props) != length(args$add_contam) || length(args$rg_props) != length(args$add_faunal) ) {
  cat('Must provide contamination and faunal proportions for each read group:\n')
  cat('RG:', args$rg_props, '\n')
  cat('contam:', args$add_contam, '\n')
  cat('faunal:', args$add_faunal, '\n')
  exit()
}


source(here('R/estim_branchpoints_fns.R'))

registerDoParallel(cores=args$ncores)
getDoParWorkers()


# sims.dat <- readRDS('~/Google Drive/branch_point_esimates/sims.dat.RDS')
sims.dat <- readRDS(args$sims.dat)
## this 0.004 doesn't matter here, because it's modified in a different function!
sims.dat <- add_linear_p_given_b_t_arcs(sims.dat, fixed_anc_p = 0.004)




# dt.sed.og <- fread('~/Google Drive/branch_point_esimates/all_simple_gts.tsv.gz')
# dt.sed.og <- fread('~/Downloads/all_simple_gts.tsv.gz')
dt.sed.og <- fread(args$simple_gts)

####
## confirm that all of the correct columns are present

req_columns <- c('sed_gt', 'v_gt', 'c_gt', 'a_gt', 'd_gt', 'f_mh', 'lib')
if (sum(!req_columns %in% colnames(dt.sed.og)) > 0) {
  cat('Not all required columns are present:\n')
  cat('Required: ', req_columns, '\n')
  cat('Missing: ', req_columns[!req_columns %in% colnames(dt.sed.og)], '\n')
}

if (!is.null(args$libs)) {
  dt.sed <- dt.sed.og[lib %in% args$libs]
} else {
  ## not really necessary, takes up a lot more space...
  dt.sed <- data.table(dt.sed.og)
}

if (args$merge_libs) {
  dt.sed[, lib := 'merged_libs']
}

# ## set up simulated data, because I didn't previously fill this in
# dt.sed[, lib := 'sim009']
# setnames(dt.sed, c('sed', 'mh_f'), c('sed_gt', 'f_mh'))

dt.sed.poly.full <- dt.sed[!(v_gt == c_gt & v_gt == a_gt & v_gt == d_gt)]
# dt.sed.poly.full[, deam53 := rep(c(T,T,F,F,F), length.out = .N)]
# dt.sed.poly.full[, f_mh := f_mh / 99]
# dt.sed.poly.full[, pos := NULL]

dt.sed.mh.full <- dt.sed[v_gt == c_gt & v_gt == a_gt & v_gt == d_gt & v_gt == 0 & f_mh > 0]
# dt.sed.mh.full[, deam53 := rep(c(T,T,F,F,F), length.out = .N)]
# dt.sed.mh.full[, f_mh := f_mh / 99]
# dt.sed.mh.full[, pos := NULL]

dt.sed.qc.full <- foreach(my.lib = dt.sed.poly.full[, unique(lib)], .combine = rbind) %do% {
  ## just duplicate the first row of dt.sed.poly.full the correct number of times
  dt.sed.qc.full <- dt.sed.poly.full[lib == my.lib][rep(1, args$n_qc0 + args$n_qc1)]
  
  qc_fill_freq_or_hap <- c(rep(0,args$n_qc0), rep(1,args$n_qc1))
  qc_fill_gt <- c(rep(0,args$n_qc0), rep(2,args$n_qc1))
  
  dt.sed.qc.full[, f_mh := qc_fill_freq_or_hap]
  dt.sed.qc.full[, mh := qc_fill_freq_or_hap]
  dt.sed.qc.full[, v_gt := qc_fill_gt]
  dt.sed.qc.full[, c_gt := qc_fill_gt]
  dt.sed.qc.full[, a_gt := qc_fill_gt]
  dt.sed.qc.full[, d_gt := qc_fill_gt]
  dt.sed.qc.full[, sed_gt := qc_fill_freq_or_hap]
  dt.sed.qc.full
}

## not really necessary, takes up a lot more space...
dt.sed.poly <- data.table(dt.sed.poly.full)
dt.sed.mh <- data.table(dt.sed.mh.full)
dt.sed.qc <- data.table(dt.sed.qc.full)
rbind(dt.sed.mh, dt.sed.qc, dt.sed.poly)


####################
## simulate contamination

## four read groups
dt.sed.analysis.mh_all.contam1.x4 = rbind(dt.sed.qc,dt.sed.poly,dt.sed.mh)
# dt.sed.analysis = rbind(dt.sed.qc,dt.sed.poly[!(v_gt == c_gt & v_gt == a_gt)])
dt.sed.analysis.mh_all.contam1.x4[, rg := paste0(lib, '_rg_', 1:4)]
all_rg = dt.sed.analysis.mh_all.contam1.x4[, unique(rg)]
all_rg_contam = c(0,0.025,.1,.15)
for (my_rg.idx in 1:length(all_rg)) {
  my_rg = all_rg[my_rg.idx]
  my_rg.contam = all_rg_contam[my_rg.idx]
  my_rg.contam.sites = dt.sed.analysis.mh_all.contam1.x4[, rg == my_rg]
  sum(my_rg.contam.sites)
  sum(my_rg.contam.sites) / length(my_rg.contam.sites)
  my_rg.contam.sites[my_rg.contam.sites] <- sample(c(T,F), sum(my_rg.contam.sites), replace = T,
                                                      prob = c(my_rg.contam,1-my_rg.contam))
  sum(my_rg.contam.sites)
  sum(my_rg.contam.sites) / length(my_rg.contam.sites)
  sum(my_rg.contam.sites) / dt.sed.analysis.mh_all.contam1.x4[rg == my_rg, .N]
  
  dt.sed.analysis.mh_all.contam1.x4[my_rg.contam.sites, sed_gt := mh]
}
dt.sed.analysis.mh_all.contam1.x4.simple.fixed_anc_p_0.004 <- eval_sed_t_and_mh(dt.sed.analysis.mh_all.contam1.x4, 'simple', max.iter = 40)
plot_eval_sed_t_and_mh(dt.sed.analysis.mh_all.contam1.x4.simple.fixed_anc_p_0.004, true_branchtime = 0.7601, true_mh_contam = all_rg_contam)
ggsave('~/Google Drive/branch_point_esimates/sims009__with_mh_contam_x4.fixed_anc_p_0.004.mh_all.pdf', width=7, height=5)
# ggsave('~/Google Drive/branch_point_esimates/sims009__with_mh_contam_x4.fixed_anc_p_0.004.mh_all.pdf', width=7, height=5)

## this doesn't go negative until very small numbers -16493.28 vs -16493.28, diff: -5.070935e-06
## but with RG it does go negative quite a bit earlier.
## check to see if the optim of the grid matches optim of the EM?  Hard to eval that for the multiple RG thing,
## since we can't do the grid for many RG.
## double check that there isn't some issue with rg?  ugh how to find this?  maybe gamma updates weirdly?
dt.sed.analysis.mh_all.contam1.x4.one_rg <- data.table(dt.sed.analysis.mh_all.contam1.x4)
dt.sed.analysis.mh_all.contam1.x4.one_rg[, rg := 'hey']
dt.sed.analysis.mh_all.contam1.x4.one_rg.simple.fixed_anc_p_0.004 <- eval_sed_t_and_mh(dt.sed.analysis.mh_all.contam1.x4.one_rg, 'simple', max.iter = 40)
plot_eval_sed_t_and_mh(dt.sed.analysis.mh_all.contam1.x4.one_rg.simple.fixed_anc_p_0.004, true_branchtime = 0.7601, true_mh_contam = mean(all_rg_contam))
#
## ok, this is really weird, this goes negative much sooner, and with much larger diffs. -33177.13 -33227.22 -50.08743 
dt.sed.analysis.mh_all.contam1.x4.one_rg2 <- data.table(dt.sed.analysis.mh_all.contam1.x4)
dt.sed.analysis.mh_all.contam1.x4.one_rg2[, rg := 'hey']
dt.sed.analysis.mh_all.contam1.x4.one_rg2[sample(.N, 10), rg := 'what']
dt.sed.analysis.mh_all.contam1.x4.one_rg2.simple.fixed_anc_p_0.004 <- eval_sed_t_and_mh(dt.sed.analysis.mh_all.contam1.x4.one_rg2, 'simple', max.iter = 40)
plot_eval_sed_t_and_mh(dt.sed.analysis.mh_all.contam1.x4.one_rg2.simple.fixed_anc_p_0.004, true_branchtime = 0.7601, true_mh_contam = mean(all_rg_contam))
#
dt.sed.analysis.mh_all.contam1.x4.one_rg2.simple.fixed_anc_p_0.004_I_fixed_the_rg_bug_in_calc_man_lik <- eval_sed_t_and_mh(dt.sed.analysis.mh_all.contam1.x4.one_rg2, 'simple', max.iter = 40)
plot_eval_sed_t_and_mh(dt.sed.analysis.mh_all.contam1.x4.one_rg2.simple.fixed_anc_p_0.004_I_fixed_the_rg_bug_in_calc_man_lik, true_branchtime = 0.7601, true_mh_contam = mean(all_rg_contam))

dt.sed.analysis.mh_all.contam1.x4.one_rg2.simple.fixed_anc_p_0.004_faunal_also <- eval_sed_t_and_mh(dt.sed.analysis.mh_all.contam1.x4.one_rg2, 'simple', max.iter = 40, set.faunal_prop = 'estim')
plot_eval_sed_t_and_mh(dt.sed.analysis.mh_all.contam1.x4.one_rg2.simple.fixed_anc_p_0.004_faunal_also, true_branchtime = 0.7601, true_mh_contam = mean(all_rg_contam))



sims.dat <- add_linear_p_given_b_t_arcs(sims.dat, fixed_anc_p = 0.03)
dt.sed.analysis.mh_all.contam1.x4.simple.fixed_anc_p_0.03 <- eval_sed_t_and_mh(dt.sed.analysis.mh_all.contam1.x4, 'simple', max.iter = 40)
plot_eval_sed_t_and_mh(dt.sed.analysis.mh_all.contam1.x4.simple.fixed_anc_p_0.03, true_branchtime = 0.7601, true_mh_contam = all_rg_contam)
ggsave('~/Google Drive/branch_point_esimates/sims009__with_mh_contam_x4.fixed_anc_p_0.03.mh_all.pdf', width=7, height=5)

sims.dat <- add_linear_p_given_b_t_arcs(sims.dat, fixed_anc_p = 0.004) 


save.image('~/Google Drive/branch_point_esimates/sims009__evals.Rdata') 

setkey(dt.sed, v_gt,c_gt,a_gt,d_gt)



x.em.all_params.n100 = sed_EM(dt.sed.analysis.mh_n5k.contam1, sims.dat, my.branch = 'v', err_rate = 0.001, 
                              max.iter = 40, ll.converge = 1e-10,
                              set.faunal_prop = 0, p_h_method = 'simple')




#################
### look at four different sets of sims


## four read groups
dt.sed.analysis.mh_all.contam1.x4.4_libs = rbind(dt.sed.qc.full,dt.sed.poly.full,dt.sed.mh.full)
# dt.sed.analysis = rbind(dt.sed.qc,dt.sed.poly[!(v_gt == c_gt & v_gt == a_gt)])
dt.sed.analysis.mh_all.contam1.x4.4_libs[, rg := paste0(lib, '_rg_', 1:4)]
all_rg = dt.sed.analysis.mh_all.contam1.x4.4_libs[, unique(rg)]
all_rg_contam = rep(c(0,0.025,.1,.15), length.out = length(all_rg))
dt.sed.analysis.mh_all.contam1.x4.4_libs[, is_contam := F]
for (my_rg.idx in 1:length(all_rg)) {
  my_rg = all_rg[my_rg.idx]
  cat(my_rg)
  my_rg.contam = all_rg_contam[my_rg.idx]
  my_rg.contam.sites = dt.sed.analysis.mh_all.contam1.x4.4_libs[, rg == my_rg]
  sum(my_rg.contam.sites)
  sum(my_rg.contam.sites) / length(my_rg.contam.sites)
  my_rg.contam.sites[my_rg.contam.sites] <- sample(c(T,F), sum(my_rg.contam.sites), replace = T,
                                                   prob = c(my_rg.contam,1-my_rg.contam))
  sum(my_rg.contam.sites)
  sum(my_rg.contam.sites) / length(my_rg.contam.sites)
  sum(my_rg.contam.sites) / dt.sed.analysis.mh_all.contam1.x4.4_libs[rg == my_rg, .N]
  
  dt.sed.analysis.mh_all.contam1.x4.4_libs[my_rg.contam.sites, sed_gt := mh]
  dt.sed.analysis.mh_all.contam1.x4.4_libs[my_rg.contam.sites, is_contam := T]
}
dt.sed.analysis.mh_all.contam1.x4.4_libs.results = list()
for (my.lib in dt.sed.analysis.mh_all.contam1.x4.4_libs[, unique(lib)]) {
  dt.sed.analysis.mh_all.contam1.x4.4_libs.results[[my.lib]] <- 
    eval_sed_t_and_mh(dt.sed.analysis.mh_all.contam1.x4.4_libs[lib == my.lib], 'simple', max.iter = 40)
  plot_eval_sed_t_and_mh(dt.sed.analysis.mh_all.contam1.x4.4_libs.results[[my.lib]], 
                         true_branchtime = 0.7601, true_mh_contam = all_rg_contam)
}

dt.sed.analysis.mh_all.contam1.x4.4_libs.one_rg.3 <- dt.sed.analysis.mh_all.contam1.x4.4_libs[lib == 'sims009.v.0.7601.3']
dt.sed.analysis.mh_all.contam1.x4.4_libs.one_rg.3[, rg := 'hey']
dt.sed.analysis.mh_all.contam1.x4.4_libs.results[['sims009.v.0.7601.3_one_rg']] <- 
  eval_sed_t_and_mh(dt.sed.analysis.mh_all.contam1.x4.4_libs.one_rg.3, 'simple', max.iter = 40)
plot_eval_sed_t_and_mh(dt.sed.analysis.mh_all.contam1.x4.4_libs.results[['sims009.v.0.7601.3_one_rg']], 
                       true_branchtime = 0.7601, true_mh_contam = mean(all_rg_contam))



for (my.lib in dt.sed.analysis.mh_all.contam1.x4.4_libs[, unique(lib)]) {
  my.lib.tag <- paste0(my.lib, '_full')
  dt.sed.analysis.mh_all.contam1.x4.4_libs.results[[my.lib.tag]] <- 
    eval_sed_t_and_mh(dt.sed.analysis.mh_all.contam1.x4.4_libs[lib == my.lib], 'full', max.iter = 40)
  plot_eval_sed_t_and_mh(dt.sed.analysis.mh_all.contam1.x4.4_libs.results[[my.lib.tag]],
                         true_branchtime = 0.7601, true_mh_contam = all_rg_contam)
}


pdf('~/Google Drive/branch_point_esimates/sims009__with_mh_contam_x4.fixed_anc_p_0.004.mh_all.4_libs.pdf', width=7, height=5)
for (my.lib in sort(names(dt.sed.analysis.mh_all.contam1.x4.4_libs.results))) {
  plot_eval_sed_t_and_mh(dt.sed.analysis.mh_all.contam1.x4.4_libs.results[[my.lib]], 
                         true_branchtime = 0.7601, true_mh_contam = all_rg_contam,
                         plot.title = my.lib)
}
dev.off()

'sims009.v.0.7601.3'
  


###########
## there is a major problem if the genotype data.table has a branch column!
## I should subset it to just the necessary columns first, in the EM.






#######################
###
"
Trying to debug this issue where the likelihood changes in a negative way sometimes.
I can get it to happen almost all the time with the following example, with just 6 sites.
But it only seems to happen when I *only* optimize mh_contam or faunal_prop separately,
which uses 'optimize'.  When I optimize both concurrently, using constraint optimization,
the likelihood doesn't go negative.
"

nsites.from.cats = 10
dt.sed.analysis.mh_n1k.contam2 = rbind(dt.sed.qc[sample(.N,nsites.from.cats)],
                                       dt.sed.poly[sample(.N,nsites.from.cats)],
                                       dt.sed.mh[sample(.N,nsites.from.cats)]
                                       )
# dt.sed.analysis = rbind(dt.sed.qc,dt.sed.poly[!(v_gt == c_gt & v_gt == a_gt)])
dt.sed.analysis.mh_n1k.contam2[, rg := paste0(lib, '_deam_', deam53)]
dt.sed.analysis.mh_n1k.contam2[, rg := 'hey']
# dt.sed.analysis.mh_n1k.contam2[deam53 == F & sample(c(T,F), .N, replace = T, prob = c(.1,.9)),
#                                .N, .(mh, sed_gt)]
# dt.sed.analysis.mh_n1k.contam2[deam53 == F & sample(c(T,F), .N, replace = T, prob = c(.1,.9)),
#                                sed_gt := mh]

x.em.all_params.n100 = sed_EM(dt.sed.analysis.mh_n1k.contam2, sims.dat, my.branch = 'v', err_rate = 0.001,
                              max.iter = 40, ll.converge = 1e-10,
                              set.faunal_prop = 0, p_h_method = 'simple', fail_on_neg_change = T)
x.em.all_params.n100 = sed_EM(dt.sed.analysis.mh_n1k.contam2, sims.dat, my.branch = 'v', err_rate = 0.001,
                              max.iter = 40, ll.converge = 1e-10,
                              set.faunal_prop = 0, p_h_method = 'simple', fail_on_neg_change_q = T)
x.em.all_params.n100 = sed_EM(dt.sed.analysis.mh_n1k.contam2, sims.dat, my.branch = 'v', err_rate = 0.001,
                              max.iter = 40, ll.converge = 1e-10,
                              set.faunal_prop = 0.5053096, p_h_method = 'simple', fail_on_neg_change = T)
# x.em.all_params.n100 = sed_EM(dt.sed.analysis.mh_n1k.contam2, sims.dat, my.branch = 'v', err_rate = 0.001,
#                               max.iter = 40, ll.converge = 1e-10,
#                               set.mh_contam = 0, p_h_method = 'simple', fail_on_neg_change = T)
x.em.all_params.n100 = sed_EM(dt.sed.analysis.mh_n1k.contam2, sims.dat, my.branch = 'v', err_rate = 0.001,
                              max.iter = 40, ll.converge = 1e-10, p_h_method = 'simple')

$max.ll
[1] -0.04951786

$man.max.ll
[1] -0.00180017

## simple case that gives an error:
## one or two sites, f_mh and archaics and sed all zeros, and p(der) == branch time. goes negative after 2 iterations
## both manual and iter/qval lik go negative
dt.sed.analysis.mh_n1k.contam2 = data.table(dt.sed.qc[sample(.N,nsites.from.cats)], rg = 'hey')
## starts with gamma_h_0 = .9
## and mh contam starts at .1, faunal is fixed at 0
x.em.all_params.n100 = sed_EM(dt.sed.analysis.mh_n1k.contam2, sims.dat, my.branch = 'v', err_rate = 0.001,
                              max.iter = 40, ll.converge = 1e-10,
                              set.faunal_prop = 0, p_h_method = 'bt', fail_on_neg_change = T)
x.em.all_params.n100.both = sed_EM(dt.sed.analysis.mh_n1k.contam2, sims.dat, my.branch = 'v', err_rate = 0.001,
                              max.iter = 40, ll.converge = 1e-10,
                              p_h_method = 'bt', fail_on_neg_change = T)
x.em.all_params.n100.fix_fp_.5 = sed_EM(dt.sed.analysis.mh_n1k.contam2, sims.dat, my.branch = 'v', err_rate = 0.001,
                              max.iter = 40, ll.converge = 1e-10, set.mh_contam = 'estim2',
                              set.faunal_prop = 0.5053096, p_h_method = 'bt',
                              fail_on_neg_change = T)
## do a grid search here, try to find the optimal value?

###########
## THIS SHOULD NOT HAPPEN!
|rg  | mh_contam| faunal_prop|
  |:---|---------:|-----------:|
  |hey | 0.9002624|   0.5053096|
  



ITER_q  : 2 -0.01531858738 -0.04951654664 0.03419795926 v 0.889291
ITER_q  : 3 -0.0495172 -0.01531859 -0.03419861 v 0.8892908 
ITER_q  : 4 -0.01531762 -0.0495172 0.03419959 v 0.8892908 
ITER_q  : 5 -0.04951786 -0.01531762 -0.03420024 v 0.8892908 



|rg  | mh_contam| faunal_prop|
  |:---|---------:|-----------:|
  |hey | 0.4946904|   0.5053096|
  
  |rg  | mh_contam| faunal_prop|
  |:---|---------:|-----------:|
  |hey | 0.4946904|   0.5053096|
  



constrOptim(.1, function(params) {
  # cat(params, '\n')
  - params[1] + .2},
  grad = NULL,
  ui=rbind(-1,  # the -x-y > -.8
           1 ),  # the y > 0
  ci=c(-.8,0))



date()

# x.em.all_params.n100 = sed_EM(dt.sed.analysis.mh_n1k, sims.dat, my.branch = 'v', err_rate = 0.001, 
#                               max.iter = 10, ll.converge = 1e-10,
#                               set.faunal_prop = 0, p_h_method = 'simple')




## full has the problem of having fit a model where p_der given fixed archaics is wildly inaccurate
eval.ret.analysis.full <- eval_sed_t_and_mh(dt.sed.analysis, 'full')
plot_eval_sed_t_and_mh(eval.ret.analysis.full)

date()

# eval.ret.analysis.simple.i500 <- eval_sed_t_and_mh(dt.sed.analysis, 'simple', max.iter = 500)
# plot_eval_sed_t_and_mh(eval.ret.analysis.simple.i500)
# 
# eval.ret.analysis.full.i500 <- eval_sed_t_and_mh(dt.sed.analysis, 'full', max.iter = 500)
# plot_eval_sed_t_and_mh(eval.ret.analysis.full.i500)
# 
# date()



#################
## test method by "simulating" data, simply from the frequencies
if (F) {
  
  
  simulate_data_sed('bt', .75, sims.dat)
  simulate_data_sed('bt2', .75, sims.dat)
  simulate_data_sed('bt3', .75, sims.dat)
  simulate_data_sed('simple', .75, sims.dat)[v_gt == 0 & c_gt == 0 & a_gt == 0]
  simulate_data_sed('full', .75, sims.dat)[v_gt == 0 & c_gt == 0 & a_gt == 0]
  

  ## very very simple fake data, using p_der=bt
  dt.sed.poly.tst <- simulate_data_sed('bt', .75, sims.dat)
  eval.ret.bt <- eval_sed_t_and_mh(dt.sed.poly.tst, 'bt')
  plot_eval_sed_t_and_mh(eval.ret.bt)
  
  ## slightly more complicated fake data, using bt2, p_der=bt/(gt_c+1)
  dt.sed.poly.tst.bt2 <- simulate_data_sed('bt2', .75, sims.dat)
  eval.ret.bt2 <- eval_sed_t_and_mh(dt.sed.poly.tst.bt2, p_h_method = 'bt2')
  plot_eval_sed_t_and_mh(eval.ret.bt2)

  ## same scenario, but randomly shuffle the sediment data
  dt.sed.poly.tst.bt2.shuf <- data.table(dt.sed.poly.tst.bt2)
  dt.sed.poly.tst.bt2.shuf[, sed_gt := sample(sed_gt)]
  eval.ret.bt2.shuf <- eval_sed_t_and_mh(dt.sed.poly.tst.bt2.shuf, p_h_method = 'bt2')
  plot_eval_sed_t_and_mh(eval.ret.bt2.shuf)

  ## slightly more complicated fake data, using bt3, p_der=bt/(gt_c+1) and p_der=0 if altai=0
  dt.sed.poly.tst.bt3 <- simulate_data_sed('bt3', .75, sims.dat)
  eval.ret.bt3 <- eval_sed_t_and_mh(dt.sed.poly.tst.bt3, p_h_method = 'bt3')
  plot_eval_sed_t_and_mh(eval.ret.bt3)
  
  ## use the real simulated model, but have p_der change linearly along a branch
  ## (rather than with a fitted spline)
  dt.sed.poly.tst.simple <- simulate_data_sed('simple', .75, sims.dat)
  eval.ret.simple <- eval_sed_t_and_mh(dt.sed.poly.tst.simple, p_h_method = 'simple')
  plot_eval_sed_t_and_mh(eval.ret.simple)
  
  ## use the real simulated model, but have p_der modeled with a fitted spline
  dt.sed.poly.tst.full <- simulate_data_sed('full', .75, sims.dat)
  eval.ret.full <- eval_sed_t_and_mh(dt.sed.poly.tst.full, p_h_method = 'full', max.iter = 20)
  plot_eval_sed_t_and_mh(eval.ret.full)
  ##
  
}
  

