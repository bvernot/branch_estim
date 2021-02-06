dump_and_quit <- function() {
    ## Save debugging info to file last.dump.rda
    traceback()
    ## Quit R with error status
    q(status = 1)
}
if (interactive()) {
  options(width=100)
} else {
  options(error = dump_and_quit, width=200)
}


c.args <- commandArgs()
cat('Command:\n\n', c.args, '\n\n', sep = " ")


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
                    help="Number of cores to use. [only currently used with --ll-surface]")

parser$add_argument("-sims", "--sims.dat", required = T,
                    help="")
parser$add_argument("-gts", "--simple-gts", required=T,
                    help="all_simple_gts.tsv.gz")

parser$add_argument("-libs", "--libs", required=F, nargs = '+',
                    help="One or more libraries to group together. The EM still treats them as separate read groups.  To treat these as a single read group, use --merge-libs.")
parser$add_argument("-libs-downsample", "--libs-downsample", required=F, nargs = '+', type='integer',
                    help="Sample N reads from each of the libraries in --libs. Requires --libs. Sampling is performed after reads are filtered e.g. to be polymorphic, ti/tv, etc. If N is larger than the number of reads in the libary, an error is given.")
parser$add_argument("-libs-add-deam", "--libs-add-deam", required=F, nargs = '+', type='double',
                    help="Randomly mark as deaminated x% of non-deaminated reads from each of the libraries in --libs. Requires --libs. Modification is performed after reads are filtered e.g. to be polymorphic, ti/tv, etc.")
parser$add_argument("-merge-libs", "--merge-libs", action="store_true", default=F,
                    help="Merge all requested libraries into a single read group (or entire file, if --libs not given).")

parser$add_argument("-table", "--output-table", required=F,
                    help="Save the results to <output-table> in tsv form.  Also saves an RDS of the em results.")
parser$add_argument("-debug-gts", "--debug-gts-at-time", required=F, type='double', nargs = '+',
                    help="Given a branch and a branch time, report expected and observed p(der) for each genotype category. Requires --branches also.")

parser$add_argument("-ll-surface", "--ll-surface", action="store_true", default=F,
                    help="Get a true maximum likelihood surface for branchtime, over a semi-random grid of times. The surface w/ ll ~ less than ll.thresh away from the max is explored more extensively with each step (nsteps)")
parser$add_argument("-nsteps", "--nsteps", required=F, default=3, type='integer',
                    help="Number of steps to run the ll-surface search [default 3]")
parser$add_argument("-nbins", "--nbins", required=F, default=10, type='integer',
                    help="Number of bins to split each branch into for ll-surface search [default 10]")
parser$add_argument("-tbreaks", "--time-breaks", required=F, default=NULL, type='double',
                    help="Force the first search in ll-surface to do a grid with breaks T over the whole tree DOES NOT CURRENTLY WORK [default not used]")
parser$add_argument("-ll-surface-thresh", "--ll-surface-thresh", required=F, default=5, type='double',
                    help="With --ll-surface, search extensively within ll-thresh of the maximum likelihood [default 5]")


parser$add_argument("-branches", "--branches", required=F, default=c('c', 'v', 'anc_1', 'a', 'anc_2', 'd', 'anc_3'), nargs='+',
                    help="Only run the EM on a single branch. e.g. --branch v.  Required for --debug-gts-at-time.")

## parser$add_argument("-true-branch", "--true-branch", required=F, default='NA',
##                     help="The simulated branch, used in plots and tables.")
## parser$add_argument("-true-branchtime", "--true-branchtime", required=F, default='NA',
##                     help="The simulated branchtime, used in plots and tables.")

parser$add_argument("-set-mh-contam", "--set-mh-contam", required=F, default='estim',
                    help="If provided, constrain MH contamination to this value. [currently only works with a single value, across all RG]")
parser$add_argument("-set-faunal-prop", "--set-faunal-prop", required=F, default='estim',
                    help="If provided, constrain faunal proportion to this value. [currently only works with a single value, across all RG]")

parser$add_argument("-sample-mh", "--sample-mh-from-freqs", default=F, action='store_true',
                    help="Sample the mh 'contamination' from f_mh frequencies rather than a separate column mh")
parser$add_argument("-f-mh", "--f-mh", default='f_mh',
                    help="Use this column ID for modern human allele frquencies.  Default is f_mh")

parser$add_argument("-tag-labels", "--tag-labels", required=F, default='tag', nargs='+',
                    help="One (or more) tag labels for this analysis.")
parser$add_argument("-tags", "--tags", required=F, default='none', nargs='+',
                    help="One (or more) tags for this analysis.")

parser$add_argument("-num-em-iters", "--num-em-iters", type='integer', default=100,
                    help="Maximum number of EM iterations")
parser$add_argument("-ll-converge", "--ll-converge", required=F, default=1e-6,
                    help="Stop the EM search when the LL changes by less than ll-converge [default 1e-6]")


parser$add_argument("-n-qc1", "--n-qc1", type='integer', default=0,
                    help="Artificially add N QC sites that are DERIVED in all hominins. These are used for calculating faunal proportions, and have to be artificially added to simulated data.")
parser$add_argument("-n-qc0", "--n-qc0", type='integer', default=0,
                    help="Artificially add N QC sites that are ANCESTRAL in all hominins. These are not present/useful in real data, so this should mostly be used for debugging.")

parser$add_argument("-downsample", "--downsample", type='integer', default=0,
                    help="Downsample data to N reads. This samples the entire dataset, so e.g. if originally each read group was 25% of the data, these proportions may change.")
parser$add_argument("-block-bootstrap", "--block-bootstrap", type='integer', default=0,
                    help="Split data randomly into N blocks (does not currently use location information), and then resample from these blocks.")
parser$add_argument("-drop-one", "--drop-one-block", type='integer', default=0,
                    help="Instead of resampling, drop block N. Requires --block-bootstrap.")

parser$add_argument("-faunal-der-rate", "--faunal-der-rate", type='double', default=0,
                    help="Assume that any faunal contamination has the derived allele with this probability")
parser$add_argument("-error-rate", "--error-rate", type='double', default=0.001,
                    help="Assume sequencing errors with this probability")


parser$add_argument("-add-contam", "--add-contam", type='double', default=0, nargs='+',
                    help="Artificially add contamination in these proportions")
parser$add_argument("-add-faunal", "--add-faunal", type='double', default=0, nargs='+',
                    help='Artificially add faunal "contamination" in these proportions')
parser$add_argument("-rg-props", "--rg-props", type='double', default=1, nargs='+',
                    help="Randomly split the simulations into read groups with these proportions")

parser$add_argument("-sites", "--sites", required=F, default = 'all', nargs='+', dest = 'site_cat',
                    help='Site categories to use. Can be a list. Options: all, poly_archaic, poly_neand, mh_seg_arc_fixed0, sed_qc_hominin')
parser$add_argument("-method", "--sim-method", required=F, default = 'simple',
                    help='Site categories to use [not currently implemented]')
parser$add_argument("-include-ti", "--include-ti", required=F, action='store_true',
                    help='By default, the script removes transitions. This includes them. Will probably change later')
parser$add_argument("-deam-only", "--deam-only", required=F, action='store_true',
                    help='By default, the script uses all reads. This drops reads with no deamination (as seen in column "deam53"). May change later.')
parser$add_argument("-ag-cols", "--aggregate-gt-columns", required=F, default = c('v_gt', 'c_gt', 'a_gt', 'd_gt'), 
                    nargs='+', dest = 'aggregate_gt_columns',
                    help='Use these columns from the simulations to model p(test=der).  Default is to use v_gt,c_gt,a_gt,d_gt [all archaics]')



parser$add_argument("-script-path", "--script-path", required=F, default = NULL,
                    help='A hack to let R find the path for scripts to source.')
parser$add_argument("-prefix", "--prefix", required=F,
                    help="Prefix for output files [not currently used].")

if (interactive()) {
  # args <- parser$parse_args(strsplit('-gts ~/Downloads/all_simple_gts.tsv.gz --sims ~/Documents/soil_dna_capture/sims.dat.RDS -libs A17273 --prefix what -nc 2', split = ' ')[[1]])
  # args <- parser$parse_args(strsplit('-gts "~/GoogleDrive/branch_point_estimates/all_simple_gts.tsv.gz" --sims "~/GoogleDrive/branch_point_estimates/sims.dat.RDS" -libs A20281 --prefix what -nc 2 -sites all -tag hey', split = ' ')[[1]])
  args <- parser$parse_args(strsplit('-gts ~/Downloads/all_simple_gts.small.tsv.gz --sims ~/Downloads/sims.dat.RDS -libs A20281 --prefix what -nc 2 -sites all -tags hey', split = ' ')[[1]])
  args <- parser$parse_args(strsplit('-gts ~/Downloads/hey.gt.txt.gz --sims ~/Downloads/sims.dat.RDS --prefix what -nc 2 -sites all -tags hey', split = ' ')[[1]])
  args <- parser$parse_args(strsplit('-gts ~/Downloads/test_sims_v_0.7601_ALL.gt.txt.gz --sims ~/Downloads/sims.dat.RDS --prefix what -nc 2 -sites all -tags hey', split = ' ')[[1]])
  args <- parser$parse_args(strsplit('-gts ~/GoogleDrive/branch_point_esimates/data/test_sims_v_0.7601_ALL.gt.txt.gz --sims ~/GoogleDrive/branch_point_esimates/data/sims.dat.RDS --debug-gts-at-time 0.7601 --prefix what -nc 2 -sites all -tags hey --sim-method simple --libs sims009.v.0.7601.1 --rg-props .2 .8 --add-contam 0 .1 --downsample 10000 --n-qc1 1000 --branch v --num-em-iters 5 -table hey.txt', split = ' ')[[1]])
  args <- parser$parse_args(strsplit('-gts ~/GoogleDrive/branch_point_esimates/data/all_simple_gts.mez.tsv.gz --sims ~/GoogleDrive/branch_point_esimates/data/sims_og_newrun.txt --prefix what -nc 2 -sites all -tags hey --sim-method simple --libs Mez1_R5661 -f-mh f_mh.yri --n-qc1 0 --branch v --num-em-iters 10 -table hey.txt', split = ' ')[[1]])
  args <- parser$parse_args(strsplit('-gts ~/GoogleDrive/branch_point_esimates/data/all_simple_gts.mez.tsv.gz --sims ~/GoogleDrive/branch_point_esimates/data/sims_og_newrun.txt --prefix what -nc 2 -sites all -tags hey --sim-method simple --libs Mez1_R5661 -f-mh f_mh.yri --n-qc1 0 --num-em-iters 2 --ll-surface --nsteps 2 -table hey.txt', split = ' ')[[1]])
  args <- parser$parse_args(strsplit('-gts data/all_simple_gts.mez.tsv.gz --sims data/sims_og_newrun.txt --prefix what -nc 2 -sites all -tags hey --sim-method simple --libs Mez1_R5661 -f-mh f_mh.yri --n-qc1 0 --num-em-iters 2 --ll-surface --nsteps 2 -table hey.txt', split = ' ')[[1]])
  #  Mez2_A9180
  # args <- parser$parse_args(strsplit('--splits ~/Documents/index_cross_contam/data/ludovic/splitstats_ludovic_orlando_001.myformat3.txt -nhits 100 --prefix splitstats_ludovic_orlando_001 -nc 2 --sources 150', split = ' ')[[1]])
} else {
  args <- parser$parse_args()
}

print('Arguments:')
print(args)

if (length(args$add_contam) == 1) args$add_contam <- rep(args$add_contam, length(args$rg_props))
if (length(args$add_faunal) == 1) args$add_faunal <- rep(args$add_faunal, length(args$rg_props))

if ( length(args$rg_props) != length(args$add_contam) || length(args$rg_props) != length(args$add_faunal) ) {
  cat('Must provide contamination and faunal proportions for each read group:\n')
  cat('RG:', args$rg_props, '\n')
  cat('contam:', args$add_contam, '\n')
  cat('faunal:', args$add_faunal, '\n')
  q(save='no', status=1)
}

if ( length(args$tags) != length(args$tag_labels) ) {
  cat('Tag labels and tags must be same length:\n')
  cat('tags:', args$tags, '\n')
  cat('tag-labels:', args$tag_labels, '\n')
  q(save='no', status=1)
}


if (is.null(args$script_path)) {
  source(here('R/estim_branchpoints_fns.R'))
  source(here('R/read_data_fns.R'))
  source(here('R/read_sims_fns.R'))
} else {
  source(paste0(args$script_path, '/estim_branchpoints_fns.R'))
  source(paste0(args$script_path, '/read_data_fns.R'))
  source(paste0(args$script_path, '/read_sims_fns.R'))
}

if (args$ncores > 1) {
  cat(sprintf('Registering %d cores\n', args$ncores))
  library(doParallel)
  registerDoParallel(cores=args$ncores)
  getDoParWorkers()
} else {
  cat('Running in single-core mode.\n')
}


# # sims.dat <- readRDS('~/Google Drive/branch_point_esimates/sims.dat.RDS')
# sims.dat <- readRDS(args$sims.dat)
# ## this 0.004 doesn't matter here, because it's modified in a different function!
# sims.dat <- add_linear_p_given_b_t_arcs(sims.dat, fixed_anc_p = 0.004)

sims.dat <- generate_sims_dat(args$sims.dat)

####
## read genotypes

dt.sed.analysis <- read_and_process_genos(args$simple_gts, f_mh.col = args$f_mh, agCols = args$aggregate_gt_columns,
                                          keep.libs = args$libs, keep.libs.downsample = args$libs_downsample, keep.libs.add_deam = args$libs_add_deam,
                                          sample_mh_from_freq = args$sample_mh_from_freqs,
                                          include_ti = args$include_ti,
                                          deam_only = args$deam_only,
                                          merge_libs = args$merge_libs,
                                          site.cats = args$site_cat,
                                          n_qc0 = args$n_qc0,
                                          n_qc1 = args$n_qc1,
                                          downsample = args$downsample, rg_props = args$rg_props,
                                          block_bootstrap = args$block_bootstrap,
                                          drop_one = args$drop_one_block)



####################
## simulate contamination

all_rg = dt.sed.analysis[, unique(rg)]

if (length(args$add_contam) == 1) args$add_contam <- rep(args$add_contam, length(all_rg))
if (length(args$add_faunal) == 1) args$add_faunal <- rep(args$add_faunal, length(all_rg))

for (my_rg.idx in 1:length(all_rg)) {
  my_rg = all_rg[my_rg.idx]
  my_rg.contam = args$add_contam[my_rg.idx]
  my_rg.contam.sites = dt.sed.analysis[, rg == my_rg]
  sum(my_rg.contam.sites)
  sum(my_rg.contam.sites) / length(my_rg.contam.sites)
  my_rg.contam.sites[my_rg.contam.sites] <- sample(c(T,F), sum(my_rg.contam.sites), replace = T,
                                                      prob = c(my_rg.contam,1-my_rg.contam))
  sum(my_rg.contam.sites)
  sum(my_rg.contam.sites) / length(my_rg.contam.sites)
  sum(my_rg.contam.sites) / dt.sed.analysis[rg == my_rg, .N]
  
  dt.sed.analysis[my_rg.contam.sites, sed_gt := mh]
}

setkey(dt.sed.analysis, v_gt, c_gt, a_gt, d_gt)




################
## debug genotypes, if requested

if (!is.null(args$debug_gts_at_time)) {
  
  # args$debug_gts_at_time <- c(.56, .65, .75, .85)
  if (length(args$branches) != 1) {
    cat('debug gts requires a single branch\n', args$branches)
    if(!interactive()) q()
  }
  
  all.gt <- sims.dat$dt.simple_p_given_b_t_arcs[branch == args$branches]
  debug.gts <- foreach(this.t = args$debug_gts_at_time, .combine = rbind) %:%
    foreach (gt.idx = 1:nrow(all.gt), .combine = rbind) %:% 
      foreach (this.rg = all_rg, .combine = rbind) %do% {
        
        this.gt <- all.gt[gt.idx, .(v_gt, c_gt, a_gt, d_gt)]
        this.p.bounds <- sims.dat$dt.simple_p_given_b_t_arcs[this.gt][branch == args$branches]
        
        p.der.exp <- sims.dat$simple_p_given_b_t_arcs(args$branches, this.t, 
                                                      my.gt = this.gt, 
                                                      sims.dat = sims.dat)
        if (dt.sed.analysis[rg == this.rg][this.gt, .N] <= 1) return(data.table())
        cbind(dt.sed.analysis[this.gt][rg == this.rg,
                                       .(p.der.obs = sum(sed_gt)/.N, .N, p.der.exp),
                                       .(v_gt, c_gt, a_gt, d_gt, gt.cat)],
              this.p.bounds[, .(p.low, p.high)],
              this.t, this.rg)
  }
  
  debug.gts[, this.lib := tstrsplit(this.rg,'_',keep=1), this.rg]
  debug.gts[, best.mod := abs(p.der.exp-p.der.obs) == min(abs(p.der.exp-p.der.obs)), .(this.rg, gt.cat)]
  
  ggplot(debug.gts, aes(x=p.der.exp, p.der.obs, size=N)) + geom_point() +
    geom_text_repel(data=debug.gts[abs(p.der.exp-p.der.obs) > .1], 
                    aes(label=paste0(gt.cat,'_',N), color=N), size=3) +
    geom_abline(slope=1) + facet_wrap(this.rg~this.t)
  
  ggplot(debug.gts[N > 2], aes(x=(p.der.exp-p.der.obs) / (p.low-p.high), y=N)) + geom_point() +
    geom_text_repel(data=debug.gts, 
                    aes(label=paste0(gt.cat,'_',N,'_',this.lib), color=N), size=3) +
    geom_abline(slope=1) + 
    scale_y_log10() + facet_wrap(this.rg~this.t)

  ggplot(debug.gts[N > 20], aes(x=(p.der.exp-p.der.obs) / (p.low-p.high), y=N, 
                               color=as.factor(this.t))) +
    geom_point() +
    geom_text_repel(data=debug.gts[N > 20 & this.t == min(args$debug_gts_at_time)],
                    aes(label=paste0(gt.cat,'_',N,'_',this.lib)), size=3, color='black') +
    geom_vline(xintercept = 0, lty=3) +
    scale_y_log10() + facet_wrap(~this.rg)
  
  ggplot(debug.gts, aes(x=(p.der.obs-p.der.exp), y=N, fill=this.lib)) + geom_point() +
    # geom_text_repel(data=debug.gts[N>30 & v_gt == 1], 
    #                 aes(label=paste0(gt.cat,'_',N,'_',this.lib), color=abs(p.high - p.low)), size=3) +
    geom_line(aes(group=gt.cat, color=abs(p.high - p.low))) +
    # geom_abline(slope=1) + 
    scale_y_log10() + facet_wrap(~this.t)
  # + scale_x_log10()

  ggplot(debug.gts, aes(x=(p.der.obs-p.der.exp), y=N, color=this.t)) + 
    geom_point() +
    geom_text_repel(data=debug.gts[N>20 & best.mod == T & gt.cat != '0000'],
                    aes(label=paste0(gt.cat,'_',N), color=this.t), size=3) +
    geom_line(aes(group=gt.cat)) +
    # geom_abline(slope=1) + 
    scale_y_log10() + facet_wrap(~this.rg)
  
  ggplot(debug.gts, aes(x=(p.der.obs-p.der.exp), y=N, color=this.t)) + 
    geom_point() +
    geom_text_repel(data=debug.gts[best.mod == T & lengths(regmatches(gt.cat, gregexpr("0", gt.cat))) == 2],
                    aes(label=paste0(gt.cat,'_',N), color=this.t), size=3) +
    geom_line(aes(group=gt.cat)) +
    # geom_abline(slope=1) + 
    scale_y_log10() + facet_wrap(~this.rg)
  
  # x.em = sed_EM(dt.sed.analysis[rg %like% 'Mez1' & !gt.cat %like% '^222' & !gt.cat %like% '^222'], 
  #               sims.dat, my.branch = args$branches, err_rate = 0.001,
  #               max.iter = 30, ll.converge = 1e-10,
  #               set.faunal_prop = args$set_faunal_prop,
  #               set.mh_contam = args$set_mh_contam,
  #               p_h_method = args$sim_method)
  
  if (!interactive()) q()
}


# sims.dat$dt.sims.p[, gt.cat := paste0(v, c, a, d)]
# sims.dat$dt.sims.p[lengths(regmatches(gt.cat, gregexpr("0", gt.cat))) == 2 & time == 0.6001000, 
#                    gt.cat, 
#                    keyby=.(n_tot=log(n_tot))]
# 
# lengths(regmatches(gt.cat, gregexpr("0", gt.cat))) == 2




if (args$ll_surface) {
  cat('Computing LL surface for branchtime\n')
  if (is.null(args$output_table)) {
    cat('--ll-surface option requires --output-table\n')
    if(!interactive()) q()
  }
  ## cat('args$error_rate', args$error_rate, '\n')
  ## cat('args$faunal_der_rate', args$faunal_der_rate, '\n')

  #print(dt.sed.analysis[, .N, rg])
  
  dt.x.new = grid_t_em_theta(dt.sed.analysis,
                             sims.dat,
                             my.branches = args$branches,
                             err_rate = args$error_rate,
                             faunal_der_rate = args$faunal_der_rate,
                             bins.t = args$nbins,
                             t.breaks= args$time_breaks,
                             max.iter = args$num_em_iters, ll.converge = args$ll_converge, 
                             nsteps = args$nsteps, ll.thresh = args$ll_surface_thresh)
  for (i in 1:length(args$tags)) {
    dt.x.new[, (args$tag_labels[i]) := args$tags[i]]
  }
  fwrite(dt.x.new, args$output_table, sep='\t')
  if(!interactive()) q()
}





################
## run the EM

## ISSUE - not sure that ll.converge works?
if (length(args$branches) == 1) {
  cat('Running EM on single branch\n')
  x.em = sed_EM(dt.sed.analysis,
                sims.dat,
                my.branch = args$branches,
                err_rate = args$error_rate,
                faunal_der_rate = args$faunal_der_rate,
                max.iter = args$num_em_iters,
                ll.converge = args$ll_converge,
                set.faunal_prop = args$set_faunal_prop,
                set.mh_contam = args$set_mh_contam,
                p_h_method = args$sim_method)

  if (!is.null(args$output_table)) {

      dt.x.new = ll_ret_to_dt(x.em, args, dt.sed.analysis)
      
      for (i in 1:length(args$tags)) {
          dt.x.new[, (args$tag_labels[i]) := args$tags[i]]
      }
      
      fwrite(dt.x.new, args$output_table, sep='\t')
      x.em$args <- args
      saveRDS(x.em, paste0(args$output_table,'.RDS'))
  }
} else {
  cat('Running EM on multiple branches, with simple grid after\n')
  dt.sed.analysis.ret <- run_simple_analysis(dt.sed.analysis, sims.dat, branches = args$branches,
                                             err_rate = args$error_rate,
                                             faunal_der_rate = args$faunal_der_rate,
                                             blocks = 10, nbootstraps = 10, max.iter = 40)
  dt.sed.analysis.ret$tag <- args$tag
  saveRDS(dt.sed.analysis.ret, paste0(args$output_table,'.full.RDS'))
}


