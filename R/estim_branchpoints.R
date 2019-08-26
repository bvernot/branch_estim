options(width=200)
# install.packages('R.utils')
# install.packages('argparse')
library(argparse)

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
                    help="One (or more, in the future) librarie to group together.")
parser$add_argument("-tag", "--tag", required=F, default='none',
                    help="One (or more, in the future) tags for this analysis.")
parser$add_argument("-niter", "--num-iters", type='integer', default=100,
                    help="Number of EM iterations")
parser$add_argument("-sites", "--site-cat", required=F, default = 'all',
                    help='Site categories to use [not currently implemented]')
parser$add_argument("-prefix", "--prefix", required=T,
                    help="Prefix for output files.")

if (interactive()) {
  # args <- parser$parse_args(strsplit('-gts ~/Downloads/all_simple_gts.tsv.gz --sims ~/Documents/soil_dna_capture/sims.dat.RDS -libs A17273 --prefix what -nc 2', split = ' ')[[1]])
  # args <- parser$parse_args(strsplit('-gts "~/Google Drive/branch_point_estimates/all_simple_gts.tsv.gz" --sims "~/Google Drive/branch_point_estimates/sims.dat.RDS" -libs A20281 --prefix what -nc 2 -sites all -tag hey', split = ' ')[[1]])
  args <- parser$parse_args(strsplit('-gts ~/Downloads/all_simple_gts.small.tsv.gz --sims ~/Downloads/sims.dat.RDS -libs A20281 --prefix what -nc 2 -sites all -tag hey', split = ' ')[[1]])
  args <- parser$parse_args(strsplit('-gts ~/Downloads/hey.gt.txt.gz --sims ~/Downloads/sims.dat.RDS --prefix what -nc 2 -sites all -tag hey', split = ' ')[[1]])
  # args <- parser$parse_args(strsplit('--splits ~/Documents/index_cross_contam/data/ludovic/splitstats_ludovic_orlando_001.myformat3.txt -nhits 100 --prefix splitstats_ludovic_orlando_001 -nc 2 --sources 150', split = ' ')[[1]])
} else {
  args <- parser$parse_args()
}

# library(ggplot2)
# library(data.table)
# library(dplyr)
# install.packages('cobs')
# library(cobs)
# library(foreach)
# install.packages('doParallel')
# library(doParallel)
# 
# setDTthreads(1)
# 
source('~/Google Drive/branch_point_esimates/estim_branchpoints_fns.R')

registerDoParallel(cores=args$ncores)
getDoParWorkers()


# sims.dat <- readRDS('~/Google Drive/branch_point_esimates/sims.dat.RDS')
sims.dat <- readRDS(args$sims.dat)

if (F) {
  x = sims.dat$branch.bounds[branch == 'v', seq(t.low, t.high, .01)]
  plot(x, sapply(x, function(x) sims.dat$p_gt_given_b_t_arcs('v', x, list(0,0,0,0), sims.dat)), ylim=c(0,.5), pch='0', col='red')
  points(x, sapply(x, function(x) sims.dat$p_gt_given_b_t_arcs('v', x, list(1,0,0,0), sims.dat)), pch='1')
  points(x, sapply(x, function(x) sims.dat$p_gt_given_b_t_arcs('v', x, list(2,0,0,0), sims.dat)), pch='2')
  plot(x, sapply(x, function(x) sims.dat$p_gt_given_b_t_arcs('v', x, list(0,0,0,1), sims.dat)), ylim=c(0,.5))
  points(x, sapply(x, function(x) sims.dat$p_gt_given_b_t_arcs('v', x, list(0,0,0,2), sims.dat)))
  points(x, sapply(x, function(x) sims.dat$p_gt_given_b_t_arcs('v', x, list(0,0,0,0), sims.dat)))
}

sims.dat <- add_linear_p_given_b_t_arcs(sims.dat)


# #######
# ## make a new function that gives you a 'simple' p_der
# sims.dat$dt.simple_p_given_b_t_arcs <- data.table(expand.grid(v_gt=0:2,c_gt=0:2,a_gt=0:2,d_gt=0:2,branch=sims.dat$branches))
# sims.dat$dt.simple_p_given_b_t_arcs <- merge(sims.dat$dt.simple_p_given_b_t_arcs, sims.dat$branch.bounds, by='branch')
# sims.dat$dt.simple_p_given_b_t_arcs[, p.low := sims.dat$p_gt_given_b_t_arcs(branch, t.low, list(v_gt,c_gt,a_gt,d_gt), sims.dat),
#                                     keyby=.(v_gt,c_gt,a_gt,d_gt,branch)]
# sims.dat$dt.simple_p_given_b_t_arcs[, p.high := sims.dat$p_gt_given_b_t_arcs(branch, t.high, list(v_gt,c_gt,a_gt,d_gt), sims.dat),
#                                     keyby=.(v_gt,c_gt,a_gt,d_gt,branch)]
# setkey(sims.dat$dt.simple_p_given_b_t_arcs, v_gt,c_gt,a_gt,d_gt)
# sims.dat$simple_p_given_b_t_arcs <- function(my.b, my.t, my.gt, sims.dat) {
#   sims.dat$dt.simple_p_given_b_t_arcs[my.gt][my.b == branch][, p.low + (my.t-t.low)/(t.high-t.low) * (p.high-p.low)]
# }

if (F) {
  sims.dat$dt.simple_p_given_b_t_arcs[list(2,2,1,2)]['v' == branch]
  sims.dat$simple_p_given_b_t_arcs('v', .8, list(2,2,1,2), sims.dat)
  x.range = seq(0.5585724, 0.8893428, .01)
  plot(x.range, sapply(x.range, function(my.t) sims.dat$simple_p_given_b_t_arcs('v', my.t, list(2,2,1,2), sims.dat)))
  sims.dat$dt.simple_p_given_b_t_arcs[list(0,2,1,2)]['v' == branch]
  sims.dat$simple_p_given_b_t_arcs('v', .8, list(0,2,1,2), sims.dat)
  plot(x.range, sapply(x.range, function(my.t) sims.dat$simple_p_given_b_t_arcs('v', my.t, list(0,2,1,2), sims.dat)))
}



# dt.sed.og <- fread('~/Google Drive/branch_point_esimates/all_simple_gts.tsv.gz')
# dt.sed.og <- fread('~/Downloads/all_simple_gts.tsv.gz')
dt.sed.og <- fread(args$simple_gts)
dt.sed.og[, .N, lib]

## subset to the desired library [eventually will be able to provide a list of libraries]
if (!is.null(args$libs)) {
  dt.sed <- dt.sed.og[lib == args$libs]
} else {
  dt.sed <- data.table(dt.sed.og)
}
# dt.sed <- dt.sed.og[lib == 'A20398']

## set up simulated data, because I didn't previously fill this in
dt.sed[, lib := 'sim001']
dt.sed.poly <- dt.sed[!(v_gt == c_gt & v_gt == a_gt & v_gt == d_gt)]
dt.sed.poly[, deam53 := c(T,T,F,F,F)]
setnames(dt.sed.poly, c('sed', 'mh_f'), c('sed_gt', 'f_mh'))
dt.sed.poly[, f_mh := f_mh / 99]
dt.sed.poly[, pos := NULL]
dt.sed.qc <- data.table(sim_count = 0, site_count = 1:2000, f_mh = 0, mh = 0, v_gt = 0, c_gt = 0, a_gt = 0, d_gt = 0, sed_gt = 0, 
                        lib = 'sim001', deam53 = c(T,F,F,F))

## look at the data a bit
dt.sed[!(v_gt == c_gt & v_gt == a_gt & v_gt == d_gt), .N, freqs.FLAG]
dt.sed[!(v_gt == c_gt & v_gt == a_gt & v_gt == d_gt), .N, snp.cat.general]
dt.sed[snp.cat.general %like% 'qc', .N, freqs.FLAG]
dt.sed[snp.cat.general %like% 'qc', .N, snp.cat.general]
dt.sed[snp.cat.general %like% 'qc', .N, f_mh]
dt.sed[!snp.cat.general %like% 'qc', .N, f_mh]

## split the dataset into sites that are polymorphic in archaics,
## and qc sites.  require commmmmplete cases.
dt.sed.poly <- dt.sed[!(v_gt == c_gt & v_gt == a_gt & v_gt == d_gt)]
dt.sed.poly <- dt.sed.poly[complete.cases(v_gt,c_gt,a_gt,d_gt,sed_gt,f_mh)]
dt.sed.qc <- dt.sed[snp.cat.general %like% 'qc' & freqs.FLAG == 'invar']
dt.sed.qc[, .N, .(v_gt,c_gt,a_gt,d_gt,sed_gt,f_mh)]
dt.sed.qc <- dt.sed.qc[complete.cases(v_gt,c_gt,a_gt,d_gt,sed_gt,f_mh)]
dt.sed.qc[, .N, .(v_gt,c_gt,a_gt,d_gt,sed_gt,f_mh)]


##############
## just some early fucking around with LL, making sure basic parts work
## this gives pretty similar results as the original composite likelihood stuff, I think?
if (F) {
  #
  my.b = 'anc_1'
  my.t = 1.1
  get_gt_ll <- function(my.b, my.t, dt.sed.poly, sims.dat) {
    cat(my.b, my.t, '..')
    og.t <- my.t
    t.high <- sims.dat$bounds.for.branch(my.b, 'high', sims.dat)
    t.low <- sims.dat$bounds.for.branch(my.b, 'low', sims.dat)
    if (my.t > t.high) my.t <- t.high
    if (my.t < t.low) my.t <- t.low
    dt.sed.poly.p <- dt.sed.poly[, .(p_h_der = sims.dat$p_gt_given_b_t_arcs(my.b, my.t,
                                                                            list(v_gt,c_gt,a_gt,d_gt),
                                                                            sims.dat),
                                     sed_gt, f_mh),
                                 .(v_gt,c_gt,a_gt,d_gt)]
    my.ll <- dt.sed.poly.p[sed_gt == 0, sum(log(1-p_h_der))] + 
      dt.sed.poly.p[sed_gt == 1, sum(log(p_h_der))]
    
    if (my.t == t.high) my.ll <- my.ll - (og.t-t.high)
    if (my.t == t.low)  my.ll <- my.ll - (t.low-og.t)
    return(my.ll)
  }
  get_gt_ll(my.b = 'v', 2, dt.sed.poly, sims.dat)
  get_gt_ll(my.b = 'v', .71, dt.sed.poly, sims.dat)
  
  # optimize(function(my.t) get_gt_ll(my.b = 'v', my.t, dt.sed.poly, sims.dat), c(0,3), maximum = T)
  # a <- sapply(seq(0,3,.1), function(x) get_gt_ll(my.b = 'v', x, dt.sed.poly, sims.dat))
  # a
  # plot(a)
  # for (my.b in sims.dat$branches) {
  #   x <- optimize(function(my.t) get_gt_ll(my.b, my.t, 
  #                                          dt.sed.poly[deam53 == T], 
  #                                          sims.dat), c(0,3), maximum = T)
  #   cat(my.b, x$maximum, x$objective, '\n')
  # }
  
  dt.branch_grid <- foreach(my.branch = sims.dat$branches, .combine = rbind) %do% {
    data.table(my.branch, branch_time = seq(sims.dat$bounds.for.branch(my.branch, 'low', sims.dat),
                                            sims.dat$bounds.for.branch(my.branch, 'high', sims.dat),
                                            .05))
  }
  dt.branch_grid <- dt.branch_grid[branch_time < 2 & my.branch != 'd']
  
  dt.branch_grid[, ll.deam := get_gt_ll(my.branch, branch_time, dt.sed.poly[deam53 == T], sims.dat), .(my.branch, branch_time)]
  dt.branch_grid[, ll.deam.polyN := get_gt_ll(my.branch, branch_time, dt.sed.poly[deam53 == T & !(v_gt == c_gt & v_gt == a_gt)], sims.dat), .(my.branch, branch_time)]
  dt.resample <- dt.sed.poly[deam53 == T & !(v_gt == c_gt & v_gt == a_gt)][sample(.N, replace = T)]
  dt.branch_grid[, ll.deam.polyN.sam0 := get_gt_ll(my.branch, branch_time, dt.resample, sims.dat), .(my.branch, branch_time)]
  dt.resample <- dt.sed.poly[deam53 == T & !(v_gt == c_gt & v_gt == a_gt)][sample(.N, replace = T)]
  dt.branch_grid[, ll.deam.polyN.sam1 := get_gt_ll(my.branch, branch_time, dt.resample, sims.dat), .(my.branch, branch_time)]
  dt.branch_grid[, ll.all := get_gt_ll(my.branch, branch_time, dt.sed.poly, sims.dat), .(my.branch, branch_time)]
  
  ggplot(dt.branch_grid[branch_time < 2 & my.branch != 'd'], aes(x=branch_time, y=ll.deam, color=my.branch)) + geom_point()
  ggplot(dt.branch_grid[branch_time < 2 & my.branch != 'd'], aes(x=branch_time, y=ll.deam.polyN, color=my.branch)) + geom_point()
  ggplot(dt.branch_grid[branch_time < 2 & my.branch != 'd'], aes(x=branch_time, y=ll.deam.polyN.sam0, color=my.branch)) + geom_point()
  ggplot(dt.branch_grid[branch_time < 2 & my.branch != 'd'], aes(x=branch_time, y=ll.deam.polyN.sam1, color=my.branch)) + geom_point()
  ggplot(dt.branch_grid[branch_time < 2 & my.branch != 'd'], aes(x=branch_time, y=ll.all, color=my.branch)) + geom_point()
  
  
  grid_basic_ll <- function(dt.data, sims.dat, n.iter = 10, tag = '') {
    cat(tag, '... \n')
    dt.branch_grid <- foreach(my.branch = sims.dat$branches, .combine = rbind) %do% {
      dt.ret <- data.table(my.branch, branch_time = seq(sims.dat$bounds.for.branch(my.branch, 'low', sims.dat),
                                                        sims.dat$bounds.for.branch(my.branch, 'high', sims.dat),
                                                        .05))
      dt.ret <- dt.ret[branch_time < 2 & my.branch != 'd']
      if (dt.ret[, .N] == 0) return(data.table())
      
      # dt.data <- dt.sed.poly[deam53 == T & !(v_gt == c_gt & v_gt == a_gt)]
      # dt.data <- dt.sed.poly[!(v_gt == c_gt & v_gt == a_gt)]
      
      dt.ret.sam <- foreach(my.iter = 1:n.iter, .combine = rbind) %do% {
        if (my.iter == 0) {
          dt.resample <- dt.data
        } else {
          set.seed(my.iter)
          dt.resample <- dt.data[sample(.N, replace = T)]
        }
        dt.ret[, .(ll = get_gt_ll(my.branch, branch_time, dt.resample, sims.dat),
                   my.iter = as.character(my.iter)), 
               .(my.branch, branch_time)]
      }
      dt.ret.sam[, tag := tag]
      dt.ret.sam
    }
  }
  
  
  dt.branch_grid <- grid_basic_ll(dt.sed.poly[deam53 == T], sims.dat, n.iter = 10)
  dt.branch_grid <- 
    rbind(grid_basic_ll(dt.sed.poly, sims.dat, n.iter = 10, tag = 'all'),
          grid_basic_ll(dt.sed.poly[!(v_gt == c_gt & v_gt == a_gt)], sims.dat, n.iter = 10, tag = 'all_polyN'),
          grid_basic_ll(dt.sed.poly[deam53 == T], sims.dat, n.iter = 10, tag = 'deam'),
          grid_basic_ll(dt.sed.poly[deam53 == T & !(v_gt == c_gt & v_gt == a_gt)], sims.dat, n.iter = 10, tag = 'deam_polyN'))
  dt.branch_grid <- 
    rbind(dt.branch_grid,
          grid_basic_ll(dt.sed.poly[!(v_gt == c_gt & v_gt == a_gt) | (v_gt == c_gt & v_gt == a_gt & v_gt == 2 & d_gt == 0)],
                        sims.dat, n.iter = 10, tag = 'all_poly_N_2_D_0'),
          grid_basic_ll(dt.sed.poly[!(v_gt == c_gt & v_gt == a_gt) | (v_gt == c_gt & v_gt == a_gt & v_gt == 0 & d_gt == 2)],
                        sims.dat, n.iter = 10, tag = 'all_poly_N_0_D_2'),
          grid_basic_ll(dt.sed.poly[!(v_gt == c_gt & v_gt == a_gt) | (v_gt == c_gt & v_gt == a_gt & v_gt == 0 & d_gt == 1)],
                        sims.dat, n.iter = 10, tag = 'all_poly_N_0_D_1'), fill = T)
  dt.branch_grid <- 
    rbind(dt.branch_grid,
          grid_basic_ll(dt.sed.poly[deam53 == T][!(v_gt == c_gt & v_gt == a_gt) | (v_gt == c_gt & v_gt == a_gt & v_gt == 2 & d_gt == 0)],
                        sims.dat, n.iter = 10, tag = 'deam_poly_N_2_D_0'),
          grid_basic_ll(dt.sed.poly[deam53 == T][!(v_gt == c_gt & v_gt == a_gt) | (v_gt == c_gt & v_gt == a_gt & v_gt == 0 & d_gt == 2)],
                        sims.dat, n.iter = 10, tag = 'deam_poly_N_0_D_2'),
          grid_basic_ll(dt.sed.poly[deam53 == T][!(v_gt == c_gt & v_gt == a_gt) | (v_gt == c_gt & v_gt == a_gt & v_gt == 0 & d_gt == 1)],
                        sims.dat, n.iter = 10, tag = 'deam_poly_N_0_D_1'), fill = T)
  dt.sed.poly[(v_gt == c_gt & v_gt == a_gt), .N, keyby=.(v_gt, d_gt, sed_gt)]
  
  dt.branch_grid[, ll.norm := -ll/max(ll), .(my.iter, tag)]
  dt.branch_grid.max <- dt.branch_grid[, .SD[ll == max(ll)], .(my.iter, tag)]
  ggplot(dt.branch_grid, aes(x=branch_time, y=ll.norm,
                             color=my.branch,
                             group=interaction(my.branch,my.iter))) +
    geom_line() +
    geom_point(data = dt.branch_grid.max) + facet_wrap(~tag)
}

###############
###############
## it looks like there's an effect where when you add sites that are fixed in N and diff btwn D & N,
## it shifts the timing down.  particularly for fixed derived sites in N, you have more ancestral alleles than for e.g. fixed ancestral alleles in N
## want to see if that's expected
if (F) {
  my.branch.dt.p = data.table(my.branch = sims.dat$branches)
  my.branch.dt.p = my.branch.dt.p[, expand.grid(branch_time = seq(sims.dat$bounds.for.branch(my.branch, 'low', sims.dat),
                                                                  sims.dat$bounds.for.branch(my.branch, 'high', sims.dat),
                                                                  .01),
                                                v_gt = c(0,1,2),
                                                c_gt = c(0,1,2),
                                                a_gt = c(0,1,2),
                                                d_gt = c(0,1,2)),
                                  my.branch]
  setkey(my.branch.dt.p, v_gt, c_gt, a_gt, d_gt)
  my.branch.dt.p[, p := sims.dat$p_gt_given_b_t_arcs(as.character(my.branch), branch_time,
                                                     list(v_gt, c_gt, a_gt, d_gt),
                                                     sims.dat),
                 .(v_gt, c_gt, a_gt, d_gt, branch_time, my.branch)]
  my.branch.dt.p[, p.norm := (p-min(p)) / (max(p)-min(p)), by=.(v_gt, c_gt, a_gt, d_gt)]
  
  dt.sed.poly.prop = dt.sed.poly[, .(p.sed=sum(sed_gt==1)/.N, .N), keyby=.(v_gt, c_gt, a_gt, d_gt)]
  dt.sed.poly.prop <- dt.sed.poly.prop[N > 10]
  my.branch.dt.p <- my.branch.dt.p[dt.sed.poly.prop]
  my.branch.dt.p[, tag := sprintf('V%d_C%d_A%d_D%d__%dsites', v_gt, c_gt, a_gt, d_gt, N)]
  my.branch.dt.p[, base.tag := sprintf('V%d_C%d_A%d_D%d', v_gt, c_gt, a_gt, d_gt)]
  tag.levs = unique(my.branch.dt.p[, .('x'), keyby=.(N, tag)])[,tag]
  my.branch.dt.p[, tag.fac := factor(tag, levels = tag.levs)]
  my.branch.dt.p <- my.branch.dt.p[my.branch != 'd' & branch_time < 2]
  
  my.branch.dt.p[, p.sed.norm := (p.sed-min(p)) / (max(p)-min(p)), by=.(v_gt, c_gt, a_gt, d_gt)]
  ## this should be binom conf intervals, and plotted w/ a vertical line
  my.branch.dt.p[, p.optim := abs(p-p.sed) == min(abs(p-p.sed)), by=.(v_gt, c_gt, a_gt, d_gt)]
  
  basic.sites = c('V0_C0_A0_D2',
                  'V2_C2_A2_D0',
                  'V2_C0_A0_D0',
                  'V1_C0_A0_D0')
  
  ggplot(my.branch.dt.p, aes(x=branch_time, y=p, color=my.branch)) + geom_line() +
    facet_wrap(~tag)
  ggplot(my.branch.dt.p[N>400], aes(x=branch_time, y=p, color=my.branch)) + geom_line() +
    facet_wrap(~tag.fac)
  ggplot(my.branch.dt.p[N>400], aes(x=branch_time, y=p, color=my.branch)) + 
    geom_hline(data=my.branch.dt.p[p.optim == T & N>400], aes(yintercept=p.sed), color='grey') +
    geom_point(data=my.branch.dt.p[p.optim == T & N>400], aes(y=p.sed)) +
    geom_line() + facet_wrap(~tag.fac)
  pdf('~/Google Drive/branch_point_esimates/mez2_estim_site_cats_debug.pdf', width=7, height=5)
  ggplot(my.branch.dt.p[base.tag %in% basic.sites], aes(x=branch_time, y=p, color=my.branch)) + 
    geom_hline(data=my.branch.dt.p[p.optim == T & base.tag %in% basic.sites], aes(yintercept=p.sed), color='grey') +
    geom_point(data=my.branch.dt.p[p.optim == T & base.tag %in% basic.sites], aes(y=p.sed)) +
    geom_line() + facet_wrap(~tag.fac) +
    ylab('Probability of der allele in test sample')
  dev.off()
  ggplot(my.branch.dt.p, aes(x=branch_time, y=p.norm, color=my.branch)) + geom_line() +
    geom_point(data=my.branch.dt.p[p.optim == T], aes(y=p.sed.norm)) +
    facet_wrap(~tag.fac)
}


##########
## 





#################################
## the original code for calculating a composite likelihood
#################################

# test.gts <- fread('~/Google Drive/soil_dna_capture/all_simple_gts.deam.tsv.gz')
# test.gts <- fread('output_v4/all_simple_gts.deam.tsv.gz')
# test.gts <- fread('zcat output_v4/all_simple_gts.deam.tsv.gz | head -n1000')

test.gts

dt.sims.p = fread('~/Google Drive/soil_dna_capture/simulate_demog/dt.sims.p.simfiles004.txt')
dt.sims.p.coarse = dt.sims.p[!endsWith(paste(time), '01') | endsWith(paste(time), '001')]


# demog_composite_liks <- function(dt.sims.p, test.gts, model.contam = NULL) {
  
  ## dt.sims.p should be of this format:
  # > dt.sims.p
  #    v c a sed    time branch      p
  # 1: 0 0 0   0 0.89999      a 0.9928
  # 2: 0 0 0   0 0.89999      c 0.9887
  # 3: 0 0 0   0 0.89999      v 0.9891
  # 4: 0 0 0   0 0.91000      a 0.9920
  
  ## test.gts should be
  # test.gts <- fread('~/Google Drive/soil_dna_capture/all_simple_gts.deam.tsv.gz')
  #       libname v c a sed deam53              lib  pct_ref pct_ref.low
  # 1:     A16036 0 0 0   0   TRUE           A16036 1.000000   0.6456696
  # 2:     A16036 2 2 2   0   TRUE           A16036 1.000000   0.6456696
  # 3:     A16036 0 2 0   0   TRUE           A16036 1.000000   0.6456696
  # 4:     A16036 0 0 0   0   TRUE           A16036 1.000000   0.6456696
  ## although I think only v,c,a,libname,lib,sed are used?
  

if (F) {
  libs.150 <- test.gts[, .N, lib][N > 150, lib]
  cat('a.s.lib\n')
  a.s.lib <- test.gts[lib %in% libs.150, demog_liks(dt.sims.p, .SD, model.contam = T)$full.tree, .(lib,libname)]
  cat('a.s.libname\n')
  a.s.libname <- test.gts[lib %in% libs.150, demog_liks(dt.sims.p, .SD, model.contam = T)$full.tree, libname]
  date()
}
 



dt.sed.analysis = rbind(dt.sed.qc,dt.sed.poly)
# dt.sed.analysis = rbind(dt.sed.qc,dt.sed.poly[!(v_gt == c_gt & v_gt == a_gt)])
dt.sed.analysis[, rg := paste0(lib, '_deam_', deam53)]

if (args$site_cat == 'tv') dt.sed.analysis <- dt.sed.analysis[freqs.FLAG != 'transi']
if (args$site_cat == 'polyN') dt.sed.analysis <- dt.sed.analysis[!(v_gt == c_gt & v_gt == a_gt)]

if (F) {
  my.ret = sed_EM(dt.sed.analysis[freqs.FLAG != 'transi', .SD[sample(.N, 100)]], sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 10, ll.converge = 1000)
  my.ret = sed_EM_allbranch(dt.sed.analysis[freqs.FLAG != 'transi', .SD[sample(.N, 100)]], sims.dat, err_rate = 0.001, max.iter = 2)
  sed_EM(dt.sed.analysis[freqs.FLAG != 'transi'], sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 10)
  # set branchtime
  sed_EM(dt.sed.analysis[freqs.FLAG != 'transi'], sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 10, set.branchtime = .8)
  sed_EM(dt.sed.analysis[freqs.FLAG != 'transi' & !(v_gt == c_gt & v_gt == a_gt)], sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 10)
  sed_EM(dt.sed.analysis, sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 10)
}

## actually estimate full model - contam, faunal, and branch point
my.ret = sed_EM_allbranch(dt.sed.analysis, sims.dat, err_rate = 0.001, max.iter = args$num_iters, ll.converge = 1)

dt.ret <- rbind(ll_ret_to_dt(my.ret$a, args$site_cat),
                ll_ret_to_dt(my.ret$v, args$site_cat),
                ll_ret_to_dt(my.ret$c, args$site_cat),
                ll_ret_to_dt(my.ret$anc_1, args$site_cat),
                ll_ret_to_dt(my.ret$anc_2, args$site_cat),
                ll_ret_to_dt(my.ret$anc_3, args$site_cat),
                ll_ret_to_dt(my.ret$d, args$site_cat))
print(dt.ret)
fwrite(dt.ret, sprintf('%s_estim_branch.ll.txt', args$prefix), sep = '\t')
saveRDS(my.ret, sprintf('%s_estim_branch.ll.RDS', args$prefix))


################################
##############
## temporarily have this function close by, so I can change it to test things.

#####
## called in q_t and calc_gamma, currently
# p_h_given_b_t <- function(dt.sed, sims.dat, gt, branch, branch_time) {
#   ## p_gt_given_b_t_arcs gives p(H==der | etc)
#   
#   # dt.p <- dt.sed[, .(p = sims.dat$p_gt_given_b_t_arcs(branch, branch_time, list(v_gt,c_gt,a_gt,d_gt), sims.dat)),
#   #                keyby=.(v_gt,c_gt,a_gt,d_gt)]
#   # dt.p <- dt.sed[, .(p = sims.dat$simple_p_given_b_t_arcs(branch, branch_time, list(v_gt,c_gt,a_gt,d_gt), sims.dat)),
#   #                keyby=.(v_gt,c_gt,a_gt,d_gt)]
#   dt.p <- dt.sed[, .(p = branch_time),
#                  keyby=.(v_gt,c_gt,a_gt,d_gt)]
#   ## this [CURRENTLY] only happens for QC sites. 
#   ## the other sites we are considering are polymorphic in archaics
#   ## WOULD HAVE TO CHANGE THIS if I add sites that are fixed in archaics but poly in MH
#   ## - this is actually a lot of sites, so it would add more information
#   dt.p[v_gt == 2 & c_gt == 2 & a_gt == 2 & d_gt == 2, p := .999]
#   ## at this point, p is p(H==der).  We want p(H==gt).  Flip if necessary.
#   if (gt == 0) dt.p[, p := 1-p]
#   
#   # print(plyr::join(dt.sed[, .(v_gt,c_gt,a_gt,d_gt)], dt.p, by=c('v_gt','c_gt','a_gt','d_gt'))[, p])
#   
#   ## unlike merge, plyr::join maintains the order of dt.sed
#   plyr::join(dt.sed[, .(v_gt,c_gt,a_gt,d_gt)], dt.p, by=c('v_gt','c_gt','a_gt','d_gt'))[, p]
# }


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
  
  ## do this on the 'real' simulated data?
  ## time bash scrm9.mod.v.sh 10000 .058 0.7601
  ## I think the simulations for this come from scrm1, but those bootstraps also look suspiciously similar?
  dt.sed.analysis = rbind(dt.sed.qc,dt.sed.poly)
  # dt.sed.analysis = rbind(dt.sed.qc,dt.sed.poly[!(v_gt == c_gt & v_gt == a_gt)])
  dt.sed.analysis[, rg := paste0(lib, '_deam_', deam53)]

  eval.ret.analysis.simple <- eval_sed_t_and_mh(dt.sed.analysis, 'simple')
  plot_eval_sed_t_and_mh(eval.ret.analysis.simple)

  eval.ret.analysis.full <- eval_sed_t_and_mh(dt.sed.analysis, 'full')
  plot_eval_sed_t_and_mh(eval.ret.analysis.full)
  
  date()
  
  eval.ret.analysis.simple.i500 <- eval_sed_t_and_mh(dt.sed.analysis, 'simple', max.iter = 500)
  plot_eval_sed_t_and_mh(eval.ret.analysis.simple.i500)
  
  eval.ret.analysis.full.i500 <- eval_sed_t_and_mh(dt.sed.analysis, 'full', max.iter = 500)
  plot_eval_sed_t_and_mh(eval.ret.analysis.full.i500)
  
  date()
  
  
  ####################
  ####################
  ####################
  ####################
  ####################
  ####################
  
  
  
  
  
  
  
    
  ## constrain mh and faunal = 0, just estimate time with EM
  x.em = sed_EM(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 10, set.faunal_prop = 0, set.mh_contam = 0, p_h_method = 'bt2')
  
  ## the ll goes up and down a bit in the trace.. just not 100% stable?
  # plot(tail(x.em$man.ll.trace,9))
  # plot(tail(x.em$ll.trace,9))
  
  ## do "grid" search, just estimating time  
  dt.sed.poly.tst.gridll = sed_grid_search(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, p_h_method = 'bt2',
                                           bins.mh_contam = 1, bins.faunal_prop = 1, bins.t = 10, nsteps = 5,
                                           range.mh_contam = c(0,0), range.faunal_prop = c(0,0), range.t = c(.55,.85))
  ## plot them, they match!
  ggplot(dt.sed.poly.tst.gridll, aes(x=my.t, y=ll)) + geom_point() +
    geom_point(aes(x=x.em$branchtime, y=x.em$man.max.ll), color='red', pch='x', size=10) +
    xlab('branch time estimate grid vs EM (red x)')
  
  
  x.em.all_params.n10 = sed_EM(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 10, set.faunal_prop = 0, p_h_method = 'bt')
  x.em.all_params.n20 = sed_EM(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 20, set.faunal_prop = 0, p_h_method = 'bt')
  x.em.all_params.n100 = sed_EM(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 100, set.faunal_prop = 0, p_h_method = 'bt')
  
  # plot(tail(x.em.all_params$man.ll.trace,9))
  # plot(tail(x.em.all_params$ll.trace,9))
  
  dt.sed.poly.tst.gridll.all_params = sed_grid_search(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, p_h_method = 'bt',
                                                      bins.mh_contam = 10, bins.faunal_prop = 1, bins.t = 10, nsteps = 4,
                                                      range.mh_contam = c(0,.2), range.faunal_prop = c(0,0), range.t = c(.55,.85))
  em.mismatch = dt.sed.poly.tst.gridll.all_params[is.max == T, ll] - x.em.all_params.n100$man.max.ll
  ggplot(dt.sed.poly.tst.gridll.all_params[step.x > 0],
         aes(x=my.t, y=mh_contam, color=ll)) +
    geom_point(aes(x=.75, y=0), color='black', pch='x', size=10) +
    geom_point() +
    geom_point(data=dt.sed.poly.tst.gridll.all_params[max(ll) - ll < em.mismatch & step.x > 0], color='red') +
    geom_point(data=dt.sed.poly.tst.gridll.all_params[max(ll) - ll < em.mismatch/10 & step.x > 0], color='green') +
    geom_point(aes(x=x.em.all_params.n10$branchtime, y=x.em.all_params.n10$dt.theta$mh_contam), color='red', pch='x', size=5) +
    geom_point(aes(x=x.em.all_params.n20$branchtime, y=x.em.all_params.n20$dt.theta$mh_contam), color='red', pch='x', size=7) +
    geom_point(aes(x=x.em.all_params.n100$branchtime, y=x.em.all_params.n100$dt.theta$mh_contam), color='red', pch='x', size=10) +
    xlab('branch time estimate grid (green/red dots) vs EM (red x)') + ylab('mh contam estimate') +
    ggtitle(sprintf('ll of red points within %g of EM: grid and EM effectively match, but EM doesn\'t quite optimize', em.mismatch))
  #
  
  
  
  
  
  
  
  
    
  
  
  ## p_h_der is equal to bt
  x.em.f_mh = foreach(my.f_mh = seq(0,.9,.1), .combine = rbind) %:%
    foreach(my.p_h_der = seq(0.1,.9,.1), .combine = rbind) %do% {
      dt.sed.poly.tst <- data.table(v_gt=2, c_gt=2, a_gt=0, d_gt=0, f_mh=my.f_mh, p_h_der = my.p_h_der, rg = 'hey_rg', x = 1:1000)
    dt.sed.poly.tst[, sed_gt := sample(c(0,1), .N, prob = c(1-p_h_der,p_h_der), replace=T), p_h_der]
    
    err.sites = dt.sed.poly.tst[, sample(.N, .N*0.001)]
    dt.sed.poly.tst[err.sites, .N, sed_gt]
    dt.sed.poly.tst[err.sites, sed_gt := sample(c(0,1), .N, prob=c(.5,.5), replace=T)]
    dt.sed.poly.tst[err.sites, .N, sed_gt]
    #
    
    # x.em = sed_EM(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 6, set.mh_contam = 0, set.faunal_prop = 0, p_h_method = 'bt')
    x.em = sed_EM(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 3, p_h_method = 'bt')
    data.table(est.bt = x.em$branchtime, x.em$dt.theta, my.f_mh, my.p_h_der)
  }
  ggplot(x.em.f_mh, aes(x=my.f_mh, y=mh_contam, color=factor(my.p_h_der))) + geom_line() + ylim(0,NA)
  ggplot(x.em.f_mh, aes(x=my.p_h_der, y=est.bt, color=factor(my.f_mh))) + geom_line() + geom_abline(slope=1)
  #
  
  ## p_h_der is equal to bt
  x.em.f_mh = foreach(my.f_mh_s1 = c(.2,.5,1), .combine = rbind) %:%
    foreach(my.t = seq(0.55,.85,.1), .combine = rbind) %:%
    foreach(my.mh_contam = seq(0.0,0.2,0.05), .combine = rbind) %do% {
      nsites = 1002
      dt.sed.poly.tst <- data.table(v_gt=2, c_gt=0:2, a_gt=0, d_gt=0, f_mh=rbeta(nsites, my.f_mh_s1, 5), rg = 'hey_rg', x = 1:nsites)
      dt.sed.poly.tst[, p_h_der := my.t]
      dt.sed.poly.tst[, sed_gt := sample(c(0,1), .N, prob = c(1-p_h_der,p_h_der), replace=T), p_h_der]
      
      dt.sed.poly.tst[, .(sum(sed_gt) / .N, .N), .(p_h_der)]
      
      mh.sites = dt.sed.poly.tst[, sample(.N, .N*my.mh_contam)]
      dt.sed.poly.tst[mh.sites, .N, keyby=sed_gt]
      dt.sed.poly.tst[mh.sites, sed_gt := sample(c(0,1), .N, prob=c(1-f_mh,f_mh), replace=T), f_mh]
      dt.sed.poly.tst[mh.sites, .(sum(sed_gt)/.N, mean(f_mh))]
      #
      err.sites = dt.sed.poly.tst[, sample(.N, .N*0.001)]
      dt.sed.poly.tst[err.sites, .N, sed_gt]
      dt.sed.poly.tst[err.sites, sed_gt := sample(c(0,1), .N, prob=c(.5,.5), replace=T)]
      dt.sed.poly.tst[err.sites, .N, sed_gt]
      #
      
      # x.em = sed_EM(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 6, set.mh_contam = 0, set.faunal_prop = 0, p_h_method = 'bt')
      x.em = sed_EM(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 3, p_h_method = 'bt')
      data.table(est.bt = x.em$branchtime, x.em$dt.theta, my.f_mh_s1, my.t, my.mh_contam)
    }
  ggplot(x.em.f_mh, aes(x=my.mh_contam, y=mh_contam, color=factor(my.t))) + geom_line() + ylim(0,NA) + geom_abline(slope=1) + facet_wrap(~my.f_mh_s1)
  ggplot(x.em.f_mh, aes(x=my.t, y=est.bt, color=factor(my.mh_contam))) + geom_line() + geom_abline(slope=1) + facet_wrap(~my.f_mh_s1)
  ggplot(x.em.f_mh, aes(x=my.t, y=mh_contam, color=factor(my.f_mh_s1))) + geom_line() + geom_abline(slope=1) + facet_wrap(~my.mh_contam)
  
  
  ## p_h_der is equal to bt / (c_gt+1)
  ## that way the derived frequency depends somewhat on the genotype, and hopefully can be better separated from mh_contam
  x.em.f_mh2 = foreach(my.f_mh_s1 = c(.2,.5,1), .combine = rbind) %:%
    foreach(my.t = seq(0.55,.85,.1), .combine = rbind) %:%
    foreach(my.mh_contam = seq(0.0,0.2,0.05), .combine = rbind) %do% {
      nsites = 1002
      dt.sed.poly.tst <- data.table(v_gt=2, c_gt=0:2, a_gt=0, d_gt=0, f_mh=rbeta(nsites, my.f_mh_s1, 5), rg = 'hey_rg', x = 1:nsites)
      dt.sed.poly.tst[, p_h_der := my.t / (c_gt+1)]
      dt.sed.poly.tst[, sed_gt := sample(c(0,1), .N, prob = c(1-p_h_der,p_h_der), replace=T), p_h_der]
      
      dt.sed.poly.tst[, .(sum(sed_gt) / .N, .N), .(p_h_der)]
      
      mh.sites = dt.sed.poly.tst[, sample(.N, .N*my.mh_contam)]
      dt.sed.poly.tst[mh.sites, .N, keyby=sed_gt]
      dt.sed.poly.tst[mh.sites, sed_gt := sample(c(0,1), .N, prob=c(1-f_mh,f_mh), replace=T), f_mh]
      dt.sed.poly.tst[mh.sites, .(sum(sed_gt)/.N, mean(f_mh))]
      #
      err.sites = dt.sed.poly.tst[, sample(.N, .N*0.001)]
      dt.sed.poly.tst[err.sites, .N, sed_gt]
      dt.sed.poly.tst[err.sites, sed_gt := sample(c(0,1), .N, prob=c(.5,.5), replace=T)]
      dt.sed.poly.tst[err.sites, .N, sed_gt]
      #

      # x.em = sed_EM(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 6, set.mh_contam = 0, set.faunal_prop = 0, p_h_method = 'bt')
      x.em = sed_EM(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 3, p_h_method = 'bt2')
      
      dt.sed.poly.tst.gridll = sed_grid_search(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, p_h_method = 'bt2',
                                               bins.mh_contam = 10, bins.faunal_prop = 1, bins.t = 10,
                                               range.mh_contam = c(0,.2), range.faunal_prop = c(0,0), range.t = c(.55,.85))
      setnames(dt.sed.poly.tst.gridll, paste0('grid.est.',names(dt.sed.poly.tst.gridll)))
      setnames(x.em$dt.theta, paste0('em.est.',names(x.em$dt.theta)))
      
      data.table(em.est.ll = x.em$max.ll, em.est.bt = x.em$branchtime, x.em$dt.theta, my.f_mh_s1, my.t, my.mh_contam, dt.sed.poly.tst.gridll)
    }
  x.em.f_mh2[, iter.tag := paste(my.f_mh_s1, my.t, my.mh_contam)]
  ggplot(x.em.f_mh2, aes(x=my.mh_contam, y=mh_contam, color=factor(my.t))) + geom_line() + ylim(0,NA) + geom_abline(slope=1) + facet_wrap(~my.f_mh_s1)
  ggplot(x.em.f_mh2, aes(x=my.t, y=est.bt, color=factor(my.mh_contam))) + geom_line() + geom_abline(slope=1) + facet_wrap(~my.f_mh_s1)
  ggplot(x.em.f_mh2, aes(x=my.t, y=mh_contam, color=factor(my.f_mh_s1))) + geom_line() + geom_abline(slope=1)

  ggplot(x.em.f_mh2[grid.est.is.max == T & faunal_prop < .001], aes(mh_contam, grid.est.mh_contam)) + geom_point() + geom_abline(slope=1)
  my.iter.tag = x.em.f_mh2[grid.est.is.max == T & faunal_prop < .001, iter.tag][1]
  ggplot(x.em.f_mh2[iter.tag == my.iter.tag], aes(grid.est.my.t, grid.est.mh_contam, color=grid.est.ll)) + geom_point() +
    geom_point(aes(est.bt, mh_contam), color='red') +
    geom_point(aes(my.t, my.mh_contam), color='green')
    

  ## shuffle freqs
  ## p_h_der is equal to bt / (c_gt+1)
  x.em.f_mh2.shuf = foreach(my.f_mh_s1 = c(.2,.5,1), .combine = rbind) %:%
    foreach(my.t = seq(0.55,.85,.1), .combine = rbind) %:%
    foreach(my.mh_contam = seq(0.0,0.2,0.05), .combine = rbind) %do% {
      nsites = 1002
      dt.sed.poly.tst <- data.table(v_gt=2, c_gt=0:2, a_gt=0, d_gt=0, f_mh=rbeta(nsites, my.f_mh_s1, 5), rg = 'hey_rg', x = 1:nsites)
      dt.sed.poly.tst[, p_h_der := my.t / (c_gt+1)]
      dt.sed.poly.tst[, sed_gt := sample(c(0,1), .N, prob = c(1-p_h_der,p_h_der), replace=T), p_h_der]
      
      dt.sed.poly.tst[, .(sum(sed_gt) / .N, .N), .(p_h_der)]
      
      mh.sites = dt.sed.poly.tst[, sample(.N, .N*my.mh_contam)]
      dt.sed.poly.tst[mh.sites, .N, keyby=sed_gt]
      dt.sed.poly.tst[mh.sites, sed_gt := sample(c(0,1), .N, prob=c(1-f_mh,f_mh), replace=T), f_mh]
      dt.sed.poly.tst[mh.sites, .(sum(sed_gt)/.N, mean(f_mh))]
      #
      err.sites = dt.sed.poly.tst[, sample(.N, .N*0.001)]
      dt.sed.poly.tst[err.sites, .N, sed_gt]
      dt.sed.poly.tst[err.sites, sed_gt := sample(c(0,1), .N, prob=c(.5,.5), replace=T)]
      dt.sed.poly.tst[err.sites, .N, sed_gt]
      #
      
      dt.sed.poly.tst[, sed_gt := sample(sed_gt)]
      
      # x.em = sed_EM(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 6, set.mh_contam = 0, set.faunal_prop = 0, p_h_method = 'bt')
      x.em = sed_EM(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 3, p_h_method = 'bt2')
      data.table(est.bt = x.em$branchtime, x.em$dt.theta, my.f_mh_s1, my.t, my.mh_contam)
    }
  ggplot(x.em.f_mh2.shuf, aes(x=my.f_mh, y=mh_contam, color=factor(my.t))) + geom_line() + ylim(0,NA) + xlab('MH allele frequency (all sites)') + ylab('Estimated MH contam [actual=0]')
  ggplot(x.em.f_mh2.shuf, aes(x=my.t, y=est.bt, color=factor(my.f_mh))) + geom_line() + geom_abline(slope=1)
  ggplot(x.em.f_mh2.shuf, aes(x=my.t, y=mh_contam, color=factor(my.f_mh))) + geom_line() + geom_abline(slope=1)
  #
  
  
  
  x.manll.dt = foreach(mh_contam = seq(0,.2,.02), .combine = rbind) %:%
    foreach(faunal_prop = seq(0,.2,.02), .combine = rbind) %:%
    foreach(my.t = seq(0.5, .9, 0.02), .combine = rbind) %do% {
      cat(mh_contam, faunal_prop, my.t)
      ll = calc_manual_lik(dt.sed.poly.tst,
                           mh_contam = mh_contam,
                           faunal_prop = faunal_prop,
                           err_rate = 0.001,
                           sims.dat = sims.dat, 
                           my.branch = 'v', 
                           branch_time = my.t,
                           method = 'bt')
      data.table(my.t, mh_contam, faunal_prop, ll)
    }
  x.manll.dt[, is.max := ll == max(ll), .(mh_contam, faunal_prop)]
  x.manll.dt[, is.max := ll == max(ll), .(mh_contam, faunal_prop)]
  # ggplot(x.manll2.dt, aes(x=my.t, y=ll, color=err_rate, group=err_rate)) + geom_line()
  ggplot(x.manll.dt, aes(x=my.t, y=ll, color=mh_contam, group=mh_contam)) + geom_line() +
    geom_point(data=x.manll.dt[is.max==T]) + 
    geom_point(data=x.manll.dt[ll == max(ll)], color='red') + 
    facet_wrap(~faunal_prop) +
    geom_vline(xintercept = .75) +
    ylim(-57000,-56344)
  ggplot(x.manll.dt[faunal_prop == 0], aes(x=my.t, y=ll, color=factor(mh_contam), group=mh_contam)) + geom_line() +
    geom_point(data=x.manll.dt[is.max==T & faunal_prop == 0]) + 
    # geom_point(data=x.manll.dt[ll == max(ll)], color='red') + 
    facet_wrap(~faunal_prop) +
    geom_vline(xintercept = .75) +
    ylim(-57000,-56344)
  #

  
  
  #
  
  
  
  
  
  
  
  
  
  ### older tests
  
  dt.sed.poly.tst <- dt.sed.poly[, .(v_gt, c_gt, a_gt, d_gt, f_mh, deam53)]
  
  dt.sed.poly.tst[, .(.N),
                  keyby=.(v_gt,c_gt,a_gt,d_gt)] %>% setkey(N) %>% print

  dt.sed.poly.tst <- data.table(v_gt=2, c_gt=2, a_gt=0, d_gt=0, f_mh=.5, x=seq(100))
  dt.sed.poly.tst <- rbind(data.table(v_gt=2, c_gt=2, a_gt=0, d_gt=0, f_mh=.5, x=seq(100000)),
                           data.table(v_gt=2, c_gt=2, a_gt=1, d_gt=0, f_mh=.5, x=seq(100000)))
  dt.sed.poly.tst <- dt.sed.poly.tst[sample(.N)]
  dt.sed.poly.tst[, rg := paste0('tst_deam_', F)]

  
  
  # dt.sed.poly.tst[, p_h_der := sims.dat$p_gt_given_b_t_arcs('v', .6, 
  #                                                           list(v_gt,c_gt,a_gt,d_gt), 
  #                                                           sims.dat),
  #                 .(v_gt,c_gt,a_gt,d_gt)]
  dt.sed.poly.tst[, p_h_der := sims.dat$simple_p_given_b_t_arcs('v', .7,
                                                            list(v_gt,c_gt,a_gt,d_gt),
                                                            sims.dat),
                  .(v_gt,c_gt,a_gt,d_gt)]
  dt.sed.poly.grid <- data.table(v_gt=2, c_gt=2, a_gt=0, d_gt=0,
                                 my.t = seq(.56,.88,.03))
  dt.sed.poly.grid[, .(p_h_der = sims.dat$simple_p_given_b_t_arcs('v', my.t,
                                                                list(v_gt,c_gt,a_gt,d_gt),
                                                                sims.dat)),
                  .(v_gt,c_gt,a_gt,d_gt,my.t)]
  dt.sed.poly.tst[, p_h_der := .75,
                  .(v_gt,c_gt,a_gt,d_gt)]
  dt.sed.poly.tst[, sed_gt := sample(c(0,1), .N, prob = c(1-p_h_der,p_h_der), replace=T),
                  .(v_gt,c_gt,a_gt,d_gt,p_h_der)]
  
  dt.sed.poly.tst <- data.table(v_gt=2, c_gt=2, a_gt=0, d_gt=0, f_mh=.5, 
                                sed_gt=c(rep(0,1), rep(1,3)), p_h_der = .75)
  dt.sed.poly.tst[, rg := paste0('tst_deam_', F)]
  
  
  # setkey(dt.sed.poly.tst, sed_gt)
  dt.sed.poly.tst[, .(sum(sed_gt)/.N, .N),
                  keyby=.(v_gt,c_gt,a_gt,d_gt,p_h_der)]
  # ggplot(dt.sed.poly.tst[, .(sum(sed_gt)/.N, .N),
  #                 .(v_gt,c_gt,a_gt,d_gt,p_h_der)],
  #        aes(x=N, y=p_h_der - V1)) + geom_point()
  # dt.sed.poly.tst <-
  #   dt.sed.poly.tst[v_gt==2 & c_gt==2 & a_gt==2 & d_gt==0]

  # sed_EM(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 1, set.mh_contam = 0, set.faunal_prop = 0, set.branchtime = .8)
  # sapply(sims.dat$branch.bounds[branch == 'v', seq(t.low, t.high, .1)],
  #        function(my.t) sed_EM(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 1, 
  #                              set.mh_contam = 0, set.faunal_prop = 0, set.branchtime = my.t)$max.ll)
  # x.em = sed_EM(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 2, set.mh_contam = 0, set.faunal_prop = 0, set.branchtime = 'grid')
  # 
  x.em = sed_EM(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 6, set.mh_contam = 0, set.faunal_prop = 0, p_h_method = 'bt')
  # plot(tail(x.em$ll.trace,4))
  # x.range = seq(x.em$branchtime-.05, min(x.em$branchtime+.05,0.8893428), .01)
  x.range = seq(0.5585724, 0.8893428, .03)
  # x.range = seq(0.5585724, .62, .003)
  x.range = seq(.74, .76, .001)
  x = sapply(x.range,
         function(my.t) sed_EM(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 1, 
                               set.mh_contam = 0, set.faunal_prop = 0, set.branchtime = my.t)$max.ll)
  true.der = dt.sed.poly.tst[, .(sum(sed_gt)/.N, .N),
                  keyby=.(v_gt,c_gt,a_gt,d_gt,p_h_der)][, V1]
  dt.sed.poly.grid[, .(true.der, p_h_der = sims.dat$simple_p_given_b_t_arcs('v', my.t,
                                                                  list(v_gt,c_gt,a_gt,d_gt),
                                                                  sims.dat)),
                   .(v_gt,c_gt,a_gt,d_gt,my.t)][, .(my.t, x.em$branchtime, abs(true.der-p_h_der) == min(abs(true.der-p_h_der)),
                                                    true.der, p_h_der, v_gt,c_gt,a_gt,d_gt)]
  plot(x.range, x)
  abline(v=x.em$branchtime, col='blue')
  abline(v=x.range[x == max(x)], col='red')
  
  
  
  
  #########
  
  all.rg = dt.sed.poly.tst[, unique(rg)]
  o.t = list(maximum = 0.77, objective = 0)
  
  t.last = o.t$maximum
  iter.ll = update_gamma(dt.sed.poly.tst, all.rg = all.rg,
                         sims.dat = sims.dat,
                         dt.theta = data.table(rg = all.rg,
                                               mh_contam = 0,
                                               faunal_prop = 0),
                         err_rate = 0.001,
                         my.branch = 'v', branch_time = o.t$maximum)

  dt.sed.poly.tst[, .N, gamma_h_0]
  dt.sed.poly.tst[, .N, gamma_h_1]
  
  o.t = optimize(function(my.t) {
    x = q_theta(dt.sed.poly.tst,
                mh_contam = 0,
                faunal_prop = 0,
                err_rate = 0.001) +
      q_t(dt.sed.poly.tst, sims.dat = sims.dat,
          branch_time = my.t, my.branch = 'v')
    x},
    c(0.5585724, 0.8893428),
    maximum = T)
  print(t.last)
  print(o.t$maximum)
  print(o.t$maximum-t.last)
  
  
  x.range = seq(.747, .752, .0002)
  x = sapply(x.range,
             function(my.t) sed_EM(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 1,
                                   set.mh_contam = 0, set.faunal_prop = 0, set.branchtime = my.t)$max.ll)
  plot(x.range, x)
  abline(v=o.t$maximum, col='blue')
  abline(v=x.range[x == max(x)], col='red')
  
  
  
  x.manll = sapply(x.range,
             function(my.t) {
               all.rg = dt.sed.poly.tst[, unique(rg)]
               iter.ll = update_gamma(dt.sed.poly.tst, all.rg = all.rg,
                                      sims.dat = sims.dat,
                                      dt.theta = data.table(rg = all.rg,
                                                            mh_contam = 0,
                                                            faunal_prop = 0),
                                      err_rate = 0.001,
                                      my.branch = 'v', branch_time = my.t)
               iter.ll
               # only works if theta is the same for all rg
               q_theta(dt.sed.poly.tst,
                       mh_contam = 0,
                       faunal_prop = 0,
                       err_rate = 0.001) +
                 q_t(dt.sed.poly.tst, sims.dat = sims.dat,
                     branch_time = my.t, my.branch = 'v')
               # iter.ll
             })
  plot(x.range, x.manll)
  abline(v=o.t$maximum, col='blue')
  abline(v=x.range[x.manll == max(x.manll)], col='red')

  
  ## this gets the optimal q_t for the most recent gamma
  x.manqt = sapply(x.range,
                   function(my.t) {
                       q_t(dt.sed.poly.tst, sims.dat = sims.dat,
                           branch_time = my.t, my.branch = 'v')
                   })
  plot(x.range, x.manqt)
  abline(v=o.t$maximum, col='blue')
  abline(v=x.range[x.manqt == max(x.manqt)], col='red')

  
  x.manll2 <- sapply(x.range,
                     function(my.t) {
                       
                       err_rate = 0.001
                       g0 <- calc_gamma_num(dt.sed.poly.tst, gt=0,
                                            mh_contam = 0,
                                            faunal_prop = 0,
                                            err_rate = err_rate,
                                            sims.dat = sims.dat, 
                                            my.branch = 'v', 
                                            branch_time = my.t)
                       g1 <- calc_gamma_num(dt.sed.poly.tst, gt=1,
                                            mh_contam = 0,
                                            faunal_prop = 0,
                                            err_rate = err_rate,
                                            sims.dat = sims.dat, 
                                            my.branch = 'v', 
                                            branch_time = my.t)
                       sum(log(g0+g1))
                     })
  plot(x.range, x.manll2)
  abline(v=o.t$maximum, col='blue')
  abline(v=x.range[x.manll2 == max(x.manll2)], col='red')
  
  
  
  dt.sed.poly.tst <- data.table(v_gt=2, c_gt=2, a_gt=0, d_gt=0, f_mh=.5, 
                                sed_gt=c(rep(0,1), rep(1,3)), p_h_der = .75)
  dt.sed.poly.tst <- data.table(v_gt=2, c_gt=2, a_gt=0, d_gt=0, f_mh=.5, 
                                sed_gt=c(rep(0,2), rep(1,2)), p_h_der = .75)
  dt.sed.poly.tst[, rg := paste0('tst_deam_', F)]
  
  
  # x.range = seq(.747, .752, .0002)
  # x.range = seq(0.1, 0.8893428, .1)
  x.manll2.dt = foreach(err_rate = seq(0,1,.1), .combine = rbind) %:%
    foreach(my.t = seq(0.4, 1, 0.01), .combine = rbind) %do% {
      # foreach(my.t = seq(0.5585724, 0.8893428, 0.005), .combine = rbind) %do% {
        # err_rate = 0.001
      g0 <- calc_gamma_num(dt.sed.poly.tst, gt=0,
                           mh_contam = 0,
                           faunal_prop = 0,
                           err_rate = err_rate,
                           sims.dat = sims.dat, 
                           my.branch = 'v', 
                           branch_time = my.t)
      g1 <- calc_gamma_num(dt.sed.poly.tst, gt=1,
                           mh_contam = 0,
                           faunal_prop = 0,
                           err_rate = err_rate,
                           sims.dat = sims.dat, 
                           my.branch = 'v', 
                           branch_time = my.t)
      data.table(my.t, err_rate, ll=sum(log(g0+g1)))
    }
  x.manll2.dt[, is.max := ll == max(ll), err_rate]
  # ggplot(x.manll2.dt, aes(x=my.t, y=ll, color=err_rate, group=err_rate)) + geom_line()
  ggplot(x.manll2.dt, aes(x=my.t, y=ll, color=err_rate, group=err_rate)) + geom_line() +
    geom_point(data=x.manll2.dt[err_rate < 1 & is.max==T])
  ggplot(x.manll2.dt[ll > -3], aes(x=my.t, y=ll, color=err_rate, group=err_rate)) + geom_line() +
    geom_point(data=x.manll2.dt[err_rate < 1 & is.max==T]) +
    geom_vline(xintercept = .75, lty=3)
  

  e=0.001
  bt=.77
  ll0 = dt.sed.poly.tst[sed_gt == 0, .N] * log( e/2     * bt + (1-e/2) * (1-bt) )
  ll1 = dt.sed.poly.tst[sed_gt == 1, .N] * log( (1-e/2) * bt + e/2     * (1-bt) )
  ll0 + ll1

  x.manll3.dt = foreach(err_rate = seq(0,1,.1), .combine = rbind) %:%
    foreach(my.t = seq(0.4, 1, 0.01), .combine = rbind) %do% {
      e=err_rate
      bt=my.t
      ll0 = dt.sed.poly.tst[sed_gt == 0, .N] * log( e/2     * bt + (1-e/2) * (1-bt) )
      ll1 = dt.sed.poly.tst[sed_gt == 1, .N] * log( (1-e/2) * bt + e/2     * (1-bt) )
      data.table(my.t, err_rate, ll=(ll0+ll1))
    }
  x.manll3.dt[, is.max := ll == max(ll), err_rate]
  ggplot(x.manll3.dt[ll > -3], aes(x=my.t, y=ll, color=err_rate, group=err_rate)) + geom_line() +
    geom_point(data=x.manll3.dt[err_rate < 1 & is.max==T]) +
    geom_vline(xintercept = .75, lty=3)
  
  
  
}

points(
c(0.69598252416964945244,
0.69594182369453505554 ,
0.69602322464476384933 ,
0.69598252416964945244 ),

c(-3308.34103426361934907618,
 -3308.34103435032784545911,
 -3308.34103467753675431595,
 -3308.34103426361934907618), col='red'
)

# Goyet_A9229
# Les_Cottes_A9230
# Mez2_A9180
# SpyR5556
# Vindija_G1_A9348

###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################
###################################################




###################################################
###################################################


