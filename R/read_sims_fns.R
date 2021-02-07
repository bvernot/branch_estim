
####
## essentially reconstruct a sims.dat from dt.sims.p table
## eventually allow any combination of genotypes?


add_gt_categories <- function(dt.sims.p, agCols) {
  dt.sims.p[, gt.category := apply(.SD, 1, paste0, collapse=''), .SDcols=agCols]
}

get_fixed_gt_categories <- function(dt.sims.p, agCols, category = 'anc') {
  ## i.e. 2222
  if (category == 'der') return(paste0(dt.sims.p[, apply(.SD, 2, max), .SDcols=agCols], collapse=''))
  if (category == 'anc') return(paste0(dt.sims.p[, apply(.SD, 2, min), .SDcols=agCols], collapse=''))
  cat('must provide anc or der to get_fixed_gt_categories()\n')
  stop(1)
}

read_sims <- function(simsfile, agCols, keep.gt.cols = F) {
  dt.sims.p = fread(simsfile)
  if ('c' %in% colnames(dt.sims.p)) setnames(dt.sims.p, c('v', 'c', 'a', 'd'), c('v_gt', 'c_gt', 'a_gt', 'd_gt'))
  
  add_gt_categories(dt.sims.p, agCols)
  
  ### ISSUE : have to eventually allow this to be flexible, and to *recalculate* p(der) based on the columns selected
  ### ISSUE : also have to make sure that everything is polymorphic when taking the probabilities from the sims
  ###  i.e., if the user requests v/c/a, but d/mh were in the sims, then p(der) at 000 sites represents p(der) when d/mh are polymorphic
  # dt.sims.p[, gt.category := paste0(v_gt, c_gt, a_gt, d_gt)]
  
  setkey(dt.sims.p, gt.category)
  
  
  if (keep.gt.cols) {
    cat('simply adding p(der) to dt.sims.p, which can leave extra rows in the dt, which might cause problems later on [not currently supported]\n')
    dt.sims.p[, p := sum(n_der)/sum(n_tot), .(gt.category, branch, time)]
    dt.sims.p[, n_der := sum(n_der), .(gt.category, branch, time)]
    dt.sims.p[, n_tot := sum(n_tot), .(gt.category, branch, time)]
  } else {
    dt.sims.p <- dt.sims.p[, .SD, .SDcols = c(agCols, 'n_der', 'n_tot', 'branch', 'time', 'gt.category')]
    cat('debug\n')
    print(dt.sims.p[gt.category %like% '^221' & branch == 'v' & time == .9])
    dt.sims.p[, p := sum(n_der)/sum(n_tot), .(gt.category, branch, time)]
    dt.sims.p[, n_der := sum(n_der), .(gt.category, branch, time)]
    dt.sims.p[, n_tot := sum(n_tot), .(gt.category, branch, time)]
    dt.sims.p <- unique(dt.sims.p) ## does this always work with floating points?  it should?
  }
  
  
  
  all.branches = dt.sims.p[, .(time = max(time)), branch]
  setkey(all.branches, time)
  trunk.branch <- all.branches[.N, branch]
  x.branch = ''
  for(x in 1:all.branches[, .N-1]) {
    t1 = all.branches[x, time]
    t2 = all.branches[x+1, time]
    dt.sims.p.tmp <- dt.sims.p[time == t1 & (branch == trunk.branch | branch == x.branch)]
    x.branch = sprintf('anc_%d', x)
    ## print(dt.sims.p.tmp)
    dt.sims.p.tmp[, branch := x.branch]
    dt.sims.p[branch == trunk.branch & time > t1 & time <= t2, branch := x.branch]
    dt.sims.p <- rbind(dt.sims.p, dt.sims.p.tmp)
  }
  setkey(dt.sims.p, gt.category)
  branch.bounds <- dt.sims.p[, .(t.low = min(time), t.high = max(time)), branch]
  branch.bounds[, constr := 'high']
  branch.bounds[branch %like% 'anc', constr := 'both']
  branch.bounds[branch == x.branch, constr := 'low']
  inner.node.times <- branch.bounds[branch %like% 'anc', sort(t.low)]
  ret <- list()
  ret$dt.sims.p <- dt.sims.p
  ret$branch.bounds <- branch.bounds
  ret$inner.node.times <- inner.node.times
  ret$branches <- branch.bounds$branch
  ret$gt.fixed_der <- get_fixed_gt_categories(dt.sims.p, agCols, 'der') ## i.e. 2222 or 2221 if last gt is haploid
  ret$gt.fixed_anc <- get_fixed_gt_categories(dt.sims.p, agCols, 'anc') ## i.e. 0000
  ret
}
read_sims('data/sims_og_newrun.txt', agCols = c('v_gt', 'c_gt', 'a_gt', 'd_gt'))$dt.sims.p[gt.category %like% '^221' & branch == 'v' & time == .9]
read_sims('data/sims_og_newrun.txt', agCols = c('v_gt', 'c_gt', 'a_gt'))


generate_sims_dat <- function(simsfile, fixed_anc_p = 0.003, fixed_der_p = 0.999, agCols = c('v_gt', 'c_gt', 'a_gt', 'd_gt')) {
  
  sims.dat <- read_sims(simsfile, agCols)
  
  ## bounds for branch
  bounds.for.branch <- function(b, which.bound, sims.dat) {
    bc <- sims.dat$branch.bounds[branch == b]
    if (which.bound == 'high') {
      return(bc$t.high)
    }
    if (which.bound == 'low') {
      return(bc$t.low)
    }
    if (which.bound == 'mid') {
      return(mean(c(bc$t.high, bc$t.low)))
    }
    return(NA)
  }
  
  sims.dat$bounds.for.branch <- bounds.for.branch
  
  # ## have to make:
  # sims.dat$bounds.for.branch()
  # sims.dat$branch.bounds
  # sims.dat$branches
  sims.dat
  
  
  # sims.dat$dt.sims.p[gt.category == '1210', .(t.low = min(time), p.low = min(p),
  #                                             t.high = min(time), p.high = max(p)), .(gt.category, branch)]
  
  x = merge(sims.dat$dt.sims.p, sims.dat$branch.bounds, by='branch')
  sims.dat$dt.simple_p_given_b_t_arcs <- merge(x[t.low == time, .(t.low, p.low = p), keyby=.(gt.category, branch)],
                                               x[t.high == time, .(t.high, p.high = p), keyby=.(gt.category, branch)])
  
  setkey(sims.dat$dt.simple_p_given_b_t_arcs, gt.category)
  
  ## sites where all archaics are fixed aren't a good match for the simulations
  ## both because they can be 'QC' sites, and thus
  sims.dat$dt.simple_p_given_b_t_arcs[sims.dat$gt.fixed_anc, p.low := fixed_anc_p]
  sims.dat$dt.simple_p_given_b_t_arcs[sims.dat$gt.fixed_anc, p.high := fixed_anc_p]
  sims.dat$dt.simple_p_given_b_t_arcs[sims.dat$gt.fixed_der, p.low := fixed_der_p]
  sims.dat$dt.simple_p_given_b_t_arcs[sims.dat$gt.fixed_der, p.high := fixed_der_p]
  
  sims.dat$simple_p_given_b_t_arcs <- function(my.b, my.t, my.gt, sims.dat) {
    sims.dat$dt.simple_p_given_b_t_arcs[my.gt][my.b == branch][, p.low + (my.t-t.low)/(t.high-t.low) * (p.high-p.low)]
  }
  sims.dat
  
}
if (F) {
  generate_sims_dat('data/sims_og_newrun.txt')$bounds.for.branch('v', 'low', generate_sims_dat('sims_og_newrun.txt'))
  generate_sims_dat('data/sims_og_newrun.txt')$simple_p_given_b_t_arcs('v', 0.7801, '1000', generate_sims_dat('sims_og_newrun.txt'))
  generate_sims_dat('data/sims_og_newrun.txt')$simple_p_given_b_t_arcs('v', 0.7811, '1000', generate_sims_dat('sims_og_newrun.txt'))
  generate_sims_dat('data/sims_og_newrun.txt')$simple_p_given_b_t_arcs('v', 0.7821, '1000', generate_sims_dat('sims_og_newrun.txt'))
  generate_sims_dat('data/sims_og_newrun.txt')
  generate_sims_dat('data/sims_og_newrun.txt', agCols = c('v_gt', 'c_gt', 'a_gt'))$dt.simple_p_given_b_t_arcs['100']
}

