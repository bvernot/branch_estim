options(width=200)
library(argparse)


library(ggplot2)
library(data.table)
library(dplyr)
library(cobs)
library(foreach)
library(doParallel)
library(ggrepel)

setDTthreads(1)

print_tables = F





# #######
# ## make a new function that gives you a 'simple' p_der, with a linear change between the estimated endpoints on a branch
# 
# ## this really needs to be rewritten, because it should also re-compute the data from a simulation table?
# 
# add_linear_p_given_b_t_arcs <- function(sims.dat, fixed_anc_p = 0.01, fixed_der_p = .99) {
#   sims.dat$dt.simple_p_given_b_t_arcs <- data.table(expand.grid(v_gt=0:2,c_gt=0:2,a_gt=0:2,d_gt=0:2,branch=sims.dat$branches))
#   sims.dat$dt.simple_p_given_b_t_arcs <- merge(sims.dat$dt.simple_p_given_b_t_arcs, sims.dat$branch.bounds, by='branch')
#   sims.dat$dt.simple_p_given_b_t_arcs[, p.low := sims.dat$p_gt_given_b_t_arcs(branch, t.low, gt.category, sims.dat),
#                                       keyby=.(v_gt,c_gt,a_gt,d_gt,branch)]
#   sims.dat$dt.simple_p_given_b_t_arcs[, p.high := sims.dat$p_gt_given_b_t_arcs(branch, t.high, gt.category, sims.dat),
#                                       keyby=.(v_gt,c_gt,a_gt,d_gt,branch)]
#   setkey(sims.dat$dt.simple_p_given_b_t_arcs, v_gt,c_gt,a_gt,d_gt)
# 
#   ## sites where all archaics are fixed aren't a good match for the simulations
#   ## both because they can be 'QC' sites, and thus
#   sims.dat$dt.simple_p_given_b_t_arcs[.(0,0,0,0), p.low := fixed_anc_p]
#   sims.dat$dt.simple_p_given_b_t_arcs[.(0,0,0,0), p.high := fixed_anc_p]
#   sims.dat$dt.simple_p_given_b_t_arcs[.(2,2,2,2), p.low := fixed_der_p]
#   sims.dat$dt.simple_p_given_b_t_arcs[.(2,2,2,2), p.high := fixed_der_p]
# 
#   sims.dat$simple_p_given_b_t_arcs <- function(my.b, my.t, my.gt, sims.dat) {
#     sims.dat$dt.simple_p_given_b_t_arcs[my.gt][my.b == branch][, p.low + (my.t-t.low)/(t.high-t.low) * (p.high-p.low)]
#   }
#   sims.dat
# }


# add_linear_p_given_b_t_arcs_even_for_fixed <- function(sims.dat) {
#   sims.dat$dt.simple_p_given_b_t_arcs_even_for_fixed <- data.table(expand.grid(v_gt=0:2,c_gt=0:2,a_gt=0:2,d_gt=0:2,branch=sims.dat$branches))
#   sims.dat$dt.simple_p_given_b_t_arcs_even_for_fixed <- merge(sims.dat$dt.simple_p_given_b_t_arcs_even_for_fixed, sims.dat$branch.bounds, by='branch')
#   sims.dat$dt.simple_p_given_b_t_arcs_even_for_fixed[, p.low := sims.dat$p_gt_given_b_t_arcs(branch, t.low, gt.category, sims.dat),
#                                       keyby=.(v_gt,c_gt,a_gt,d_gt,branch)]
#   sims.dat$dt.simple_p_given_b_t_arcs_even_for_fixed[, p.high := sims.dat$p_gt_given_b_t_arcs(branch, t.high, gt.category, sims.dat),
#                                       keyby=.(v_gt,c_gt,a_gt,d_gt,branch)]
#   setkey(sims.dat$dt.simple_p_given_b_t_arcs_even_for_fixed, v_gt,c_gt,a_gt,d_gt)
#   sims.dat$dt.simple_p_given_b_t_arcs_even_for_fixed[.(0,0,0,0), p.low := 0.001]
#   sims.dat$dt.simple_p_given_b_t_arcs_even_for_fixed[.(0,0,0,0), p.high := 0.001]
#   sims.dat$dt.simple_p_given_b_t_arcs_even_for_fixed[.(1,1,1,1), p.low := 0.999]
#   sims.dat$dt.simple_p_given_b_t_arcs_even_for_fixed[.(1,1,1,1), p.high := 0.999]
# 
#   sims.dat$simple_p_given_b_t_arcs_even_p_for_fixed <- function(my.b, my.t, my.gt, sims.dat) {
#     sims.dat$dt.simple_p_given_b_t_arcs_even_for_fixed[my.gt][my.b == branch][, p.low + (my.t-t.low)/(t.high-t.low) * (p.high-p.low)]
#   }
#   sims.dat
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






##########
## 

p_der_given_h_theta <- function(dt.sed, mh_contam, faunal_prop, gt, err_rate) {
  ## calculate the probability of sampling a derived allele given an underlying (haploid) genotype gt,
  ## and contamination rates from fauna and MH, and a sequencing error rate
  ## e * (1-p') + (1-e) * p'
  ## p' = mh_contam * freq_mh + faunal_prop * 0 + (1-faunal_prop-mh_contam) * gt
  
  # p_prime = mh_contam * dt.sed[,f_mh] + (1-faunal_prop-mh_contam) * gt
  # p_prime = mh_contam * dt.sed[,f_mh] + faunal_prop * (pfd) + max(1-faunal_prop-mh_contam,0) * gt
  
  p_prime = mh_contam * dt.sed[,f_mh] + faunal_prop * 0 + (1-faunal_prop-mh_contam) * gt
  ## why should the errors be tied to the estimated der allele frequency in the data? [err_rate * (1-p_prime)]
  # err_rate * (1-p_prime) + (1-err_rate) * p_prime
  ## here it's: 1/2 the time an error makes a derived allele, regardless of p_prime. 1/2 the time it makes ancestral, and the rest of the time you sample from the underlying data
  1*err_rate/2 + 0*err_rate/2 + (1-err_rate) * p_prime
}
if (F) {
  p_der_given_h_theta(dt.sed.qc, mh_contam = 0.1, faunal_prop = 0.1, gt = 1, err_rate = 0.001)
  p_der_given_h_theta(dt.sed.poly, mh_contam = 0.1, faunal_prop = 0.1, gt = 1, err_rate = 0.001)
  p_der_given_h_theta(dt.sed.poly, mh_contam = 0.1, faunal_prop = .9, gt = 1, err_rate = 0.001)
  p_der_given_h_theta(dt.sed.poly, mh_contam = 0.1, faunal_prop = .1, gt = 1, err_rate = 0.00) %>% unique()
}


##########
## This is the probability of observing the reads from the sediment data, given a "true" underlying genotype, and contam rates [theta]
##########

p_O_given_h_theta <- function(dt.sed, mh_contam, faunal_prop, gt, err_rate) {
  p <- p_der_given_h_theta(dt.sed, mh_contam = mh_contam,
                           faunal_prop = faunal_prop,
                           gt=gt, err_rate=err_rate)
  # print(head(p))
  sed_gt <- dt.sed[, sed_gt]
  ## if the sediment state is ancestral, then we need the probability of observing the 
  ##  ancestral allele (not the derived allele, which is p)
  p[sed_gt == 0] <- 1-p[sed_gt == 0]
  p
}
if (F) {
  p_O_given_h_theta(dt.sed.poly, mh_contam = 0.1, faunal_prop = 0.1, gt = 1, err_rate = 0.001)
  p_O_given_h_theta(dt.sed.poly, mh_contam = 0.1, faunal_prop = .1, gt = 0, err_rate = 0.001)
  p_O_given_h_theta(dt.sed.poly, mh_contam = 0.0, faunal_prop = 1, gt = 1, err_rate = 0.001)
  
  dt.p_O <- data.table(expand.grid(mh_contam = seq(0,1,0.1), faunal_prop = seq(0,1,0.1)))
  dt.p_O <- dt.p_O[mh_contam + faunal_prop <= 1]
  dt.p_O[, lik := sum(log(p_O_given_h_theta(rbind(dt.sed.qc,dt.sed.poly), mh_contam, faunal_prop, gt = 1, err_rate = 0.001) +
                            p_O_given_h_theta(rbind(dt.sed.qc,dt.sed.poly), mh_contam, faunal_prop, gt = 0, err_rate = 0.001))),
         .(mh_contam, faunal_prop)]
  ## always gives pretty high estimate of MH contam - probably because this is the only source of any "data" in the calculation
  ## because h = 0 and h = 1 are equally weighted
  ggplot(dt.p_O, aes(x=mh_contam, y=faunal_prop, fill=lik)) + geom_tile() + 
    scale_fill_distiller(type = 'div')
}

##########
## This is the probability of having a particular genotype, given a branching point and time of branching [theta]
##########


#####
## called in q_t and calc_gamma, currently
p_h_given_b_t <- function(dt.sed, sims.dat, gt, branch, branch_time, method) {
  ## p_gt_given_b_t_arcs gives p(H==der | etc)
  
  dt.sed
  
  if (method == 'full') {
    # cat('full')
    dt.p <- dt.sed[, .(p = sims.dat$p_gt_given_b_t_arcs(branch, branch_time, gt.category, sims.dat)),
                   keyby=gt.category]
  } else if (method == 'simple') {
    # cat('simple')
    dt.p <- dt.sed[, .(p = sims.dat$simple_p_given_b_t_arcs(branch, branch_time, gt.category, sims.dat)),
                   keyby=gt.category]
  } else if (method == 'bt') {
    # cat('bt')
    dt.p <- dt.sed[, .(p = branch_time),
                   keyby=gt.category]
  } else if (method == 'bt2') {
    # cat('bt2')
    dt.p <- dt.sed[, .(p = branch_time/(c_gt+1)),
                   keyby=gt.category]
  } else if (method == 'bt3') {
    # cat('bt3')
    dt.p <- dt.sed[, .(p = ifelse(a_gt == 0, 0.001, branch_time / (c_gt+1))),
                   keyby=gt.category]
  }
  # dt.p <- dt.sed[, .(p = branch_time),
  #                keyby=gt.category]
  ## this [CURRENTLY] only happens for QC sites. 
  ## the other sites we are considering are polymorphic in archaics
  ## WOULD HAVE TO CHANGE THIS if I add sites that are fixed in archaics but poly in MH
  ## - this is actually a lot of sites, so it would add more information
  ##################
  ###### THIS SHOULDN'T BE DONE HERE, RIGHT?  SHOULD BE DONE IN THE FUNCTIONS?  OR MAYBE NOT, THIS CAN BE A FLAG
  ###### EITHER WAY, REMEMBER THIS, PUT A FLAG ON IT - IT TOTALLY OBVIATES CODE IN add_linear_p_given_b_t_arcs!!
  ##################
  ###### this also needs to be adjusted based on genotype category - I guess it should be set up above, in the main 'simple' function and not here
  ###### ISSUE
  ##################
  # dt.p[v_gt == 2 & c_gt == 2 & a_gt == 2 & d_gt == 2, p := .999]
  # dt.p[v_gt == 0 & c_gt == 0 & a_gt == 0 & d_gt == 0, p := .004]
  ## ISSUE - I JUST COMMENTED THESE OUT, AFTER RUNNING (AND SAVING) SOME DATA USING JUST NEAND POLYMORPHISMS
  # dt.p['2222', p := .999]
  # dt.p['0000', p := .004]
  ## at this point, p is p(H==der).  We want p(H==gt).  Flip if necessary.
  if (gt == 0) dt.p[, p := 1-p]
  ## unlike merge, plyr::join maintains the order of dt.sed
  # dt.sed.p = plyr::join(dt.sed[, .(v_gt,c_gt,a_gt,d_gt)], dt.p, by=c('v_gt','c_gt','a_gt','d_gt'))[, p]
  # dt.sed.p = plyr::join(dt.sed[, gt.category], dt.p, by='gt.category')
  dt.sed.p = dt.p[dt.sed$gt.category]
  if (sum((dt.sed[, gt.category] == dt.sed.p[, gt.category]) == F) > 0) {
    cat('join in p_h_given_b_t gives wrong order??\n')
    print(knitr::kable(dt.sed))
    print(knitr::kable(dt.sed.p))
    return('hey')
  }
  dt.sed.p[, p]
}
if (F) {
  p_h_given_b_t(dt.sed.poly, sims.dat, gt = 1, branch = 'v', branch_time = .7) %>% head
  p_h_given_b_t(dt.sed.poly, sims.dat, gt = 0, branch = 'v', branch_time = .7) %>% head
  p_h_given_b_t(dt.sed.qc, sims.dat, gt = 0, branch = 'v', branch_time = .7)
  
  
  dt.p_0.sed <- data.table(dt.sed.poly)
  dt.p_0.sed <- data.table(dt.sed.qc)
  dt.p_0.sed <- rbind(dt.sed.qc,dt.sed.poly)
  dt.p_O <- data.table(expand.grid(mh_contam = seq(0,1,0.05), faunal_prop = seq(0,1,0.05)))
  dt.p_O <- dt.p_O[mh_contam + faunal_prop <= 1]
  p_h0 = p_h_given_b_t(dt.p_0.sed, sims.dat, gt = 0, branch = 'anc_1', branch_time = 1.1)
  p_h1 = p_h_given_b_t(dt.p_0.sed, sims.dat, gt = 1, branch = 'anc_1', branch_time = 1.1)
  dt.p_O[, lik := sum(log(p_O_given_h_theta(dt.p_0.sed, mh_contam, faunal_prop, 
                                            gt = 1, err_rate = 0.001) * p_h1 +
                            p_O_given_h_theta(dt.p_0.sed, mh_contam, faunal_prop,
                                              gt = 0, err_rate = 0.001) * p_h0)),
         .(mh_contam, faunal_prop)]
  ggplot(dt.p_O, aes(x=mh_contam, y=faunal_prop, fill=lik)) + geom_tile() + 
    scale_fill_distiller(type = 'div')
  
  
  my.branch = 'c'
  dt.p_0.sed <- rbind(dt.sed.qc,dt.sed.poly)
  # dt.p_0.sed <- rbind(dt.sed.qc)
  dt.p_O <- data.table(expand.grid(mh_contam = seq(0,1,0.1), 
                                   faunal_prop = seq(0,1,0.1),
                                   branch_time = seq(sims.dat$bounds.for.branch(my.branch, 'low', sims.dat),
                                                     sims.dat$bounds.for.branch(my.branch, 'high', sims.dat),
                                                     .05)))
  dt.p_O <- dt.p_O[mh_contam + faunal_prop <= 1]
  
  dt.p_O.res <- foreach(bt = unique(dt.p_O$branch_time), .combine=rbind) %do% {
    cat(bt, ' ')
    p_h0 = p_h_given_b_t(dt.p_0.sed, sims.dat, gt = 0, branch = my.branch, branch_time = bt)
    p_h1 = p_h_given_b_t(dt.p_0.sed, sims.dat, gt = 1, branch = my.branch, branch_time = bt)
    ## not correct?
    dt.p_O[branch_time == bt, .(lik = sum(log(p_O_given_h_theta(dt.p_0.sed, mh_contam, faunal_prop, 
                                                                gt = 1, err_rate = 0.001) * p_h1 +
                                                p_O_given_h_theta(dt.p_0.sed, mh_contam, faunal_prop,
                                                                  gt = 0, err_rate = 0.001) * p_h0))),
           .(mh_contam, faunal_prop, branch_time)]
  }
  ## always gives pretty high estimate of MH contam - probably because this is the only source of any "data" in the calculation
  ## because h = 0 and h = 1 are equally weighted
  ggplot(dt.p_O.res, aes(x=mh_contam, y=faunal_prop, color=lik)) + geom_jitter() + 
    scale_color_distiller(type = 'div')
  ggplot(dt.p_O.res[, .(lik=max(lik)), branch_time], aes(x=branch_time, y=lik)) + geom_line()
}


###################################################
###################################################
###################################################
###################################################

# for each branch separately

if (F) {
  # data:
  dt.sed.analysis = rbind(dt.sed.qc,dt.sed.poly)
  dt.sed.analysis = rbind(dt.sed.qc,dt.sed.poly[!(v_gt == c_gt & v_gt == a_gt)])
  dt.sed.analysis[, rg := paste0('deam_', deam53)]
  all.rg = dt.sed.analysis[, unique(rg)]
  # initialized parameters. theta is different for each rg, but there's only one branch_time
  dt.theta = data.table(rg = all.rg,
                        mh_contam = .5,
                        faunal_prop = .5)
  branch_time = 1
  err_rate = 0.001
  
  # # initialized genotype likelihoods - random? .5?
  # dt.sed.analysis[, gamma_h_0 := .5]
  # dt.sed.analysis[, gamma_h_1 := .5]
}


########################
## for each read group, split dt.sed.analysis, and optimize theta using q_theta
########################
q_theta <- function(dt.sed, mh_contam, faunal_prop, err_rate) {
  # cat('q_theta', mh_contam, faunal_prop, '\n')
  q_h0 <- p_O_given_h_theta(dt.sed, mh_contam, faunal_prop, 0, err_rate)
  q_h1 <- p_O_given_h_theta(dt.sed, mh_contam, faunal_prop, 1, err_rate)
  if (print_tables) print(knitr::kable(dt.sed[, .(q_theta_h0 = log(q_h0) * gamma_h_0, q_theta_h1 = log(q_h1) * gamma_h_1)]))
  dt.sed[, sum(log(q_h0) * gamma_h_0 + log(q_h1) * gamma_h_1)]
}
if (F) {
  q_theta(dt.sed.analysis, mh_contam = 0.05, faunal_prop = 0.05, err_rate = 0.001)
  q_theta(dt.sed.analysis, mh_contam = 0.05, faunal_prop = 0.05, err_rate = 0.00) # -Inf, can't give 0 error rate
  my.rg = 'deam_FALSE'
  my.rg = 'deam_TRUE'
  q_params = c(mh_contam = dt.theta[rg == my.rg, mh_contam],
               faunal_prop = dt.theta[rg == my.rg, faunal_prop],
               err_rate = err_rate) ## doesn't work all the time?
  q_theta(dt.sed.analysis[rg == my.rg], 
          mh_contam = q_params[1], 
          faunal_prop = q_params[2], 
          err_rate = q_params[3])
  optim(q_params, function(params) {
    cat(params, '\n')
    -q_theta(dt.sed.analysis[rg == my.rg], 
             mh_contam = params[1], 
             faunal_prop = params[2], 
             err_rate = params[3])},
    control = list(), lower = 0.0000001, upper = 0.9999999, method = 'L-BFGS-B')
  q_params = c(mh_contam = .1,
               faunal_prop = .1)
  constrOptim(q_params, function(params) {
    cat(params, '\n')
    -q_theta(dt.sed.analysis[rg == my.rg],
             mh_contam = params[1],
             faunal_prop = params[2],
             err_rate = err_rate)}, grad = NULL,
    ui=rbind(c(-1,-1),  # the -x-y > -1
             c(1,0),    # the x > 0
             c(0,1) ),  # the y > 0
    ci=c(-1,0, 0))
}

########################
## and optimize t using q_t and all of dt.sed.analysis
########################
q_t <- function(dt.sed, sims.dat, branch_time, my.branch, method) {
  q_h0 = p_h_given_b_t(dt.sed, sims.dat = sims.dat, gt = 0,
                       branch = my.branch, branch_time = branch_time, method = method)
  q_h1 = p_h_given_b_t(dt.sed, sims.dat = sims.dat, gt = 1,
                       branch = my.branch, branch_time = branch_time, method = method)
  qt = dt.sed[, sum(log(q_h0) * gamma_h_0 + log(q_h1) * gamma_h_1)]
  cat(sprintf('q_t %.20g %0.20f\n', branch_time, qt))
  qt
}
if (F) {
  q_t(dt.sed.analysis, sims.dat = sims.dat, branch_time = 1, my.branch = 'anc_1')
}

########################
## and recalculate genotype likelihoods - have to do this separately 
########################
calc_gamma_num <- function(dt.sed, gt, mh_contam, faunal_prop, err_rate, sims.dat, my.branch, branch_time, method) {
  if (gt == 0) {
    p_O_h0 <- p_O_given_h_theta(dt.sed, mh_contam = mh_contam, faunal_prop = faunal_prop, gt = 0, err_rate = err_rate)
    p_H_h0 <- p_h_given_b_t(dt.sed, sims.dat = sims.dat, gt = 0, branch = my.branch, branch_time = branch_time, method = method)
    if (print_tables) {
      cat('\ncalc_gamma_num0:')
      print(knitr::kable(data.table(p_O_h0 = p_O_h0, p_H_h0 = p_H_h0, p_O_h0_x_p_H_h0 = p_O_h0 * p_H_h0)))
    }
    return(p_O_h0 * p_H_h0)
  }
  if (gt == 1) {
    p_O_h1 <- p_O_given_h_theta(dt.sed, mh_contam = mh_contam, faunal_prop = faunal_prop, gt = 1, err_rate = err_rate)
    p_H_h1 <- p_h_given_b_t(dt.sed, sims.dat = sims.dat, gt = 1, branch = my.branch, branch_time = branch_time, method = method)
    if (print_tables) {
      cat('\ncalc_gamma_num1:')
      print(knitr::kable(data.table(p_O_h1 = p_O_h1, p_H_h1 = p_H_h1, p_O_h1_x_p_H_h1 = p_O_h1 * p_H_h1)))
    }
    return(p_O_h1 * p_H_h1)
  }
}
calc_manual_lik <- function(dt.sed, all.rg, sims.dat, dt.theta, err_rate, my.branch, branch_time, method) {
  iter.ll = 0
  
  for (my.rg in all.rg) {
    mh_contam = dt.theta[rg == my.rg, mh_contam]
    faunal_prop = dt.theta[rg == my.rg, faunal_prop]
    dt.sed.rg <- dt.sed[rg == my.rg]
    
    #### ISSUE
    #### should change calc_gamma_num to return a dt w/ g0 and g1 in it, to reduce time of calculation
    #### would cut running time [for calc_manual_lik] in half
    g0 <- calc_gamma_num(dt.sed.rg,
                         gt=0,
                         mh_contam = mh_contam,
                         faunal_prop = faunal_prop,
                         err_rate = err_rate,
                         sims.dat = sims.dat, 
                         my.branch = my.branch, 
                         branch_time = branch_time,
                         method = method)
    g1 <- calc_gamma_num(dt.sed.rg, ## wtf, this used to be just dt.sed
                         gt=1,
                         mh_contam = mh_contam,
                         faunal_prop = faunal_prop,
                         err_rate = err_rate,
                         sims.dat = sims.dat, 
                         my.branch = my.branch, 
                         branch_time = branch_time,
                         method = method)
    iter.ll = iter.ll + sum(log(g0+g1))
  }
  iter.ll
}

calc_gamma <- function(dt.sed, gt, mh_contam, faunal_prop, err_rate, sims.dat, my.branch, branch_time, method = 'full') {
  p_O_h0 <- p_O_given_h_theta(dt.sed, mh_contam = mh_contam, faunal_prop = faunal_prop, gt = 0, err_rate = err_rate)
  p_O_h1 <- p_O_given_h_theta(dt.sed, mh_contam = mh_contam, faunal_prop = faunal_prop, gt = 1, err_rate = err_rate)
  p_H_h0 <- p_h_given_b_t(dt.sed, sims.dat = sims.dat, gt = 0, branch = my.branch, branch_time = branch_time, method = method)
  p_H_h1 <- p_h_given_b_t(dt.sed, sims.dat = sims.dat, gt = 1, branch = my.branch, branch_time = branch_time, method = method)
  if (gt == 0) return(p_O_h0 * p_H_h0 / (p_O_h0 * p_H_h0 + p_O_h1 * p_H_h1))
  if (gt == 1) return(p_O_h1 * p_H_h1 / (p_O_h0 * p_H_h0 + p_O_h1 * p_H_h1))
}

## I don't think this works correctly - not sure why!
update_gamma <- function(dt.sed.analysis, all.rg, sims.dat,
                         dt.theta, err_rate, my.branch, branch_time, method = 'full') {
  iter.ll = 0
  
  for (my.rg in all.rg) {
    q_params = c(mh_contam = dt.theta[rg == my.rg, mh_contam],
                 faunal_prop = dt.theta[rg == my.rg, faunal_prop])
    dt.sed.analysis.rg <- dt.sed.analysis[rg == my.rg]
    g0 <- calc_gamma(dt.sed.analysis.rg, 0,
                     mh_contam = q_params['mh_contam'],
                     faunal_prop = q_params['faunal_prop'],
                     err_rate = err_rate, sims.dat = sims.dat, 
                     my.branch = my.branch, 
                     branch_time = branch_time, 
                     method = method)
    g1 <- calc_gamma(dt.sed.analysis.rg, 1, 
                     mh_contam = q_params['mh_contam'],
                     faunal_prop = q_params['faunal_prop'],
                     err_rate = err_rate, sims.dat = sims.dat, 
                     my.branch = my.branch, 
                     branch_time = branch_time, 
                     method = method)
    dt.sed.analysis[rg == my.rg, gamma_h_0 := g0]
    dt.sed.analysis[rg == my.rg, gamma_h_1 := g1]
    
    ## now get lik for this round
    iter.ll = iter.ll +
      q_theta(dt.sed.analysis[rg == my.rg],
              mh_contam = q_params['mh_contam'],
              faunal_prop = q_params['faunal_prop'],
              err_rate = err_rate) +
      q_t(dt.sed.analysis[rg == my.rg], sims.dat = sims.dat,
          branch_time = branch_time, my.branch = my.branch, 
          method = method)
  }
  iter.ll
}
if (F) {
  all.rg = dt.sed.analysis[, unique(rg)]
  iter.ll = update_gamma(dt.sed.analysis, all.rg = all.rg,
                         sims.dat = sims.dat,
                         dt.theta = data.table(rg = all.rg,
                                               mh_contam = 0,
                                               faunal_prop = 0),
                         err_rate = 0.001,
                         my.branch = 'v', branch_time = .75)

  o.t = optimize(function(my.t) q_t(dt.sed.analysis, sims.dat = sims.dat,
                                    branch_time = my.t, my.branch = 'v'),
                 c(0.5585724, 0.8893428),
                 maximum = T)
}


## if blocks > 0, bootstrap over blocks rather than single sites
## if drop.one > 0, then drop just that block
bootstrap_gt_data <- function(dt.sed.analysis, blocks = 0, drop.one = 0, shuffle_first = T) {
  dt <- data.table(dt.sed.analysis[sample(.N)])
  if (blocks == 0) return(dt[sample(.N,.N,replace=T)])

  block.labels <- sort(rep(1:blocks, length.out=dt[, .N]))
  dt[, block.label := as.character(block.labels)]
  setkey(dt, block.label)
  if (drop.one == 0) {
    these.blocks = as.character(sample(blocks, blocks, replace = T))
    dt[these.blocks]
  } else {
    dt[!as.character(drop.one)]
  }
}
# bootstrap_gt_data(dt.sed.analysis, 100)[, .N, block.label]
# bootstrap_gt_data(dt.sed.analysis, 10, 3)[, .N, block.label]




get_expanded_range <- function(vals.x, range.x) {
  ## ensure that there are only two values here
  range.x = range(range.x)
  if (!range.x[1] %in% vals.x) {
    cat('error in get_expanded_range [1]:', range.x[1], 'not in', vals.x, '\n')
    return(NULL)
  }
  if (!range.x[2] %in% vals.x) {
    cat('error in get_expanded_range [2]:', range.x[2], 'not in', vals.x, '\n')
    return(NULL)
  }
  idx.x = which(vals.x %in% range.x) + c(-1,1)
  idx.x[1] = max(idx.x[1], 1)
  idx.x[2] = min(idx.x[2], length(vals.x))
  vals.x[idx.x]
}
if (F) {
  get_expanded_range(2:10, 3:4)
  get_expanded_range(2:10, 2:4)
  get_expanded_range(2:10, 4:11)
  get_expanded_range(2:2, 2)
}


sed_grid_search <- function(dt.sed.analysis, sims.dat, my.branch, err_rate = 0.001, p_h_method = 'simple',
                            range.mh_contam = c(0,.3), range.faunal_prop = c(0,.2), range.t = NULL,
                            bins.mh_contam = 10, bins.faunal_prop = 5, bins.t = 10,
                            nsteps = 3, print.debug = F) {
  
  ######
  ## ISSUE - match.call does lazy-evaluation at its worst, so if e.g. dt.sed.analysis is:
  ##   dt.sed.analysis[sample(.N, .N, replace=T)]
  ## then it would randomly sample it for 
  ######
  ## NEW ISSUE : any time sed_grid_search_single_branch changes, this also has to change
  
  dt.ret <- foreach (my.b = my.branch, .combine = rbind) %do% {
    # mc = match.call()
    # mc$dt.sed.analysis <- dt.sed.analysis
    # mc$my.branch <- my.b
    # mc[[1]] <- as.name('sed_grid_search_single_branch')
    # eval(mc)
    sed_grid_search_single_branch(dt.sed.analysis = dt.sed.analysis, 
                                  sims.dat = sims.dat, 
                                  my.branch = my.b, ## only difference
                                  err_rate = err_rate,
                                  p_h_method = p_h_method,
                                  range.mh_contam = range.mh_contam,
                                  range.faunal_prop = range.faunal_prop, 
                                  range.t = range.t,
                                  bins.mh_contam = bins.mh_contam,
                                  bins.faunal_prop = bins.faunal_prop,
                                  bins.t = bins.t,
                                  nsteps = nsteps,
                                  print.debug = print.debug)
  }
  dt.ret[, is.max := ll == max(ll)]
  dt.ret[, is.top1pct := ll > quantile(ll,.99)]
  dt.ret
}
if (F) {
  gt.tmp = dt.sed.analysis[lib %like% 'Mez1'][sample(.N, .N, replace=T)]
  x.grid.mez1 = sed_grid_search(gt.tmp, sims.dat, my.branch = c('v','c', 'anc_1'),
                           p_h_method = 'simple', range.mh_contam = .0, range.faunal_prop = 0)
  
  x.grid.mez1 = sed_grid_search(dt.sed.analysis[lib %like% 'Mez1'][sample(.N, .N, replace=T)],
                                sims.dat, my.branch = c('v','c', 'anc_1'),
                                p_h_method = 'simple', range.mh_contam = .02, range.faunal_prop = 0, nsteps = 2)
  x.grid.mez1
  ggplot(x.grid.mez1, aes(y=ll, x=my.t, color=my.branch)) + geom_line() + geom_point(data=x.grid.mez1[is.max == T])
  
}

sed_grid_search_single_branch <- function(dt.sed.analysis, sims.dat, my.branch, err_rate = 0.001, p_h_method = 'simple',
                            range.mh_contam = c(0,.3), range.faunal_prop = c(0,.2), range.t = NULL,
                            bins.mh_contam = 10, bins.faunal_prop = 5, bins.t = 10,
                            nsteps = 3, print.debug = F) {
  cat('Grid search on', my.branch)
  ## if the user sets only a single value for the range of values to search, then set it up so the code still works [i.e. so it's still a range]
  if (length(range.mh_contam) == 1) {
    range.mh_contam <- rep(range.mh_contam, 2)
    bins.mh_contam <- 1
  }
  if (length(range.faunal_prop) == 1) {
    range.faunal_prop <- rep(range.faunal_prop, 2)
    bins.faunal_prop <- 1
  }
  if (length(range.t) == 1) {
    range.t <- rep(range.t, 2)
    bins.t <- 1
  }
  
  if (is.null(range.t)) {
    range.t <- sims.dat$branch.bounds[my.branch, c(t.low, t.high), on='branch']
    # cat(range.t, '\n')
    # cat(bins.t, '\n')
  }
  
  
  dt.ret = data.table()
  for (step.x in 1:nsteps) {
    ## do a grid search
    vals.mh_contam   = unique(seq(range.mh_contam[1], range.mh_contam[2], length.out = bins.mh_contam))
    vals.faunal_prop = unique(seq(range.faunal_prop[1], range.faunal_prop[2], length.out = bins.faunal_prop))
    vals.t           = unique(seq(range.t[1], range.t[2], length.out = bins.t))
    x.manll.dt = foreach(mh_contam = vals.mh_contam, .combine = rbind) %:%
      foreach(faunal_prop = vals.faunal_prop, .combine = rbind) %:%
      foreach(my.t = vals.t, .combine = rbind) %do% {
        if (print.debug) cat(mh_contam, faunal_prop, my.t, '\n')
        all.rg = dt.sed.analysis[, unique(rg)]
        ll = calc_manual_lik(dt.sed.analysis,
                             all.rg = all.rg,
                             dt.theta = data.table(rg = all.rg, mh_contam = mh_contam, faunal_prop = faunal_prop),
                             err_rate = err_rate,
                             sims.dat = sims.dat, 
                             my.branch = my.branch, 
                             branch_time = my.t,
                             method = p_h_method)
        data.table(my.t, mh_contam, faunal_prop, ll)
      }
    ## get range of top X ll per parameter
    range.faunal_prop = get_expanded_range(vals.faunal_prop, x.manll.dt[ll > quantile(ll,.99), range(faunal_prop)])
    range.mh_contam = get_expanded_range(vals.mh_contam, x.manll.dt[ll > quantile(ll,.99), range(mh_contam)])
    range.t = get_expanded_range(vals.t, x.manll.dt[ll > quantile(ll,.99), range(my.t)])
    cat('\n')
    cat('faunal_prop range:', step.x, range.faunal_prop, '\n')
    cat('mh_contam   range:', step.x, range.mh_contam, '\n')
    cat('branch time range:', step.x, range.t, '\n')
    dt.ret = rbind(dt.ret, data.table(x.manll.dt, step.x))
  }
  dt.ret[, is.max := ll == max(ll)]
  dt.ret[, is.top1pct := ll > quantile(ll,.99)]
  dt.ret[, my.branch := my.branch]
  dt.ret
}
if (F) {
  dt.sed.poly.tst.gridll = sed_grid_search(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, p_h_method = 'bt2',
                                           bins.mh_contam = 10, bins.faunal_prop = 6, bins.t = 10,
                                           range.mh_contam = c(0,.2), range.faunal_prop = c(0,.2), range.t = c(.55,.85))
  dt.sed.poly.tst.gridll = sed_grid_search(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, p_h_method = 'bt2',
                                           bins.mh_contam = 10, bins.faunal_prop = 1, bins.t = 10,
                                           range.mh_contam = c(0,.2), range.faunal_prop = c(0,0), range.t = c(.55,.85))
  ggplot(dt.sed.poly.tst.gridll[, .(ll = max(ll)), by=.(my.t, mh_contam)], aes(x=my.t, y=mh_contam, color = ll)) + geom_point()
}



################
# procedure:
################

sed_EM <- function(dt.sed.analysis, sims.dat, my.branch, err_rate = 0.001, max.iter = 20, ll.converge = 0, 
                   set.branchtime = 'estim', set.mh_contam = 'estim', set.faunal_prop = 'estim', p_h_method = 'full',
                   fail_on_neg_change = F, fail_on_neg_change_q = F) {


  # hey  
  ## save data to return
  ret = list()
  ret$ll.trace = c()
  ret$man.ll.trace = c()
  last.ll = -1e200
  last.iter.ll = -1e200

  # because you can't be sure that this isn't running in parallel..
  dt.sed.analysis <- data.table(dt.sed.analysis)
  
  # all readgroups:
  all.rg = dt.sed.analysis[, unique(rg)]
  
  branch.low <- sims.dat$bounds.for.branch(my.branch, 'low', sims.dat = sims.dat)
  branch.high <- sims.dat$bounds.for.branch(my.branch, 'high', sims.dat = sims.dat)
  
  # if ((set.faunal_prop != 'estim' || set.mh_contam != 'estim') && (set.faunal_prop == 'estim' || set.mh_contam == 'estim')) {
  #   cat('set.faunal_prop and set.mh_contam must both be provided, or neither.:', set.faunal_prop, set.mh_contam, '\n')
  #   stop(paste('set.faunal_prop and set.mh_contam must both be provided, or neither:', set.faunal_prop, set.mh_contam, '\n'))
  # }
  
  # initialized parameters. theta is different for each rg, but there's only one branch_time
  dt.theta = data.table(rg = all.rg,
                        mh_contam = ifelse(is.numeric(set.mh_contam), set.mh_contam, .1),
                        faunal_prop = ifelse(is.numeric(set.faunal_prop), set.faunal_prop, .1))
  if (set.branchtime == 'estim' | set.branchtime == 'grid') {
    branch_time = sims.dat$bounds.for.branch(my.branch, which.bound = 'mid', sims.dat = sims.dat)
  } else if (branch.low <= set.branchtime &
             branch.high >= set.branchtime) {
    branch_time = set.branchtime
  } else {
    cat('set.branchtime in sed_EM must be either "estim" or between low and high bounds for branch:', set.branchtime, '\n')
    stop(paste('set.branchtime in sed_EM must be either "estim" or between low and high bounds for branch:', set.branchtime, '\n'))
  }
  
  # initialized genotype likelihoods - random? .5?
  dt.sed.analysis[, gamma_h_0 := ifelse(sed_gt == 0, .9, .1)]
  dt.sed.analysis[, gamma_h_1 := ifelse(sed_gt == 1, .9, .1)]
  
  for (iter in 1:max.iter) {
    
    # iter.ll = 0
    
    ## for each read group, split dt.sed.analysis, and optimize theta using q_theta
    if (set.faunal_prop == 'estim' && set.mh_contam == 'estim') {
      for (my.rg in all.rg) {
        dt.sed.analysis.rg <- dt.sed.analysis[rg == my.rg]
        q_params = c(mh_contam = dt.theta[rg == my.rg, mh_contam],
                     faunal_prop = dt.theta[rg == my.rg, faunal_prop])
        # o.theta = optim(q_params, function(params) -q_theta(dt.sed.analysis[rg == my.rg],
        #                                                     mh_contam = params['mh_contam'],
        #                                                     faunal_prop = params['faunal_prop'],
        #                                                     err_rate = err_rate),
        #                 control = list(), lower = 0.000001, upper = 0.999999, method = 'L-BFGS-B')
        o.theta = constrOptim(q_params, function(params) {
          # cat(params, '\n')
          -q_theta(dt.sed.analysis.rg,
                   mh_contam = params[1],
                   faunal_prop = params[2],
                   err_rate = err_rate)}, grad = NULL,
          ui=rbind(c(-1,-1),  # the -x-y > -1
                   c(1,0),    # the x > 0
                   c(0,1) ),  # the y > 0
          ci=c(-1,0, 0))
        dt.theta[rg == my.rg, mh_contam := o.theta$par['mh_contam']]
        dt.theta[rg == my.rg, faunal_prop := o.theta$par['faunal_prop']]
        # print(o.theta)
        cat('theta convergence:', o.theta$convergence, '\n')
        # iter.ll = iter.ll - o.theta$value
      }
    } else if (set.mh_contam == 'estim') {
      for (my.rg in all.rg) {
        dt.sed.analysis.rg <- dt.sed.analysis[rg == my.rg]
        q_params = c(mh_contam = dt.theta[rg == my.rg, mh_contam],
                     faunal_prop = dt.theta[rg == my.rg, faunal_prop])
        o.theta = optimize(function(my.mc) {
          q_theta(dt.sed.analysis.rg,
                  mh_contam = my.mc,
                  faunal_prop = q_params[2],
                  err_rate = err_rate)},
          c(1e-200, 1-q_params[2]-1e-200), maximum = T)
        dt.theta[rg == my.rg, mh_contam := o.theta$maximum]
        dt.theta[rg == my.rg, faunal_prop := q_params[2]]
      }
    } else if (set.mh_contam == 'estim2') {
      for (my.rg in all.rg) {
        dt.sed.analysis.rg <- dt.sed.analysis[rg == my.rg]
        q_params = c(mh_contam = dt.theta[rg == my.rg, mh_contam],
                     faunal_prop = dt.theta[rg == my.rg, faunal_prop])
        o.theta = optimize(function(my.mc) {
          q_theta(dt.sed.analysis.rg,
                  mh_contam = my.mc,
                  faunal_prop = q_params[2],
                  err_rate = err_rate)},
          c(1e-200, 1-q_params[2]-1e-200), maximum = T)
        dt.theta[rg == my.rg, mh_contam := o.theta$maximum]
        dt.theta[rg == my.rg, faunal_prop := q_params[2]]
      }
    } else if (set.faunal_prop == 'estim') {
      for (my.rg in all.rg) {
        dt.sed.analysis.rg <- dt.sed.analysis[rg == my.rg]
        q_params = c(mh_contam = dt.theta[rg == my.rg, mh_contam],
                     faunal_prop = dt.theta[rg == my.rg, faunal_prop])
        o.theta = optimize(function(my.fp) {
          q_theta(dt.sed.analysis.rg,
                  mh_contam = q_params[1],
                  faunal_prop = my.fp,
                  err_rate = err_rate)},
          c(1e-200, 1-q_params[1]-1e-200), maximum = T)
        dt.theta[rg == my.rg, mh_contam := q_params[1]]
        dt.theta[rg == my.rg, faunal_prop := o.theta$maximum]
      }
    }
    # cat(iter, dt.sed.analysis[, sum(log(gamma_h_1) + log(gamma_h_1))], my.branch, o.t$maximum)
    
    if (set.branchtime == 'estim') {
      ## and optimize t using q_t and all of dt.sed.analysis
      o.t = optimize(function(my.t) q_t(dt.sed.analysis, sims.dat = sims.dat,
                                        branch_time = my.t, my.branch = my.branch, method = p_h_method),
                     c(branch.low, branch.high),
                     maximum = T)
    } else if (set.branchtime == 'grid') {
      ## and optimize t using q_t and all of dt.sed.analysis
      o.t = NULL
      for (my.t in seq(branch.low, branch.high, .01)) {
        q_t_obj = q_t(dt.sed.analysis, sims.dat = sims.dat,
                      branch_time = my.t, my.branch = my.branch, method = p_h_method)
        if (my.t == branch.low || o.t$objective < q_t_obj) {
          o.t = list(objective = q_t_obj,
                     maximum = my.t)
        }
      }
      ## could use sed_grid_search, but that does a manual likelihood calculation, and does not use gamma - shouldn't I be using gamma and optimizing q_t?  does it make a difference?
      # dt.sed.poly.tst.gridll.all_params = sed_grid_search(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, p_h_method = p_h_method,
      #                                                     bins.mh_contam = 10, bins.faunal_prop = 1, bins.t = 10, nsteps = 4,
      #                                                     range.mh_contam = c(0,.2), range.faunal_prop = c(0,0),
      #                                                     range.t = sims.dat$branch.bounds[branch == 'v', c(t.low, t.high)])
      
    } else {
      ## if we're not supposed to estimate branch time, just get the likelihood of the given branch time
      ## TRUE? we don't actually have to calculate q_t, we don't use it, and the likelihood shouldn't change 
      ## [it does change, because gamma changes - keeping it in for now] <- obv removed the q_t call
      o.t = list(objective = 0,
                 maximum = branch_time)
    }
    
    # iter.ll = iter.ll + o.t$objective
    
    ## now calculate new p(h) [gamma], using learned theta and t, and save
    iter.ll = update_gamma(dt.sed.analysis, all.rg = all.rg, sims.dat = sims.dat,
                           dt.theta = dt.theta, err_rate = err_rate,
                           my.branch = my.branch, branch_time = o.t$maximum, method = p_h_method)
    
    manual.ll = calc_manual_lik(dt.sed.analysis,
                                all.rg = all.rg,
                                sims.dat = sims.dat, 
                                dt.theta = dt.theta,
                                err_rate = err_rate,
                                my.branch = my.branch, 
                                branch_time = o.t$maximum,
                                method = p_h_method)
    

    ## report results of this iteration
    ret$ll.trace <- c(ret$ll.trace, iter.ll)
    ret$man.ll.trace <- c(ret$man.ll.trace, manual.ll)
    ret$dt.theta.trace <- rbind(ret$dt.theta.trace, 
                                data.table(dt.theta, iter = length(ret$man.ll.trace)))
    ret$branchtime.trace <- c(ret$branchtime.trace, o.t$maximum)

    
    if (nrow(dt.sed.analysis) < 10) {
      print(knitr::kable(dt.sed.analysis))
    }
    
    print(knitr::kable(dt.theta))
    print(sprintf('%g', dt.theta$mh_contam))
    print(sprintf('%g', dt.theta$faunal_prop))
    print(sprintf('%g', manual.ll-last.ll))
    cat('ITER_man:', iter, manual.ll, last.ll, manual.ll-last.ll, my.branch, o.t$maximum, '\n')
    cat('ITER_q  :', iter, iter.ll, last.iter.ll, iter.ll-last.iter.ll, my.branch, o.t$maximum, '\n')
    # cat(sprintf('ITER_man: %d %10.10g %10.10g %10.10g %s %f\n', iter, manual.ll, last.ll, manual.ll-last.ll, my.branch, o.t$maximum))
    # cat(sprintf('ITER_q  : %d %10.10g %10.10g %10.10g %s %f\n', iter, iter.ll, last.iter.ll, iter.ll-last.iter.ll, my.branch, o.t$maximum))
    # cat('ITER_q  :', iter, iter.ll, last.iter.ll, iter.ll-last.iter.ll, my.branch, o.t$maximum, '\n')
    cat('\n')
    
    if (manual.ll-last.ll < 0 | iter.ll-last.iter.ll < 0 & fail_on_neg_change_q) {
      cat('WARNING ERROR ETC ETC: likelihood went down in this iteration!', last.ll, manual.ll, manual.ll-last.ll, '\n')
      cat(' -what about iter lik:                                        ', last.iter.ll, iter.ll, iter.ll-last.iter.ll, '\n')
    }
    if (manual.ll-last.ll < 0 & fail_on_neg_change) {
      break()
    }
    if (iter.ll-last.iter.ll < 0 & fail_on_neg_change_q) {
      break()
    }
    if (ll.converge != 0 & ll.converge > abs(manual.ll-last.ll)) break()
    # if (ll.converge != 0 & ll.converge > manual.ll-last.ll) break()
    last.ll <- manual.ll
    last.iter.ll <- iter.ll
  }
  
  ret$dt.theta = dt.theta
  ret$branchtime = o.t$maximum
  ret$branch = my.branch
  ret$max.ll = iter.ll
  ret$man.max.ll = manual.ll
  ret
}

ll_ret_to_dt_sims <- function(my.ret, args) {
  dt.theta = data.table(my.ret$dt.theta)
  dt.theta[, true_mh_contam := args$add_contam]
  dt.theta[, true_faunal_prop := args$add_faunal]
  
  dt = data.table(branch = my.ret$branch,
             true_branch = args$true_branch,
             branchtime = my.ret$branchtime,
             true_branchtime = args$true_branchtime,
             branchtime.rel = (my.ret$branchtime - sims.dat$bounds.for.branch(my.ret$branch, 'low', sims.dat)) / 
               (sims.dat$bounds.for.branch(my.ret$branch, 'high', sims.dat) - sims.dat$bounds.for.branch(my.ret$branch, 'low', sims.dat)),
             max.ll.iter = length(my.ret$ll.trace),
             max.ll = my.ret$max.ll,
             next.ll = tail(my.ret$ll.trace,2)[1],
             lib = paste0(args$libs, collapse = '_'),
             downsample = args$downsample,
             nreads = args$nreads,
             tag = args$tag)
  
  cbind(dt, dt.theta)
}
# ll_ret_to_dt_sims(x.em, args)

ll_ret_to_dt <- function(my.ret, sites.cat) {
  data.table(max.ll = my.ret$max.ll,
             branch = my.ret$branch,
             branchtime = my.ret$branchtime,
             branchtime.rel = (my.ret$branchtime - sims.dat$bounds.for.branch(my.ret$branch, 'low', sims.dat)) / 
               (sims.dat$bounds.for.branch(my.ret$branch, 'high', sims.dat) - sims.dat$bounds.for.branch(my.ret$branch, 'low', sims.dat)),
             mh_contam_deamT = my.ret$dt.theta[rg %like% 'deam_TRUE', mh_contam],
             faunal_prop_deamT = my.ret$dt.theta[rg %like% 'deam_TRUE', faunal_prop],
             mh_contam_deamF = my.ret$dt.theta[rg %like% 'deam_FALSE', mh_contam],
             faunal_prop_deamF = my.ret$dt.theta[rg %like% 'deam_FALSE', faunal_prop],
             max.ll.iter = length(my.ret$ll.trace),
             next.ll = tail(my.ret$ll.trace,2)[1],
             sites.cat = sites.cat, 
             lib = args$libs,
             tag = args$tag)
}
# ll_ret_to_dt(my.ret$a, args$site_cat)

# function(dt.sed.analysis, sims.dat, my.branch, err_rate = 0.001, max.iter = 20, ll.converge = 0, 
#          set.branchtime = 'estim', set.mh_contam = 'estim', set.faunal_prop = 'estim', p_h_method = 'full',
#          fail_on_neg_change = F, fail_on_neg_change_q = F)

sed_EM_allbranch <- function(dt.sed.analysis, sims.dat, branches = NULL, ...) {
  ret = list()
  max.ll = -1e200
  if (is.null(branches)) branches <- sims.dat$branches
  for (my.branch in branches) {
    ret[[my.branch]] <- sed_EM(dt.sed.analysis = dt.sed.analysis, sims.dat = sims.dat, my.branch = my.branch, ...)
    ### ISSUE - manual vs q-ll
    ## here I use the manual ll, because the q-based ll was giving inconsistent results across the branches
    ## i.e., the manual ll may be highest for 'v', and the q-ll highest for anc_1.  I think this is because
    ## the q-ll "remembers" the path it took to get there [in the form of gamma], and so it's not comparable
    ## across branches. I'm not 100% on this, though.
    if (ret[[my.branch]]$man.max.ll > max.ll) ret$max <- ret[[my.branch]]
    max.ll <- ret[[my.branch]]$man.max.ll
  }
  ret
}

# sed_EM_allbranch <- function(dt.sed.analysis, sims.dat, err_rate = 0.001, max.iter = 20) {
#   max.ll = -1e200
#   my.ret <- foreach (my.branch = sims.dat$branches, .combine = c) %dopar% {
#     cat(my.branch, '\n')
#     ret = list()
#     ret[[my.branch]] <- sed_EM(dt.sed.analysis = dt.sed.analysis, sims.dat = sims.dat, my.branch = my.branch, err_rate = err_rate, max.iter = max.iter)
#     ret[[my.branch]]
#   }
#   # if (ret$max.ll > max.ll) ret$max <- ret[[my.branch]]
#   my.ret
# }






######## get a true maximum likelihood surface for branchtime, over a semi-random grid
## the surface w/ ll ~ less than ll.thresh away from the max is explored more extensively with each step (nsteps)
## returns a data.table

grid_t_em_theta <- function(dt.sed.analysis, sims.dat, my.branches, err_rate = 0.001, p_h_method = 'simple',
                            bins.t = 10,
                            max.iter = 30, ll.converge = 1e-6,
                            nsteps = 3, print.debug = F, ll.thresh = 4) {
  
  ## first get the true max likelihood estimate of theta and t
  hey.em <- sed_EM_allbranch(dt.sed.analysis, sims.dat,
                             branches = my.branches,
                             err_rate = 0.001, set.branchtime = 'estim',
                             max.iter = max.iter, ll.converge = ll.converge,
                             set.faunal_prop = 'estim',
                             set.mh_contam = 'estim',
                             p_h_method = 'simple')
  hey.em <- hey.em$max
  
  dt.ret <- data.table(hey.em$dt.theta)
  dt.ret[, branchtime := hey.em$branchtime]
  dt.ret[, branch := hey.em$branch]
  dt.ret[, man.max.ll := hey.em$man.max.ll]
  dt.ret[, man.max.ll.last := tail(hey.em$man.ll.trace,1)]
  dt.ret[, n.iter := length(hey.em$man.ll.trace)]
  dt.ret[, my.t.idx := 0]
  dt.ret[, max.ll := man.max.ll]
  dt.ret[, step.x := 0]
  
  for (step.x in 1:nsteps) {
    cat(sprintf('\nstep %d/%d\n', step.x, nsteps))
    
    x.manll.dt <- foreach (my.b = my.branches, .combine = rbind) %:%
      foreach(my.t.idx = 1:bins.t, .combine = rbind) %dopar% {
        
        if (step.x == 1) {
          vals.t = sims.dat$branch.bounds[my.b, unique(seq(t.low, t.high, length.out = bins.t)), on='branch']
        } else {
          # print()
          cat(sprintf('\nstep %d; branch %s; idx %d\n', step.x, my.b, my.t.idx))
          # print()
          ## sometimes the ll are very slightly different, so uniq doesn't work.  so I take the mean per branchtime, thus guaranteeing that each branchtime is represented only once (but preserving the ll)
          dt.x.branch = setorder(unique(dt.ret[branch == my.b, .(branchtime, man.max.ll = mean(man.max.ll), max.ll), branchtime]), branchtime)
          print(dt.x.branch)
          a = dt.x.branch[, max.ll-man.max.ll < ll.thresh]
          print(a)
          if (sum(a) == 0) {
            print(sprintf('branch has no in-bound points: %s', my.b))
            return(data.table()) ## there are no points on this branch that are within the bounds
          }
          a.new <- a
          a.new[-1] <- a.new[-1] | a[-length(a)]
          a.new[-length(a)] <- a.new[-length(a)] | a[-1]
          a
          a.new
          print(a.new)
          print(dt.x.branch[a.new])
          # vals.t = dt.x.branch[a.new, unique(seq(min(branchtime), max(branchtime), length.out = bins.t+1))] ## the endpoints are included, but we have already sampled them
          vals.t = dt.x.branch[a.new, runif(bins.t, min(branchtime), max(branchtime))] ## the endpoints are included, but we have already sampled them
          print(vals.t)
          vals.t <- vals.t[!vals.t %in% dt.x.branch$branchtime]
          print(vals.t)
        }
        if (my.t.idx > length(vals.t)) {
          print(sprintf('idx is larger than length of vals.t: %s, %d > %d', my.b, my.t.idx, length(vals.t)))
          return(data.table())
        }
        
        my.t = vals.t[my.t.idx]
        if (print.debug) cat(mh_contam, faunal_prop, my.t, '\n')
        hey.em <- sed_EM(dt.sed.analysis, sims.dat,
                         my.branch = my.b,
                         err_rate = 0.001, set.branchtime = my.t,
                         max.iter = max.iter, ll.converge = ll.converge,
                         set.faunal_prop = 'estim',
                         set.mh_contam = 'estim',
                         p_h_method = 'simple')
        
        dt.tmp <- data.table(hey.em$dt.theta)
        dt.tmp[, branchtime := hey.em$branchtime]
        dt.tmp[, branch := hey.em$branch]
        dt.tmp[, man.max.ll := hey.em$man.max.ll]
        dt.tmp[, man.max.ll.last := tail(hey.em$man.ll.trace,1)]
        dt.tmp[, n.iter := length(hey.em$man.ll.trace)]
        dt.tmp[, my.t.idx := my.t.idx]
        dt.tmp[, max.ll := NA_real_]
        dt.tmp[, step.x := step.x]
        
      }
    
    dt.ret <- rbind(x.manll.dt, dt.ret)
    dt.ret[, max.ll := max(man.max.ll)]
  }
  
  dt.ret
}

if (F) {
  dt.x.new = grid_t_em_theta(dt.sed.analysis.archaics[lib == 'Mez1_R5661'],
                             sims.dat.archaics,
                             # my.branches = c('v'), bins.t = 5,
                             my.branches = c('v','anc_1','c'), bins.t = 5,
                             max.iter = 30, ll.converge = 1e-6, nsteps = 6, ll.thresh = 10)
  ggplot(dt.x.new, aes(x=branchtime, y=man.max.ll, color=branch)) + geom_point() + geom_line() + geom_hline(aes(yintercept = max(man.max.ll)-2)) #+ facet_wrap(~step.x)
  ggplot(dt.x.new, aes(x=branchtime, y=man.max.ll, color=branch)) + geom_point() + geom_line() + geom_hline(aes(yintercept = max(man.max.ll)-10)) + facet_wrap(~step.x)
}










###########################33
## this just runs some analyses that approach what we would want to do with each sample


run_simple_analysis <- function(dt.sed.analysis, sims.dat, branches = c('c','v','anc_1'), 
                                blocks = 10, nbootstraps = 0, max.iter = 30,
                                set.faunal_prop = 'estim',
                                set.mh_contam = 'estim') {
  dt.sed.analysis.em <- sed_EM_allbranch(dt.sed.analysis,
                                         sims.dat, 
                                         branches = branches,
                                         err_rate = 0.001,
                                         max.iter = max.iter, ll.converge = 1e-6,
                                         set.faunal_prop = set.faunal_prop,
                                         set.mh_contam = set.mh_contam,
                                         p_h_method = 'simple')
  
  dt.sed.analysis.em.theta <- merge(dt.sed.analysis[, .N, rg], dt.sed.analysis.em$max$dt.theta)
  dt.sed.analysis.em.theta[, sum(mh_contam * N / sum(N))]
  
  dt.sed.analysis.gridll = sed_grid_search(dt.sed.analysis,
                                           sims.dat, my.branch = branches,
                                           err_rate = 0.001, p_h_method = 'simple', nsteps = 2,
                                           range.mh_contam = dt.sed.analysis.em.theta[, sum(mh_contam * N / sum(N))],
                                           range.faunal_prop = dt.sed.analysis.em.theta[, sum(faunal_prop * N / sum(N))])
  
  p1 <- ggplot(dt.sed.analysis.gridll, aes(x=my.t, y=ll, color=my.branch)) + geom_line() + 
    geom_point(data=dt.sed.analysis.gridll[is.max == T]) + ggtitle('grid search full likelihood')
  p1
  
  ret <- list()
  ret$em <- dt.sed.analysis.em
  ret$grid <- dt.sed.analysis.gridll
  ret$em.theta <- dt.sed.analysis.em.theta
  ret$p1 <- p1
  
  if (nbootstraps > 0) {
    dt.sed.analysis.gridll.bootstraps <- foreach(iter = 1:nbootstraps, .combine = rbind) %do% {
      cat('\n\nBootstrap', iter, '\n')
      dt = sed_grid_search(bootstrap_gt_data(dt.sed.analysis, blocks),
                           sims.dat, my.branch = branches,
                           err_rate = 0.001, p_h_method = 'simple', nsteps = 2,
                           range.mh_contam = dt.sed.analysis.em.theta[, sum(mh_contam * N / sum(N))],
                           range.faunal_prop = dt.sed.analysis.em.theta[, sum(faunal_prop * N / sum(N))])
      dt[, bs.iter := iter]
      dt
    }
    p1.boots <- ggplot(dt.sed.analysis.gridll, aes(x=my.t, y=ll, color=my.branch)) + geom_line() + 
      geom_point(data=dt.sed.analysis.gridll.bootstraps[is.max == T], 
                 y=dt.sed.analysis.gridll[is.max == T, ll], alpha=.5) +
      ggtitle('grid search full likelihood')
    ret$p1.boots <- p1.boots
    ret$grid.bootstraps <- dt.sed.analysis.gridll.bootstraps
  }
  ret
}




#################################
## the original code for calculating a composite likelihood
#################################
# 
# # test.gts <- fread('~/Google Drive/soil_dna_capture/all_simple_gts.deam.tsv.gz')
# # test.gts <- fread('output_v4/all_simple_gts.deam.tsv.gz')
# # test.gts <- fread('zcat output_v4/all_simple_gts.deam.tsv.gz | head -n1000')
# 
# 
# dt.sims.p = fread('~/GoogleDrive/soil_dna_capture/simulate_demog/dt.sims.p.simfiles004.txt')
# setnames(dt.sims.p, c('v', 'c', 'a', 'sed'), c('v_gt', 'c_gt', 'a_gt', 'sed_gt'))
# dt.sims.p.coarse = dt.sims.p[!endsWith(paste(time), '01') | endsWith(paste(time), '001')]

apply_contam <- function(dt, contam=0.05) {
  dt = data.table(dt)
  ## contam % of the time, when we should see a 1, we'll see a 0
  dt[sed == 1, p := p * (1-contam)]
  ## also increase 0 p by the same amount
  dt[sed == 0, p := p + (1-p) * contam]
  dt[, model.contam := contam]
  return(dt)
}
# a = apply_contam(dt.sims.p)

demog_composite_liks <- function(dt.sims.p, test.gts, model.contam = NULL) {
  
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
  
  test.gts <- data.table(test.gts)
  dt.sims.p <- data.table(dt.sims.p)
  
  if (!'v_gt' %in% colnames(test.gts)) {
    # also save a1/a2/pantro for pairing w/ MH freqs?
    test.gts <- test.gts[!is.na(read_base.gt) &
                           !is.na(vind.gt) &
                           !is.na(chag.gt) &
                           !is.na(altai.gt),
                         .(v_gt = vind.gt, c_gt = chag.gt, a_gt = altai.gt, 
                           sed_gt = read_base.gt/2,
                           lib, libname)]
  }
  
  test.gts = test.gts[!(v_gt == c_gt & v_gt == a_gt)]
  
  if (!is.null(model.contam)) {
    if (!(length(model.contam) > 0 && is.numeric(model.contam))) model.contam = seq(0,.3,.05)
    dt.sims.p = do.call(rbind, lapply(model.contam, function(x) apply_contam(dt.sims.p, x)))
  } else {
    dt.sims.p[, model.contam := 0]
  }
  
  setkey(dt.sims.p, time, branch, model.contam)
  sim.vals = dt.sims.p[, .(time,branch,model.contam)] %>% unique
  test.gts.t = foreach(x = sim.vals[,.I], .combine = rbind) %do% {
    x = sim.vals[x]
    print(x)
    x = merge(test.gts, dt.sims.p[x], by=c('v_gt','c_gt','a_gt','sed_gt'))
    x[, .(ll.og=sum(log(p)),
          ll=log(choose(.N,24))+sum(log(p))), .(time,branch,model.contam,lib)]
  }
  
  test.gts.t[, ll.max := max(ll), lib]
  test.gts.t[, ll.k := exp(ll.max-ll)]
  
  ret = list()
  ret[['full.tree']] <- test.gts.t
  ret[['max.model']] <- test.gts.t[ll.max == ll]
  ll.max.contam = ret[['max.model']][, model.contam]
  # ret[['k20.ci']] <- test.gts.t[ll.k < 20 & model.contam == ll.max.contam, range(time), branch]
  # ret[['k20.contam.ci']] <- test.gts.t[ll.k < 20, range(model.contam)]
  return(ret)
}

if (F) {
  libs.150 <- test.gts[, .N, lib][N > 150, lib]
  cat('a.s.lib\n')
  a.s.lib <- test.gts[lib %in% libs.150, demog_liks(dt.sims.p, .SD, model.contam = T)$full.tree, .(lib,libname)]
  cat('a.s.libname\n')
  a.s.libname <- test.gts[lib %in% libs.150, demog_liks(dt.sims.p, .SD, model.contam = T)$full.tree, libname]
  date()
  
  
  dt.sims.p = fread('~/GoogleDrive/soil_dna_capture/simulate_demog/dt.sims.p.simfiles004.txt')
  setnames(dt.sims.p, c('v', 'c', 'a', 'sed'), c('v_gt', 'c_gt', 'a_gt', 'sed_gt'))
  lik.res <- demog_composite_liks(dt.sims.p, dt.sed.poly[, .(v_gt, c_gt, a_gt, sed_gt, lib)])
  ggplot(lik.res$full.tree, aes(x=time, y=ll, color=branch, lty=lib)) + geom_line()
  
  # $max.model
  # time branch model.contam        lib     ll.og        ll    ll.max ll.k
  # 1: 0.6301      v            0 Mez2_A9180 -672.0432 -554.3336 -554.3336    1
  # 2: 0.7801      c            0 Mez1_R5661 -634.6606 -518.4100 -518.4100    1
  # 0.6301 is the lowest "real" point on the vindija lineage
  # ggplot(dt.sims.p[sed_gt == 1 & v_gt == 1 & c_gt == 2 & a_gt == 2], aes(time, p, color=branch)) + geom_line() + xlim(.5,.7) + geom_vline(xintercept = 0.6301)

  sims.dat.dt.sims.p.i10 <- rbind(sims.dat$dt.sims.p[, .(sed = 0, p = 1-sum(n_der)/sum(n_tot)), .(v, c, a, branch, time)],
                                  sims.dat$dt.sims.p[, .(sed = 1, p =   sum(n_der)/sum(n_tot)), .(v, c, a, branch, time)])
  
  setnames(sims.dat.dt.sims.p.i10, c('v', 'c', 'a', 'sed'), c('v_gt', 'c_gt', 'a_gt', 'sed_gt'))
  lik.res <- demog_composite_liks(sims.dat.dt.sims.p.i10, dt.sed.poly[sample(.N, 5000, replace = T), .(v_gt, c_gt, a_gt, sed_gt, lib)])
  lik.res <- demog_composite_liks(sims.dat.dt.sims.p.i10, dt.sed.poly[, .(v_gt, c_gt, a_gt, sed_gt, lib)])
  ggplot(lik.res$full.tree, aes(x=time, y=ll, color=branch, lty=lib)) + geom_line()
  
  # $max.model
  # time branch model.contam        lib     ll.og        ll    ll.max ll.k
  # 1: 0.5585724      v            0 Mez2_A9180 -678.2588 -560.5491 -560.5491    1
  # 2: 0.7801000      c            0 Mez1_R5661 -637.0174 -520.7667 -520.7667    1
  ## with the newer sims [only looking at polymorphic neandertal sites, completely disregarding denisova]
  ## I get similar point estimates, namely that the best branch time of mez1 is 0.7801 on chag,
  ##  and the best branch time of Mez2 is ~ the tip of vindija
  ## BUT the likelihood surface looks much "better" for older sims, in that Mez1 is pretty much at the branchpoint,
  ## and Mez2 is close to the tip of Vindija
  ## soooo maybe differences in simulations are not so much of an issue?
  
  
  sims.dat.dt.sims.p.i10.forhets <- rbind(sims.dat$dt.sims.p[, .(sed = 0), .(v, c, a, branch, time)],
                                  sims.dat$dt.sims.p[, .(sed = 1), .(v, c, a, branch, time)])
  test.v_hets_only <- rbind(sims.dat$dt.sims.p[, .(sed = 0, p = 1 - sum(n_der)/sum(n_tot)), .(v, branch, time)],
                            sims.dat$dt.sims.p[, .(sed = 1, p =     sum(n_der)/sum(n_tot)), .(v, branch, time)])
  sims.dat.dt.sims.p.i10.forhets <- merge(sims.dat.dt.sims.p.i10.forhets, test.v_hets_only, by=c('v', 'branch', 'time', 'sed'))
  setnames(sims.dat.dt.sims.p.i10.forhets, c('v', 'c', 'a', 'sed'), c('v_gt', 'c_gt', 'a_gt', 'sed_gt'))
  lik.res <- demog_composite_liks(sims.dat.dt.sims.p.i10.forhets, dt.sed.poly[, .(v_gt, c_gt, a_gt, sed_gt, lib)])
  # $max.model
  # time branch model.contam        lib     ll.og        ll    ll.max ll.k
  # 1: 0.7401      v            0 Mez2_A9180 -721.5131 -603.8035 -603.8035    1
  # 2: 1.3201  anc_2            0 Mez1_R5661 -742.7500 -626.4994 -626.4994    1
  ## hets-only gives quite old estimates for Mez1
  
  ggplot(lik.res$full.tree, aes(x=time, y=ll, color=branch, lty=lib)) + geom_line()
  # ggplot(lik.res$full.tree, aes(x=time, y=ll.og, color=branch, lty=lib)) + geom_line() ## doesn't change at all
  
  dt.sed.poly.tst <- data.table(dt.sed.poly)
  dt.sed.poly.tst <- dt.sed.poly.tst[!(v_gt == c_gt & v_gt == a_gt)]
  dt.sed.poly.tst[, rg := 'hey']
  x.em = sed_EM(dt.sed.poly.tst[lib %like% 'Mez1'], sims.dat, my.branch = 'v', err_rate = 0.001,
                max.iter = 10, ll.converge = 1e-10,
                set.faunal_prop = 0,
                set.mh_contam = 0,
                p_h_method = 'simple')
  # Mez1: ITER_man: 4 -640.8327 -640.8327 4.547474e-13 v 0.7625362 
  
  x.em = sed_EM(dt.sed.poly.tst[lib %like% 'Mez1'], sims.dat, my.branch = 'c', err_rate = 0.001,
                max.iter = 10, ll.converge = 1e-10,
                set.faunal_prop = 0,
                set.mh_contam = 0,
                p_h_method = 'simple')
  # Mez1: ITER_man: 5 -634.0518 -634.0518 -2.569323e-11 c 0.7796873

  
  
  ###
  ## compare new simulations from old simulation script with og simultions
  sims.dat.newog <- generate_sims_dat('data/sims_og_newrun.txt') ###

  ########  THIS DOESN'T WORK BECAUSE SIMPLE_P DOESN'T "REMEMBER" N_DER AND N_TOT, SO YOU CAN'T RECOMPUTE P WITHOUT DEN
  ## WOULD HAVE TO MAKE THE GENERATE FUNCTION SO IT CAN TAKE A SET OF COLUMNS TO USE FOR GENERATING THE DATA


  
  dt.new_og_sims.raw <- fread('data/sims_og_newrun.txt')
  dt.new_og_sims.raw[, gt.category := paste0(v, c, a)]
  dt.new_og_sims.raw.polyn <- rbind(dt.new_og_sims.raw[, .(sed = 0, p = 1-sum(n_der)/sum(n_tot)), .(v, c, a, gt.category, branch, time)],
                                    dt.new_og_sims.raw[, .(sed = 1, p =   sum(n_der)/sum(n_tot)), .(v, c, a, gt.category, branch, time)])
  
  setnames(dt.new_og_sims.raw.polyn, c('v', 'c', 'a', 'sed'), c('v_gt', 'c_gt', 'a_gt', 'sed_gt'))
  # lik.res <- demog_composite_liks(dt.new_og_sims.raw.polyn, dt.sed.poly[sample(.N, 5000, replace = T), .(v_gt, c_gt, a_gt, sed_gt, lib)])
  lik.res <- demog_composite_liks(dt.new_og_sims.raw.polyn, dt.sed.poly[, .(v_gt, c_gt, a_gt, sed_gt, lib)])
  ggplot(lik.res$full.tree, aes(x=time, y=ll, color=branch, lty=lib)) + geom_line()
  ## this is just *very* coarse, but it does seem to point that my new simulations accurately represent 
  ## the old simulations, and give the "correct" answers for Mez1 and Mez2
  
  sims.dat.newog.onlyn <- generate_sims_dat('data/sims_og_newrun.txt', agCols = c('v_gt', 'c_gt', 'a_gt'))
  
  dt.new_og_sims <- sims.dat.newog.onlyn$branch.bounds[, .(time = c(seq(t.low, t.high, .02), t.high)), branch]
  dt.new_og_sims <- merge(unique(sims.dat.newog.onlyn$dt.sims.p[, .(v_gt, c_gt, a_gt, gt.category, branch)]), dt.new_og_sims, by='branch', allow.cartesian = T)

  dt.new_og_sims <- rbind(dt.new_og_sims[, .(sed_gt = 1, p = sims.dat.newog.onlyn$simple_p_given_b_t_arcs(branch, time, gt.category, sims.dat.newog.onlyn)),
                                         by = .(branch, time, v_gt, c_gt, a_gt, gt.category)],
                          dt.new_og_sims[, .(sed_gt = 0, p = 1-sims.dat.newog.onlyn$simple_p_given_b_t_arcs(branch, time, gt.category, sims.dat.newog.onlyn)),
                                         by = .(branch, time, v_gt, c_gt, a_gt, gt.category)])
  lik.res <- demog_composite_liks(dt.new_og_sims, dt.sed.poly[, .(v_gt, c_gt, a_gt, sed_gt, lib)])
  p1 <- ggplot(lik.res$full.tree, aes(x=time, y=ll, color=branch, lty=lib)) + geom_line() + 
    geom_point(data=lik.res$full.tree[ll == ll.max]) + ggtitle('composite likelihood')
  p1
  # lik.res <- demog_composite_liks(dt.new_og_sims, bootstrap_gt_data(dt.sed.poly[, .(v_gt, c_gt, a_gt, sed_gt, lib)], blocks = 10, drop.one = 5))
  # ggplot(lik.res$full.tree, aes(x=time, y=ll, color=branch, lty=lib)) + geom_line()
  
  dt.sed.poly.tst <- data.table(dt.sed.poly)
  dt.sed.poly.tst[, rg := 'hey']
  dt.sed.poly.tst[, gt.category := paste0(v_gt, c_gt, a_gt)]
  dt.sed.poly.tst.gridll = sed_grid_search(dt.sed.poly.tst, sims.dat.newog.onlyn, my.branch = c('v','c','anc_1', 'anc_2', 'anc_3', 'a', 'd'),
                                           err_rate = 0.001, p_h_method = 'simple', nsteps = 2,
                                           range.mh_contam = 0.025, range.faunal_prop = 0)
  dt.sed.poly.tst.gridll
  p2 <- ggplot(dt.sed.poly.tst.gridll, aes(x=my.t, y=ll, color=my.branch)) + geom_line() + 
    geom_point(data=dt.sed.poly.tst.gridll[is.max == T]) + ggtitle('grid search full likelihood contam=0.025')
  p2
  
  # ggsave('figures/demog_newog_composite_lik_mez1.pdf')
  
  
  
  ###############################333
  ###############################333 debugging ^^^
  ###############################333
  
  
  
  
  ## maybe I continue with my programming project?  or first just write code to modify these into "simple"
  ## linear models for the existing EM code?
  
  #
  
  
  sims.dat.newog <- generate_sims_dat('data/sims_og_newrun.txt')
  dt.sed.poly.tst <- data.table(dt.sed.poly)
  dt.sed.poly.tst[, gt.category := paste0(v_gt, c_gt, a_gt, d_gt)]
  # dt.sed.poly.tst <- dt.sed.poly.tst[!(v_gt == c_gt & v_gt == a_gt)]
  # dt.sed.poly.tst[lib %like% 'Mez1']
  dt.sed.poly.tst[, rg := 'hey']
  # dt.sed.analysis[, gt.category := gt.cat]
  x.em = sed_EM(dt.sed.analysis[lib %like% 'Mez1'], sims.dat.newog, my.branch = 'v', err_rate = 0.001,
                max.iter = 40, ll.converge = 1e-10,
                set.faunal_prop = 0,
                set.mh_contam = 'estim',
                p_h_method = 'simple')
  
  
}





#################
## test method by "simulating" data
  

simulate_data_sed <- function(method, my.t, sims.dat, my.mh_contam = 0, my.faunal_prop = 0, err_rate = 0.001, nsites = 1002) {
  
  dt.sed.poly.tst <- data.table(v_gt=0:2, c_gt=0:2, a_gt=0:2, d_gt=0, f_mh=rbeta(nsites, .2, 5), rg = 'hey_rg', x = 1:nsites)
  dt.sed.poly.tst[, v_gt := sample(v_gt)]
  dt.sed.poly.tst[, c_gt := sample(c_gt)]
  dt.sed.poly.tst[, a_gt := sample(a_gt)]
  
  if (method == 'bt')  dt.sed.poly.tst[, p_h_der := my.t]
  if (method == 'bt2') dt.sed.poly.tst[, p_h_der := my.t / (c_gt+1)]
  if (method == 'bt3') dt.sed.poly.tst[, p_h_der := ifelse(a_gt == 0, 0.001, my.t / (c_gt+1))]
  
  if (method == 'simple') {
    dt.sed.poly.tst[, p_h_der := sims.dat$simple_p_given_b_t_arcs('v', my.t,
                                                                  gt.category,
                                                                  sims.dat),
                    gt.category]
  }
  if (method == 'full') {
    dt.sed.poly.tst[, p_h_der := sims.dat$p_gt_given_b_t_arcs('v', my.t,
                                                              gt.category,
                                                              sims.dat),
                    gt.category]
  }
  
  
  dt.sed.poly.tst[, sed_gt := sample(c(0,1), .N, prob = c(1-p_h_der,p_h_der), replace=T), p_h_der]
  
  dt.sed.poly.tst[, .(sum(sed_gt) / .N, .N), .(p_h_der)]
  
  mh.sites = dt.sed.poly.tst[, sample(.N, .N*my.mh_contam)]
  dt.sed.poly.tst[mh.sites, .N, keyby=sed_gt]
  dt.sed.poly.tst[mh.sites, sed_gt := sample(c(0,1), .N, prob=c(1-f_mh,f_mh), replace=T), f_mh]
  dt.sed.poly.tst[mh.sites, .(sum(sed_gt)/.N, mean(f_mh))]
  #
  err.sites = dt.sed.poly.tst[, sample(.N, .N*err_rate)]
  dt.sed.poly.tst[err.sites, .N, sed_gt]
  dt.sed.poly.tst[err.sites, sed_gt := sample(c(0,1), .N, prob=c(.5,.5), replace=T)]
  dt.sed.poly.tst[err.sites, .N, sed_gt]
  dt.sed.poly.tst
}
  
if (F) {
  simulate_data_sed('bt', .75, sims.dat)
  simulate_data_sed('bt2', .75, sims.dat)
  simulate_data_sed('bt3', .75, sims.dat)
  simulate_data_sed('simple', .75, sims.dat)[v_gt == 0 & c_gt == 0 & a_gt == 0]
  simulate_data_sed('full', .75, sims.dat)[v_gt == 0 & c_gt == 0 & a_gt == 0]
}
  


eval_sed_t_and_mh <- function(dt.sed.poly.tst, p_h_method, max.iter = 100, set.faunal_prop = 0) {
  ## constrain mh and faunal = 0, just estimate time with EM
  x.em = sed_EM(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, max.iter = 4, set.faunal_prop = 0, set.mh_contam = 0, p_h_method = p_h_method)
  
  ## the ll goes up and down a bit in the trace sometimes?  just not 100% stable?
  # plot(tail(x.em$man.ll.trace,9))
  # plot(tail(x.em$ll.trace,9))
  
  ## do "grid" search, just estimating time  
  dt.sed.poly.tst.gridll = sed_grid_search(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, p_h_method = p_h_method,
                                           bins.mh_contam = 1, bins.faunal_prop = 1, bins.t = 10, nsteps = 2,
                                           range.mh_contam = c(0,0), range.faunal_prop = c(0,0), 
                                           range.t = sims.dat$branch.bounds[branch == 'v', c(t.low, t.high)])
  
  ## do EM search, estimating time and MH contam
  x.em.all_params.n100 = sed_EM(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, 
                                max.iter = max.iter, ll.converge = 1e-10, 
                                set.faunal_prop = set.faunal_prop, p_h_method = p_h_method)
  
  dt.sed.poly.tst.gridll.all_params = sed_grid_search(dt.sed.poly.tst, sims.dat, my.branch = 'v', err_rate = 0.001, p_h_method = p_h_method,
                                                      bins.mh_contam = 10, bins.faunal_prop = 1, bins.t = 10, nsteps = 3,
                                                      range.mh_contam = c(0,.2), range.faunal_prop = c(0,0),
                                                      range.t = sims.dat$branch.bounds[branch == 'v', c(t.low, t.high)])
  em.mismatch = dt.sed.poly.tst.gridll.all_params[is.max == T, ll] - x.em.all_params.n100$man.max.ll
  return(list(x.em = x.em,
              dt.sed.poly.tst.gridll = dt.sed.poly.tst.gridll,
              x.em.all_params.n100 = x.em.all_params.n100,
              dt.sed.poly.tst.gridll.all_params = dt.sed.poly.tst.gridll.all_params,
              em.mismatch = em.mismatch))
}

plot_eval_sed_t_and_mh <- function(eval.ret, true_branchtime = 0.75, true_mh_contam = 0, plot.0 = F, plot.title = NULL) {
  ## constrain mh and faunal = 0, just estimate time with EM
  
  if (plot.0) {
    ## plot them, they match!
    p1 <- ggplot(eval.ret$dt.sed.poly.tst.gridll, aes(x=my.t, y=ll)) + geom_point() +
      geom_point(aes(x=eval.ret$x.em$branchtime, y=eval.ret$x.em$man.max.ll), color='red', pch='x', size=10) +
      xlab('branch time estimate grid vs EM (red x)') +
      geom_vline(xintercept = eval.ret$dt.sed.poly.tst.gridll[is.max == T, my.t], color='green', lty=3) +
      geom_vline(xintercept = true_branchtime)
    if (!is.null(plot.title)) {
      p1 <- p1 + ggtitle(plot.title)
    }
    print(p1)
  }
  
  p2 <-   ggplot(eval.ret$dt.sed.poly.tst.gridll.all_params[step.x > 0],
                 aes(x=my.t, y=mh_contam, color=ll)) +
    geom_point() +
    geom_point(data=eval.ret$dt.sed.poly.tst.gridll.all_params[max(ll) - ll < eval.ret$em.mismatch & step.x > 0], color='red') +
    geom_point(data=eval.ret$dt.sed.poly.tst.gridll.all_params[max(ll) - ll < eval.ret$em.mismatch/10 & step.x > 0], color='green') +
    geom_path(data=data.table(my.t = rep(eval.ret$x.em.all_params.n100$branchtime.trace,
                                         each = length(unique(eval.ret$x.em.all_params.n100$dt.theta.trace$rg))),
                              eval.ret$x.em.all_params.n100$dt.theta.trace),
              aes(group = rg),
              color='red') +
    geom_point(data = data.table(x=rep(true_branchtime, length(true_mh_contam)), y=true_mh_contam), 
               aes(x=x,y=y), color='black', pch='x', size=10) +
    geom_point(data = data.table(x=rep(eval.ret$x.em.all_params.n100$branchtime, 
                                       length(eval.ret$x.em.all_params.n100$dt.theta$mh_contam)), 
                                 y=eval.ret$x.em.all_params.n100$dt.theta$mh_contam),
               aes(x=x,y=y), color='slateblue1', pch='x', size=10) +
    geom_point(data=eval.ret$dt.sed.poly.tst.gridll.all_params[max(ll) == ll], color='black', pch='*', size=4) +
    xlab('branch time estimate grid (green/red dots) vs EM (red x)') + ylab('mh contam estimate') +
    coord_cartesian(ylim = c(0,.3)) +
    ggtitle(sprintf('ll of red points within %g of EM ll; red line is EM trace', eval.ret$em.mismatch))
  
  if (!is.null(plot.title)) {
    p2 <- p2 + ggtitle(plot.title)
  }
  
  # dt.sed.analysis.mh_all.contam1.x4.one_rg2.simple.fixed_anc_p_0.004$dt.sed.poly.tst.gridll.all_params[ll == max(ll)]
  #
  print(p2)
}



  
if (F) {
    
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



