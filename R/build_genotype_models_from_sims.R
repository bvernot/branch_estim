library(ggplot2)
library(data.table)
library(dplyr)
library(cobs)
library(foreach)
library(bit64)


dt.sims.p.all <- data.table(file = sprintf('~/Google Drive/soil_dna_capture/simulate_demog/scrm%d.summary.txt', 1:2))
dt.sims.p.all <- dt.sims.p.all[, fread(file), file]
dt.sims.p.all[, file := NULL]
setkey(dt.sims.p.all, v, c, a, d)
ggplot(dt.sims.p.all[list(2,1,1,1)], aes(x=time, y=p, color=branch, lty=factor(scrm_iter))) +
  geom_line()# + geom_smooth(se = F)
#




dt.sims.p = fread('~/Google Drive/soil_dna_capture/simulate_demog/dt.sims.p.simfiles004.txt')
dt.sims.p = fread('~/Google Drive/soil_dna_capture/simulate_demog/sims_summary.txt')

simsfile <- '~/Google Drive/soil_dna_capture/simulate_demog/scrm1.sims_summary.txt'


dt.sims.p = fread(simsfile)
setkey(dt.sims.p, v, c, a, d)

all.branches = dt.sims.p[, .(time = max(time)), branch]
setkey(all.branches, time)
trunk.branch <- all.branches[.N, branch]
x.branch = ''
for(x in 1:all.branches[, .N-1]) {
  t1 = all.branches[x, time]
  t2 = all.branches[x+1, time]
  dt.sims.p.tmp <- dt.sims.p[time == t1 & (branch == trunk.branch | branch == x.branch)]
  x.branch = sprintf('anc_%d', x)
  print(dt.sims.p.tmp)
  dt.sims.p.tmp[, branch := x.branch]
  dt.sims.p[branch == trunk.branch & time > t1 & time <= t2, branch := x.branch]
  dt.sims.p <- rbind(dt.sims.p, dt.sims.p.tmp)
}
setkey(dt.sims.p, v, c, a, d)
branch.bounds <- dt.sims.p[, .(t.low = min(time), t.high = max(time)), branch]
branch.bounds[, constr := 'high']
branch.bounds[branch %like% 'anc', constr := 'both']
branch.bounds[branch == x.branch, constr := 'low']
inner.node.times <- branch.bounds[branch %like% 'anc', sort(t.low)]



ggplot(dt.sims.p[list(2,1,1,1)], aes(x=time, y=p, color=branch)) + geom_line()# + geom_smooth(se = F)
ggplot(dt.sims.p[list(1,2,1,1)], aes(x=time, y=p, color=branch)) + geom_line()# + geom_smooth(se = F)
ggplot(dt.sims.p[list(2,2,1,0)], aes(x=time, y=p, color=branch)) + geom_line()# + geom_smooth(se = F)
ggplot(dt.sims.p[list(2,2,2,0)], aes(x=time, y=p, color=branch)) + geom_line()# + geom_smooth(se = F)


check.range <- function(r, x) {
  x >= r[1] & x <= r[2]
}
check.range(c(1,2), seq(0,5))

constraints.for.branch <- function(init.constraints, b) {
  if (length(init.constraints) == 1) init.constraints <- rep(init.constraints, length(inner.node.times))
  bc <- branch.bounds[branch == b]
  if (bc[, constr] == 'high') {
    return(matrix(c(0, bc$t.high, init.constraints[bc$t.high == inner.node.times]), ncol=3))
  }
  if (bc[, constr] == 'low') {
    return(matrix(c(0, bc$t.low, init.constraints[bc$t.low == inner.node.times]), ncol=3))
  }
  if (bc[, constr] == 'both') {
    return(matrix(c(0, bc$t.low, init.constraints[bc$t.low == inner.node.times],
                       0, bc$t.high, init.constraints[bc$t.high == inner.node.times]),
                    ncol = 3, byrow = T))
  }
}
constraints.for.branch(.5, 'v')
constraints.for.branch(.5, 'anc_1')
constraints.for.branch(.5, 'anc_2')
constraints.for.branch(.5, 'anc_3')


dt <- dt.sims.p[list(2,1,1,1)]
dt <- dt.sims.p[list(1,2,1,1)]


## set number of knots based on length of branch?

optim.fn2 <- function(init.constraints, dt, do.plots=0, nknots = 5, gen.fn = F, ...) {

  cat(init.constraints)
  if (do.plots == 2) dt.predict <- data.table()
  if (gen.fn) models.env <- new.env(parent = emptyenv())

  x <- foreach (b = branch.bounds$branch, .combine = rbind) %do% {
    dt.b <- dt[branch == b]
    # b.nknots = max(2,as.integer(dt.b[, diff(range(time))] / dt[, diff(range(time))] * nknots + 1))
    # b.nknots = 5
    constraints.b <- constraints.for.branch(init.constraints, b)
    # cat('modeling branch', b, 'with nknots=', b.nknots, '\n')
    # print(constraints.b)
    model <- with(dt.b, cobs(time, p, pointwise=constraints.b, print.warn = F, print.mesg = F, nknots = nknots, ...))
    if (do.plots == 1) plot(model, main=b)
    if (do.plots == 1) points(constraints.b[, 2], constraints.b[, 3])
    if (do.plots == 2) dt.predict <- 
      rbind(dt.predict, data.table(time=model$x, p=model$fitted, branch = b))

    ## save the model for each branch
    if (gen.fn) {
      models.env[[b]] <- model
    }
    sum(model$resid^2)
  }
  if (gen.fn) {
    fn <- function(my.b, my.t) {
      if (branch.bounds[branch == my.b, my.t < t.low | my.t > t.high]) return(0)
      model <- models.env[[my.b]]
      predict(model, my.t)[,2]
    }
    return(fn)
  }
  
  if (do.plots == 2) {
    p <- ggplot(dt.predict, aes(x=time, y=p, color=branch)) + 
      geom_point(data=dt, alpha=.5) +
      geom_line() +
      geom_point(data=data.table(p=init.constraints, time=inner.node.times), color='red')
      NULL
    print(p)
  }
  cat(':', sum(x), '\n')
  sum(x)
}

# optim.fn2(rep(1,length(inner.node.times)), do.plots = T)
optim.fn2(rep(1,length(inner.node.times)), dt, do.plots = 2)
fn <- optim.fn2(rep(1,length(inner.node.times)), dt, gen.fn = T)
fn('v', .7)

# my.constraints <- optim(rep(.5,length(inner.node.times)), optim.fn2, control = list(reltol = 1e-10), nknots = 10, maxiter = 200)
my.constraints2 <- optim(rep(.5,length(inner.node.times)), optim.fn2, dt = dt)
optim.fn2(my.constraints2$par, dt, do.plots = 2)
fn <- optim.fn2(my.constraints2$par, dt, gen.fn = T)
sapply(seq(0,2,.1), function(x) fn('anc_1', x))

dt.all_models = dt.sims.p[, .N, .(v,c,a,d)]
setkey(dt.all_models, v, c, a, d)
my.gt = list(1,2,1,1)
dt.all_models[my.gt, fn := .(list(fn=fn))]
dt.all_models[my.gt, fn[[1]]]('v', .7)

p_gt_given_b_t_arcs <- function(my.b, my.t, my.gt) {
  dt.all_models[my.gt, fn[[1]]](my.b, my,t)
}
