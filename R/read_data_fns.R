


####
## confirm that all of the correct columns are present
## set the f_mh column
## sample a mh from f_mh if mh is not present
## filter down only to required columns
## require 
## ISSUE  - this modifies dt.sed.og a little, but not fully, and then returns a new dt.
##          should probably not modify dt.sed.og, keeping it, in fact, OG.
check_and_add_required_cols <- function(dt.sed.og, f_mh.col, agCols, sample_mh_from_freq) {
  
  ## is it an issue that this is a string, and the column is usually T/F?
  if (!'deam53' %in% colnames(dt.sed.og)) dt.sed.og[, deam53 := 'unknown']
  if (!'freqs.FLAG' %in% colnames(dt.sed.og)) dt.sed.og[, freqs.FLAG := '.']
  
  req_columns <- c('sed_gt', agCols, f_mh.col, 'lib', 'deam53')
  if (sum(!req_columns %in% colnames(dt.sed.og)) > 0) {
    cat('Not all required columns are present:\n')
    cat('Required: ', req_columns, '\n')
    cat('Missing: ', req_columns[!req_columns %in% colnames(dt.sed.og)], '\n')
    stop()
  }
  
  ## setting requested f_mh column
  cat('Using', f_mh.col, 'as modern human frequency column.\n')
  setnames(dt.sed.og, f_mh.col, 'f_mh')
  
  ## if there's no mh column, then sample from f_mh
  ## anyway, we require complete.cases for f_mh
  dt.sed.og <- dt.sed.og[!is.na(dt.sed.og$f_mh)]
  if (!'mh' %in% colnames(dt.sed.og) || sample_mh_from_freq) {
    cat('Sampling mh [haploid] state from modern human frequencies\n')
    dt.sed.og$mh <- sapply(dt.sed.og[, f_mh], function(f_mh) sample(0:1, 1, prob = c(1-f_mh,f_mh)))
    # ggplot(dt.sed.og[, .(p = sum(mh)/.N, .N), f_mh], aes(x=f_mh, y=p, size=N)) + geom_point()
  }
  
  ## now remove unneeded columns
  ### ISSUE - this should be dynamic, and take into account agCols
  dt.sed <- dt.sed.og[, .(sed_gt, v_gt, c_gt, a_gt, d_gt, f_mh, mh, lib, freqs.FLAG, deam53)]
  dt.sed
}


basic_filtering_gts <- function(dt.sed, keep.libs, keep.libs.downsample, keep.libs.add_deam, sample_mh_from_freq, include_ti, merge_libs) {

    dt.sed.prefilter_counts <- dt.sed[, .(N.prefilter = .N), keyby=lib]
    
  cat('Filtering rows with NA in required columns:', sum(!complete.cases(dt.sed)), '/', dt.sed[, .N], '\n')
  dt.sed <- dt.sed[complete.cases(dt.sed)]
  
  if (!include_ti) {
    cat('Filtering transversions:', dt.sed[!freqs.FLAG %like% 'transi', .N], '/', dt.sed[, .N], 'sites removed\n')
    dt.sed <- dt.sed[!freqs.FLAG %like% 'transi']
  } else {
    cat('Keeping transversions! Should implement strand filtering!\n')
  }

    dt.sed.postfilter_counts <- dt.sed[, .(N.postfilter = .N), keyby=lib]

    cat('Read counts pre and post filter:\n')
    print(merge(dt.sed.prefilter_counts, dt.sed.postfilter_counts))

  if (!is.null(keep.libs)) {
    if (sum(!keep.libs %in% dt.sed[, unique(lib)]) > 0) {
      cat('Requested library is not in dataset:', keep.libs[!keep.libs %in% dt.sed[, unique(lib)]], '\n')
      q(save='no', status=1)
    }
    cat('Restricting to requested libraries: ', keep.libs, '\n')
    dt.sed <- dt.sed[lib %in% keep.libs]

    if (!is.null(keep.libs.downsample)) {
        cat('Downsampling requested libraries: ', keep.libs, '\n')
        cat('To proportions/counts: ', keep.libs.downsample, '\n')

        ## check that downsample counts are less than .N
        libs.counts <- sapply(keep.libs, function(x) dt.sed[lib == x, .N])
        cat('Pre-filter lib counts:', libs.counts, '\n')
        if (sum(libs.counts < keep.libs.downsample) > 0) {
            cat('Requested lib downsample numbers are larger than number of reads in lib:', keep.libs[libs.counts < keep.libs.downsample], '\n')
            cat(' - requested:', keep.libs.downsample[libs.counts < keep.libs.downsample], '\n')
            cat(' - actual:', libs.counts[libs.counts < keep.libs.downsample], '\n')
            q(save='no', status=1)
        }

        dt.sed <- foreach (my.idx = 1:length(keep.libs), .combine = rbind) %do% {
            my.lib <- keep.libs[my.idx]
            dt.tmp <- dt.sed[lib == my.lib]
            my.downsample <- keep.libs.downsample[my.idx]
            
            if (my.downsample > 1) return(dt.tmp[sample(.N, my.downsample)])
            ## if (my.downsample < 1) return(dt.tmp[sample(.N, .N * downsample)])
        }
        libs.counts <- sapply(keep.libs, function(x) dt.sed[lib == x, .N])
        cat('Post-filter lib counts:', libs.counts, '\n')
            
    }
    cat('Keeping: ', dt.sed[, .N], 'sites\n')

    if (!is.null(keep.libs.add_deam)) {
        cat('Adding "deamination" to requested libraries: ', keep.libs, '\n')
        cat(' - specifically to this proportion of non-deaminated reads: ', keep.libs.add_deam, '\n')
        libs.deam <- sapply(keep.libs, function(x) dt.sed[lib == x, sum(deam53 == T)/.N])
        cat('Deam rates before modification: ', libs.deam, '\n')
        dt.sed <- foreach (my.idx = 1:length(keep.libs), .combine = rbind) %do% {
            my.lib <- keep.libs[my.idx]
            dt.tmp <- dt.sed[lib == my.lib]
            my.deam <- keep.libs.add_deam[my.idx]

            if (my.deam <= 1) {
                dt.tmp[deam53 == F, deam53 := sample(c(T,F), .N, replace = T, prob = c(my.deam,1-my.deam))]
                return(dt.tmp)
            }
        }
        libs.deam <- sapply(keep.libs, function(x) dt.sed[lib == x, sum(deam53 == T)/.N])
        cat('Deam rates after modification: ', libs.deam, '\n')
    }
    
  } else {
    ## not really necessary, takes up a lot more space...
    dt.sed <- data.table(dt.sed)
    cat('Using all libraries\n')
  }
  
  if (merge_libs) {
    dt.sed[, lib := 'merged_libs']
    cat('Merging libraries\n')
  }
  dt.sed
}


subset_gt_states <- function(dt.sed, site_cats) {

  ### ISSUE - this requires that v_gt, c_gt, a_gt, d_gt are kept/present in sed.gt
  if ('poly_arc' %in% site_cats || 'all' %in% site_cats) {
    dt.sed.poly.full <- dt.sed[!(v_gt == c_gt & v_gt == a_gt & v_gt == d_gt)]
    cat('Polymorphic in archaic: ', dt.sed.poly.full[, .N], 'sites\n')
  } else if ('poly_neand' %in% site_cats) {
    dt.sed.poly.full <- dt.sed[!(v_gt == c_gt & v_gt == a_gt)]
    cat('Polymorphic in neands: ', dt.sed.poly.full[, .N], 'sites\n')
  } else {
    dt.sed.poly.full <- data.table()
  }
  
  ### ISSUE - this requires that v_gt, c_gt, a_gt, d_gt are kept/present in sed.gt
  if ('mh_seg_arc_fixed0' %in% site_cats || 'all' %in% site_cats) {
    dt.sed.mh.full <- dt.sed[v_gt == c_gt & v_gt == a_gt & v_gt == d_gt & v_gt == 0 & f_mh > 0]
    cat('Ancestral in archaics, seg in MH: ', dt.sed.mh.full[, .N], 'sites\n')
  } else {
    dt.sed.mh.full <- data.table()
  }

  ## ISSUE - this uses freqs.FLAG, which might not be present for most people [unless they merge w/ steffi's gts using my script]
  ## ISSUE - this also doesn't check if they are the hominin specific QC or pan_gor specific QC - this might cause some issue for really
  ## early branching hominins / denisovans have an elevated faunal contam proportion?
  if ('sed_qc_hominin' %in% site_cats || 'all' %in% site_cats) {
    dt.sed.qc.full <- dt.sed[freqs.FLAG == 'invar' & c_gt == 2]
    cat('Fixed derived in archaics and all [surveyed] MH: ', dt.sed.qc.full[, .N], 'sites\n')
  } else {
    dt.sed.qc.full <- data.table()
  }
  
  ### ISSUE - we end up with a dt.sed that is not sorted by position, which causes issues for the block sampling / jackknifing that I do later  
  rbind(dt.sed.poly.full, dt.sed.mh.full, dt.sed.qc.full)
}



add_sim_qc_sites <- function(dt.sed, n_qc0, n_qc1) {
  
  ## ISSUE - this assumes certain columns are used / present, not a flexible way of doing this at all.
  if (n_qc0 + n_qc1 == 0) return(dt.sed)
  
  ## simulate QC sites - only really useful for simulated data
  dt.sed.qc.full <- foreach(my.lib = dt.sed[, unique(lib)], .combine = rbind) %do% {
    ## just duplicate the first row of dt.sed the correct number of times
    dt.sed.qc.full <- dt.sed[lib == my.lib][rep(1, n_qc0 + n_qc1)]
    
    qc_fill_freq_or_hap <- c(rep(0,n_qc0), rep(1,n_qc1))
    qc_fill_gt <- c(rep(0,n_qc0), rep(2,n_qc1))

    ## ISSUE - this assumes certain columns are used / present, not a flexible way of doing this at all.
    ## should really drop this altogether, and require that the user simulate the sites they want [by adding an outgroup]
    dt.sed.qc.full[, f_mh := qc_fill_freq_or_hap]
    dt.sed.qc.full[, mh := qc_fill_freq_or_hap]
    dt.sed.qc.full[, v_gt := qc_fill_gt]
    dt.sed.qc.full[, c_gt := qc_fill_gt]
    dt.sed.qc.full[, a_gt := qc_fill_gt]
    dt.sed.qc.full[, d_gt := qc_fill_gt]
    dt.sed.qc.full[, sed_gt := qc_fill_freq_or_hap]
    dt.sed.qc.full
  }
  rbind(dt.sed.qc.full, dt.sed)
}

downsample_gt_sites <- function(dt.sed, downsample) {

  if (downsample == 0) return(dt.sed)
  
  ### ISSUE - should check that downsample is less than .N
  if (downsample > 1) return(dt.sed[sample(.N, downsample)])
  if (downsample < 1) return(dt.sed[sample(.N, .N * downsample)])

}


add_rg <- function(dt.sed, rg_props) {
  ####################
  ## create read groups
  
  ## make sure that proportions sum to 1
  rg_props <- rg_props / sum(rg_props)
  
  ## ISSUE - this doesn't work correctly if there are multple libraries and multiple read groups,
  ## because you end up with more read groups than contam proportions, etc
  dt.sed[, rg := paste0(lib, '_rg_', sample(1:length(rg_props), .N, replace = T, prob = rg_props), '_', deam53)]
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


read_and_process_genos <- function(gts_file, f_mh.col = 'f_mh', agCols = c('v_gt', 'c_gt', 'a_gt', 'd_gt'),
                                   keep.libs = NULL, keep.libs.downsample = NULL, keep.libs.add_deam = NULL,
                                   sample_mh_from_freq = T, include_ti = F, merge_libs = F, 
                                   site.cats = 'all', n_qc0 = 0, n_qc1 = 0, downsample = 0, rg_props = 1,
                                   block_bootstrap = 0) {
  cat('Reading genotypes..', gts_file, '\n')
  dt.sed.og <- fread(gts_file)
  cat('  Read:', dt.sed.og[, .N], 'sites\n')
  cat('  Libraries in file:', dt.sed.og[, unique(lib)], '\n')
  
  dt.sed <- check_and_add_required_cols(dt.sed.og, f_mh.col = f_mh.col, agCols = agCols, sample_mh_from_freq = sample_mh_from_freq)
  
  dt.sed <- basic_filtering_gts(dt.sed, 
                                keep.libs = keep.libs, 
                                keep.libs.downsample = keep.libs.downsample, 
                                keep.libs.add_deam = keep.libs.add_deam, 
                                sample_mh_from_freq = sample_mh_from_freq,
                                include_ti = include_ti,
                                merge_libs = merge_libs)
  
  add_gt_categories(dt.sed, agCols = agCols)
  
  dt.sed.analysis <- subset_gt_states(dt.sed, site_cats = site.cats)
  dt.sed.analysis <- add_sim_qc_sites(dt.sed.analysis, n_qc0 = n_qc0, n_qc1 = n_qc1)
  
  dt.sed.analysis <- downsample_gt_sites(dt.sed.analysis, downsample = downsample)
  
  if (block_bootstrap > 1) dt.sed.analysis <- bootstrap_gt_data(dt.sed.analysis, blocks = block_bootstrap)
  
  add_rg(dt.sed.analysis, rg_props = rg_props)
  
}












