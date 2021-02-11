# branch_estim
This R script implements an EM algorithm for estimating the branching point of a sample _X_ from a larger phylogeny, given a fixed demographic model learned from high-quality samples. In this case, the phylogeny is composed of four high-coverage archaic genomes: three Neandertals (Altai Neandertal, Vindija 33.19, Chagyrskaya 8) and one Denisovan (Denisova 3). For a demographic model with tips _T_, one can calculate the probability of observing a derived allele in _X_, given that _X_ branches from a specific part of the model, and given the tip genotypes. These probabilities are then informative for the branching point of the sample _X_. An example using Neandertals and Denisovans is given in the figure below, where _P(X = 1)_ is the probability of observing a derived allele in _X_.


given joint frequency spectra for _X_ and _T_ at arbitrary branching points of _X_ from the tree, estimate the most likely branching point of _X_ from the model. I have used this for estimating the split time of low-quality archaic hominin samples. for branch-point modeling of a low-quality sample, . Primarily used for sediment libraries.

## Input Data



### Genotype file

The genotypes are read from a tab-delimited file, with one line per read. This represents a typical file:

| chrom | pos       | freqs.FLAG | v_gt | c_gt | a_gt | d_gt | sed_gt | f_mh.yri   | deam53 | lib    |
|-------|-----------|------------|------|------|------|------|--------|------------|--------|--------|
| 11    | 115953451 | invar      | 2    | 2    | 2    | 2    | 1      | 1          | FALSE  | A16112 |
| 1     | 52606336  | transi     | 2    | 2    | 2    | 0    | 1      | 0          | TRUE   | A16112 |
| 2     | 125524971 | .          | 0    | 0    | 0    | 0    | 0      | 0.2037037  | FALSE  | A16112 |
| 8     | 133156471 | .          | 1    | 0    | 0    | 0    | 0      | 0          | FALSE  | A16112 |
| 16    | 50532171  | transi     | 0    | 0    | 0    | 0    | 0      | 0.09259259 | FALSE  | A16112 |
| 15    | 65654559  | invar      | 2    | 2    | 2    | 2    | 1      | 1          | FALSE  | A16112 |

The following columns are required:

| COLUMN     | DESCRIPTION                                                                                                                                                                                                                   |
|------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| chrom      | Chromosome                                                                                                                                                                                                                    |
| pos        | Positionm                                                                                                                                                                                                                     |
| freqs.FLAG | A flag indicating certain categories of site. 'invar' indicates sites   that are derived in all known hominins, and is used to calculate faunal   contamination. 'transi' indicates transition sites. '.' is everything else. |
| v_gt       | Derived diploid genotype in Vindija 33.19 (0/1/2).                                                                                                                                                                            |
| c_gt       | Derived diploid genotype in Chagyrskaya 8 (0/1/2).                                                                                                                                                                            |
| a_gt       | Derived diploid genotype in the Altai Neandertal (0/1/2).                                                                                                                                                                     |
| d_gt       | Derived diploid genotype in Denisova 5 (the high coverage Denisovan)   (0/1/2).                                                                                                                                               |
| sed_gt     | Derived haploid genotype in the sample of interest (0/1).                                                                                                                                                                     |
| f_mh.yri   | Derived allele frequency in a modern human population - here, YRI.                                                                                                                                                            |
| deam53     | Is this read deaminated (TRUE/FALSE)                                                                                                                                                                                          |
| lib        | Library ID                                                                                                                                                                                                                    |

## Example commands

    sims=data/simulated_sfs.txt
    genos=data/all_simple_gts.A16112.short3b.tsv.gz

    time Rscript R/estim_branchpoints_from_sims.R \
        -gts $genos  \
        --sims $sims \
        -nc 1 \
        --tag-labels mylib \
        --tags A16112 \
        -sites all \
        --sim-method simple \
        --libs A16112 \
        -f-mh f_mh.yri \
        --n-qc1 0 \
        -table em_output.A16112.tsv \
        --faunal-der-rate 0.0 \
        --num-em-iters 100 \
        --ll-surface \
        --nsteps 0

This command creates the file `em_output.A16112.tsv`:

| rg                | mh_contam  | faunal_prop | nsnps | branchtime | branch | man.max.ll | man.max.ll.last | n.iter | my.t.idx | max.ll     | step.x | mylib  |
|-------------------|------------|-------------|-------|------------|--------|------------|-----------------|--------|----------|------------|--------|--------|
| A16112_rg_1_FALSE | 0.01930535 | 4.00E-10    | 1616  | 0.82130775 | v      | -261.27706 | -261.27706      | 13     | 0        | -261.27706 | 0      | A16112 |
| A16112_rg_1_TRUE  | 6.37E-10   | 9.07E-11    | 886   | 0.82130775 | v      | -261.27706 | -261.27706      | 13     | 0        | -261.27706 | 0      | A16112 |

With the following columns:

| COLUMN          | DESCRIPTION                                                                                                                                                                |
|-----------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| rg              | Read Group - this is a group of reads which   are analyzed together. Typically, this is all deaminated (TRUE) or   non-deaminated (FALSE) reads in one or more libraries.  |
| mh_contam       | Maximum likelihood estimate (MLE) of modern human contamination in this   read group.                                                                                      |
| faunal_prop     | MLE of faunal contamination in this read group.                                                                                                                            |
| nsnps           | Number of SNPs in this read group.                                                                                                                                         |
| branchtime      | MLE branchtime - this is estimated on all read groups together   (typically, deam and non-deam reads of one or more libraries).                                            |
| branch          | MLE branch.                                                                                                                                                                |
| man.max.ll      | Log-Likelihood of MLE solutions.                                                                                                                                           |
| man.max.ll.last | Next-to-last log-likelihood (can be used to identify the change in LL in   the last step).                                                                                 |
| n.iter          | Number of EM iterations performed.                                                                                                                                         |
| my.t.idx        | When performing a likelihood surface calculation, which gridpoint is   this?                                                                                               |
| max.ll          | When performing a likelihood surface calculation, the global maximum LL.                                                                                                   |
| step.x          | When performing a likelihood surface calculation, the step in the grid   search.                                                                                           |
| mylib           | Library ID (from column 'lib' in genotypes file).                                                                                                                          |
