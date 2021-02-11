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

    time Rscript R/estim_branchpoints_from_sims.R -gts GENOTYPE_FILE.tsv \
        --sims data/simulated_sfs.tsv \
        -nc 1 \
        --tag-labels SampleID \
        --tags my_sample \
        -sites all \
        --sim-method simple \
        --libs my_lib1 my_lib2 \
        -f-mh f_mh.yri \
        --n-qc1 0 \
        -table OUTPUT_FILE.tsv \
        --faunal-der-rate 0.0 \
        --num-em-iters 100 \
        --ll-surface --nbins 30 --ll-surface-thresh 5
        
        
    --deam-only
