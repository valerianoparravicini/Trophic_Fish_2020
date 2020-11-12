
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Global gut content data synthesis and phylogeny delineate reef fish trophic guilds

This repository contains code and data needed to reproduce the figures
and result of the manuscript:

FULL REFERENCE AND DOI

## Content

The directory contains:  
\- [:page\_facing\_up: DESCRIPTION](/DESCRIPTION):  
\- [:page\_facing\_up: make.R](/make.R): Script to reproduce all parts
of the project.  
\- [:file\_folder: data](/data): Folder containing all data.  
\- [:file\_folder: scripts](/scripts): Folder containing scripts to
perform multiple parts of the analysis to produce results and figures.  
\- [:file\_folder: output](/output): Folder containing all results and
figures.

### scripts

  - **01\_guilds\_compare\_experts.R**: We conducted a systematic review
    in the literature to collect information about how different authors
    have classified reef fish species into different trophic guilds. We
    could obtain 30 distinct classifications and this information is
    contained in the file “data/original\_experts\_classification.csv”.
    These data were then made comparable (see
    “data/converted\_experts\_classification.csv”. The script compares
    pairs of experts and quantify the agreement/disagreement among
    them.  
  - **02\_global\_network\_guilds\_site\_specific.R**:  
    This script defines the diet for each combination of reef fish
    species and site. Then we used the method introduced by Becket
    (2016) to identify trophic guilds according to a maximization of the
    modularity of the weighted bipartite network. Once we identified the
    trophic guilds for each species x site combination we realized that
    in the vast majority of cases the same species was classified in the
    same way regardless the location of sampling.

Beckett, S.J. 2016 Improved community detection in weighted bipartite
networks. Royal Society open science 3, 140536.

  - **03\_global\_network\_guilds.R**:  
    Once we realized that, at our resolution, the trophic guilds for the
    same species in different locations were exactly the same, we
    decided to perform the same analysis at global level (i.e. combining
    ditar data for the same species coming from different sites) to
    define the final trophic guilds. This is what this script does.

  - **04\_trophic\_interactions\_calibration\_dataset.R**  
    Once the guilds were defined we tested whether it was possible to
    directly predict the probability of interaction between reef fishes
    and the prey categories used. We did that using machine learning to
    predict trophic interactions according to the phylogenetic position
    of reef fishes and their body size. This script perform a series of
    preliminary analyses (i.e. definition of the most likely
    phylogenetic tree for reef fishes, extraction of phylogenetic
    eigenvector maps for reef fishes) that are needed to obtain the
    calibration dataset.

  - **05\_h20\_ensemble\_model\_calibration.R**  
    This script sets up a h2o cluster (<https://www.h2o.ai>) to predict
    the probability of trophic interaction between reef fish species and
    our prey categories. Then it uses the calibration dataset to perform
    an ensemble model composed of: boosted regression trees, extreme
    gradient boosting, random forest and a glm used as a super-learner.

  - **06\_h20\_model\_check\_and\_extrapolation.R** This script contains
    the functions needed to check the model performance (e.g. AUC, TSS)
    and to extrapolate the likely trophic interaction for all reef fish
    species.

## Working environment

``` r
sessionInfo()
#> R version 3.6.3 (2020-02-29)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 18.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
#> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
#> 
#> locale:
#>  [1] LC_CTYPE=fr_FR.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=fr_FR.UTF-8        LC_COLLATE=fr_FR.UTF-8    
#>  [5] LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=fr_FR.UTF-8   
#>  [7] LC_PAPER=fr_FR.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] compiler_3.6.3  magrittr_1.5    tools_3.6.3     htmltools_0.4.0
#>  [5] yaml_2.2.0      Rcpp_1.0.5      stringi_1.4.3   rmarkdown_2.1  
#>  [9] knitr_1.25      stringr_1.4.0   xfun_0.10       digest_0.6.22  
#> [13] rlang_0.4.7     evaluate_0.14
```
