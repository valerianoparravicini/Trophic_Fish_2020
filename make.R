# Executing this file will reproduce all analyses and produce results and figures for the manuscript 
# "Global gut content data synthesis and phylogeny delineate reef fish trophic guilds"

# Install all needed packages
installed.packages("devtools")
devtools::install_deps()

# Load needed packages and functions
devtools::load_all()

# Run scripts one by one

# CAUTION: Rerunning these scripts requires a supercomputer 
# with at least 50 processors and may take about a week

scripts <- list.files("scripts/")

lapply(scripts, source)


