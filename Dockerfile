# Instructions for Running
# 1. sudo docker build -t mutspot-env .
# 2. sudo docker run -it --rm -v "$(pwd)":/mutspot-analysis mutspot-env

FROM rocker/r-ver:3.6.3

WORKDIR /mutspot-analysis

# Download from archives
RUN sed -i 's/deb.debian.org/archive.debian.org/g' /etc/apt/sources.list && \
    sed -i 's|security.debian.org/debian-security|archive.debian.org/debian-security|g' /etc/apt/sources.list && \
    sed -i '/stretch-updates/d' /etc/apt/sources.list && \
    sed -i '/buster-updates/d' /etc/apt/sources.list

# Download libraries that R needs
RUN apt-get update && apt-get install -y \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    libxt-dev \
    libglu1-mesa-dev \
    libfreetype6-dev \
    libpng-dev \
    && rm -rf /var/lib/apt/lists/*

# Install everything via BiocManager to ensure version compatibility
RUN R -e "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/focal/2020-12-31')); \
    install.packages('BiocManager'); \
    BiocManager::install(version = '3.10', ask = FALSE, update = FALSE); \
    BiocManager::install(c( \
    'arrow', 'askpass', 'assertthat', 'backports', 'base64enc', 'BH', 'bit', 'bit64', 'bitops', \
    'blob', 'brew', 'brio', 'cachem', 'callr', 'cellranger', 'checkmate', 'cli', 'clipr', \
    'colorspace', 'commonmark', 'covr', 'cpp11', 'crayon', 'crosstalk', 'curl', 'data.table', \
    'DBI', 'desc', 'devtools', 'diffobj', 'digest', 'dplyr', 'DT', 'ellipsis', 'evaluate', \
    'fansi', 'farver', 'fastmap', 'fastshap', 'fontawesome', 'foreach', 'formatR', 'Formula', \
    'fs', 'futile.logger', 'futile.options', 'future', 'future.apply', 'generics', 'ggfittext', \
    'gggenes', 'ggplot2', 'ggrepel', 'gh', 'git2r', 'gitcreds', 'glmnet', 'globals', 'glue', \
    'gridtext', 'gtable', 'highr', 'htmltools', 'htmlwidgets', 'httr', 'ids', 'iml', 'ini', \
    'isoband', 'iterators', 'jpeg', 'jquerylib', 'jsonlite', 'knitr', 'labeling', 'lambda.r', \
    'later', 'lazyeval', 'lgr', 'lifecycle', 'listenv', 'litedown', 'magrittr', 'markdown', \
    'matrixStats', 'memoise', 'Metrics', 'mime', 'mirai', 'mlbench', 'mlr3', 'mlr3measures', \
    'mlr3misc', 'munsell', 'nanonext', 'openssl', 'palmerpenguins', 'paradox', 'parallelly', \
    'patchwork', 'pillar', 'pkgbuild', 'pkgconfig', 'pkgload', 'plyr', 'png', 'poibin', \
    'praise', 'prettyunits', 'processx', 'promises', 'PRROC', 'ps', 'pscl', 'purrr', 'R6', \
    'rcmdcheck', 'RColorBrewer', 'Rcpp', 'RcppArmadillo', 'RcppEigen', 'RCurl', 'rematch', \
    'rematch2', 'remotes', 'renv', 'rex', 'rlang', 'roxygen2', 'rprojroot', 'rstudioapi', \
    'rversions', 'scales', 'sessioninfo', 'shades', 'shape', 'shapviz', 'snow', 'stringi', \
    'stringr', 'sys', 'testthat', 'tibble', 'tidyselect', 'tinytex', 'usethis', 'utf8', \
    'uuid', 'vctrs', 'viridisLite', 'waldo', 'whisker', 'withr', 'xfun', 'xgboost', 'XML', \
    'xml2', 'xopen', 'yaml', \
    'Biobase', 'BiocGenerics', 'BiocParallel', 'BiocVersion', 'Biostrings', 'BSgenome', \
    'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38', 'DelayedArray', \
    'GenomeInfoDb', 'GenomeInfoDbData', 'GenomicAlignments', 'GenomicRanges', 'IRanges', \
    'Rhtslib', 'Rsamtools', 'rtracklayer', 'S4Vectors', 'SummarizedExperiment', 'XVector', 'zlibbioc'))"

# Install MutSpot from Github
RUN R -e "remotes::install_github('skandlab/MutSpot', subdir='MutSpot_Rpackage', upgrade='never')"

# Run R
CMD ["R"]
