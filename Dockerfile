FROM docker.io/rocker/tidyverse:4.5.1

LABEL maintainer="you@example.com"

# Make ALL tar extractions rootless-friendly:
ENV TAR_OPTIONS="--no-same-owner --no-same-permissions"
ENV RUNROOTLESS=false

# Shiny server
RUN /rocker_scripts/install_shiny_server.sh 1.5.23.1030

# 1. System dependencies
##############################################
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    wget \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libgsl-dev \
    python3 \
    zlib1g-dev \
    git \
    openssh-client \
    curl \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# Ensure the RStudio user exists and has a password
RUN set -eux; \
    if ! getent passwd rstudio >/dev/null; then \
      useradd -m -s /bin/bash rstudio; \
    fi; \
    echo "rstudio:rstudio" | chpasswd; \
    mkdir -p /home/rstudio; \
    chown -R rstudio:rstudio /home/rstudio

# Make RStudio start in the project directory automatically (rstudio login)
RUN echo 'setwd("/srv/shiny-server/ERASOR")' >> /home/rstudio/.Rprofile && \
    chown rstudio:rstudio /home/rstudio/.Rprofile

# Configure git (system-wide defaults)
##############################################
RUN git config --system user.name "rstudio" && \
    git config --system user.email "rstudio@example.com"

# ViennaRNA 2.7.0 installation
##############################################
RUN wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_7_x/ViennaRNA-2.7.0.tar.gz && \
    tar --no-same-owner --no-same-permissions -xzf ViennaRNA-2.7.0.tar.gz && \
    cd ViennaRNA-2.7.0 && \
    ./configure && \
    make -j"$(nproc)" && \
    make install && \
    ldconfig && \
    cd .. && \
    rm -rf ViennaRNA-2.7.0 ViennaRNA-2.7.0.tar.gz

# Install R packages (CRAN + Bioconductor)
##############################################
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" \
 && R -e "BiocManager::install(c( \
        'GenomicFeatures', \
        'AnnotationDbi', \
        'BSgenome.Hsapiens.NCBI.GRCh38', \
        'biomaRt', \
        'Biostrings', \
        'txdbmaker' \
     ), ask = FALSE)" \
 && R -e "install.packages(c( \
        'shinythemes', 'shiny', 'shinyjs', 'tidyverse', 'cluster', \
        'rlang', 'dplyr', 'DT', 'shinyBS', 'shinydashboard', 'openxlsx', \
        'future', 'future.apply' \
     ), repos='https://cloud.r-project.org', dependencies=TRUE)"

# Make necessary directories
##############################################
RUN mkdir -p /srv/shiny-server/ERASOR && \
    chown -R shiny:shiny /srv/shiny-server && \
    chmod -R 775 /srv/shiny-server/ERASOR

RUN mkdir -p /home/rstudio/.ssh && \
    chown -R rstudio:rstudio /home/rstudio/.ssh && \
    chmod 700 /home/rstudio/.ssh

# Ensure Shiny has its log/state directories
##############################################
RUN mkdir -p /var/log/shiny-server /var/lib/shiny-server && \
    chown -R shiny:shiny /var/log/shiny-server /var/lib/shiny-server

# Copy the whole project into the Shiny Server root
##############################################
WORKDIR /srv/shiny-server/ERASOR
COPY . .

# Ownership/permissions: Shiny needs full access; RStudio needs at least read access.
# Keep Shiny as owner, but ensure others can read/enter dirs.
RUN chown -R shiny:shiny /srv/shiny-server/ERASOR && \
    chmod -R a+rX /srv/shiny-server/ERASOR

# Run Shiny apps as the 'shiny' system user
##############################################
RUN sed -i 's/^run_as.*/run_as shiny;/' /etc/shiny-server/shiny-server.conf || \
    sed -i '1irun_as shiny;' /etc/shiny-server/shiny-server.conf

# Download and create txdb database
##############################################
RUN mkdir -p /opt/ERASOR
RUN R -e "\
  library(txdbmaker); \
  gtf <- 'Homo_sapiens.GRCh38.115.gtf.gz'; \
  gtf_url <- paste0('https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/', gtf); \
  download.file(gtf_url, gtf, mode = 'wb'); \
  txdb <- makeTxDbFromGFF(gtf); \
  saveDb(txdb, '/opt/ERASOR/txdb_hsa_biomart.db'); \
  unlink(gtf); \
"

# Make RStudio start in the project directory automatically (rstudio login mode)
##############################################
RUN echo 'setwd("/srv/shiny-server/ASOstool-v2")' >> /home/rstudio/.Rprofile && \
    chown rstudio:rstudio /home/rstudio/.Rprofile

# Expose ports
##############################################
EXPOSE 3838 8787

# Default credentials for RStudio (rstudio user)
##############################################
ENV PASSWORD=rstudio

# Start services (RStudio + Shiny)
##############################################
CMD ["/init"]
