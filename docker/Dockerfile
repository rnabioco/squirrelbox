# Install R version 4.0.0
FROM r-base:4.0.0

# Install Ubuntu packages
RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev/unstable \
    libxt-dev \
    libssl-dev

# Download and install ShinyServer (latest version)
RUN wget --no-verbose https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb

# Install R packages that are required
# TODO: add further package if you need!
RUN R -e "install.packages(c('shiny', 'remotes', 'BiocManager', 'shinydashboard', 'dplyr', 'tibble', 'purrr', 'readr', 'stringr', 'tidyr', 'feather', 'ggplot2', 'ggrepel', 'cowplot', 'DT', 'igraph', 'valr', 'plotly', 'crosstalk', 'shinyjs', 'shinythemes', 'shinycustomloader', 'shinyjqui', 'shinyWidgets', 'rintrojs', 'shinyBS', 'bsplus', 'BioCircos'), repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('transite')"
RUN installGithub.r jokergoo/ComplexHeatmap

# Copy configuration files into the Docker image
COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY /app /srv/shiny-server/

# Make the ShinyApp available at port 80
EXPOSE 80

# Copy further configuration files into the Docker image
COPY shiny-server.sh /usr/bin/shiny-server.sh

CMD ["/usr/bin/shiny-server.sh"]
