##### BASE IMAGE #####
FROM ubuntu:26.04

##### METADATA #####
LABEL base.image="ubuntu:26.04"
LABEL version="2.0"
LABEL software="Sanity"
LABEL software.version="2.0"
LABEL software.description="Sanity"
LABEL software.website="https://github.com/jmbreda/Sanity"
LABEL software.documentation="https://github.com/jmbreda/Sanity"
LABEL software.license="GNU General Public License v3.0"
LABEL software.tags="Genomics, Transcriptomics"
LABEL maintainer="mikhail.pachkov@unibas.ch"
LABEL maintainer.organisation="Biozentrum, University of Basel"
LABEL maintainer.location="Spitalstrasse 41, CH-4056 Basel, Switzerland"
LABEL maintainer.lab="Erik van Nimwegen Lab"

##### VARIABLES #####
# Use variables for convenient updates/re-usability
ENV SOFTWARE_VERSION=2.0

RUN apt-get update \
  && apt-get install -y tzdata \
  && ln -fs /usr/share/zoneinfo/Europe/Berlin /etc/localtime \
  && dpkg-reconfigure --frontend noninteractive tzdata \
  && apt-get install --yes git make g++ libgomp1 zlib1g zlib1g-dev\
  && git clone https://github.com/jmbreda/Sanity.git \
  && cd Sanity \
  && cd src \
  && make \
  && cp ../bin/Sanity /usr/bin \
  && cd ../../ \
  && rm -rf Sanity \
  && apt-get remove --purge --yes git make g++ \
  && apt-get autoremove --purge --yes

# Sanity itself is PID 1, so the container lives exactly as long as the run and exits with the
# tool's own exit code. Arguments given to `docker run` after the image name are passed straight
# through to Sanity. With no arguments the CMD default prints Sanity's own help, which cannot
# drift from the actual set of options the way a hand-written usage block does.
ENTRYPOINT ["/usr/bin/Sanity"]
CMD ["--help"]

#### USAGE ####
# Build:
# docker build -t jmbreda/sanity:2.0 .
# Show the available options:
#  docker run --rm jmbreda/sanity:2.0
# Run Sanity on data in the current directory, writing results back to it:
#  docker run --rm -v "$PWD":/mnt -w /mnt jmbreda/sanity:2.0 -f [data.tsv] -d sanity_results -e 1 -n 1
#### ####
