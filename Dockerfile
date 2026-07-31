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
ENV SOFTWARE_VERSION 2.0

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

CMD ["/bin/sh", "-c", "cat <<'EOF'\n##### USAGE #####\n# Build:\n# docker build -t jmbreda/sanity:2.0 .\n# Run:\n# - Start image with mounting working directory to /mnt\n#   docker run -d --name sanity -v [path to working dir]:mnt -t jmbreda/sanity:2.0 > /dev/null\n# - Run Sanity:\n#  docker exec -w mnt -t sanity bash -c \"/usr/bin/Sanity -f [data.tsv] -d sanity_results -e1 -n 1\"\n# - Stop and remove container:\n#  docker stop sanity && docker rm sanity\n#### ####\nEOF"]

#### USAGE ####
# Build:
# docker build -t jmbreda/sanity:2.0 .
# Run:
# - Start image with mounting working directory to /mnt
#   docker run -d --name sanity -v [path to working dir]:mnt -t jmbreda/sanity:2.0 > /dev/null
# - Run Sanity:
#  docker exec -w mnt -t sanity bash -c "/usr/bin/Sanity -f [data.tsv] -d sanity_results -e1 -n 1"
# - Stop and remove container:
#  docker stop sanity && docker rm sanity
#### ####
