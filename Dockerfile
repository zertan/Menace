FROM alpine:3.4

MAINTAINER CognitiveScale <devops@cognitivescale.com>

RUN apk --update  --repository http://dl-4.alpinelinux.org/alpine/edge/community add \
    nano \
    bash \
    git \
	  gcc \
	  musl-dev \
	  zlib \
    zlib-dev \
	  linux-headers \
    parallel \
    openssh-client \
    curl \
    ca-certificates \
    bzip2 \
    unzip \
    sudo \
    libstdc++ \
    glib \
    libxext \
    libxrender \
    tini \
    && curl -L "https://github.com/andyshinn/alpine-pkg-glibc/releases/download/2.23-r1/glibc-2.23-r1.apk" -o /tmp/glibc.apk \
    && curl -L "https://github.com/andyshinn/alpine-pkg-glibc/releases/download/2.23-r1/glibc-bin-2.23-r1.apk" -o /tmp/glibc-bin.apk \
    && curl -L "https://github.com/andyshinn/alpine-pkg-glibc/releases/download/2.23-r1/glibc-i18n-2.23-r1.apk" -o /tmp/glibc-i18n.apk \
    && apk add --allow-untrusted /tmp/glibc*.apk \
    && /usr/glibc-compat/sbin/ldconfig /lib /usr/glibc-compat/lib \
    && /usr/glibc-compat/bin/localedef -i en_US -f UTF-8 en_US.UTF-8 \
    && rm -rf /tmp/glibc*apk /var/cache/apk/*

SHELL ["/bin/bash", "-c"]

# Configure environment
ENV CONDA_DIR=/opt/conda MENACE_DIR=/opt/PTR-Pipeline BITSEQ_DIR=/opt/BitSeq CONDA_VER=4.2.12
ENV PATH=$CONDA_DIR/bin:$MENACE_DIR:$PATH LANG=C SHELL=/bin/bash

# Install conda
RUN mkdir -p $CONDA_DIR && \
    echo export PATH=$CONDA_DIR/bin:'$PATH' > /etc/profile.d/conda.sh && \
    curl https://repo.continuum.io/miniconda/Miniconda2-${CONDA_VER}-Linux-x86_64.sh -o mconda.sh && \
    /bin/bash mconda.sh -f -b -p $CONDA_DIR && \
    rm mconda.sh && \
    $CONDA_DIR/bin/conda install --yes conda==${CONDA_VER}

# install packages
RUN conda config --add channels conda-forge \
 && conda config --add channels r \
 && conda config --add channels bioconda \
 && conda install blas \
 && conda install samtools \
 && conda install bamtools \
 && conda install bowtie2 \
 && conda clean -a --yes \
 && pip install --no-cache-dir --no-cache menace

# bitseq
#RUN mkdir -p $BITSEQ_DIR && \
#    cd /opt && \
#    git clone https://github.com/BitSeq/BitSeq.git \
# && cd BitSeq
# && make

# menace pipeline
#RUN mkdir -p $MENACE_DIR && \
#    cd /opt && \
#    git clone https://github.com/zertan/PTR-Pipeline.git \
# && pip install --no-cache-dir --no-cache -r PTR-Pipeline/requirements.txt
