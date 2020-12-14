FROM nfcore/base:1.10.2
LABEL authors="Barry Digby" \
      description="MA5112 week1 tutorial"

WORKDIR ./
COPY week1.yml ./
RUN conda env create -f week1.yml && conda clean -a
ENV PATH /opt/conda/envs/week1/bin:$PATH
