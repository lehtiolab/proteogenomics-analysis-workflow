FROM nfcore/base
LABEL description="Docker image containing all requirements for IPAW pipeline"

COPY envs /envs/
RUN conda env create -f /envs/environment.yml && conda clean -a

RUN git clone https://github.com/yafeng/SpectrumAI /SpectrumAI
RUN cd /SpectrumAI && git pull && git reset --hard b8e7001807d834db633c30d265ef6e8361cdcb3c

ENV PATH /opt/conda/envs/ipaw-0.4/bin:$PATH
