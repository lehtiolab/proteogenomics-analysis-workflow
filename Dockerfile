FROM nfcore/base
LABEL description="Docker image containing all requirements for IPAW pipeline"

COPY environment.yml /
COPY tools /tools/

RUN conda env create -f /environment.yml && conda clean -a
RUN conda env create -f /tools/openms/environment.yml && conda clean -a

RUN git clone https://github.com/yafeng/SpectrumAI /SpectrumAI
RUN cd /SpectrumAI && git pull && git reset --hard d9fc290cd76a5ec09aa17c03a380ad09cbce2387

ENV PATH /opt/conda/envs/ipaw-1.0/bin:$PATH
