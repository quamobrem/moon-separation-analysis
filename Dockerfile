FROM continuumio/miniconda3:main

# COPY environment.yml environment.yml
COPY orekit-data.zip orekit-data.zip

RUN conda config --system --set remote_connect_timeout_secs 30.0
RUN conda config --system --set remote_max_retries 10
RUN conda config --system --set remote_read_timeout_secs 300.0

RUN conda create --name separation-analysis-env python=3.13.1

RUN conda install -c conda-forge -n separation-analysis-env orekit=12.2
RUN conda install -n separation-analysis-env numpy=2.2.2
RUN conda install -n separation-analysis-env scipy=1.15.1
RUN conda install -n separation-analysis-env matplotlib=3.10.0


COPY src src
COPY main.py main.py

ENTRYPOINT [ "conda", "run", "-n",  "separation-analysis-env", "python", "-m", "main" ]
CMD [ "-h" ]
