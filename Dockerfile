FROM continuumio/miniconda3:main

COPY environment.yml environment.yml
COPY orekit-data.zip orekit-data.zip

RUN conda env create -f environment.yml

COPY src src
COPY main.py main.py

ENTRYPOINT [ "conda run -n separation-analysis-env python -m main" ]
CMD [ "-h" ]
