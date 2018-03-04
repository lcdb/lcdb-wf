FROM bioconda/bioconda-utils-build-env
COPY requirements.txt requirements.txt
COPY .circleci/setup.sh /tmp/setup.sh
RUN bash /tmp/setup.sh
ENV PATH="/miniconda/bin:${PATH}"
