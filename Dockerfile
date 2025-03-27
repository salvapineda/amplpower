# Use any image as base image
FROM ubuntu:22.04
# Install curl in order to download the modules necessary
RUN apt-get update && apt-get install -y curl python3 python3-pip python3-venv libc6-dev

# Install full bundle (includes AMPL and all solvers)
RUN cd /opt/ && curl -OL https://ampl.com/dl/amplce/ampl.linux64.tgz && \
    tar oxzvf ampl.linux64.tgz && rm ampl.linux64.tgz

# Add installation directory to the environment variable PATH
ENV PATH="/opt/ampl.linux-intel64/:${PATH}"
