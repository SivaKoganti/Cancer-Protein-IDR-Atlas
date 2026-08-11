FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    build-essential python3 python3-pip wget git default-jre perl gzip \
    && rm -rf /var/lib/apt/lists/*

# Install Python deps
COPY requirements.txt /tmp/requirements.txt
RUN pip3 install --no-cache-dir -r /tmp/requirements.txt

# Add workflow scripts
WORKDIR /opt/pipeline
COPY . /opt/pipeline

ENV PATH="/opt/pipeline/scripts:${PATH}"

CMD ["snakemake", "--cores", "4"]
