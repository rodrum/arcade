# Use an official Ubuntu 20.04 base image
FROM mambaorg/micromamba:1-focal

USER root

ARG DEBIAN_FRONTEND=noninteractive

# Update and install system packages
RUN apt-get update && \
    apt-get install -y gfortran-8 libfftw3-dev python3-dev curl git tar wget build-essential cmake time nano && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Get ECcodes
WORKDIR /

RUN wget -q https://confluence.ecmwf.int/download/attachments/45757960/eccodes-2.30.0-Source.tar.gz \
    && tar -xzf eccodes-2.30.0-Source.tar.gz \
    && mkdir build

WORKDIR /build

RUN cmake ../eccodes-2.30.0-Source -DENABLE_AEC=OFF -DENABLE_FORTRAN=OFF \
    && make -j4 \
    && ctest -j4 \
    && make install

WORKDIR /

RUN rm eccodes-2.30.0-Source.tar.gz \
    && rm -rf eccodes-2.30.0-Source \
    && rm -rf build

# Download git repos
WORKDIR /app/arcade

COPY . .

WORKDIR /app/arcade/repos

# Get NRLMSIS2.0
ADD https://map.nrl.navy.mil/map/pub/nrl/NRLMSIS/NRLMSIS2.0/NRLMSIS2.0.tar.gz .

RUN mkdir NRLMSIS2.0 && \
    tar -xzf NRLMSIS2.0.tar.gz -C./NRLMSIS2.0 && \
    rm NRLMSIS2.0.tar.gz

# Get HMW14
ADD https://map.nrl.navy.mil/map/pub/nrl/HWM/HWM14/HWM14_ess224-sup-0002-supinfo.tgz .
   
RUN tar -xf HWM14_ess224-sup-0002-supinfo.tgz &&\
    rm HWM14_ess224-sup-0002-supinfo.tgz ._HWM14
    
## Modify makefile for gfortran-8

RUN sed -i "s|	ifort checkhwm14.f90 hwm14.f90 -o check.opt.exe|#	ifort checkhwm14.f90 hwm14.f90 -o check.opt.exe|g" HWM14/makefile && \
    sed -i "s|#	gfortran checkhwm14.f90 hwm14.f90 -o check.gcc.exe|	gfortran-8 checkhwm14.f90 hwm14.f90 -o check.gcc.exe|g" HWM14/makefile && \
    sed -i "s|	./check.opt.exe > Check/result.txt|	./check.gcc.exe > Check/result.txt|g" HWM14/makefile && \
    sed -i "s|	diff Check/result.txt Check/ifort.opt.txt |	diff Check/result.txt Check/gfortran.txt |g" HWM14/makefile
    
# Get infraGA
RUN git clone -b master https://github.com/LANL-Seismoacoustics/infraGA.git && \
    mv infraGA infraGA-master
    
# Initialize the repos

WORKDIR /app/arcade/repos/NRLMSIS2.0

RUN gfortran-8 -O3 -cpp -o msis2.0_test.exe alt2gph.F90 msis_constants.F90 msis_init.F90 msis_gfn.F90 msis_tfn.F90 msis_dfn.F90 msis_calc.F90 msis_gtd8d.F90 msis2.0_test.F90

WORKDIR /app/arcade/repos/HWM14

RUN make

WORKDIR /app/arcade/repos/infraGA-master

RUN make

WORKDIR /app/arcade

RUN make all-prep

# Set up python env with mamba
USER $MAMBA_USER

RUN micromamba install --yes --name base -f environment_cross_platform.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN micromamba install --yes --name base pip

RUN python3 -m pip install cdsapi

USER root

RUN printf "url: https://cds.climate.copernicus.eu/api/v2\nkey: {uid}:{api-key}\n" > /root/.cdsapirc

# Set the entrypoint
ENTRYPOINT ["/bin/bash"]
