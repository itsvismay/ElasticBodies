FROM ubuntu:16.04
MAINTAINER zack441@mac.com

##ADD ubuntu-xenial-core-cloudimg-amd64-root.tar.gz /

# a few minor docker-specific tweaks
# see https://github.com/docker/docker/blob/9a9fc01af8fb5d98b8eec0740716226fadb3735c/contrib/mkimage/debootstrap
##RUN set -xe \
	\
# https://github.com/docker/docker/blob/9a9fc01af8fb5d98b8eec0740716226fadb3735c/contrib/mkimage/debootstrap#L40-L48
##	&& echo '#!/bin/sh' > /usr/sbin/policy-rc.d \
##	&& echo 'exit 101' >> /usr/sbin/policy-rc.d \
##	&& chmod +x /usr/sbin/policy-rc.d \
	\
# https://github.com/docker/docker/blob/9a9fc01af8fb5d98b8eec0740716226fadb3735c/contrib/mkimage/debootstrap#L54-L56
##	&& dpkg-divert --local --rename --add /sbin/initctl \
##	&& cp -a /usr/sbin/policy-rc.d /sbin/initctl \
##	&& sed -i 's/^exit.*/exit 0/' /sbin/initctl \
##	\
# https://github.com/docker/docker/blob/9a9fc01af8fb5d98b8eec0740716226fadb3735c/contrib/mkimage/debootstrap#L71-L78
##	&& echo 'force-unsafe-io' > /etc/dpkg/dpkg.cfg.d/docker-apt-speedup \
	\
# https://github.com/docker/docker/blob/9a9fc01af8fb5d98b8eec0740716226fadb3735c/contrib/mkimage/debootstrap#L85-L105
##	&& echo 'DPkg::Post-Invoke { "rm -f /var/cache/apt/archives/*.deb /var/cache/apt/archives/partial/*.deb /var/cache/apt/*.bin || true"; };' > /etc/apt/apt.conf.d/docker-clean \
##	&& echo 'APT::Update::Post-Invoke { "rm -f /var/cache/apt/archives/*.deb /var/cache/apt/archives/partial/*.deb /var/cache/apt/*.bin || true"; };' >> /etc/apt/apt.conf.d/docker-clean \
##	&& echo 'Dir::Cache::pkgcache ""; Dir::Cache::srcpkgcache "";' >> /etc/apt/apt.conf.d/docker-clean \
##	\
# https://github.com/docker/docker/blob/9a9fc01af8fb5d98b8eec0740716226fadb3735c/contrib/mkimage/debootstrap#L109-L115
##	&& echo 'Acquire::Languages "none";' > /etc/apt/apt.conf.d/docker-no-languages \
##	\
# https://github.com/docker/docker/blob/9a9fc01af8fb5d98b8eec0740716226fadb3735c/contrib/mkimage/debootstrap#L118-L130
##	&& echo 'Acquire::GzipIndexes "true"; Acquire::CompressionTypes::Order:: "gz";' > /etc/apt/apt.conf.d/docker-gzip-indexes \
##	\
# https://github.com/docker/docker/blob/9a9fc01af8fb5d98b8eec0740716226fadb3735c/contrib/mkimage/debootstrap#L134-L151
##	&& echo 'Apt::AutoRemove::SuggestsImportant "false";' > /etc/apt/apt.conf.d/docker-autoremove-suggests

# delete all the apt list files since they're big and get stale quickly
##RUN rm -rf /var/lib/apt/lists/*
# this forces "apt-get update" in dependent images, which is also good

# enable the universe
##RUN sed -i 's/^#\s*\(deb.*universe\)$/\1/g' /etc/apt/sources.list

# make systemd-detect-virt return "docker"
# See: https://github.com/systemd/systemd/blob/aa0c34279ee40bce2f9681b496922dedbadfca19/src/basic/virt.c#L434
##RUN mkdir -p /run/systemd && echo 'docker' > /run/systemd/container

# overwrite this with 'CMD []' in a dependent Dockerfile
# CMD ["/bin/bash"]

# required libraries
RUN echo 'apt-get update'
RUN apt-get update
RUN echo 'install apt-utils'
RUN apt-get install -y apt-utils
RUN echo 'install aptitude'
RUN apt-get install -y aptitude
RUN echo 'install vim'
RUN apt-get install -y vim
RUN echo 'install git'
RUN apt-get install -y git
RUN echo 'install build-essential'
RUN apt-get install -y build-essential
RUN echo 'install libgl1-mesa-dev'
RUN apt-get install -y libgl1-mesa-dev
RUN echo 'install freeglut3-dev'
RUN apt-get install -y freeglut3-dev
RUN echo 'install cmake'
RUN apt-get install -y cmake
RUN echo 'install libblas-dev'
RUN apt-get install -y libblas-dev
RUN echo 'install libboost-dev'
RUN apt-get install -y libboost-dev
RUN echo 'install libmpfr4 and libmpfr4-dbg'
RUN apt-get install -y libmpfr4 libmpfr4-dbg
RUN echo 'install libopenlbas-dev'
RUN apt-get install -y libopenblas-dev
RUN echo 'install libpng12-dec'
RUN apt-get install -y libpng12-dev
RUN echo 'install python-matplotlib'
RUN apt-get install -y python-matplotlib
RUN echo 'install python-scipy'
RUN apt-get install -y python-scipy
RUN echo 'install python3-dev'
RUN apt-get install python3-dev
RUN echo 'install swig'
RUN apt-get install -y swig
RUN echo 'install zlib1g-dev'
RUN apt-get install -y zlib1g-dev
RUN echo 'install libxcursor-dev'
RUN apt-get install -y libxcursor-dev
RUN echo 'install libxi-dev'
RUN apt-get install -y libxi-dev
RUN echo 'install libxinerama-dev'
RUN apt-get install -y libxinerama-dev
RUN echo 'install libxrandr-dev'
RUN apt-get install -y libxrandr-dev
RUN echo 'install libcgal-dev'
RUN apt-get install -y libcgal-dev
RUN echo 'purge libcgal-dev'
RUN apt-get purge -y libcgal-dev
RUN echo 'install libeigen3-dev'
RUN apt-get install -y libeigen3-dev
RUN echo 'install openscad and openscad-dbg'
RUN apt-get install -y openscad openscad-dbg
RUN echo 'install slic3r'
RUN apt-get install -y slic3r
RUN echo 'install meshlab'
RUN apt-get install -y meshlab
RUN echo 'install clang'
RUN apt-get install -y clang
RUN echo 'install libtool'
RUN apt-get install -y libtool
RUN echo 'install libtbb-dev'
RUN apt-get install -y libtbb-dev
RUN echo 'install libqt4-dev'
RUN apt-get install -y libqt4-dev
RUN echo 'install wget'
RUN apt-get install -y wget
#RUN echo 'install qt4-dev-tool'
#RUN apt-get install -y qt4-dev-tool
RUN echo 'INSTALLED EVERYTHING'
#RUN echo 'install libcgal-qt5-dev'
#RUN apt install libcgal-qt5-dev
RUN apt install libgmp3-dev
RUN apt install libmpfr-dev


# git clones
RUN git clone --recursive https://github.com/qnzhou/PyMesh.git
RUN git clone https://github.com/chokkan/liblbfgs.git
RUN git clone --recursive https://github.com/libigl/libigl.git
RUN git clone https://github.com/jluttine/suitesparse.git
RUN git clone https://github.com/itsvismay/ElasticBodies.git
RUN git clone --recursive https://github.com/CGAL/cgal.git

RUN ls

# required downloads - use wget
# Anaconda
RUN wget https://repo.continuum.io/archive/Anaconda2-4.2.0-Linux-x86_64.sh
RUN bash Anaconda2-4.2.0-Linux-x86_64.sh -b
# CGAL
#RUN wget https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.9/CGAL-4.9.tar.xz
#RUN tar xf CGAL-4.9.tar.xz

RUN ls

# Set up CGAL
RUN cd cgal
RUN cmake .
RUN make
RUN make install
RUN cd ..

# Set up Suitesparse
RUN cd suitesparse
RUN make
RUN make install
RUN cd ..

# Set up libigl
RUN cd libigl
RUN cd external
RUN cd nanogui
RUN mkdir build
RUN cd build
RUN cmake ..
RUN make
RUN cd ..
RUN cd ..
RUN cd ..
RUN cd tutorial
# TODO :: IMPORT CMAKELISTS.TXT FROM ELSEWHERE
# TODO :: MOVE 3dUnion TO THIS FOLDER
RUN mkdir build
RUN cd build
RUN cmake ..
RUN make
RUN cd ..
RUN cd ..
RUN cd ..

# Set up PyMesh
RUN cd PyMesh
RUN cd third_party
RUN mkdir build
RUN cd build
RUN cmake ..
RUN make
RUN make install
RUN cd ..
RUN cd ..
RUN mkdir build
RUN cd build
RUN cmake ..
RUN make
RUN make src_tests
RUN make tools
RUN make tools_tests
