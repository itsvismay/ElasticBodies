FROM ubuntu:16.04
MAINTAINER zack441@mac.com

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
RUN echo 'INSTALLED EVERYTHING'
RUN apt install libgmp3-dev
RUN apt install libmpfr-dev


# git clones
RUN git clone --recursive https://github.com/qnzhou/PyMesh.git
RUN git clone https://github.com/chokkan/liblbfgs.git
RUN git clone --recursive https://github.com/libigl/libigl.git
RUN git clone https://github.com/jluttine/suitesparse.git
RUN git clone https://github.com/itsvismay/ElasticBodies.git
RUN git clone --recursive https://github.com/CGAL/cgal.git

# Anaconda
RUN echo ‘wget https://repo.continuum.io/archive/Anaconda2-4.2.0-Linux-x86_64.sh’
RUN wget https://repo.continuum.io/archive/Anaconda2-4.2.0-Linux-x86_64.sh
RUN echo ‘Anaconda2-4.2.0-Linux-x86_64.sh -b’
RUN bash Anaconda2-4.2.0-Linux-x86_64.sh -b

# Needed for PyMesh
RUN echo ‘wget https://bootstrap.pypa.io/ez_setup.py -O - | python’
RUN wget https://bootstrap.pypa.io/ez_setup.py -O - | python

# Moving Files Around
RUN echo ‘cp -r ElasticBodies/Pipeline/3dUnion libigl/tutorial’
cp -r ElasticBodies/Pipeline/3dUnion libigl/tutorial
# TODO :: IMPORT CMAKELISTS.TXT FROM ELSEWHERE

RUN ls

# Set up CGAL
RUN echo ‘cd cgal’
RUN cd cgal
RUN echo ‘cmake .’
RUN cmake .
RUN echo ‘make’
RUN make
RUN echo ‘make install’
RUN make install
RUN echo ‘cd ..’
RUN cd ..

# Set up Suitesparse
RUN echo ‘cd suitesparse’
RUN cd suitesparse
RUN echo ‘make’
RUN make
RUN echo ‘make install’
RUN make install
RUN echo ‘cd ..’
RUN cd ..

# Set up libigl
RUN echo ‘cd libigl’
RUN cd libigl
RUN echo ‘cd external’
RUN cd external
RUN echo ‘cd nanogui’
RUN cd nanogui
RUN echo ‘mkdir build’
RUN mkdir build
RUN echo ‘cd build’
RUN cd build
RUN echo ‘cmake ..’
RUN cmake ..
RUN echo ‘make’
RUN make
RUN echo ‘cd ..’
RUN cd ..
RUN echo ‘cd ..’
RUN cd ..
RUN echo ‘cd ..’
RUN cd ..
RUN echo ‘cd tutorial’
RUN cd tutorial
RUN echo ‘mkdir build’
RUN mkdir build
RUN echo ‘cd build’
RUN cd build
RUN echo ‘cmake ..’
RUN cmake ..
RUN echo ‘make’
RUN make
RUN echo ‘cd ..’
RUN cd ..
RUN echo ‘cd ..’
RUN cd ..
RUN echo ‘cd ..’
RUN cd ..

# Set up PyMesh
RUN echo ‘cd PyMesh’
RUN cd PyMesh
RUN echo ‘export PYMESH_PATH=$(pwd)’
RUN export PYMESH_PATH=$(pwd)
RUN echo ‘cd third_party’
RUN cd third_party
RUN echo ‘mkdir build’
RUN mkdir build
RUN echo ‘cd build’
RUN cd build
RUN echo ‘cmake ..’
RUN cmake ..
RUN echo ‘make’
RUN make
RUN echo ‘make install’
RUN make install
RUN echo ‘cd ..’
RUN cd ..
RUN echo ‘cd ..’
RUN cd ..
RUN echo ‘cmake ..’
RUN cmake ..
RUN echo ‘make’
RUN make
RUN echo ‘make src_tests’
RUN make src_tests
RUN echo ‘make tools’
RUN make tools
RUN echo ‘make tools_tests’
RUN make tools_tests
RUN echo ‘./setup.py build’
RUN ./setup.py build
RUN echo ‘./setup.py install’
RUN ./setup.py install
RUN echo ‘cd ..’
RUN cd ..

# Set up liblbfgs

# Set up ElasticBodies
