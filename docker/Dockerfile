FROM ubuntu:16.04
MAINTAINER zack441@mac.com

# required libraries
RUN echo 'apt-get update' && \
    apt-get update && \
    echo 'install apt-utils' && \
    apt-get install -y apt-utils && \
    echo 'install aptitude' && \
    apt-get install -y aptitude && \
    echo 'install vim' && \
    apt-get install -y vim && \
    echo 'install git' && \
    apt-get install -y git && \
    echo 'install build-essential' && \
    apt-get install -y build-essential && \
    echo 'install libgl1-mesa-dev' && \
    apt-get install -y libgl1-mesa-dev && \
    echo 'install freeglut3-dev' && \
    apt-get install -y freeglut3-dev && \
    echo 'install cmake' && \
    apt-get install -y cmake && \
    echo 'install libblas-dev' && \
    apt-get install -y libblas-dev && \
    echo 'install libboost-dev' && \
    apt-get install -y libboost-dev && \
    echo 'install libmpfr4 and libmpfr4-dbg' && \
    apt-get install -y libmpfr4 libmpfr4-dbg && \
    echo 'install libopenlbas-dev' && \
    apt-get install -y libopenblas-dev && \
    echo 'install libpng12-dev' && \
    apt-get install -y libpng12-dev && \
    echo 'install python-matplotlib' && \
    apt-get install -y python-matplotlib && \
    echo 'install python-scipy' && \
    apt-get install -y python-scipy && \
    echo 'install python3-dev' && \
    apt-get install -y python3-dev && \
    echo 'install swig' && \
    apt-get install -y swig && \
    echo 'install zlib1g-dev' && \
    apt-get install -y zlib1g-dev && \
    echo 'install libxcursor-dev' && \
    apt-get install -y libxcursor-dev && \
    echo 'install libxi-dev' && \
    apt-get install -y libxi-dev && \
    echo 'install libxinerama-dev' && \
    apt-get install -y libxinerama-dev && \
    echo 'install libxrandr-dev' && \
    apt-get install -y libxrandr-dev && \
    echo 'install libcgal-dev' && \
    apt-get install -y libcgal-dev && \
    echo 'purge libcgal-dev' && \
    apt-get purge -y libcgal-dev && \
    echo 'install libeigen3-dev' && \
    apt-get install -y libeigen3-dev && \
    echo 'install openscad and openscad-dbg' && \
    apt-get install -y openscad openscad-dbg && \
    echo 'install slic3r' && \
    apt-get install -y slic3r && \
    echo 'install meshlab' && \
    apt-get install -y meshlab && \
    echo 'install clang' && \
    apt-get install -y clang && \
    echo 'install libtool' && \
    apt-get install -y libtool && \
    echo 'install libtbb-dev' && \
    apt-get install -y libtbb-dev && \
    echo 'install libqt4-dev' && \
    apt-get install -y libqt4-dev && \
    echo 'install wget' && \
    apt-get install -y wget && \
    echo ‘install libgmp3-dev’ && \
    apt-get install -y libgmp3-dev && \
    echo ‘install libmpfr-dev’ && \
    apt-get install -y libmpfr-dev && \
    echo ‘install m4’ && \
    apt-get install -y m4 && \
    echo ‘install autoconf’ && \
    apt-get install -y autoconf && \


# git clones
RUN git clone --recursive https://github.com/qnzhou/PyMesh.git && \
    git clone https://github.com/chokkan/liblbfgs.git && \
    git clone --recursive https://github.com/libigl/libigl.git && \
    git clone https://github.com/jluttine/suitesparse.git && \
    git clone https://github.com/itsvismay/ElasticBodies.git && \
    git clone --recursive https://github.com/CGAL/cgal.git

# Anaconda
RUN echo ‘wget https://repo.continuum.io/archive/Anaconda2-4.2.0-Linux-x86_64.sh’ && \
    wget https://repo.continuum.io/archive/Anaconda2-4.2.0-Linux-x86_64.sh && \
    echo ‘Anaconda2-4.2.0-Linux-x86_64.sh -b’ && \
    bash Anaconda2-4.2.0-Linux-x86_64.sh -b

# Needed for PyMesh
RUN python --version && \
    echo ‘wget https://bootstrap.pypa.io/ez_setup.py -O - python’ && \
    wget https://bootstrap.pypa.io/ez_setup.py -O - | python

# Moving Files Around
RUN echo ‘cp -r ElasticBodies/Pipeline/3dUnion libigl/tutorial’ && \
    cp -r ElasticBodies/Pipeline/3dUnion libigl/tutorial && \
    cp ElasticBodies/Pipeline/CMakeLists.txt libigl/tutorial/CMakeLists.txt

RUN ls

# Set up CGAL
RUN echo ‘cd cgal’ && \
    cd cgal && \
    echo ‘cmake .’ && \
    cmake . && \
    echo ‘make’ && \
    make && \
    echo ‘make install’ && \
    make install && \
    echo ‘cd ..’ && \
    cd ..

# Set up Suitesparse
RUN echo ‘cd suitesparse’ && \
    cd suitesparse && \
    echo ‘make’ && \
    make && \
    echo ‘make install’ && \
    make install && \
    echo ‘cd ..’ && \
    cd ..

# Set up libigl
RUN echo ‘cd libigl’ && \
    cd libigl && \
    echo ‘cd external’ && \
    cd external && \
    echo ‘cd nanogui’ && \
    cd nanogui && \
    echo ‘mkdir build’ && \
    mkdir build && \
    echo ‘cd build’ && \
    cd build && \
    echo ‘cmake ..’ && \
    cmake .. && \
    echo ‘make’ && \
    make && \
    echo ‘cd ..’ && \
    cd .. && \
    echo ‘cd ..’ && \
    cd .. && \
    echo ‘cd ..’ && \
    cd .. && \
    echo ‘cd tutorial’ && \
    cd tutorial && \
    echo ‘mkdir build’ && \
    mkdir build && \
    echo ‘cd build’ && \
    cd build && \
    echo ‘cmake ..’ && \
    cmake .. && \
    echo ‘make’ && \
    make && \
    echo ‘cd ..’ && \
    cd .. && \
    echo ‘cd ..’ && \
    cd .. && \
    echo ‘cd ..’ && \
    cd ..

# Set up PyMesh
RUN echo ‘cd PyMesh’ && \
    cd PyMesh && \
    echo ‘export PYMESH_PATH=$(pwd)’ && \
    export PYMESH_PATH=$(pwd) && \
    echo ‘cd third_party’ && \
    cd third_party && \
    echo ‘mkdir build’ && \
    mkdir build && \
    echo ‘cd build’ && \
    cd build && \
    echo ‘cmake ..’ && \
    cmake .. && \
    echo ‘make’ && \
    make && \
    echo ‘make install’ && \
    make install && \
    echo ‘cd ..’ && \
    cd .. && \
    echo ‘cd ..’ && \
    cd .. && \
    echo ‘cmake ..’ && \
    cmake .. && \
    echo ‘make’ && \
    make && \
    echo ‘make src_tests’ && \
    make src_tests && \
    echo ‘make tools’ && \
    make tools && \
    echo ‘make tools_tests’ && \
    make tools_tests && \
    echo ‘./setup.py build’ && \
    ./setup.py build && \
    echo ‘./setup.py install’ && \
    ./setup.py install && \
    echo ‘cd ..’ && \
    cd ..

# Set up liblbfgs
RUN echo ‘cd liblbfgs’ && \
    cd liblbfgs && \
    echo ‘cd ..’ && \
    cd ..

# Set up ElasticBodies
RUN echo ‘cd ElasticBodies’ && \
    cd ElasticBodies && \
    echo ’cmake .’ && \
    cmake . && \
    echo ‘make’ && \
    make && \
    echo ‘cd Pipeline/SimPrep’ && \
    cd Pipeline/SimPrep && \
    echo ‘make’ && \
    make && \
    echo ‘cd ..’ && \
    cd .. && \
    echo ‘cd ..’ && \
    cd .. && \
    echo ‘cd ..’ && \
    cd ..

# Run the Simulation

