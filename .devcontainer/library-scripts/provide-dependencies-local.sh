# Stop the program at the first error
set -e

# Install the necessary programs for developement
# The programs listed in provide-dependencies-ci.sh are omitted here
# therefore both scripts must be run in order to set up the local developement environment
apt-get install -y \
    clang-format \
    git \
    doxygen \
    graphviz \
    gdb \
    libboost-all-dev


# For generating latex documentation with doxygen, texlive and ghostscript also need to be installed
# See https://www.doxygen.nl/manual/install.html for more info