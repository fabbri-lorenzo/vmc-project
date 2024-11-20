# Stop the program at the first error
set -e

# Install the necessary programs for CI
apt-get install -y \
    g++ \
    cmake \
    libtbb-dev \
    libboost-all-dev
