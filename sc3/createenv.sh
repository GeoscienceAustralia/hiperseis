#!/usr/bin/env bash
set -e
set -o pipefail
# Setup virtualenv
echo '' >> $HOME/.bashrc
echo "source /usr/bin/virtualenvwrapper.sh" >> $HOME/.bashrc
source $HOME/.bashrc

# manage env variables, if ilocroot exists then skip
if [ -z ${ILOCROOT+x} ]; then
    cat $HOME/passive-seismic/sc3/sc3_envs.txt >> $HOME/.bashrc;
    source $HOME/.bashrc
fi

# clone passive-seismic in home directory after git install
if [ ! -d "$HOME/passive-seismic" ]; then
	CWD=`pwd`
	cd $HOME
	git clone https://github.com/GeoscienceAustralia/passive-seismic.git
	cd $CWD
fi

echo "Installing passive seismic software....."


if [ -z ${CIRCLECI+x} ];
    then
    source $HOME/.bashrc && \
    cd $HOME/passive-seismic && \
    git submodule update --init --recursive --remote && \
    mkvirtualenv --system-site-packages seismic && \
    workon seismic && \
    pip install --ignore-installed ipython && \
    export PATH=$HOME/.virtualenvs/seismic/bin/:$PATH && \
    export ELLIPCORR=$PWD/ellip-corr && \
    env GEOS_DIR=/usr/include/ pip install --process-dependency-links -e .[dev] && \
    $HOME/passive-seismic/iloc_rstt/install_iloc_rstt.sh;  # install iloc and rstt
fi
