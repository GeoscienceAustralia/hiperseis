# Created by Fei Zhang on 2019-12-13 in VDI

# Step-1 set up install ellip-corr
conda create -n hiperseispy37 python=3.7
conda activate hiperseispy37
python -V
pip list

export ELLIPCORR=/g/data/ha3/fxz547/Githubz/hiperseis/ellip-corr

cd $ELLIPCORR
python setup.py install
# this will install numpy into the env hiperseispy37


python ../seismic/traveltime/test_ellipcorr.py
###################################################################################################
#$ conda list
## packages in environment at /g/data1a/ha3/fxz547/miniconda3/envs/hiperseispy37:
##
## Name                    Version                   Build  Channel
#_libgcc_mutex             0.1                        main    conda-forge
#bzip2                     1.0.8                h516909a_2    conda-forge
#ca-certificates           2019.11.28           hecc5488_0    conda-forge
#certifi                   2019.11.28               py37_0    conda-forge
#cython                    0.29.14                  pypi_0    pypi
#ellip-corr                0.0.1                    pypi_0    pypi
#ellipcorr                 0.0.0                    pypi_0    pypi
#ld_impl_linux-64          2.33.1               h53a641e_7    conda-forge
#libffi                    3.2.1             he1b5a44_1006    conda-forge
#libgcc-ng                 9.2.0                hdf63c60_0    conda-forge
#libstdcxx-ng              9.2.0                hdf63c60_0    conda-forge
#ncurses                   6.1               hf484d3e_1002    conda-forge
#openssl                   1.1.1d               h516909a_0    conda-forge
#pip                       19.3.1                   py37_0    conda-forge
#python                    3.7.3                h357f687_2    conda-forge
#readline                  8.0                  hf8c457e_0    conda-forge
#setuptools                42.0.2                   py37_0    conda-forge
#sqlite                    3.30.1               hcee41ef_0    conda-forge
#tk                        8.6.10               hed695b0_0    conda-forge
#wheel                     0.33.6                   py37_0    conda-forge
#xz                        5.2.4             h14c3975_1001    conda-forge
#zlib                      1.2.11            h516909a_1006    conda-forge

#-----------------------------------------------------------------------------------------------
#(hiperseispy37) fxz547@vdi-n23 /g/data/ha3/fxz547/Githubz/hiperseis/ellip-corr ((no branch))
#$ pip list
#Package    Version
#---------- -------------------
#certifi    2019.11.28
#Cython     0.29.14
#ellip-corr 0.0.1
#ellipcorr  0.0.0
#numpy      1.17.4
#pip        19.3.1
#setuptools 42.0.2.post20191201
#wheel      0.33.6
########################################################################################################

# After the ellip-corr is installed and tested working, then install other python packages

conda install --file requirements_conda_py37.txt

# then run the quick test script, which must be successful.

bash seismic/traveltime/test_sort_rays.sh