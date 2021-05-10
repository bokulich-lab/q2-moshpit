#!/bin/sh

# first arg should be Q2 version to test against (e.g. 2021.8)
# second arg should be OS version (osx or linux)
Q2V="$1"
OSV="$2"
COL="\033[0;32m"
NC="\033[0m" # NoColor

if test $Q2V = ""
  then Q2V="2021.8"
fi

if test $OSV =  ""
  then OSV="osx"
fi

CONDA_ENVS=$(conda info | grep "envs directories" | sed "s/envs directories \: //" | xargs)
CONDA_MAIN=$(echo $CONDA_ENVS | sed "s/\/envs//")
echo "${COL}Detected conda path: ${CONDA_MAIN}${NC}"

echo "${COL}Preparing a new conda environment (qiime2-${Q2V}-buildtest)...${NC}"
conda create -y -n "qiime2-${Q2V}-buildtest" conda-build conda-verify
source "${CONDA_MAIN}/etc/profile.d/conda.sh"
conda activate "qiime2-${Q2V}-buildtest"
conda info

# build new package
echo "${COL}Starting build for QIIME2 ${Q2V} (${OSV})...${NC}"
conda build -c "https://packages.qiime2.org/qiime2/${Q2V}/staged" -c conda-forge -c bioconda -c defaults --override-channels --no-anaconda-upload --cache-dir conda_cache .

# test the build
echo "${COL}Testing the build...${NC}"
wget -O env.yml "https://raw.githubusercontent.com/qiime2/environment-files/master/${Q2V}/staging/qiime2-${Q2V}-py38-${OSV}-conda.yml"
conda env create -q -p "./testing-${Q2V}" --file env.yml
conda install -p "./testing-${Q2V}" -q -y -c "$CONDA_ENVS/qiime2-${Q2V}-buildtest/conda-bld/${OSV}-64" -c conda-forge -c bioconda -c defaults --override-channels --strict-channel-priority q2-plugin-name

conda activate "./testing-${Q2V}"
pytest --pyargs q2_plugin_name
conda deactivate

# clean up previous build tests
echo "${COL}Cleaning up...${NC}"
rm -rf "./testing-${Q2V}/*"
rm env.yml
conda build purge

conda deactivate
conda env remove -n "qiime2-${Q2V}-buildtest"

echo "${COL}All done!${NC}"
