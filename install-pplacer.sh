#!/usr/bin/env bash

if [[ "$OSTYPE" == "darwin"* ]]; then
  DOWNLOAD_URL="https://github.com/matsen/pplacer/releases/download/v1.1.alpha17/pplacer-Darwin-v1.1.alpha17.zip"
elif [[ "$OSTYPE" == "linux"* ]]; then
  DOWNLOAD_URL="https://github.com/matsen/pplacer/releases/download/v1.1.alpha17/pplacer-Linux-v1.1.alpha17.zip"
else
  echo "Detected OS version (${OSTYPE}) is not supported. Aborting."
  exit 1
fi

echo "Fetching pplacer binaries from ${DOWNLOAD_URL}..."
curl -L "${DOWNLOAD_URL}" > pplacer.zip

echo "Extracting..."
unzip pplacer.zip
rm pplacer.zip

if [[ "$PREFIX" == "" ]]; then
  echo "Setting PREFIX=$CONDA_PREFIX"
  PREFIX="$CONDA_PREFIX"
fi

echo "Installing pplacer in $PREFIX..."
if [[ ! -d "$PREFIX/bin/" ]]; then
  mkdir $PREFIX/bin/
fi
mv pplacer*/guppy $PREFIX/bin/
mv pplacer*/pplacer $PREFIX/bin/
mv pplacer*/rppr $PREFIX/bin/

mkdir $PREFIX/bin/scripts/
mv pplacer*/scripts/* $PREFIX/bin/scripts/
rm -r pplacer*

echo "Testing installation..."
if [[ $(which pplacer) == "$PREFIX/bin"* ]]; then
  echo "Success!"
# TODO: make sure later that this really is not necessary (try to install with conda on macOS and Ubuntu)
#  pplacer --version
#  retVal=$?
#  if [[ $retVal -ne 0 ]]; then
#    echo "pplacer was installed in the correct location but there was a problem with the installation."
#    echo "See below for the potential error message:"
#    pplacer --version
#    exit 1
#  else
#    echo "Success!"
#  fi
else
  echo "Installation failed."
  exit 1
fi
