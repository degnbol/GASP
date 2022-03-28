#!/usr/bin/env zsh
if [[ `uname` == "Darwin" ]]; then
    OS=MacOSX
elif [[ `uname` == "Linux" ]]; then
    OS=Linux
else
    OS=Windows
fi
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-$OS-x86_64.sh
bash Miniconda3-latest-$OS-x86_64.sh
rm Miniconda3-latest-$OS-x86_64.sh
