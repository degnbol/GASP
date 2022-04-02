#!/usr/bin/env zsh
ROOT=`git root`
echo "
# Copy-paste these two lines to your ~/.zshrc to make sure all scripts can be found for GT
export PATH=\"\$PATH:$ROOT/src:$ROOT/tools/degnlib\"
export PYTHONPATH=\"\$PYTHONPATH:$ROOT:$ROOT/tools/degnlib\"
"
