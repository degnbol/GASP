#!/usr/bin/env zsh
ROOT=`git root`
echo "
# Putting GT in your paths
export PATH=\"\$PATH:$ROOT/src:$ROOT/tools/degnlib\"
export PYTHONPATH=\"\$PYTHONPATH:$ROOT:$ROOT/tools/degnlib\"
"
