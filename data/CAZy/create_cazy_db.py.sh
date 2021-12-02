#!/usr/bin/env zsh
# The create_cazy_db.py script was installed from pip and then modified to restrict parsing to GT1.
pip install --upgrade cazy-parser
# fix import bug so we get urllib.request loaded properly:
sed -i '' 's/^import urllib$/import urllib.request/' ~/miniconda3/envs/gt/lib/python3*/site-packages/cazy_parser/create_cazy_db.py
# I also commented out families that are not Glycosyltransferase since we don't care about those and moved the file here.
