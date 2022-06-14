#!/usr/bin/env zsh
# "Functional and informatics analysis enables glycosyltransferase activity prediction" by Yang et al 2018
# 
# https://www.nature.com/articles/s41589-018-0154-9
# data made available at https://ora.ox.ac.uk/objects/uuid:1b174bc0-4058-4057-8db4-59872c2b6d99
wget -O GTPredict_All_Data.zip 'https://ora.ox.ac.uk/objects/uuid:1b174bc0-4058-4057-8db4-59872c2b6d99/download_file?file_format=zip&safe_filename=GTPredict_All_Data.zip&type_of_work=Dataset'
unzip GTPredict_All_Data.zip
rm -r __MACOSX
# selected data files have then been moved to this folder.
