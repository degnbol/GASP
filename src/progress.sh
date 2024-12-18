#!/usr/bin/env zsh
# Instead of printing stdout, pipe it into this script and see progress.
# USAGE: CMD | progress.sh N [SKIP] 1>&2 
# N = total number of lines
# SKIP = number of initial lines to skip
# 1>&2 to print to stderr instead of default stdout
N=$1
SKIP=$2 # may be empty

function progress_bar () {
    i=$1
    char=$2
    # print progress as a fraction that replaces the line (\r)
    # and i is width adjusted base on N
    printf "\r%${#N}d/$N " $i
    # calculate number of columns of progress symbol to draw.
    # Subtract 2 for '/' and a space. Subtract 2*char length of N for the two printed numbers.
    let cols='(COLUMNS - 2 - 2*'$#N')*i/N'
    # a trick to print a character ($char) n times, where n=$cols.
    # Works for integers except 0
    if [[ "$i" > 0 ]]; then
        printf '%0.1s' $char{1..$cols}
    fi
}

if [ -n "$SKIP" ]; then
    # skip SKIP lines
    for i in {1..$SKIP}; do
        read line
    done
fi

progress_bar 0 '#'
read line
let i=1

# exit on eof (altho also on any other empty line)
until [ -z "$line" ]; do
    # allow lines starting with '#' through
    if [ "${line::1}" = '#' ]; then
        # first clear progress bar by overwriting spaces
        progress_bar $i ' '
        echo "\r$line"
    else
        let i+=1
    fi

    progress_bar $i '#'
    read line
done

# clear
printf '\33[2K\r'

