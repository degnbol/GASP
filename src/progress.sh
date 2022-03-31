#!/usr/bin/env zsh
# Instead of printing stdout, pipe it into this script and see progress.
# USAGE: CMD | progress.sh N [SKIP]
# N = total number of lines
# SKIP = number of initial lines to skip
N=$1
SKIP=$2 # may be empty

echo -n "0/$N"

if [ -n "$SKIP" ]; then
    # skip SKIP lines
    for i in {1..$SKIP}; do
        read line
    done
fi

let i=0
while true; do
    let i++
    
    read line
    # exit on eof (altho also on any other empty line)
    if [ -z "$line" ]; then
        break
    fi
    
    # print progress as a fraction that replaces the line (\r)
    # and i is width adjusted base on N
    printf "\r%${#N}d/$N " $i
    # calculate number of columns of progress symbol to draw
    let cols='(COLUMNS - 2 - 2*'$#N')*i/N'
    # a hacky method to repeat print a character, here "#"
    printf '%0.1s' "#"{1..$cols}
done

# clear
echo -n "\r"

