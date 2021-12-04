cat selected*.txt | tr ' ' '\n' | sort -u | tr '\n' ' ' > all_selected.arg
