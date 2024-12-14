#!/usr/bin/env python3
from pymol import cmd

# after pymol 2.5 pymol started ignoring atoms with flag ignore,
# which was on all atoms
cmd.flag("ignore", "all", "clear")

print("cid\tSESA")
for name in cmd.get_names():
    area = cmd.get_area(name)
    print(f"{name}\t{area}")


