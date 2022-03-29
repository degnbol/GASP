#!/usr/bin/env python3
import cirpy
import sys

for line in sys.stdin:
    print(cirpy.resolve(line.strip(), "smiles"))

