#!/usr/bin/env zsh
LIB="`git root`/tools/degnlib/subfunctions"
hmmalign --trim --amino --outformat A2M gtpred.muscle.hmm `git root`/data/inhouse/hts_ugt/hts_ugt_constructs.fa |
    $LIB/fasta_delete_lower.py > hts_ugt-gtpred.muscle.hmm.fa
