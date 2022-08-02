#!/usr/bin/env zsh
mlr -t --from ../10-experimentalValidation/experimentYield.tsv put '$reaction = $Yield >= 25 ? 1 : ($Yield < 12.5 ? 0 : 0.5)' + \
    filter '$reaction != 0.5' + uniq -f enzyme,reaction,cid |
    mlr --from /dev/stdin --from ../08-features/train.tsv cat > DON+betanidin-trainset.tsv
