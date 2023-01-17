#!/usr/bin/env zsh
{
    echo 'AUC on dataset 1 predicting with ESP:'
    auc ESP_output.tsv -c reaction -p pred_ESP | mlr --t2d cat
    echo 'When filtering by nMetabInESP>0:'
    mlr -t filter '$nMetabInESP > 0' ESP_output.tsv | auc -c reaction -p pred_ESP | mlr --t2d cat
} > aucs.txt
