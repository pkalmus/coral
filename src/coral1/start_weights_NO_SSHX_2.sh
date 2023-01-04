#!/bin/bash

model=64

while [ $model -le 127 ]
do
 echo "model: $model"

 screen -d -m -S mod$model
 screen -m -S mod$model -X stuff 'geo \n'
 screen -m -S mod$model -X stuff 'python weighting.py params '"$model"' \n'
 ((model++))
done
