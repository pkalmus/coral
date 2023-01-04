#!/bin/bash

model=0

while [ $model -le 127 ]
do
 echo "model: $model"

 screen -d -m -S mod$model
 screen -m -S mod$model -X stuff 'python weight_trend.py params '"$model"' \n'
 ((model++))
done
