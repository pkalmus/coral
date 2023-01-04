#!/bin/bash

model=0

while [ $model -le 126 ]
do
 echo "model: $model"

 screen -d -m -S mod$model
 screen -m -S mod$model -X stuff 'geo \n'
 screen -m -S mod$model -X stuff 'python weight_dhw.py params '"$model"' \n'
 ((model++))
done
