#!/bin/bash
year=$1
month=1

while [ $month -le 12 ]
do
 echo "year/month: $year/$month"

 screen -d -m -S $year$month
 screen -m -S $year$month -X stuff 'python dailyISDmatchup.py '"$year"' '"$month"' \n'
 ((month++))
done

