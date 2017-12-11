#!/bin/sh


for i in $(seq 1 925)
do

	n=$((i%10))
        if [ $n -eq 0 ]
	then
		python main.py "$i"

	else
		python main.py "$i"
	fi
	echo "$i"

done
