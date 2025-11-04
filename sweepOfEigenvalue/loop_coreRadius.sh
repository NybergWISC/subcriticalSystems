#!/bin/bash

paramArray='100 125 150 175 200 225 240 241 242 250 275 300 325 350'
paramArray='125 150 175 200 225 240 241 242 250 275 300 325 350'
paramArray='350 375 400'
# paramArray='0.075 0.2'

FILE_END="_outTemp.txt"
mkdir outputs_temp

for value in $paramArray
do
	echo $value
	python3 sweepRadisu.py $value
	python3 readTally.py $value
	filename="$value$FILE_END"
	cp tallies.out outputs_temp/$filename
done

# for value in $paramArray
# do
#         python3 eigen_loop.py $value
# done
