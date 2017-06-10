#!/bin/bash


export frame=0
export final=9900
export distance="9a5" 

cd confined_wat

while [ $frame -le $final ]; do
#	echo "frame: $frame"
#	sander -O -i NVT.in -o wat\_c\_$distance\_$frame.out -p wat\_c.parm7 -c wat\_c\_$distance\_$frame.rst7 
#	sander -O -i NVT.in -o wat\_graphene\_$distance\_$frame.out -p density\_dis$distance.parm7 -c wat\_graphene\_$distance\_$frame.rst7
	echo "frame: $frame"
	sander -O -i NVT.in -o wat\_confined\_$distance\_$frame.out -p wat\_confined.parm7 -c wat\_confined\_$distance\_$frame.rst7 
	grep "VDWAALS" wat\_confined\_$distance\_$frame.out | awk '{print $11}' >> ../wat_confined_VDWAALS.dat
	let frame="$frame +100" 
done

cd ..





