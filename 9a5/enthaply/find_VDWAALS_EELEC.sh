#!/bin/bash


export frame=0
export final=9900
export distance="9a5" 

cd 333

while [ $frame -le $final ]; do
#	echo "frame: $frame"
#	sander -O -i NVT.in -o wat\_c\_$distance\_$frame.out -p wat\_c.parm7 -c wat\_c\_$distance\_$frame.rst7 
#	sander -O -i NVT.in -o wat\_graphene\_$distance\_$frame.out -p density\_dis$distance.parm7 -c wat\_graphene\_$distance\_$frame.rst7
	grep "VDWAALS" wat\_c\_$distance\_$frame.out | awk '{print $11}' >> ../wat_c_VDWAALS.dat
	grep "VDWAALS" wat\_graphene\_$distance\_$frame.out | awk '{print $11}' >> ../wat_graphene_VDWAALS.dat
	grep "EELEC" wat\_c\_$distance\_$frame.out | awk '{print $3}' >> ../wat_c_EELEC.dat
	grep "EELEC" wat\_graphene\_$distance\_$frame.out | awk '{print $3}' >> ../wat_graphene_EELEC.dat
	let frame="$frame +100" 
done

cd ..

paste wat_c_VDWAALS.dat wat_graphene_VDWAALS.dat | awk '{print $2-$1}' > VDWAALS_energy.dat
paste wat_c_EELEC.dat wat_graphene_EELEC.dat | awk '{print $2-$1}' > ELLEC_energy.dat



