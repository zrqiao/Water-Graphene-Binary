#!/bin/bash

#export fln="water_hydrophobicball_mix300_test"
#export parm="prmt7"
#export rst='rst7'
#export in_rst="n0"
#export suffix="n"
#export in_1st="nmd.in"
#export in_rest="nmd.in"
#export out=1
#export end=30

export frame=0
export final=9900
export distance="9a5" 

cd 333

while [ $frame -le $final ]; do
	echo "frame: $frame"
	sander -O -i NVT.in -o wat\_c\_$distance\_$frame.out -p wat\_c.parm7 -c wat\_c\_$distance\_$frame.rst7 
	sander -O -i NVT.in -o wat\_graphene\_$distance\_$frame.out -p density\_dis$distance.parm7 -c wat\_graphene\_$distance\_$frame.rst7
	grep "EPtot" wat\_c\_$distance\_$frame.out | awk '{print $9}' >> ../wat_c.dat
	grep "EPtot" wat\_graphene\_$distance\_$frame.out | awk '{print $9}' >> ../wat_graphene.dat
	let frame="$frame +100" 
done

cd ..

paste wat_c.dat wat_graphene.dat | awk '{print $2-$1}' > interaction_energy.dat


#sander -O -i min.in -o min.out -c water_ion_graphene.rst7 -p water_ion_graphene.parm7 -r min.rst7 
#sander -O -i heat360.in -o heat360.out -c min.rst7 -p water_ion_graphene.parm7 -r heat360.rst7 
#sander -O -i eq360.in -o eq360.out -c heat360.rst7 -p water_ion_graphene.parm7 -r eq360.rst7 
#sander -O -i cool300.in -o cool300.out -c eq360.rst7 -p water_ion_graphene.parm7 -r cool300.rst7
#sander -O -i eq300.in -o eq300.out -c cool300.rst7 -p water_ion_graphene.parm7 -r eq300.rst7
#sander -O -i NVT.in -o water_ion_result_1.out -p water_ion_graphene.parm7 -c eq300.rst7 -r water_ion_result_1.rst7 -x water_ion_result_1.nc


###############################################

#cp eq300.rst7 $fln\_n0.rst7

#let in="$out -1"

#  while [ $out -le $end ]; do
#	    ~/amber14/bin/pmemd.MPI -O -i $in_rest -p $fln.$parm -c $fln\_$suffix$in.rst7 -o $fln\_$suffix$out.out -r $fln\_$suffix$out.$rst -x $fln\_$suffix$out.nc 

    
#    let in="$in +1"
#    let out="$in +1"
#  done


###############################################

#let in="$out -1"
#  while [ $out -le $end ]; do
#             sander -O -i $in_rest -p $fln.$parm -c $fln\_$suffix$in.rst7 -o $fln\_$suffix$out.out -r $fln\_$suffix$out.$rst -x $fln\_$suffix$out.nc
#  let in="$in +1"
#  let out="$in +1"
#done

#sander -O -i NPT.in -o $fln.out -p $fln.parm7 -c $fln\_$in.rst7 -r water_tip3p_result.rst7 -x water_tip3p.nc
