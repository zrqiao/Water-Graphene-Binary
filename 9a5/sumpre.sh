#!/bin/bash  

export filename="density_dis9a5"
export MDCURRENTJOB=14

#/share/home/yanglj/amber14/bin/pmemd.MPI -O -i min.in -p $filename.parm7 -c $filename.rst7 -o $filename\_min.out -r $filename\_min.rst7 -ref $filename.rst7
#/share/home/yanglj/amber14/bin/pmemd.MPI -O -i heat360.in -p $filename.parm7 -c $filename\_min.rst7 -o $filename\_360.out -r $filename\_360.rst7 -x $filename\_360.nc -ref $filename\_min.rst7
#/share/home/yanglj/amber14/bin/pmemd.MPI -O -i eq360.in -p $filename.parm7 -c $filename\_360.rst7 -o $filename\_360md.out -r $filename\_360md.rst7 -x $filename\_360md.nc -ref $filename\_360.rst7
#/share/home/yanglj/amber14/bin/pmemd.MPI -O -i cool280.in -p $filename.parm7 -c $filename\_360md.rst7 -o $filename\_cool280.out -r $filename\_cool280.rst7 -x $filename\_cool280.nc  -ref $filename\_360md.rst7
#/share/home/yanglj/amber14/bin/pmemd.MPI -O -i NPT.in -p $filename.parm7 -c $filename\_cool280.rst7 -o $filename\_0.out -r $filename\_0.rst7 -x $filename\_0.nc -ref $filename\_cool280.rst7

#/share/home/yanglj/amber14/bin/pmemd.MPI -O -i NPT.in -p $filename.parm7 -c $filename\_0.rst7 -o $filename\_1.out -r $filename\_1.rst7 -x $filename\_1.nc -ref $filename\_0.rst7


while [ $MDCURRENTJOB -le 60 ]; do
	let MDINPUT="$MDCURRENTJOB -1"
        mpirun -np 6 sander.MPI -O -i NVT_entropy.in -o $filename\_$MDCURRENTJOB.out -p $filename.parm7 -c $filename\_$MDINPUT.rst7  -r $filename\_$MDCURRENTJOB.rst7 -x $filename\_$MDCURRENTJOB.nc -ref $filename\_$MDINPUT.rst7
        let MDCURRENTJOB="$MDCURRENTJOB +1"
	done
