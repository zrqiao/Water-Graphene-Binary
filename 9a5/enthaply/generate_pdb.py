#usr/bin/python
#This file is generating crystal
import sys
import random

#ofile1 = open('wat_pdb','r')
ofile2 = open('wat_confined_pdb','r')
#ofile3 = open('wat.pdb','w')
ofile4 = open('wat_confined.pdb','w')


#distance = 6.5  
#graphene1_num = 1400
#graphnen2_num = 1400

	

#for line in ofile1:
#	line_string = line.rsplit('n')
#	list_string = line_string[0].split()
#	if list_string[0] == 'ATOM':
#		ofile3.write("%4s%7d%5s%4s%6d%12.3f%8.3f%8.3f\n"%(list_string[0],int(list_string[1]),list_string[2],list_string[3],int(list_string[4]),float(list_string[5]),float(list_string[6]),float(list_string[7])))
#	elif list_string[0] == 'TER':
#               ofile3.write('TER\n')

for line in ofile2:
	line_string = line.rsplit('n')
	list_string = line_string[0].split()
	if list_string[0] == 'ATOM':
		ofile4.write("%4s%7d%5s%4s%6d%12.3f%8.3f%8.3f\n"%(list_string[0],int(list_string[1]),list_string[2],list_string[3],int(list_string[4]),float(list_string[5]),float(list_string[6]),float(list_string[7])))
	elif list_string[0] == 'TER':
                ofile4.write('TER\n')

#ofile1.close()
ofile2.close()
#ofile3.close()
ofile4.close()
