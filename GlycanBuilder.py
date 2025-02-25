#!/usr/bin/python

import sys
import math
from math import *
import argparse

# Argument Parser
parser = argparse.ArgumentParser(description='Create a custom coordinate and topology file for coarse-grained linear glycan')
parser.add_argument( "-l", "--chainlength", type=int, default=5, help='Chain length (Corresponds to number of monomers per glycan chain)')
parser.add_argument( "-o", "--output",   type=str,   default=None,   help='Name of output file' )
args = parser.parse_args()

l = args.chainlength

if args.output == None:
	output = "Glycan-"+str(l)+"monomers"
else:
        output = args.output

# Total number of atoms 
natoms = l*3

print( "Generated Martini model for glycan with %s monomers" % (l))

# Function to define atom names
# Code obtained from https://stackoverflow.com/questions/60039572/how-to-increment-alphanumeric-number-in-python

def excel_format(num):
    res = ""
    while num:
        mod = (num - 1) % 26
        res = chr(65 + mod) + res
        num = (num - mod) // 26
    return res

def full_format(num, d=3):
    chars = num // (10**d-1) + 1 # this becomes   A..ZZZ
    digit = num %  (10**d-1) + 1 # this becomes 001..999
    return excel_format(chars) + "{:0{}d}".format(digit, d)


# CREATE GRO FILE

# Opens the gro file for writing

structure_file = open(output+".gro", 'w')


# Writes header for gro file

structure_file.write( "Glycan-%smonomers\n" % (l) )
structure_file.write( "  %d\n" % (natoms) )


# Writes coordinates for beads

bl1 = 0.376
bl2 = 0.329
bl3 = 0.276
bl4 = 0.220
bl = 0.356

for i in range(1,l+1):
	j=(i-1)*3+1
	x=bl*(i-1)
	at_name1 = full_format(j, d=3)
	at_name2 = full_format(j+1, d=3)
	at_name3 = full_format(j+2, d=3)
	if i%2==1:
		y1=bl1
		y2=-bl2
	else:
		y1=bl3
		y2=-bl4
	structure_file.write(  "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (1, "POL", at_name1, j, x, 0, 0, 0, 0, 0) )	
	structure_file.write(  "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (1, "POL", at_name2, j+1, x, y1, 0, 0, 0, 0) )	
	structure_file.write(  "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % (1, "POL", at_name3, j+2, x, y2, 0, 0, 0, 0) )	

# Box dimensions

xdim = l*bl + 1
ydim = bl1+bl2+1
zdim = 2

structure_file.write(  "   %u  %u  %u\n" % (xdim, ydim, zdim) )

structure_file.close()


# CREATE TOPOLOGY FILE

# Opens the itp file for writing

topology_file = open(output+".itp", 'w')

# Writes header for the itp file

topology_file.write( "; \n; Martini topology for glycan with %s monomers\n" % (l))
topology_file.write( "; \n; Topology generated using GlycanBuilder.py \n" )
topology_file.write( "; Written by Raman Preet Singh \n\n")

topology_file.write( "[ moleculetype ]\n" )
topology_file.write( "; Name	 nrexcl\n" )
topology_file.write( "POL 1\n")

# Atoms

topology_file.write( "\n[ atoms ]\n" )
topology_file.write( "; id	 type	 resnr	 residue	 atom	 cgnr	 charge	 mass\n" )

for i in range(1,l+1):
	j=(i-1)*3+1
	at_name1 = full_format(j, d=3)
	at_name2 = full_format(j+1, d=3)
	at_name3 = full_format(j+2, d=3)
	if i==1:
		m=59.0448
	else:
		m=43.0454
	topology_file.write( "%3d    %4s   %6d   POL    %5s     %3d     0       %2.4f\n" % (j, "P2", 1, at_name1, i, m) )
	topology_file.write( "%3d    %4s   %6d   POL    %5s     %3d     0       60.0528\n" % (j+1, "P4", 1, at_name2, i) )
	topology_file.write( "%3d    %4s   %6d   POL    %5s     %3d     0       60.0528\n" % (j+2, "P1", 1, at_name3, i) )


# Constraints

topology_file.write( "\n[ constraints ]\n" )
topology_file.write( "; i	 j	  funct	 length	 force\n" )

for i in range(1,l+1):
	j=(i-1)*3+1
	if i%2==1:
		y1=bl1
		y2=bl2
	else:
		y1=bl3
		y2=bl4
	topology_file.write( "     %3d     %3d       1   %4.3f    30000; P2-P4\n" % (j, j+1, y1) )
	topology_file.write( "     %3d     %3d       1   %4.3f    30000; P2-P1\n" % (j, j+2, y2) )
	if i<l:
		topology_file.write( "     %3d     %3d       1   %4.3f    30000; P2-P2\n" % (j, j+3, bl) )

# Angles

topology_file.write( "\n[ angles ]\n" )
topology_file.write( "; i	 j	 k	 funct	 angle	 force\n" )

for i in range(1,l+1):
	j=(i-1)*3+1
	if i<l:
		topology_file.write( "     %3d     %3d     %3d       2     54.0     80; P1-P2-P2\n" % (j+2, j, j+3) )
		topology_file.write( "     %3d     %3d     %3d       2    124.0     200; P4-P2-P2\n" % (j+1, j, j+3) )
		topology_file.write( "     %3d     %3d     %3d       2     44.0     500; P2-P2-P4\n" % (j, j+3, j+4) )
		topology_file.write( "     %3d     %3d     %3d       2     67.0     800; P2-P2-P1\n" % (j, j+3, j+5) )
	if i<l-1:
		topology_file.write( "     %3d     %3d     %3d       2    136.0     500; P2-P2-P2\n" % (j, j+3, j+6) )

# Dihedrals
topology_file.write( "\n[ dihedrals ]\n" )
topology_file.write( ";  ai    aj    ak    al funct    angle fc    mult \n" )

for i in range(2,l+1):
	j=(i-1)*3+1
	if i<l and i%2==0:
		topology_file.write( "     %3d     %3d     %3d     %3d       1     20.0     15       1; P1-P2-P2-P4\n" % (j+2, j, j+3, j+4) )
		topology_file.write( "     %3d     %3d     %3d     %3d       1     55.0      5       1; P1-P2-P2-P1\n" % (j+2, j, j+3, j+5) )
		topology_file.write( "     %3d     %3d     %3d     %3d       1     42.0      5       1; P4-P2-P2-P4\n" % (j+1, j, j+3, j+4) )

