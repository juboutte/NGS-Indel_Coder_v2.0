# -*- coding: us-ascii -*-

import sys

### Usage ###

#1-parsing_Samtools_depth-files.py - Parse a Samtools outputfile generated using
# depth option. Generate 1 output file per sample for each locus
#For example, a depth file generated using 20 species and 300 locus will generate
#until 6,000 files.
#each output file contains read depth for each position for one species/one locus

#By Julien Boutte, April 2019
#Copyright (c) 2019 Julien Boutte.
#Version 1.0.0

#This program is free software: you can redistribute it and/or modify it under
#the terms of the GNU General Public License as published by the Free Software
#Foundation, either version 3 of the License, or (at your option) any later
#version. A copy of this license is available at <http://www.gnu.org/licenses/>.
#Great effort has been taken to make this software perform its said
#task, however, this software comes with ABSOLUTELY NO WARRANTY,
#not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

inputfile=open(sys.argv[1],'r')

test=0
for line in inputfile:
    if test==0:
        temp=line.split('\t')
        tempname=temp[0]
        output=open(tempname+'.txt','w')
        output.write(line)
        test=1
    else:
        temp=line.split('\t')
        if temp[0]==tempname:
            output.write(line)
        else:
            output.close()
            tempname=temp[0]
            output=open(tempname+'.txt','w')            
            output.write(line)
output.close()

