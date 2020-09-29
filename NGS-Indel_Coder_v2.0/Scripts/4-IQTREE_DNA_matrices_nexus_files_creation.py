# -*- coding: us-ascii -*-

### Usage ###

#By Julien Boutte, September 2020
#Copyright (c) 2020 Julien Boutte.
#Version 1.0.0

#This program is free software: you can redistribute it and/or modify it under
#the terms of the GNU General Public License as published by the Free Software
#Foundation, either version 3 of the License, or (at your option) any later
#version. A copy of this license is available at <http://www.gnu.org/licenses/>.
#Great effort has been taken to make this software perform its said
#task, however, this software comes with ABSOLUTELY NO WARRANTY,
#not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

import sys
import os
from optparse import OptionParser

opts=OptionParser()
opts.add_option('-y',type='str',dest='phy',default='')
opts.add_option('-p',type='str',dest='part',default='')
opts.add_option('-f',type='str',dest='folder',default='')
opts.add_option('-i',type='str',dest='indel_phy',default='')
opts.add_option('-d',type='str',dest='deletion',default='no')

(options, args)=opts.parse_args(sys.argv[1:])

phy=options.phy
part=options.part
folder=options.folder
indel_phy=options.indel_phy
deletion=options.deletion

#PART 1

#using files created by 2MATRIX create 2 .phy for IQ-tree
input_phy=open(phy,'r')
input_part=open(part,'r')

output_phy1=open(phy.replace('.phy','')+'_dna.phy','w')
#output_phy2=open(sys.argv[1].replace('.phy','')+'_indel.phy','w')

start_dna=0
end_dna=0

start_indel=0
end_indel=0

test=0
for line in input_part:
    if test==0:
        temp=line.replace('\n','').split('=')
        temp2=temp[1].split('-')
        start_dna=int(temp2[0])
        end_dna=int(temp2[1])        
        test=1
    else:
        temp=line.replace('\n','').split('=')
        temp2=temp[1].split('-')
        start_indel=int(temp2[0])
        end_indel=int(temp2[1])

input_part.close()

list_supp=[]

test=0
for line in input_phy:
    if test==0:
        test=1
    else:
        temp=line.split(' ')
        name=temp[0].split('_')
        #if temp[-1][start_indel-1:].count('0')==0 and temp[-1][start_indel-1:].count('1')==0 and temp[-1][start_indel-1:].count('-')==0:
        if len(temp[-1][start_indel-1:].replace('\n',''))==temp[-1][start_indel-1:].replace('\n','').count('?'):
            list_supp.append(name[0])

print list_supp

input_phy.close()
input_phy=open(phy,'r')

test=0
for line in input_phy:
    if test==0:
        test=1
        temp=line.split(' ')
        output_phy1.write(str(int(temp[0])-len(list_supp))+' '+str(end_dna-start_dna+1)+'\n')
        #output_phy2.write(str(int(temp[0])-len(list_supp))+' '+str(end_indel-start_indel+1)+'\n')
    else:
        temp=line.split(' ')
        name=temp[0].split('_')
        if name[0] not in list_supp:
            output_phy1.write(name[0]+' ')
            output_phy1.write(temp[-1][:start_indel-1].replace('\n','')+'\n')
            #output_phy2.write(name[0]+' ')
            #output_phy2.write(temp[-1][start_indel-1:].replace('\n','')+'\n')

input_phy.close()
output_phy1.close()
#output_phy2.close()

#PART 2

inputfile1=open(phy.replace('.phy','')+'_dna.phy','r') # _dna.phy
inputfile2=open(indel_phy,'r') # _indel.phy

output=open(phy.replace('.phy','')+'.nex','w')

try:
    path=folder
except:
    path=''
    
end_dna=0
end_indel=0

test=0
for line in inputfile1:
    if test==0:
        temp=line.replace('\n','').split(' ')
        end_dna=int(temp[1])        
        break
inputfile1.close()

test=0
for line in inputfile2:
    if test==0:
        temp=line.replace('\n','').split(' ')
        end_indel=int(temp[1])        
        break
inputfile2.close()

tempo=phy.replace('.phy','')+'_dna.phy'

temp=tempo.split('/')
name1=temp[-1]
temp=indel_phy.split('/')
name2=temp[-1]

#output write
output.write('#nexus'+'\n'+'\n')
output.write('begin sets;'+'\n')
output.write('	charset part1 = '+str(path)+name1+': '+str(1)+'-'+str(end_dna)+';'+'\n')
output.write('	charset part2 = '+str(path)+name2+': '+str(1)+'-'+str(end_indel)+';'+'\n'+'\n')
#output.write('	charpartition mine = GTR+G:part1, MK+ASC:part2;'+'\n')
output.write('end;'+'\n')
output.close()

if deletion!='no':

    os.remove(name1.replace('_dna.phy','.part'))
    os.remove(name1.replace('_dna.phy','.phy'))

    os.remove(name2.replace('_indel.phy','.part'))
    os.remove(name2.replace('_indel.phy','.phy'))
