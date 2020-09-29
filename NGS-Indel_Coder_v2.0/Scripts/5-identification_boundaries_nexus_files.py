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
#Goal of this program:

import sys
import os
from optparse import OptionParser

opts=OptionParser()
opts.add_option('-b',type='str',dest='blastn',default='')
opts.add_option('-f',type='str',dest='fasta',default='')
opts.add_option('-n',type='str',dest='nexus',default='')
opts.add_option('-d',type='str',dest='deletion',default='no')

(options, args)=opts.parse_args(sys.argv[1:])

blastn=options.blastn
fasta=options.fasta
nexus=options.nexus
deletion=options.deletion

#PART 1

inputfile=open(blastn,'r')
output=open('exon_positions.txt','w')

list_introns=[]

for line in inputfile:
    if line[0]!='#':
        temp=line.split('\t')
        tempo=[]
        if int(temp[6])<int(temp[7]):
            tempo.append(int(temp[6]))
            tempo.append(int(temp[7]))
        else:
            tempo.append(int(temp[7]))
            tempo.append(int(temp[6]))
        list_introns.append(tempo)

#print list_introns

test=0
while test==0:
    test=1
    i=0
    while i<=len(list_introns)-2:
        j=0
        while j<=len(list_introns)-1:
            if i!=j:
                #print i,j, len(list_introns)
                #print list_introns[i],list_introns[j]
                if list_introns[i][0]<=list_introns[j][0] and list_introns[i][1]>=list_introns[j][1]:
                    del(list_introns[j])
                    test=0
                elif list_introns[j][0]<=list_introns[i][0] and list_introns[j][1]>=list_introns[i][1]:
                    del(list_introns[i])
                    test=0
                    j=0
                elif list_introns[i][0]<=list_introns[j][0] and list_introns[j][0]<=list_introns[i][1]:
                    tempo=[]
                    tempo.append(min(list_introns[i][0],list_introns[j][0]))
                    tempo.append(max(list_introns[i][1],list_introns[j][1]))
                    #print list_introns[i], list_introns[j], tempo
                    if i<j:
                        del(list_introns[j])
                        del(list_introns[i])
                    else:
                        del(list_introns[i])
                        del(list_introns[j])
                    j=0
                    i=0
                    list_introns.append(tempo)
                    test=0
                elif list_introns[j][0]<=list_introns[i][0] and list_introns[i][0]<=list_introns[j][1]:
                    tempo=[]
                    tempo.append(min(list_introns[i][0],list_introns[j][0]))
                    tempo.append(max(list_introns[i][1],list_introns[j][1]))
                    #print list_introns[i], list_introns[j], tempo
                    if i<j:
                        del(list_introns[j])
                        del(list_introns[i])
                    else:
                        del(list_introns[i])
                        del(list_introns[j])
                    j=0
                    i=0
                    list_introns.append(tempo)
                    test=0
                else:
                    j=j+1
            else: j=j+1
        i=i+1


#reorder list_introns
list_temp=[]
for i in list_introns:
    list_temp.append(i[0])
list_temp.sort()

#print list_temp

list_introns_finale=[]

i=0
while i!=len(list_temp):
    j=0
    while list_introns[j][0]!=list_temp[i]:
        j=j+1
    list_introns_finale.append(list_introns[j])
    i=i+1

temp=fasta.split('/')
temp=temp[-1].split('.')

output.write(str(temp[0])+'\n')
output.write(str(list_introns_finale)+'\n')
output.close()

#PART 2

#Inputfile 1: exons_positions.txt

input1=open('exon_positions.txt',"r")
dico_exon={}

test=0
for line in input1:
    if test==0:
        test=1
        temp=line.replace('\n','')
        temp=temp.replace('.txt','')
        dico_exon[temp]=[]
        name1=temp
    else:
        test=0
        temp2=line.replace('\n','')
        temp2=temp2.replace('[','')
        temp2=temp2.replace(']','')
        temp2=temp2.replace(' ','')
        temp2=temp2.split(',')
        list_pos=[]
        plop=[]
        j=0
        test2=0
        while j!=len(temp2):
            if len(plop)==2:
                list_pos.append(plop)
                plop=[]
            plop.append(int(temp2[j]))
            j=j+1
        list_pos.append(plop)
        dico_exon[temp]=list_pos
        
input1.close()

#for i in dico_exon:
#    print i, dico_exon[i]
#exit()

#Open each fasta, and identify exon and intron

input_file=open(fasta,'r')
list_name=[]
list_seq=[]

seq_temp=''

for line in input_file:
    if line[0]!='>':
        temp=line.split()
        if len(temp)>0:
            seq_temp=seq_temp+temp[0]
    else:
        list_name.append(line[1:].replace('\n','').replace('\r',''))
        if seq_temp!='':
            list_seq.append(seq_temp)
            seq_temp=''

list_seq.append(seq_temp) 

list_intron=[]
list_exon=[]
if dico_exon[name1][0][0]!=1:
    temp=[]
    temp.append(1)
    temp.append(dico_exon[name1][0][0]-1)
    list_intron.append(temp)
temp2=[]
j=0
while j!=len(dico_exon[name1]):
    temp=[]
    temp.append(dico_exon[name1][j][0])
    temp.append(dico_exon[name1][j][1])
    list_exon.append(temp)
    j=j+1
    
end=dico_exon[name1][-1][1]

if len(dico_exon[name1])>1:
    j=0
    tempo=[]
    while j!=len(dico_exon[name1]):
        tempo.append(dico_exon[name1][j][0])
        tempo.append(dico_exon[name1][j][1])
        j=j+1
    del(tempo[-1])
    del(tempo[0])
    j=0
    while j!=len(tempo):
        temp=[]
        temp.append(tempo[j]+1)
        temp.append(tempo[j+1]-1)
        list_intron.append(temp)
        j=j+2
if end<len(list_seq[0]):
    temp=[]
    temp.append(end+1)
    temp.append(len(list_seq[0]))
    list_intron.append(temp)

###

input3=open(nexus,"r")
temp=''
for line in input3:
    if "part1" in line and "charset" in line:
        temp=line.split(': ')
        #print line
#print temp
input3.close()

output=open(name1+"_partition.nex","w")
output.write("#nexus"+'\n'+'\n')
output.write("begin sets;"+'\n')
output.write(temp[0]+': ')
#DNA-Exon
j=0
while j!=len(list_exon)-1:
    output.write(str(list_exon[j][0])+'-'+str(list_exon[j][1])+'\\3,')
    output.write(str(list_exon[j][0]+1)+'-'+str(list_exon[j][1])+'\\3,')
    output.write(str(list_exon[j][0]+2)+'-'+str(list_exon[j][1])+'\\3,')
    j=j+1
output.write(str(list_exon[j][0])+'-'+str(list_exon[j][1])+'\\3,')
output.write(str(list_exon[j][0]+1)+'-'+str(list_exon[j][1])+'\\3,')
output.write(str(list_exon[j][0]+2)+'-'+str(list_exon[j][1])+'\\3;'+'\n')
#DNA-intron
intron=0
if len(list_intron)>0:
    output.write(temp[0].replace("part1","part2")+': ')
    j=0
    while j!=len(list_intron)-1:
        output.write(str(list_intron[j][0])+'-'+str(list_intron[j][1])+',')
        j=j+1
    output.write(str(list_intron[j][0])+'-'+str(list_intron[j][1])+';'+'\n')
    intron=1
#Binary
input3=open(nexus,"r")
temp=''
for line in input3:
    if "indel" in line and "charset" in line:
        temp=line.split(': ')
        #print line
#print temp
input3.close()
if len(temp)==2:
    if intron==1:
        output.write(temp[0].replace("part2","part3")+': ')
    else:
        output.write(temp[0]+': ')
    output.write(temp[1].replace('\n','')+'\n'+'end;'+'\n')
output.close()
#exit()
#create output file with name of output files created

if deletion!='no':
    os.remove('exon_positions.txt')
    
