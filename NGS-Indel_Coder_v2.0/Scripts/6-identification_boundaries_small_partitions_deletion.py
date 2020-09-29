# -*- coding: us-ascii -*-

### Usage ###

#By Julien Boutte, Septembe 2020
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
opts.add_option('-t',type='int',dest='threshold',default='100')
opts.add_option('-o',type='str',dest='outputfile_fasta_aligned',default='')
opts.add_option('-d',type='str',dest='deletion',default='no')

(options, args)=opts.parse_args(sys.argv[1:])

blastn=options.blastn
fasta=options.fasta
threshold=options.threshold
outputfile_fasta_aligned=options.outputfile_fasta_aligned
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

# PART 2

#1- Inputfile 1
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
    
lg=len(list_seq[0])

#2- Identification list_intron and list_exon
        
#check start

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
if end<lg:
    temp=[]
    temp.append(end+1)
    temp.append(lg)
    list_intron.append(temp)

#3- identify small introns and small exons

list_tot=[]
ordo=[]
for i in list_exon:
    if int(i[1])<int(i[0]): plop=1
    else:
        temp=[]
        temp.append(int(i[0]))
        temp.append(int(i[1]))
        ordo.append(int(i[0]))
        temp.append('E')
        list_tot.append(temp)
    
for i in list_intron:
    if int(i[1])<int(i[0]): plop=1
    else:
        temp=[]
        temp.append(int(i[0]))
        temp.append(int(i[1]))
        ordo.append(int(i[0]))
        temp.append('I')
        list_tot.append(temp)

#Ordonate list_tot
ordo_ordo=sorted(ordo)

list_tot_ordo=[]
for i in ordo_ordo:
    for j in list_tot:
        if j[0]==i:
            list_tot_ordo.append(j)

if len(list_tot)!=len(list_tot_ordo):
    print 'error length, exit'
    print list_tot
    print list_tot_ordo
    exit()

#delete small regions in list_tot_ordo
list_delete_fasta=[] #for fasta file

save=[]
for i in list_tot_ordo:
    save.append(i)

length=int(threshold)

i=0
while i!=len(list_tot_ordo):
    if list_tot_ordo[i][1]-list_tot_ordo[i][0]+1<length:
        #for fasta file
        temp=[]
        temp.append(list_tot_ordo[i][0]-1) ##### ATTENTION -1 #####
        temp.append(list_tot_ordo[i][1]-1) ##### ATTENTION -1 #####
        list_delete_fasta.append(temp) 
        #change values
        supp=list_tot_ordo[i][1]-list_tot_ordo[i][0]+1
        j=i+1
        while j<len(list_tot_ordo):
            temp=[]
            temp.append(list_tot_ordo[j][0]-supp)
            temp.append(list_tot_ordo[j][1]-supp)
            temp.append(list_tot_ordo[j][2])
            list_tot_ordo[j]=temp
            j=j+1
        del(list_tot_ordo[i])
    else:
        i=i+1

#create list_final_exon containing new positions of exons
list_final_exon=[]

for i in list_tot_ordo:
    if i[2]=='E':
        temp=[]
        temp.append(i[0])
        temp.append(i[1])
        list_final_exon.append(temp)

#Using list_delete_fasta, list_name and list_seq
#create a new fasta file

list_retained=[]

if len(list_delete_fasta)>0:
    if list_delete_fasta[0][0]!=0:
        temp=[]
        temp.append(0)
        temp.append(list_delete_fasta[0][0]-1)
        list_retained.append(temp)
    i=0
    while i!=len(list_delete_fasta)-1:
        temp=[]
        temp.append(list_delete_fasta[i][1]+1)
        temp.append(list_delete_fasta[i+1][0]-1)
        list_retained.append(temp)
        i=i+1

    if list_delete_fasta[i][1]!=lg-1:
        temp=[]
        temp.append(list_delete_fasta[i][1]+1)
        temp.append(lg-1)
        list_retained.append(temp)

    output=open(outputfile_fasta_aligned,'w')

    i=0
    while i!=len(list_name):
        output.write('>'+list_name[i]+'\n')
        j=0
        while j!=len(list_retained):
            output.write(list_seq[i][list_retained[j][0]:list_retained[j][1]])
            j=j+1
        output.write('\n')
        i=i+1
    output.close()

else:
    
    output=open(outputfile_fasta_aligned,'w')

    i=0
    while i!=len(list_name):
        output.write('>'+list_name[i]+'\n')
        output.write(list_seq[i]+'\n')
        i=i+1
    output.close()

if deletion!='no':
    os.remove('exon_positions.txt')












