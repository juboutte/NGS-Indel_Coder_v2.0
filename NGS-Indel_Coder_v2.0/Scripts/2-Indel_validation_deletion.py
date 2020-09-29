# -*- coding: us-ascii -*-

import sys
import os
from optparse import OptionParser

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

def mean(liste):
    total=0
    for i in liste:
        total=total+float(i)
    return (float(float(total)/len(liste)))
    
#start program

opts=OptionParser()
opts.add_option('-f',type='str',dest='fasta_file',default='')
opts.add_option('-r',type='str',dest='read_path',default='')
opts.add_option('-t',type='int',dest='threshold1',default='0')
opts.add_option('-d',type='str',dest='deletion',default='no')
opts.add_option('-l',type='int',dest='length',default='2')

(options, args)=opts.parse_args(sys.argv[1:])

fasta_file=options.fasta_file
read_path=options.read_path
threshold1=options.threshold1
deletion=options.deletion
length=options.length

print 'fasta file used:',fasta_file
print 'read depth path used:',read_path
print 'threshold selected:',threshold1
print 'indel length selected:',length

#Part 1

#Parsing of the fasta file
input_file=open(fasta_file,'r')

liste_name=[]
liste_seq=[]

seq_temp=''

for line in input_file:
    if line[0]!='>':
        temp=line.split()
        if len(temp)>0:
            seq_temp=seq_temp+temp[0]
    else:
        liste_name.append(line[1:].replace('\n','').replace('\r',''))
        if seq_temp!='':
            liste_seq.append(seq_temp)
            seq_temp=''

liste_seq.append(seq_temp) 

liste_deb=[]
liste_fin=[]

#print liste_name
#print liste_seq

#beginning and end of sequences
i=0
while i!=len(liste_seq):
    j=0
    while liste_seq[i][j]=='-' or liste_seq[i][j:j+9].count('-')!=0:
        j=j+1
    liste_deb.append(j)
    j=len(liste_seq[0])-1
    if j!='-':
        liste_fin.append(j+1)
    else:
        while liste_seq[i][j]=='-' or liste_seq[i][j:j-9].count('-')!=0:
            j=j-1
        liste_fin.append(j+1)
    i=i+1

#print liste_deb
#print liste_fin

nb_SNPs_un=0
nb_SNPs_clean=0

liste_SNPs_clean=[]


if len(liste_seq)==0:
   print "error with this sample:", fasta_file
   print "script exit"
   exit()
else:
    i=0
    while i!=len(liste_seq[0]): #each position
        temp=[]
        j=0
        while j!=len(liste_seq): #each sample
            if liste_deb[j]<=i and liste_fin[j]>=i:
                #print liste_name[j], i
                temp.append(liste_seq[j][i]) 
            j=j+1
        #print temp
        if len(list(set(temp)))>1: #SNPs
            proba=list(set(temp))
            temp2=[]
            for n in proba:
                if temp.count(n)>1:
                    temp2.append(n)
            if len(temp2)>1:
                nb_SNPs_clean=nb_SNPs_clean+1 #SNPs with nucleotide present in more than 1 seq
                if temp.count('-')>=1:
                    liste_SNPs_clean.append(i)
            else:
                nb_SNPs_un=nb_SNPs_un+1 # SNP specific to one sequence
        i=i+1

    #print liste_SNPs_clean

    mini=len(liste_seq[0])
    maxi=0

    i=1
    boundaries=[]
    while i!=len(liste_SNPs_clean):
        deb=liste_SNPs_clean[i-1]
        while liste_SNPs_clean[i]-1==liste_SNPs_clean[i-1] and i!=len(liste_SNPs_clean)-1:
            i=i+1
        fin=liste_SNPs_clean[i-1]
        temp=[]
        if fin-deb>=length: #>1 try >0 to see if we include indel=1 or not #old:if fin-deb>1:
            temp.append(deb)
            temp.append(fin)
            boundaries.append(temp)
            if fin-deb+1<mini:
                mini=fin-deb+1
            if fin-deb+1>maxi:
                maxi=fin-deb+1
        #print i, liste_SNPs_clean[i], liste_SNPs_clean[i-1]
        i=i+1
    #print boundaries

    #PART 2: look read depth
    output=open('temp_part1.txt','w')
    i=0
    while i!=len(boundaries): #for each insertion/deletion
        output.write('>'+fasta_file+'\t'+str(boundaries[i][0])+'\t'+str(boundaries[i][1])+'\n')
        j=0
        while j!=len(liste_seq): #each sample
            if liste_deb[j]<=boundaries[i][0] and liste_fin[j]>=boundaries[i][1]: #the sequence present the insertion/deletion
                #open read depth speicific to the sequence
                file_depth=open(str(read_path)+str(liste_name[j])+'.txt','r')
                liste_depth=[]
                h=0
                for line in file_depth:
                    #print liste_name[i], h, len(liste_seq[i])
                    temp=line.replace('\n','').split('\t')
                    if liste_seq[j][h]!='N' and liste_seq[j][h]!='n':
                        liste_depth.append(temp[2])
                    else:
                        liste_depth.append('NA')
                        h=h+1
                file_depth.close()

                ###
                
                #identification of the real position of begining and end of the insertion/deletion
                pos=0
                k=0
                while k!=boundaries[i][0]: #identification of the real start position pos
                    if liste_seq[j][k]!='-':
                        pos=pos+1
                    k=k+1
                    pos2=0
                k=0
                while k!=boundaries[i][1]+1: #identification of the real end position pos2
                    if liste_seq[j][k]!='-':
                        pos2=pos2+1
                    k=k+1                
                output.write(liste_name[j]+'\t'+str(pos)+'\t'+str(pos2)+'\t'+str(liste_depth[pos:pos2])+'\t') #test to validate position
                #next step: indicate N or 0 if deletion and start 10 bp before and end 10 bp after
                if pos-10>0 and pos2+10<len(liste_seq[j].replace('-','')):
                    output.write(str(pos-10)+'\t'+str(pos2+10)+'\t'+str(liste_depth[pos-10:pos2+10])+'\n') 
                elif pos-10>0:
                    #output.write(str(pos-10)+'\t'+str(len(liste_seq[j].replace('-','')))+'\t'+str(liste_depth[pos-10:len(liste_seq[j].replace('-',''))])+'\n')
                    output.write(str(pos-10)+'\t'+str(len(liste_seq[j].replace('-',''))-1)+'\t'+str(liste_depth[pos-10:])+'\n')
                else:
                    output.write(str(0)+'\t'+str(pos2+10)+'\t'+str(liste_depth[0:pos2+10])+'\n')
                ###
                
            j=j+1
        i=i+1
    output.close()

#Part 2

inputfile=open('temp_part1.txt','r')
output=open('temp_part2.txt','w')

pb=0

for line in inputfile:
    if line[0]!='>':
        temp=line.replace('\n','').split('\t')
        if len(temp)==7:
            output.write(temp[0]+'\t')
            if temp[3]=='[]':
                output.write('GAP'+'\t')
                output.write('NA'+'\t'+'NA'+'\t'+'NA'+'\t')
            else:
                output.write('No-GAP'+'\t')
                depth=temp[3].replace('[','')
                depth=depth.replace(']','')
                depth=depth.replace("'",'')
                depth=depth.split(',')
                liste_depth=[]
                for i in depth:
                    liste_depth.append(int(i))
                output.write(str(min(liste_depth))+'\t'+str(max(liste_depth))+'\t'+str(mean(liste_depth))+'\t')
            ###
            if temp[6]=='['']' or temp[6]=='[]':
                output.write(str('NA')+'\t'+str('NA')+'\t'+str('NA')+'\n')
                pb=pb+1
                #print line
                #print temporaire
            else:
                depth=temp[6].replace('[','')
                depth=depth.replace(']','')
                depth=depth.replace("'",'')
                depth=depth.split(',')
                liste_depth=[]
                for i in depth:
                    #print i
                    if 'NA' not in i:
                        liste_depth.append(int(i))
                output.write(str(min(liste_depth))+'\t'+str(max(liste_depth))+'\t'+str(mean(liste_depth))+'\n')
        else:
            print 'error'
            print line
            print 'program exit'
            exit()
        ###
    else:
        output.write(line)

output.close()           

print "Number of insertion/deletion per sample for which not information were available:",pb
print "a low number of this value is not a serious problem (example generated 1)"

#Part 3

threshold=int(threshold1)

print "threshold selected:", threshold

#Creation of 1 output files

inputfile=open('temp_part2.txt','r')

#D1= region of GAP
#D2= region of GAP + 10bp before+10 bp after

liste_mean_D2_GAP=[]
liste_mean_D2_NoGAP=[]

output4=open('MRD2_MRD3_T'+str(threshold)+'.txt','w')

for line in inputfile:
    if line[0]=='>':
        if len(liste_mean_D2_GAP)!=0 and len(liste_mean_D2_NoGAP)!=0:
            #print 'step1'
            liste_mean_D2_GAP=sorted(liste_mean_D2_GAP)
            liste_mean_D2_NoGAP=sorted(liste_mean_D2_NoGAP)
            if max(liste_mean_D2_NoGAP)>=max(liste_mean_D2_GAP):
                toto=max(liste_mean_D2_NoGAP)
                plop=0
                for i in liste_mean_D2_NoGAP:
                    if abs(max(liste_mean_D2_GAP)-i)<toto:
                        toto=abs(max(liste_mean_D2_GAP)-i)
                        tata=i
                        plop=1
                if plop==1:
                    plup=tata
                else:
                    plup=max(liste_mean_D2_NoGAP)
                ### HERE
                if plup<=threshold or max(liste_mean_D2_GAP)<=threshold:
                    output4.write(name+'\t'+str(max(liste_mean_D2_GAP))+'\t'+str(plup)+'\n')
            else:
                toto=max(liste_mean_D2_GAP)
                plop=0
                for i in liste_mean_D2_GAP:
                    if abs(max(liste_mean_D2_NoGAP)-i)<toto:
                        toto=abs(max(liste_mean_D2_NoGAP)-i)
                        tata=i
                        plop=1
                if plop==1:
                    plup=tata
                else:
                    plup=max(liste_mean_D2_GAP)
                ### HERE
                if plup<=threshold or max(liste_mean_D2_NoGAP)<=threshold:
                    output4.write(name+'\t'+str(plup)+'\t'+str(max(liste_mean_D2_NoGAP))+'\n')
            liste_mean_D2_GAP=[]
            liste_mean_D2_NoGAP=[]
        name=line[1:].replace('\n','')
    else:
        temp=line.split('\t')
        if temp[-1].replace('\n','')!='NA':
            if temp[1]=='No-GAP' and temp[7]!='NA':
                #print 'step2'
                liste_mean_D2_NoGAP.append(float(temp[7]))
            elif temp[7]!='NA':
                #print 'step3'
                liste_mean_D2_GAP.append(float(temp[7].replace('\n','')))

if len(liste_mean_D2_GAP)!=0 and len(liste_mean_D2_NoGAP)!=0:
    liste_mean_D2_GAP=sorted(liste_mean_D2_GAP)
    liste_mean_D2_NoGAP=sorted(liste_mean_D2_NoGAP)
    if max(liste_mean_D2_NoGAP)>=max(liste_mean_D2_GAP):
        toto=max(liste_mean_D2_NoGAP)
        plop=0
        for i in liste_mean_D2_NoGAP:
            if abs(max(liste_mean_D2_GAP)-i)<toto:
                toto=abs(max(liste_mean_D2_GAP)-i)
                tata=i
                plop=1
        if plop==1:
            plup=tata
        else:
            plup=max(liste_mean_D2_NoGAP)
        ### HERE
        if plup<=threshold or max(liste_mean_D2_GAP)<=threshold:
            output4.write(name+'\t'+str(max(liste_mean_D2_GAP))+'\t'+str(plup)+'\n')
    else:
        toto=max(liste_mean_D2_GAP)
        plop=0
        for i in liste_mean_D2_GAP:
            if abs(max(liste_mean_D2_NoGAP)-i)<toto:
                toto=abs(max(liste_mean_D2_NoGAP)-i)
                tata=i
                plop=1
        if plop==1:
            plup=tata
        else:
            plup=max(liste_mean_D2_GAP)
        ### HERE
        if plup<=threshold or max(liste_mean_D2_NoGAP)<=threshold:
            output4.write(name+'\t'+str(plup)+'\t'+str(max(liste_mean_D2_NoGAP))+'\n')
output4.close()

#Part 4

#Function to ordonate windows for each alignment fasta file
def tri_dico(liste):
    if len(liste)==1:
        return liste
    else:
        tempo=[]
        for i in liste:
            tempo.append(i[0])
        order=sorted(tempo, reverse=True) ###
        if len(list(sorted(order)))!=len(liste):
            print 'ERROR FUNCTION 1 - Part 4'
        final=[]
        i=0
        while i!=len(order):
            for j in liste:
                if j[0]==order[i]:
                    final.append(j)
            i=i+1
        if len(final)==len(liste):
            return final
        else:
            print 'ERROR FUNCTION 2 - Part 4'
#Part 1: extraction in .txt file of the different windows
       
inputfile=open('MRD2_MRD3_T'+str(threshold)+'.txt','r')
dico_seq={}

for line in inputfile:
    temp=line.split('\t')
    name=temp[0].replace('\n','').replace('\r','').split('/')
    if name[-1] not in dico_seq:
        dico_seq[name[-1]]=[]
    liste_tempo=[]
    liste_tempo.append(int(temp[1]))
    liste_tempo.append(int(temp[2]))
    dico_seq[name[-1]].append(liste_tempo)
    liste_tempo=[]

inputfile.close()

#print len(dico_seq)
#print dico_seq

#Part 2: for each element of dico_seq:
#open fasta file and modify it

#Warning: necessary to ordonate windows by descendant values

#print dico_seq

for i in dico_seq:
    i2=i.replace('\n','').replace('\r','').split('/')
    #organise dico_seq[i]
    tempo=[]
    for j in dico_seq[i]:
        tempo.append(j)
    tempo2=tri_dico(tempo)
    #tempo2 is the descending order of dico_seq
    output=open(i.replace('.fasta','_temp.fasta'),'w')
    liste_name=[]
    liste_seq=[]
    seq_temp=''
    #inputfile=open(str(fasta_file)+str(i2[-1]),'r')
    inputfile=open(str(fasta_file),'r')
    for line in inputfile:
        if line[0]!='>':
            temp=line.split()
            if len(temp)>0:
                seq_temp=seq_temp+temp[0]
        else:
            liste_name.append(line[1:].replace('\n','').replace('\r',''))
            if seq_temp!='':
                liste_seq.append(seq_temp)
                seq_temp=''

    liste_seq.append(seq_temp)    
    #outputfile
    j=0
    while j!=len(liste_seq):
        seq_temp=liste_seq[j]
        for k in tempo2:
            seq_temp=seq_temp[:k[0]]+seq_temp[k[1]+1:]
        output.write('>'+liste_name[j]+'\n'+seq_temp+'\n')
        seq_temp=''
        j=j+1
    output.close()

if deletion!='no':
    os.remove('temp_part1.txt')
    os.remove('temp_part2.txt')
    os.remove('MRD2_MRD3_T'+str(threshold)+'.txt')

