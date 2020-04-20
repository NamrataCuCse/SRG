from Bio import Seq, SeqIO, SeqRecord
import random 
import sys
import os

filename=sys.argv[1] #Reference
depth=int(sys.argv[2]) #depth %
#USAGE : python simulated_read_PE.py reference.fa 10


read1file = open("simulatedR_1.fastq","w")
read2file = open("simulatedR_2.fastq","w")

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'N':'N'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)
    def reverse_complement(s):
        return complement(s[::-1])

def reverse(s): 
  str = "" 
  for i in s: 
    str = i + str
  return str
k=0
for record in SeqIO.parse(filename, "fasta"):
    read=record.seq
    rd_len= len(read)
    if rd_len> 160 :
	    depth_cover= (rd_len*depth/100)/2
	    i=1
	    while(i <= depth_cover):
	       frag_start = random.randrange(0,rd_len+1,1)       
	       insert = random.randrange(160,221,1)
	       frag_end = frag_start + insert 
	       if frag_end <= rd_len :
		    print( k, frag_end , rd_len)
		    #i = i + 1
		    flag = random.randrange(0,2) 
		    if flag == 1:
		        fragment = complement(reverse(read[frag_start-1:frag_end]))
		        if 'N' not in fragment :
                           i = i + 1
                           k = k + 1  
		           read1 = fragment[:100]
		           read2 = complement(reverse(fragment[-100:]))
                           #print read1
		           read1file.write("@simulated_read"+"|-|"+str(k)+"\n"+str(read1.upper())+"\n"+str("+")+"\n"+str("I"*100)+"\n")
		           read2file.write("@simulated_read"+"|-|"+str(k)+"\n"+str(read2.upper())+"\n"+str("+")+"\n"+str("I"*100)+"\n")
		    else: 
			fragment = read[frag_start-1:frag_end]
		        if 'N' not in fragment :
                           i = i + 1
                           k= k + 1
		           read1 = fragment[:100]
			   read2 = complement(reverse(fragment[-100:]))
                           #print read1
		           read1file.write("@simulated_read"+"|+|"+str(k)+"\n"+str(read1.upper())+"\n"+str("+")+"\n"+str("I"*100)+"\n")
			   read2file.write("@simulated_read"+"|+|"+str(k)+"\n"+str(read2.upper())+"\n"+str("+")+"\n"+str("I"*100)+"\n")
	       

