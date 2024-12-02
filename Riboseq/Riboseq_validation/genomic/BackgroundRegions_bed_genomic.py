##############################################################################
### BackgroundRegions_adapted_genomic.py takes a bedfile and a fasta as input and
### generates random regions of the fasta file with the length distribution as 
### given in the bed file and returns them in bed format
### for the genmic regions the selection of the random regions is made by the
### name, not by the chromsome (the chromsome was the name for transcriptomic
### alignment)
##############################################################################

#usage: python BackgroundRegions_adapted.py in.bed in.fasta out.bed
import sys
import random
from pybedtools import BedTool


def get_length_dist(bed_file):
    #read in bedfile and obtain length distribution
    lengthdistribution=[]
    for line in bed_file:
        elems = line.split('\t')
        length=int(elems[2])-int(elems[1])
        lengthdistribution.append(length)
    return lengthdistribution


def get_random_bed_ids(lengthdistribution, UTR_bed_file):
    #select random keys for the fasta file, as many as there are 
        #entries in the bed file
    #do not allow for duplicate keys
    three_prime_UTRs = BedTool(UTR_bed_file)
    bed_names = [bed.name for bed in three_prime_UTRs]
    randomlist=[]
    for i in range(len(lengthdistribution)):
        random_key = random.choices(bed_names, k = 1)[0]
        #if key has already been sampled, randomly choose a new key
        while random_key in randomlist:
            random_key = random.choices(bed_names, k = 1)[0]
        randomlist.append(random_key)
    return randomlist, three_prime_UTRs, bed_names


def write_bed_output(outname, randomlist, three_prime_UTRs, lengthdistribution, bed_names):
    #for each length from the length distribution
        #get a random sample of the same length from the fasta sequence
    random_list_reference = randomlist.copy()
    with open(outname, 'w') as out:
        i = len(randomlist) - 1
        while i > 0:
            length = lengthdistribution.pop()
            random_key = randomlist.pop()
            random_seq = three_prime_UTRs.filter(lambda interval: interval.name == random_key)[0]
            UTR_length = random_seq.end - random_seq.start
            if length < UTR_length:
                start = random.choice(range(random_seq.start, UTR_length - length + random_seq.start))
                end = start + length
                if random_seq.strand == '1':
                    out.write(random_seq.chrom + '\t' + str(start) + '\t' + str(end) + '\t' + str(random_key)
                    + '\t' + str(0) + '\t' + '+' + '\n')
                elif random_seq.strand == '-1':
                    out.write(random_seq.chrom + '\t' + str(start) + '\t' + str(end) + '\t' + str(random_key)
                    + '\t' + str(0) + '\t' + '-' + '\n')
                i -= 1 
            else:
                #if the sampled fasta sequence is too short, sample a new random fasta entry
                #and reappend the length to resample
                random_key = random.choices(bed_ids, k = 1)[0]
                while random_key in random_list_reference:
                    random_key = random.choices(bed_ids, k = 1)[0]
                randomlist.append(random_key)
                random_list_reference.append(random_key)
                lengthdistribution.append(length)



if len(sys.argv) < 3:
        print('usage python BackgroundRegions_adapted_genomic.py in.bed UTR.bed out.bed')
else :
    file = open(sys.argv[1],'r')
    lengthdistribution = get_length_dist(file)  
    
    randomlist, three_prime_UTRs, bed_ids = get_random_bed_ids(lengthdistribution, sys.argv[2])
    
    write_bed_output(sys.argv[3], randomlist, three_prime_UTRs, lengthdistribution, bed_ids)
    
        
