#start in trimmedFastqs directory

from Bio import SeqIO
from Bio.Seq import Seq
import os
import subprocess
import sys
import pandas as pd
from collections import Counter

# read in fastqs and make dictionary
recordsDict = {}
for filename in os.listdir():
    if ".fastq" in filename:
        record = list(SeqIO.parse(filename,'fastq'))
        recordsDict[filename] = record

# Filter for correct primer sequence and plasmid att subsequence
print('filtering by primer...')

primer = str(sys.argv[1]) #'GCTCCTCGCCCTTGCTCACG' # requires the primer sequence as an argument when running the script. Change to sys.argv in final code
att_sub = str(sys.argv[2]) #'CTCCGTCGTCAGGATCAT'# requires the primer sequence as an argument when running the script. Change to sys.argv in final code

primerLen = len(primer)
primerRev = Seq(primer).reverse_complement()
primerRead = str(primerRev) # when the primer appears reverse complemented in reads

#make a dictionary with keys as filenames and values as lists: correct primer reads, wrong primer reads, total reads
recordStats = {}
correctPrimerDict = {}
for records in recordsDict:
    
    correctPrimerReads = []
    wrongPrimer = 0
    for x in recordsDict[records]:
        sequence = str(x.seq)
        if primerRead in sequence and att_sub in sequence: #require primer and att
            #remove everything after the primer in the read
            primer_index = sequence.find(primerRead)
            end_bound = primer_index + len(primerRead)

            new_entry = SeqIO.SeqRecord(Seq(sequence[:end_bound]), id=x.id, description=x.description)
            new_entry.letter_annotations["phred_quality"] = [40] * len(sequence[:end_bound])
            correctPrimerReads.append(new_entry) 
        else:
            wrongPrimer += 1
    
    correctPrimerDict[records]=correctPrimerReads
    recordStats[records]=[len(correctPrimerReads),wrongPrimer,len(recordsDict[records])]

#filter for unique records
uniqueDict = {}
umi_unique_dict = {}
umiReport = {}

print('analyzing UMIs...')

for y in correctPrimerDict:
    
    #uniqueRecords = {}
    umiDict = {}
    #umi_unique_records = {}
    dups = 0
    umi_dups = 0

    
    for x in correctPrimerDict[y]:
        sequence = str(x.seq)
        
        #get the UMI
        primer_index = sequence.find(primerRead)
        umi_end = primer_index - 1
        umi_start = primer_index - 10
        umi = sequence[umi_start:umi_end+1]
        
        #check the umi and its context
        #print(umi)
        
        #check duplicates even though filtering by UMI just to see how they compare
       # if sequence in uniqueRecords:
        #    dups += 1
        #else:
         #   uniqueRecords[sequence] = x
        
        #count the umi if present or create new umi entry and keep the read if new umi
        #use umiDict to hold list of all sequences with the UMI 
        if umi in umiDict:
            umiDict[umi].append(x)
            umi_dups += 1
        else: 
            umiDict[umi] = [x]
            #umi_unique_records[sequence] = x

    #make a dictionary for the sequences for each UMI (like umiDict but sequences instead of full fastq entry)
    bins = {}
    for umi in list(umiDict):
        sequences = []
        for entry in umiDict[umi]:
            sequences.append(str(entry.seq))
        bins[umi] = sequences
    
    
    
    #get most common sequence and add corresponding fastq entry to the final dictionary
    #don't use this filter this time
    seq_artifacts = []
    #seq_artifacts = ['GGCTTGTCGACGACGGCGGTCTCCGTCGTCAGGATCAT','AGATAGAACCGCGGCCCCCCACCGCCAGGT','GGTTTGTCTGGTCAACCACCGCGGTCTCAGTGGTGTACGGTACAAACC','TAGATAGAACCGCGGATCACTTGGG','ATCAACTTGAAAAAGTGGCACCGAGTCG','GCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGT',
                     #'TAGTTGGTTTAACGCGTAACTAGATAGAACCGCG','ACGCGTAACTAGATAGAACCGCG']

    umi_unique_records = {}
    for binned_umi in list(bins):

        #check if the bin has any artifacts
        has_artifacts = False
        for umi_ex in bins[binned_umi]:
            for artifact in seq_artifacts:
                if artifact in umi_ex:
                    has_artifacts = True
                    break

        most_common_seq = Counter(bins[binned_umi]).most_common(1)[0][0]

        if has_artifacts == False:
            for x in umiDict[binned_umi]:
                

                if str(x.seq) == most_common_seq: #requires that no artifacts were found to retain the umi
                    umi_unique_records[binned_umi] = x
                    break
                else:
                    pass
            
    #uniqueDict[y] = uniqueRecords
    umi_unique_dict[y] = umi_unique_records
    
    #fill the dictionary for summary statistics for each sample: unique reads, umi duplicates, total reads from primer filtering
    #umiReport[y] = [len(umi_unique_records), umi_dups, len(correctPrimerDict[y])]

#Write the filtered fastq reads to a file that will be used for alignment
   
for name in umi_unique_dict:
    prefix = 'UMI'
    SeqIO.write(umi_unique_dict[name].values(),prefix + name,'fastq')

# write summary stats to csv
#statsDF = pd.DataFrame.from_dict(umiReport)
#cols = ['UMI-unique reads','UMI-duplicate reads','input reads with primer and att present']
#statsDF=statsDF.T
#statsDF.columns=cols
#statsDF.to_csv('UMI_Stats.csv')




