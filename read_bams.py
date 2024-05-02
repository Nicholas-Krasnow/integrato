import os
import pysam
import pandas as pd
import subprocess
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2

#Start with bam file alignments to the human genome in one folder and alignments to the plasmid in another folder. 
#Alignments from the same fastq will have the same name in the two different directories. For each bam, make a list of the alignment names in the plasmid bam. 
#For each alignment in the genome bam, if the name is mapped in the plasmid bam list, write the alignment to a new outfile (filtered). 
#Else, add 1 to the count of reads unmapped.

# assign directory
#initialize list to hold the dictionary for each bam
# iterate over files in that directory

cwd = os.getcwd()
dir1='bams_human'
dir2 = 'bams_Pmid'

humanDicts = {}

for filename in os.listdir(dir1):
    print(filename)

    #import the fastq file that has the reads aligned in this bam, make into dictionary
    
    name = 'trimmedFastqs/' + filename[:-4] + '.fastq'
    file_fastq = list(SeqIO.parse(name,'fastq'))
    fastq_dict = {}
    for x in list(file_fastq):
        fastq_dict[x.id] = x.seq

    #make a dictionary to hold data for all alignments in the bam, with read names as keys
    dictr = {} 
    bamfile = pysam.AlignmentFile(os.path.join(dir1,filename),"rb")  

    for alignment in bamfile.fetch():
        if alignment.is_unmapped == False:
            dictr[alignment.query_name]=[alignment.is_unmapped,alignment.query_sequence,alignment.reference_name,alignment.reference_start,alignment.reference_end,str(fastq_dict[alignment.query_name])]
        else:
            dictr[alignment.query_name]=[alignment.is_unmapped]
    humanDicts[filename]=dictr
        
pmidDicts = {}

for filename in os.listdir(dir2):

    #make a dictionary to hold data for all alignments in the bam, with read names as keys
    dictPmid = {} 
    bamfile = pysam.AlignmentFile(os.path.join(dir2,filename),"rb") 

    for alignment in bamfile.fetch():
        if alignment.is_unmapped == False:
            dictPmid[alignment.query_name]=[alignment.is_unmapped,alignment.query_sequence,alignment.reference_name,alignment.reference_start,alignment.reference_end]
        else:
            dictPmid[alignment.query_name]=[alignment.is_unmapped]
    pmidDicts[filename]=dictPmid

pmidDictsmapped = {}

for x in pmidDicts:
    unmapped = 0
    mappedDictPmid = {}
    for y in pmidDicts[x]:
        if pmidDicts[x][y][0]==False:
            mappedDictPmid[y]=pmidDicts[x][y]
        else:
            unmapped += 1
            
    pmidDictsmapped[x]=mappedDictPmid
    
    print('number of reads aligned to donor = %s' %len(list(mappedDictPmid)))
    print('Unmapped = %s' %unmapped)

humanDictsmapped = {}

for x in humanDicts:
    unmapped = 0
    mappedDictHuman = {}
    for y in humanDicts[x]:
        if humanDicts[x][y][0]==False:
            mappedDictHuman[y]=humanDicts[x][y]
        else:
            unmapped += 1
            
    humanDictsmapped[x]=mappedDictHuman
    print('number of reads aligned to genome = %s' %len(list(mappedDictHuman)))
    print('Unmapped = %s' %unmapped)

#count the reads that align to plasmid only
    
mappedStats = {}
mappedToBothDict = {}
for file in list(pmidDictsmapped):
    
    mappedToBoth = {}
    unmapped = 0
    pmidOnly = {}
    for x in pmidDictsmapped[file]:
        if x in humanDictsmapped[file]:
            mappedToBoth[x]=humanDictsmapped[file][x] #keep the human genome alignment records
        else:
            unmapped += 1
            
    
    mappedToBothDict[file]=mappedToBoth
    mappedStats[file] = [len(list(humanDictsmapped[file])),unmapped,len(list(humanDicts[file]))] 

cols = ['reads mapped to genome','reads mapped to donor only','total filtered reads']

#Write the mapping stats to csv file

mappedStatsdf = pd.DataFrame.from_dict(mappedStats)
mappedStatsdf = mappedStatsdf.T
mappedStatsdf.columns = cols
mappedStatsdf.to_csv('mapping_stats.csv')

#Count the mapped genomic sites

os.chdir('../genomeref')
reference_genome = pysam.FastaFile('GRCh38.p14.fna')
genome = list(SeqIO.parse('GRCh38.p14.fna','fasta'))

chrs_dict = {}
digits = ['0','1','2','3','4','5','6','7','8','9']
for x in genome:
    index1 = x.description.find('chromosome') + len('chromosome') + 1
    index2 = x.description.find('chromosome') + len('chromosome') + 2
    if x.description[index2] in digits:
        chrs_dict[x.id] = x.description[index1:index2+1]
    else:
        if x.description[index1] in digits:
            chrs_dict[x.id] = '0' + x.description[index1]
        else:
            chrs_dict[x.id] = x.description[index1]

os.chdir(cwd)

primer = str(sys.argv[1])

primerRev = str(Seq(primer).reverse_complement())

attDonor = str(sys.argv[2])

#count the sites mapped
sitesDictDict = {}

#define the function to search the genomic locus for the putative genomic att subsequence from the read
def check_with_mismatch(subsequence, larger_sequence, max_mismatch):
    subsequence = Seq(subsequence)
    larger_sequence = Seq(larger_sequence)

    alignments = pairwise2.align.localms(subsequence, larger_sequence, 2, -1, -2, -1)
    for alignment in alignments:
        if alignment[2] >= len(subsequence) - max_mismatch:
            start_index = alignment[3]
            return True, start_index
    return False, -1

for z in humanDictsmapped:
    
    ind = z.find('_')
    concise_name = z[ind+1:]
    
    sitesDict = {}
    for x in humanDictsmapped[z]:
        query_seq = humanDictsmapped[z][x][1]
        name = str(humanDictsmapped[z][x][2])
        chr_num = 'chr'+chrs_dict[name]
        start= int(humanDictsmapped[z][x][3])
        end= int(humanDictsmapped[z][x][4])

        title = name + '_'+str(start) +'-'+ str(end)

        #check if there exists a site with the same start or end on the same chromosome
        found = False
        for y in sitesDict:
            if sitesDict[y][0] == chr_num:
                if sitesDict[y][1] ==start or sitesDict[y][2] ==end:
                    sitesDict[y][3] += 1
                    found = True
                    if sitesDict[y][6] == '': #check for att site if there isn't one found in previous reads
                        read = humanDictsmapped[z][x][-1]
                        if primer in humanDictsmapped[z][x][1]:
                            preAtt = str(Seq(reference_genome.fetch(name,humanDictsmapped[z][x][3]-50,humanDictsmapped[z][x][3]+50)).reverse_complement()).upper()
                            coord = str(read).find(attDonor)
                            att_left = str(read)[coord-25:coord]

                            coord2 = check_with_mismatch(att_left,preAtt,1)[1]

                            exact_att = str(preAtt)[coord2:coord2+48]
                            sitesDict[y][6] = exact_att
                        if primerRev in humanDictsmapped[z][x][1]:
                            preAtt = str(Seq(reference_genome.fetch(name,humanDictsmapped[z][x][4]-50,humanDictsmapped[z][x][4]+50))).upper()
                            coord = str(read).find(attDonor)
                            att_left = str(read)[coord-25:coord]

                            coord2 = check_with_mismatch(att_left,preAtt,1)[1]
                           
                            exact_att = str(preAtt)[coord2:coord2+48]
                            sitesDict[y][6] = exact_att

        #create a new site entry if none was found with the same start or end
        #extract the full attP site based on orientation of read relative to reference boundaries
        #try extracting the post-ntegration att site too
        if found == False:
            read = humanDictsmapped[z][x][-1]
            if primer in humanDictsmapped[z][x][1]:
                preAtt = str(Seq(reference_genome.fetch(name,humanDictsmapped[z][x][3]-50,humanDictsmapped[z][x][3]+50)).reverse_complement()).upper()
                coord = str(read).find(attDonor)
                att_left = str(read)[coord-25:coord]

                coord2 = check_with_mismatch(att_left,preAtt,1)[1]
                
                exact_att = str(preAtt)[coord2:coord2+48]
                
                sitesDict[title] = [chr_num,start,end,1,reference_genome.fetch(name,humanDictsmapped[z][x][3],humanDictsmapped[z][x][4]+1),preAtt,exact_att,'start',humanDictsmapped[z][x][3]-24,humanDictsmapped[z][x][3]+23,humanDictsmapped[z][x][-1]]
            if primerRev in humanDictsmapped[z][x][1]:
                preAtt = str(Seq(reference_genome.fetch(name,humanDictsmapped[z][x][4]-50,humanDictsmapped[z][x][4]+50))).upper()
                coord = str(read).find(attDonor)
                att_left = str(read)[coord-25:coord]

                coord2 = check_with_mismatch(att_left,preAtt,1)[1]
                
                exact_att = str(preAtt)[coord2:coord2+48]
                
                sitesDict[title] = [chr_num,start,end,1,reference_genome.fetch(name,humanDictsmapped[z][x][3],humanDictsmapped[z][x][4]+1),preAtt,exact_att,'end',humanDictsmapped[z][x][4]-24,humanDictsmapped[z][x][4]+23,humanDictsmapped[z][x][-1]]



    sitesDictDict[concise_name]=sitesDict


cols = ['Chr number','align start','align end','read count','reference sequence','integration site','pre-integration att', 'Integration at start or end of aligned sequence','att start coordinate','att end coordinate','example read']
cols2 = ['Chr number','align start','align end','read count','reference sequence','integration site','pre-integration att', 'Integration at start or end of aligned sequence','att start coordinate','att end coordinate','example read','nbrd']

for z in sitesDictDict:
    old_ext = z.find('.bam')
    df = pd.DataFrame.from_dict(sitesDictDict[z])

    if sitesDictDict[z] != {}:
        df = df.T
        df.columns = cols
        df = df.sort_values(by=['Chr number','align start'])
        df.reset_index(drop=True,inplace=True)
    
        #write the csv output
        df.to_csv(z[:old_ext]+'.csv')
        #write the .bed output for bedtools site clustering
        df.to_csv(z[:old_ext]+'.bed', sep="\t",index=False,header=False)
    
        #run bedtools on command line to cluster the sites into neighborhoods
        with open(z[:old_ext]+'_clustered.bed', "w") as file:
            res = subprocess.run(["bedtools", "cluster", "-d", "500", "-i", z[:old_ext]+'.bed'],stdout=file)
    
        #read in bedtools cluster output
        df_bed = pd.read_csv(z[:old_ext]+'_clustered.bed', sep='\t',header=None)
    
        dict_out = {}
        dict_out2 = {}
        
        nbrds = []
        for x in range(len(df_bed)):
            nbrd = df_bed.iloc[x][0] + str(df_bed.iloc[x][11])
            if nbrd in nbrds:
                dict_out[nbrd][3] += df_bed.iloc[x][3] # add the read count to the neihgborhood
            else:
                dict_out[nbrd] = list(df_bed.iloc[x])
                nbrds.append(nbrd)
                
        #make a new dictionary with the single-read sites removed
        dn_core = str(sys.argv[3])
        for entry in list(dict_out):
            if dict_out[entry][3] > 1:
                if dn_core in dict_out[entry][5]: #require dinucleotide core
                    dict_out2[entry] = dict_out[entry]
        
        #format dataframes for export
        
        dfout = pd.DataFrame.from_dict(dict_out)
        dfout = dfout.T
        dfout.columns = cols2
        dfout.reset_index(drop=True,inplace=True)
        

        dfout2 = pd.DataFrame.from_dict(dict_out2)
        dfout2 = dfout2.T
        
        if dict_out2 != {}: #need to check if dictionary 2 is empty
            dfout2.columns = cols2
            dfout2.reset_index(drop=True,inplace=True)
            dfout2.to_csv(z[:old_ext]+'_clustered_multiple'+'.csv')

        #write the csv output
        dfout.to_csv(z[:old_ext]+'_clustered'+'.csv')

#sort the output files
    
subprocess.run(['mkdir', 'output_clustered'])
#subprocess.run(['mkdir', 'intermediate_files'])
for filename in os.listdir():
    if '_clustered.csv' in filename:
        subprocess.run(['mv',filename,'output_clustered'])

subprocess.run(['mkdir', 'output_clustered_multiple_reads'])

for filename in os.listdir():
    if '_clustered_multiple.csv' in filename:
        subprocess.run(['mv',filename,'output_clustered_multiple_reads'])
        












