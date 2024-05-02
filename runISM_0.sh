#conda activate integrato_env

dn_core="GT" #Dinucleotide core for the recombinase, change to the particular core of your recombinase if different

mDir=$(pwd)
# sort files into by target site

sites=("AAVS1" "CCR5")
for site in "${sites[@]}"; do
    mkdir "$site"
    mkdir "$site"/fastqs
done

for sample in $(ls *.fastq.gz); do
    file_name=$(basename "$sample")
    echo "$file_name"

    for site in "${sites[@]}"; do
        if [[ "$file_name" == *"$site"* ]]; then
            mv "$file_name" "$site"/fastqs
        fi
     done
done

# iterate over the sites

#get the donor-specific primer to analyze with

primer=$(python3 "$mDir"/get_fasta.py primer.fasta)

echo "$primer"
#echo "$att_sub"

for site in "${sites[@]}"; do

    #define the att subsequence specific to the donor for this site
    site_att="$site"
    site_att+="_"
    site_att+="att_sub.fasta"
    att_sub=$(python3 "$mDir"/get_fasta.py "$site_att")

    # unzip
    cd "$site"/fastqs
    for sample in $(ls *.fastq.gz); do
        
        
        gunzip $sample
        
    done
    cd ..
    # adapter trimming and quality filtering
    sh "$mDir"/trim.sh

    #filtering by primer presence and UMI depduplication
    cd trimmedFastqs
    python3 "$mDir"/primer_donorUMI_filter.py "$primer" "$att_sub"
    cd ..

    #run BWA alignment
    sh "$mDir"/bwa.sh 

    #Bam analysis
    python3 "$mDir"/read_bams.py "$primer" "$att_sub" "$dn_core"

    cd .. #back to the main folder


done

