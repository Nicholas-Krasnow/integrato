mkdir trimmedFastqs

adapterDir=$(dirname "$(pwd)")
adapter=adapterSeqs.fasta
adapterFile="$adapterDir/$adapter"
echo "$adapterFile"

cd fastqs

for sample in $(ls *.fastq); do
	prefix="trim_"
	name=$(basename "$sample")
	fastp -i "$sample" -o "trim_$name" --adapter_fasta "$adapterFile"
	mv "$prefix$name" ../trimmedFastqs
	
done

cd ..