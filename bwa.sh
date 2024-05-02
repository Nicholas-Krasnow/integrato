indDir=$(pwd)/..
mkdir bams_human
mkdir bams_pmid

dir='trimmedFastqs'
prefix="UMItrim_"

for sample in "$dir/$prefix"*; do

	if [ -f "$sample" ]; then
		name="$(cut -d'.' -f1 <<<"$sample")"
		
		bwa mem "$indDir"/genomeref/GRCh38.p14.fna ${sample} >"${name}".bam
		mv "${name}".bam bams_human

		bwa mem "$indDir"/plasmidref/CCR5donorPlasmid.fasta ${sample} >"${name}".bam
		mv "${name}".bam bams_pmid

	fi

done