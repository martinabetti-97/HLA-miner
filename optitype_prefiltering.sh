while getopts a:i:o:r option
do
case "${option}"
in
a) bwa=${OPTARG} ;;
i) in=${OPTARG} ;;
o) out=${OPTARG} ;;
r) rm=$(echo "on") ;;
esac
done

if [ ! -d "$out" ]; then mkdir $out
fi

echo "writing output in $out"
reference_sequence=$bwa"/hla_rna.fasta"
for i in $in/*.fastq.gz;
do 
        file="$(basename $i)"
        sample=${file%.fastq.gz*}
        echo "Running filtering step on: $file"
        start=$(date +%s.%N)
        bwa mem -t 4 $reference_sequence $i > $out$sample".sam"
        echo "Alignment of $file on reference genome completed. 
        Creating new fastq file."
        samtools fastq -F4 $out$sample".sam" > $out$sample".fastq"
        end=$(date +%s.%N)
        time_elapsed=$(echo "scale=1; $end - $start" | bc)
        aligned_reads=$(echo $(cat $out$sample".fastq"|wc -l)/4|bc)
        echo "$file
        Number of aligned reads: $aligned_reads 
        Time required for filtering and conversion: $time_elapsed"
        echo "Removing temporary alignmnet file and compressing filtered fastq file"
        gzip $out$sample".fastq" 
        #if [[ $e == "on" ]] ; then rm $out$sample".sam" ; fi
done 
exit 0

