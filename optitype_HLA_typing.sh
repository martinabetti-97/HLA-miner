#!/bin/bash
for dir in $in/*R1.fastq.gz;
do
        file=$(basename $dir)
        sample_id=${file%_*}
        echo "running analysis on:
        file 1: $in/$sample_id"_R1.fastq.gz"
        file 2: $in/$sample_id"_R2.fastq.gz" "
        start=$(date +%s.%N)
        docker run -v "$in":/data/ -t fred2/optitype -i $sample_id"_R1.fastq.gz"
 $sample_id"_R2.fastq.gz" -r -o $out -p $sample_id -v
        end=$(date +%s.%N)
        time_elapsed=$(echo "scale=1; $end - $start" | bc)
	FILE=$in/$out$sample_id"_result.tsv"
	if test -f "$FILE"; then
		echo "Optitype analysis on $sample_id completed."
    		echo "$sample_id $time_elapsed"  >> $in/$out"time.txt"
	else 
    		echo "Warning: Optiytpe analysis on $sample_id was not succesful
. It could be a memory issue, try to operate the pre-filtering step on this samp
le."
		echo "$sample_id failed"  >> $in/$out"time.txt"
	fi
done

echo "Summarizing results..."
abs=$(echo $PWD)
cd $in/$out
python3 $abs"/optitype_results.py"
echo "Analysis completed."

exit 0
