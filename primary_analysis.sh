#!/bin/bash
echo "Starting optitype analysis"
abs=$(echo $PWD)

while getopts i:o: option
do case "${option}" in
i) input=${OPTARG};;
o) output=${OPTARG};;
esac
done

echo "Output subdirectory is $output" 
echo "Input directory is $input"

in=$input out=$output bash $abs/optitype_HLA_typing.sh 

exit 0
