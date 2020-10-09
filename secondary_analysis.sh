abs=$(echo $PWD)
wd1=$abs"/output" 
wd2=$abs"/files"  
e=$(echo 'off') #default variable for ethnicity analysis

while getopts i:f:o:eh option
do case "${option}" in
i) Info="${OPTARG}" ;;
f) Hla="${OPTARG}" ;;
e) e=$(echo "on") ;;
o) Out=${OPTARG}  
	wd1=$Out ;;
h)	echo "Usage:"
	echo "		-h	Display this help message."
	echo "		-i 	Path to input info file."		 
	echo "		-f	Path to hla results table."
	echo "		-e	This parameter activates the ethnicity analysis."
	echo "		-o	Path where results will be saved."      
	exit 0 
	;;
esac
done

#------------------------------------------------------------------------

if [ ! -d "$Out" ] ; then mkdir $Out ; fi

cd $wd1

echo 'Collecting files...'

if [ ! -f "$Hla" ] ; then echo "File $Hla not found." ; exit 0 ;
else cp $Hla . ; fname=$(basename $Hla) ; if [[ $fname != 'hla.tsv' ]] ; then mv $fname hla.tsv ; fi ;  
fi

if [ ! -f "$Info" ] ; then echo "File $Info not found." ; exit 0 ;
else cp $Info . ; fname=$(basename $Info) ; if [[ $fname != 'info.tsv' ]] ; then mv $fname info.tsv ; fi ; 
fi

cp -n $wd2'/supertypes.xls' . 
cp -n $wd2'/sequence.fa' . 
cp -n $wd2'/distance.tsv' . 

#------------------------------------------------------------------------
python3 $abs/data_preparation.py

echo "
Dividing samples for condition.. and calculating group-specific allele frequency"

if [[ $e == "on" ]] ; then Rscript $abs/ethnicity.R ; fi 

python3 $abs/calc_freq.py $e

#-----------------------------------------------------------------------
echo "
Calculating Grantham distance"
python3 $abs/hed.py
#-----------------------------------------------------------------------
echo "
Plotting graphs "

echo "."
#python3 $abs/matrix.py
R CMD BATCH $abs/circlize.R 

echo "."
R CMD BATCH $abs/heatmap.R 

echo "."
if [[ $e == "on" ]] ; then R CMD BATCH $abs/ethnicity_heatmap.R ; fi 

echo "Done"
#------------------------------------------------------------------------
echo "
Running statistical analysis"
R CMD BATCH $abs/stat.R
#------------------------------------------------------------------------
echo "
Removing intermediate files..."

rm freq_hla* hla_* hlaST_* 
#rm sequence.fa distance.tsv info.tsv supertypes.xls
rm *.Rout 

echo "Analysis completed!"
