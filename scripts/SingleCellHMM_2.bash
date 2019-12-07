# bash SingleCellHMM.bash  Path_to_bam_file/PREFIX.bam numberOfThread Path_to_SingleCellHMM.R

INPUT_BAM=$1 #pbmc4k_possorted_genome_bam.bam
CORE=$2
MINCOV=$3
MERGEBP=$4
PL=$5

${CORE:=5}
${MINCOV:=5}
${MERGEBP:=500}
${PL:=./scripts}

PREFIX=`echo ${INPUT_BAM} | rev | cut -d / -f 1 |cut -d . -f 2- |rev`
tmp="HMM_features"
TMPDIR=${PREFIX}_${tmp}
mkdir ${TMPDIR}

exec > >(tee SingleCellHMM_Run_${TMPDIR}.log)
exec 2>&1
echo "Path to SingleCellHMM.R   $PL" 
echo "INPUT_BAM                 $INPUT_BAM"
echo "temp folder               $TMPDIR"
echo "number Of thread          $CORE"
echo "minimum coverage		$MINCOV"
echo ""
echo "Reads spanning over splicing junction will join HMM blocks"
echo "To avoid that, split reads into small blocks before input to groHMM"
echo "Spliting and sorting reads..."
bedtools bamtobed -i ${INPUT_BAM} -split |LC_ALL=C sort -k1,1V -k2,2n --parallel=30| awk '{print $0}' | gzip > ${TMPDIR}/${PREFIX}_split.sorted.bed.gz 

cd ${TMPDIR}
zcat ${PREFIX}_split.sorted.bed.gz  |awk '{print $0 >> "chr"$1".bed"}' 
find -name "chr*.bed" -size -1024k -delete
#wc chr*.bed -l > chr_read_count.txt

echo ""
echo "Start to run groHMM in each individual chromosome..."


wait_a_second() {
	joblist=($(jobs -p))
    while (( ${#joblist[*]} >= ${CORE} ))
	    do
	    sleep 1
	    joblist=($(jobs -p))
	done
}


for f in chr*.bed
do 
wait_a_second
R --vanilla --slave --args $(pwd) ${f}  < ${PL}/SingleCellHMM.R  > ${f}.log 2>&1 &
done
#R --vanilla --slave --args $(pwd) ${PREFIX}_split.sorted.bed.gz  < ${PL}/SingleCellHMM.R 
wait


echo ""
echo "Merging HMM blocks within ${MERGEBP}bp..."
for f in chr*_HMM.bed
do	
  LC_ALL=C sort -k1,1V -k2,2n --parallel=30 ${f} > ${f}.sorted.bed
  cat ${f}.sorted.bed | grep + > ${f}_plus
  cat ${f}.sorted.bed | grep - > ${f}_minus
  bedtools merge -s -d ${MERGEBP} -i ${f}_plus > ${f}_plus_merge${MERGEBP} &
  bedtools merge -s -d ${MERGEBP} -i ${f}_minus > ${f}_minus_merge${MERGEBP} &
  wait_a_second
done

wait 


#f=${PREFIX}_split.sorted_HMM
#gzip ${f}.bed &

cat chr*_HMM.bed_plus_merge${MERGEBP} | awk 'BEGIN{OFS="\t"} {print $0, ".", ".", "+"}' > ${PREFIX}_merge${MERGEBP}
cat chr*_HMM.bed_minus_merge${MERGEBP} | awk 'BEGIN{OFS="\t"} {print $0, ".", ".", "-"}' >> ${PREFIX}_merge${MERGEBP}

mkdir toremove
for f in chr*_HMM.bed
do	
mv ${f}.sorted.bed ${f}_plus ${f}_minus ${f}_plus_merge${MERGEBP} ${f}_minus_merge${MERGEBP} toremove/.
done


echo ""
echo "Calculating the coverage..." 
f=${PREFIX}
LC_ALL=C sort -k1,1V -k2,2n ${f}_merge${MERGEBP} --parallel=30 > ${f}_merge${MERGEBP}.sorted.bed
rm ${f}_merge${MERGEBP}


bedtools coverage -a ${f}_merge${MERGEBP}.sorted.bed -b <(zcat ${PREFIX}_split.sorted.bed.gz) -s -counts -split -sorted > ${f}_merge${MERGEBP}.sorted.bed_count

echo ""
echo "Filtering the HMM blocks by coverage..." 
#cat ${f}_merge${MERGEBP}.sorted.bed_count | awk 'BEGIN{OFS="\t"} ($7 >= 2){print $1, $2, $3, $4, $5, $6}' | gzip > ${f}_merge${MERGEBP}_2reads.bed.gz
cat ${f}_merge${MERGEBP}.sorted.bed_count | awk 'BEGIN{OFS="\t"} ($7 >= '$MINCOV'){print $1, $2, $3, $4, $5, $6, $7}' | gzip > ${f}_merge${MERGEBP}_${MINCOV}reads.bed.gz

echo "" 
echo "#### Please examine if major chromosomes are all present in the final PREFIX_merge${MERGEBP}_${MINCOV}reads.bed.gz file ####"
zcat ${f}_merge${MERGEBP}_${MINCOV}reads.bed.gz |cut -f 1 |uniq

echo "" 
echo "Link the final PREFIX_merge${MERGEBP}_${MINCOV}reads.bed.gz file to the working directory"
cd ..
ln -s ${TMPDIR}/${f}_merge${MERGEBP}_${MINCOV}reads.bed.gz .


echo ""
echo "Move intermediate files to  ${TMPDIR}/toremove ..." 
echo ""
echo "${TMPDIR}/toremove can be deleted if no error message in SingleCellHMM_Run log file and "
echo "all major chromosomes are present in the final PREFIX_merge${MERGEBP}_${MINCOV}reads.bed.gz file"
echo ""
echo ""

cd ${TMPDIR}
mv chr* toremove/.

for f in *
do gzip ${f} &
done

cd toremove
for f in *
do gzip ${f} &
done

cd ../..

echo ""
echo "done!"

