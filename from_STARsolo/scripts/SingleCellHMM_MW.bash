# bash SingleCellHMM.bash  Path_to_bam_file/PREFIX.bam numberOfThread Path_to_SingleCellHMM.R

INPUT_BAM=$1 #path to .bam file
CORE=$2 # number of cores for parallelization
MEM=$3
MERGEBP=$4
THRESH=$5
OUTDIR=$6 #DWM; STARsolo output directory path
PL=$7 #path to SingleCellHMM.R

CORE="${CORE:-5}"
#MINCOV="${MINCOV:-5}"
MERGEBP="${MERGEBP:-500}"
THRESH="${THRESH:-10000000}"
CURDIR=`pwd` #snakemake directory
PL="${PL:-${CURDIR}/scripts}"

#can this samtools call be parallelized with -@ ?
reads=`samtools view -q 255 ${INPUT_BAM} | wc -l`
echo "Number of aligned reads is ${reads}"
minCovReads=`expr ${reads} / ${THRESH}`
MINCOV=${minCovReads}

# PREFIX=`echo ${INPUT_BAM} | rev | cut -d / -f 1 | cut -d . -f 2- | rev`
PREFIX="HMMout"

TMPDIR=${OUTDIR}/${PREFIX}_HMM_features
mkdir -p ${TMPDIR}

exec > >(tee SingleCellHMM_Run_${TMPDIR}.log)
exec 2>&1
echo "Path to SingleCellHMM.R:  	${PL}"
echo "input .bam                	${INPUT_BAM}"
echo "STARsolo output directory: 	${OUTDIR}"
echo "tmp folder:                	${TMPDIR}"
echo "number of threads:         	${CORE}"
echo "memory usage:              	${MEM}"
echo "minimum coverage:						${MINCOV}"
echo "thresholded at 1 in ${THRESH} reads"
echo ""
echo "Reads spanning over splicing junction will join HMM blocks"
echo "To avoid that, split reads into small blocks before input to groHMM"
echo "Spliting and sorting reads..."
bedtools bamtobed -i ${INPUT_BAM} -split | LC_ALL=C sort -k1,1V -k2,2n --buffer-size=${MEM} --parallel=${CORE} | awk '{print $0}' | gzip > ${TMPDIR}/${PREFIX}_split.sorted.bed.gz

cd ${TMPDIR}
zcat ${PREFIX}_split.sorted.bed.gz  | awk '{print $0 >> $1".bed"}'
find -name "*.bed" -size -1024k -delete
wc chr*.bed -l > chr_read_count.txt

echo ""
echo "Start to run groHMM on each individual chromosome..."


wait_a_second() {
	joblist=($(jobs -p))
    while (( ${#joblist[*]} >= ${CORE} ))
	    do
	    sleep 1
	    joblist=($(jobs -p))
	done
}

for f in *.bed
do
  wait_a_second
  echo "    " ${f}
  R --vanilla --slave --args $(pwd) ${f}  < ${PL}/SingleCellHMM.R  > ${f}.log 2>&1 & pids+=($!)
done
wait "${pids[@]}"


echo ""
echo "Merging HMM blocks within ${MERGEBP}bp..."
for f in *_HMM.bed
do
  LC_ALL=C sort -k1,1V -k2,2n --parallel=${CORE} ${f} > ${f}.sorted.bed
  cat ${f}.sorted.bed | grep + > ${f}_plus
  cat ${f}.sorted.bed | grep - > ${f}_minus
  bedtools merge -s -d ${MERGEBP} -i ${f}_plus > ${f}_plus_merge${MERGEBP} & pids2+=($!)
  bedtools merge -s -d ${MERGEBP} -i ${f}_minus > ${f}_minus_merge${MERGEBP} & pids2+=($!)
  wait_a_second
done
wait "${pids2[@]}"

echo "Combining HMM output from all chromosomes..."
cat *_HMM.bed_plus_merge${MERGEBP}  | awk 'BEGIN{OFS="\t"} {print $0, ".", ".", "+"}' > ${PREFIX}_merge${MERGEBP}
cat *_HMM.bed_minus_merge${MERGEBP} | awk 'BEGIN{OFS="\t"} {print $0, ".", ".", "-"}' >> ${PREFIX}_merge${MERGEBP}

mkdir toremove
for f in *_HMM.bed
do
	mv ${f}.sorted.bed ${f}_plus ${f}_minus ${f}_plus_merge${MERGEBP} ${f}_minus_merge${MERGEBP} toremove/.
done


echo ""
echo "Sorting combined .bed file..."
f=${PREFIX}
LC_ALL=C sort -k1,1V -k2,2n ${f}_merge${MERGEBP} --parallel=${CORE} > ${f}_merge${MERGEBP}.sorted.bed
#-k1,1V
# rm ${f}_merge${MERGEBP} #TODO
# make a toy genome file for bedtools to specify sort order (chromosome lengths are all 42, not needed here)
zcat ${f}_split.sorted.bed.gz | awk {'print $1'} | sort -k1,1V | uniq | sed s/$/'\t42'/ > tmp.genome

echo ""
echo "Calculating the coverage..."
bedtools coverage -nonamecheck -a ${f}_merge${MERGEBP}.sorted.bed -b <(zcat ${f}_split.sorted.bed.gz) -s -counts -split -sorted -g tmp.genome > ${f}_merge${MERGEBP}.sorted.bed_count
rm tmp.genome

echo ""
echo "Filtering the HMM blocks by coverage..."
cat ${f}_merge${MERGEBP}.sorted.bed_count | awk 'BEGIN{OFS="\t"} ($7 >= '$MINCOV'){print $1, $2, $3, $4, $5, $6, $7}' | gzip > TAR_reads.bed.gz

echo ""
echo "#### Please examine if major chromosomes are all present in the final TAR_reads.bed.gz file ####"
zcat TAR_reads.bed.gz | cut -f 1 | uniq

echo ""
echo "Link the final TAR_reads.bed.gz file to the working directory"
cd ..
ln -s ${TMPDIR}/TAR_reads.bed.gz .


echo ""
echo "Move intermediate files to  ${TMPDIR}/toremove ..."
echo ""
echo "${TMPDIR}/toremove can be deleted if no error message in SingleCellHMM_Run log file and all major chromosomes are present in the final TAR_reads.bed.gz file"
echo ""

cd ${TMPDIR}
mv chr* toremove/.

for f in *
do gzip --quiet ${f} &
done

# cd toremove
# for f in *
# do rm ${f} &
# done

cd ${CURDIR}

echo "Done!"
