# Convert .bam file to .bed format, and split across chromosomes. Files saved in ./beds

INPUT_BAM=$1 #path to .bam file
CORE=$2 # number of cores for parallelization
OUTDIR=$3 #DWM; cellranger .../count_directory/multicr_TAR/beds

PREFIX=`echo ${INPUT_BAM} | rev | cut -d / -f 1 |cut -d . -f 2- |rev` #this is the same for all cellranger pipeline, could just directly name it here
CURDIR=`pwd` #snakemake directory

mkdir ${OUTDIR}/beds

echo "Spliting and sorting reads for ${PREFIX}..."
bedtools bamtobed -i ${INPUT_BAM} -split | LC_ALL=C sort -k1,1V -k2,2n --parallel=30 | awk '{print $0}' | gzip > ${OUTDIR}/${PREFIX}_split.sorted.bed.gz

cd ${OUTDIR}
# TODO - should this step be gzipped? coul
zcat ${PREFIX}_split.sorted.bed.gz  | awk '{print $0 >> "chr"$1".bed"}' | gzip
# find -name "chr*.bed" -size -1024k -delete

cd ${CURDIR}
