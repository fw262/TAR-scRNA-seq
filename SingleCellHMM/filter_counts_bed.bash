input=$1
read=$2
output=$3

cat ${input} | awk -v r=${read} 'BEGIN{OFS="\t"} ($7 >= r){print $1, $2, $3, $4, $5, $6}' | gzip > ${output}.bed.gz


