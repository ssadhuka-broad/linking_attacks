cat datasets/GSE68951_family.soft | grep "SAMPLE = GSM" | cut -d "=" -f2 > datasets/GSMs.txt

while read p; do
	TIMEPT=$(cat datasets/GSE68951_family.soft | grep "SAMPLE = ${p}" -A 14 | grep timepoint | cut -d ":" -f2)
	PAT_ID=$(cat datasets/GSE68951_family.soft | grep "SAMPLE = ${p}" -A 14 | grep "patient id" | cut -d ":" -f2)
	echo "${p} ${PAT_ID} ${TIMEPT}" >> GSM_conversion.txt
done <datasets/GSMs.txt
#cat datasets/GSE68951_family.soft | grep "SAMPLE = GSM1688368" -A 14 | grep timepoint | cut -d ":" -f2