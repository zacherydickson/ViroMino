#!/bin/bash
set -o errexit
set -o pipefail

ExecDir="$(dirname "$(readlink -f "$0")")"
Positioner="$ExecDir/getReadPositionsByAllelle.pl"
RPBCalculator="$ExecDir/readPosByAlleleTest.R"


if [ "$#" -lt 2 ]; then
	>&2 echo -e "Usage: $(basename $0) bam inVCF.gz [ outVCF.gz = inVCF.gz ]\n" \
		"\tNote: outVCF is erased before processing unless out==in\n" \
        ;
	exit 1;
fi

Bam=$1;shift;
inVCF=$1;shift;
outVCF=$1;
[ -z $outVCF ] && outVCF=$inVCF;

baseVCF=${outVCF/.gz/}

rm -f $baseVCF;
if [ $inVCF != $outVCF ]; then 
	 rm -f $outVCF;
fi

[ -f $Bam.bai ] || samtools index $Bam;

zcat $inVCF |
	while read -r line; do
		if grep -q '^##' <<< $line; then
			echo "$line" >> $baseVCF;
			continue;
		fi
		if grep -q '^#' <<< $line; then
			echo "##INFO=<ID=RPB,Number=A,Type=Integer,Description=\"Phred-scaled read position bias at this position\">" >> $baseVCF;
			echo "$line" >> $baseVCF;
			continue;
		fi
		read -r chrom pos id ref altstr null <<< $line;
		readarray -td, altArr <<< $altstr;
		#>&2 echo "$Positioner $Bam $chrom $pos $ref ${altArr[@]} | $RPBCalculator /dev/stdin"
		RPB=$($Positioner $Bam $chrom $pos $ref ${altArr[@]} | $RPBCalculator /dev/stdin)
		echo "$line" | awk -v RPB=$RPB -F '\\t' 'BEGIN{OFS="\t"}{$8=$8";RPB="RPB}1' >> $baseVCF;
	done

bgzip -f $baseVCF;
bcftools index $outVCF;
