#!/bin/bash

outdir="Results/"
SAMPLES=$(find Data/*.txt -type f -printf '%f\n')

for s in $SAMPLES; do
	samplename="$(echo $s | cut -d. -f1-1)"
	cp Data/$s Data/Input_Data.txt
	echo "Beginning LDA 1-feature selection for $samplename."
	./OneFeature
	mv Results/LDA_output_1_feature.txt $outdir/$samplename-1_feature.txt
	echo "Finished with LDA 1-feature selection for $samplename., beginning 2-feature selection"
	./TwoFeature
	mv Results/LDA_output_2_feature.txt $outdir/$samplename-2_feature.txt
	echo "Finished with LDA 2-feature selection for $samplename, beginning 3-feature selection"
	./ThreeFeature
	mv Results/LDA_output_3_feature.txt $outdir/$samplename-3_feature.txt
	echo "Finished with LDA 1, 2 & 3 feature selection for $samplename"
	rm Data/Input_Data.txt
done

##