#!/bin/bash

beds=$(echo "$1" | tr ";" " ")
chrom_sizes=$2
segmentation_file=$3
TSS_pc_lincRNA_antiS=$4
odir=$5

#Merge all regions peaks in the same file
cat $beds | bedtools sort -i stdin | bedtools merge -i stdin > "${odir}pm.bed"
#Filter out regions that would be out of the genome range
bedtools intersect -a "${odir}pm.bed" -b $chrom_sizes > "${odir}pm.classic.bed"
#Sort
bedtools sort -i "${odir}pm.classic.bed" > "${odir}pm.classic.sorted.bed"

#Intersect with regions of the input matrix(ces)
bedtools intersect -a $segmentation_file -b "${odir}pm.classic.sorted.bed" > "${odir}pm.cutbywindow.bed"
bedtools sort -i "${odir}pm.cutbywindow.bed" > "${odir}pm.cutbywindow.sorted.bed" 
bedtools closest -a "${odir}pm.cutbywindow.sorted.bed" -b $TSS_pc_lincRNA_antiS -wa -wb -d > "${odir}pm.annot.gene.bed"
bedtools closest -a "${odir}pm.cutbywindow.sorted.bed" -b $segmentation_file -wa -d > "${odir}pm.annot.window.bed"

rm "${odir}pm.bed"
rm "${odir}pm.classic.bed"
rm "${odir}pm.classic.sorted.bed"
rm "${odir}pm.cutbywindow.bed"
rm "${odir}pm.cutbywindow.sorted.bed"
rm  $segmentation_file
# awk '{ if ($8 <= 1000) { print } }' pm.annot.gene.bed > closer1000.bed
