cd /home/juan/Desktop/juan/bio/svevo-zavitan
wew=/home/juan/Desktop/juan/bio/svevo-zavitan/disk/151210_Emmer_zavitan_v2_pseudomolecules.fasta.fa
svevo=/home/juan/Desktop/juan/bio/svevo-zavitan/disk/160802_Svevo_v2_pseudomolecules.fasta
path_hisat2=/home/juan/Desktop/juan/bio/sw/hisat2-2.1.0

#MERGED_TD1_2_3__accepted_hits_SVEVO.bam
#MERGED_TL1_2_3__accepted_hits_SVEVO.bam

#Program: bcftools (Tools for variant calling and manipulating VCFs and BCFs)
#Version: 1.7 (using htslib 1.7-2)

#TD
bcftools mpileup -f data/160802_Svevo_v2_pseudomolecules.fasta data/MERGED_TD1_2_3__accepted_hits_SVEVO.bam | bcftools call --multiallelic-caller --variants-only -O v -o data/variant_td_svevo_all.bcf
bcftools view data/variant_td_svevo_all.bcf | bcftools filter -e "QUAL < 50" | bcftools filter -i "AVG(INFO/DP)>=4 & (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>0.75 "  > data/variant_td_svevo_all.filtered.vcf

#TL
bcftools mpileup -f data/160802_Svevo_v2_pseudomolecules.fasta data/MERGED_TL1_2_3__accepted_hits_SVEVO.bam | bcftools call --multiallelic-caller --variants-only -O v -o data/variant_tl_svevo_all.bcf
bcftools view data/variant_tl_svevo_all.bcf | bcftools filter -e "QUAL < 50" | bcftools filter -i "AVG(INFO/DP)>=4 & (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>0.75 "  > data/variant_tl_svevo_all.vcf

#run notebook
bcftools view data/diff.vcf | bcftools filter -e "QUAL < 50" | bcftools filter -i "AVG(INFO/DP)>=4 & (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>0.75 "  > data/diff.filtered.vcf

#SNPEFF
java -Xmx4g -jar data/snpEff/snpEff.jar svevo data/variant_td_svevo_all.vcf > data/variant_td_svevo_all.ann.vcf
java -Xmx4g -jar data/snpEff/snpEff.jar svevo data/variant_tl_svevo_all.vcf > data/variant_tl_svevo_all.ann.vcf

#counts

samtools view -h -o data/td.sam data/MERGED_TD1_2_3__accepted_hits_SVEVO.bam
samtools view -h -o data/tl.sam data/MERGED_TL1_2_3__accepted_hits_SVEVO.bam

htseq-count -m intersection-nonempty --stranded=yes data/td.sam data/PGSB_TRITD_Jan2017_all.gff3 -t CDS -i Parent > data/td.counts.csv
htseq-count -m intersection-nonempty --stranded=yes data/tl.sam data/PGSB_TRITD_Jan2017_all.gff3 -t CDS -i Parent > data/tl.counts.csv




#demos
    

bcftools mpileup -f --per-sample-mF --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \ data/160802_Svevo_v2_pseudomolecules.fasta -b data/bams.txt | bcftools call --multiallelic-caller --variants-only -O v -o data/all.2.bcf
bcftools view data/all.2.bcf | bcftools filter -e "QUAL < 50" | bcftools filter -i "AVG(INFO/DP)>=4 & (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>0.75 "  > data/all.filtered.2.vcf