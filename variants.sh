#!/usr/bin/env bash

#TL - Langdon (S) (Durum) (Svevo is durum)
#TD - Dicoccoides (R) (Diccocoides) (Zavitan is diccocoides)

#tophat --no-coverage-search --num-threads 6 --min-intron-length 40 --max-intron-length 50000 --output-dir /media/paolo/MYBACK_B/AAA_SORESI/TOPHAT_OUT_DURUM_NRGENE_gen18//TD1/ /media/paolo/TOSHIBA_EXT/WHEAT_DURUM_13gen17/DURUM_GENOMA_SVEVO/1

#tophat --no-coverage-search --num-threads 6 --min-intron-length 40 --max-intron-length 50000 --output-dir /media/paolo/MYBACK_B/AAA_SORESI/TOPHAT_OUT_DI_FAR_ZAVI_14gen18//TD1/ /media/paolo/TOSHIBA_EXT/WHEAT_EMMER_DICOCC_FARRO_ZAVI/151210_Emmer_zavitan_v2_pseudomolecules.fasta /media/paolo/MYBACK_B/AAA_SORESI/FASTQ/TD1.R1.fastq.gz

cd /Users/juan/Documents/manu/dev/vms/bio/daniela
ln -s /Volumes/toshiba/bio/daniela 
source /Users/juan/Documents/venvs/bio/bin/activate

#python biopyutils/gffSlicer.py -a PGSB_TRITD_Jan2017_all.gff3 -i chr3A -s 54100000 -e 78600000 -o data/svevo.region.gff3
#python biopyutils/gffSlicer.py -a SORTED_WEW_v2_HC_e_LC_GFF3_CATTATI_gff_PAO_26feb18.gff3 -i chr3A -s 54600000 -e 79100000 -o data/wew_v2.region.gff3

samtools faidx 151210_Emmer_zavitan_v2_pseudomolecules.fasta.fa
samtools faidx 160802_Svevo_v2_pseudomolecules.fasta

samtools view -h data/MERGED_TD1_2_3__accepted_hits_SVEVO.bam chr3A:54100000-78700000 > data/dic_svevo.region.sam
samtools view -h data/MERGED_TL1_2_3__accepted_hits_SVEVO.bam chr3A:54100000-78700000 > data/lan_svevo.region.sam
samtools view -h data/MERGED_TD1_2_3__accepted_hits_WEW.bam chr3A:54600000-79100000 > data/dic_wew.region.sam
samtools view -h data/MERGED_TL1_2_3__accepted_hits_WEW.bam chr3A:54600000-79100000 > data/lan_wew.region.sam

samtools view -Sb data/dic_svevo.region.sam > data/dic_svevo.region.bam
samtools view -Sb data/lan_svevo.region.sam > data/lan_svevo.region.bam
samtools view -Sb data/dic_wew.region.sam > data/dic_wew.region.bam
samtools view -Sb data/lan_wew.region.sam > data/lan_wew.region.bam

bcftools mpileup -f data/160802_Svevo_v2_pseudomolecules.fasta data/lan_svevo.region.bam | bcftools call --multiallelic-caller --variants-only -O v -o data/variant_lan_svevo.bcf
bcftools view -i 'AVG(INFO/DP)>=4 & AVG(INFO/MQ)>=50 & (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>0.75' data/variant_lan_svevo.bcf > data/variant_lan_svevo.vcf

bcftools mpileup -f data/160802_Svevo_v2_pseudomolecules.fasta data/dic_svevo.region.bam | bcftools call --multiallelic-caller --variants-only -O v -o data/variant_dic_svevo.bcf
bcftools view -i 'AVG(INFO/DP)>=4 & AVG(INFO/MQ)>=50 & (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>0.75' data/variant_dic_svevo.bcf > data/variant_dic_svevo.vcf

bcftools mpileup -f  data/151210_Emmer_zavitan_v2_pseudomolecules.fasta.fa data/lan_wew.region.bam | bcftools call --multiallelic-caller --variants-only -O v -o data/variant_lan_wew.bcf
bcftools view -i 'AVG(INFO/DP)>=4 & AVG(INFO/MQ)>=50 & (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>0.75' data/variant_lan_wew.bcf > data/variant_lan_wew.vcf

bcftools mpileup -f  data/151210_Emmer_zavitan_v2_pseudomolecules.fasta.fa data/dic_wew.region.bam | bcftools call --multiallelic-caller --variants-only -O v -o data/variant_dic_wew.bcf
bcftools view -i 'AVG(INFO/DP)>=4 & AVG(INFO/MQ)>=50 & (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>0.75' data/variant_dic_wew.bcf > data/variant_dic_wew.vcf



#check vcfs
wc -l data/*.vcf


####BEGIN snpEff databases

#prepare SnpEff databases files
mkdir snpEff/data/svevo
mkdir snpEff/data/zavitan
mkdir snpEff/data/zavitan_v2

#svevo
ln -s 160802_Svevo_v2_pseudomolecules.fasta snpEff/data/svevo/sequences.fa
ln -s PGSB_TRITD_Jan2017_HC_CDS.fasta /snpEff/data/svevo/cds.fa
ln -s PGSB_TRITD_Jan2017_HC_protein.fasta /snpEff/data/svevo/protein.fa
ln -s PGSB_TRITD_Jan2017_all.gff3 snpEff/data/svevo/genes.gff

#zavitan
ln -s 151210_Emmer_zavitan_v2_pseudomolecules.fasta.fa snpEff/data/zavitan/sequences.fa
ln -s TRIDC_WEWseq_PGSB_20160501_CDS.fasta snpEff/data/zavitan/cds.fa
ln -s TRIDC_WEWseq_PGSB_20160501.gtf snpEff/data/zavitan/genes.gtf
ln -s TRIDC_WEWseq_PGSB_20160501_Proteins.fasta snpEff/data/zavitan/protein.fa

#zavitan_v2
ln -s 151210_Emmer_zavitan_v2_pseudomolecules.fasta.fa snpEff/data/zavitan_v2/sequences.fa
ln -s SORTED_WEW_v2_HC_e_LC_GFF3_CATTATI_gff_PAO_26feb18.gff3 snpEff/data/zavitan_v2/genes.gff
ln -s Zavitan_HC_and_LC_predict_protein_version2/transcripts.genes.combined.collapsed_cocla_newid_HC_CDS_protein_OK_2mar17.fasta snpEff/data/zavitan_v2/protein.fa


#create SnpEff databases
cd snpEff
java -Xms5000m -Xmx8000m -jar snpEff.jar build -gff3 -v svevo
java -Xms5000m -Xmx8000m -jar snpEff.jar build -gtf22 -v zavitan
java -Xms5000m -Xmx8000m -jar snpEff.jar build -gff3 -v zavitan_v2

####END snpEff databases

#run SnpEff
cd snpEff
java -Xmx4g -jar snpEff.jar svevo /home/juan/Desktop/juan/bio/svevo-zavitan/data/variant_lan_svevo.vcf > /home/juan/Desktop/juan/bio/svevo-zavitan/data/variant_lan_svevo.ann.vcf
java -Xmx4g -jar snpEff.jar svevo /home/juan/Desktop/juan/bio/svevo-zavitan/data/variant_dic_svevo.vcf > /home/juan/Desktop/juan/bio/svevo-zavitan/data/variant_dic_svevo.ann.vcf
java -Xmx4g -jar snpEff.jar zavitan_v2 /home/juan/Desktop/juan/bio/svevo-zavitan/data/variant_lan_wew.vcf > /home/juan/Desktop/juan/bio/svevo-zavitan/data/variant_lan_wew.ann.vcf
java -Xmx4g -jar snpEff.jar zavitan_v2 /home/juan/Desktop/juan/bio/svevo-zavitan/data/variant_dic_wew.vcf > /home/juan/Desktop/juan/bio/svevo-zavitan/data/variant_dic_wew.ann.vcf

#extract interesting genes
python biopyutils/gff2FA.py -s /Volumes/toshiba/bio/daniela/160802_Svevo_v2_pseudomolecules.fasta -a /Volumes/toshiba/bio/daniela/PGSB_TRITD_Jan2017_all.gff3 -i TRITD3Av1G031440

#Run DEG

cd ../../

htseq-count -m intersection-nonempty --stranded=yes data/lan_svevo.region.sam data/svevo.region.gff3 -t CDS -i Parent > data/lan_svevo.region.counts.sam
htseq-count -m intersection-nonempty --stranded=yes data/dic_svevo.region.sam  data/svevo.region.gff3 -t CDS -i Parent > data/dic_svevo.region.counts.sam
htseq-count -m intersection-nonempty --stranded=yes data/dic_wew.region.sam  data/wew_v2.region.gff3 -t CDS -i Parent > data/lan_wew_v2.region.counts.sam
