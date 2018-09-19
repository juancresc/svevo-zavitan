#!/usr/bin/env bash

#TL - Langdon (S) (Durum) (Svevo is durum)
#TD - Dicoccoides (R) (Diccocoides) (Zavitan is diccocoides)

#tophat --no-coverage-search --num-threads 6 --min-intron-length 40 --max-intron-length 50000 --output-dir /media/paolo/MYBACK_B/AAA_SORESI/TOPHAT_OUT_DURUM_NRGENE_gen18//TD1/ /media/paolo/TOSHIBA_EXT/WHEAT_DURUM_13gen17/DURUM_GENOMA_SVEVO/1

#tophat --no-coverage-search --num-threads 6 --min-intron-length 40 --max-intron-length 50000 --output-dir /media/paolo/MYBACK_B/AAA_SORESI/TOPHAT_OUT_DI_FAR_ZAVI_14gen18//TD1/ /media/paolo/TOSHIBA_EXT/WHEAT_EMMER_DICOCC_FARRO_ZAVI/151210_Emmer_zavitan_v2_pseudomolecules.fasta /media/paolo/MYBACK_B/AAA_SORESI/FASTQ/TD1.R1.fastq.gz

cd /Users/juan/Documents/manu/dev/vms/bio/daniela
ln -s /Volumes/toshiba/bio/daniela disk

source /Users/juan/Documents/venvs/bio/bin/activate

#python biopyutils/gffSlicer.py -a disk/PGSB_TRITD_Jan2017_all.gff3 -i chr3A -s 54100000 -e 78600000 -o disk/data/svevo.region.gff3
#python biopyutils/gffSlicer.py -a disk/SORTED_WEW_v2_HC_e_LC_GFF3_CATTATI_gff_PAO_26feb18.gff3 -i chr3A -s 54600000 -e 79100000 -o disk/data/wew_v2.region.gff3

samtools faidx 151210_Emmer_zavitan_v2_pseudomolecules.fasta.fa
samtools faidx 160802_Svevo_v2_pseudomolecules.fasta

samtools view -h disk/MERGED_TD1_2_3__accepted_hits_SVEVO.bam chr3A:54100000-78600000 > disk/data/dic_svevo.region.sam
samtools view -h disk/MERGED_TL1_2_3__accepted_hits_SVEVO.bam chr3A:54100000-78600000 > disk/data/lan_svevo.region.sam
samtools view -h disk/MERGED_TD1_2_3__accepted_hits_WEW.bam chr3A:54600000-79100000 > disk/data/dic_wew.region.sam
samtools view -h disk/MERGED_TL1_2_3__accepted_hits_WEW.bam chr3A:54600000-79100000 > disk/data/lan_wew.region.sam

samtools view -Sb disk/data/dic_svevo.region.sam > disk/data/dic_svevo.region.bam
samtools view -Sb disk/data/lan_svevo.region.sam > disk/data/lan_svevo.region.bam
samtools view -Sb disk/data/dic_wew.region.sam > disk/data/dic_wew.region.bam
samtools view -Sb disk/data/lan_wew.region.sam > disk/data/lan_wew.region.bam

bcftools mpileup -f disk/160802_Svevo_v2_pseudomolecules.fasta disk/data/lan_svevo.region.bam | bcftools call --multiallelic-caller --variants-only -O v -o disk/data/variant_lan_svevo.bcf
bcftools view disk/data/variant_lan_svevo.bcf > disk/data/variant_lan_svevo.vcf

bcftools mpileup -f disk/160802_Svevo_v2_pseudomolecules.fasta disk/data/dic_svevo.region.bam | bcftools call --multiallelic-caller --variants-only -O v -o disk/data/variant_dic_svevo.bcf
bcftools view disk/data/variant_dic_svevo.bcf > disk/data/variant_dic_svevo.vcf

bcftools mpileup -f disk/151210_Emmer_zavitan_v2_pseudomolecules.fasta.fa disk/data/lan_wew.region.bam | bcftools call --multiallelic-caller --variants-only -O v -o disk/data/variant_lan_wew.bcf
bcftools view  disk/data/variant_lan_wew.bcf > disk/data/variant_lan_wew.vcf

bcftools mpileup -f disk/151210_Emmer_zavitan_v2_pseudomolecules.fasta.fa disk/data/dic_wew.region.bam | bcftools call --multiallelic-caller --variants-only -O v -o disk/data/variant_dic_wew.bcf
bcftools view disk/data/variant_dic_wew.bcf > disk/data/variant_dic_wew.vcf

####BEGIN snpEff databases

#prepare SnpEff databases files
mkdir disk/snpEff/data/svevo
mkdir disk/snpEff/data/zavitan
mkdir disk/snpEff/data/zavitan_v2

#svevo
ln -s disk/160802_Svevo_v2_pseudomolecules.fasta disk/snpEff/data/svevo/sequences.fa
ln -s disk/PGSB_TRITD_Jan2017_HC_CDS.fasta /snpEff/data/svevo/cds.fa
ln -s disk/PGSB_TRITD_Jan2017_HC_protein.fasta /snpEff/data/svevo/protein.fa
ln -s   /snpEff/data/svevo/genes.gff

#zavitan
ln -s disk/151210_Emmer_zavitan_v2_pseudomolecules.fasta.fa disk/snpEff/data/zavitan/sequences.fa
ln -s disk/TRIDC_WEWseq_PGSB_20160501_CDS.fasta disk/snpEff/data/zavitan/cds.fa
ln -s disk/TRIDC_WEWseq_PGSB_20160501.gtf disk/snpEff/data/zavitan/genes.gtf
ln -s disk/TRIDC_WEWseq_PGSB_20160501_Proteins.fasta disk/snpEff/data/zavitan/protein.fa

#zavitan_v2
ln -s disk/151210_Emmer_zavitan_v2_pseudomolecules.fasta.fa disk/snpEff/data/zavitan_v2/sequences.fa
ln -s disk/SORTED_WEW_v2_HC_e_LC_GFF3_CATTATI_gff_PAO_26feb18.gff3 disk/snpEff/data/zavitan_v2/genes.gff
ln -s disk/Zavitan_HC_and_LC_predict_protein_version2/transcripts.genes.combined.collapsed_cocla_newid_HC_CDS_protein_OK_2mar17.fasta disk/snpEff/data/zavitan_v2/protein.fa


#create SnpEff databases
java -Xms5000m -Xmx8000m -jar snpEff.jar build -gff3 -v svevo
java -Xms5000m -Xmx8000m -jar snpEff.jar build -gtf22 -v zavitan
java -Xms5000m -Xmx8000m -jar snpEff.jar build -gff3 -v zavitan_v2

####END snpEff databases

#run SnpEff
cd disk/snpEff
java -Xmx4g -jar snpEff.jar svevo ../data/variant_lan_svevo.vcf > ../data/variant_lan_svevo.ann.vcf
java -Xmx4g -jar snpEff.jar svevo ../data/variant_dic_svevo.vcf > ../data/variant_dic_svevo.ann.vcf
java -Xmx4g -jar snpEff.jar zavitan_v2 ../data/variant_lan_wew.vcf > ../data/variant_lan_wew.ann.vcf
java -Xmx4g -jar snpEff.jar zavitan_v2 ../data/variant_dic_wew.vcf > ../data/variant_dic_wew.ann.vcf

#Run DEG

cd ../../

htseq-count -m intersection-nonempty --stranded=yes disk/data/lan_svevo.region.sam disk/data/svevo.region.gff3 -t CDS -i Parent > disk/data/lan_svevo.region.counts.sam
htseq-count -m intersection-nonempty --stranded=yes disk/data/dic_svevo.region.sam  disk/data/svevo.region.gff3 -t CDS -i Parent > disk/data/dic_svevo.region.counts.sam
htseq-count -m intersection-nonempty --stranded=yes disk/data/dic_wew.region.sam  disk/data/wew_v2.region.gff3 -t CDS -i Parent > disk/data/lan_wew_v2.region.counts.sam
