#!/bin/bash

# This script performs the analysis of publicly available groundwater metagenomes from SRA
#Required Tools: 
# METAXA2 version 2.2.3 
# DIAMOND 
# cutadapt 
# The script Needs the ResFinder Database Version 4, translated from nucleotides into amino-acids. 
# 

mkdir DLPBL
cd    DLPBL
mkdir logs
mkdir ARG-hits
mkdir PhyHits

export  PATH=$PATH:"/home/ioannis/ShotgunSequencingAnalysis/sratoolkit/bin" 

echo "Staring downloads"

export  PATH=$PATH:"/home/ioannis/bin"  

echo "Path exported"




for i in   	ERR2191546 ERR3281021 ERR3281022 ERR3281023 ERR3281024 ERR3281025 ERR3281026 ERR3281027 ERR3281028 ERR3281029 ERR3281030 ERR3281031	\
ERR3281034ERR3281035 ERR3281036 ERR3281037  \
ERR3281038 ERR3281039 ERR3281040 ERR3281041 ERR4388077 ERR4388078 ERR4388079 ERR4388080 SRR10419845 SRR10992209 SRR10992210 SRR3546449 SRR3546450 SRR3546451 SRR3546452 SRR3546453	\
SRR3546454 SRR3546456 SRR3546457 SRR6208705 SRR6209370 SRR6211631 SRR6212576 SRR6213089 SRR6213120 SRR6213885 SRR6222115 SRR8797411 SRR8797412	\
SRR8797413 SRR8797414 SRR9644739 SRR9644740 SRR9644741 SRR9644742 SRR2638363 SRR2638365 SRR2638367 SRR2638368 SRR2638369 SRR2638370 SRR2638366 SRR2638364 	\
SRR2090165 SRR2090167 SRR2090168 SRR2090172 SRR2090164 SRR2090166 SRR2090169 SRR2090171 SRR2090173 SRR2090175 SRR4388699 SRR4389122 SRR4389352 \
SRR2090174 SRR8512245 SRR8512246 SRR8512247	\
SRR8512248 SRR8512249 SRR8512250 SRR8512251	\
SRR12666050 SRR12918151 SRR12918152 SRR12926096	\
SRR12927484 SRR11844736 SRR11844737 SRR11844738 \
SRR11844739 SRR11844740 SRR11844784 SRR11844785 	\
SRR13662870 SRR13662871 SRR11600161 SRR11600162   \  
SRR11600163 SRR11614986 SRR11614987 SRR11614988
	                         
	 
	do prefetch "$i"
	echo "Download Complete"
	mv "$i"/"$i".sra  "$i"/..
	echo "Moving Complete"
        fastq-dump --split-files "$i".sra
	echo "FASTQdump Complete"
	rm "$i" -r
	rm "$i".sra
	echo "SRA file removal Complete"
	~/.local/bin/cutadapt --cores=10        --minimum-length 90 --max-n 0 --max-ee 1  	  \
     	-o trimmed_"$i"_1.fastq \
      	"$i"_1.fastq  > logs/"$i"_1.fastq.log
	~/.local/bin/cutadapt --cores=10        --minimum-length 90 --max-n 0 --max-ee 1  	  \
     	-o trimmed_"$i"_2.fastq \
      	"$i"_2.fastq  > logs/"$i"_2.fastq.log
	echo "Allinging ARGs to sequencing data with Diamond"
	diamond blastx -d   resfinderdb -q trimmed_"$i"_1.fastq -o ARG-hits/"$i"_1.txt   --min-orf 30 --id 99    -e 1e-10 --more-sensitive -f 6  --max-target-seqs 1 --threads 10
	diamond blastx -d   resfinderdb -q trimmed_"$i"_2.fastq -o ARG-hits/"$i"_2.txt  --min-orf 30 --id 99    -e 1e-10 --more-sensitive -f 6  --max-target-seqs 1 --threads 10
	metaxa2 -1 trimmed_"$i"_1.fastq -2 trimmed_"$i"_2.fastq -f q -o PhyHits/"$i"_sample  --cpu 10  
         cd    PhyHits
	 rm *.fasta
	 rm *.graph
	 rm *.extraction.results
	 cd ../
	 rm *.fastq
	echo "Cycle complete"

done
      
echo "All complete"








