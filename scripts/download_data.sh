#!/bin/bash
# Script to process NCBI fastq sequences
# to then use for a BLAST search within Genbank

# Sanam Sajjadi, ssajjadi@dons.usfca.edu
# November 22, 2018

export PATH=$PATH:/Users/sanamsajjadi/Desktop/Sanamsajjadi-final-project/sratoolkit.2.9.2-mac64/bin
# Make sure to have sra-tool kit in library
# Download the list of files in the SRA run table
# to the raw data directory
# the pipe and tail -n+2 is a good way to exclude the first line
for SRA_number in $(cut -f 8 data/metadata/SraRunTable.txt | tail -n +2)
do
  fastq-dump -v "$SRA_number" -O data/raw_data
 
done
echo "List of files downloaded."