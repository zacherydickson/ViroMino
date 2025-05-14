# ViroMino

## Description

A Minor Variant analysis Pipeline optimized for finding minor variants in viral populations.
Given a set of samples with common origin, generates a set of sites which have been confidently called as minor variant sites in each sample.
Additionally extracts the information at those sites in all samples for comparison.

The pipeline runs from raw fastq files, preprocessing with fastp, aligning to the reference with bowtie, calling minor variant sites  with both lofreq and freebays
Filtering is done with custom scripts.

## Dependencies

The current implementation was developed with:
- fastp V0.23.4
- bowtie2 V2.2.5
- lofreq V2.1.5
- freebayes V1.3.6
- samtools V1.21
- bcftools V1.21
- bedtool V2.31.1
- perl V5.32.1
  - Bio::DB::HTS V3.01
  - Bio V1.7.8
- R V4.3.3

## Installation

No compilation required

Clone the repo and its submodules, then place a link in your path
```
    git clone --recurse-submodules https://github.com/zacherydickson/ViroMino.git ViroMino
    ln -s $(readlink -f ViroMino) /somewhere/in/path
```

## Running 

The most basic call is:
```
  ViroMino metaInfo referenceFasta referenceBowtieIndex
```
Use `ViroMino -h` for more options and information

The metaInfo file defines sample names and paths to the raw reads for each sample.

Each line in the file will look like `Sample1:PathtoR1.fq.gz:PathtR2.fq.gz`

## Output

ViroMino will retain the preprocessed fastq files, alignments, as well as raw and filtered calls with each variant caller.

The `CommonCalls.tab` file will lay out each variant site confidently identified in each sample.

The `expanded.vcf.gz` file will include information for each confidently called site in each sample, regardless of whether the site was called in that sample.

The `haplotyped.vcf.gz` file contains the same information, but variants near enough to be supported by the same reads are combined together.
