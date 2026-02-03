# RNA SHAPE Analysis Pipeline
Scripts for running RNA SHAPE (Selective 2'-Hydroxyl Acylation analyzed by Primer Extension) analysis on tRNA data.

## Overview
This pipeline processes high-throughput sequencing data from tRNA SHAPE experiments through quality control, alignment, UMI deduplication, reactivity quantification, and secondary structure prediction.

## tRNA Pipeline

This contains scripts for running SHAPE analysis on tRNA data for Nuno Santos and Arianna Di Fazio.

### Pipeline Execution Order

#### Step 1: Prepare Samplesheet

**01_samplesheet.sh** - Generates a samplesheet catalog for downstream analysis

Creates a CSV file listing all samples from the input directory.

Usage:
```bash
bash tRNA/01_samplesheet.sh <input_directory> <output_samplesheet.csv>
```

Example:
```bash
bash tRNA/01_samplesheet.sh /path/to/fastq/ samplesheet.csv
```

---

#### Step 2: Preprocess Reads

**02_preprocess_reads.sh** - Prepares raw sequencing reads for alignment

Performs the following preprocessing steps:
1. Adapter trimming (removes sequencing adapters)
2. UMI (Unique Molecular Identifier) extraction
3. Hard-clipping to a specified read length
4. Collapse paired-end reads across R1 and R2
5. Combine with R1 singleton reads
6. Reverse complement strands as needed
7. Rearrange FASTQ headers to preserve UMI at end (separated by underscore)

Output: Preprocessed FASTQ files with standardized headers

---

#### Step 3: Merge Sequencing Lanes (Optional)

**merge_lanes.sh** - Combines reads from multiple sequencing lanes

If samples were sequenced across multiple lanes, this script merges R1 and R2 reads of the same sample from different lanes into unified files.

Output: Lane-merged FASTQ files

---

#### Step 4: Concatenate Sample Batches (Optional)

**concat_reads.sh** - Merges reads from multiple sequencing runs

Concatenates reads from different sequencing batches or experiments into a single set for unified downstream analysis.

Output: Concatenated FASTQ files

---

#### Step 5: Quality Control

**03_fastqc.sh** - Performs quality assessment of preprocessed reads

Runs FastQC on all preprocessed FASTQ files to assess read quality metrics:
- Per-base quality scores
- GC content distribution
- Adapter presence
- Sequence duplication levels

Generates MultiQC report for summary visualization.

Output: FastQC HTML reports and MultiQC summary

---

#### Step 6: Build Alignment Index

**04a_bowtie2_index.sh** - Creates Bowtie2 alignment index for reference sequences

Builds searchable indices from reference tRNA sequences using Bowtie2.

Input: Reference FASTA files (e.g., tRNA_Ala_AGC_2_1.fa, tRNA_Pro_TGG_3_5.fa)

Output: Bowtie2 index files (.bt2)

---

#### Step 7: Align Reads to References

**04b_bowtie2.sh** - Aligns preprocessed reads to reference tRNA sequences

Maps preprocessed reads to tRNA reference genomes using Bowtie2.
- Aligns against multiple indices (Ala and Pro tRNAs)
- Converts alignments to sorted BAM format
- Generates alignment statistics

Output: Sorted BAM files (.bam) with alignment statistics

---

#### Step 8: UMI Deduplication

**umidedup.sh** - Removes PCR duplicates using UMI information

Deduplicates aligned reads by collapsing reads with identical UMIs to single consensus reads, reducing PCR amplification bias.

Input: Sorted BAM files
Output: UMI-deduplicated BAM files

---

#### Step 9: Gather UMI Statistics

**umigather.sh** - Collects UMI deduplication metrics

Extracts and summarizes UMI deduplication statistics from log files, creating a summary CSV of input/output read counts.

Output: CSV file with UMI deduplication statistics

---

#### Step 10: Quantify Reactivity

**05_rf-count.sh** - Counts nucleotide-level SHAPE reactivity using RFcount

Calculates per-nucleotide SHAPE reactivities from aligned BAM files using the rf-count tool from the Leonore singularity container. Processes both DMS and 5NIA treated samples.

Input: UMI-deduplicated BAM files and reference FASTA files
Output: Reactivity count files (.rc format)

---

#### Step 11: Normalize Reactivity

**07_rf-norm.sh** - Normalizes SHAPE reactivity values

Normalizes reactivity values across samples and conditions using rf-normalize:
- Groups samples by tRNA type (Ala, Pro, Pro_Dic) and reagent (DMS, 5NIA, DMSO)
- Applies inter-sample normalization
- Scales values to 0-1 range based on population average

Input: Reactivity count files
Output: Normalized reactivity XML files

---

#### Step 12: Generate Wiggle Tracks

**08_rf-wiggle.sh** - Converts reactivity data to genome browser format

Creates Wiggle format files for visualization in genome browsers:
- One track per sample per nucleotide position
- Allows visualization of per-nucleotide reactivity across the tRNA

Input: Normalized reactivity XML files
Output: Wiggle (.wig) files for genome browser visualization

---

#### Step 13: Predict Secondary Structure with Reactivity Constraints

**09_rf-fold.sh** - Predicts tRNA secondary structure using SHAPE data

Folds tRNA sequences using RNAfold with SHAPE reactivity constraints:
- Uses normalized reactivities to guide structure prediction
- Generates structure prediction in dot-bracket notation
- Produces multiple output formats (CT, DotBracket, ViennaRNA)

Input: Normalized reactivity XML files and tRNA sequences
Output: Folded structure files in multiple formats

---

#### Step 14: Generate Coverage Analysis and Combined Folds

**10_rftools_cov.sh** - Calculates per-nucleotide coverage statistics

Computes coverage metrics from reactivity count files to assess sequencing depth per position.

Output: Coverage files for analysis quality assessment

---

**11_rf_comb_fold.sh** - Combines reactivities and generates consensus structures

Merges reactivity data from multiple replicates per condition and generates consensus secondary structures:
- Combines replicates for each treatment condition
- Predicts consensus structure with averaged reactivities
- Outputs consensus structures for downstream analysis

Input: Reactivity count files from multiple replicates
Output: Combined reactivity files and consensus structure predictions

---

### Supporting Analysis Scripts

**CovPlots.R** - Generates coverage plots for visualization

**deltaShape.R** - Calculates and plots differential SHAPE reactivity between conditions

**shapeplotr.R** - R package functions for SHAPE data visualization

**shapeplotr.sh** - Bash wrapper for R-based SHAPE plotting functions