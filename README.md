# Introduction
Celero is a fast detection method of microsatellite status in RNA-seq data that leverages k-mer methods. 

## Overview of Workflow
<img width="458" alt="Screenshot 2025-02-18 at 11 17 18 AM" src="https://github.com/user-attachments/assets/44640da3-be84-4ba7-9e89-fcec0e9962c7" />

## Additional Data
- Synthetic references (includes l = 7, 9, 11, 13)

# Running Celero
- Predict MSI status with already downloaded tumor and normal FASTQ sample files
```
python3 celero_program.py -td tumor.fastq -nd normal.fastq -l 9 -p /path/to/kmc -o output
```
- Predict MSI status and download tumor and normal FASTQ sample files
```
python3 celero_program.py -ds -t tumorSRR -n normalSRR -l 9 -p /path/to/kmc -o output
```

# Input and output
## Input
### FASTQ Files
- One tumor and one normal FASTQ file are required as input to run `Celero`
  - These files can be downloaded by specifying the `-ds` parameter and inputting the SRR number for both tumor and normal FASTQ samples
### Synthetic Reference
- One synthetic reference is required to run `Celero`
- The reference used is specified using the `-l` parameter, with the default set to `l=9`
 
## Output
- Predicted MSI status is output to the terminal unless the `-r` parameter is used, generating an output_file_results.txt file with all results
- All generated files except the BAM and VCF files are deleted unless the `-d` parameter is specified 

# Installation guide
## Dependencies
All dependencies must be downloaded and installed to run `Celero`. All dependencies except `KMC` are assumed to be in the current working directory when running `Celero`.

- SRA-Toolkit (https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)
- seqtk (https://github.com/lh3/seqtk)
- KMC (https://github.com/refresh-bio/KMC)
- STAR (https://github.com/alexdobin/STAR)
- Bcftools (https://github.com/samtools/bcftools)
## Installing the package
To install `Celero`, download and install all dependencies and run the following:
```
git clone https://github.com/jdavilal/celero
```

# Resources

* The analysis and R code for our [paper](analysis_final.md). 
