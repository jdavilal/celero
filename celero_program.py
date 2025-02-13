import argparse
import os
import subprocess
import sys 
import shutil
import time
import math

script_path_list = os.path.normpath(os.path.abspath(__file__)).split(os.sep)
home_dir = os.path.join("/", script_path_list[1], script_path_list[2])
parser = argparse.ArgumentParser()
requiredparser = parser.add_argument_group('required arguments')
parser.add_argument("-ds", "--download_samples", action = 'store_false', help = "True if tumor/normal FASTQ files need to be downloaded")
parser.add_argument("-td", "--tumor_downloaded", type = str, help = "Name of the tumor sample in FASTQ format")
parser.add_argument("-nd", "--normal_downloaded", type = str, help = "Name of the normal sample in FASTQ format")

parser.add_argument("-t", "--tumor_SRR", type = str, help = "Name of the tumor sample from SRA")
parser.add_argument("-n", "--normal_SRR", type = str, help = "Name of the normal sample from SRA")
parser.add_argument("-l", "--synthetic_reference_length", type = int, default = 9, help = "Synthetic reference used")
parser.add_argument("-p", "--kmc_path", type = str, default = str(os.getcwd()), help = "Path to KMC from current directory. Default is your current working directory.")
parser.add_argument("-d", "--delete_files", action = 'store_false', help = "Delete all produced files except the final VCF and BAM files")
parser.add_argument("-o", "--output_name", type = str, default = "output_file", help = "Name of the output file")
parser.add_argument("-r", "--output_results", action = 'store_true', help = "Output all file statistics to a .txt file")
args = parser.parse_args()

download_start = time.time()

if(args.download_samples == True):
  try:
    os.chdir("sratoolkit.3.0.7-ubuntu64")
    os.chdir("bin")
  except FileNotFoundError:
    print("FileNotFoundError")
  except NotADirectoryError:
    print(NotADirectoryError)
  
  SRA_accession_ids = [args.tumor_SRR, args.normal_SRR]
  for id in SRA_accession_ids:
    try:
      print("Currently downloading: " + id)
      prefetch = "./prefetch " + id
      subprocess.call(prefetch, shell=True)

      print("Generating fastq for: " + id)
      fastq_dump = "./fastq-dump " + id
      subprocess.call(fastq_dump, shell=True)
    except subprocess.CalledProcessError as e:
      print(f"Error downloading {accession_number}: {e}")
    except FileNotFoundError:
      print("Error: prefetch command not found. Ensure SRA Toolkit is installed and in your PATH.")

  for id in SRA_accession_ids:
    try:
      parent_dir = os.path.dirname(os.getcwd())
      destination = os.path.join(parent_dir, id + ".fastq")
      shutil.move(id + ".fastq", destination)
      os.chdir("..")
  
      parent_dir = os.path.dirname(os.getcwd())
      destination = os.path.join(parent_dir, id + ".fastq")
      shutil.move(id + ".fastq", destination)
      os.chdir("..")
      os.chdir("sratoolkit.3.0.7-ubuntu64")
      os.chdir("bin")
    except Exception as e:
      print(f"{e}")
  
  os.chdir("..")
  os.chdir("..")
else:
  SRA_accession_ids = [args.tumor_downloaded, args.normal_downloaded]
  
  
for id in SRA_accession_ids:
  try:
    seqtk = "seqtk/seqtk trimfq -b 6 -e 45 " + id + ".fastq > " + id + "_trimfq.fastq"
    subprocess.call(seqtk, shell=True)
    if(args.delete_files == True and args.download_samples == True):
      os.remove(id + ".fastq")
  except Exception as e:
    print(f"{e}")

download_end = time.time()
download_time = download_end - download_start

for id in SRA_accession_ids:
  kmc = str(args.kmc_path) + "/kmc -k27 " + id + "_trimfq.fastq " + id + "_kmer_num ."
  subprocess.call(kmc, shell=True)

if(args.download_samples == True):
  tumor_sample = args.tumor_SRR
  normal_sample = args.normal_SRR
else:
  tumor_sample = args.tumor_downloaded
  normal_sample = args.normal_downloaded
  
try: 
  kmc_simple = str(args.kmc_path) + "/kmc_tools simple " + tumor_sample + "_kmer_num " + normal_sample + "_kmer_num -ci10 kmers_subtract final_complement"
  subprocess.call(kmc_simple, shell=True)
except Exception as e:
  print(f"{e}")
  
try:
  kmc_simple2 = str(args.kmc_path) + "/kmc_tools simple final_complement complete_ref_kmer_" + str(args.synthetic_reference_length) + " -ci5 intersect intersect_combine_" + str(args.synthetic_reference_length)
  subprocess.call(kmc_simple2, shell=True)
  kmc_transform = str(args.kmc_path) + "/kmc_tools transform intersect_combine_" + str(args.synthetic_reference_length) + " dump intersect_combine_kmer" + str(args.synthetic_reference_length) + ".txt"
  subprocess.call(kmc_transform, shell=True)
except Exception as e:
  print(f"{e}")

input_file = "intersect_combine_kmer" + str(args.synthetic_reference_length) + ".txt"
output_file = "kmer" + str(args.synthetic_reference_length) + ".txt"

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
  for line in infile:
    fields = line.split()
    if len(fields) > 1 and fields[1].isdigit() and int(fields[1]) > 5:
      outfile.write(line)

try:
  to_fasta = "python3 to_fasta.py final_complement.fasta kmer" + str(args.synthetic_reference_length) + ".txt"
  subprocess.call(to_fasta, shell=True)
except Exception as e:
  print(f"{e}")

try:
  star = "STAR --runThreadN 8 --genomeDir " + str(os.getcwd()) + "/data/GenomeDir \
  --readFilesIn " + str(os.getcwd()) + "/final_complement.fasta \
  --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 39540834309"
  subprocess.call(star, shell=True)
except Exception as e:
  print(f"{e}") 

try:
  bcftools1 = "bcftools-1.18/bcftools mpileup -f data/ref_genome.fa -o final_complement_indels.bcf Aligned.sortedByCoord.out.bam"
  subprocess.call(bcftools1, shell=True)
  bcftools2 = "bcftools-1.18/bcftools index final_complement_indels.bcf"
  subprocess.call(bcftools2, shell=True)
  bcftools3 = "bcftools-1.18/bcftools call --skip-variants snps --multiallelic-caller \
  --variants-only -O v final_complement_indels.bcf -o final_complement_indels.vcf"
  subprocess.call(bcftools3, shell=True)
except Exception as e:
  print(f"{e}") 
  
process_end = time.time()
process_time = process_end - download_end

os.rename("final_complement_indels.vcf", str(args.output_name) + ".vcf")  
os.rename("Aligned.sortedByCoord.out.bam", str(args.output_name) + ".bam")

print("download time: " + str(download_time))
print("process time: " + str(process_time))

with open("kmer" + str(args.synthetic_reference_length) + ".txt") as file:
  kmer_count = sum(1 for _ in file)
    
with open(str(args.output_name) + ".vcf", "r") as file:
    indel_count = sum(1 for line in file if not line.startswith("#"))
  
deletion_count = 0
with open(str(args.output_name) + ".vcf", "r") as file:
    for line in file:
        if not line.startswith("#"):
            columns = line.split()
            if len(columns) > 4 and len(columns[3]) > len(columns[4]):
                deletion_count += 1

indel_density = (indel_count / kmer_count)*1000
deletion_density = (deletion_count / kmer_count)*1000

log_odds = -3.38 + 0.83*(deletion_density) - 0.45*(indel_density)
odds = math.exp(log_odds)
probability = odds / (1+odds)
if probability > .5:
  status_prediction = "MSI-H"
else:
  status_prediction = "MSS"
  
if(args.output_results == True):
  with open(str(args.output_name)+"_results.txt", "w") as file:
    file.write(str(args.output_name) + "\n***********RESULTS**********\n" + "Number of k-mers: " + str(kmer_count) + "\nNumber of indels: " + str(indel_count) + "\nNumber of deletions: " + str(deletion_count)+"\n")
    file.write("Indel density: " + str(indel_density) + "\nDeletion density: " + str(deletion_density) + "\n*****************************\n")
    file.write("Download time: " + str(download_time) + "\nProcess time: " + str(process_time) + "\n*****************************\n")
    file.write("Log odds: " + str(log_odds) + "\nOdds: " + str(odds) + "\nProbability: " + str(probability) + "\n")
    file.write("MSI STATUS PREDICTION: " + str(status_prediction))

if(args.delete_files == True):
  os.remove("final_complement_indels.bcf")
  os.remove("final_complement_indels.bcf.csi")
  os.remove("Log.final.out")
  os.remove("Log.out")
  os.remove("Log.progress.out")
  os.remove("SJ.out.tab")
  os.remove("final_complement.fasta")
  os.remove("intersect_combine_" + str(args.synthetic_reference_length) + ".kmc_pre")
  os.remove("intersect_combine_" + str(args.synthetic_reference_length) + ".kmc_suf")
  os.remove("intersect_combine_kmer" + str(args.synthetic_reference_length) + ".txt")
  os.remove("kmer" + str(args.synthetic_reference_length) + ".txt")
  os.remove("final_complement.kmc_pre")
  os.remove("final_complement.kmc_suf")
  os.remove(tumor_sample + "_kmer_num.kmc_pre")
  os.remove(tumor_sample + "_kmer_num.kmc_suf")
  os.remove(normal_sample + "_kmer_num.kmc_pre")
  os.remove(normal_sample + "_kmer_num.kmc_suf")
  os.remove(tumor_sample + "_trimfq.fastq")
  os.remove(normal_sample + "_trimfq.fastq")
  print("Excess files deleted")
