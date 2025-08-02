# ITMfinder: Intratumoral Microbiome Finder Tool

ITMfinder is a Python software tool designed for identifying and analyzing intratumoral microbiomes. It processes FASTQ data through multiple steps, including quality control, host gene removal, taxonomic analysis, and integration of results, enabling researchers to gain insights into the tumor microenvironment.

## Features
- Quality control of input FASTQ data using `fastp`
- Host gene removal using `Bowtie2`
- Taxonomic analysis using `Kraken2`, with abundance estimation using `Bracken`
- Integration of results and diversity analysis

## Installation

### Prerequisites
- Python 3.6+
- [fastp](https://github.com/OpenGene/fastp)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [Kraken2](https://ccb.jhu.edu/software/kraken2/)
- [Bracken](https://ccb.jhu.edu/software/bracken/)

To clone the repository and install required Python packages:
```sh
# Clone the repository
git clone https://github.com/your-username/itmfinder.git
cd itmfinder

# Install required Python packages
pip install -r requirements.txt
```

## Usage
ITMfinder is a command-line tool that executes different analysis steps by specifying the desired step.

### Basic Usage
```sh
itmfinder --step [1-8] [additional arguments]
```

### Arguments
- `--step`
  - Specifies the step to execute (1-8).
- `--num_threads`
  - Number of threads to use for processing (default: `8`).
- `--batch_size`
  - Number of samples to process in parallel (default: `1`).
- `--memory_mapping`
  - Use memory mapping for Kraken2 to reduce memory usage during analysis (Step 6). This reduces memory consumption by using memory-mapped files instead of loading everything into memory, which is helpful if system memory is limited.

### Analysis Steps
ITMfinder comprises eight main steps, each representing a different phase of the analysis:

1. **Step 1: Quality Control (FASTQ Data)**
   Performs quality control on input FASTQ data using `fastp`.
   ```sh
   itmfinder --step 1 \
     --path1_fastq /path/to/raw_fastq \
     --path2_fastp /path/to/cleaned_fastq \
     --num_threads 8 \
     --suffix1 _1.fastq.gz \
     --length_required 50 \
     --batch_size 1
   ```

2. **Step 2: Host Gene Removal**
   Removes host gene sequences using `Bowtie2`.
   ```sh
   itmfinder --step 2 \
     --path2_fastp /path/to/cleaned_fastq \
     --path3_hr /path/to/host_removed \
     --db_bowtie2 /path/to/bowtie2_db \
     --num_threads 8 \
     --batch_size 1
   ```

3. **Step 3: (Optional) Second Round of Host Gene Removal**
   Performs an additional round of host gene removal if required.
   ```sh
   itmfinder --step 3 \
     --path3_hr /path/to/host_removed \
     --path3_hr2 /path/to/host_removed_round2 \
     --db_bowtie2 /path/to/bowtie2_db \
     --num_threads 8 \
     --batch_size 1
   ```

4. **Step 4: Taxonomic Analysis (Kraken2)**
   Performs taxonomic classification using `Kraken2`.
   ```sh
   itmfinder --step 4 \
     --path3_hr /path/to/host_removed \
     --path4_ku1 /path/to/kraken2_output \
     --db_ku /path/to/kraken2_db \
     --num_threads 8 \
     --memory_mapping
   ```

5. **Step 5: Extract Microbiome Sequences and Remove Contaminants**
   Extracts microbiome sequences and removes potential contaminants.
   ```sh
   itmfinder --step 5 \
     --path3_hr /path/to/host_removed \
     --path4_ku1 /path/to/kraken2_output \
     --path5_mr /path/to/microbiome_reads \
     --TID 2 10239 4751 2157 \
     --batch_size 1
   ```

6. **Step 6: Kraken2 and Bracken Analysis**
   Performs Kraken2 analysis again and uses `Bracken` for abundance estimation.
   ```sh
   itmfinder --step 6 \
     --path5_mr /path/to/microbiome_reads \
     --path6_ku2 /path/to/kraken2_output2 \
     --path7_bracken /path/to/bracken_output \
     --db_ku /path/to/kraken2_db \
     --num_threads 8 \
     --memory_mapping \
     --force
   ```

7. **Step 7: Integrate Analysis Data**
   Merges all sample data into a matrix for integrated analysis.
   ```sh
   itmfinder --step 7 \
     --path6_ku2 /path/to/kraken2_output2 \
     --path7_bracken /path/to/bracken_output \
     --path8_mpa /path/to/merged_mpa \
     --mpa_suffix .kraken.mpa.std.txt \
     --report_suffix kraken.report.txt \
     --output_file 1-combine_mpa_std.txt
   ```

8. **Step 8: Data Summary**
   Performs file size statistics to identify potential outliers in the data.
   ```sh
   itmfinder --step 8 \
     --path2_fastp /path/to/cleaned_fastq \
     --path3_hr /path/to/host_removed \
     --path6_ku2 /path/to/kraken2_output2 \
     --path5_mr /path/to/microbiome_reads \
     --out_dir /path/to/data_summary
   ```

## Troubleshooting

### 1. Memory Issues
`Kraken2` and `Bracken` are memory-intensive programs. If you experience memory-related issues:
- Reduce the `--num_threads` value to lower memory usage.
- Use the `--memory_mapping` flag to reduce Kraken2 memory usage.

### 2. Long File Paths
If file paths are too long, errors may occur. Consider shortening directory names or paths.

## Contributions
We welcome contributions and suggestions! Please submit any issues or feature requests on [GitHub Issues](https://github.com/your-username/itmfinder/issues).

## License
This project is licensed under the MIT License. See [LICENSE](LICENSE) for more details.

## Download
To download ITMfinder, visit the [GitHub repository](https://github.com/your-username/itmfinder) and follow the installation instructions.

