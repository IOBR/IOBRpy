#!/usr/bin/env python3
import os
import random
import subprocess
import argparse

try:
    from tqdm.auto import tqdm
except Exception:  # tqdm not installed; degrade gracefully
    tqdm = None


def set_ulimit():
    """
    Set the maximum number of open files to 65535 for the current process.
    This helps avoid issues with too many open files during large-scale processing.
    """
    try:
        # Attempt to set the ulimit value
        subprocess.run("ulimit -n 65535", shell=True, check=True)
        # print("ulimit set to 65535 successfully.")
    except subprocess.CalledProcessError:
        print("Warning: ulimit raise failed.")

def process_sample(f1, path_out, index, suffix1, num_threads):
    """
    Process a single sample using STAR for alignment.
    Args:
        f1 (str): Path to the first FASTQ file.
        path_out (str): Output directory for STAR results.
        index (str): Path to STAR index file.
        suffix1 (str): Suffix for the first FASTQ file.
        num_threads (int): Number of threads to use for STAR processing.
    """
    # Set ulimit before processing
    set_ulimit()
    
    # Build corresponding second FASTQ filename by replacing '1' with '2' in the suffix
    suffix2 = suffix1.replace("1", "2")
    f2 = f1.replace(suffix1, suffix2)

    # Extract sample ID assuming no other underscores in the filename
    sample_id = os.path.basename(f1).replace(suffix1, "")

    # Check if task.complete and BAM files already exist and are not empty
    if os.path.exists(os.path.join(path_out, f"{sample_id}.task.complete")) and \
            os.path.exists(os.path.join(path_out, f"{sample_id}_Aligned.sortedByCoord.out.bam")) and \
            os.path.getsize(os.path.join(path_out, f"{sample_id}_Aligned.sortedByCoord.out.bam")) > 0:
        # Output progress message
        print(f"Skipping {sample_id} - task.complete and BAM files exist and are not empty.")
    else:
        # Output progress message
        print(f">>> Processing {sample_id}...")

        # Construct STAR command (splitted into multiple lines)
        command = (
            f"STAR --genomeDir {index} "
            f"--readFilesIn {f1} {f2} "
            f"--outFileNamePrefix {os.path.join(path_out, sample_id)}_ "
            f"--twopassMode Basic "
            f"--runThreadN {num_threads} "
            f"--outBAMsortingThreadN {num_threads} "
            f"--outTmpDir {os.path.join(path_out, sample_id + '_STARtmp')} "
            f"--readFilesCommand zcat "
            f"--outSAMtype BAM SortedByCoordinate "
            f"--quantMode GeneCounts "
            f"--limitBAMsortRAM 137438953472"
        )
        
        # Run STAR command using subprocess
        subprocess.run(command, shell=True, check=True)

        # Create task.complete file
        os.makedirs(os.path.join(path_out, sample_id), exist_ok=True)
        open(os.path.join(path_out, f"{sample_id}.task.complete"), 'w').close()
        print(f">>> {sample_id} has completed the STAR process. ")
        print(f" ")

def _print_outputs_summary(files_1, path_out, suffix1):
    """
    Print the final output file locations for each sample.
    This does not change computation, only reports paths.
    """
    if not files_1:
        return
    print("\n===== Output files summary =====")
    for f1 in sorted(files_1):
        sample_id = os.path.basename(f1).replace(suffix1, "")
        bam = os.path.join(path_out, f"{sample_id}_Aligned.sortedByCoord.out.bam")
        gene = os.path.join(path_out, f"{sample_id}_ReadsPerGene.out.tab")
        # Print whether file exists
        bam_status = "OK" if os.path.exists(bam) and os.path.getsize(bam) > 0 else "missing"
        gene_status = "OK" if os.path.exists(gene) and os.path.getsize(gene) > 0 else "missing"
        print(f"{sample_id}\n  BAM:        {bam} [{bam_status}]\n  GeneCounts: {gene} [{gene_status}]")

def process_samples(path_fq, path_out, index, suffix1="_1.fastq.gz", batch_size=1, num_threads=8):
    """
    Process FASTQ files using STAR for alignment.
    Args:
        path_fq (str): Path to directory containing FASTQ files.
        path_out (str): Output directory for STAR results.
        index (str): Path to STAR index file.
        suffix1 (str): Suffix for the first FASTQ file.
        batch_size (int): Number of samples to process concurrently.
        num_threads (int): Number of threads to use for STAR processing.
    """
    # Create output directory if it does not exist
    os.makedirs(path_out, exist_ok=True)

    # Get all _1.fastq.gz files and shuffle them
    files_1 = sorted([os.path.join(path_fq, file) for file in os.listdir(path_fq) if file.endswith(suffix1)])
    random.shuffle(files_1)

    # Progress setup
    total = len(files_1)
    if total == 0:
        print(f"No files with suffix '{suffix1}' found under: {path_fq}")
        return
    pbar = tqdm(total=total, desc="STAR (GeneCounts)", unit="sample") if tqdm else None
    processed = 0

    # Initialize batch index
    batch_index = 0
    while batch_index < len(files_1):
        # Get batch of files
        batch_files = files_1[batch_index:batch_index + batch_size]

        # Process batch of files
        for f1 in batch_files:
            process_sample(f1, path_out, index, suffix1, num_threads)
            if pbar:
                pbar.update(1)
            else:
                processed += 1
                if (processed % 10 == 0) or (processed == total):
                    print(f"[Progress] {processed}/{total} samples")

        # Increment batch index
        batch_index += batch_size

    if pbar:
        pbar.close()

    # Final summary of output files
    _print_outputs_summary(files_1, path_out, suffix1)

# Main execution block
def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description="Process FASTQ files using STAR for alignment.")
    parser.add_argument("--index", type=str, required=True, help="Path to STAR index file.")
    parser.add_argument("--path_fq", type=str, required=True, help="Path to directory containing FASTQ files.")
    parser.add_argument("--path_out", type=str, required=True, help="Output directory for STAR results.")
    parser.add_argument("--suffix1", type=str, default="_1.fastq.gz", help="Suffix for the first FASTQ file.")
    parser.add_argument("--batch_size", type=int, default=1, help="Number of samples to process concurrently.")
    parser.add_argument("--num_threads", type=int, default=8, help="Number of threads to use for STAR processing.")
    args = parser.parse_args()

    # Call the function to process FASTQ files
    process_samples(args.path_fq, args.path_out, args.index, args.suffix1, args.batch_size, args.num_threads)

if __name__ == '__main__':
    main()