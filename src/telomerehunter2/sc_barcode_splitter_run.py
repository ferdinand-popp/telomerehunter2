#!/usr/bin/env python3
"""
sc-ATAC Barcode Splitter and Telomerehunter2 Runner
Splits sc-ATAC BAM files by cell barcodes (Sinto) and processes them with telomerehunter2
"""

import argparse
import datetime
import logging
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pysam

# Setup logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def sanitize_barcode(barcode):
    """Convert barcode to safe filename"""
    # Replace problematic characters with underscores
    safe_barcode = re.sub(r"[^a-zA-Z0-9_-]", "_", barcode)
    return safe_barcode


def get_barcodes_from_bam(
    input_bam,
    max_barcodes=None,
    barcode_txt_path=None,
    skip_barcode_file=False,
    min_reads_per_barcode=10000,
):
    logger.info(f"Barcodes from {input_bam}")

    barcode_counts = defaultdict(int)
    processed_reads = 0

    with pysam.AlignmentFile(input_bam, "rb") as bamfile:
        for read in bamfile.fetch(until_eof=True):
            processed_reads += 1
            if processed_reads % 100000 == 0:
                logger.info(
                    f"Processed {processed_reads} reads, found {len(barcode_counts)} barcodes"
                )
            bc = read.get_tag("CB") if read.has_tag("CB") else None
            if bc:
                barcode_counts[bc] += 1

    logger.info(
        f"Found: {len(barcode_counts)} unique barcodes in {processed_reads} reads"
    )

    # Filter barcodes by min_reads_per_barcode
    filtered_barcodes = {
        bc: count
        for bc, count in barcode_counts.items()
        if count >= min_reads_per_barcode
    }
    logger.info(
        f"{len(filtered_barcodes)} barcodes with at least {min_reads_per_barcode} reads"
    )

    # Optionally limit number of barcodes
    sorted_barcodes = sorted(
        filtered_barcodes.items(), key=lambda x: x[1], reverse=True
    )
    if max_barcodes and len(sorted_barcodes) > max_barcodes:
        logger.info(f"Limiting to top {max_barcodes} barcodes")
        sorted_barcodes = sorted_barcodes[:max_barcodes]

    # Write barcodes as: bc<TAB>bc<TAB>count
    if barcode_txt_path and not skip_barcode_file:
        with open(barcode_txt_path, "w") as f:
            for bc, count in sorted_barcodes:
                f.write(f"{bc}\t{bc}\t{count}\n")
        logger.info(f"Barcode list written to {barcode_txt_path}")
        logger.info(
            f"{len(sorted_barcodes)} barcodes written to file after filtering and limiting."
        )

    # Return dict: barcode -> read count
    return dict(sorted_barcodes)


def read_barcodes_from_file(barcode_txt_path, max_barcodes=None):
    """Read barcodes from an existing barcode file (strip fake group if present)"""
    logger.info(f"Reading barcodes from existing file: {barcode_txt_path}")
    barcodes = []
    with open(barcode_txt_path, "r") as f:
        for line in f:
            bc = line.strip()
            if bc:
                bc = bc.split("\t")[0]  # Remove fake group if present
                barcodes.append(bc)
    if max_barcodes and len(barcodes) > max_barcodes:
        logger.info(f"Limiting to top {max_barcodes} barcodes from file")
        barcodes = barcodes[:max_barcodes]
    logger.info(f"Loaded {len(barcodes)} barcodes from file")
    # No read counts available from file, so return as set
    return set(barcodes)


def split_with_sinto(
    input_bam, output_dir, patient_id, barcode_dict, barcode_txt_path, max_parallel=4
):
    """Use Sinto for efficient BAM splitting"""
    logger.info("Using Sinto for BAM splitting")

    # Prepare BAM output subfolder
    bam_out_dir = output_dir / "bam_files_by_barcode"
    bam_out_dir.mkdir(parents=True, exist_ok=True)

    # Check if BAM index exists, if not, create it
    bam_index = str(input_bam) + ".bai"
    if not os.path.exists(bam_index):
        logger.info(f"BAM index {bam_index} not found. Creating index with pysam...")
        pysam.index(str(input_bam))
        logger.info("BAM index created.")

    cmd = [
        "sinto",
        "filterbarcodes",
        "-b",
        str(input_bam),
        "-c",
        str(barcode_txt_path),
        "--outdir",
        str(bam_out_dir),
        "-p",
        str(max_parallel),
    ]

    logger.info(f"Running Sinto: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        logger.info("Sinto completed successfully")

        # Collect output files
        barcode_files = {}
        for barcode in barcode_dict:
            safe_bc = sanitize_barcode(barcode)
            bam_file = bam_out_dir / f"{safe_bc}.bam"

            if bam_file.exists():
                # Index BAM
                pysam.index(str(bam_file))
                # Rename to include patient ID
                new_name = bam_out_dir / f"{patient_id}_{safe_bc}.bam"
                shutil.move(str(bam_file), str(new_name))
                shutil.move(str(bam_file) + ".bai", str(new_name) + ".bai")
                barcode_files[barcode] = {
                    "file": str(new_name),
                    "safe_barcode": safe_bc,
                    "read_count": barcode_dict[barcode],
                }

        return barcode_files

    except subprocess.CalledProcessError as e:
        logger.error(f"Sinto failed: {e}")
        logger.error(f"Stdout: {e.stdout}")
        logger.error(f"Stderr: {e.stderr}")
        return None


def run_telomerehunter2(
    patient_barcode_id, bam_file, cytobanding_file, output_dir, th2_threads=1
):
    """Run telomerehunter2 on a single BAM file, mit Thread-Anzahl."""
    cmd = [
        "telomerehunter2",
        "-p",
        patient_barcode_id,
        "-ibt",
        str(bam_file),
        "-b",
        str(cytobanding_file),
        "-o",
        str(output_dir),
        "-pno",  # no plotting
        "-c",
        str(th2_threads),
    ]

    logger.info(
        f"Running telomerehunter2 for {patient_barcode_id} with {th2_threads} threads"
    )
    logger.debug(f"Command: {' '.join(cmd)}")

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(th2_threads)

    try:
        result = subprocess.run(
            cmd, cwd=output_dir, capture_output=True, text=True, check=True, env=env
        )
        logger.info(f"telomerehunter2 completed successfully for {patient_barcode_id}")
        return True, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        logger.error(f"telomerehunter2 timed out for {patient_barcode_id}")
        return False, "", "Process timed out"
    except subprocess.CalledProcessError as e:
        logger.error(
            f"telomerehunter2 failed for {patient_barcode_id} (exit code: {e.returncode})"
        )
        return False, e.stdout, e.stderr
    except Exception as e:
        logger.error(
            f"Unexpected error running telomerehunter2 for {patient_barcode_id}: {e}"
        )
        return False, "", str(e)


def th2_wrapper(args):
    """Wrapper for running telomerehunter2, for multiprocessing."""
    (
        barcode,
        file_info,
        patient_id,
        cytobanding_file,
        th2_results_dir,
        cleanup_bams,
        th2_threads,
    ) = args
    patient_barcode_id = f"{patient_id}_{file_info['safe_barcode']}"
    bam_file = file_info["file"]

    try:
        success, stdout, stderr = run_telomerehunter2(
            patient_barcode_id,
            bam_file,
            cytobanding_file,
            th2_results_dir,
            th2_threads=th2_threads,
        )
    except Exception as e:
        return barcode, {
            "success": False,
            "stdout": "",
            "stderr": str(e),
            "patient_barcode_id": patient_barcode_id,
        }

    # Clean up
    if cleanup_bams:
        try:
            os.unlink(bam_file)
            os.unlink(bam_file + ".bai")
        except Exception as e:
            logger.warning(f"Failed to clean up {bam_file}: {e}")

    return barcode, {
        "success": success,
        "stdout": stdout,
        "stderr": stderr,
        "patient_barcode_id": patient_barcode_id,
    }


def process_barcodes_with_telomerehunter2(
    barcode_files,
    patient_id,
    cytobanding_file,
    th2_results_dir,
    cleanup_bams=True,
    max_parallel=4,
):
    """Process all barcode BAM files with telomerehunter2, limiting parallelism and cores per job"""

    total_cores = multiprocessing.cpu_count()
    max_parallel = max_parallel if max_parallel < 6 else 6
    th2_threads = max(1, total_cores // max_parallel)
    logger.info(
        f"Processing {len(barcode_files)} barcodes with max {max_parallel} concurrent processes, each using {th2_threads} threads"
    )

    results = {}

    # Prepare argument tuples for each barcode
    args_list = [
        (
            barcode,
            file_info,
            patient_id,
            cytobanding_file,
            th2_results_dir,
            cleanup_bams,
            th2_threads,
        )
        for barcode, file_info in barcode_files.items()
    ]

    with ProcessPoolExecutor(max_workers=max_parallel) as executor:
        future_to_barcode = {
            executor.submit(th2_wrapper, args): args[0] for args in args_list
        }

        for i, future in enumerate(as_completed(future_to_barcode), 1):
            barcode = future_to_barcode[future]
            try:
                bc, result = future.result()
                results[bc] = result
                status = "SUCCESS" if result["success"] else "FAILED"
                logger.info(f"[{i}/{len(future_to_barcode)}] {bc} -> {status}")
            except Exception as e:
                logger.error(f"Exception processing {barcode}: {e}")
                results[barcode] = {"success": False, "error": str(e)}

    return results


def write_summary_report(
    results, patient_id, input_bam, th2_results_dir, cytobanding_file
):
    """Write detailed summary report"""

    summary_file = th2_results_dir / "telomerehunter2_summary.txt"

    successful = sum(1 for r in results.values() if r["success"])

    with open(summary_file, "w") as f:
        f.write("=== sc-ATAC Telomerehunter2 Processing Summary ===\n")
        f.write(f"Patient ID: {patient_id}\n")
        f.write(f"Input BAM: {input_bam}\n")
        f.write(f"Cytobanding file: {cytobanding_file}\n")
        f.write(f"Output directory: {th2_results_dir}\n")
        f.write(f"Total barcodes: {len(results)}\n")
        f.write(f"Successfully processed: {successful}\n")
        f.write(f"Failed: {len(results) - successful}\n")
        f.write(f"Success rate: {successful / len(results) * 100:.1f}%\n")
        f.write("\n=== Per-Barcode Results ===\n")
        f.write("Barcode\tStatus\tPatient_Barcode_ID\tError\n")

        for barcode, result in results.items():
            status = "SUCCESS" if result["success"] else "FAILED"
            error = (
                result.get("stderr", "").replace("\n", " ").strip()
                if not result["success"]
                else ""
            )
            patient_barcode_id = result.get("patient_barcode_id", "")
            f.write(f"{barcode}\t{status}\t{patient_barcode_id}\t{error}\n")

    logger.info(f"Summary report written to {summary_file}")
    return summary_file


def combine_summary_tsvs(output_dir, date_str):
    """
    Combine all *_summary.tsv files from telomerehunter2_results_YYYYMMDD/*/ subfolders into one file.
    """
    th2_results_dir = output_dir / f"telomerehunter2_results_{date_str}"
    combined_path = output_dir / f"combined_results_th2_{date_str}.tsv"

    # Look for summary files one level deeper: telomerehunter2_results_YYYYMMDD/*/*_summary.tsv
    summary_files = sorted(th2_results_dir.glob("*/*_summary.tsv"))
    if not summary_files:
        logger.warning(f"No summary TSV files found in {th2_results_dir} subfolders")
        return

    logger.info(
        f"Combining {len(summary_files)} summary TSV files into {combined_path}"
    )

    header_written = False
    with open(combined_path, "w") as outfile:
        for summary_file in summary_files:
            with open(summary_file, "r") as infile:
                lines = infile.readlines()
                if not lines:
                    continue
                # Write header only once
                if not header_written:
                    outfile.write(lines[0])
                    header_written = True
                # Write data lines (skip header)
                for line in lines[1:]:
                    outfile.write(line)
    logger.info(f"Combined summary written to {combined_path}")


def collect_existing_barcode_bams(bam_out_dir, patient_id, expected_barcodes=None):
    """
    Collects already split BAM files in the barcode folder and creates the barcode_files structure.
    Supports both {patient_id}_{barcode}.bam and {barcode}.bam file patterns.
    If expected_barcodes is given, checks that all are present and raises an error if not.
    """
    barcode_files = {}
    # Collect all BAM files matching both patterns
    all_bam_files = list(bam_out_dir.glob("*.bam"))
    for bam_file in all_bam_files:
        stem = bam_file.stem
        if stem.startswith(f"{patient_id}_"):
            safe_bc = stem[len(patient_id) + 1 :]
        else:
            safe_bc = stem
        barcode = (
            safe_bc  # Assumption: safe_barcode = barcode, adjust if mapping is needed
        )
        # Avoid duplicates (prefer patient_id_ prefix if both exist)
        if barcode not in barcode_files or stem.startswith(f"{patient_id}_"):
            barcode_files[barcode] = {
                "file": str(bam_file),
                "safe_barcode": safe_bc,
                "read_count": None,  # Read count unknown
            }
    # If expected_barcodes is provided, check completeness
    if expected_barcodes is not None:
        found = set(barcode_files.keys())
        expected = set(expected_barcodes)
        missing = expected - found
        logger.info(
            f"Found BAMs for {len(found & expected)}/{len(expected)} expected barcodes."
        )
        if missing:
            logger.error(
                f"Missing BAM files for {len(missing)} barcodes: {', '.join(sorted(missing))}"
            )
            raise RuntimeError(
                f"Not all expected barcode BAMs are present in {bam_out_dir}."
            )
    return barcode_files


def main():
    parser = argparse.ArgumentParser(
        description="Split sc-ATAC BAM by barcodes and run telomerehunter2",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run all steps (default)
  python %(prog)s input.bam ./output Patient1 cytobanding_hg38.txt

  # Only run telomerehunter2 on existing barcode BAMs
  python %(prog)s input.bam ./output Patient1 cytobanding_hg38.txt --steps th2

  # Only split and run th2 (use existing barcode file)
  python %(prog)s input.bam ./output Patient1 cytobanding_hg38.txt --steps split,th2

  # Only generate barcode file
  python %(prog)s input.bam ./output Patient1 cytobanding_hg38.txt --steps barcode
        """,
    )
    parser.add_argument(
        "-ibt", "--inputBamTumor", dest="input_bam", help="Input sc-ATAC BAM/CRAM file"
    )
    parser.add_argument("-o", "--output_dir", help="Output directory")
    parser.add_argument("-p", "--patient_id", help="Patient ID prefix")
    parser.add_argument(
        "-b", "--cytobanding_file", help="Cytobanding file for telomerehunter2"
    )
    parser.add_argument(
        "--max-barcodes",
        type=int,
        default=None,
        help="Maximum number of barcodes to process",
    )
    parser.add_argument(
        "--keep-bams", action="store_true", help="Keep individual barcode BAM files"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Enable verbose logging"
    )
    parser.add_argument(
        "--max-parallel",
        type=int,
        default=os.cpu_count(),
        help="Maximum number of parallel runs, default max number of cores",
    )
    parser.add_argument(
        "--min-reads-per-barcode",
        type=int,
        default=10000,
        help="Minimum number of reads per barcode to keep (default: 10000)",
    )
    parser.add_argument(
        "--steps",
        type=str,
        default="barcode,split,th2",
        help="Comma-separated list of steps to run: barcode,split,th2 (default: all)",
    )

    args = parser.parse_args()

    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Print max_parallel and available cores
    available_cores = multiprocessing.cpu_count()
    logger.info(
        f"Max parallel: {args.max_parallel} | Available cores: {available_cores}"
    )

    # Parse steps
    steps = [s.strip().lower() for s in args.steps.split(",") if s.strip()]
    run_barcode = "barcode" in steps
    run_split = "split" in steps
    run_th2 = "th2" in steps

    # Validate inputs
    if not os.path.exists(args.input_bam):
        logger.error(f"Input BAM file not found: {args.input_bam}")
        sys.exit(1)

    if not os.path.exists(args.cytobanding_file):
        logger.error(f"Cytobanding file not found: {args.cytobanding_file}")
        sys.exit(1)

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    barcode_txt_path = output_dir / "barcodes.txt"

    # Prepare subfolders
    bam_out_dir = output_dir / "bam_files_by_barcode"
    bam_out_dir.mkdir(parents=True, exist_ok=True)
    today_str = datetime.datetime.now().strftime("%Y%m%d")
    th2_results_dir = output_dir / f"telomerehunter2_results_{today_str}"
    th2_results_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Extract barcodes and write to text file
    if run_barcode:
        logger.info("=== Extracting barcodes and writing barcode list ===")
        barcode_dict = get_barcodes_from_bam(
            args.input_bam,
            args.max_barcodes,
            barcode_txt_path=barcode_txt_path,
            skip_barcode_file=False,
            min_reads_per_barcode=args.min_reads_per_barcode,
        )
        if not barcode_dict:
            logger.error("No barcodes found in input BAM")
            sys.exit(1)
        logger.info(f"{len(barcode_dict)} barcodes extracted")
        # If only barcode step is requested, exit after writing barcode file
        if steps == ["barcode"]:
            logger.info(
                "Only barcode step requested, exiting after barcode file generation."
            )
            sys.exit(0)
    else:
        if not barcode_txt_path.exists():
            logger.error(
                f"Barcode file {barcode_txt_path} does not exist but barcode step was skipped."
            )
            sys.exit(1)
        barcode_set = read_barcodes_from_file(barcode_txt_path, args.max_barcodes)
        if not barcode_set:
            logger.error("No barcodes found in existing barcode file")
            sys.exit(1)
        logger.info(f"{len(barcode_set)} barcodes loaded from existing barcode file")
        barcode_dict = {bc: None for bc in barcode_set}

    # Step 2: Sinto splitting
    if run_split:
        logger.info("=== Starting Sinto splitting ===")
        barcode_files = split_with_sinto(
            args.input_bam,
            output_dir,
            args.patient_id,
            barcode_dict,
            barcode_txt_path,
            max_parallel=args.max_parallel,
        )
        if not barcode_files:
            logger.error("Sinto splitting failed or no BAM files produced")
            sys.exit(1)
        logger.info(
            f"BAM successfully split into {len(barcode_files)} barcode-specific files"
        )
        # If only split step is requested, exit after splitting
        if steps == ["split"]:
            logger.info("Only split step requested, exiting after BAM splitting.")
            sys.exit(0)
    else:
        logger.info("Skipping BAM splitting, using existing barcode BAMs")
        if run_th2:
            try:
                barcode_files = collect_existing_barcode_bams(
                    bam_out_dir, args.patient_id, barcode_dict.keys()
                )
            except RuntimeError as e:
                logger.error(str(e))
                sys.exit(1)
            logger.info(f"{len(barcode_files)} barcode BAMs found")

    # Step 3: Telomerehunter2
    if run_th2:
        logger.info("=== Starting telomerehunter2 processing ===")
        results = process_barcodes_with_telomerehunter2(
            barcode_files,
            args.patient_id,
            args.cytobanding_file,
            th2_results_dir,
            cleanup_bams=not args.keep_bams,
            max_parallel=args.max_parallel,
        )
        summary_file = write_summary_report(
            results,
            args.patient_id,
            args.input_bam,
            th2_results_dir,
            args.cytobanding_file,
        )
        successful = sum(1 for r in results.values() if r["success"])
        logger.info("=== Processing finished ===")
        logger.info(f"Successfully processed: {successful}/{len(results)} barcodes")
        logger.info(f"Details: {summary_file}")
        if successful < len(results):
            logger.warning(f"{len(results) - successful} barcodes failed")
            logger.info("See summary file for error details")
        combine_summary_tsvs(output_dir, today_str)
        # If only th2 step is requested, exit after processing
        if steps == ["th2"]:
            logger.info(
                "Only th2 step requested, exiting after telomerehunter2 processing."
            )
            sys.exit(0)
    else:
        logger.info("Skipping telomerehunter2 step.")
        split_summary = output_dir / "split_summary.txt"
        with open(split_summary, "w") as f:
            f.write("=== BAM Splitting Summary ===\n")
            f.write(f"Input BAM: {args.input_bam}\n")
            f.write(f"Patient ID: {args.patient_id}\n")
            f.write(f"Total barcodes: {len(barcode_files)}\n")
            f.write(f"Output directory: {output_dir}\n\n")
            f.write("Barcode\tRead_Count\tOutput_File\n")
            for barcode, info in barcode_files.items():
                read_count = info.get("read_count", "NA")
                f.write(f"{barcode}\t{read_count}\t{info['file']}\n")
        logger.info(f"Split summary written to {split_summary}")


if __name__ == "__main__":
    main()
