import os
import re
import shutil
import subprocess
from ftplib import FTP
from urllib.parse import urlparse

import pysam


def download_file(url, file_name, local_path):
    """
    Download a bam file from the specified URL to the local path.

    Args:
    url (str): The URL of the file to download.
    file_name (str): The name of the file to be saved locally.
    local_path (str): The local directory where the file will be saved.

    Returns:
    str: The path to the downloaded BAM file.
    """

    os.makedirs(local_path, exist_ok=True)

    parsed_url = urlparse(url)
    ftp = FTP(
        parsed_url.netloc
    )  # https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam
    ftp.login()
    ftp.cwd(parsed_url.path)
    with open(os.path.join(local_path, file_name), "wb") as f:
        ftp.retrbinary("RETR " + file_name, f.write)
    ftp.quit()
    print("Download done")

    return os.path.join(local_path, file_name)


def get_file_type(filename):
    if filename.endswith(".bam"):
        return "BAM"
    elif filename.endswith(".cram"):
        return "CRAM"
    else:
        return None


def run_telomerehunter_package(
    bam_file_path,
    results_path,
    patient_name,
    bam_file_path_control=None,
    banding_file=None,
):
    """
    Run Telomerehunter2 on the specified BAM file via the commandline.

    Args:
    bam_file_path (str): The path to the input BAM file.
    results_path (str): The path to the directory where results will be saved.
    patient_name (str): The name of the patient or sample.
    bam_file_path_control (str) optional: The path to the input control BAM file.


    Returns:
    None
    """
    command = [
        "telomerehunter2" "-ibt",
        bam_file_path,
        "-o",
        results_path,
        "-p",
        patient_name,
    ]
    if bam_file_path_control:
        command += ["-ibc", bam_file_path_control]
    if banding_file:
        command += ["-b", banding_file]

    print(*command)

    # Run the command and capture its output
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, universal_newlines=True)

    # Iterate over the output line by line and print it in real-time
    for line in proc.stdout:
        print(
            line, end=""
        )  # The end='' argument is used to prevent adding extra newlines

    # Wait for the command to complete and capture the exit code
    exit_code = proc.wait()

    assert exit_code == 0, f"Command finished with error. Exit code: {exit_code}"


def run_telomerehunter_live(bam_file_path, results_path, patient_name, parameters=None):
    """
    Run Telomerehunter2 on the specified BAM file via the python interpreter (for debugging).

    Args:
    bam_file_path (str): The path to the input BAM file.
    results_path (str): The path to the directory where results will be saved.
    patient_name (str): The name of the patient or sample.
    bam_file_path_control (str) optional: The path to the input control BAM file.

    Returns:
    None
    """
    telomerehunter_path = os.path.join(
        os.getcwd(), "..", "src", "telomerehunter2", "telomerehunter2_main.py"
    )
    command = [
        "python",
        telomerehunter_path,
        "-ibt",
        bam_file_path,
        "-o",
        results_path,
        "-p",
        patient_name,
    ]
    if parameters:
        command += parameters

    print(*command)

    # Run the command and capture its output
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, universal_newlines=True)

    # Iterate over the output line by line and print it in real-time
    for line in proc.stdout:
        print(
            line, end=""
        )  # The end='' argument is used to prevent adding extra newlines

    # Wait for the command to complete
    proc.wait()


def run_telomerehunter_live_complete(parameters):
    telomerehunter_path = os.path.join(
        os.getcwd(), "telomerehunter2", "telomerehunter2_main.py"
    )
    command = [
        "python",
        telomerehunter_path,
    ]
    if parameters:
        command += parameters

    print(*command)

    # Run the command and capture its output
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, universal_newlines=True)

    # Iterate over the output line by line and print it in real-time
    for line in proc.stdout:
        print(
            line, end=""
        )  # The end='' argument is used to prevent adding extra newlines

    # Wait for the command to complete
    proc.wait()


def subsample_bam(input_path, format, fraction=0.2):
    """
    Subsample the input BAM/CRAM file to the specified number of reads.

    Args:
    input_path (str): The path to the input BAM file.
    fraction (float): The fraction of reads to keep

    Returns:
    str: The path to the subsampled BAM file.
    """
    if format == "BAM":
        subsampled_filename = re.sub(r"\.bam$", "_subsampled.bam", input_path)

        # Use pysam.view with '-s' option to subsample by fraction
        pysam.view(
            "-b",
            "-s",
            str(fraction),
            input_path,
            "-o",
            subsampled_filename,
            catch_stdout=False,
        )
    elif format == "CRAM":
        subsampled_filename = re.sub(r"\.cram$", "_subsampled.cram", input_path)

        # Use pysam.view with '-s' option to subsample by fraction
        pysam.view(
            "-C",
            "-s",
            str(fraction),
            input_path,
            "-o",
            subsampled_filename,
            catch_stdout=False,
        )
    else:
        raise TypeError("No subsetting possible due to format")

    print("Done subsampling")

    return subsampled_filename


def get_unmapped_bam(input_path, format):
    """
    Get unmapped reads file from BAM or CRAM input file.

    Args:
    input_path (str): The path to the input BAM file.

    Returns:
    str: The path to the subsampled BAM file.
    """
    if format == "BAM":
        unmapped_filename = re.sub(r"\.bam$", "_unmapped.bam", input_path)

        # Use pysam.view with '-s' option to subsample by fraction
        pysam.view(
            "-b", "-f 4", input_path, "-o", unmapped_filename, catch_stdout=False
        )
    elif format == "CRAM":
        unmapped_filename = re.sub(r"\.cram$", "_unmapped.cram", input_path)

        # Use pysam.view with '-s' option to subsample by fraction
        pysam.view(
            "-C", "-f 4", input_path, "-o", unmapped_filename, catch_stdout=False
        )
    else:
        raise TypeError("No unmapped filtering possible due to format")

    print("Done filtering unmapped reads.")

    return unmapped_filename


def create_example_bam_file(data_folder):
    # Define the sequence data
    sequences = [
        (
            "read1",
            "CGATGCTATCGGTTAGGGTTAGGGATAGGGGTAGTAGATTAGGGGTTAGGCGTAGTGCTAGCTTAGGGCGTACTTAGGGTTAGGG",
        ),
        (
            "read2",
            "CGATTACGATTACGTTAGGGTTAGGGTTAGGGTTGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTAGCGATCGAT",
        ),
        ("read3", "GATCGATCGATCGA"),
        ("read4", "TACGTACGTACGTA"),
    ]

    # Create a new BAM file
    bam_file_cre_path = f"{data_folder}/singletons_example.bam"
    bam_file_cre = pysam.AlignmentFile(
        bam_file_cre_path, "wb", reference_names=["chr1"], reference_lengths=[100]
    )

    # Add the reads to the BAM file
    for name, seq in sequences:
        read = pysam.AlignedSegment()
        read.query_name = name
        read.query_sequence = seq
        read.flag = 0
        read.reference_id = 0
        read.reference_start = 100
        read.mapping_quality = 60
        read.cigarstring = "{}M".format(len(seq))
        bam_file_cre.write(read)

    bam_file_cre.close()

    return bam_file_cre_path

def inspect_bam_file(count_reads=True):
    bam_path = "//results/tumor_banding/tumor_TelomerCnt_tumor_banding/tumor_banding_filtered_intrachromosomal.bam"

    with pysam.AlignmentFile(bam_path, "rb") as bam_file:
        if count_reads:

            # Initialize a counter
            read_count = 0

            # Iterate through all reads
            for _ in bam_file:
                read_count += 1

            # Print the total number of reads
            print(f"Total number of reads: {read_count}")

        else:
            for read in bam_file:
                print(f"Read name: {read.query_name}")
                print(f"Reference name: {read.reference_name}")
                print(f"Mapping quality: {read.mapping_quality}")
                print(f"Aligned position: {read.reference_start}")
                print(f"Is mapped: {read.is_mapped}")
                print(f"Is reverse strand: {read.is_reverse}")
                print(f"Is paired: {read.is_paired}")
                print(f"Is proper_paired: {read.is_proper_pair}")
                print(f"Sequence: {read.query_sequence}")
                print(f"Qualities: {read.query_qualities}")
                print(f"CIGAR string: {read.cigarstring}")
                print(f"Flags: {read.flag}")

                # Get optional tags
                for tag, value in read.get_tags():
                    print(f"Tag {tag}: {value}")
                print("------------------------------")


def run_telomerehunter2_sc(bam_file_path, extra_params=None):
    script_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
        "src",
        "telomerehunter2",
        "telomerehunter2_sc.py"
    )
    command = ["python", script_path, "-ibt", bam_file_path]
    if extra_params:
        command += extra_params
    print("Running:", " ".join(command))
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, universal_newlines=True)
    for line in proc.stdout:
        print(line, end="")
    exit_code = proc.wait()
    assert exit_code == 0, f"Command finished with error. Exit code: {exit_code}"


if __name__ == "__main__":
    # Get the directory path of the current script or file
    script_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    # Combine the script directory and data folder using os.path.join
    data_folder = os.path.join(script_dir, "data")

    # Specify testfile
    # file_name = "HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam"
    # file_name = "HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam"
    # file_name = "HG00096.unmapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam"
    # file_name = "HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam.cram"
    # file_name = "HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522_subsampled.bam"
    # file_name = "HG00096.combination11.bam"
    # file_name = "test.sorted.bam"
    # file_name = "HG00096.combination11.cram"
    # file_name = "atac_hgmm_1k_nextgem_possorted_bam_subsampled_10pct.bam"
    # file_name = "atac_pbmc_500_nextgem_possorted_bam.bam"
    file_name = "HG00097.chrom11.ILLUMINA.bwa.GBR.low_coverage.20130415.bam"

    bam_file_path = os.path.join(data_folder, file_name)
    file_format = get_file_type(file_name)

    # optional files
    control_bam = bam_file_path  # os.path.join(data_folder, "HG00097.chrom11.ILLUMINA.bwa.GBR.low_coverage.20130415.bam")
    banding_file = os.path.join(
        script_dir, "src", "telomerehunter2", "cytoband_files", "hg19_cytoBand.txt"
    )
    banding_file_38 = os.path.join(
        script_dir, "src", "telomerehunter2", "cytoband_files", "hg38_cytoBand.txt"
    )

    print("---------------------------")
    print(f"Executing on {os.getcwd()}")

    assert (
        shutil.which("telomerehunter2") is not None
    ), "Error: telomerehunter2 not found in PATH or is not executable"

    # Download if not exist yet # alternatively run wget
    # https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment/HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam
    if not os.path.exists(bam_file_path):
        print("Downloading testfile")
        download_file(
            "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/alignment",
            file_name,
            data_folder,
        )
    else:
        print("Testfile already present")

    # Subsample if fast testing is needed
    preprocessing = "None"
    print(f"Preprocessing: {preprocessing}")
    switch = {
        "Subsample": lambda: subsample_bam(bam_file_path, file_format, fraction=0.3),
        "Unmapped": lambda: get_unmapped_bam(bam_file_path, file_format),
        "None": lambda: bam_file_path,
    }
    bam_file_path_sub = switch.get(preprocessing, lambda: print("Invalid choice"))()

    assert os.path.exists(bam_file_path_sub), "Subsampled BAM file not found."

    results_path = os.path.join(script_dir, "results")
    os.makedirs(results_path, exist_ok=True)
    assert os.path.exists(results_path), "Results directory not created."

    print("Running Telomerehunter2")
    # testing the package version
    # run_telomerehunter_package(bam_file_path_sub, results_path, "package_test", banding_file=banding_file)

    # DEPRECATED DEVELOPER: testing the live dev version from python directly

    # only tumor
    # run_telomerehunter_live(bam_file_path_sub, results_path, "tumor_banding_parallel", parameters=["-b", banding_file])
    run_telomerehunter_live(bam_file_path_sub, results_path, "tumor_banding_parallel", parameters=["-b", banding_file, "-pno"])
    # run_telomerehunter_live(bam_file_path_sub, results_path, "tumor_banding_parallel_subsample", parameters=["-b", banding_file, "-pno", "--subsample", "0.2"])

    # tumor and control banding
    # run_telomerehunter_live(bam_file_path_sub, results_path, "tumor_control_banding",
    #                         parameters=["-b", banding_file, "-ibc", control_bam])

    # tumor and control banding subsampled
    # run_telomerehunter_live(
    #     bam_file_path_sub,
    #     results_path,
    #     "tumor_control_banding_subsampled",
    #     parameters=[
    #         "-b",
    #         banding_file,
    #         "-ibc",
    #         control_bam,
    #         "--subsample",
    #         "0.2",
    #         "-pno",
    #     ],
    # )

    # tumor banding fast mode
    # run_telomerehunter_live(bam_file_path_sub, results_path, "tumor_banding_fast", parameters=["-b", banding_file, "--fast"])

    # tumor without banding file
    # run_telomerehunter_live(bam_file_path_sub, results_path, "tumor_no_banding")

    # tumor flexible input repeats and hexamers
    # run_telomerehunter_live(bam_file_path_sub, results_path, "tumor_heptamers", parameters=["-r", "TTAGGGG", "TGAGGGG", "TCAGGGG", "TTGGGGG", "-bp", "21", "-rc", "TCAGGGG", "TGAGGGG", "TTGGGGG", "TTCGGGG", "TTTGGGG", "ATAGGGG", "CATGGGG", "CTAGGGG", "GTAGGGG", "TAAGGGG"])

    # run_telomerehunter2_sc(bam_file_path, extra_params=["-o", results_path, "-p", "sc_th2", "-b", banding_file, "--min-reads-per-barcode", "10000"])

    print("Done testing")
