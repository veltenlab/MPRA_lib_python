"""This file contains functions that are used in the GRE-BC mapping script"""
import gzip
import pysam
import re


def count_reads(file_path):
    """Counts the number of reads in a BAM or FASTQ file, whether gzipped or not
    
    Args:
        file_path (str): path to the BAM or FASTQ file.
    
    Returns:
        (int): the number of reads.
    """
    # Determine file type and if it is gzipped
    is_bam = file_path.endswith('.bam') or file_path.endswith('.bam.gz')
    is_gzipped = file_path.endswith('.gz')
    
    if is_bam:
        mode = 'rb' if is_gzipped else 'r'
        with pysam.AlignmentFile(file_path, mode) as bam_file:
            return sum(1 for _ in bam_file)
    else:
        read_count = 0
        if is_gzipped:
            with gzip.open(file_path, 'rt') as fastq_file:
                for line in fastq_file:
                    if line.startswith('@') and not line.startswith('@+'):
                        read_count += 1
        else:
            with open(file_path, 'r') as fastq_file:
                for line in fastq_file:
                    if line.startswith('@') and not line.startswith('@+'):
                        read_count += 1
        return read_count

# Old cunt_reads function

# def count_reads(bam_filename, bc_filename, guide_filename):
#     bam_file_read = pysam.AlignmentFile(bam_filename, 'rb')
#     alignment_count = sum(1 for _ in bam_file_read)
#     print("Bam alignments count: ", alignment_count)
    
#     with gzip.open(bc_filename, 'rt') if bc_filename.endswith('.gz') else open(bc_filename, 'r') as f:
#         record_count = 0
#         for line in f:
#             if line.startswith("@"):
#                 record_count += 1
#         print("BC fastq reads count: ", record_count)

#     with gzip.open(guide_filename, 'rt') if guide_filename.endswith('.gz') else open(guide_filename, 'r') as f:
#         record_count = 0
#         for line in f:
#             if line.startswith("@"):
#                 record_count += 1
#         print("Guide fastq reads count: ", record_count)

# Defining a function to convert cigar string to score

def cigar_to_score(cigar_tuples):
    """Extract numbers from tuples where the first element is 0 = Match"""
    match_numbers = [num for op, num in cigar_tuples if op == 0]
    total_match_numbers = sum(match_numbers)
    return total_match_numbers

def filter_bam_file(input_bam_path, output_bam_path):
    """Function to filter the BAM file to keep only the best alignment per read."""
    with pysam.AlignmentFile(input_bam_path, "rb") as input_bam, pysam.AlignmentFile(output_bam_path, "wb", template=input_bam) as output_bam:
        best_alignments = {}
        
        for alignment in input_bam:
            readname = alignment.query_name
            if readname not in best_alignments:
                best_alignments[alignment.query_name] = alignment
            elif cigar_to_score(alignment.cigar) > cigar_to_score(best_alignments[readname].cigar):
                best_alignments[alignment.query_name] = alignment

        for alignment in best_alignments.values():
            output_bam.write(alignment)

def get_next_barcode(guide_file, bc_file = None, is_bam = False):
    """Function to generate next barcode

    Args:
        bc_file (fastq): fastq file with GFP barcodes seq
        guide_file (fastq): fastq file with gRNA barcodes seq

    Returns:
        out (str): concantenated GFP barcode with the first 15 bp of the gRNA barcode
    """
    readname = None
    barcode = None
    guide = None

    if is_bam:
        bamfile = pysam.AlignmentFile(guide_file, "rb")
        try:
            alignment = next(bamfile)
            readname = alignment.query_name
            barcode = alignment.query_sequence
        except StopIteration:
            return None  # Indicate that EOF was reached
        bamfile.close()

        for i in range(1, 5):
            bcline = bc_file.readline().rstrip('\n')
            guideline = guide_file.readline().rstrip('\n')

            # Check for EOF in either file
            if not bcline or not guideline:
                return None  # Indicate that EOF was reached or one file is shorter than the other

            if i == 1:
                readname = re.sub(r'^@(\S+)\s.+', r'\1', bcline)
            elif i == 2:
                barcode = bcline
                guide = guideline

        out = [readname, barcode + guide]
    
    if not is_bam:
        for i in range(1, 5):
            bcline = bc_file.readline().rstrip('\n')
            guideline = guide_file.readline().rstrip('\n')

            # Check for EOF in either file
            if not bcline or not guideline:
                return None  # Indicate that EOF was reached or one file is shorter than the other

            if i == 1:
                readname = re.sub(r'^@(\S+)\s.+', r'\1', bcline)
            elif i == 2:
                barcode = bcline
                guide = guideline

        out = [readname, barcode + guide]
    # print("new_fastq_barcode: ", out)
    return out

def concat_sequences(fastq_file_path, second_file_path, is_bam=False):
    """Concatenate sequences from a FastQ file and a second file which can be either FastQ or BAM.

    Args:
        fastq_file_path (str): Path to the FastQ file.
        second_file_path (str): Path to the second file (either FastQ or BAM).
        is_bam (bool): Indicates whether the second file is a BAM file.

    Returns:
        list of tuple: List of tuples where each tuple contains the read name and the concatenated sequence.
    """
    results = []

    with open(fastq_file_path, 'r') as fastq_file:
        if is_bam:
            second_file = pysam.AlignmentFile(second_file_path, "rb")
        else:
            second_file = open(second_file_path, 'r')

        try:
            while True:
                # Read the next sequence block from the FastQ file
                fastq_file.readline().strip()
                fastq_seq = fastq_file.readline().strip()
                fastq_file.readline()  # skip '+' line
                fastq_file.readline()  # skip quality score line

                if not fastq_seq:
                    break  # End of file or empty line

                if is_bam:
                    try:
                        alignment = next(second_file)
                        second_seq = alignment.query_sequence
                    except StopIteration:
                        break  # No more reads in BAM
                else:
                    second_file.readline().strip()
                    second_seq = second_file.readline().strip()
                    second_file.readline()  # skip '+' line
                    second_file.readline()  # skip quality score line

                    if not second_seq:
                        break  # End of file or empty line

                # Concatenate the sequences
                results.append(fastq_seq + second_seq)

        finally:
            if is_bam:
                second_file.close()
            else:
                second_file.close()

    return results



def format_number(value):
    try:
        float_value = float(value)
        if float_value.is_integer():
            return int(float_value)
        else:
            return "%.12f" % float_value
    except ValueError:
        # Return the original value if it's not a number
        return value
    



    
################################################################
### Optional functions:
def hamming_distance(string1, string2):
    if len(string1) != 30 or len(string2) != 30:
        return 100  # Return a high number for unequal lengths
    
    distance = 0
    for ch1, ch2 in zip(string1, string2):
        if ch1 != ch2:
            distance += 1
            if distance > 3:
                break  # Stop calculation once distance exceeds 3
    return distance

def safe_concatenate(str1, str2):
    # Convert None to an empty string for both inputs before concatenating
    str1 = str1 if str1 is not None else ""
    str2 = str2 if str2 is not None else ""
    return str1 + str2

def concat_sequences(fastq_file_path, second_file_path, is_bam=False):
    """Concatenate sequences from a FastQ file and a second file which can be either FastQ or BAM.

    Args:
        fastq_file_path (str): Path to the FastQ file.
        second_file_path (str): Path to the second file (either FastQ or BAM).
        is_bam (bool): Indicates whether the second file is a BAM file.

    Returns:
        list of tuple: List of tuples where each tuple contains the read name and the concatenated sequence.
    """
    results = []

    line_count = 0

    with gzip.open(fastq_file_path, 'rt') as fastq_file:
        if is_bam:
            second_file = pysam.AlignmentFile(second_file_path, "rb")
        else:
            second_file = gzip.open(second_file_path, 'rt')

        try:
            while True:
                # Read the next sequence block from the FastQ file
                fastq_file.readline()
                fastq_seq = fastq_file.readline().strip()
                fastq_file.readline()  # skip '+' line
                fastq_file.readline()  # skip quality score line

                # line_count += 1
            
                # # Break out of the loop if line count reaches 110
                # if line_count >= 10:
                #     break

                if not fastq_seq:
                    break  # End of file or empty line


                if is_bam:
                    try:
                        alignment = next(second_file)
                        second_seq = second_file.get_reference_name(alignment.reference_id)
                    except StopIteration:
                        break  # No more reads in BAM
                else:
                    second_file.readline()
                    second_seq = second_file.readline().strip()
                    second_file.readline()  # skip '+' line
                    second_file.readline()  # skip quality score line

                    if not second_seq:
                        break  # End of file or empty line

                # Concatenate the sequences
                results.append(safe_concatenate(fastq_seq,second_seq))

        finally:
            if is_bam:
                second_file.close()
            else:
                second_file.close()

    return results



