"""This file contains functions that are used in the GRE-BC mapping script"""
import gzip
import pysam
import re
import Levenshtein

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

def cigar_to_score(cigar_tuples):
    """Extract numbers from tuples where the first element is 0 = Match"""
    match_numbers = [num for op, num in cigar_tuples if op == 0]
    total_match_numbers = sum(match_numbers)
    return total_match_numbers

def get_next_barcode(bc_file, guide_file = None, mode = "sc"):
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

    if mode == "trans":
        try:
            alignment = next(guide_file)
            readname = alignment.query_name
            barcode = alignment.query_sequence
            guide = guide_file.get_reference_name(alignment.reference_id)
        except StopIteration:
            return None  # Indicate that EOF was reached

        for i in range(1, 5):
            bcline = bc_file.readline().rstrip('\n')
            # guideline = guide_file.readline().rstrip('\n')

            # Check for EOF in either file
            if not bcline:
                return None  # Indicate that EOF was reached or one file is shorter than the other

            if i == 1:
                readname = re.sub(r'^@(\S+)\s.+', r'\1', bcline)
            elif i == 2:
                barcode = bcline
        # Barcode and guide sequences are separated by the vertical line to be able to separate them later
        out = [readname, f"{barcode}|{guide}"]
    
    elif mode == "sc":
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
        out = [readname, f"{barcode}|{guide}"]
    else:
        for i in range(1, 5):
            bcline = bc_file.readline().rstrip('\n')
            # Check for EOF in either file
            if not bcline:
                return None  # Indicate that EOF was reached or one file is shorter than the other
            if i == 1:
                readname = re.sub(r'^@(\S+)\s.+', r'\1', bcline)
            elif i == 2:
                barcode = bcline
                out = [readname, f"{barcode}|"] 
    return out

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

def correct_barcodes(BC_CRS, CRS_BC, threshold = 1, filter_ones = True):
    """Function to merge low-presented barcodes with the similar ones, which are better presented

    Args:
        BC_CRS (dic): dictionary of barcodes and corresponded CRSs with read counts
        CRS_BC (dic): dictionary of CRSs with corresponded matched barcodes
        threshold (int): the threshold value for the max accepted Levenstein distance for merge, default 1
        filter_ones (bool): specifies, whether the filtering step for deleting all one-reads for CRSs with more than 1000 barcodes should be applied

    Returns:
        BC_CRS_fixed (dic): shorter barcodes dictionary after merging step
        total_mapped_reads (int): summ statistics about total mapped reads
    """
    BC_CRS_fixed = {}
    nstep = 0
    ncrs = len(CRS_BC)
    removed_oneread = 0
    total_mapped_reads = 0

    for crs, barcodes in CRS_BC.items():
        nstep += 1
        BARCODES = list(barcodes)  # All barcodes mapping to a given CRS

        # Calculate total reads
        totreads = sum(BC_CRS[bc][2] for bc in BARCODES)
        total_mapped_reads += totreads

        # Filter (optional)
        if filter_ones == True:
            if len(BARCODES) > 1000:
                removed_oneread += len([bc for bc in BARCODES if BC_CRS[bc][2] == 1])
                BARCODES = [bc for bc in BARCODES if BC_CRS[bc][2] > 1]

        # Sort BARCODES by number of reads
        BARCODES.sort(key=lambda bc: BC_CRS[bc][2])

        while BARCODES:
            this_bc = BARCODES.pop(0)  # Start with the barcode with the least reads

            matched = False
            for target_bc in reversed(BARCODES):  # Compare with the barcode with the most reads first
                if Levenshtein.distance(this_bc, target_bc, weights=(10,10,1)) <= threshold:  # Adjust the threshold as needed
                    if target_bc in BC_CRS_fixed:
                        BC_CRS_fixed[target_bc][2] += BC_CRS[this_bc][2]  # Add read counts
                        BC_CRS_fixed[target_bc][3] += BC_CRS[this_bc][3]  # Add deviant read counts
                        BC_CRS_fixed[target_bc][1].extend(BC_CRS[this_bc][1])  # Append mapping statistics
                    else:
                        BC_CRS_fixed[target_bc] = BC_CRS[this_bc][:]  # Copy it over if not matched
                    matched = True
                    break

            if not matched:
                if this_bc in BC_CRS_fixed:
                    BC_CRS_fixed[this_bc][2] += BC_CRS[this_bc][2]  # Add read counts
                    BC_CRS_fixed[this_bc][3] += BC_CRS[this_bc][3]  # Add deviant read counts
                    BC_CRS_fixed[this_bc][1].extend(BC_CRS[this_bc][1])  # Append mapping statistics
                else:
                    BC_CRS_fixed[this_bc] = BC_CRS[this_bc][:]  # Copy it over if not matched

    return BC_CRS_fixed, total_mapped_reads