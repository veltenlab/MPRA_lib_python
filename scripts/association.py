"""
This is a script that takes sequencing and aligned files and aligns matches them together, creating an association library
This library will be later used to identify which barcode corresponds to which tested CRS
"""

import gzip
import pysam
import sys

from collections import defaultdict

from utils import *

# Check if at least one argument is provided
if len(sys.argv) != 8:
    print("Usage: association.py <mode> <bc_file> <guide_file> <crs_file> <outfile> <bc_corr_threshold> <bc_corr_filter>", file=sys.stderr)
    sys.exit(1)

# Save the arguments given from argv
mode, crs_file_path, bc_file_path, guide_file_path, outfile, bc_corr_threshold, bc_corr_filter = sys.argv[1:]
print(bc_corr_threshold, type(bc_corr_threshold))
print(bc_corr_filter, type(bc_corr_filter))

BC_CRS_raw = defaultdict(dict)
counter = 0
maps_score = {}
maps_count_fwd = {}
maps_count_rev = {}


#################################################
####      Mapping barcodes and alignments    ####
#################################################
# Create the first BC_CRS_raw dictionary by reading bc and crs files, then guide file mode-specific
# BC_CRS_raw = {
    
#     'barcode1': {'CRS1': [[score1, score2, ...], read_counts],
#                 'CRS10': [[score1, score2, ...], read_counts]},
                
#     'barcode2': {'CRS2': [[score1, score2, ...], read_counts]},

#     ...}

with gzip.open(bc_file_path, 'rt') as bc_file, \
        pysam.AlignmentFile(crs_file_path, 'rb') as crs_file:
        if mode == "trans":
            guide_file = pysam.AlignmentFile(guide_file_path, 'rb')
        elif mode == "sc":
            guide_file = gzip.open(guide_file_path, 'rt')               
        try:
            # Generate the bc sequence
            readname_fastq, bcline_fastq = get_next_barcode(guide_file, bc_file, mode)

            if readname_fastq is None:
                print("No barcodes found", file=sys.stderr)
                sys.exit(1)

            for bamline in crs_file:

                bits = format(bamline.flag, '05b')[::-1]  # Parse flags
                bam_ref = crs_file.get_reference_name(bamline.reference_id)

                # If there is more than one alignment for a sequence, check for orientation and amount of reads
                if bamline.query_name != readname_fastq:
                    # For bulk mode, accept only reads with 1 fow and 1 rev aligments
                    if mode == "bulk":
                        considered = [k for k in maps_score if maps_count_fwd.get(k, 0) == 1 and maps_count_rev.get(k, 0) == 1]
                    # For other modes, accept reads with only 1 fow alignment and 0 rev
                    else:
                        considered = [k for k in maps_score if maps_count_fwd.get(k, 0) == 0 and maps_count_rev.get(k, 0) == 1]
                
                    if len(considered) == 1:
                        selectedmap = considered[0]
                        score = maps_score[selectedmap]
                    elif len(considered) > 1:
                        sorted_considered = sorted(considered, key=lambda x: maps_score[x], reverse=True)
                        selectedmap = sorted_considered[0]
                        score = maps_score[selectedmap]

                    if considered and selectedmap != "*":
                        if bcline_fastq in BC_CRS_raw:
                            if selectedmap in BC_CRS_raw[bcline_fastq]:
                                BC_CRS_raw[bcline_fastq][selectedmap][0].append(score)
                                BC_CRS_raw[bcline_fastq][selectedmap][1] += 1
                            else:
                                BC_CRS_raw[bcline_fastq][selectedmap] = [[score], 1]
                        else:
                            BC_CRS_raw[bcline_fastq] = {selectedmap: [[score], 1]}

                    maps_score = {}
                    maps_count_fwd = {}
                    maps_count_rev = {}

                    readname_fastq, bcline_fastq = get_next_barcode(guide_file, bc_file, mode)

                    if readname_fastq is None:
                        print("No barcodes found", file=sys.stderr)
                        sys.exit(1)

                    counter += 1
                    
                    if bamline.query_name != readname_fastq:
                        print("Something is wrong", file=sys.stderr)
                        sys.exit(1)

                if bam_ref in maps_score:
                    maps_score[bam_ref] += cigar_to_score(bamline.cigar)
                    if bits[4] == "1":
                        maps_count_rev[bam_ref] += 1
                    else:
                        maps_count_fwd[bam_ref] += 1
                else:
                    maps_score[bam_ref] = cigar_to_score(bamline.cigar)
                    maps_count_rev[bam_ref] = int(bits[4])
                    maps_count_fwd[bam_ref] = 1 - int(bits[4])
        finally:
        # Ensure the guide_file is properly closed
            guide_file.close()

#############################################
####       Clean up multiple assignments#####
#############################################
#There are some barcodes that get assigned to multiple CRS. actually this happens a lot. In these cases count how often this happens but assign the barcode to the
# CRS that is supported with more reads, buty report the number of reads supporting a different assignment. This does not remove anything 
# BC_CRS = {
    
#     'barcode1': ['CRS1', [score1, score2, ...], read_counts, other_CRSs_read_counts],

#     'barcode2': ['CRS2', [score3, score4, ...], read_counts, other_CRSs_read_counts],

#     ...}

total_reads = 0
BC_CRS = defaultdict(list)
for barcode, crs in BC_CRS_raw.items():
    total_reads += sum(count for _, (score, count) in crs.items())
    sorted_crs = sorted(crs, key=lambda x: crs[x][1], reverse=True)
    current = sorted_crs[0]

    BC_CRS[barcode] = [current, crs[current][0], crs[current][1], 0]
    for c in sorted_crs[1:]:
        BC_CRS[barcode][3] += crs[c][1]
print(f"Cleaned dictionary BC_CRS with the length {len(BC_CRS)} was created")

# Revert dictionary to the crs-bc association from bc-crs
CRS_BC = defaultdict(list)

# Create a hash that maps each CRS to a list of barcodes (from BC:CRS zu CRS:BC)
for barcode, crs in BC_CRS.items():
    CRS_BC[crs[0]].append(barcode)
print("Completed creating CRS_BC")
print("CRS_BC length: :", len(CRS_BC))


#############################################
####        Run barcode correction      #####
#############################################
# If two barcodes map to the same CRS and are similar, add the less abundant barcode to the more abundant one
# Rank barcodes per CRS by the read counts and start with comparison of barcode with the lowest to one with the highest read counts

BC_CRS_fixed, total_mapped_reads = correct_barcodes(BC_CRS, CRS_BC, int(bc_corr_threshold), bc_corr_filter)

print(f"Total mapped reads: {total_mapped_reads}")
print(f"Length of BC_CRS: {len(BC_CRS)}")
print(f"Length of BC_CRS_fixed: {len(BC_CRS_fixed)}")


#############################################
####         Output results to csv       ####
#############################################
# Save the information about barcodes and corresponding crs sequences in a csv.gz file
with gzip.open(outfile, 'wt') as out:
    # Write the header
    out.write("BARCODE\tCRS\tREADS\tDEVIANTREADS\tMEANMATCHES\tMINMATCHES\tMAXMATCHES\n")
    count = 0
    # Iterate over each barcode in BC_CRS_fixed dictionary
    for barcode, crs in BC_CRS_fixed.items():
        # Calculate mean, min, and max matches
        mean_matches = format_number(sum(crs[1]) / len(crs[1])) if crs[1] else 0
        min_matches = min(crs[1]) if crs[1] else 0
        max_matches = max(crs[1]) if crs[1] else 0
        
        # Create the line to write to the output file
        line = f"{barcode}\t{crs[0]}\t{crs[2]}\t{crs[3]}\t{mean_matches}\t{min_matches}\t{max_matches}\n"
        # Write the line to the output file
        out.write(line)
