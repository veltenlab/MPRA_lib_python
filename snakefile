"""
This is a snakemake pipeline file with rules for the sequencing processing
It consists of 3 main steps: alignment, association and vizualisation
"""

from datetime import datetime

mode = config["mode"]
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

# Rule for the final target, I specified here both outputs from alignment and association steps
rule all:
    input:
        f"results/{mode}/alignment_crs.bam" if config["mode"] != "trans" else [f"results/{mode}/alignment_crs.bam", f"results/{mode}/alignment_guide.bam"],
        f"results/{mode}/counts_matrix_{timestamp}.csv.gz",
        # f"results/{mode}/report_{timestamp}.html"

# Rule to check if reference exists and create .fa if necessary
rule check_and_convert_reference:
    input:
        reference=lambda wildcards: config['references'][wildcards.sample]
    output:
        reference="data/{sample}.fa"
    run:
        reference = input.reference
        output = output.reference
        if reference.endswith(".csv"):
            shell(r"""awk -F, '{{print ">"$1"\n"$2}}' {input} > {output}""")
        else:
            shell("echo 'Preparing reference file'")
            shell(f"cp -f {reference} {output}")


# Rule to align guide
rule align_trans:
    input:
        fastq_guide = config["input_files"]["guide"],
        reference_guide = "data/reference_guide.fa",
        
        fastq_crs = config["input_files"]["crs"],
        reference_crs = "data/reference_crs.fa"
    output:
        bam_guide = "results/trans/alignment_guide.bam",
        bam = "results/trans/alignment_crs.bam"
    params:
        threads = config["threads_bwa"]
    shell:
        """
        echo "Running trans alignment"

        bwa index {input.reference_guide}
        bwa aln -n 0 -o 0 -e 0 {input.reference_guide} {input.fastq_guide} | bwa samse {input.reference_guide} - {input.fastq_guide} | samtools view -b > {output.bam_guide}


        bwa index {input.reference_crs}
        bwa mem -t {params.threads} {input.reference_crs} {input.fastq_crs} | samtools view -b > {output.bam}
        """

# Rule to align crs
rule align_sc:
    input:
        fastq = config["input_files"]["crs"],
        reference = "data/reference_crs.fa"
    output:
        bam = "results/sc/alignment_crs.bam"
    params:
        threads = config["threads_bwa"]
    shell:
        """
        echo "Running sc alignment"

        bwa index {input.reference}
        bwa mem -t {params.threads} {input.reference} {input.fastq} | samtools view -b > {output.bam}
        """

# Rule to align crs with paired-end reads
rule align_bulk:
    input:
        reference = "data/reference_crs.fa",
        FWD = config["input_files"]["crs"],
        REV = config["input_files"]["crs_paired"]
    output:
        bam = "results/bulk/alignment_crs.bam"
    params:
        threads = config["threads_bwa"]
    shell:
        """
        echo "Running bulk alignment"

        bwa index {input.reference}
        bwa mem -t {params.threads} {input.reference} {input.FWD} {input.REV} | samtools view -b > {output.bam}
        """

def get_guide_flag(wildcards, input):
    if input.get("guide"):
        return f"--guide_file {input.guide}"
    return ""

# Rule to process files with Python script and create an association library
rule process_files:
    input:
        crs = f"results/{mode}/alignment_crs.bam",
        bc = config["input_files"]["bc"],
        guide = f"results/{mode}/alignment_guide.bam" if mode == "trans" else (config["input_files"]["guide"] if mode == "sc" else [])
        # guide = f"results/{mode}/alignment_guide.bam" if mode == "trans" else config["input_files"]["guide"]
    output:
        csv = f"results/{mode}/counts_matrix_{timestamp}.csv.gz"
    params:
        threshold = config["bc_corr_threshold"],
        filter_bc = config["bc_corr_filtering"],
        # guide_flag = f"--guide_file {input.guide}" if input.guide else ""
        guide_flag = get_guide_flag
    shell:
        """
        # If mode = bulk, create a temporary empty fastq.gz file
        if [ "{mode}" == "bulk" ]; then
            temp_guide=$(mktemp --suffix=.fastq.gz)
            gzip -c /dev/null > "$temp_guide"
        else
            temp_guide={input.guide}
        fi

        python scripts/association.py --mode {mode} \
                          --crs_file {input.crs} \
                          --bc_file {input.bc} \
                          --guide_file "$temp_guide" \
                          --outfile {output.csv} \
                          --bc_corr_threshold {params.threshold} \
                          --bc_corr_filter {params.filter_bc}

        """

# # Rule to generate the HTML report using R Markdown
# rule generate_report:
#     input:
#         rmd = "scripts/generate_report.Rmd",
#         csv = f"results/{mode}/counts_matrix_{timestamp}.csv.gz"
#     output:
#         html = f"results/{mode}/report_{timestamp}.html"
#     shell:
#         """
#         Rscript -e "rmarkdown::render(input = '{input.rmd}', output_file = '{output.html}', params = list(data_file = '{input.csv}'))"
#         """