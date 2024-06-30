import os

configfile: "config_test.yaml"
# Get the mode from the config file

mode = config["mode"]

# Rule for the final target
rule all:
    input:
        f"results/{mode}/alignment_crs.bam" if config["mode"] != "trans" else [f"results/{mode}/alignment_crs.bam", f"results/{mode}/alignment_guide.bam"]

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
        bwa index {input.reference_guide}
        bwa mem -B 100 -O 100 -E 100 -t {params.threads} {input.reference_guide} {input.fastq_guide} | samtools view -b > {output.bam_guide}
        
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
        bwa index {input.reference}
        bwa mem -t {params.threads} {input.reference} {input.FWD} {input.REV} | samtools view -b > {output.bam}
        """
         
# Rule to process files with a Python script
rule process_files:
    input:
        crs = f"results/{mode}/alignment_crs.bam",
        bc = config["input_files"]["bc"]
        guide = f"results/{mode}/alignment_guide.bam" if config["mode"] == "trans" else (config["input_files"]["guide"] if config["mode"] == "sc" else None),
    output:
        csv_gz = f"results/{mode}/counts_matrix.csv.gz"
    script:
        "scripts/process_files.py input.crs input.bc input."

# # # Rule to generate the HTML report using R Markdown
# # rule generate_report:
# #     input:
# #         csv_gz = "results/processed_data.csv.gz"
# #     output:
# #         html = "results/final_report.html"
# #     shell:
# #         """
# #         {load_modules()}
# #         Rscript -e "rmarkdown::render('scripts/generate_report.Rmd', params=list(input='{input.csv_gz}'), output_file='{output.html}')"
# #         """