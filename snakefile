import os

configfile: "config_test.yaml"

# Rule for the final target
rule all:
    input:
        # "results/final_report.html"
        "results/alignment_crs.bam"

#Add here a rule to check, whether reference exists and create .fa if necessary
rule check_and_convert_reference:
    input:
        reference = lambda wildcards: config['references'][wildcards.sample]
    output:
        reference = "data/{sample}.fa"
    always_run:
        True
    run:
        reference = input.reference
        output = output.reference
        if reference.endswith(".csv"):
            shell(f"awk -F, '{{print \">{wildcards.sample}\"$1\"\\n\"$2}}' {reference} > {output}")
        else:
            shell(f"cp -f {reference} {output}")

rule align_guide:
    input:
        fastq = config["input_files"]["guide"],
        reference = "data/reference_guide.fa"
    output:
        bam = "results/alignment_guide.bam"
    params:
        threads=config["threads_bwa"]
    shell:
        """
        bwa index {input.reference}
        bwa mem -B 100 -O 100 -E 100 -t {params.threads} {input.reference} {input.fastq} | samtools view -b > {output.bam}
        """

rule align_crs:
    input:
        fastq = config["input_files"]["crs"],
        reference = "data/reference_crs.fa"   
    output:
        bam = "results/alignment_crs.bam"
    params:
        threads = config["threads_bwa"]
    shell:
        """
        bwa index {input.reference}
        bwa mem -t {params.threads} {input.reference} {input.fastq} | samtools view -b > {output.bam}
        """

rule align_crs_paired:
    input:
        reference = reference = "data/reference_crs.fa"

        FWD = config["input_files"]["crs"],
        REV = config["input_files"]["crs_rev"]
    output:
        bam = "results/alignment_crs.bam" 
    params:
        threads=config["threads_bwa"]   
    log:
        OutDir+"/log_files/bwa_map/{Design}_{FastqDir}_{LibraryName}_bwa_map.log"
    shell:
        "bwa mem -a -t {params.threads} {input.design} {input.FWD} {input.REV} | samtools view -b > {output.bam} " 
        "2> {log}" 

# # Rule to process files with a Python script
# rule process_files:
#     input:
#         preprocess_guide = "results/alignment_guide.bam" if config["mode"] == "trans" else config["input_files"]["guide"],
#         preprocess_crs = "results/alignment_crs.bam" 
#     output:
#         csv_gz = "results/processed_data.csv.gz"
#     script:
#         "scripts/process_files.py"

# # Rule to generate the HTML report using R Markdown
# rule generate_report:
#     input:
#         csv_gz = "results/processed_data.csv.gz"
#     output:
#         html = "results/final_report.html"
#     shell:
#         """
#         {load_modules()}
#         Rscript -e "rmarkdown::render('scripts/generate_report.Rmd', params=list(input='{input.csv_gz}'), output_file='{output.html}')"
#         """

# Function to dynamically select the rules for step 1 based on mode
def get_alignment_rules(mode):
    if mode == "trans":
        return ["align_guide", "align_crs"]
    elif mode == "ss":
        return ["align_crs"]
    elif mode == "bulk":
        return ["align_crs_paired"]
    else:
        raise ValueError(f"Unknown mode: {mode}")

# Get the mode from the config file
mode = config["mode"]

# Dynamic selection of step 1 rules
rule alignment:
    output:
        ["results/alignment_guide.bam", "results/alignment_crs.bam"] if mode == "trans" else ["results/alignment_crs.bam"]
    run:
        selected_rules = get_alignment_rules(mode)
        for rule in selected_rules:
            shell(f"snakemake {rule}")

# Final rule to run the entire workflow
rule all:
    input:
        "final_output.txt"