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
    run:
        reference = input.reference
        output = output.reference
        if reference.endswith(".csv"):
            shell(f"awk -F, '{{print \">{wildcards.sample}\"$1\"\\n\"$2}}' {reference} > {output}")
        else:
            shell(f"cp {reference} {output}")

# Rule to align the A FASTQ file using BWA mem
rule align_guide:
    input:
        fastq = config["input_files"]["guide"],
        reference = "data/reference_guide.fa"
    output:
        bam = "results/alignment_guide.bam"
    # conda:
    #     "MPRA_env.yml"
    params:
        threads=config["threads_bwa"]
    shell:
        """
        bwa index {input.reference}
        bwa mem -B 100 -O 100 -E 100 -t {params.threads} {input.reference} {input.fastq} | samtools view -b > {output.bam}
        """

# Rule to align the B FASTQ file using BWA mem with specific parameters, only if mode is X
rule align_crs:
    input:
        fastq = config["input_files"]["crs"],
        reference = "data/reference_crs.fa"   
    output:
        bam = "results/alignment_crs.bam"
    # conda:
    #     "MPRA_env.yml"
    params:
        threads = config["threads_bwa"]
    shell:
        """
        bwa index {input.reference}
        bwa mem -t {params.threads} {input.reference} {input.fastq} | samtools view -b > {output.bam}
        """

# rule align_crs_paired:
#     input:
#         reference = reference = "data/reference_crs.fa"

#         FWD = config["input_files"]["crs"],
#         REV = config["input_files"]["crs_rev"]
#     output:
#         bam = "results/alignment_crs.bam"
#     conda:
#        "MPRA_processing.yml"   
#     params:
#         threads=config["threads_bwa"]   
#     log:
#         OutDir+"/log_files/bwa_map/{Design}_{FastqDir}_{LibraryName}_bwa_map.log"
#     shell:
#         "bwa mem -a -t {params.threads} {input.design} {input.FWD} {input.REV} | samtools view -b > {output.bam} " 
#         "2> {log}" 

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

# # Conditional execution of align_B and subsequent steps based on the mode
# if config["mode"] == "trans":
#     include: ...
# else:
#     include: ...
