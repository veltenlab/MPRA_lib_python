import os

configfile: "config.yaml"

# Rule for the final target
rule all:
    input:
        "results/final_report.html"

#Add here a rule to check, whether reference exists and create .fa if necessary
# rule check_and_convert_reference:
#     input:
#         reference = lambda wildcards: config[f'reference_{sample}']
#     output:
#         reference_fa = temp("data/reference_{sample}.fa")
#     run:
#         input_ref = input.reference
#         output_ref = output.reference_fa
#         if input_ref.endswith(".csv"):
#             shell(f"awk -F, '{{print \">{wildcards.sample}\"$1\"\\n\"$2}}' {input_ref} > {output_ref}")
#         else:
#             shell(f"cp {input_ref} {output_ref}")

# Rule to align the A FASTQ file using BWA mem
rule align_guide:
    input:
        fastq = config["input_files"]["guide"],
        reference = config["reference_guide"]
    output:
        bam = "results/alignment_guide.bam"
    shell:
        """
        module use /software/as/el7.2/EasyBuild/CRG/modules/all
        module load BWA/0.7.17
        module load SAMtools/1.10-GCC-9.3.0

        bwa index {input.reference}
        bwa mem -B 100 -O 100 -E 100 -t 1 {input.reference} {input.fastq} | samtools view -b > {output.bam}
        """

# Rule to align the B FASTQ file using BWA mem with specific parameters, only if mode is X
rule align_crs:
    input:
        fastq = config["input_files"]["crs"],
        reference = config["reference_crs"]
    output:
        bam = "results/sample_crs.bam"
    shell:
        """
        module use /software/as/el7.2/EasyBuild/CRG/modules/all
        module load BWA/0.7.17
        module load SAMtools/1.10-GCC-9.3.0

        bwa index {input.reference}
        bwa mem -t 1 {input.reference} {input.fastq} | samtools view -b > {output.bam}
        """
# Rule to process files with a Python script
rule process_files:
    input:
        preprocess_guide = "results/alignment_guide.bam" if config["mode"] == "trans" else config["input_files"]["guide"],
        preprocess_crs = "results/alignment_crs.bam" 
    output:
        csv_gz = "results/processed_data.csv.gz"
    script:
        "scripts/process_files.py"

# Rule to generate the HTML report using R Markdown
rule generate_report:
    input:
        csv_gz = "results/processed_data.csv.gz"
    output:
        html = "results/final_report.html"
    shell:
        """
        {load_modules()}
        Rscript -e "rmarkdown::render('scripts/generate_report.Rmd', params=list(input='{input.csv_gz}'), output_file='{output.html}')"
        """

# Conditional execution of align_B and subsequent steps based on the mode
if config["mode"] == "trans":
    include: ...
else:
    include: ...
