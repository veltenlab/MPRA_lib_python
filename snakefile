
# logic: all needs crs.bam > 
 
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
            shell(r"""awk -F, '{{print ">"$1"\\n"$2}}' {input} > {output}""")
        else:
            shell(f"cp -f {reference} {output}")
            
            
            # awk -F, '{print ">"$1"\n"$2}' {wildcards.sample} > {output}
            # awk -F, '{{print \">{wildcards.sample}\\n\"$2}}' {reference} > {output}

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
        REV = config["input_files"]["crs_rev"]
    output:
        bam = "results/bulk/alignment_crs.bam"
    params:
        threads = config["threads_bwa"]
    shell:
        """
        bwa index {input.reference}
        bwa mem -t {params.threads} {input.reference} {input.FWD} {input.REV} | samtools view -b > {output.bam}
        """

# # Function to dynamically select the rules for alignment based on mode
# def get_alignment_rules(mode):
#     if mode == "trans":
#         return ["align_guide", "align_crs"]
#     elif mode == "ss":
#         return ["align_crs"]
#     elif mode == "bulk":
#         return ["align_crs_paired"]
#     else:
#         raise ValueError(f"Unknown mode: {mode}")



# Dynamic selection of alignment rules

# Checkpoint for dynamic rule selection
# checkpoint select_alignment_rules:
#     output:
#         touch("checkpoint/rules_selected")
#     run:
#         selected_rules = get_alignment_rules(mode)
#         with open(output[0], "w") as f:
#             f.write("\n".join(selected_rules))

# rule alignment:
#     input:
#         rules_selected="checkpoint/rules_selected"
#     output:
#         bams=expand("results/alignment_{sample}.bam", sample=["crs", "guide"] if config["mode"] == "trans" else ["crs"])
#     run:
#         with open(input.rules_selected) as f:
#             selected_rules = f.read().splitlines()
#         for rule in selected_rules:
#             shell(f"echo test")
#             shell(f"echo {selected_rules}")
#             shell(f"echo {rule}")
#             shell(f"snakemake -f --cores 1 {rule}")
            

# # Rule to execute after checkpoint
# rule alignment:
#     input:
#         rules_selected="checkpoint/rules_selected",
#         guide_bam="results/alignment_guide.bam" if mode == "trans" else None,
#         crs_bam="results/alignment_crs.bam"
#     output:
#         bams=["results/alignment_guide.bam", "results/alignment_crs.bam"] if mode == "trans" else ["results/alignment_crs.bam"]
#     run:
#         with open(input.rules_selected) as f:
#             selected_rules = f.read().splitlines()
#         for rule in selected_rules:
#             shell(f"snakemake -f {rule}")

# import os

# configfile: "config_test.yaml"

# # Rule for the final target
# rule all:
#     input:
#         # "results/final_report.html"
#         "results/alignment_crs.bam",

# #Add here a rule to check, whether reference exists and create .fa if necessary
# rule check_and_convert_reference:
#     input:
#         reference = lambda wildcards: config['references'][wildcards.sample]
#         # reference = config['references'][wildcards.sample]
#     output:
#         reference = "data/{sample}.fa"
#     run:
#         reference = input.reference
#         output = output.reference
#         if reference.endswith(".csv"):
#             shell(f"awk -F, '{{print \">{wildcards.sample}\"$1\"\\n\"$2}}' {reference} > {output}")
#         else:
#             shell(f"cp -f {reference} {output}")

# rule align_guide:
#     input:
#         fastq = config["input_files"]["guide"],
#         reference = "data/reference_guide.fa"
#     output:
#         # bam = "results/alignment_guide.bam"
#         touch("data/reference_guide.fa")
#     params:
#         threads=config["threads_bwa"]
#     shell:
#         """
#         bwa index {input.reference}
#         bwa mem -B 100 -O 100 -E 100 -t {params.threads} {input.reference} {input.fastq} | samtools view -b > results/alignment_guide.bam
#         """

# rule align_crs:
#     input:
#         fastq = config["input_files"]["crs"],
#         reference = "data/reference_crs.fa"   
#     output:
#         # bam = "results/alignment_crs.bam"
#         touch("data/reference_crs.fa")
#     params:
#         threads = config["threads_bwa"]
#     shell:
#         """
#         bwa index {input.reference}
#         bwa mem -t {params.threads} {input.reference} {input.fastq} | samtools view -b > results/alignment_crs.bam
#         """

# rule align_crs_paired:
#     input:
#         reference = "data/reference_crs.fa",
#         FWD = config["input_files"]["crs"],
#         REV = config["input_files"]["crs_rev"]
#     output:
#         # bam = "results/alignment_crs.bam" 
#         touch("data/reference_crs.fa")
#     params:
#         threads = config["threads_bwa"]   
#     shell:
#         "bwa mem -a -t {params.threads} {input.design} {input.FWD} {input.REV} | samtools view -b > {output.bam} "

# # Function to dynamically select the rules for step 1 based on mode
# def get_alignment_rules(mode):
#     if mode == "trans":
#         return ["align_guide", "align_crs"]
#     elif mode == "ss":
#         return ["align_crs"]
#     elif mode == "bulk":
#         return ["align_crs_paired"]
#     else:
#         raise ValueError(f"Unknown mode: {mode}")

# # Get the mode from the config file
# mode = config["mode"]

# # Dynamic selection of alignment rules
# selected_rules = get_alignment_rules(mode)

# # Rule alignment that dynamically includes the selected rules
# rule alignment:
#     input:
#         guide_ref = "data/reference_guide.fa" if "align_guide" in selected_rules else None,
#         crs_ref = "data/reference_crs.fa"
#     output:
#         bam = "results/alignment_crs.bam",     
#     run:
#         for rule in selected_rules:
#             shell(f"snakemake -f {rule}")
            
            
            
            
            
# # # Rule to process files with a Python script
# # rule process_files:
# #     input:
# #         preprocess_guide = "results/alignment_guide.bam" if config["mode"] == "trans" else config["input_files"]["guide"],
# #         preprocess_crs = "results/alignment_crs.bam" 
# #     output:
# #         csv_gz = "results/processed_data.csv.gz"
# #     script:
# #         "scripts/process_files.py"

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

# # # Function to dynamically select the rules for step 1 based on mode
# # def get_alignment_rules(mode):
# #     if mode == "trans":
# #         return ["align_guide", "align_crs"]
# #     elif mode == "ss":
# #         return ["align_crs"]
# #     elif mode == "bulk":
# #         return ["align_crs_paired"]
# #     else:
# #         raise ValueError(f"Unknown mode: {mode}")

# # # Get the mode from the config file
# # mode = config["mode"]

# # # Dynamic selection of step 1 rules
# # rule alignment:
# #     output:
# #         ["results/alignment_guide.bam", "results/alignment_crs.bam"] if mode == "trans" else ["results/alignment_crs.bam"]
# #     run:
# #         selected_rules = get_alignment_rules(mode)
# #         for rule in selected_rules:
# #             shell(f"snakemake {rule}")
