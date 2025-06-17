import os
import re
from snakemake.io import expand

#Change

# The config file is loaded that specifies the input and output directory
configfile: "config/config.yaml"

# All files in the input directory are listed
input_dir = config["input_dir"]
all_files = os.listdir(input_dir)

# We pull the fastq.gz files specifically
input_files = [
    f for f in all_files 
    if re.match(r".*_R[12]_.*\.fastq\.gz$", f)
]

# Unique subject IDs are extracted
input_subjects = {
    re.match(r"([A-Za-z0-9]+)_.*_R[12]_.*", f).group(1)
    for f in input_files if re.match(r"([A-Za-z0-9]+)_.*_R[12]_.*", f)
}

#METADATA = pd.read_csv(config["METADATA"].strip())
#SAMPLES = METADATA["Sample"].tolist()
#RAW_FWD_READS = METADATA[config["fwd_reads_path"]]
#RAW_REV_READS = METADATA[config["rev_reads_path"]]

#print("Input samples:", input_files)




# Debug: The input files and input subjects are printed for verification
print("Input files:", input_files)
print("Input subjects:", input_subjects)

# Full filenames for R1 and R2 based on subject are pulled for future outputs
input_files_dict = {
    subject: [os.path.join(input_dir, f) for f in input_files if re.match(rf"{subject}_.*_R[12]_.*\.fastq\.gz$", f)]
    for subject in input_subjects
}

#The rull all specifies all the final outputs of the snakemake workflow
rule all:
    input:
        expand(os.path.join(config["output_dir"], "unzipped/{subject}_R1.fastq"), subject=input_subjects) + \
        expand(os.path.join(config["output_dir"], "unzipped/{subject}_R2.fastq"), subject=input_subjects) + \
        expand(os.path.join(config["output_dir"], "trimmed/{subject}_R1_trimmed.fastq"), subject=input_subjects) + \
        expand(os.path.join(config["output_dir"], "trimmed/{subject}_R2_trimmed.fastq"), subject=input_subjects) + \
        expand(os.path.join(config["output_dir"], "summary/{subject}_trimmomatic_summary.txt"), subject=input_subjects) + \
        expand(os.path.join(config["output_dir"], "trimlog/{subject}_trimmomatic_log.txt"), subject=input_subjects) + \
        expand(os.path.join(config["output_dir"], "truncated/{subject}_trimmed_truncated.R1.fastq"), subject=input_subjects) + \
        expand(os.path.join(config["output_dir"], "truncated/{subject}_trimmed_truncated.R2.fastq"), subject=input_subjects) + \
        expand(os.path.join(config["output_dir"], "repaired/{subject}_repaired.R1.fastq"), subject=input_subjects) + \
        expand(os.path.join(config["output_dir"], "repaired/{subject}_repaired.R2.fastq"), subject=input_subjects) + \
        expand(os.path.join(config["output_dir"], "singletons/{subject}_singletons.fastq"), subject=input_subjects) + \
        pj(config["output_dir"], "FASTQC_reports") + \
        pj(config["output_dir"], "FASTQC_reports/multiqc_report")

#The fastq.gz files are unzipped
rule unzip:
    input:
        gz_file=lambda wildcards: [f for f in input_files_dict[wildcards.subject] if f"_R{wildcards.read}_" in f]
    output:
        unzipped=os.path.join(config["output_dir"], "unzipped/{subject}_R{read}.fastq")
    shell:
        """
        gunzip -c {input.gz_file[0]} > {output.unzipped}
        """

#The unzipped files are trimmed remvoing the adapters
rule trim:
    input:
        r1=os.path.join(config["output_dir"], "unzipped/{subject}_R1.fastq"),
        r2=os.path.join(config["output_dir"], "unzipped/{subject}_R2.fastq")
    output:
        r1_trimmed=os.path.join(config["output_dir"], "trimmed/{subject}_R1_trimmed.fastq"),
        r2_trimmed=os.path.join(config["output_dir"], "trimmed/{subject}_R2_trimmed.fastq"),
        summary=os.path.join(config["output_dir"], "summary/{subject}_trimmomatic_summary.txt"),
        trimlog=os.path.join(config["output_dir"], "trimlog/{subject}_trimmomatic_log.txt")
    conda:
        "envs/trimmomatic.yaml"
    shell:
        """
        trimmomatic PE {input.r1} {input.r2} \
        {output.r1_trimmed} {output.r1_trimmed} \
        {output.r2_trimmed} {output.r2_trimmed} \
        -trimlog {output.trimlog} \
        -summary {output.summary} \
        ILLUMINACLIP:src/Trimmomatic-0.39/adapters/adapters_all.fa:2:30:10 LEADING:20 TRAILING:20 MINLEN:50 SLIDINGWINDOW:4:20 \
        -phred33
        > {output.summary}
        """

#Trimmed reads are not truncated 
rule truncate_reads:
    input:
        r1_read_trim=os.path.join(config["output_dir"], "trimmed/{subject}_R1_trimmed.fastq"),
        r2_read_trim=os.path.join(config["output_dir"], "trimmed/{subject}_R2_trimmed.fastq")
    output:
        r1_truncated=os.path.join(config["output_dir"], "truncated/{subject}_trimmed_truncated.R1.fastq"),
        r2_truncated=os.path.join(config["output_dir"], "truncated/{subject}_trimmed_truncated.R2.fastq")
    conda:
        "envs/seqtk.yaml"  # Path to your seqtk conda environment
    shell:
        """
        # Truncate forward reads (5 bases from the beginning)
        seqtk trimfq -b 5 -e 0 {input.r1_read_trim} > {output.r1_truncated}
        echo "Forward read truncated for {wildcards.subject}"

        # Truncate reverse reads (10 bases from the beginning)
        seqtk trimfq -b 10 -e 0 {input.r2_read_trim} > {output.r2_truncated}
        echo "Reverse read truncated for {wildcards.subject}"
        """
    
#We repair the reads if they were differentially truncated or trimmed
rule repair:
    input:
        r1_read_truncated=os.path.join(config["output_dir"], "truncated/{subject}_trimmed_truncated.R1.fastq"),
        r2_read_truncated=os.path.join(config["output_dir"], "truncated/{subject}_trimmed_truncated.R2.fastq")
    output:
        r1_repaired=os.path.join(config["output_dir"], "repaired/{subject}_repaired.R1.fastq"),
        r2_repaired=os.path.join(config["output_dir"], "repaired/{subject}_repaired.R2.fastq"),
        singleton=os.path.join(config["output_dir"], "singletons/{subject}_singletons.fastq"),
    conda:
        "envs/bbmap.yaml"  # A conda environment containing Bowtie2
    shell:
        """
        repair.sh in={input.r1_read_truncated} in2={input.r2_read_truncated} out={output.r1_repaired} out2={output.r2_repaired} outs={output.singleton}
        """

#FastQC is used for quality control 
rule fastQC:
    input: 
        r1_repaired=os.path.join(config["output_dir"], "repaired/{subject}_repaired.R1.fastq"),
        r2_repaired=os.path.join(config["output_dir"], "repaired/{subject}_repaired.R2.fastq")
    output:
        fastqc_dir=directory(os.path.join(config["output_dir"], "FastQC_reports"))
    conda:
        "envs/fastqc.yaml"
    shell:
        """
        mkdir -p {output.fastqc_dir} 
        fastqc {input.r1_repaired} {input.r2_repaired} -o {output.fastqc_dir} 
        """ 

#MultiQC    
rule multiQC:
    input: 
        fastqc_dir=directory(os.path.join(config["output_dir"], "FastQC_reports"))
    output:
        multiqc_report=os.path.join(config["output_dir"], "FASTQC_reports/multiqc_report")
    conda:
        "envs/fastqc.yaml"
    shell:
        """
        multiqc {input.fastqc_dir} {output.multiqc_report} 
        """    

