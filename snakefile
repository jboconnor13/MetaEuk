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
        expand(os.path.join(config["output_dir"], "singletons/{subject}_singletons.fastq"), subject=input_subjects) + \
        expand(os.path.join(config["output_dir"], "FastQC_reports/{subject}"), subject=input_subjects) + \
        [os.path.join(config["output_dir"], "FASTQC_reports", "multiqc_report.html")] + \
        expand(os.path.join(config["output_dir"], "nonhost/{subject}.hostile_report.html"), subject=input_subjects)

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
        fastqc_output_dir=directory(os.path.join(config["output_dir"], "FastQC_reports/{subject}")),
        html_report_R1=os.path.join(config["output_dir"], "FastQC_reports/{subject}", "{subject}_repaired.R1_fastqc.html"), # FastQC names this based on input file name
        html_report_R2=os.path.join(config["output_dir"], "FastQC_reports/{subject}", "{subject}_repaired.R2_fastqc.html") # FastQC names this based on input file name
    conda:
        "envs/fastqc.yaml"
    shell:
        """
        # MultiQC will later look for these individual FastQC outputs
        # Ensure the output directory exists
        mkdir -p {output.fastqc_output_dir}
        # Run fastqc, outputting to the subject-specific directory
        fastqc {input.r1_repaired} {input.r2_repaired} -o {output.fastqc_output_dir}
        """

#MultiqC compiles results
rule multiQC:
    input:
        fastqc_output_dir=directory(os.path.join(config["output_dir"], "FastQC_reports")) 
    output:
        multiqc_report_name=[os.path.join(config["output_dir"], "FASTQC_reports", "multiqc_report.html")]
    params:
        fastqc_search_dir=directory(os.path.join(config["output_dir"], "FastQC_reports"))
    conda:
        "envs/fastqc.yaml"
    shell:
        """
        multiqc {params.fastqc_search_dir} --outdir {input.fastqc_output_dir} --filename {output.multiqc_report_name} --dirs-depth 2
        """


#HOstile is used for host filtering
rule hostile:
    input:
        r1_repaired=os.path.join(config["output_dir"], "repaired/{subject}_repaired.R1.fastq"),
        r2_repaired=os.path.join(config["output_dir"], "repaired/{subject}_repaired.R2.fastq"),
        hostile_dir=directory(config["hostile_dir"])
    output:
        r1_nonhost=os.path.join(config["output_dir"], "nonhost/{subject}.R1.fq.gz"),
        r2_nonhost=os.path.join(config["output_dir"], "nonhost/{subject}.R2.fq.gz"),
        hostile_report=os.path.join(config["output_dir"], "nonhost/{subject}.hostile_report.html")
    params:
        hostile_out_dir=directory(os.path.join(config["output_dir"], "nonhost"))
    conda:
        "envs/hostile.yaml"
    shell:
        """
        hostile clean \
            --fastq1 {input.r1_repaired} --fastq2 {input.r2_repaired} \
            -o {params.hostile_out_dir} \
            --threads 1 \
            --index {input.hostile_dir} \
            --debug \
            --aligner bowtie2 \
            > {output.hostile_report}
        """




 hostile clean \
 --fastq1 /Users/johnoconnor/MetaEuk_Remote/output/repaired/ZIM089_repaired.R1.fastq \
 --fastq2 /Users/johnoconnor/MetaEuk_Remote/output/repaired/ZIM089_repaired.R2.fastq \
 -o /Users/johnoconnor/MetaEuk_Remote/output/nonhost \
 --threads 1 \
 --index /Users/johnoconnor/MetaEuk_Remote/ref_databases/hostile \
 --debug \
 --aligner bowtie2 \
 > /Users/johnoconnor/MetaEuk_Remote/output/nonhost/ZIM089.hostile_report.html
 
 taEuk %  hostile clean \
 --fastq1 /Users/johnoconnor/MetaEuk_Remote/output/repaired/ZIM089_repaired.R1.fastq \
 --fastq2 /Users/johnoconnor/MetaEuk_Remote/output/repaired/ZIM089_repaired.R2.fastq \
 -o /Users/johnoconnor/MetaEuk_Remote/output/nonhost \
 --threads 1 \
 --index /Users/johnoconnor/MetaEuk_Remote/ref_databases/hostile \
 --debug \
 --aligner bowtie2 \
 > /Users/johnoconnor/MetaEuk_Remote/output/nonhost/ZIM089.hostile_report.html
15:05:28 DEBUG: clean_paired_fastqs() threads=1 aligner_threads=1 compression_threads=0 util.CACHE_DIR=PosixPath('/Users/johnoconnor/Library/Application Support/hostile') util.INDEX_REPOSITORY_URL='https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o'
15:05:28 INFO: Hostile v2.0.1. Mode: paired short read (Bowtie2)
15:05:30 DEBUG: Fetching bucket contents
15:05:30 DEBUG: connect_tcp.started host='objectstorage.uk-london-1.oraclecloud.com' port=443 local_address=None timeout=5.0 socket_options=None
15:05:30 DEBUG: connect_tcp.failed exception=ConnectError(gaierror(8, 'nodename nor servname provided, or not known'))
Traceback (most recent call last):
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/httpx/_transports/default.py", line 101, in map_httpcore_exceptions
    yield
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/httpx/_transports/default.py", line 250, in handle_request
    resp = self._pool.handle_request(req)
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/httpcore/_sync/connection_pool.py", line 256, in handle_request
    raise exc from None
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/httpcore/_sync/connection_pool.py", line 236, in handle_request
    response = connection.handle_request(
        pool_request.request
    )
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/httpcore/_sync/connection.py", line 101, in handle_request
    raise exc
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/httpcore/_sync/connection.py", line 78, in handle_request
    stream = self._connect(request)
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/httpcore/_sync/connection.py", line 124, in _connect
    stream = self._network_backend.connect_tcp(**kwargs)
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/httpcore/_backends/sync.py", line 207, in connect_tcp
    with map_exceptions(exc_map):
         ~~~~~~~~~~~~~~^^^^^^^^^
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/contextlib.py", line 162, in __exit__
    self.gen.throw(value)
    ~~~~~~~~~~~~~~^^^^^^^
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/httpcore/_exceptions.py", line 14, in map_exceptions
    raise to_exc(exc) from exc
httpcore.ConnectError: [Errno 8] nodename nor servname provided, or not known

The above exception was the direct cause of the following exception:

Traceback (most recent call last):
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/hostile/util.py", line 157, in fetch_manifest
    r = httpx.get(f"{url}/manifest.json")
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/httpx/_api.py", line 195, in get
    return request(
        "GET",
    ...<9 lines>...
        trust_env=trust_env,
    )
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/httpx/_api.py", line 109, in request
    return client.request(
           ~~~~~~~~~~~~~~^
        method=method,
        ^^^^^^^^^^^^^^
    ...<8 lines>...
        follow_redirects=follow_redirects,
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    )
    ^
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/httpx/_client.py", line 825, in request
    return self.send(request, auth=auth, follow_redirects=follow_redirects)
           ~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/httpx/_client.py", line 914, in send
    response = self._send_handling_auth(
        request,
    ...<2 lines>...
        history=[],
    )
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/httpx/_client.py", line 942, in _send_handling_auth
    response = self._send_handling_redirects(
        request,
        follow_redirects=follow_redirects,
        history=history,
    )
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/httpx/_client.py", line 979, in _send_handling_redirects
    response = self._send_single_request(request)
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/httpx/_client.py", line 1014, in _send_single_request
    response = transport.handle_request(request)
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/httpx/_transports/default.py", line 249, in handle_request
    with map_httpcore_exceptions():
         ~~~~~~~~~~~~~~~~~~~~~~~^^
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/contextlib.py", line 162, in __exit__
    self.gen.throw(value)
    ~~~~~~~~~~~~~~^^^^^^^
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/httpx/_transports/default.py", line 118, in map_httpcore_exceptions
    raise mapped_exc(message) from exc
httpx.ConnectError: [Errno 8] nodename nor servname provided, or not known

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/bin/hostile", line 10, in <module>
    sys.exit(main())
             ~~~~^^
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/hostile/cli.py", line 170, in main
    defopt.run(
    ~~~~~~~~~~^
        {
        ^
    ...<9 lines>...
        cli_options="has_default",
        ^^^^^^^^^^^^^^^^^^^^^^^^^^
    )
    ^
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/defopt.py", line 308, in run
    return bind(
           ~~~~~
    ...<2 lines>...
        no_negated_flags=no_negated_flags, version=version,
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        argparse_kwargs=argparse_kwargs, intermixed=intermixed, argv=argv)()
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~^^
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/hostile/cli.py", line 71, in clean
    stats = lib.clean_paired_fastqs(
        [(fastq1, fastq2)],
    ...<10 lines>...
        airplane=airplane,
    )
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/hostile/lib.py", line 278, in clean_paired_fastqs
    index_path = aligner.value.check_index(index, airplane=airplane)
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/hostile/aligner.py", line 56, in check_index
    elif not airplane and util.fetch_manifest(util.INDEX_REPOSITORY_URL).get(
                          ~~~~~~~~~~~~~~~~~~~^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/Users/johnoconnor/MetaEuk/.snakemake/conda/1370aea1fc2311860281e28af47228cb_/lib/python3.13/site-packages/hostile/util.py", line 160, in fetch_manifest
    raise httpx.HTTPError(
    ...<3 lines>...
    )
httpx.HTTPError: Failed to fetch manifest.json from object storage. Ensure you are connected to the internet, or provide a valid path to a local index