###
#
###

import csv
import itertools

configfile: "./config.yaml"

def sniff_delimiter(file_path, sample_size=2048):
    with open(file_path, 'r') as fp:
        sample = fp.read(sample_size)
        return csv.Sniffer().sniff(sample).delimiter

SAMPLES = []
SYMBOLS = []
max_rows = 5 if config.get("test") else None
sample_delim = sniff_delimiter(config['sample_metadata'])
symbol_delim = sniff_delimiter(config['gene_metadata'])

with open(config['sample_metadata'], mode='r', newline='') as fp:
    reader = csv.DictReader(fp, delimiter=sample_delim)
    for row in itertools.islice(reader, max_rows):
        SAMPLES.append(row['run'])

with open(config['gene_metadata'], mode='r', newline='') as fp:
    reader = csv.DictReader(fp, delimiter=symbol_delim)
    for row in itertools.islice(reader, max_rows):
        SYMBOLS.append(row['symbol'])

OUTDIR = [ config.get('output_directory') if config.get('output_directory') is not None else './workflow-results' ]

# Dictionary for correct slurm batch allocations
HICPUS_JOBS = {1: ['low2', 1], 2: ['low2', 1], 3: ['med2', 25], 4: ['med2', 25], 5: ['high2', 50]}
BIGMEM_JOBS = {1: ['bml', 1], 2: ['bml', 1], 3: ['bmm', 25], 4: ['bmm', 25], 5: ['bmh', 50]}

wildcard_constraints:
    sample = "(SRR|ERR|DRR)\d+",
    symbol = "(?!SRR|ERR|DRR)[A-Z0-9_]+(?:\.[0-9]+)?"


# how does snakemake work if you need to activate a conda environment or request resources?
# change to general paths so this is usable for Colton 
rule all:
  input:
    # how do the wildcards work when you get rid of an input in the rule all because you make a new rule? 
    # look into the wildcards thing
   # expand("{o}/fastp_output/trim/{sample}_1_trim.fastq", sample = ALL_NAMES),
   # expand("{o}/fastp_output/trim/{sample}_2_trim.fastq", sample = ALL_NAMES),
   # expand("{o}/fastp_output/trim/reports/{sample}_trim.json", sample = ALL_NAMES),
   # expand("{o}/fastp_output/trim/reports/{sample}_trim.html", sample = ALL_NAMES),
    expand("{o}/mapped_reads/stats/{sample}.x.{symbol}_mapping_stats.tsv", o=OUTDIR, sample=SAMPLES, symbol=SYMBOLS),
   # expand("{o}/mapped_reads/stats/{sample}.sorted.bam.bai", sample=ALL_NAMES),
   # expand("{o}/blastp_results/{sample}_blast_filtered.csv", sample=ALL_NAMES),
   # "{o}/prodigal_output/combined_prodigal.faa",
   # "{o}/blast_dbs/prodigal_db.pin",
   # "{o}/blast_dbs/prodigal_db.phr",
   # "{o}/blast_dbs/prodigal_db.psq",
  #  expand("{o}/mapped_reads/{sample}.sorted.bam", sample = ALL_NAMES),
  
# Add R1 and R2 dynamically 
# Ask Colton about no phred score qc, detecting adapter seq, minimum length req, see other snkfile

#make rules for getting ibd fastq samples from SRA 

rule download_sra:
    output:
        sra_file = temp("{o}/sra/{sample}.sra"),
        checksum = "{o}/validate/{sample}-checksum.log",
        stats = '{o}/stats/{sample}-counts.xml',
    threads:
        10
    resources:
        mem_mb = lambda wildcards, attempt: 8 * 1024 * attempt,
        disk_mb = lambda wildcards, attempt: 8 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs= lambda wildcards, attempt: HICPUS_JOBS[attempt][1],
        partition= lambda wildcards, attempt: HICPUS_JOBS[attempt][0],
    conda:
        "envs/sra.yaml"
    shell:
       """
           (
               aws s3 cp --no-sign-request \
                   s3://sra-pub-run-odp/sra/{wildcards.sample}/{wildcards.sample} \
                   {output.sra_file} \
               || prefetch {wildcards.sample} -o {output.sra_file}
           )
           # Validate the sra file against the checksums provided
           vdb-validate {output.sra_file} -I no &> {output.checksum}

           if grep -q 'err' {output.checksum};
             then
               echo 'Validation of {wildcards.sample} failed. Removing {output.sra_file}...'
               rm {output.sra_file}
               echo 'Validation failed, exiting.'
               exit 1
           else
               echo 'Sample {wildcards.sample} validated without error'
           fi

           # Stats for the sequence contained in the sra file
           sra-stat --alignment off --quick -x {output.sra_file} &> {output.stats}
       """

rule dump_fastq:
    input:
        xml_script = 'scripts/process_xml_file.py',
        seq_script = 'scripts/process_fastq_files.sh',
        sra_file = "{o}/sra/{sample}.sra",
        stats = '{o}/stats/{sample}-counts.xml',
    output:
        r1 = temp("{o}/fastq/{sample}_1.fastq.gz"),
        r2 = temp("{o}/fastq/{sample}_2.fastq.gz"),
    conda:
        "envs/sra.yaml",
    resources:
        mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
        disk_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs= lambda wildcards, attempt: HICPUS_JOBS[attempt][1],
        partition= lambda wildcards, attempt: HICPUS_JOBS[attempt][0],
    threads:
        16
    params:
        sra_counts = 'all-sra-counts.csv',
        fastq_counts = 'all-fastq-counts.csv',
    shell:
        """
        # job ID returning 0 in log files when running on cluster (use SLURM_JOB_ID if available, otherwise use Snakemake jobid)
        job=${{SLURM_JOB_ID:-{jobid}}}
        # make a directory specific to user and job
        export MYTMP="/scratch/${{USER}}/slurm_{rule}.${{job}}"
        mkdir -v -p "$MYTMP"/scripts

        # force clean it up after job script ends
        function cleanup() {{ rm -rf "$MYTMP"; }}
        trap cleanup EXIT

        # change to it!
        cp {input.xml_script} "$MYTMP"/scripts
        cp {input.seq_script} "$MYTMP"/scripts
        cd "$MYTMP"

        # run your code, telling it to use that dir:
        echo "Changing to $(pwd)"

        echo "Running fasterq-dump for {wildcards.sample}"
        echo ''

        # Get fastqc file
        # To verify any errors due to threads -> https://github.com/sourmash-bio/sourmash/issues/2941
        # Quality check, if failed check use Worts conservative approach
        # https://github.com/sourmash-bio/wort/blob/09d6ee552c85360fdc60e9d59e2fc757c2c6dd99/wort/blueprints/compute/tasks.py#L38-L61
        fasterq-dump {input.sra_file} --skip-technical --split-files --progress \
                     --threads {threads} --bufsize 1000MB --curcache 10000MB --mem {resources.mem_mb}

        echo ''
        echo "Fastq files generated by fasterq-dump:"
        ls -h "$MYTMP"
        echo ''

        # For checking the status of fastq files because the md5checksums are difficult to use with this file type
        {input.xml_script} {input.stats} -o {params.sra_counts}
        {input.xml_script} {input.stats} -o {wildcards.o}/stats/{params.sra_counts}

        # capture the awk variables in the 3rd and 4th column of the counts file for each sample
        read sra_read_count sra_base_count < <(awk -F',' -v sample="{wildcards.sample}" '
            $0 ~ sample {{
                sra_read_count = $3
                sra_base_count = $4
                print sra_read_count, sra_base_count
                exit
            }}
        ' {params.sra_counts})

        sra_base_count=$(echo "$sra_base_count" | tr -d '\r')

        R1=$MYTMP/{wildcards.sample}_1.fastq
        R2=$MYTMP/{wildcards.sample}_2.fastq
        R1_T=$MYTMP/{wildcards.sample}_1_trim.fastq
        R2_T=$MYTMP/{wildcards.sample}_2_trim.fastq

        {input.seq_script} {threads} ./{params.fastq_counts} "$R1" "$R2"
        # Only add the output if it is not in the output...
        if ! grep -q "{wildcards.sample}_1" ./{params.fastq_counts} && ! grep -q "{wildcards.sample}_2" ./{params.fastq_counts}; then
            cat ./{params.fastq_counts} >> {wildcards.o}/stats/{params.fastq_counts}
        fi

        # Create the count variables from the file
        read fq1_read_count fq1_base_count < <(awk -F',' -v sample="{wildcards.sample}_1" '
            $0 ~ sample {{
                fq1_read_count = $2
                fq1_base_count = $3
                print fq1_read_count, fq1_base_count
                exit
            }}
        ' {params.fastq_counts})
        read fq2_read_count fq2_base_count < <(awk -F',' -v sample="{wildcards.sample}_2" '
            $0 ~ sample {{
                fq2_read_count = $2
                fq2_base_count = $3
                print fq2_read_count, fq2_base_count
                exit
            }}
        ' {params.fastq_counts})

        # if there is a comma, remove the comma
        fq1_read_count=$(echo "$fq1_read_count" | tr -d ',')
        fq2_read_count=$(echo "$fq2_read_count" | tr -d ',')
        fq1_base_count=$(echo "$fq1_base_count" | tr -d ',')
        fq2_base_count=$(echo "$fq2_base_count" | tr -d ',')

        echo "$fq1_read_count reads in {wildcards.sample}_1"
        echo "$fq2_read_count reads in {wildcards.sample}_2"
        echo "$sra_read_count reads in sra file"
        echo "$(( fq1_base_count + fq2_base_count )) bases in {wildcards.sample}_1 and {wildcards.sample}_2"
        echo "$sra_base_count in sra file"

        # Compare the SRA file counts and the FASTQ file counts
        if [[ $fq1_read_count == $sra_read_count && $fq2_read_count == $sra_read_count && $(( fq1_base_count + fq2_base_count )) == $sra_base_count ]]; then
            echo "Read and base counts match!!!"

       # Remove and use the unthreaded, tried-and-true approach
        else
            echo "Read or base counts do not match!!!"
            rm "$R1" "$R2"
            fastq-dump --disable-multithreading --skip-technical --split-files --clip {input.sra_file} -O "$MYTMP"
        fi

        echo ''
        echo "Processing R1..."
        seqtk seq -C "$R1" | \
            perl -ne 's/\.([12])$/\/$1/; print $_' | \
            gzip -9 -c > {output.r1} &

        echo "Processing R2..."
        seqtk seq -C "$R2" | \
            perl -ne 's/\.([12])$/\/$1/; print $_' | \
            gzip -9 -c > {output.r2} &

        wait
        echo "Finished downloading raw reads"
        """


rule fastp:
  input:
    fq1 = "{o}/fastq/{sample}_1.fastq.gz",
    fq2 = "{o}/fastq/{sample}_2.fastq.gz"
  output: 
    trim1 = "{o}/fastp_output/trim/{sample}_1_trim.fastq",
    trim2 = "{o}/fastp_output/trim/{sample}_2_trim.fastq",
    json = "{o}/fastp_output/trim/reports/{sample}_trim.json",  
    html = "{o}/fastp_output/trim/reports/{sample}_trim.html"
  resources:
    mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
    disk_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
    time = lambda wildcards, attempt: 1.5 * 60 * attempt,
    runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
    allowed_jobs= lambda wildcards, attempt: BIGMEM_JOBS[attempt][1],
    partition= lambda wildcards, attempt: BIGMEM_JOBS[attempt][0],
  conda: "envs/sebastian_env.yml"
  shell:
    """
    fastp \
      --in1 {input.fq1} \
      -I {input.fq2} \
      --out1 {output.trim1} \
      -O {output.trim2} \
      --json {output.json} \
      --html {output.html}
    """

rule fastqc:
  input:
    trim1 = "{o}/fastp_output/trim/{sample}_1_trim.fastq",
    trim2 = "{o}/fastp_output/trim/{sample}_2_trim.fastq",
  output:
    html1 = "{o}/fastqc_reports/{sample}_1_fastqc.html",
    html2 = "{o}/fastqc_reports/{sample}_2_fastqc.html",
    zip1 = "{o}/fastqc_reports/{sample}_1_fastqc.zip",
    zip2 = "{o}/fastqc_reports/{sample}_2_fastqc.zip",
  resources:
    mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
    disk_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
    time = lambda wildcards, attempt: 1.5 * 60 * attempt,
    runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
    allowed_jobs= lambda wildcards, attempt: BIGMEM_JOBS[attempt][1],
    partition= lambda wildcards, attempt: BIGMEM_JOBS[attempt][0],
  conda: "envs/sebastian_env.yml"
  shell:
    """
    fastqc \
      {input.trim1} \
      {input.trim2} \
      --outdir {FASTQC_DIR}
    """

rule megahit:
  input:
    trim1 = "{o}/fastp_output/trim/{sample}_1_trim.fastq",
    trim2 = "{o}/fastp_output/trim/{sample}_2_trim.fastq",
  output:
    contigs = "{o}/megahit_output/{sample}/final.contigs.fa",
  conda: "envs/sebastian_env.yml"
  threads: 16
  resources:
    mem_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
    disk_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
    time = lambda wildcards, attempt: 8 * 60 * attempt,
    runtime = lambda wildcards, attempt: 8 * 60 * attempt,
    allowed_jobs= lambda wildcards, attempt: HICPUS_JOBS[attempt][1],
    partition= lambda wildcards, attempt: HICPUS_JOBS[attempt][0],
  shell:
    """
    # Define expected and temporary output directories
    expected_outdir="{wildcards.o}/megahit_output/{wildcards.sample}"
    timestamp=$(date +"%Y%m%d_%H%M%S")
    tmp_outdir="${{expected_outdir}}_${{timestamp}}"

    # Run MEGAHIT with a unique directory
    megahit \
      -1 {input.trim1} \
      -2 {input.trim2} \
      -o "$tmp_outdir" \
      -t {threads} \
      --continue

    # Move contents to expected output directory
    mkdir -p "$expected_outdir"
    mv "$tmp_outdir"/* "$expected_outdir"/

    # Clean up the temp directory
    rmdir "$tmp_outdir"
    """    

rule prodigal:
  input:
    contigs = "{o}/megahit_output/{sample}/final.contigs.fa",
  output:
    fna = "{o}/prodigal_output/{sample}.prod.fna",
    faa = "{o}/prodigal_output/{sample}.prod.faa",
    gbk = "{o}/prodigal_output/{sample}.prod.gbk",
  conda: "envs/sebastian_env.yml"
  resources:
    mem_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
    disk_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
    time = lambda wildcards, attempt: 8 * 60 * attempt,
    runtime = lambda wildcards, attempt: 8 * 60 * attempt,
    allowed_jobs= lambda wildcards, attempt: HICPUS_JOBS[attempt][1],
    partition= lambda wildcards, attempt: HICPUS_JOBS[attempt][0],
  shell:
    """
    prodigal \
      -p meta \
      -m \
      -i {input.contigs} \
      -d {output.fna} \
      -a {output.faa} \
      -o {output.gbk}
    """

rule combine_prodigal: #Review the shell script because I need to change the header in the files
  input:
    faa = expand("{o}/prodigal_output/{sample}.prod.faa", o=OUTDIR, sample = SAMPLES),
  output:
    combined_faa = "{o}/prodigal_output/combined_prodigal.faa"
  conda: "envs/sebastian_env.yml"
  shell:
    """
    for file in {input.faa}; do
      awk -v filename="$(basename $file)" '/^>/{{print ">" filename "_" substr($0, 2)}} !/^>/' "$file" >> {output.combined_faa}
    done
    """

rule blastp_database:
  input:
    combined_faa = "{o}/prodigal_output/combined_prodigal.faa"
  output:
    pin = "{o}/blast_dbs/prodigal_db.pin",
    phr = "{o}/blast_dbs/prodigal_db.phr",
    psq = "{o}/blast_dbs/prodigal_db.psq"
  params:
    out_prefix = "{o}/blast_dbs/prodigal_db"
  resources:
    mem_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
    disk_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
    time = lambda wildcards, attempt: 8 * 60 * attempt,
    runtime = lambda wildcards, attempt: 8 * 60 * attempt,
    allowed_jobs= lambda wildcards, attempt: BIGMEM_JOBS[attempt][1],
    partition= lambda wildcards, attempt: BIGMEM_JOBS[attempt][0],
  conda: "envs/sebastian_env.yml"
  shell:
    """
    makeblastdb -in {input.combined_faa} -dbtype prot -out {params.out_prefix}
    """

rule get_gene_seq:
    input:
        csv = 'metadata/genes-to-get.csv',
        seq_script = 'scripts/get-gene-seq.py',
    output:
        pgenes = '{o}/gene-seqs/{symbol}_data/protein.faa',
    conda: "envs/requests.yaml",
    shell:"""
        awk -F',' -v sym="{wildcards.symbol}" 'NR==1 {{next}} $1 == sym {{
            printf "Running for symbol=%s, organism=%s, db=%s\\n", $1, $3, $5;
            printf "python {input.seq_script} --genes %s --organism %s --extract -db %s --all-annotation --ext-dir {wildcards.o}/gene-seqs/", $1, $3, $5;
            cmd = sprintf("python {input.seq_script} --genes %s --organism %s --extract -db %s --all-annotation --ext-dir {wildcards.o}/gene-seqs/", $1, $3, $5);
            system(cmd);
        }}' {input.csv}

        grep "{wildcards.symbol}" {output.pgenes}
    """

rule blastp: # should we make a bastn rule and a blastn_db?
  input:
    genes = "{o}/gene-seqs/{symbol}_data/protein.faa",
    db = "{o}/blast_dbs/prodigal_db.pin"
  output:
     temp = temporary("{o}/blastp_results/temp_{symbol}_blast.csv"),
     csv = "{o}/blastp_results/{symbol}_blast_output.csv",
  params:
     db = "{o}/blast_dbs/prodigal_db"
  resources:
    mem_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
    disk_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
    time = lambda wildcards, attempt: 8 * 60 * attempt,
    runtime = lambda wildcards, attempt: 8 * 60 * attempt,
    allowed_jobs= lambda wildcards, attempt: BIGMEM_JOBS[attempt][1],
    partition= lambda wildcards, attempt: BIGMEM_JOBS[attempt][0],
  conda: "envs/sebastian_env.yml"
  shell:
    """
    blastp \
      -query {input.genes} \
      -db {params.db} \
      -outfmt "10 qseqid sseqid pident length evalue bitscore" \
      -out {output.temp}

    cat <(printf "qseqid,sseqid,pident,length,e_value,bit_score\\n") \
        {output.temp} > {output.csv}
    """
 
rule parse_blast_output:
  input:
    csv = "{o}/blastp_results/{symbol}_blast_output.csv",
    script = "scripts/parse_blast_output.py",
  output:
    filtered = "{o}/blastp_results/{symbol}_blast_filtered.csv"
  params:
    cutoff = "1e-07"
  resources:
    mem_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
    disk_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
    time = lambda wildcards, attempt: 8 * 60 * attempt,
    runtime = lambda wildcards, attempt: 8 * 60 * attempt,
    allowed_jobs= lambda wildcards, attempt: BIGMEM_JOBS[attempt][1],
    partition= lambda wildcards, attempt: BIGMEM_JOBS[attempt][0],
  conda: "envs/sebastian_env.yml"
  shell:
    """
    {input.script} {input.csv} \
      -o {output.filtered} \
      --cutoff {params.cutoff} \
      --verbose
    """
    # look at this py script to see what output files. I need a contigs list fna file for downstream rules
rule make_contig_list: 
  input:
    blastp_hits = "{o}/blastp_results/{symbol}_blast_filtered.csv",
    script = "scripts/extract_prodigal_contigs.py",
  output:
    contig_list_fna = "{o}/contiglists/{sample}.x.{symbol}.contiglist.fna"
  params:
    contig_list_csv = "{o}/contiglist.csv",
    fna_dir = "{o}/prodigal_output/", 
    fna_ext = "fna"
  resources:
    mem_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
    disk_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
    time = lambda wildcards, attempt: 8 * 60 * attempt,
    runtime = lambda wildcards, attempt: 8 * 60 * attempt,
    allowed_jobs= lambda wildcards, attempt: BIGMEM_JOBS[attempt][1],
    partition= lambda wildcards, attempt: BIGMEM_JOBS[attempt][0],
  conda: "envs/sebastian_env.yml"
  shell:
    """
    {input.script} {input.blastp_hits} \
      -v \
      -o {params.contig_list_csv} \
      -d {params.fna_dir} \
      -f {params.fna_ext} \
      --outdir $(dirname {output.contig_list_fna})
    """

rule minimap2:
  input:
    contig_list_fna = "{o}/contiglists/{sample}.x.{symbol}.contiglist.fna",
    fq1 = "{o}/fastq/{sample}_1.fastq.gz",
    fq2 = "{o}/fastq/{sample}_2.fastq.gz"
  output:
    sorted_bam = "{o}/mapped_reads/{sample}.x.{symbol}.sorted.bam"
  threads: 16
  resources:
    mem_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
    disk_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
    time = lambda wildcards, attempt: 8 * 60 * attempt,
    runtime = lambda wildcards, attempt: 8 * 60 * attempt,
    allowed_jobs= lambda wildcards, attempt: HICPUS_JOBS[attempt][1],
    partition= lambda wildcards, attempt: HICPUS_JOBS[attempt][0],
  conda: "envs/sebastian_env.yml"
  shell:
    """
    minimap2 -ax sr -t {threads} \
      {input.contig_list_fna} \
      {input.fq1} {input.fq2} \
    | samtools sort -@ {threads} -o {output.sorted_bam}
    """
   
rule index_bam:
  input:
    sorted_bam = "{o}/mapped_reads/{sample}.x.{symbol}.sorted.bam"
  output:
    sorted_bai = "{o}/mapped_reads/{sample}.x.{symbol}.sorted.bam.bai"
  conda: "envs/sebastian_env.yml"
  resources:
    mem_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
    disk_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
    time = lambda wildcards, attempt: 8 * 60 * attempt,
    runtime = lambda wildcards, attempt: 8 * 60 * attempt,
    allowed_jobs= lambda wildcards, attempt: BIGMEM_JOBS[attempt][1],
    partition= lambda wildcards, attempt: BIGMEM_JOBS[attempt][0],
  shell:
    """
    samtools index {input.sorted_bam}
    """

rule reads_mapped:
  input:
    sorted_bam = "{o}/mapped_reads/{sample}.x.{symbol}.sorted.bam.bai",
    bam = "{o}/mapped_reads/{sample}.x.{symbol}.sorted.bam",
    script = "scripts/aggregating_depth.py",
  output:
    map = "{o}/mapped_reads/stats/{sample}.x.{symbol}_mapping_stats.tsv",
    coverage = "{o}/mapped_reads/stats/{sample}.x.{symbol}_coverage_stats.tsv",
    ave_cov = "{o}/mapped_reads/stats/{sample}.x.{symbol}_average_coverage_stats.tsv",
  resources:
    mem_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
    disk_mb = lambda wildcards, attempt: 24 * 1024 * attempt,
    time = lambda wildcards, attempt: 8 * 60 * attempt,
    runtime = lambda wildcards, attempt: 8 * 60 * attempt,
    allowed_jobs= lambda wildcards, attempt: BIGMEM_JOBS[attempt][1],
    partition= lambda wildcards, attempt: BIGMEM_JOBS[attempt][0],
  conda: "envs/sebastian_env.yml"
  shell:
    """
    (echo -e "contig\tlength\tmapped\tunmapped\tnet_mapped"; \
    samtools idxstats {input.bam} | \
    awk -F'\t' '{{print $0 "\t" ($3 - $4)}}') > {output.map}

    samtools depth -aa {input.bam} > {output.coverage}

    {input.script} {output.coverage} {output.ave_cov}
    """
