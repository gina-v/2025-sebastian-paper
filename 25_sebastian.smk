ALL_NAMES = ["ERR13295987", "ERR13295988", "ERR13295990", "ERR13295993", "ERR13295994"]

FASTQC_DIR = "/home/gmvaz/2025-sebastian-paper/fastqc_reports"
MEGAHIT_DIR = "/home/gmvaz/2025-sebastian-paper/megahit_output"

# how does snakemake work if you need to activate a conda environment or request resources?
# change to general paths so this is usable for Colton 
rule all:
  input:
    # how do the wildcards work when you get rid of an input in the rule all because you make a new rule? 
    # look into the wildcards thing
   # expand("/home/gmvaz/2025-sebastian-paper/fastp_output/trim/{sample}_1_trim.fastq", sample = ALL_NAMES),
   # expand("/home/gmvaz/2025-sebastian-paper/fastp_output/trim/{sample}_2_trim.fastq", sample = ALL_NAMES),
   # expand("/home/gmvaz/2025-sebastian-paper/fastp_output/trim/reports/{sample}_trim.json", sample = ALL_NAMES),
   # expand("/home/gmvaz/2025-sebastian-paper/fastp_output/trim/reports/{sample}_trim.html", sample = ALL_NAMES),
    expand("/home/gmvaz/2025-sebastian-paper/mapped_reads/stats/{sample}_mapping_stats.tsv", sample=ALL_NAMES),
   # expand("/home/gmvaz/2025-sebastian-paper/mapped_reads/stats/{sample}.sorted.bam.bai", sample=ALL_NAMES),
   # expand("/home/gmvaz/2025-sebastian-paper/blastp_results/{sample}_blast_filtered.csv", sample=ALL_NAMES),
   # "/home/gmvaz/2025-sebastian-paper/prodigal_output/combined_prodigal.faa",
   # "/home/gmvaz/2025-sebastian-paper/blast_dbs/prodigal_db.pin",
   # "/home/gmvaz/2025-sebastian-paper/blast_dbs/prodigal_db.phr",
   # "/home/gmvaz/2025-sebastian-paper/blast_dbs/prodigal_db.psq",
  #  expand("/home/gmvaz/2025-sebastian-paper/mapped_reads/{sample}.sorted.bam", sample = ALL_NAMES),
  
# Add R1 and R2 dynamically 
# Ask Colton about no phred score qc, detecting adapter seq, minimum length req, see other snkfile

#make rules for getting ibd fastq samples from SRA 
rule fastp:
  input:
    fq1 = "/group/ctbrowngrp4/2025-ccbaumler-ldtartrate/test_data/{sample}_1.fastq.gz",
    fq2 = "/group/ctbrowngrp4/2025-ccbaumler-ldtartrate/test_data/{sample}_2.fastq.gz"
  output: 
    trim1 = "/home/gmvaz/2025-sebastian-paper/fastp_output/trim/{sample}_1_trim.fastq",
    trim2 = "/home/gmvaz/2025-sebastian-paper/fastp_output/trim/{sample}_2_trim.fastq",
    json = "/home/gmvaz/2025-sebastian-paper/fastp_output/trim/reports/{sample}_trim.json",  
    html = "/home/gmvaz/2025-sebastian-paper/fastp_output/trim/reports/{sample}_trim.html"
  conda: "/home/gmvaz/2025-sebastian-paper/sebastian_env.yml"
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
    trim1 = "/home/gmvaz/2025-sebastian-paper/fastp_output/trim/{sample}_1_trim.fastq",
    trim2 = "/home/gmvaz/2025-sebastian-paper/fastp_output/trim/{sample}_2_trim.fastq",
  output:
    html1 = "/home/gmvaz/2025-sebastian-paper/fastqc_reports/{sample}_1_fastqc.html",
    html2 = "/home/gmvaz/2025-sebastian-paper/fastqc_reports/{sample}_2_fastqc.html",
    zip1 = "/home/gmvaz/2025-sebastian-paper/fastqc_reports/{sample}_1_fastqc.zip",
    zip2 = "/home/gmvaz/2025-sebastian-paper/fastqc_reports/{sample}_2_fastqc.zip",
  conda: "/home/gmvaz/2025-sebastian-paper/sebastian_env.yml"
  shell:
    """
    fastqc \
      {input.trim1} \
      {input.trim2} \
      --outdir {FASTQC_DIR}
    """

rule megahit:
  input:
    trim1 = "/home/gmvaz/2025-sebastian-paper/fastp_output/trim/{sample}_1_trim.fastq",
    trim2 = "/home/gmvaz/2025-sebastian-paper/fastp_output/trim/{sample}_2_trim.fastq",
  output:
    contigs = "/home/gmvaz/2025-sebastian-paper/megahit_output/{sample}/final_contigs.fa",
  conda: "/home/gmvaz/2025-sebastian-paper/sebastian_env.yml"
  threads: 16
  shell:
    """
    megahit \
    -1 {input.trim1} \
    -2 {input.trim2} \
    -o {MEGAHIT_DIR}/{wildcards.sample} \
    -t {threads}
    """    

rule prodigal:
  input:
    contigs = "/home/gmvaz/2025-sebastian-paper/megahit_output/{sample}/final_contigs.fa",
  output:
    fna = "/home/gmvaz/2025-sebastian-paper/prodigal_output/{sample}_prod.fna",
    faa = "/home/gmvaz/2025-sebastian-paper/prodigal_output/{sample}_prod.faa",
    gbk = "/home/gmvaz/2025-sebastian-paper/prodigal_output/{sample}_prod.gbk",
  conda: "/home/gmvaz/2025-sebastian-paper/sebastian_env.yml"
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
    faa = expand("/home/gmvaz/2025-sebastian-paper/prodigal_output/{sample}_prod.faa", sample = ALL_NAMES),
  output:
    combined_faa = "/home/gmvaz/2025-sebastian-paper/prodigal_output/combined_prodigal.faa"
  conda: "/home/gmvaz/2025-sebastian-paper/sebastian_env.yml"
  shell:
    """
    for file in {input.faa}; do
    awk -v filename="$(basename $file)" '/^>/{{print ">" filename "_" substr($0, 2)}} !/^>/' "$file" >> {output.combined_faa}
done
    """

rule blastp_database:
  input:
    combined_faa = "/home/gmvaz/2025-sebastian-paper/prodigal_output/combined_prodigal.faa"
  output:
    pin = "/home/gmvaz/2025-sebastian-paper/blast_dbs/prodigal_db.pin",
    phr = "/home/gmvaz/2025-sebastian-paper/blast_dbs/prodigal_db.phr",
    psq = "/home/gmvaz/2025-sebastian-paper/blast_dbs/prodigal_db.psq"
  params:
    out_prefix = "blast_dbs/prodigal_db"
  conda: "/home/gmvaz/2025-sebastian-paper/sebastian_env.yml"
  shell:
    """
    makeblastdb -in {input.combined_faa} -dbtype prot -out {params.out_prefix}
    """
    
rule blastp: # should we make a bastn rule and a blastn_db?
  input:
    faa = "/home/gmvaz/2025-sebastian-paper/prodigal_output/{sample}_prod.faa",
    db = "/home/gmvaz/2025-sebastian-paper/blast_dbs/prodigal_db.psq"
  output:
     csv = "/home/gmvaz/2025-sebastian-paper/blastp_results/{sample}_blast_output.csv"
  conda: "/home/gmvaz/2025-sebastian-paper/sebastian_env.yml"
  shell:
    """
    blastp \
      -query {input.faa} \
      -db /home/gmvaz/2025-sebastian-paper/blast_dbs/prodigal_db \
      -outfmt "10 qseqid sseqid pident length evalue bitscore" \
      -out temp_{wildcards.sample}_blast.txt

    cat <(printf "qseqid,sseqid,pident,length,e_value,bit_score\\n") \
        temp_{wildcards.sample}_blast.txt > {output.csv}
    """
 
rule parse_blast_output:
  input:
    csv = "/home/gmvaz/2025-sebastian-paper/blastp_results/{sample}_blast_output.csv"
  output:
    filtered = "/home/gmvaz/2025-sebastian-paper/blastp_results/{sample}_blast_filtered.csv"
  params:
    cutoff = "1e-100"
  conda: "/home/gmvaz/2025-sebastian-paper/sebastian_env.yml"
  shell:
    """
    ./parse_blast_output.py {input.csv} \
      -o {output.filtered} \
      --cutoff {params.cutoff} \
      --verbose
    """
    # look at this py script to see what output files. I need a contigs list fna file for downstream rules
rule make_contig_list: 
  input:
    blastp_hits = "/home/gmvaz/2025-sebastian-paper/blastp_results/{sample}_blast_filtered.csv"
  output:
    contig_list_csv = "/home/gmvaz/2025-sebastian-paper/contiglists/{sample}_contigs.csv",
    contig_list_fna = "/home/gmvaz/2025-sebastian-paper/contiglists/{sample}_contigs.fna"
  params:
    fna_dir = "/home/gmvaz/2025-sebastian-paper/prodigal_output/", 
    fna_ext = "fna"
  conda: "/home/gmvaz/2025-sebastian-paper/sebastian_env.yml"
  shell:
    """
    ./extract_prodigal_contigs.py {input.blastp_hits} \
      -v \
      -o {output.contig_list_csv} \
      -d {params.fna_dir} \
      -f {params.fna_ext} \
      --outdir $(dirname {output.contig_list_csv})
    """

rule minimap2:
  input:
    contig_list_fna = "/home/gmvaz/2025-sebastian-paper/contiglists/{sample}_contigs.fna",
    fq1 = "/group/ctbrowngrp4/2025-ccbaumler-ldtartrate/test_data/{sample}_1.fastq.gz",
    fq2 = "/group/ctbrowngrp4/2025-ccbaumler-ldtartrate/test_data/{sample}_2.fastq.gz"
  output:
    sorted_bam = "/home/gmvaz/2025-sebastian-paper/mapped_reads/{sample}.sorted.bam"
  threads: 8
  conda: "/home/gmvaz/2025-sebastian-paper/sebastian_env.yml"
  shell:
    """
    minimap2 -ax sr -t {threads} \
      {input.contig_list_fna} \
      {input.fq1} {input.fq2} \
    | samtools sort -@ {threads} -o {output.sorted_bam}
    """
   
rule index_bam:
  input:
    sorted_bam = "/home/gmvaz/2025-sebastian-paper/mapped_reads/{sample}.sorted.bam"
  output:
    sorted_bai = "/home/gmvaz/2025-sebastian-paper/mapped_reads/stats/{sample}.sorted.bam.bai"
  conda: "/home/gmvaz/2025-sebastian-paper/sebastian_env.yml"
  shell:
    """
    samtools \
    index {input.sorted_bam}
    """

rule reads_mapped:
  input:
    sorted_bam = "/home/gmvaz/2025-sebastian-paper/mapped_reads/{sample}.sorted.bam"
  output:
    stats = "/home/gmvaz/2025-sebastian-paper/mapped_reads/stats/{sample}_mapping_stats.tsv"
  conda: "/home/gmvaz/2025-sebastian-paper/sebastian_env.yml"
  shell:
    """
      (echo -e "contig\tlength\tmapped\tunmapped\tnet_mapped"; \
    samtools idxstats {input.sorted_bam} | \
    awk -F'\t' '{{print $0 "\t" ($3 - $4)}}')
    """