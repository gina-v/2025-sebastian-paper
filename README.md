# Workflow to help Sebastian with Adherent  E. Coli and S. Typhiomerium

Here we will annotate the workflow steps that must be done to create a group of contigs for blasting purposes
 
But, before we do anything lets look at conda....

go to your `.condarc` file that stores your personal configuration for your conda install:

```
nano ~/.condarc
```

```
channels:
  - conda-forge
  - bioconda
  - defaults
envs_dirs:
  - /share/apps/conda/environments
pkgs_dirs:
  - /share/apps/conda/pkgs
  - ~/.conda/pkgs

channel_priority: strict
```

I started by downloading `megahit`, `fastqc`, and `fastp` into the `sebastian.yml` env.

The code below creates a yaml config file for the env name sebastian:
```
conda env export --from-history --name sebastian > "sebastian_env.yml"
```

To install:
```
conda env create --file sebastian_env.yml
```

# Metadata manipulation

To count all the variables in a tsv file for the diagnosis column:
```
awk -F'\t' '{count[$8]++} END {for (val in count) print val, count[val]}' metadata/working_metadata.tsv
```
which returns:
```
 1
diagnosis 1
CD 566
UC 122
nonIBD 372
```

To see the unique values in a tsv file for the projects column:
```
awk -F'\t' '{print $1}' metadata/working_metadata.tsv | sort | uniq
```
The values returned are:
```
PRJNA237362
PRJNA237362
PRJEB2054
PRJNA237362
PRJNA385949
PRJNA400072
SRP057027
study_accession
```

Notice there are a few repeated values... that's unexpected.
```
awk -F'\t' '{print $1}' metadata/working_metadata.tsv | sort | uniq -c | sort -nr
```
The repeated unique values must contain a whitespace character (or some other invisible unicode character) somewhere in their string:
```
    448 PRJNA237362
    263 PRJEB2054
    218 PRJNA400072
    112 SRP057027
     17 PRJNA385949
      1 study_accession
      1 PRJNA237362
      1 PRJNA237362
      1
```


# Some helpful bash commands

Find all the files that are greater that 64 bits in size. Count them all.
```
find /group/ctbrowngrp/2025-ccbaumler-team-tartrate/mapped_reads/stats/ -size +64 -print | wc -l
```

List all the files in a human readable format by the time modified in reverse order:
```
ls -lathr /group/ctbrowngrp/2025-ccbaumler-team-tartrate/mapped_reads/stats/
```

Return how much diskspace is currently available:
```
df -h /group/ctbrowngrp 
```

Update the timestamp on all the files in a directory and its sub-directories:
```
find /group/ctbrowngrp/2025-ccbaumler-team-tartrate/contiglists/ -exec touch {} +
```

Look at the queue of all jobs running in my group but exclude my running jobs:
```
(squeue -A ctbrowngrp -o "%.10i %.10u %.10P %.20j %.2t %.10M %.4C %.6D %.10m %R" | head -n 1 && squeue -A ctbrowngrp -o "%.10i %.10u %.10P %.20j %.2t %.10M %.4C %.6D %.10m %R"| awk -v user="$USER" 'NR>1 && $2 != user')
```

Watch the changes in my squeue:
```
watch 'squeue --me -o "%.10i %.10P %.20j %.2t %.10M %.4C %.6D %.10m %R"'
```

# Gina, here is a workflow!

1. Download the samples
I think you did this manually, but you could try to automate this by parsing a metadata file of accessions

2. Run `fastp` to clean up metagenomes

```
fastp \
    -i "$R1" \
    -I "$R2" \
    -o "$OUTDIR/${base}_R1_fastp.fq.gz" \
    -O "$OUTDIR/${base}_R2_fastp.fq.gz" \
    -h "$OUTDIR/${base}_fastp.html" \
    -j "$OUTDIR/${base}_fastp.json"
```

3. Run `megahit` and `fastqc` (should probably add `multiqc`)

    a.
        ```
        megahit \
            -1 "$R1" \
            -2 "$R2" \
            -o "$OUTDIR/${base}" \
            -t $SLURM_CPUS_PER_TASK
        ```

     b. 
        ```
        fastqc \
            "$R1" \
            "$R2" \
            --outdir "$OUTDIR${base}"
        ```

4. Run `prodigal`

```
# Metagenome mode, # Treat runs of N as masked sequence; dont build genes across them.
prodigal \
    -p meta \
    -m \
    -i "${CONTIGS[$sample]}" \
    -d "$OUTDIR${base}.prod.fna" \
    -a "$OUTDIR${base}.prod.faa" \
    -o "$OUTDIR${base}.prod.gbk"
```

5. Collate the prodigal output with simple bash cmd

This collects all the prodigal sequences for all the samples we are using (currently 5 sequences)
```
# Loop over each Prodigal output file
for file in /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/prodigal-output/*.faa; do
    # Add the file name as a prefix to each header
    awk -v filename="$(basename $file)" '/^>/{print ">" filename "_" substr($0, 2)} !/^>/' "$file" >> combined_prodigal_with_fp.faa
done
```

6. Create a BLAST DB from all the prodigal sequences

Make an sample type (ibd) database
```
makeblastdb -in combined_prodigal_with_fp.faa -dbtype prot -out blast_dbs/prodigal_ibd_db
```


To use a better output format and add a header:
```
blastp -query seqs/STM3354_1254877_data/protein.faa -out blast_output -db blast_dbs/prodigal_ibd_db -outfmt "10 qseqid sseqid pident length evalue bitscore" && cat <(printf "qseqid,sseqid,pident,length,e_value,bit_score\n") blast_output > STM3354_blast_output.csv
```

7. Parse the BLAST results to trim the unwanted and curate the seqfile and contig ids

```
./parse_blast_output.py blast_output_header -o test.csv --cutoff 1e-100 --verbose
```
> `--cutoff` is the minimum allowed e-value in the output. I think 1e-7 is a decent, general cutoff to indicate that this result is not random. The e-value, or estimate value, is an estimate of the number of alignments with a better score that would be expected to occur by chance. A measurement of 1e-7 means only 1 alignment in 10 million random comparisons.
> `-o` will save a csv file separating the seq_filename and sseqid

8. To extract the contigs from the appropriate prodigal file and filetype.

```
./extract_prodigal_contigs.py test.csv -v -o trial.csv -d /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/prodigal-output/ -f fna --outdir /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/contiglists/
```
> `-v` for verbose output
> `-o` for the csv file containing the dictionary like object
> `-d` for the directory path to the prodigal output that has already been generated
> `-f` for the prodigal filetype to extract the contig to store in a similar filetype

9. Map the contig list to the pair of raw fasta files containing the reads

```
minimap2 -ax sr -t 8 /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/contiglists/ERR13295994.x.NP_462264.1.contiglist.fna /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/test_data/ERR13295994_1.fastq.gz /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/test_data/ERR13295994_2.fastq.gz | samtools sort -@ 8 -o /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/mapped_reads/ERR13295994.x.NP462264.1_mapped_reads.sorted.bam
```

10. Index and return the info we came for

index the binary alignment file
```
samtools index /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/mapped_reads/ERR13295994.x.NP462264.1_mapped_reads.sorted.bam 
```

The reads mapped to each contig returned by the BLAST search
```
(echo -e "contig\tlength\tmapped\tunmapped\tnet_mapped"; \
 samtools idxstats /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/mapped_reads/ERR13295994.x.NP462264.1_mapped_reads.sorted.bam | \
 awk -F'\t' '{print $0 "\t" ($3 - $4)}')
```

Per-base read depth (coverage?) across contig
```
samtools depth -aa /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/mapped_reads/ERR13295994.x.NP462264.1_mapped_reads.sorted.bam >> /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/ERR13295994.x.NP462264.1_coverage.tsv
```

Find the average coverage across the individual contig sequence
```
./aggregating_depth.py /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/ERR13295994.x.NP462264.1_coverage.tsv /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/ERR13295994.x.NP462264.1_agg_coverage.tsv
```

---

# Prodigal to BLAST note

I needed to install both the biopython module and the blast program (to make the local db)
```
conda install biopython blast
```

I use the `cat` cli tool to combine all the `faa` files into one for the blastdb.
```
cat /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/prodigal-output/*.prod.faa > all_prod.ibd.faa
```
or
```
# Loop over each Prodigal output file
for file in /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/prodigal-output/*.faa; do
    # Add the file name as a prefix to each header
    awk -v filename="$(basename $file)" '/^>/{print ">" filename "_" substr($0, 2)} !/^>/' "$file" >> combined_proteins.faa
done
```

```
makeblastdb -in all_prod.ibd.faa -dbtype prot -out combined_prodigal_db
```

```
blastp -query seqs/STM3354_1254877_data/protein.faa -out blast_test -db blast_dbs/prodigal_db_with_filenames
```

To use a better output format and add a header:
```
blastp -query seqs/STM3354_1254877_data/protein.faa -out blast_output -db blast_dbs/prodigal_db_with_filenames -outfmt "10 qseqid sseqid pident length evalue bitscore" && cat <(printf "qseqid,sseqid,pident,length,e_value,bit_score\n") blast_output > blast_output_header
```


# Notes on mapping reads to contigs

Makes sure that the contig is nucleotide and not amino acid. The DNA (or fna file) contig is the same header as the protein (or faa file) contig!

minimap is by far faster than bwa, but bwa seems to find more reads. Are bwa's read erroneous?

I used 8 cores in an srun to get results quickly but maybe 16 will be even faster.

```
minimap2 -ax sr -t 8 contigs.fasta raw_reads.fastq \
  | samtools sort -@ 8 -o mapped_reads.sorted.bam
```

A real example looks like:
```
minimap2 -ax sr -t 8 /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/contiglists/ERR13295994.contiglist.fna /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/test_data/ERR13295994_1.fastq.gz /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/test_data/ERR13295994_2.fastq.gz | samtools sort -@ 8 -o /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/ERR13295994_mapped_reads.sorted.bam
```


I think indexing and sorting speeds things up
```
samtools index mapped_reads.sorted.bam
```

samtools has multiple stat functions:

total mapping stats:

```
samtools flagstat mapped_reads.sorted.bam
```

All the extra stats like average read length, insert size, GC content:
```
samtools stats mapped_reads.sorted.bam
```

For us, we only need to know per-contig read counts:
```
samtools idxstats mapped_read.sorted.bam
```

For per-base depth
```
samtools depth -aa mapped_read.sorted.bam
```

# Scripts to make BLAST output and contig mapping easier

1. Make the blast output easier to parse 

NOTE! Make sure to add the filename to the prodigal results
```
# Loop over each Prodigal output file
for file in /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/prodigal-output/*.faa; do
    # Add the file name as a prefix to each header
    awk -v filename="$(basename $file)" '/^>/{print ">" filename "_" substr($0, 2)} !/^>/' "$file" >> combined_proteins.faa
done
```

2. Parse the BLAST results to trim the unwanted and curate the seqfile and contig ids

```
./parse_blast_output.py blast_output_header -o test.csv --cutoff 1e-100 --verbose
```
> `--cutoff` is the minimum allowed e-value in the output. I think 1e-7 is a decent, general cutoff to indicate that this result is not random. The e-value, or estimate value, is an estimate of the number of alignments with a better score that would be expected to occur by chance. A measurement of 1e-7 means only 1 alignment in 10 million random comparisons.
> `-o` will save a csv file separating the seq_filename and sseqid

3. To extract the contigs from the appropriate prodigal file and filetype.

```
./extract_prodigal_contigs.py test.csv -v -o trial.csv -d /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/prodigal-output/ -f fna --outdir /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/contiglists/
```
> `-v` for verbose output
> `-o` for the csv file containing the dictionary like object
> `-d` for the directory path to the prodigal output that has already been generated
> `-f` for the prodigal filetype to extract the contig to store in a similar filetype

4. Map the contig list to the pair of raw fasta files containing the reads

```
minimap2 -ax sr -t 8 /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/contiglists/ERR13295994.x.NP_462264.1.contiglist.fna /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/test_data/ERR13295994_1.fastq.gz /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/test_data/ERR13295994_2.fastq.gz | samtools sort -@ 8 -o /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/mapped_reads/ERR13295994.x.NP462264.1_mapped_reads.sorted.bam
```

5. Index and return the info we came for

index the binary alignment file
```
samtools index /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/mapped_reads/ERR13295994.x.NP462264.1_mapped_reads.sorted.bam 
```

for the full breakdown of mapped reads per contig returned by the BLAST search
```
(echo -e "contig\tlength\tmapped\tunmapped\tnet_mapped"; \
 samtools idxstats /group/ctbrowngrp4/2025-ccbaumler-ldtartrate/mapped_reads/ERR13295994.x.NP462264.1_mapped_reads.sorted.bam | \
 awk -F'\t' '{print $0 "\t" ($3 - $4)}')
```

# After everything is mapped this may help

[For an explanation of the output of `samtools idxstats`](https://toolshed.g2.bx.psu.edu/repository/display_tool?repository_id=612bfd208a70f3c3&render_repository_actions_for=tool_shed&tool_config=%2Fsrv%2Ftoolshed%2Fmain%2Fvar%2Fdata%2Frepos%2F000%2Frepo_692%2Ftools%2Fsamtools_idxstats%2Fsamtools_idxstats.xml&changeset_revision=71afa65f444a)

`samtools idxstats` tab-delimited column output descriptions:
| Column | Name       | Meaning                                                                                                       |
| ------ | ---------- | ------------------------------------------------------------------------------------------------------------- |
| **1**  | `contig`   | The reference sequence name (from the FASTA used to align the reads) — e.g. a contig, chromosome, or scaffold |
| **2**  | `length`   | The length of the reference sequence in bases                                                                 |
| **3**  | `mapped`   | The number of reads **mapped** to this contig                                                                 |
| **4**  | `unmapped` | The number of **unmapped** reads **whose mate** mapped to this contig (in paired-end data)                    |

Adding header and subtracting the mapped reads from the reads without a mapped pair:
```
(echo -e "contig\tlength\tmapped\tunmapped\tnet_mapped"; \
 samtools idxstats mapped_reads/ERR13295994_mapped_reads.sorted.bam | \
 awk -F'\t' '{print $0 "\t" ($3 - $4)}')
```
> add output with `> outputfile.tsv` 
