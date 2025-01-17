ssh inf-51-2023@bioinf-serv2.cob.lu.se
# Snakemake (https://snakemake.readthedocs.io) 
# build an index for bowtie2
bowtie2-build yeastGenome.fa yeastGenome
# map reads using bowtie2
bowtie2 -x yeastGenome -U yeastReads.fastq -S readsAligned.sam
# sort the file and convert to bam
samtools sort -o readsAligned.sortedPositions.bam readsAligned.sam
# index the sorted bam file
samtools index readsAligned.sortedPositions.bam
# pileup variants
bcftools mpileup -Ob -f yeastGenome.fa \
readsAligned.sortedPositions.bam \ > variants.bcf
# call variants
bcftools call -m -v variants.bcf > variants.vcf
# Beyond this, you also spent some time installing the correct versions of different
# softwares, etc. For basic reproducibility, you could essentially just provide a file
# containing the above script, together with the necessary data files and software
# versions, and somebody could rerun this and get the same results you did. How-
# ever, if you wanted to align this against a different genome, or align different
# read sets to the same genome, you would need to edit the script file manually,
# both to avoid rerunning the indexing, but also to update the file names in each
# of the steps. Snakemake helps us with this as well.
# Snakemake
# Install and activate Snakemake by running
conda create -n snakemake -c conda-forge -c bioconda snakemake  # -c means channel
conda activate snakemake
# 2.2 Snakemake rules
# A rule may look like this:
rule sam_to_bam: # sam_to_bam is the rule name
input:
sam = "readsAligned.sam"
output:
bam = "readsAligned.sortedPositions.bam"
shell:
"samtools sort "
"-o {output.bam} "
"{input} "
# input and output files are referenced within the curly brackets.
# introducing wildcards into the names
# of input and output files
rule sam_to_bam:
input:
sam = "{sample}.sam"
output:
bam = "{sample}.sortedPositions.bam"
shell:
"samtools sort "
"-o {output.bam} "
"{input} "
#the rule above could be used to turn any file ending in .sam into a sorted BAM file
# file name would be identical to the input file name, but replacing the .sam file ending with 
# a .sortedPositions.bam file ending
rule index_bam:
input:
bam = "{sample}.sortedPositions.bam"
output:
index = "{sample}.sortedPositions.bam.bai"
shell:
"samtools index {input} "
# 2.3 Snakemake target files
# first rule of the Snakefile is named all
rule all:
input:
'variants.vcf'
# rule all should specify one or a list of files to be generated, as opposed to
# a wildcard, as a wildcard on its own only provides a file name pattern.
# If we are dealing with multiple samples, we may want to use the expand function
# The expand function takes a python string and performs substitutions
# using the python formatting mini-language: https://docs.python.org/3/librar
# y/string.html#formatspec. 
# The expand function returns a list of strings
SAMPLE = ['sample_1', 'sample_2']
rule all:
input:
expand('{sample}.vcf', sample=SAMPLE) 
# curly brackets in the first argument to expand is not a Snakemake wildcard, but simply a placeholder
# gets replaced by whatever is provided as other named arguments to expand, in this case SAMPLE
# the input to rule all is one or several names of actual files that are to be generated by the workflow
# The SAMPLE variable here can be created using any valid python code.
# providing a config tsv file containing sample names 
configfile: 'config/config.yaml' 
sample_df = (pd.read_csv(config['samples'],
sep='\t',
dtype={'sample_name':str, 'fastq':str})
.set_index('sample_name', drop=False))
rule all:
input:
expand('{sample}.vcf', sample=sample_df.sample_name)
# have a config file located at config/config.yaml. 
samples: 'config/samples.tsv'
# the samples property may be accessed using config['samples'] within the Snakefile. 
# access this file and parse it using the pandas.read_csv function, saving the resulting dataframe to the sample_df variable
# reading the tsv file, we see that it must be tab separated (sep='\t'),
# and have the columns sample_name and fastq, each of which contains string data
sample_name fastq
sample_1 resources/s1.fq
# expand function reads the sample_name column of this dataframe, and creates the list of target files.
# Running the workflow for additional datasets would simply mean adding rows with sample names and fastq file locations to
# the tsv file for each additional sample we want to run
# specifying the samples in a tsv file
# glob_wildcards function to identify wildcards based on file names
SAMPLE, = glob_wildcards('{sample}.fastq')
rule all:
input:
expand('{sample}.vcf', sample=SAMPLE)
# extracts the first part of all fastq file names and saves them to the SAMPLE
# variable as a list, which we may then use with the expand function
# aggregate the output of several different samples, expand may be very useful
mkdir example_workflow
cd example_workflow
# create a Snakefile
touch Snakefile
# create a config directory
snakemake --cores 1 # run the workflow
Add the missing file by touching it. What happens if you run snakemake
-n? Try re-running the workflow. What happens?
# The -n flag is used to perform a dry run, which means that Snakemake will
SAMPLE, = glob_wildcards('{sample}.fastq')

rule all:
    input:
        expand('{sample}.vcf', sample=SAMPLE)

# A full Snakemake pipeline
mkdir snake_pipeline
cd snake_pipeline 
import pandas as pd
import re
configfile: 'config/config.yaml'

sample_df = (pd.read_csv(config['samples'],
sep='\t',
dtype={'sample_name':str, 'fastq':str})
.set_index('sample_name', drop=False))
rule all:
input:
expand('results/01_called_variants/{sample}.vcf',
sample=sample_df.sample_name)
rule bowtie_index:
input:
genome = config['genome']
output:
multiext(config['genome'],
".1.bt2", ".2.bt2", ".3.bt2",
".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
conda: 'envs/bowtie2.yaml'
shell:
"bowtie2-build {input.genome} {input.genome} "
rule map_reads:
input:
idx = rules.bowtie_index.output,
reads = lambda wildcards: sample_df.loc[wildcards.sample,
'fastq']
output:
temp('results/00_mapped_reads/{sample}.unsorted.sam')
params:
idx = config['genome']
conda: 'envs/bowtie2.yaml'
shell:
'bowtie2 -x {params.idx} '
'-U {input.reads} '
'-S {output} '
rule sam_to_bam:
input:
rules.map_reads.output
output:
'results/00_mapped_reads/{sample}.bam'
threads: 2
conda: 'envs/htslib.yaml'
shell:
'samtools sort '
'-@ {threads} '
'-o {output} {input} '
rule index_bam:
input:
rules.sam_to_bam.output
output:
'results/00_mapped_reads/{sample}.bam.bai'
conda: 'envs/htslib.yaml'
shell:
'samtools index {input} '

rule call_variants:
input:
rules.index_bam.output,
aligned_reads = rules.sam_to_bam.output,
genome = config['genome']
output:
'results/01_called_variants/{sample}.vcf'
conda: 'envs/htslib.yaml'
shell:
'bcftools mpileup -Ou '
'-f {input.genome} '
'{input.aligned_reads} '
'| bcftools call -m -v '
'> {output} '
snakemake -p --use-conda -j2
# -j controls the concurrency (number of parallel jobs), and -p enables the printing of shell commands before they are executed. These options provide 
# flexibility and transparency when running and debugging Snakemake workflows.
# set up the yaml files in the workflow/envs directory
rule NAME:
    input:
        "table.txt"
    output:
        "plots/myplot.pdf"
    conda:
        "envs/ggplot.yaml"
    script:
        "scripts/plot-stuff.R"
scp -r /Users/med-snt/snake_pipeline/Snakefile inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/snake_pipelineS
cp /home/inf-51-2023/Variant_calling/yeastGenome.fa .
cp /home/inf-51-2023/Variant_calling/yeastReads.fastq .
/home/inf-51-2023/snake_pipelineS/resources/yeastReads.fastq
conda install mamba
tree
snakemake -p --use-conda -j2
/home/inf-51-2023/snake_pipelineS/results/01_called_variants/sample_1.vcf
less sample1.vcf
# actual variants one per line
cat sample1.vcf | grep "^#CHROM" | tr "\t" "\n" | cat -n
# The VCF specification is available here:
# https://samtools.github.io/hts-specs/VCFv4.2.pdf
# question 10
snakemake --cores 1 --forceall # before adding the new sample
snakemake -n --forceall # adding the new sample
# question 11:
# Download fastq (in hooi environment)
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp
# how to check if fastp is installed
which fastp
# move fastp to bin
mv fastp ~/bin
# install fastqc
conda install -c bioconda fastqc  
# You can add a directory to your PATH by adding a line like this to your `.bashrc` or `.bash_profile` file: 
export PATH=$PATH:/path/to/directory
# Replace `/path/to/directory` with the actual directory where FastQC is installed. 
# After adding this line, you'll need to restart your terminal or run `source ~/.bashrc` or `source ~/.bash_profile` for the changes to take effect.
# 11. Add rules for running fastqc and fastp (https://github.com/OpenGene/fastp). 
# Ensure that the fastqc step is run for reads both before and
# after trimming. The trimming step can be done using the default fastp
# settings, but make sure that reads are trimmed prior to mapping

fastqc /home/inf-51-2023/snake_pipelineS/resources/yeastReads.fastq > qc_b4.fastq
cp /home/inf-51-2023/snake_pipelineS/resources/yeastReads.fastq .
# create another copy of the yeastReads.fastq file
cp ./yeastReads.fastq ./raw_reads.fastq 
fastp -i raw_reads.fastq -o trimmed_reads.fastq
fastqc trimmed_reads.fastq
# copy file to the server
scp -r /Users/med-snt/Downloads/Snakefile inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/fastp_snakemake

# the orinal file of yeastReads.fastq is in the resources folder
/home/inf-51-2023/snake_pipelineS/resources/yeastReads.fastq 
head raw_reads.fastq
head trimmed_reads.fastq
# To search for installed software from multiple environments
#!/bin/bash

# Get a list of all conda environments
environments=$(conda env list | awk '{print $1}' | tail -n +3)

for env in $environments
do
  # Activate the environment
  conda activate $env

  # Check if fastqc is installed in this environment
  if which fastqc > /dev/null; then
    echo "fastqc is installed in the $env environment"
  else
    echo "fastqc is not installed in the $env environment"
  fi

  # Deactivate the environment before moving on to the next one
  conda deactivate
done
# 
rule raw_fastqc:
    input:
        reads = lambda wildcards: sample_df.loc[wildcards.sample, "fastq"]
    output:
        "resources/fastqc_{sample}.html"
    conda: "envs/fastqc.yaml"
    shell:
        "fastqc {input.reads} "

rule trim:
    input:
        raw = rules.raw_fastqc.input.reads
    output:
        "resources/{sample}-trim.fastq"
    conda: "envs/fastp.yaml"
    shell:
        "fastp -i {input.raw} -o {output}"

rule trim_fastqc:
    input:
        trim = rules.trim.output
    output:
        "resources/fastqc_{sample}-trim.html"
    conda: "envs/fastqc.yaml"
    shell:
        "fastqc {input.trim}"
# copy the fastp.yaml file to the envs folder
scp -r /Users/med-snt/Downloads/fastp.yaml inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/fastp_snakemake/workflow/envs
# copy the fastqc.yaml file to the envs folder 
scp -r /Users/med-snt/Downloads/fastqc.yaml inf-51-2023@bioinf-serv2.cob.lu.se:/home/inf-51-2023/fastp_snakemake/workflow/envs
# copy bowtie2.yaml and htslib file to the envs folder (fastp_snakemake/workflow/envs)
cp /home/inf-51-2023/snake_pipelineS/workflow/envs/bowtie2.yaml .
cp /home/inf-51-2023/snake_pipelineS/workflow/envs/htslib.yaml .
# copy the config.yaml and samples.tsv files to the config folder (fastp_snakemake/config)
cp /home/inf-51-2023/snake_pipelineS/config/config.yaml .
cp /home/inf-51-2023/snake_pipelineS/config/samples.tsv .
# run snakemake
snakemake -p --use-conda -j2
Complete log: .snakemake/log/2024-01-28T190225.994472.snakemake.log
.
├── config
│   ├── config.yaml
│   └── samples.tsv
├── fastp.html
├── fastp.json
├── resources
│   ├── genome
│   │   ├── yeastGenome.fa
│   │   ├── yeastGenome.fa.1.bt2
│   │   ├── yeastGenome.fa.2.bt2
│   │   ├── yeastGenome.fa.3.bt2
│   │   ├── yeastGenome.fa.4.bt2
│   │   ├── yeastGenome.fa.fai
│   │   ├── yeastGenome.fa.rev.1.bt2
│   │   └── yeastGenome.fa.rev.2.bt2
│   ├── samp1-trim.fastq
│   ├── sample1-trim.fastq
│   └── yeastReads.fastq
├── results
│   ├── 00_mapped_reads
│   │   ├── samp1.bam
│   │   ├── samp1.bam.bai
│   │   ├── sample1.bam
│   │   └── sample1.bam.bai
│   └── 01_called_variants
│       ├── samp1.vcf
│       └── sample1.vcf
└── workflow
    ├── envs
    │   ├── bowtie2.yaml
    │   ├── fastp.yaml
    │   ├── fastqc.yaml
    │   └── htslib.yaml
    └── Snakefile
# dry run to check how many jobs will be run
snakemake -n --forceall 
# check sample1.vcf file after trimming (quality control and word count)

