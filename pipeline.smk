from datetime import datetime
import os

shell.prefix("source ~jgomez/init_shell.sh;")

if not os.path.exists("logs"):
   os.makedirs("logs")
   
date = datetime.now().strftime('%Y%m%d.%H%M%S')

#IDS = "1 2 3 ...".split() # the list of desired ids
#també podem fer una wildcard global IDS, = glob_wildcards("thedir/{id}.fastq")
IDS, = glob_wildcards("/scratch/devel/imoscardo/snakemake/reads/{sample}.fastq.gz")

rule all:
    input:
        "/scratch/devel/imoscardo/snakemake/annotation/assembly/assembly_annotation/cmsearch.out"
    log:
        "logs/" + str(date) + ".j%j.all.out",
        "logs/" + str(date) + ".j%j.all.err"
	
rule trim_galore:
    input:
        read = "/scratch/devel/imoscardo/snakemake/reads/{sample}.fastq.gz"
    output:
        trimmed_read = "/scratch/devel/imoscardo/snakemake/alignment/trim_galore/{sample}_trimmed.fq.gz"
    params:
        min_length = 16,
        max_length = 200,
        trimming = "/scratch/devel/imoscardo/snakemake/alignment/trim_galore"
    log:
        "logs/" + str(date) + ".j%j.trim_galore.{sample}.out",
        "logs/" + str(date) + ".j%j.trim_galore.{sample}.err"
    threads: 1
    shell:
        "conda activate ~jgomez/conda_environments/preprocess_illumina;"
	"trim_galore --length {params.min_length} --max_length {params.max_length} -o {params.trimming} {input.read}"
        
rule index_STAR:
    input:
        genome = "/scratch/devel/imoscardo/data/genome/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa.gz"
    output:
        index = directory("/scratch/devel/imoscardo/snakemake/alignment/genome"),
	genome = "/scratch/devel/imoscardo/snakemake/genome/genome.fa" 
    log:
        "logs/" + str(date) + ".j%j.index_STAR.out",
        "logs/" + str(date) + ".j%j.index_STAR.err"
    threads: 4
    shell:
        "mkdir /scratch/devel/imoscardo/snakemake/alignment/genome;"
        "gunzip -c {input.genome} > {output.genome};"
        "module purge;"
        "module load gcc/6.3.0;"
        "module load STAR/2.7.2a;"
        "STAR --runMode genomeGenerate --runThreadN {threads} --genomeDir {output.index} --genomeFastaFiles {output.genome}"

rule STAR:
    input:
        trimmed_read = "/scratch/devel/imoscardo/snakemake/alignment/trim_galore/{sample}_trimmed.fq.gz",
        index = "/scratch/devel/imoscardo/snakemake/alignment/genome"
    output:
        aligned_read = "/scratch/devel/imoscardo/snakemake/alignment/STAR/{sample}_Aligned.sortedByCoord.out.bam"
    params:
        prefix = "{sample}_",
        max_intron = 1        
    log:
        "logs/" + str(date) + ".j%j.STAR.{sample}.out",
        "logs/" + str(date) + ".j%j.STAR.{sample}.err"
    threads: 4
    shell:
        "cd /scratch/devel/imoscardo/snakemake/alignment/STAR;"
        "module purge;"
        "module load gcc/6.3.0;"
        "module load STAR/2.7.2a;"
        "STAR --genomeDir {input.index}  --readFilesCommand gunzip -c  --readFilesIn {input.trimmed_read} --outFileNamePrefix {params.prefix} --runThreadN {threads} --alignEndsType EndToEnd --outFilterMismatchNmax 1 --outFilterMultimapScoreRange 0 --outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 10 --outSAMunmapped Within --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 16 --alignSJDBoverhangMin 1000 --alignIntronMax {params.max_intron}"

rule StringTie:
     input:
        aligned_read = "/scratch/devel/imoscardo/snakemake/alignment/STAR/{sample}_Aligned.sortedByCoord.out.bam",
        genome = "/scratch/devel/imoscardo/snakemake/genome/genome.fa"
     output:
        assembled_read_stringtie = "/scratch/devel/imoscardo/snakemake/assembly/stringtie/{sample}.gtf"
     params:
        min_lenght = 30
     log:
        "logs/" + str(date) + ".j%j.StringTie.{sample}.out",
        "logs/" + str(date) + ".j%j.StringTie.{sample}.err"
     threads: 4
     shell:
        "/scratch/project/devel/aateam/src/STRINGTIE/stringtie_v2_2_0/stringtie -o {output.assembled_read_stringtie} --ref {input.genome} -m {params.min_lenght} {input.aligned_read}"

rule gtf_list_stringtie:
     input:
        list_files = "/scratch/devel/imoscardo/snakemake/gtf_list/gtf_list_stringtie.py",
        gtfs = lambda wildcards: expand(rules.StringTie.output.assembled_read_stringtie, sample = IDS)
     output:
        gtf_list = "/scratch/devel/imoscardo/snakemake/assembly/stringtie/gtf_list_stringtie.txt"
     log:
        "logs/" + str(date) + ".j%j.gtf_list_stringtie.out",
        "logs/" + str(date) + ".j%j.gtf_list_stringtie.err"
     shell:
        "python {input.list_files}"

rule stringtie_merge:
     input:
        gtf_list = "/scratch/devel/imoscardo/snakemake/assembly/stringtie/gtf_list_stringtie.txt",
	gtfs = lambda wildcards: expand(rules.gtf_list_stringtie.output.gtf_list, sample = IDS)
     output:
        merged_stringtie = "/scratch/devel/imoscardo/snakemake/assembly/stringtie/merged_stringtie.gtf"
     log:
        "logs/" + str(date) + ".j%j.stringtie_merge.out",
        "logs/" + str(date) + ".j%j.stringtie_merge.err"
     threads: 1
     shell:
        "/scratch/project/devel/aateam/src/STRINGTIE/stringtie_v2_2_0/stringtie --merge -o {output.merged_stringtie} {input.gtf_list}"

rule Cufflinks:
     input: 
     	aligned_read = "/scratch/devel/imoscardo/snakemake/alignment/STAR/{sample}_Aligned.sortedByCoord.out.bam",
        genome = "/scratch/devel/imoscardo/snakemake/genome/genome.fa"
     output: 
        assembled_read_cufflinks = directory("/scratch/devel/imoscardo/snakemake/assembly/cufflinks/{sample}")
     params:
        min_lenght = 16
     log:
        "logs/" + str(date) + ".j%j.Cufflinks.{sample}.out",
        "logs/" + str(date) + ".j%j.Cufflinks.{sample}.err"
     threads: 4
     shell:
        "module load gcc/4.9.3-gold;"
        "module load boost/1.55.0;"
        "module load Cufflinks;"
        "cufflinks  -o {output.assembled_read_cufflinks} -b {input.genome} -m {params.min_lenght} {input.aligned_read}"

rule gtf_list_cufflinks:
     input:
        list_files = "/scratch/devel/imoscardo/snakemake/gtf_list/gtf_list_cufflinks.py",
        gtfs = lambda wildcards: expand(rules.Cufflinks.output.assembled_read_cufflinks, sample = IDS)
     output:
        gtf_list = "/scratch/devel/imoscardo/snakemake/assembly/cufflinks/gtf_list_cufflinks.txt"
     log:
        "logs/" + str(date) + ".j%j.gtf_list_cufflinks.out",
        "logs/" + str(date) + ".j%j.gtf_list_cufflinks.err"
     shell:
        "python {input.list_files}"

rule cuffmerge:
     input:
        gtf_list = "/scratch/devel/imoscardo/snakemake/assembly/cufflinks/gtf_list_cufflinks.txt",
        gtfs = lambda wildcards: expand(rules.gtf_list_cufflinks.output.gtf_list, sample = IDS)
     output:
        merged_cufflinks = directory("/scratch/devel/imoscardo/snakemake/assembly/cufflinks/merged_cufflinks.gtf")
     log:
        "logs/" + str(date) + ".j%j.cuffmerge.out",
        "logs/" + str(date) + ".j%j.cuffmerge.err"	
     threads: 1
     shell:
        "cd /scratch/devel/imoscardo/snakemake/assembly/cufflinks;"
        "module load gcc/4.9.3-gold;" 
        "module load boost/1.55.0;"
        "module load module load PYTHON/2.7.17;"
        "module load Cufflinks;"
        "cuffmerge -o {output.merged_cufflinks} {input.gtf_list}"

rule TransBorrow:
     input:
        assembled_read_stringtie = "/scratch/devel/imoscardo/snakemake/assembly/stringtie/{sample}.gtf",
        assembled_read_cufflinks = "/scratch/devel/imoscardo/snakemake/assembly/cufflinks/{sample}/transcripts.gtf",
        genome = "/scratch/devel/imoscardo/snakemake/genome/genome.fa",
        aligned_read = "/scratch/devel/imoscardo/snakemake/alignment/STAR/{sample}_Aligned.sortedByCoord.out.bam"
     output:
        combined_assemblies = "/scratch/devel/imoscardo/snakemake/assembly/transborrow/combined_{sample}.gtf",
	assembled_transborrow = directory("/scratch/devel/imoscardo/snakemake/assembly/transborrow/assembled_{sample}"),
        assembled_transborrow_modified = "/scratch/devel/imoscardo/snakemake/assembly/transborrow/assembled_{sample}/TransBorrow_results/TransBorrow_modified.gtf"
     params:
        strand_info = "single_unstranded",
        min_lenght = 16
     log:
        "logs/" + str(date) + ".j%j.TransBorrow.{sample}.out",
        "logs/" + str(date) + ".j%j.TransBorrow.{sample}.err"
     threads: 4
     shell:
        "module purge;"
        "module load CMAKE/3.19.4;"
        "module unload gcc;"       
        "module load gcc/10.2.0;" 
        "cat {input.assembled_read_stringtie} {input.assembled_read_cufflinks} > /scratch/devel/imoscardo/snakemake/assembly/transborrow/combined_{wildcards.sample}.gtf;"
        "mkdir /scratch/devel/imoscardo/snakemake/assembly/transborrow/assembled_{wildcards.sample};"
        "cd /scratch/devel/imoscardo/snakemake/assembly/transborrow/assembled_{wildcards.sample};"
        "TransBorrow  --min_trans_len {params.min_lenght} -r {output.combined_assemblies} -g {input.genome} -b {input.aligned_read} -s {params.strand_info};"    
        "sed \'s/TPM/FPKM/\' /scratch/devel/imoscardo/snakemake/assembly/transborrow/assembled_{wildcards.sample}/TransBorrow_results/TransBorrow.gtf > /scratch/devel/imoscardo/snakemake/assembly/transborrow/assembled_{wildcards.sample}/TransBorrow_results/TransBorrow_modified.gtf"

rule gtf_list:
     input:
        list_files = "/scratch/devel/imoscardo/snakemake/gtf_list/gtf_list_transborrow.py",
        gtfs = lambda wildcards: expand(rules.TransBorrow.output.assembled_transborrow_modified, sample = IDS)
     output:
        gtf_list = "/scratch/devel/imoscardo/snakemake/assembly/transborrow/gtf_list_transborrow.txt"
     log:
        "logs/" + str(date) + ".j%j.gtf_list.out",
        "logs/" + str(date) + ".j%j.gtf_list.err"
     shell:
        "python {input.list_files}"

rule TACO:
     input:
        gtf_list = "/scratch/devel/imoscardo/snakemake/assembly/transborrow/gtf_list_transborrow.txt",
        gtfs = lambda wildcards: expand(rules.gtf_list.output.gtf_list, sample = IDS),
        genome = "/scratch/devel/imoscardo/snakemake/genome/genome.fa"
     output:
        assembled_read_taco = directory("/scratch/devel/imoscardo/snakemake/assembly/taco")
     params:
        min_lenght = 16
     log:
        "logs/" + str(date) + ".j%j.TACO.out",
        "logs/" + str(date) + ".j%j.TACO.err"
     threads: 1
     shell:
        "/apps/TACO/SRC/taco-0.6.3/build/scripts-2.7/taco_run.py --filter-min-length {params.min_lenght} --assemble-unstranded --ref-genome-fasta {input.genome} -o {output.assembled_read_taco} {input.gtf_list}"

rule Rfam_assembling:
     input:
        genome = "/scratch/devel/imoscardo/snakemake/genome/genome.fa",
        assembled_read_taco = "/scratch/devel/imoscardo/snakemake/assembly/taco/assembly.bed"
     output:
        fa_file = "/scratch/devel/imoscardo/snakemake/annotation/assembly/allseq.fa",
        annotation = "/scratch/devel/imoscardo/snakemake/annotation/assembly/assembly_annotation/cmsearch.out"
     log:
       "logs/" + str(date) + ".j%j.Rfam_assembling.out",
       "logs/" + str(date) + ".j%j.Rfam_assembling.err"
     threads: 4
     shell: 
        "module load bedtools;"
        "bedtools getfasta -fi {input.genome} -bed {input.assembled_read_taco} -fo {output.fa_file};" 
        "cd /scratch/devel/imoscardo/snakemake/annotation/assembly/assembly_annotation;"
        "~jgomez/bin/tRNAscan-SE-2.0.8/bin/cmsearch --noali -o {output.annotation} --tblout /scratch/devel/imoscardo/snakemake/annotation/assembly/assembly_annotation/cmsearch.tblout /scratch/devel/jgomez/RFAM_db/290721/Rfam.cm {output.fa_file}"