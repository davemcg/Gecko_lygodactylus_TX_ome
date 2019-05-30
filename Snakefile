'''
Snakefile for gekko transcfiptome
-currently only goes to building gtf part, need to add blast and building annotation


'''
import string

def readSampleFile(samplefile):
    # returns a dictionary of dictionaries where first dict key is sample id and second dict key are sample  properties
    res={}
    with open(samplefile) as file:
        for line in file:
            info=line.strip('\n').split('\t')
            res[info[0]]={'path':info[1].split(','),'type':info[2].strip('\n')}
    return(res)


configfile:'config.yaml'
sample_dict=readSampleFile(config['sampleFile'])# sampleID:dict{path,paired,metadata}
sample_names=list(sample_dict.keys())
STARindex='ref/STARindex'
refGenome= config['refGenome']
refGFF= config['refGFF']
stringtie_gtf='results/gekko_st.gtf'

rule all:
    input: 'notebooks/compare_txome_builds.html',expand('re_quant_files/{sample}/quant.sf', sample=sample_names)

'''
****PART 1****  get annotation run alignment
'''
rule downloadGencode:
    output:refGenome, refGFF
    shell:
        '''
        bash scripts/downloadGencode.sh {config[refGenome_url]} {config[refGFF_url]} {output}
        '''
rule run_trinity:
    input:[sample_dict[sample]['path'] for sample in sample_names]
    output: 'trinity_out_dir/Trinity.fasta'
    shell:
        '''
        module load trinity
        left=$(ls -d fastq_files/* | grep R1_ | tr '\n' ',')
        right=$(ls -d fastq_files/* | grep R2_ | tr '\n' ',')
        Trinity --seqType fq --CPU 32 --left $left --right $right --max_memory 1000G --trimmomatic
        '''
rule run_trans_decoder:
    input:'trinity_out_dir/Trinity.fasta'
    output:'trinity_out_dir/Trinity.fasta.transdecoder.pep'
    '''
    cd trinity_out_dir
    module load TransDecoder
    TransDecoder.LongOrfs -t {input}
    TransDecoder.Predict --single_best_only -t {input}
    '''

rule cluster_raw_tx:
    input: 'trinity_out_dir/Trinity.fasta', 'trinity_out_dir/Trinity.fasta.transdecoder.pep'
    output: tab='ref/pep_info.tsv', tmp='ref/tmp.fasta', clust='results/trinity_filtered.fasta'
    shell:
        '''
        python3 scripts/selectPCtx.py {input} {output.tab} {output.tmp}
        module load cd-hit
        cd-hit-est -i {output.tmp} -o {output.clust}  -c .95 -M 0 -T 0 -d 0
        '''

rule build_STARindex:
    input:'results/trinity_filtered.fasta'
    output: directory('ref/STARindex_trinity')
    shell:
        '''
        module load STAR
        mkdir -p {output}
        STAR --runThreadN 4 --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input}
        '''
rule align_star:
    input:fastqs=lambda wildcards: sample_dict[wildcards.id]['path'], index='ref/STARindex_trinity'
    params: outfolder='STARbams/{id}/'
    output:'STARbams/{id}/Aligned.out.bam', 'STARbams/{id}/Log.final.out'
    shell:
        '''
        module load STAR
        mkdir -p {params.outfolder}
        STAR  --outSAMstrandField intronMotif --outSAMtype BAM Unsorted  \
         --alignIntronMax 299999 --runThreadN 8 --genomeDir {input.index}  \
         --readFilesIn {input.fastqs} --readFilesCommand gunzip -c --outFileNamePrefix {params.outfolder}
        '''
rule sort_bams:
    input:'STARbams/{id}/Aligned.out.bam'
    output:'STARbams/{id}/Sorted.out.bam'
    shell:
        '''
        module load samtools
        samtools sort -@4 -O bam -o {output} {input}
        samtools index -b  {output}
        '''

rule calculate_cov:
    input:'STARbams/{id}/Sorted.out.bam'
    output: 'coverage_files/{id}.per-base.bed.gz'
    shell:
        '''
        module load mosdepth
        sample={wildcards.id}
        mosdepth coverage_files/$sample {input}
        '''
rule build_salmon_index:
    input: 'results/trinity_filtered.fasta'
    output: 'ref/salmonIndex_trinity'
    shell:
        '''
        module load salmon
        salmon index -t {input} -i {output} --type quasi --perfectHash -k 31
        '''
rule run_salmon:
    input: fastqs = lambda wildcards: sample_dict[wildcards.sample]['path'],
            index='ref/salmonIndex_trinity'
            #index= 'ref/salmonIndex_gekko'
    output: 'quant_files/{sample}/quant.sf'
    shell:
        '''
        module load salmon
        salmon quant -i {input.index} -l A --gcBias --seqBias -p 8 -1 {input.fastqs[0]} -2 {input.fastqs[1]} -o quant_files/{wildcards.sample}

        '''
rule aggregate_files_make_report:
    input:expand('coverage_files/{id}.per-base.bed.gz', id=sample_names), \
    expand('quant_files/{id}/quant.sf', id=sample_names)
    output:'notebooks/compare_txome_builds.html', 'results/filtered_tx_to_gene.tsv'
    shell:
        '''
        module load R
        Rscript ~/scripts/render_rmd.R notebooks/compare_txome_builds.Rmd
        '''
rule filter_fasta_for_high_conf_tx:
    input:infasta='results/trinity_filtered.fasta',  txtab='results/filtered_tx_to_gene.tsv'
    output:'results/annoted_tx_high_conf.fasta'
    shell:
        '''
        cut -f1 {input.txtab} > /tmp/tx_to_keep.txt
        python3 scripts/filterFasta.py {input.infasta} /tmp/tx_to_keep.txt {output}
        '''
rule rebuild_salmon_index:
    input:'results/annoted_tx_high_conf.fasta'
    output:directory('ref/salmonIndex_annoTx')
    shell:
        '''
        module load salmon
        salmon index -t {input} -i {output} --type quasi --perfectHash -k 31
        '''

rule requantify_salmon:
    input: fastqs = lambda wildcards: sample_dict[wildcards.sample]['path'],
            index='ref/salmonIndex_annoTx'
            #index= 'ref/salmonIndex_gekko'
    output: 're_quant_files/{sample}/quant.sf'
    shell:
        '''
        module load salmon
        salmon quant -i {input.index} -l A --gcBias --seqBias -p 8 -1 {input.fastqs[0]} -2 {input.fastqs[1]} -o quant_files/{wildcards.sample}

        '''
