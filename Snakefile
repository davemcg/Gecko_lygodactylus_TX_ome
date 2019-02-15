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
    input: 'ref/st_gffc_gekko.stats', expand('quant_files/{sampleID}/quant.sf',sampleID=sample_names), \
    'results/all_blastn.tsv', 'results/all_blastp.tsv', expand('STARbams_realigned/{samples}/Aligned.out.bam',samples=sample_names)

'''
****PART 1****  get annotation run alignment
'''
rule downloadGencode:
    output:refGenome, refGFF
    shell:
        '''
        bash scripts/downloadGencode.sh {config[refGenome_url]} {config[refGFF_url]} {output}
        '''

rule build_STARindex:
    input: refGenome, refGFF
    output: STARindex
    shell:
        '''
        module load STAR
        mkdir -p ref/STARindex
        STAR --runThreadN 4 --runMode genomeGenerate --genomeDir {output[0]} --genomeFastaFiles {refGenome} --genomeChrBinNbits 9 \
         --sjdbGTFtagExonParentTranscript transcript_id --sjdbGTFfile {refGFF} --sjdbOverhang 100
        '''

rule run_STAR_alignment:
    input: lambda wildcards: sample_dict[wildcards.sample]['path'], STARindex
    output:temp('STARbams/{sample}/raw.Aligned.out.bam'),'STARbams/{sample}/raw.Log.final.out'
    shell:
        '''
        module load STAR
        mkdir -p STARbams/{wildcards.sample}
        STAR --runThreadN 8 --genomeDir {STARindex} --outSAMstrandField intronMotif \
            --readFilesCommand gunzip -c --outFileNamePrefix STARbams/{wildcards.sample}/raw. \
            --readFilesIn {input[0]} {input[1]} --outSAMtype BAM Unsorted
        '''

rule sort_bams:
    input:'STARbams/{id}/raw.Aligned.out.bam'
    output:'STARbams/{id}/Aligned.out.bam',
    shell:
        '''
        module load samtools
        samtools sort -o {output[0]} --threads 7 {input[0]}
        '''

rule run_stringtie:
    input: 'STARbams/{sample}/Aligned.out.bam'
    output:'st_out/{sample}.gtf'
    shell:
        '''
        module load stringtie
        stringtie {input[0]} -o {output[0]} -p 8

        '''
'''
***PART 1 - Build Transcriptome***
-gffcompare give good stats but does no filtering, bit st --merge auto filters (using transcript w min 1 tpm in 1 sample)
-remember to escape braces for non variable use, like awk
'''
rule merge_gtfs_gffcomp:
    #this could be split into multiple rules, but its a lot easier to string it together
    input: expand('st_out/{sample}.gtf', sample=sample_names )
    output: 'ref/gekko.combined.gtf', 'ref/gekko.stats'
    shell:
        '''
        module load gffcompare
        gffcompare -r {refGFF} -o ref/gekko {input}
        '''
rule merge_gtf_st:
    input: expand('st_out/{sample}.gtf', sample=sample_names ), 'ref/gekko.combined.gtf'
    output: 'results/gekko_st.gtf','results/gekko_st_TXidx.gtf'
    shell:
        '''
        module load stringtie
        stringtie --merge -G {refGFF} -o {output[0]} {input}
        module load R
        Rscript scripts/swapChromforTx.R {output}
        '''
rule compare_st_gffcompare:
    input: 'results/gekko_st.gtf','ref/gekko.combined.gtf'
    output:'ref/st_gffc_gekko.stats'
    shell:
        '''
        module load gffcompare
        gffcompare -r {refGFF} -o ref/st_gffc_gekko {input}
        '''
rule extract_tx_seqs:
    input: 'results/gekko_st.gtf',refGenome
    output: 'results/gekko_tx.fa'
    shell:
        '''
        ./gffread/gffread -w {output[0]} -g {input[1]} {input[0]}
        '''
rule run_trans_decoder:
    input:'results/gekko_tx.fa'
    output:'results/best_orfs.pep'
    shell:
        '''
        module load TransDecoder
        mkdir -p trdec_dir
        cd trdec_dir
        TransDecoder.LongOrfs -t ../{input}
        TransDecoder.Predict --single_best_only -t ../{input}
        mv gekko_tx.fa.transdecoder.pep best_orfs.pep
        mv gekko_tx.fa.transdecoder.*   ../results/
        '''

'''
***PART 2 Evaluate Transcriptome Accuraccy***
-align fasta to transcriptome and look at coverage of built transcripts
-blast transcripts aginst known nt, prot

'''


rule build_STAR_transcriptomic_index:
    input: fasta='results/gekko_tx.fa',
        gtf='results/gekko_st_TXidx.gtf'
    output:'ref/STAR_TXidx'
    shell:
        '''
        module load STAR
        mkdir -p {output}
        STAR --runThreadN 12 --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input.fasta} --genomeChrBinNbits 9 \
        --sjdbGTFfile {input.gtf} --sjdbOverhang 100
        '''

rule align_to_txome_STAR:
    input: fastqs=lambda wildcards: sample_dict[wildcards.sample]['path'], index='ref/STAR_TXidx'
    output:'STARbams_realigned/{sample}/Aligned.out.bam','STARbams_realigned/{sample}/Log.final.out'
    shell:
        '''
        module load STAR
        out_dir=STARbams_realigned/{wildcards.sample}/
        mkdir -p $out_dir
        STAR --runThreadN 8 --genomeDir {input.index} --outSAMstrandField intronMotif \
            --readFilesCommand gunzip -c --outFileNamePrefix $out_dir \
            --readFilesIn {input.fastqs} --outSAMtype BAM Unsorted
        '''

rule run_blastn:
    input:'results/gekko_tx.fa'
    output:'results/all_blastn.tsv'
    shell:
        '''
        module load blast
        blastn -query {input[0]} -db /fdb/blastdb/nt  -max_target_seqs 250 -max_hsps 3 -outfmt 6 -num_threads 8 > {output[0]}
        '''
rule getFastaLengths:
    input:'results/best_orfs.pep'
    output:'results/orf_lengths.csv'
    shell:
        '''
        python3 scripts/getFastaLengths.py {input} {output}
        '''
rule run_blastp:
    input: 'results/best_orfs.pep', 'results/orf_lengths.csv'
    output:'results/all_blastp.tsv'
    shell:
        '''
        module load blast
        blastp -query {input[0]} -db /fdb/blastdb/swissprot  -max_target_seqs 250 -max_hsps 3 -outfmt 6 -num_threads 8 > {output[0]}
        '''

rule build_salmon_index:
    input: 'results/gekko_tx.fa'
    output: 'results/salmonIndex_gekko'
    shell:
        '''
        module load salmon
        salmon index -t {input} code -i {output} --type quasi --perfectHash -k 31
        '''

rule run_salmon:
    input: fastqs = lambda wildcards: sample_dict[wildcards.sample]['path'],
            index= 'ref/salmonIndex_gekko'
    output: 'quant_files/{sample}/quant.sf'
    shell:
        '''
        module load salmon
        salmon quant -i {input.index} -l A --gcBias --seqBias -p 8 -1 {input.fastqs[0]} -2 {input.fastqs[1]} -o quant_files/{wildcards.sample}

        '''
