'''
Snakefile for gekko transcfiptome
-currently only goes to building gtf part, need to add blast and building annotation


'''
import string
from itertools import product


def readSampleFile(samplefile):
    # returns a dictionary of dictionaries where first dict key is sample id and second dict key are sample  properties
    res={}
    with open(samplefile) as file:
        for line in file:
            info=line.strip('\n').split('\t')
            res[info[0]]={'path':info[1].split(','),'type':info[2].strip('\n')}
    return(res)

def make_chunk_input(stem):
    alphbet_pfx=list(product(list(string.ascii_lowercase),list(string.ascii_lowercase)))
    paths=[stem + ''.join(i) for i in alphbet_pfx]
    paths.sort()
    return(paths[0:100])

configfile:'config.yaml'
sample_dict=readSampleFile(config['sampleFile'])# sampleID:dict{path,paired,metadata}
#subtissues_SE=["RPE_Stem.Cell.Line","RPE_Cell.Line","Retina_Adult.Tissue","RPE_Fetal.Tissue","Cornea_Adult.Tissue","Cornea_Fetal.Tissue","Cornea_Cell.Line","Retina_Stem.Cell.Line",'body']
sample_names=list(sample_dict.keys())
STARindex='ref/STARindex'
refGenome= config['refGenome']
refGFF= config['refGFF']


rule all:
    input: 'ref/gekko.combined.gtf', expand('quant_files/{sampleID}/quant.sf',sampleID=sample_names), \
    'ref/all_blastn.tsv','ref/all_blastp.tsv'
    #,'smoothed_filtered_tpms.csv'
'''
****PART 1****  get annotation run alignment
'''
rule downloadGencode:
    output:refGenome, refGFF
    shell:
        '''
        wget -O - {config[refGenome_url]} | gunzip -c > {refGenome}.tmp
        wget -O - {config[refGFF_url]} | gunzip -c > {refGFF}
        python3 scripts/clean_genome_names.py {refGenome}.tmp {refGenome}
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
        stringtie {input[0]} -o {output[0]} -p 8 -G {refGFF}
        '''
'''
---making gtfs
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
    output: 'ref/gekko_st.gtf'
    shell:
        '''
        module load stringtie
        stringtie --merge -G {refGFF} -o {output} {input}
        '''
rule compare_st_gffcompare:
    input: 'ref/gekko_st.gtf','ref/gekko.combined.gtf'
    output:'ref/st_gffc_gekko.stats'
    shell:
        '''
        module load gffcompare
        gffcompare -r {refGFF} -o ref/ref/st_gffc_gekko {input}
        '''
rule extract_tx_seqs:
    input: 'ref/gekko_st.gtf',refGenome
    output: 'ref/gekko_tx.fa'
    shell:
        '''
        ./gffread/gffread -w {output[0]} -g {input[1]} {input[0]}
        '''

rule run_trans_decoder:
    input:'ref/gekko_tx.fa'
    output:'trdec_dir/gekko_tx.fa.transdecoder_dir/best_orfs.pep'
    shell:
        '''
        module load TransDecoder
        cd trdec_dir
        TransDecoder.LongOrfs -t ../{input}
        TransDecoder.Predict --single_best_only -t ../{input}
        mv gekko_tx.fa.transdecoder.pep   gekko_tx.fa.transdecoder_dir/best_orfs.pep
        mv gekko_tx.fa.transdecoder.*   gekko_tx.fa.transdecoder_dir/ || true
        '''
rule run_blastn:
    input:'ref/gekko_tx.fa'
    output:'ref/all_blastn.tsv'
    shell:
        '''
        module load blast
        blastn -query {input[0]} -db /fdb/blastdb/nt  -max_target_seqs 250 -max_hsps 3 -outfmt 6 -num_threads 8 > {output[0]}
        '''

rule run_blastp:
    input: 'trdec_dir/gekko_tx.fa.transdecoder_dir/best_orfs.pep'
    output:'ref/all_blastp.tsv'
    shell:
        '''
        module load blast
        blastp -query {input[0]} -db /fdb/blastdb/swissprot  -max_target_seqs 250 -max_hsps 3 -outfmt 6 -num_threads 8 > {output[0]}
        '''



rule build_salmon_index:
        input: 'ref/gekko_tx.fa'
        output: 'ref/salmonIndex_gekko'
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
