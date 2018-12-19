'''
Snakefile for gekko transcfiptome
-currently only goes to building gtf part, need to add blast and building annotation


'''

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
#subtissues_SE=["RPE_Stem.Cell.Line","RPE_Cell.Line","Retina_Adult.Tissue","RPE_Fetal.Tissue","Cornea_Adult.Tissue","Cornea_Fetal.Tissue","Cornea_Cell.Line","Retina_Stem.Cell.Line",'body']
sample_names=list(sample_dict.keys())
STARindex='ref/STARindex'
refGenome= config['refGenome']
refGFF= config['refGFF']


rule all:
    input: 'ref/gekko.combined.gtf', expand('quant_files/{sampleID}/quant.sf',sampleID=sample_names)
    #,'smoothed_filtered_tpms.csv'
'''
****PART 1**** download files
-still need to add missing fastq files
-gffread needs indexed fasta
-need to add versioning of tools to yaml
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
        ./gffread/gffread -w {output[0]} -g {refGenome} {input[0]}
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
