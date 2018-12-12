'''
Snakefile for gekko transcfiptome
-currently only goes to building gtf part, need to add blast and building annotation
-need to get the refernce still

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
    input: 'ref/gekko_combined.gtf'#, expand('quant_files/{sampleID}/quant.sf',sampleID=sample_names)
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
        wget {config[refGenome_url]} | gunzip -c > {refGenome}
        wget {config[refGFF_url]} | gunzip -c > {refGFF}
        '''

rule build_STARindex:
    input: refGenome, refGFF
    output: STARindex
    shell:
        '''
        module load STAR
        mkdir -p ref/STARindex
        STAR --runThreadN 16 --runMode genomeGenerate --genomeDir {output[0]} --genomeFastaFiles {refGenome} \
         --sjdbGTFtagExonParentTranscript transcript_id --sjdbGTFfile {refGFF} --sjdbOverhang 100
        '''

rule run_STAR_alignment:
    input: lambda wildcards: sample_dict[wildcards.sample]['path'], STARindex
    output:temp('STARbams/{sample}/raw.Aligned.out.bam')
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
    output:'STARbams/{id}/Aligned.out.bam'
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
        stringtie {input[0]} -o {output[0]} -p 8 -G ref/gencodeAno_bsc.gtf
        '''

rule merge_gtfs:
    #this could be split into multiple rules, but its a lot easier to string it together
    input: expand('st_out/{sample}.gtf', sample=sample_names )
    output: 'ref/gekko_combined.gtf', 'ref/gekko.stats'
    shell:
        '''
        module load gffcompare
        gffcompare -r {refGFF} -o ref/gekko {input}
        '''
