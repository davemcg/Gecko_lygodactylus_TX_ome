import sys
import subprocess as sp
args=sys.argv
comp_fasta=args[1]
trim_fasta=args[2]
trim_table=args[3]
filt_fasta=args[4]
sp.check_output("./scripts/trinity_defline_oneliner.sh {} {} ".format(trim_fasta,trim_table) , shell=True)
with open(comp_fasta) as infasta, open(trim_table) as bad_tx, open(filt_fasta,'w+') as outfasta:
    names=set()
    for line in bad_tx:
        names.add('>'+line.strip().split(' ')[0])
    oldline=infasta.readline().strip().split(' ')[0]
    while oldline:
        if oldline in names and '>' in oldline:
            write=True
        elif oldline not in names and '>' in oldline:
            write=False
        if write:
            outfasta.write(oldline+'\n')
        oldline=infasta.readline().strip().split(' ')[0]
