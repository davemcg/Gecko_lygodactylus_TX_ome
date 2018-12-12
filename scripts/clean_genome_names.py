import sys
with open(sys.argv[1]) as infl, open(sys.argv[2],'w+') as ofl:
    for line in infl:
        if line[0]=='>':
            line='>'+line.split('|')[3]+'\n'
        ofl.write(line)
