import sys
import pandas as pd
with open(sys.argv[1]) as infasta:
    line=infasta.readline().strip()
    od=dict()
    length='length'
    val='tx'
    while line:
        if '>' in line:
            od[val]=length
            val=line[1:]
            length=0
        else:
            length+=len(line)
        line=infasta.readline().strip()
    df=pd.DataFrame.from_dict(od,orient='index')
    df.to_csv(sys.argv[2],header=False)
