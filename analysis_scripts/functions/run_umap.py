#!/home/user/miniconda3/bin/python3

import sys, getopt

from pathlib import Path
import umap
import pandas as pd
import re

def run_umap(fin, fout, nc = 2, nn=25, md=0.4, niter=1):
    df=pd.read_csv(fin, index_col=0)
    for i in range(niter):
        embedding2d = umap.UMAP(n_components=nc, metric='correlation', n_neighbors=nn, min_dist=md, verbose=True).fit_transform(df.values)
        rd = pd.DataFrame(embedding2d)
        rd.index = df.index
        if niter ==1:
            rd.to_csv(fout)
        else:
            rd.to_csv(Path(re.sub("csv", str(i)+".csv", str(fout))))
            
def run_umap_init(fin, finit, fout, nc = 2, nn=25, md=0.3, niter=1):
    df1=pd.read_csv(fin, index_col=0)
    df2=pd.read_csv(finit, index_col=0)
    for i in range(niter):
        embedding2d = umap.UMAP(n_components=nc, metric='correlation', n_neighbors=nn, min_dist=md, verbose=True, init = df2.values).fit_transform(df1.values)
        rd = pd.DataFrame(embedding2d)
        rd.index = df1.index
        if niter ==1:
            rd.to_csv(fout)
        else:
            rd.to_csv(fout+i)


def main(argv):
    fin = ''
    fout = ''
    finit = ''
    nc=2
    nn = 25
    md=0.4
    niter=1
    try:
        opts, args = getopt.getopt(argv, "hc:d:n:k:i:o:t:",["infile=", "outfile=", "init_file="])
    except getopt.GetoptError:
        print('run_umap.py -i <fin> -o <fout> -t <finit> -c <nc> -d <md> -k <nn> -n <niter>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == 'h':
            print('run_umap.py -i <fin> -o <fout> -t <finit> -c <nc> -d <md> -k <nn> -n <niter>')            
        if opt in ("-i", "--infile"):
            fin = arg
        if opt in ("-o", "--outfile"):
            fout = arg
        if opt in ("-t", "--init_outfile"):
            finit = arg
        if opt == "-c":
            nc = int(arg)
        if opt == "-k":
            nn = int(arg)
        if opt == "-d":
            md = float(arg)
        if opt in ("-n"):
            niter = int(arg)
    print("fin:",fin, "finit:", finit,"fout:", fout,"niter:", niter,"\n")
    if finit =='':
            run_umap(fin,fout, nc = nc, nn=nn, md=md, niter=niter)
    else:
        run_umap_init(fin,finit,fout, nc = nc, nn=nn, md=md, niter=niter)
        

if __name__ == "__main__":
   main(sys.argv[1:])

