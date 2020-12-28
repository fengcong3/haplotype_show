#fengcong@caas.cn
#2020-12-27 17:59:11
#usage:python __file__ -h

#function:  generate a json file for plot. 

#resovle gene structure according to SNP/SV/CNV matrix

import sys,os
import argparse
import re
import json

Rscript="/public/home/fengcong/anaconda2/envs/R/bin/Rscript"
py3="/public/home/fengcong/anaconda2/envs/py3/bin/python"

def gene_structure(gene_name,gff_file):
    '''
    "gene_structure":{
        "gene":[[1,8200]],
        "gene_name":"TraesCS1A02G000100.1",
        "upstream":[[1,2000]],
        "5p_UTR":[[2001,2050]],
        "exon":[[2051,3400],[5100,5700],[7300,8000]],
        "intron":[[3401,5099],[5701,7299]],
        "3p_UTR":[[8001,8200]]
    }
    '''
    gene_structure = {
        "gene":[],
        "gene_name":"",
        "upstream":[],
        "5p_UTR":[],
        "exon":[],
        "intron":[],
        "3p_UTR":[]
    }
    #get gene information
    info=[]
    inf = open(gff_file,"r")
    for line in  inf.readlines():
        if line.find(gene_name.split(".")[0]) != -1:
            info.append(line.strip())
    inf.close()

    #deal with the information 
    if len(info) == 0:
        sys.stderr.write("can't find this gene in gff file\n")
        exit(-1)
    else:
        for line in info:
            ls = line.split()
            if ls[2] == "gene":
                gene_structure["gene"].append([int(ls[3]),int(ls[4])])
            else:
                if ls[2] =="mRNA":
                    pass
                ##unfinished 



if __name__ == "__main__":
    cmdparser = argparse.ArgumentParser(description="haploblock pipline(plot and calc common hap)")
    cmdparser.add_argument("-g", "--gene", dest="gene", type=str, required=True,
                           help="gene name.")
    cmdparser.add_argument("-f","--gff", dest="gff",type=str, required=True,
                            help="gff file for gene annotation.")
    cmdparser.add_argument("-S","--SNP", dest="SNP",type=str, required=True,
                            help="SNP matrix")
    cmdparser.add_argument("-X", "--sv", dest="sv", type=str,
                           help="Structure variation matrix.")
    cmdparser.add_argument("-C", "--cnv", dest="cnv", type=str,
                           help="cnv matrix. ")
    cmdparser.add_argument("-i", "--info", dest="info", type=str,required=True,
                           help="sample information.")
    cmdparser.add_argument("-o", "--outputprefix", dest="outputprefix", type=str,
                           help="outputprefix")


    args = cmdparser.parse_args()


    ##deal with the args
    gene_name = args.gene if args.gene.find(".") != -1 else args.gene+".1" 
    gff_file = args.gff
    snp_file = args.SNP
    sv_file = args.sv    #can be ignore
    cnv_file = args.cnv  #can be ignore
    inf_file = args.info

    ##gene_structure
    

