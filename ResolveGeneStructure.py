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
    pass


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
    gene_name = args.gene 
    gff_file = args.gff
    snp_file = args.SNP
    sv_file = args.sv    #can be ignore
    cnv_file = args.cnv  #can be ignore
    inf_file = args.info

    ##gene_structure
