#fengcong@caas.cn
#2020-12-27 17:59:11
#usage:python __file__ -h

#function:  generate a json file for plot. 

#resovle gene structure according to SNP/SV/CNV matrix

import sys,os
import argparse
import re
import json
import subprocess
import copy

Rscript="/public/home/fengcong/anaconda2/envs/R/bin/Rscript"
py3="/public/home/fengcong/anaconda2/envs/py3/bin/python"

def get_gene_structure(chr_name,gene_name,gff_file):
    sys.stderr.write("resolve Gene structure ...\n")
    '''
    "gene_structure":{
        "gene":[[1,8200]],
        "gene_name":"TraesCS1A02G000100.1",
        "upstream":[[1,2000]],
        "5p_UTR":[[2001,2050]],
        "exon":[[2051,3400],[5100,5700],[7300,8000]],
        "3p_UTR":[[8001,8200]],
        "ori":"+"
    }
    '''
    gene_structure = {
        "chr":"",
        "gene":[],
        "gene_name":"",
        "upstream":[],
        "5p_UTR":[],
        "exon":[],
        "3p_UTR":[],
        "ori":""
    }
    gene_structure["chr"]=chr_name
    gene_structure["gene_name"]=gene_name

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
                gene_structure["ori"]=ls[6]
                if gene_structure["ori"] == "+":
                    up_start=int(ls[3])-2001
                    up_end=int(ls[3])-1
                    if up_start<0:
                        up_start=0
                    gene_structure["upstream"].append([up_start,up_end])
                else:
                    up_start=int(ls[3])+2001
                    up_end=int(ls[3])+1
                    gene_structure["upstream"].append([up_end,up_start])
                
            else:
                if line.find(gene_name) != -1:
                    if ls[2] =="mRNA":
                        pass
                    elif ls[2] == "five_prime_UTR":
                        gene_structure["5p_UTR"].append([int(ls[3]),int(ls[4])])
                    elif ls[2] == "exon":
                        gene_structure["exon"].append([int(ls[3]),int(ls[4])])
                    elif ls[2] == "three_prime_UTR":
                        gene_structure["3p_UTR"].append([int(ls[3]),int(ls[4])])
                    else:
                        pass
    
    return gene_structure
                
def deal_SNP(chr_name,snp_sample_order,gene_structure,snp_file):
    sys.stderr.write("resolve SNP ...\n")
    '''
    "SNP":
    [
        {
            "position":2500,
            "stat_in_each_sample":[
                "A",
                "A/T",
                "T",
                "T",
                "T",
                "T"
            ]
        }
    ]
    '''
    '''
    snp_inf={
        "chr":"chrUn",
        "position":2500,
        "ref":"A",
        "alt":"T",
        "ann":[[ann1],[ann2]]
        ##ann1: Annotation , Gene_Name ,Feature_ID ,HGVS.c ,HGVS.p
    }
    '''
    ##check vcf.gz file
    if snp_file.endswith(".gz"):
        if os.path.exists(snp_file+".csi"):
            pass
        else:
            sys.stderr.write("didnt find index file for the vcf.gz file\n")
            exit(-1)
    else:
        sys.stderr.write("vcf file must be bgziped and indexed using bcftools\n")
        exit(-1)
    
    ##calc the gene+upstream region
    need_region=[]
    if gene_structure["ori"] == "+":
        need_region=[gene_structure["upsteam"][0][0],gene_structure["gene"][0][1]]
    else:
        need_region=[gene_structure["gene"][0][0],gene_structure["upstream"][0][1]]

    ##get the snp in the region
    child = subprocess.Popen("bcftools view -r %s:%d-%d %s"%(chr_name,need_region[0],need_region[1],snp_file),shell=True,stdout=subprocess.PIPE)
    # output = child.communicate()

    #############################################template
    SNP=[]
    SNP_item={
        "position":0,
        "stat_in_each_sample":[
            
        ]
    }

    SNP_INF=[]
    SNP_INF_item={
        "chr":"",
        "position":0,
        "ref":"",
        "alt":"",
        "ann":[]
        ##ann1: allele,Annotation , Gene_Name ,Feature_ID ,HGVS.c ,HGVS.p
    }
    #############################################template
    ##deal the output

    for line in child.stdout.readlines():
        line = line.decode("utf-8")
        if line.startswith("##"):
            pass
        elif line.startswith("#"):
            snp_sample_order=line.strip().split()[9:]
        else:
            ls = line.strip().split()
            tmp_snp_item=copy.deepcopy(SNP_item)
            tmp_snp_inf_item=copy.deepcopy(SNP_INF_item)
            tmp_snp_item["position"] = int(ls[1])
            tmp_snp_inf_item["chr"]=chr_name
            tmp_snp_inf_item["position"] = int(ls[1])
            tmp_snp_inf_item["ref"] = ls[3]
            tmp_snp_inf_item["alt"] = ls[4]
            pattern = re.compile(r'.*ANN=(.*)')
            match = pattern.match(ls[7])
            ann = match.groups()[0]
            ann_list = ann.split(",")
            for a in ann_list:
                ann_ls = a.split("|")
                tmp_snp_inf_item["ann"].append([ann_ls[0],ann_ls[1],ann_ls[3],ann_ls[6],ann_ls[9],ann_ls[10]])
            for sample in ls[9:]:
                gt = sample.split(":")[0]
                if gt[0] == gt[-1]:
                    if gt[0] == "1":
                        tmp_snp_item["stat_in_each_sample"].append("%s/%s"%(ls[4],ls[4]))
                    elif gt[0] == "0":
                        tmp_snp_item["stat_in_each_sample"].append("%s/%s"%(ls[3],ls[3]))
                    else:
                        tmp_snp_item["stat_in_each_sample"].append("./.")
                else:
                    tmp_snp_item["stat_in_each_sample"].append("%s/%s"%(ls[3],ls[4]))
        
            SNP_INF.append(tmp_snp_inf_item)
            SNP.append(tmp_snp_item)


    return (SNP,SNP_INF)






if __name__ == "__main__":
    cmdparser = argparse.ArgumentParser(description="resolve gene structure")
    cmdparser.add_argument("-g", "--gene", dest="gene", type=str, required=True,
                           help="gene name.")
    cmdparser.add_argument("-c", "--chr", dest="chr", type=str, required=True,
                           help="chr name.")
    cmdparser.add_argument("-f","--gff", dest="gff",type=str, required=True,
                            help="gff file for gene annotation.")
    cmdparser.add_argument("-S","--SNP", dest="SNP",type=str, required=True,
                            help="SNP matrix")
    cmdparser.add_argument("-X", "--sv", dest="sv", type=str,
                           help="Structure variation matrix.")
    cmdparser.add_argument("-C", "--cnv", dest="cnv", type=str,
                           help="cnv matrix. ")
    # cmdparser.add_argument("-i", "--info", dest="info", type=str,required=True,
    #                        help="sample information.")
    cmdparser.add_argument("-o", "--outputprefix", dest="outputprefix", type=str,
                           help="outputprefix")


    args = cmdparser.parse_args()


    ##deal with the args
    gene_name = args.gene if args.gene.find(".") != -1 else args.gene+".1" 
    gff_file = args.gff
    snp_file = args.SNP
    sv_file = args.sv    #can be ignore
    cnv_file = args.cnv  #can be ignore
    # inf_file = args.info
    chr_name= args.chr
    output_file = args.outputprefix+".json"

    ##gene_structure
    gene_structure=get_gene_structure(chr_name,gene_name,gff_file)

    ##maybe they have diff order
    snp_sample_order=[]
    sv_sample_order=[]
    cnv_sample_order=[]

    ##deal SNP matrix
    SNP_and_INF=deal_SNP(chr_name,snp_sample_order,gene_structure,snp_file)






    ##test ouput 
    output_dict={
        "gene_structure":gene_structure,
        "variation":{
            "SNP":SNP_and_INF[0]
            },
        "info":{
            "SNP_INF":SNP_and_INF[1]
            }
    }
    js = json.dumps(output_dict, sort_keys=False, indent=4, separators=(',', ':'))
    ouf = open(output_file,"w")
    ouf.write(js)
    ouf.close()