#fengcong@caas.cn
#2020-12-27 17:59:11
#usage:python __file__ -h

#function:  generate a json file for plot. 
#require : bcftools , bedtools

#resovle gene structure according to SNP/SV/CNV matrix

import sys,os
import argparse
import re
import json
import subprocess
import copy
import gzip

Rscript="/public/home/fengcong/anaconda2/envs/R/bin/Rscript"
py3="/public/home/fengcong/anaconda2/envs/py3/bin/python"
bcftools="/public/agis/chengshifeng_group/fengcong/WGRS/software/bcftools1.9/bin/bcftools"
bedtools="/public/agis/chengshifeng_group/fengcong/IPK_assembly/software/bedtools2/bin/bedtools"

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
                    up_start=int(ls[4])+2001
                    up_end=int(ls[4])+1
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
    child = subprocess.Popen("%s view -r %s:%d-%d %s"%(bcftools,chr_name,need_region[0],need_region[1],snp_file),shell=True,stdout=subprocess.PIPE)
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
                break
            for sample in ls[9:]:
                gt = sample.split(":")[0]
                if gt[0] == gt[-1]:
                    if gt[0] == "1":
                        tmp_snp_item["stat_in_each_sample"].append("%s"%(ls[4]))
                    elif gt[0] == "0":
                        tmp_snp_item["stat_in_each_sample"].append("%s"%(ls[3]))
                    else:
                        tmp_snp_item["stat_in_each_sample"].append("./.")
                else:
                    tmp_snp_item["stat_in_each_sample"].append("%s/%s"%(ls[3],ls[4]))
        
            SNP_INF.append(tmp_snp_inf_item)
            SNP.append(tmp_snp_item)

    
    return (SNP,SNP_INF,snp_sample_order)


def deal_indel(chr_name,indel_sample_order,gene_structure,indel_file):
    sys.stderr.write("resolve Indel ...\n")
    '''{
            "type":"DEL",
            "position":6500,
            "length":700,
            "stat_in_each_sample":[
                0,
                0,
                0,
                0,
                0,
                1
            ]
        }
    '''
    ##check vcf.gz file
    if indel_file.endswith(".gz"):
        if os.path.exists(indel_file+".csi"):
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
    child = subprocess.Popen("%s view -r %s:%d-%d %s"%(bcftools,chr_name,need_region[0],need_region[1],indel_file),shell=True,stdout=subprocess.PIPE)
    #########################indel template
    INDEL=[]
    INDEL_item={
        "type":"",
        "position":0,
        "length":0,
        "stat_in_each_sample":[
            
        ]
    }

    INDEL_INF=[]
    INDEL_INF_item={
        "chr":"",
        "position":0,
        "length":0,
        "ref":"",
        "alt":"",
        "ann":[]
        ##ann1: allele,Annotation , Gene_Name ,Feature_ID ,HGVS.c ,HGVS.p
    }
    ###########################INDEL TEMPLATE
    for line in child.stdout.readlines():
        line = line.decode("utf-8")
        if line.startswith("##"):
            pass
        elif line.startswith("#"):
            indel_sample_order=line.strip().split()[9:]
        else:
            ls = line.strip().split()
            tmp_indel_item=copy.deepcopy(INDEL_item)
            tmp_indel_inf_item=copy.deepcopy(INDEL_INF_item)
            tmp_indel_item["position"] = int(ls[1])
            tmp_indel_inf_item["chr"]=chr_name
            tmp_indel_inf_item["position"] = int(ls[1])
            tmp_indel_inf_item["ref"] = ls[3]
            tmp_indel_inf_item["alt"] = ls[4]
            if len(ls[3]) > len(ls[4]) :
                tmp_indel_item["type"] = "DEL"
            else:
                tmp_indel_item["type"] = "INS"
            
            tmp_indel_inf_item["length"] = abs(len(ls[3]) - len(ls[4]))
            tmp_indel_item["length"] = abs(len(ls[3]) - len(ls[4]))
            
            pattern = re.compile(r'.*ANN=(.*)')
            match = pattern.match(ls[7])
            ann = match.groups()[0]
            ann_list = ann.split(",")
            for a in ann_list:
                ann_ls = a.split("|")
                tmp_indel_inf_item["ann"].append([ann_ls[0],ann_ls[1],ann_ls[3],ann_ls[6],ann_ls[9],ann_ls[10]])
                break
            for sample in ls[9:]:
                gt = sample.split(":")[0]
                if gt[0] == gt[-1]:
                    if gt[0] == "1":
                        tmp_indel_item["stat_in_each_sample"].append("1")
                    elif gt[0] == "0":
                        tmp_indel_item["stat_in_each_sample"].append("0")
                    else:
                        tmp_indel_item["stat_in_each_sample"].append("./.")
                else:
                    tmp_indel_item["stat_in_each_sample"].append("0/1")

            INDEL_INF.append(tmp_indel_inf_item)
            INDEL.append(tmp_indel_item)
        
    return (INDEL,INDEL_INF,indel_sample_order)
        
def deal_sv(chr_name,sv_sample_order,gene_structure,sv_file,part_len):
    sys.stderr.write("resolve SV ...\n")
    '''{
            "type":"DEL", #DUP INV
            "position":6500,
            "length":700,
            "stat_in_each_sample":[
                0,
                0,
                0,
                0,
                0,
                1
            ]
        }
    '''
    ##calc the gene+upstream region
    need_region=[]
    if gene_structure["ori"] == "+":
        need_region=[gene_structure["upsteam"][0][0],gene_structure["gene"][0][1]]
    else:
        need_region=[gene_structure["gene"][0][0],gene_structure["upstream"][0][1]]

    #chr pos convert to part pos
    need_chr=chr_name
    if chr_name!="chrUn":
        part_len_d={}
        inf = open(part_len,"r")
        for line in inf.readlines():
            part_len_d[line.strip().split()[0]] = int(line.strip().split()[1])
        inf.close()

        ## didnt check the region whether cross the part chr ##
        ## this can be a bug ##
        ## if u find that this isnt right, u should fix ths bug ##
        if need_region[0]  > part_len_d[chr_name+"_part1"]:
            need_chr = chr_name+"_part2"
            need_region[0] = need_region[0]-part_len_d[chr_name+"_part1"]
            need_region[1] = need_region[1]-part_len_d[chr_name+"_part1"]
        else:
            need_chr = chr_name+"_part1"
    else:
        pass
    

    ##check vcf.gz file
    if sv_file.endswith(".gz"):
        if os.path.exists(sv_file+".csi"):
            pass
        else:
            sys.stderr.write("didnt find index file for the vcf.gz file\n")
            exit(-1)
    else:
        sys.stderr.write("vcf file must be bgziped and indexed using bcftools\n")
        exit(-1)

    ##get the snp in the region
    
    child = subprocess.Popen("%s view -r %s:%d-%d %s"%(bcftools,need_chr,need_region[0],need_region[1],sv_file),shell=True,stdout=subprocess.PIPE)
    #########################SV template
    SV=[]
    SV_item={
        "type":"",
        "position":0,
        "length":0,
        "stat_in_each_sample":[
            
        ]
    }

    SV_INF=[]
    SV_INF_item={
        "chr":"",
        "position":0,
        "length":0,
        "ref":"",
        "alt":"",
        "ann":[]
        ##ann1: allele,Annotation , Gene_Name ,Feature_ID ,HGVS.c ,HGVS.p
    }
    ###########################SV TEMPLATE
    for line in child.stdout.readlines():
        line = line.decode("utf-8")
        if line.startswith("##"):
            pass
        elif line.startswith("#"):
            sv_sample_order=line.strip().split()[9:]
        else:
            
            # print(line)
            ls = line.strip().split()
            if int(ls[1]) >= need_region[0] and int(ls[1]) < need_region[1]:
                pass
            else:
                continue
            tmp_sv_item=copy.deepcopy(SV_item)
            tmp_sv_inf_item=copy.deepcopy(SV_INF_item)
            tmp_sv_item["position"] = int(ls[1])
            tmp_sv_inf_item["chr"]=chr_name
            tmp_sv_inf_item["position"] = int(ls[1]) if ls[0].endswith("_part1") else int(ls[1])+part_len_d[ls[0]]
            tmp_sv_inf_item["ref"] = ls[3]
            tmp_sv_inf_item["alt"] = ls[4]
            pattern = re.compile(r".*SVTYPE=(.{,10});.*END=(\d+);.*")
            match=pattern.match(ls[7])
            if not match:
                line = inf.readline()
                continue
            svtype=match.groups()[0]
            ed=int(match.groups()[1])

            tmp_sv_inf_item["type"] = svtype
            tmp_sv_item["type"] = svtype
            
            tmp_sv_inf_item["length"] = ed - int(ls[1])
            tmp_sv_item["length"] = ed - int(ls[1])
            
            
            for sample in ls[9:]:
                gt = sample.split(":")[0]
                if gt[0] == gt[-1]:
                    if gt[0] == "1":
                        tmp_sv_item["stat_in_each_sample"].append("1")
                    elif gt[0] == "0":
                        tmp_sv_item["stat_in_each_sample"].append("0")
                    else:
                        tmp_sv_item["stat_in_each_sample"].append("./.")
                else:
                    tmp_sv_item["stat_in_each_sample"].append("0/1")

            SV_INF.append(tmp_sv_inf_item)
            SV.append(tmp_sv_item)
        
    return (SV,SV_INF,sv_sample_order)

def filter_and_sort_sv(SV_and_INF,snp_sample_order,sv_sample_order):
    new_order_index=[]
    for i in snp_sample_order :
        new_order_index.append(sv_sample_order.index(i))
    
    ##modify each SV item
    for item in SV_and_INF[0]:
        item["stat_in_each_sample"] = [item["stat_in_each_sample"][n] for n in new_order_index]


def deal_cnv(gene_name,cnv_sample_order,cnv_file):
    sys.stderr.write("resolve cnv ...\n")
    '''
    "CN":[1,3,2,2,2,2]
    '''
    gene_name=gene_name.split(".")[0]
    CN=[] 
    inf = gzip.open(cnv_file,"rt") if cnv_file.endswith(".gz") else open(cnv_file,"r")

    line = inf.readline()
    while line.startswith("##"):
        line = inf.readline()
    
    cnv_sample_order = line.strip().split()[9:]
    line = inf.readline()
    while line:
        ls = line.strip().split()
        if ls[2] != gene_name:
            pass
        else:
            for gt in ls[9:]:
                if gt[0] == 0 :
                    CN.append("1")
                else:
                    CN.append("N") # this is the first version,dont have indeed copy number

        line = inf.readline()
    inf.close()

    if len(CN) ==0:
        CN = ["1" for n in range(len(cnv_sample_order))]

    return (CN,cnv_sample_order)

def filter_and_sort_cnv(CN,snp_sample_order,cnv_sample_order):
    new_order_index=[]
    for i in snp_sample_order :
        new_order_index.append(cnv_sample_order.index(i))
    
    CN = [CN[0][n] for  n in new_order_index]
    return CN

def hap_clustering(sample_order,SNP_item_list,INDEL_item_list,SV_item_list,CN_list):
    sys.stderr.write("haplotype clustering ...\n")
    '''
    "SNP":
    [
        {
            "position":2500,
            "stat_in_each_sample":[
                "A",                "A/T",                "T",                "T",                "T",                "T"
            ]
        }
    ]
    '''
    '''{
            "type":"DEL",
            "position":6500,
            "length":700,
            "stat_in_each_sample":[
                0,                0,                0,                0,                0,                1
            ]
        }
    '''
    '''{
            "type":"DEL", #DUP INV
            "position":6500,
            "length":700,
            "stat_in_each_sample":[
                0,                0,                0,                0,                0,                1
            ]
        }
    '''
    '''
    "CN":[1,3,2,2,2,2]
    '''
    #init sample_variants dict
    sample_variants={} #{sample1:[],sample2:[]}
    for sample in sample_order:
        sample_variants[sample] = []
    
    #append SNP/indel/sv/cnv for each sample
    for index,sample in enumerate(sample_order):
        for snp_item in SNP_item_list :
            sample_variants[sample].append(snp_item["stat_in_each_sample"][index])
        for indel_item in INDEL_item_list :
            sample_variants[sample].append(indel_item["stat_in_each_sample"][index])
        for sv_item in SV_item_list:
            sample_variants[sample].append(sv_item["stat_in_each_sample"][index])
        for cn in CN_list:
            sample_variants[sample].append(CN_list[index])

    #this list will be returned
    sample_order_and_hap_cluster_inf=[] #[(sample1,0),(sample2,1),...] sample_name and cluster index

    ##cluster information
    cluster_inf=[]  #[[variants list1],[]]

    #now, clustering 
    for sample in sample_order:
        if len(cluster_inf) == 0: #there is no hap cluster exist
            cluster_inf.append(sample_variants[sample])
            sample_order_and_hap_cluster_inf.append((sample,0))
        else:
            already_confirm_cluster=False
            for index,hap in enumerate(cluster_inf):
                equal=True #this sample == this hap
                for variants_index in range(len(sample_variants[sample])):
                    if sample_variants[sample][variants_index] == cluster_inf[index][variants_index]:
                        pass
                    else:
                        equal=False
                        break
                #judge they equal or not 
                if equal:
                    sample_order_and_hap_cluster_inf.append((sample,index))
                    already_confirm_cluster=True
            #judge it has been clustering
            if not already_confirm_cluster:
                cluster_inf.append(sample_variants[sample])
                sample_order_and_hap_cluster_inf.append((sample,len(cluster_inf)-1))
    
    #clustering complete ,now,return
    sys.stderr.write("cluster number: %d\n"%(len(cluster_inf)))
    return (sample_order_and_hap_cluster_inf,cluster_inf)

def rep_idx_str(str,idx,char):
    if 0<= idx < len(str)-1:
        return str[:idx] + char + str[idx+1:]
    elif idx ==  len(str)-1:
        return str[:idx] + char
    else:
        sys.stderr.write("SNP rep error\n")
        return str

def resolve_hap_seq(INDEL_inf_item_list,SV_inf_item_list,ref_file,chr_name,gene_structure,cluster_inf,SNP_item_list,INDEL_item_list,SV_item_list,CN_list,outputprefix):
    sys.stderr.write("haplotype sequence resolving ...\n")
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
    #generate gene origin
    fasta_origin=[]
    if gene_structure["ori"] == "+":
        fasta_origin=[gene_structure["upstream"][0][0],gene_structure["gene"][0][1]]
    else:
        fasta_origin=[gene_structure["gene"][0][0],gene_structure["upstream"][0][1]]

    #output to bed file
    ouf = open("%s.bed"%(outputprefix),"w")
    ouf.write("%s\t%d\t%d\n"%(chr_name,fasta_origin[0]-1,fasta_origin[1]))
    ouf.close()


    #get the reference sequence 
    child = subprocess.Popen("%s getfasta -fi %s -bed %s"%(bedtools,ref_file,outputprefix+".bed"),shell=True,stdout=subprocess.PIPE)
    ref_seq=""
    for line in child.stdout.readlines():
        line = line.decode("utf-8")
        if not line.startswith(">"):
            ref_seq = line.strip()

    #
    ouf = open("%s.hap.fasta"%(outputprefix),"w")
    for hap_index,sample_variants_list in enumerate(cluster_inf):
        # for ,variant in enumerate(sample_variants_list):
        hap_seq=[n for n in ref_seq]
        # print(hap_seq)
        variant_index = 0
        for snp_item in SNP_item_list :
            rep_char=sample_variants_list[variant_index][-1]  ##make three character to one char
            # print(snp_item["position"],fasta_origin[0])
            hap_seq[snp_item["position"]-fasta_origin[0]] = rep_char ## didnt check the position
            variant_index+=1 ## move forward

        for indel_index,indel_item in enumerate(INDEL_item_list) :
            if sample_variants_list[variant_index]==1:
                if indel_item["type"] == "INS":
                    hap_seq[indel_item["position"]-fasta_origin[0]] = INDEL_inf_item_list[indel_index]["alt"]
                elif indel_item["type"] == "DEL":
                    #make the position's later to null string
                    for i in range(1,indel_item["length"]):
                        hap_seq[indel_item["position"]-fasta_origin[0]+i] = ''
                else:
                    sys.stderr.write("\tresolving indel error\n")                
            variant_index+=1 ## move forward

        for sv_index,sv_item in enumerate(SV_item_list):
            if sample_variants_list[variant_index]==1:
                if sv_item["type"] == "INS":
                    hap_seq[sv_item["position"]-fasta_origin[0]] = SV_inf_item_list[sv_index]["alt"]
                elif sv_item["type"] == "DEL":
                    #make the position's later to null string
                    for i in range(1,sv_item["length"]):
                        hap_seq[sv_item["position"]-fasta_origin[0]+i] = ''
                elif sv_item["type"] == "DUP":
                    ## deal it or not ?
                    hap_seq[sv_item["position"]-fasta_origin[0]+sv_item["length"]-1] = hap_seq[sv_item["position"]-fasta_origin[0]:sv_item["position"]-fasta_origin[0]+sv_item["length"]]
                else:
                    sys.stderr.write("\tresolving SV error\n")                
            variant_index+=1 ## move forward
        ouf.write(">%s_hap%d\n"%(gene_structure["gene_name"].split(".")[0],hap_index+1))
        ouf.write("".join(hap_seq)+"\n")

    ouf.close()
            
    



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
    cmdparser.add_argument("-I","--indel", dest="indel",type=str, 
                            help="indel matrix")
    cmdparser.add_argument("-X", "--sv", dest="sv", type=str,
                           help="Structure variation matrix.")
    cmdparser.add_argument("-C", "--cnv", dest="cnv", type=str,
                           help="cnv matrix. ")
    cmdparser.add_argument("-r", "--ref", dest="ref", type=str,
                           help="reference path. ")
    # cmdparser.add_argument("-i", "--info", dest="info", type=str,required=True,
    #                        help="sample information.")
    cmdparser.add_argument("-o", "--outputprefix", dest="outputprefix", type=str,
                           help="outputprefix")


    args = cmdparser.parse_args()
    part_len="/public/agis/chengshifeng_group/fengcong/WGRS/software/Fc-code/partlen.txt"


    ##deal with the args
    gene_name = args.gene if args.gene.find(".") != -1 else args.gene+".1" 
    gff_file = args.gff
    snp_file = args.SNP
    indel_file = args.indel
    sv_file = args.sv    #can be ignore
    cnv_file = args.cnv  #can be ignore
    # inf_file = args.info
    chr_name= args.chr
    ref_file = args.ref
    output_file = args.outputprefix

    ##gene_structure
    gene_structure=get_gene_structure(chr_name,gene_name,gff_file)

    ##maybe they have diff order
    snp_sample_order=[]
    indel_sample_order=[]
    sv_sample_order=[]
    cnv_sample_order=[]

    ##deal SNP matrix
    SNP_and_INF=deal_SNP(chr_name,snp_sample_order,gene_structure,snp_file)
    snp_sample_order=SNP_and_INF[2]
    
    #deal indel matrix
    INDEL_and_INF=deal_indel(chr_name,indel_sample_order,gene_structure,indel_file)
    indel_sample_order=INDEL_and_INF[2]
    # print(len(SNP_and_INF[0]),len(SNP_and_INF[1]))

    #deal SV matrix
    SV_and_INF=deal_sv(chr_name,sv_sample_order,gene_structure,sv_file,part_len)
    sv_sample_order=SV_and_INF[2]
    #filter SV sample and sort
    filter_and_sort_sv(SV_and_INF,snp_sample_order,sv_sample_order)

    #deal cnv matrix
    #/vol3/agis/chengshifeng_group/xianwenfei/wheat-GWAS/11.CNV/zCNV.vcf
    CN= deal_cnv(gene_name,cnv_sample_order,cnv_file)
    cnv_sample_order=CN[1]
    # print(CN)
    #filter cnv sample and sort
    CN=filter_and_sort_cnv(CN,snp_sample_order,cnv_sample_order)

    #hap clustering
    sample_order_and_hap_cluster_inf , cluster_inf = hap_clustering(snp_sample_order,SNP_and_INF[0],INDEL_and_INF[0],SV_and_INF[0],CN)
    sample_order_and_hap_cluster_str = ["%s(hap%d)"%(tu[0],tu[1]+1) for tu in sample_order_and_hap_cluster_inf]
    
    #generate each hap's sequence
    resolve_hap_seq(INDEL_and_INF[1],SV_and_INF[1],ref_file,chr_name,gene_structure,cluster_inf,SNP_and_INF[0],INDEL_and_INF[0],SV_and_INF[0],CN,output_file)
    

    ##test ouput 
    output_dict={
        "gene_structure":gene_structure,
        #"sample_name":snp_sample_order,
        "sample_name":sample_order_and_hap_cluster_str,
        "variation":{
            "SNP":SNP_and_INF[0],
            "INDEL":INDEL_and_INF[0],
            "SV":SV_and_INF[0],
            "CN":CN
            }
    }
    #js = json.dumps(output_dict, sort_keys=True, indent=4, separators=(',', ':'))
    js = json.dumps(output_dict)
    ouf = open(output_file+".json","w")
    ouf.write(js)
    ouf.close()


    #out info in csv format
    ##SNP
    ouf = open(output_file+".snp.csv","w")
    ouf.write("%s,%s,%s,%s,allele,Annotation,Gene_Name,Feature_ID,HGVS.c,HGVS.p\n"%("chr","position","ref","alt"))
    keys=["chr","position","ref","alt"] #+ "ann"
    for i in range(len(SNP_and_INF[1])):
        for key in keys:
            ouf.write("%s,"%( str(SNP_and_INF[1][i][key]) ) )
        
        ouf.write("%s\n"%(",".join(SNP_and_INF[1][i]["ann"][0])))

    ouf.close()

    ##indel 
    ouf = open(output_file+".indel.csv","w")
    ouf.write("%s,%s,%s,%s,%s,type,allele,Annotation,Gene_Name,Feature_ID,HGVS.c,HGVS.p\n"%("chr","position","length","ref","alt"))
    keys=["chr","position","length","ref","alt"] #+ "ann"
    for i in range(len(INDEL_and_INF[1])):
        for key in keys:
            ouf.write("%s,"%( str(INDEL_and_INF[1][i][key]) ) )
        ouf.write("%s,"%( str(INDEL_and_INF[0][i]["type"]) ) )
        ouf.write("%s\n"%(",".join(INDEL_and_INF[1][i]["ann"][0])))

    ouf.close()

    ##sv
    ouf = open(output_file+".sv.csv","w")
    ouf.write("%s,%s,%s,%s,%s,type\n"%("chr","position","length","ref","alt"))
    keys=["chr","position","length","ref","alt"] #+ "ann"
    for i in range(len(SV_and_INF[1])):
        for key in keys:
            ouf.write("%s,"%( str(SV_and_INF[1][i][key]) ) )
        ouf.write("%s\n"%( str(SV_and_INF[0][i]["type"]) ) )
        # ouf.write("%s\n"%(",".join(INDEL_and_INF[1][i]["ann"][0])))

    ouf.close()