#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   pipe_dada2_picrust2.py
@Time    :   2024/07/16 15:48:09
@Author  :   Wanjin.Hu 
@Version :   1.0
@Contact :   wanjin.hu@outlook.com
@Description :
'''

import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser(description="DADA2 and picrust2 analysis")
parser.add_argument("-i", "--input_dir", dest="inDir", required=True, type=str, help="input dir")
parser.add_argument("-o", "--output_dir", dest="outDir", required=True, type=str, help="output dir")
parser.add_argument("-r", "--ref", dest="refDB", required=False, type=str, help="silva db")
args = parser.parse_args()

class DADA2(object):
    """
    """
    def __init__(self):
        super(DADA2,self).__init__()
        self.r_env = "/root/miniconda3/envs/R-4.2.1/bin/Rscript"
        self.bash_env = "/usr/bin/bash"
        self.pipe_env = "/root/pipe_script/dada2"
        self.cmd_dada2 = self.pipe_env + "/scripts/dada2.R"
        self.cmd_pic = self.pipe_env + "/scripts/picrust2.sh"
        self.silva = self.pipe_env + "/scripts/silva_nr99_v138.1_train_set.fa.gz"
    
    def run_dada2(self):
        """运行 dada2"""
        if args.refDB:
            ref_db = args.refDB
        else:
            ref_db = self.silva
        asv_fa = args.outDir + "/ASV.fasta"
        if os.path.exists(asv_fa):
            print("ASV.fasta 文件已存在，跳过 dada2 分析")
        else:
            print("ASV.fasta 文件不存在，开始 dada2 分析")
            cmd = "{} {} --indir {} --refdb {} --outdir {}".format(self.r_env,
                                                                   self.cmd_dada2,
                                                                   args.inDir,
                                                                   ref_db,
                                                                   args.outDir)
            os.system(command=cmd)
        # asv_table.csv 格式修改
        in_asvTable = args.outDir + "/ASV_table.csv"
        out_asvTable = args.outDir + "/ASV_table.xls"
        df_asvTable = pd.read_csv(in_asvTable, sep = ",", header = 0, index_col = 0)
        df_asvTable.T.to_csv(out_asvTable,
                             index = True, 
                             index_label = "#ASV", 
                             sep = "\t")
            
    def run_taxonomy_trim(self,tax_raw,tax_trim):
        """物种 taxonomy 格式修改"""
        if os.path.exists(tax_trim):
            print("物种 taxonomy 文件已存在，跳过格式修改")
        else:
            with open(tax_raw,"r") as f1, open(tax_trim,"w") as f2:
                sample_line = f1.readline() # 首行
                sample_line = sample_line.strip().split("\t")
                sample_title = "\t".join(sample_line)
                f2.write(sample_title + "\n")
                for line in f1:
                    line = line.strip().split("\t")
                    sample_list = line[7:]
                    sample = "\t".join(sample_list)
                    kindom = line[0]
                    phylum = line[1]
                    classs = line[2]
                    order = line[3]
                    family = line[4]
                    genus = line[5]
                    asv = line[6]
                    if genus != "g__NA":
                        if genus != "g__Incertae Sedis":
                            genus_t = genus
                            family_t = family
                            order_t = order
                            class_t = classs
                            phylum_t = phylum
                        if genus == "g__Incertae Sedis":
                            genus_t = family + "_" + genus
                            family_t = family
                            order_t = order
                            class_t = classs
                            phylum_t = phylum
                    if genus == "g__NA":
                        if family != "f__NA":
                            genus_t = family + "_g__unclassified"
                            family_t = family
                            order_t = order
                            class_t = classs
                            phylum_t = phylum
                        if family == "f__NA":
                            if order!= "o__NA":
                                genus_t = order + "_g__unclassified"
                                family_t = order + "_f__unclassified"
                                order_t = order
                                class_t = classs
                                phylum_t = phylum
                            if order == "o__NA":
                                if classs!= "c__NA":
                                    genus_t = classs + "_g__unclassified"
                                    family_t = classs + "_f__unclassified"
                                    order_t = classs + "_o__unclassified"
                                    class_t = classs
                                    phylum_t = phylum
                                if classs == "c__NA":
                                    if phylum!= "p__NA":
                                        genus_t = phylum + "_g__unclassified"
                                        family_t = phylum + "_f__unclassified"
                                        order_t = phylum + "_o__unclassified"
                                        class_t = phylum + "_c__unclassified"
                                        phylum_t = phylum
                                    if phylum == "p__NA":
                                        genus_t = kindom + "_g__unclassified"
                                        family_t = kindom + "_f__unclassified"
                                        order_t = kindom + "_o__unclassified"
                                        class_t = kindom + "_c__unclassified"
                                        phylum_t = kindom + "_p__unclassified"
                    sp_list = [kindom,phylum_t,class_t,order_t,family_t,genus_t,asv]
                    sp = "\t".join(sp_list)
                    f2.write(sp + "\t" + sample + "\n")
                f2.close()

    def run_picrust2(self):
        """运行 picrust2 功能预测"""
        pic_out = args.outDir + "/Picrust2"
        asv_fa = args.outDir + "/ASV.fasta"
        asv_table = args.outDir + "/ASV_table.xls"
        if os.path.exists(pic_out):
            print("picrust2 功能结果已存在，跳过 picrust2 分析")
        else:
            print("运行 picrust2 功能预测")
            cmd = "{} {} {} {} {}".format(self.bash_env,
                                          self.cmd_pic,
                                          asv_fa,
                                          asv_table,
                                          pic_out)
            os.system(command = cmd)

    def main(self):
        self.run_dada2()
        tax_raw = args.outDir + "/SamplesTaxonomyComposition.txt"
        tax_trim = args.outDir + "/SamplesTaxonomyComposition_trim.txt"
        self.run_taxonomy_trim(tax_raw, tax_trim)
        self.run_picrust2()

if __name__ == "__main__":
    run = DADA2()
    run.main()