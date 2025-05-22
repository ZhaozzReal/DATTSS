from argparse import ArgumentParser,ArgumentTypeError
import HTSeq,subprocess,os
from multiprocessing import Pool
import numpy as np
parser = ArgumentParser(description = "Dynamic analysis of alternative tandem TSS usage from standard RNA-seq, comparison between different conditions")
parser.add_argument("-b",dest = 'bamfiles',action = "store",type = str,help = "Input text file with all bamfiles")
parser.add_argument('-anno',dest = 'anno_txt',action = "store",type = str,help = "Input annotation file contains first exon regions and annotated tandem TSSs within them")
parser.add_argument("-p",dest = "processors",action = "store",default = 10,type = int,help = "<INT> Number of processors used [default: 10]")
parser.add_argument("-r",dest = "exonRegion",action = "store",type = str,help = "Input annotation file contains annotated exon regions")
parser.add_argument("-d",dest = "exonCountDir",action = "store",type = str,help = "The directory for storing files with exon count")
parser.add_argument("-o",dest = "outfile",action = "store",type = str,help = "Output all inferred tandem TSS events and the quantification of tandem TSS usage")
args = parser.parse_args()


def parse_cfgfile(bamfile_txt):
    for line in open(bamfile_txt,"r"):
        lines = line.strip().split("=")
        if lines[0] == "condition1":
            bamfiles_condition1 = lines[1].split(",")
        if lines[0] == "condition2":
            bamfiles_condition2 = lines[1].split(",")
    return bamfiles_condition1,bamfiles_condition2


def Get_region_cvg_list(input):
    region,bamfile = input
    chrom,start_end = region.split(":")
    start,end = start_end.split("-")
    try:
        bam_reader = HTSeq.BAM_Reader(bamfile)
        ga = HTSeq.GenomicArray("auto", stranded=False, typecode="i")
        read_seq = bam_reader.fetch(region = region)
        for a in read_seq:
            iv_seq = (cigop.ref_iv for cigop in a.cigar if cigop.type == "M" and cigop.size >0)
            for iv in iv_seq:
                ga[iv] += 1
        region_iv = HTSeq.GenomicInterval(chrom,int(start),int(end),"-")
        cvg = list(ga[region_iv])
    except:
        cvg = [0]*(int(end) - int(start))
    return cvg


def Estimation_abundance(Region_Coverage,break_point):
    downstream_cov_mean = np.mean(Region_Coverage[break_point:])
    upstream_cov_mean = np.mean(Region_Coverage[:break_point])
    Coverage_diff = Region_Coverage[:break_point]-upstream_cov_mean
    Coverage_diff = np.append(Coverage_diff,Region_Coverage[break_point:]-downstream_cov_mean)
    Mean_Squared_error = np.mean(Coverage_diff**2)
    return Mean_Squared_error


def Get_proximal_TSS(cvg,peaks,first_exon_end):
    point_list = [abs(int(i) - int(first_exon_end)) for i in peaks.split("_")]
    for point in sorted(point_list,reverse = True):
        downcov = np.mean(cvg[point:])
        upcov = np.mean(cvg[:point])
        ratio = downcov/upcov
        if len(cvg[point:]) > 30:
            downratio = len(list(filter(lambda x:x>int(downcov/3),cvg[point:-30])))/len(cvg[point:-30])
        elif len(cvg[point:]) > 10:
            downratio = len(list(filter(lambda x:x>int(downcov/3),cvg[point:-10])))/len(cvg[point:-10])
        else:
            downratio = len(list(filter(lambda x:x>int(downcov/3),cvg[point:])))/len(cvg[point:])
        if ratio < 0.05 or downratio < 0.7:
            cvg = cvg[:point]
    global_mse = np.mean((cvg - np.mean(cvg))**2)
    peaks_list = []
    mse_ratio_list = []
    min_peak = "NA"
    min_mse_ratio = "NA"
    for i in peaks.split("_"):
        point = abs(int(i) - int(first_exon_end))
        if len(cvg) > point:
            mse_ratio = Estimation_abundance(cvg, point)/global_mse
            peaks_list.append(int(i))
            mse_ratio_list.append(mse_ratio)
    if mse_ratio_list != []:
        min_mse_ratio = min(mse_ratio_list)
        if min_mse_ratio < 0.5:
            min_mse_ratio_index = mse_ratio_list.index(min_mse_ratio)
            min_peak = peaks_list[min_mse_ratio_index]
    return min_peak,min_mse_ratio,len(cvg)



def Cal_distalTSS_usage(point,all_cvg_list):
    ratio_list = []
    coverage_threshold = 20
    for cvg in all_cvg_list:
        cUTR = np.mean(sorted(cvg[:point],reverse = True)[:30])
        if cUTR > coverage_threshold:
            dis = min(max(point,100),int(len(cvg[point:])/2))
            aUTR = np.mean(sorted(cvg[point:],reverse = True)[:dis])
            ratio = round(aUTR/(cUTR + 0.1),3)
            if ratio < 1:
                ratio_list.append(ratio)
            else:
                ratio_list.append(None)
        else:
            ratio_list.append(None)
    return ratio_list



def Get_exon_count(input_tuple):
    bamfile,exon_information = input_tuple
    bam_reader = HTSeq.BAM_Reader(bamfile)
    if "distal" in exon_information:
        genename,chrom,start_end,strand,distal = exon_information.split(":")
        region_fetch = chrom + ":" + start_end
        read_seq = bam_reader.fetch(region = region_fetch)
        count = len([a for a in read_seq])
        label = genename + ":" + distal 
    else:
        genename,chrom,start_end,strand = exon_information.split(":")
        region_fetch = chrom + ":" + start_end
        read_seq = bam_reader.fetch(region = region_fetch)
        count = len([a for a in read_seq])
        label = genename + ":" + chrom + "_" + start_end
    return label,count


def output_exoncount(dir_path,all_bamfiles,exonRegion,outfile,processors):
    gene_exons_dict = {}
    for line in open(exonRegion,"r"):
        genename,exon_region = line.strip().split("\t")
        if genename not in gene_exons_dict:
            gene_exons_dict[genename] = [genename + ":" + exon_region]
        else:
            gene_exons_dict[genename].append(genename + ":" + exon_region)
    exon_region_list = []
    for line in open(outfile,"r"):
        SYMBOL,first_exon,Proximal_TSS = line.strip().split("\t")[:3]
        if "chr" not in first_exon:
            continue
        chrom,start_end,strand = first_exon.split(":")
        start,end = start_end.split("-")
        Proximal_TSS = Proximal_TSS.split(":")[1]
        if strand == '+':
            distal_exon = SYMBOL + ":" + chrom + ":" + start + "-" + Proximal_TSS + ":" + strand + ":" + "distal" + "_" + Proximal_TSS
            proximal_exon = SYMBOL + ":" + chrom + ":" + Proximal_TSS + "-" + end + ":" + strand
        else:
            distal_exon = SYMBOL + ":" + chrom + ":" + Proximal_TSS + "-" + end + ":" + strand + ":" + "distal" + "_" + Proximal_TSS
            proximal_exon = SYMBOL + ":" + chrom + ":" + start + "-" + Proximal_TSS + ":" + strand
        exon_region_list.append(distal_exon)
        exon_region_list.append(proximal_exon)
        if SYMBOL in gene_exons_dict:
            exon_region_list.extend(gene_exons_dict[SYMBOL])
    exon_region_list = list(set(exon_region_list))
    if os.path.exists(dir_path) == False:
        os.makedirs(dir_path)
    for bamfile in all_bamfiles:
        out = open( dir_path  + "/" + bamfile.split("/")[-1].split(".")[0] + "_exoncount.txt","w")
        from multiprocessing import Pool
        pool = Pool(processors)
        input_tuple = list(zip([bamfile]*len(exon_region_list),exon_region_list))
        result_list = pool.map(Get_exon_count,input_tuple)
        for exon_lst in result_list:
            out.write("{}\t{}\n".format(exon_lst[0],exon_lst[1]))
        out.close()


def DATTSS_main(outfile,anno_txt,processors,all_bamfiles):
    bamnames = [i.split("/")[-1] for i in all_bamfiles]
    out = open(outfile,"w")
    out.write("{}\t{}\t{}\t{}\t{}\n".format("genename","first_exon_region","proximal_TSS","MSE_ratio","\t".join(bamnames)))
    pool = Pool(processors)
    for line in open(anno_txt,"r"):
        SYMBOL,first_exon,strand,Annotated_TSSs = line.strip().split("\t")
        chrom = first_exon.split(":")[0]
        if strand == "+":
            first_exon_end = int(first_exon.split(":")[1].split("-")[1])
        else:
            first_exon_end = int(first_exon.split(":")[1].split("-")[0])
        exon_bamfiles_list = list(zip([first_exon]*len(all_bamfiles),all_bamfiles))
        all_cvg_list = pool.map(Get_region_cvg_list,exon_bamfiles_list)
        merged_cvg = np.sum(all_cvg_list,axis = 0).tolist()
        coverage_threshold = 30
        if max(merged_cvg) > coverage_threshold*len(all_bamfiles):
            if strand == "+":
                merged_cvg = merged_cvg[::-1]
            Proximal_TSS,min_mse_ratio,exon_length = Get_proximal_TSS(merged_cvg,Annotated_TSSs,first_exon_end)
            if Proximal_TSS != "NA":
                point = abs(Proximal_TSS - first_exon_end)
                if strand == "-":
                    all_cvg_list = [ cvg[:exon_length] for cvg in all_cvg_list]
                    first_exon = chrom + ":" + str(first_exon_end) + "-" + str(first_exon_end + exon_length) + ":" + strand
                else:
                    all_cvg_list = [ cvg[::-1][:exon_length] for cvg in all_cvg_list]
                    first_exon = chrom + ":" + str(first_exon_end - exon_length) + "-" + str(first_exon_end) + ":" + strand
                ratio_list = Cal_distalTSS_usage(point,all_cvg_list)
                Proximal_TSS = chrom + ":" + str(Proximal_TSS)
                out.write("{}\t{}\t{}\t{}\t{}\n".format(SYMBOL,first_exon,Proximal_TSS,round(min_mse_ratio,3),"\t".join(list(map(str,ratio_list)))))
                out.flush()


bamfiles_condition1,bamfiles_condition2 = parse_cfgfile(args.bamfiles)
all_bamfiles = bamfiles_condition1 + bamfiles_condition2
DATTSS_main(args.outfile,args.anno_txt,args.processors,all_bamfiles)
output_exoncount(args.exonCountDir,all_bamfiles,args.exonRegion,args.outfile,args.processors)
