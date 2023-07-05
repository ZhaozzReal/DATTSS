from argparse import ArgumentParser,ArgumentTypeError
import HTSeq
from multiprocessing import Pool
import numpy as np
parser = ArgumentParser(description = "Detect dynamic tandem TSS usage from standard RNA-seq")
parser.add_argument("-b",dest = 'bamfiles',action = "store",type = str,help = "Input text file with all bamfiles")
parser.add_argument('-anno',dest = 'anno_txt',action = "store",type = str,help = "Input annotation file contains first exon regions and annotated tandem TSSs within them")
parser.add_argument("-p",dest = "processors",action = "store",default = 10,type = int,help = "<INT> Number of processors used [default: 10]")
parser.add_argument("-o",dest = "outfile",action = "store",type = str,help = "Output all inferred tandem TSS events and the quantification of tandem TSS usage")
args = parser.parse_args()




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
            ratio = round(aUTR/cUTR,3)
            if ratio < 1:
                ratio_list.append(ratio)
            else:
                ratio_list.append(None)
        else:
            ratio_list.append(None)
    return ratio_list




def Get_sample_dict(gdcfile,path):
    tumor_dict = {}
    for line in open(gdcfile,"r"):
        if "Sample Type" not in line:
            tumor_dict[line.strip().split("\t")[6]] = path + line.strip().split("\t")[1]
    return tumor_dict


bamfile_txt = args.bamfiles
bamfiles = open(bamfile_txt,"r").readlines()[0].strip().split(',')
bamnames = [ i.split("/")[-1] for i in bamfiles]


out = open(args.outfile,"w")
out.write("{}\t{}\t{}\t{}\t{}\n".format("genename","first_exon_region","Proximal_TSS","MSE_ratio","\t".join(bamnames)))
pool = Pool(args.processors)
for line in open(args.anno_txt,"r"):
    SYMBOL,first_exon,strand,Annotated_TSSs = line.strip().split("\t")
    chrom = first_exon.split(":")[0]
    if strand == "+":
        first_exon_end = int(first_exon.split(":")[1].split("-")[1])
    else:
        first_exon_end = int(first_exon.split(":")[1].split("-")[0])
    exon_bamfiles_list = list(zip([first_exon]*len(bamfiles),bamfiles))
    all_cvg_list = pool.map(Get_region_cvg_list,exon_bamfiles_list)
    merged_cvg = np.sum(all_cvg_list,axis = 0).tolist()
    coverage_threshold = 30
    if max(merged_cvg) > coverage_threshold*len(bamfiles):
        if strand == "+":
            merged_cvg = merged_cvg[::-1]
        Proximal_TSS,min_mse_ratio,exon_length = Get_proximal_TSS(merged_cvg,Annotated_TSSs,first_exon_end)
        if Proximal_TSS != "NA":
            point = abs(Proximal_TSS - first_exon_end)
            upstream_abundance = round(np.mean(merged_cvg[point:]),3)
            downstream_abundance = round(np.mean(merged_cvg[:point]),3)
            if strand == "-":
                all_cvg_list = [ cvg[:exon_length] for cvg in all_cvg_list]
                first_exon = chrom + ":" + str(first_exon_end) + "-" + str(first_exon_end + exon_length) + ":" + strand
            else:
                all_cvg_list = [ cvg[::-1][:exon_length] for cvg in all_cvg_list]
                first_exon = chrom + ":" + str(first_exon_end - exon_length) + "-" + str(first_exon_end) + ":" + strand
            ratio_list = Cal_distalTSS_usage(point,all_cvg_list)
            Proximal_TSS = chrom + ":" + str(Proximal_TSS)
            out.write("{}\t{}\t{}\t{}\t{}\n".format(SYMBOL,first_exon,Proximal_TSS,min_mse_ratio,"\t".join(list(map(str,ratio_list)))))


out.close()
