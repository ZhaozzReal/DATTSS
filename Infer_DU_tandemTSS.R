suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(DEXSeq))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
args = commandArgs(TRUE)
option.list = list(
  make_option(c("-b","--bam"),action = "store_true",default = FALSE,help = "input file contains all filenames of bamfile [%default]"),
  make_option(c("-I","--IPUI"),action = "store_true",default = FALSE,help = "input file generated in the previous step which contains all detected tandem TSS events [%default]"),
  make_option(c("-d","--dir"),action = "store_true",default = FALSE,help = "input directory name contains all exoncount files [%default]"),
  make_option(c("-o","--output"),action = "store_true",default = FALSE,help = "output final results[%default]")
)
desc = "Infer significantly dysregulated tandem TSS usage between conditions using DEXSeq model"
parser = OptionParser(option_list = option.list, description = desc)
opt = parse_args(parser, args = args, positional_arguments = TRUE)
cfg_file = read.table(opt$args[1])
filenames = apply(cfg_file,1,function(x)strsplit(x,"=") %>% unlist() %>% .[2] %>% strsplit(",") %>% unlist()) %>% as.character() %>% 
  lapply(function(x)strsplit(x,".",fixed=T) %>% unlist() %>% .[1]) %>% unlist()
condition_length = apply(cfg_file,1,function(x)strsplit(x,"=") %>% unlist() %>% .[2] %>% strsplit(",") %>% unlist() %>% length())
sampleTable = data.frame(row.names = filenames,condition = c(rep("condition1",condition_length[1]),rep("condition2",condition_length[2])))
countFiles = list.files(opt$args[3], pattern = "exoncount.txt", full.names=TRUE)
dxd = DEXSeqDataSetFromHTSeq(countFiles,sampleData = sampleTable,design = ~ sample + exon + condition:exon)
dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd )
dxd = testForDEU( dxd )
dexresult = DEXSeqResults( dxd ) %>% as.data.frame() %>% dplyr::filter(grepl("distal",featureID)) %>% dplyr::select("groupID","featureID","padj") %>% rowwise() %>% dplyr::mutate(featureID=strsplit(featureID,"_",fixed=T) %>% unlist() %>% .[2]) %>% dplyr::mutate(anno=paste0(groupID,":",featureID)) %>% dplyr::select(anno,padj)

DTUI_table = read.table(opt$args[2],header = T)
DTUI_table = DTUI_table[rowSums(DTUI_table[,5:ncol(DTUI_table)]!="None")==length(filenames),]
DTUI_table[,5:ncol(DTUI_table)] = apply(DTUI_table[,5:ncol(DTUI_table)],2,as.numeric)
DTUI_table$DTUI_condition1 = round(rowSums(DTUI_table[,5:(5+condition_length[1]-1)])/condition_length[1],3)
DTUI_table$DTUI_condition2 = round(rowSums(DTUI_table[,(5+condition_length[1]):(5+condition_length[1]+condition_length[2]-1)])/condition_length[2],3)
DTUI_table$DTUI_diff = round(DTUI_table$DTUI_condition2 - DTUI_table$DTUI_condition1,3)
DTUI_table = DTUI_table %>% rowwise() %>% dplyr::mutate(anno=paste0(genename,":",strsplit(as.character(proximal_TSS),":",fixed=T) %>% unlist() %>% .[2])) %>% 
  inner_join(dexresult)
final = DTUI_table %>% dplyr::mutate(change=ifelse(padj<0.05&abs(DTUI_diff)>0.05,ifelse(DTUI_diff>0,"Lengthening","Shortening"),"Nochange")) %>% dplyr::distinct(genename,proximal_TSS,.keep_all=T) %>% dplyr::select(-"anno")
write.table(final,opt$args[4],quote = FALSE,sep = "\t",row.names = FALSE)
