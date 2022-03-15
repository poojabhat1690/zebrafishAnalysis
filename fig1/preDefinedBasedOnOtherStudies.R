### this script is used to define a list of high confidence list of maternal and zygotic genes based on previous 
### datasets as well as known GO terms. 
### Author : Pooja Bhat

library(reshape)
theme_ameres <- function (type) { 
  
  types_plot = c("boxplot","barplot","violin")
  if(type %in% types_plot ==T)
  {
    
    if(type == "boxplot"){
      theme(legend.title=element_blank(),axis.text.x = ggplot2::element_text(margin=ggplot2::margin(10,15,10,15,"pt"),size = 15),axis.text.y = element_text(margin=ggplot2::margin(5,15,10,5,"pt"),size = 15), axis.ticks.length = unit(-0.25 , "cm"),legend.position="bottom",axis.line = element_line(colour = "black", lineend = "round"), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))   
    }
    
    if(type == "barplot"){
      theme(legend.title=element_blank(),axis.text.x = ggplot2::element_text(margin=ggplot2::margin(10,15,10,15,"pt"),size = 15),axis.text.y = element_text(margin=ggplot2::margin(5,15,10,5,"pt"),size = 15), axis.ticks.length = unit(-0.25 , "cm"),legend.position="bottom",axis.line = element_line(colour = "black", lineend = "round"), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))   
    }
    
    
  }
  
}


#############defining a set of purely maternal genes
#### this rna seq data is taken from pauli et al 2012. 
#### 
classifiedGenes = read.table('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/packages/data/countingWindows_classified_conversions.txt',
                             sep="\t", header = T,stringsAsFactors = F)

expression_polyA = read.table('/Volumes/Macintosh HD/Users/pooja.bhat/Downloads//reads_STARmapping.txt',
                              sep="\t",stringsAsFactors = F, header = T)

expression_polyA = expression_polyA %>%  mutate_at(vars(contains('scratch')),.funs= funs(./(Length/1000))) %>%  
  mutate_at(vars(contains('scratch')),.funs= funs(./(sum(.)/1000000))) %>% 
  dplyr::mutate(cell_2 = (X.scratch.pooja.timecourse_polyA.SRR372787Aligned.sortedByCoord.out.bam + X.scratch.pooja.timecourse_polyA.SRR372788Aligned.sortedByCoord.out.bam)/2) %>%
  dplyr::mutate(hpf_28 = (X.scratch.pooja.timecourse_polyA.SRR372798Aligned.sortedByCoord.out.bam+X.scratch.pooja.timecourse_polyA.SRR372799Aligned.sortedByCoord.out.bam)/2 ) %>% 
  dplyr::mutate(oocyte_mean =  (X.scratch.pooja.timecourse_polyA.zf_oocyte_1_5Aligned.sortedByCoord.out.bam+X.scratch.pooja.timecourse_polyA.zf_oocyte_1_6Aligned.sortedByCoord.out.bam+
                                  X.scratch.pooja.timecourse_polyA.zf_oocyte_1_7Aligned.sortedByCoord.out.bam )/3) %>%
  dplyr::mutate(  bud_mean = (X.scratch.pooja.timecourse_polyA.SRR372796Aligned.sortedByCoord.out.bam+
       X.scratch.pooja.timecourse_polyA.SRR372797Aligned.sortedByCoord.out.bam)/2) %>% 
  dplyr::mutate(dpf5_mean = (X.scratch.pooja.timecourse_polyA.SRR372802Aligned.sortedByCoord.out.bam+X.scratch.pooja.timecourse_polyA.SRR372803Aligned.sortedByCoord.out.bam)/2) %>%
  dplyr::mutate(dome_mean=(X.scratch.pooja.timecourse_polyA.SRR372791Aligned.sortedByCoord.out.bam+X.scratch.pooja.timecourse_polyA.SRR372792Aligned.sortedByCoord.out.bam)/2) %>%
  dplyr::mutate(ratio_bud_2 = bud_mean/cell_2) %>%
  dplyr::mutate(ratio_28_2 = hpf_28/cell_2) %>%  dplyr::mutate(ratio_dome_2 = dome_mean/cell_2) %>% dplyr::mutate(ratio_5dpf_2 = dpf5_mean/cell_2) 
dr11_data =  read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/packages/data/figure2//allGenes_dr11.txt",
                        sep="\t",stringsAsFactors = F, header = T)
expression_polyA = expression_polyA %>% dplyr::mutate(ensembl_gene_id = Geneid)
expression_polyA = plyr::join(expression_polyA,dr11_data)

expression_polyA = expression_polyA[!duplicated(expression_polyA$external_gene_name),]


DecreasedExpression = expression_polyA %>% dplyr::filter(ratio_28_2<0.1  & cell_2 > 5)   
DecreasedExpression = DecreasedExpression[order(DecreasedExpression$cell_2,decreasing = T),] %>% dplyr::mutate(ensembl_gene_id = Geneid)
DecreasedExpression = DecreasedExpression[DecreasedExpression$external_gene_name %in% classifiedGenes$gene,]



GO_danRer = read.delim('/Volumes/Macintosh HD/Users/pooja.bhat/Downloads/GO_oogenesis_danRer.txt',stringsAsFactors = F)
GO_danRer_meiosis = read.delim('/Volumes/Macintosh HD/Users/pooja.bhat/Downloads/GO_meiosis_danRer.txt',stringsAsFactors = F)
GO_danRer_meiosis1 = read.delim('/Volumes/Macintosh HD/Users/pooja.bhat/Downloads/GO_meiosis1_danRer.txt',stringsAsFactors = F)

GO_danRer =  rbind.data.frame(GO_danRer,GO_danRer_meiosis,GO_danRer_meiosis1)

GO_danRer = GO_danRer[-which(GO_danRer$Gene==""),]
GO_danRer = GO_danRer[!duplicated(GO_danRer$Gene),]

GO_danRer$external_gene_name = GO_danRer$Gene
DecreasedExpression = plyr::join(DecreasedExpression,GO_danRer,by='external_gene_name')
DecreasedExpression = DecreasedExpression %>% dplyr::select(c(ensembl_gene_id,external_gene_name,gene_biotype,ratio_bud_2,ratio_28_2,ratio_dome_2,ratio_5dpf_2,cell_2,GO))
DecreasedExpression = DecreasedExpression[order(DecreasedExpression$cell_2,decreasing = T) ,]
oocyte = DecreasedExpression[DecreasedExpression$external_gene_name %in% GO_danRer$Gene,]

write.table(oocyte,'/Volumes/Macintosh HD/Users/pooja.bhat/Desktop//oocyte_meiosis.txt',sep="\t",
            col.names = T,quote = F,row.names = F)

write.table(DecreasedExpression,'/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data/oocyte_meiosis.txt',sep="\t",
            col.names = T,quote = F)



################## defining zygotic genes ####################

########################################## rRNA depleted libraries #####################################
riboMinus = read.table('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data/masterTable_quantSeqRNAseqRPMs.txt',
                       sep="\t",stringsAsFactors = F, header = T)
riboMinus = riboMinus %>% dplyr::mutate(ratio_ribocop = T1_riboCop/T4_riboCop)
a = riboMinus %>% dplyr::filter(T1_riboCop>0 | T4_riboCop>0) %>% dplyr::mutate(ratio_ribocop = T1_riboCop/T4_riboCop) %>% dplyr::filter(T1_riboCop<1) %>%
  dplyr::filter(ratio_ribocop<0.01) %>% dplyr::filter(T4_riboCop>5)

neuralCrestDevelopment = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Downloads/GO_neuralCrest.txt",sep="\t",stringsAsFactors = F)
neuralCrestDevelopment = neuralCrestDevelopment[neuralCrestDevelopment$V2 %in% a$external_gene_name,]

bmp_heart = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Downloads/GO_bmp_he.txt",sep="\t",stringsAsFactors = F)
bmp_heart[bmp_heart$V2 %in% a$external_gene_name,]
####################### Heyn et al 2014....


heynEtal = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/packages/data/figure2/heynEtal_s1.txt",
                      sep="\t",stringsAsFactors = F, header = T)
ensembl_geneNames_dr7 = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/packages/data/figure2/allGenes_dr7.txt",
                                   sep="\t",stringsAsFactors = F, header = T)
ensembl_geneNames_dr7 = ensembl_geneNames_dr7 %>% dplyr::distinct(ensembl_gene_id,.keep_all=T)
heynEtal = plyr::join(heynEtal,ensembl_geneNames_dr7)
heynEtal = heynEtal %>% dplyr::filter(gene_biotype=='protein_coding')



dr11 = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/packages/data/figure2//allGenes_dr11.txt",
                  sep="\t",stringsAsFactors = F, header = T)

`%!in%` = Negate(`%in%`)
heynEtal_dr11 = heynEtal[heynEtal$external_gene_name %in% dr11$external_gene_name,]
heynEtal_dr11 = heynEtal_dr11[-which(heynEtal_dr11$external_gene_name == ""),]
heynEtal_dr11 = heynEtal_dr11[-which(heynEtal_dr11$Gene_name  == ""),]


RPMs_ribominus = read.table('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data/masterTable_quantSeqRNAseqRPMs.txt',
                            stringsAsFactors = F, header = T, sep="\t")
#RPMs_ribominus %>% dplyr::mutate(ratio = (T4_riboCop+1)/(T1_riboCop+1) )%>% dplyr::filter(ratio<0.5)
RPMs_ribominus = RPMs_ribominus %>% dplyr::select(c('external_gene_name','T1_riboCop','T4_riboCop')) %>% dplyr::mutate(ration =(T4_riboCop+1)/(T1_riboCop+1) )
ZygoticGenes_heynEtal = plyr::join(heynEtal_dr11,RPMs_ribominus)

###### harvey et al., 2013

harveyEtal = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/packages/data/figure2/harveyEtal.txt",
                        sep="\t", stringsAsFactors = F, header = T)
harveyEtal_Z = harveyEtal %>% dplyr::filter(class == "Z")
harveyEtal_Z = data.frame(ensembl_gene_id = unlist(strsplit(harveyEtal_Z$Gene.ID,",",T)),class = "Z") %>% dplyr::distinct(ensembl_gene_id,.keep_all=T)
totalHarvey = rbind.data.frame(harveyEtal_Z)

dr11 = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/packages/data/figure2/allGenes_dr11.txt",
                  sep="\t",stringsAsFactors = F, header = T)

totalHarvey = plyr::join(totalHarvey,dr11) %>% dplyr::distinct(ensembl_gene_id,.keep_all=T)
totalHarvey = totalHarvey[complete.cases(totalHarvey$external_gene_name),]

####### lee et al...

leeEtal = read.delim("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/packages/data/figure2/leeEtal_complete.txt",
                     sep="\t", stringsAsFactors = F, header = T,skip = 23)
t4 = leeEtal %>% dplyr::filter(WT4_txed == 'E' | WT4_txed =='I' | WT4_txed == "EI") %>% dplyr::filter(Maternal_contr=="Z")
t6 = leeEtal %>% dplyr::filter(WT6_txed == 'E' | WT6_txed =='I' | WT6_txed == "EI") %>% dplyr::filter(Maternal_contr=="Z")
u1u2 = leeEtal %>% dplyr::filter(U1U2_txed == 'E' | U1U2_txed =='I' | U1U2_txed == "EI") %>% dplyr::filter(Maternal_contr=="Z")
leeEtal_Z = rbind.data.frame(t4,t6, u1u2)
leeEtal_Z = leeEtal_Z[!duplicated(leeEtal$Symbol),]

leeEtal_Z = leeEtal_Z %>% dplyr::mutate(external_gene_name = Symbol, ensembl_gene_id = Gene_id) %>%  
  dplyr::select(c('external_gene_name','ensembl_gene_id','Maternal_contr'))
dr11 = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/packages/data/figure2/allGenes_dr11.txt",
                  sep="\t",stringsAsFactors = F, header = T)
leeEtal_Z  = plyr::join(leeEtal_Z, dr11,by='external_gene_name') 

leeEtal_Z= leeEtal_Z[!duplicated(leeEtal_Z$external_gene_name),]
leeEtal_Z = leeEtal_Z[complete.cases(leeEtal_Z$external_gene_name),] #### using this only to get the genes which we can detect in dr11
leeEtal_Z$ensembl_gene_id =  NULL



heynEtal = ZygoticGenes_heynEtal %>% dplyr::mutate(study = 'heyn') %>% dplyr::select(external_gene_name,study)
leeEtal_Z = leeEtal_Z %>% dplyr::mutate(study = 'lee') %>% dplyr::select(external_gene_name,study)
totalHarvey = totalHarvey %>% dplyr::mutate(study = 'harvey') %>% dplyr::select(external_gene_name,study)

heyn_lee = intersect(heynEtal$external_gene_name,leeEtal_Z$external_gene_name)
harvey_lee = intersect(totalHarvey$external_gene_name,leeEtal_Z$external_gene_name)
heyn_harvey  = intersect(totalHarvey$external_gene_name,heynEtal$external_gene_name)


atLeastIn2 = unique(c(heyn_lee,harvey_lee,heyn_harvey))

atLeastIn2_Z = expression_polyA[expression_polyA$external_gene_name %in% atLeastIn2,]
atLeastIn2_Z = atLeastIn2_Z %>% 
  dplyr::select(-dplyr::contains('ratio')) %>% dplyr::mutate(ratio_cell2_hpf28 =  cell_2/hpf_28) %>% dplyr::filter(cell_2<1)    %>% dplyr::filter(dome_mean>5)  
  

allZgenes = rbind.data.frame(heynEtal,leeEtal_Z,totalHarvey)

allZgenes = allZgenes[!duplicated(allZgenes$external_gene_name),]
allZgenes = plyr::join(allZgenes,RPMs_ribominus) %>% dplyr::filter(T1_riboCop<1)


###################### defining maternal genes... 


harveyEtal_M = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/packages/data/figure2/harveyEral_Mgenes.txt",
                          sep="\t", stringsAsFactors = F, header = F)
colnames(harveyEtal_M) = c('ensembl_gene_id','class')
harveyEtal_M = data.frame(ensembl_gene_id = unlist(strsplit(harveyEtal_M$ensembl_gene_id,",",T)),class = "M") %>% dplyr::distinct(ensembl_gene_id,.keep_all=T)

dr11 = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/packages/data/figure2/allGenes_dr11.txt",
                  sep="\t",stringsAsFactors = F, header = T)

harveyEtal_M = plyr::join(harveyEtal_M,dr11) %>% dplyr::distinct(ensembl_gene_id,.keep_all=T)
harveyEtal_M = harveyEtal_M[complete.cases(harveyEtal_M$external_gene_name),]
harveyEtal_M = plyr::join(harveyEtal_M,RPMs_ribominus)



#### lee et al 
leeEtal = read.delim("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/packages/data/figure2/leeEtal_complete.txt",
                     sep="\t", stringsAsFactors = F, header = T,skip = 23)
leeEtal = leeEtal[!duplicated(leeEtal$Symbol),]

leeEtal = leeEtal %>% dplyr::mutate(external_gene_name = Symbol, ensembl_gene_id = Gene_id) %>%  
  dplyr::select(c('external_gene_name','ensembl_gene_id','Maternal_contr'))
dr11 = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/packages/data/figure2/allGenes_dr11.txt",
                  sep="\t",stringsAsFactors = F, header = T)
leeEtal  = plyr::join(leeEtal, dr11,by='external_gene_name') 

leeEtal= leeEtal[!duplicated(leeEtal$external_gene_name),]
leeEtal = leeEtal[complete.cases(leeEtal),] #### using this only to get the genes which we can detect in dr11
leeEtal$ensembl_gene_id =  NULL

leeEtal_M = leeEtal %>% dplyr::filter(Maternal_contr == "M")
leeEtal_M =  plyr::join(leeEtal_M,RPMs_ribominus) %>% dplyr::filter(ration <0.5)

leeEtal_M =leeEtal_M %>% dplyr::mutate(study='lee')%>% dplyr::select(external_gene_name,study,T1_riboCop, T4_riboCop, ration) %>% dplyr::mutate(type = "M")

##### 
allZgenes = allZgenes %>% dplyr::mutate(type = "Z")
preDefinedGenes= rbind.data.frame(leeEtal_M,allZgenes)
write.table(preDefinedGenes,"/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data//preDefGenes.txt",sep="\t",quote = F,row.names = F,col.names = T)

##### defining Zygotic and maternal genes... 

`%!in%` = Negate(`%in%`)
masterTable = read.table('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data/masterTable_quantSeqRNAseqRPMs.txt',
                         sep="\t",stringsAsFactors = F, header = T)
masterTable$ratio_expression = (masterTable$T4_riboCop+1)/(masterTable$T1_riboCop+1)
masterTable = masterTable %>% dplyr::select(c(external_gene_name,dplyr::contains('Inj_mean'),T1_riboCop:T4_riboCop,ratio_expression))
masterTable = masterTable %>% dplyr::filter(T4_riboCop>1 | T1_riboCop>1)
zygoticGenes = masterTable %>% dplyr::filter(ratio_expression>2 & T1_riboCop<1 &Inj_mean_TP9>2 ) %>% dplyr::mutate(type = 'Zygotic')
maternalGenes = masterTable %>% dplyr::filter(ratio_expression<0.5) %>% dplyr::mutate(type = 'Maternal')
nonZorM = masterTable[masterTable$external_gene_name %!in% c(zygoticGenes$external_gene_name,maternalGenes$external_gene_name),] %>% dplyr::mutate(type = 'InBetween')

totalGenes = rbind.data.frame(zygoticGenes,maternalGenes,nonZorM)
write.table(totalGenes,"/Volumes//Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data//allGenes.txt",sep="\t",quote = F)



