theme_ameres <- function (type) { 
  
  types_plot = c("boxplot","barplot","violin")
  if(type %in% types_plot ==T)
  {
    
    if(type == "boxplot"){
      theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt"),size = 15),axis.text.y = element_text(margin=margin(5,15,10,5,"pt"),size = 15), axis.ticks.length = unit(-0.25 , "cm"),legend.position="bottom",axis.line = element_line(colour = "black", lineend = "round"), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))   
    }
    
    if(type == "barplot"){
      theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt"),size = 15, hjust = 1),axis.text.y = element_text(margin=margin(5,15,10,5,"pt"),size = 15), axis.ticks.length = unit(-0.25 , "cm"),legend.position="bottom",axis.line = element_line(colour = "black", lineend = "round"), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18))   
    }
    
    
  }
  
}

sampleInfo = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data//barcodes_description.txt",
                        sep="\t",stringsAsFactors = F, header = F)
sampleInfo = sampleInfo %>% dplyr::mutate(sample = paste0("combinedFile_",V2,".fastq.gz_adapterTrimmed_slamdunk_mapped_filtered_tcperutr.csv"))

sample="TP9"
R2_sample = paste0('Inj_R2_',sample)
R3_sample = paste0('Inj_R3_',sample)

R2_TP= paste0(sampleInfo[which(sampleInfo$V3 == R2_sample),'sample'])
R3_TP= paste0(sampleInfo[which(sampleInfo$V3 == R3_sample),'sample'])

####################
### R2
reads_Z_R2 = read.table(paste0("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data//Zgenes/",R2_TP),
                   sep="\t",stringsAsFactors = F)

reads_M_R2 =  read.table(paste0("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data//Mgenes/",R2_TP),
                      sep="\t",stringsAsFactors = F)

reads_Mconversions_R2 =  read.table(paste0("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data//M_Conversiongenes//",R2_TP),
                         sep="\t",stringsAsFactors = F)

reads_MNOconversions_R2 =  read.table(paste0("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data//M_noConversiongenes///",R2_TP),
                                    sep="\t",stringsAsFactors = F)
reads_Mexampleconversions_R2 =  read.table(paste0("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data//MaternalGenesExamples_Conversiongenes///",R2_TP),
                                      sep="\t",stringsAsFactors = F)

reads_Z_example_R2 = read.table(paste0("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data//Zgenes_n1//",R2_TP),
                        sep="\t",stringsAsFactors = F)

reads_other_R2 = read.table(paste0("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data//OtherGenes_Conversiongenes//",R2_TP),
                                sep="\t",stringsAsFactors = F)
### R3

reads_Z_R3 = read.table(paste0("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data/Zgenes/",R3_TP),
                        sep="\t",stringsAsFactors = F)

reads_M_R3 =  read.table(paste0("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data//Mgenes/",R3_TP),
                         sep="\t",stringsAsFactors = F)

reads_Mconversions_R3 =  read.table(paste0("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data//M_Conversiongenes/",R3_TP),
                         sep="\t",stringsAsFactors = F)

reads_MNOconversions_R3 =  read.table(paste0("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data//M_noConversiongenes/",R3_TP),
                                      sep="\t",stringsAsFactors = F)
reads_Mexampleconversions_R3 =  read.table(paste0("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data//MaternalGenesExamples_Conversiongenes///",R3_TP),
                                           sep="\t",stringsAsFactors = F)

reads_Z_example_R3 = read.table(paste0("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data//Zgenes_n1//",R3_TP),
                                sep="\t",stringsAsFactors = F)

reads_other_R3 = read.table(paste0("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data//OtherGenes_Conversiongenes//",R3_TP),
                            sep="\t",stringsAsFactors = F)



Z_1 = reads_Z_R2 %>% dplyr::mutate(TC_minus = V4/V6, TC_plus = V3/V5,TC = TC_minus+ TC_plus, TC = TC*100) %>% dplyr::mutate(type = 'Z_1')
M_1 =  reads_M_R2 %>% dplyr::mutate(TC_minus = V4/V6, TC_plus = V3/V5,TC = TC_minus+ TC_plus, TC = TC*100) %>% dplyr::mutate(type = 'M_1')
MConv_1 = reads_Mconversions_R2 %>% dplyr::mutate(TC_minus = V4/V6, TC_plus = V3/V5,TC = TC_minus+ TC_plus, TC = TC*100) %>% dplyr::mutate(type = 'Mconv_1')
MNoConv_1 = reads_MNOconversions_R2 %>% dplyr::mutate(TC_minus = V4/V6, TC_plus = V3/V5,TC = TC_minus+ TC_plus, TC = TC*100) %>% dplyr::mutate(type = 'MNoconv_1')
Z_1_example = reads_Z_example_R2 %>% dplyr::mutate(TC_minus = V4/V6, TC_plus = V3/V5,TC = TC_plus, TC = TC*100) %>% dplyr::mutate(type = 'Z_1_example')### example gene is on minys
Mexample_1 = reads_Mexampleconversions_R2 %>% dplyr::mutate(TC_minus = V4/V6, TC_plus = V3/V5,TC = TC_minus+ TC_plus, TC = TC*100) %>% dplyr::mutate(type = 'Mexample')
Other_1 =  reads_other_R2 %>% dplyr::mutate(TC_minus = V4/V6, TC_plus = V3/V5,TC = TC_minus+ TC_plus, TC = TC*100) %>% dplyr::mutate(type = 'Other')

Z_2 = reads_Z_R3 %>% dplyr::mutate(TC_minus = V4/V6, TC_plus = V3/V5,TC = TC_minus+ TC_plus, TC = TC*100) %>% dplyr::mutate(type = 'Z_2')
M_2 =  reads_M_R3 %>% dplyr::mutate(TC_minus = V4/V6, TC_plus = V3/V5,TC = TC_minus+ TC_plus, TC = TC*100) %>% dplyr::mutate(type = 'M_2')
MConv_2 = reads_Mconversions_R3 %>% dplyr::mutate(TC_minus = V4/V6, TC_plus = V3/V5,TC = TC_minus+ TC_plus, TC = TC*100) %>% dplyr::mutate(type = 'Mconv_2')
MNoConv_2 = reads_MNOconversions_R3 %>% dplyr::mutate(TC_minus = V4/V6, TC_plus = V3/V5,TC = TC_minus+ TC_plus, TC = TC*100) %>% dplyr::mutate(type = 'MNoconv_2')
Z_2_example = reads_Z_example_R3 %>% dplyr::mutate(TC_minus = V4/V6, TC_plus = V3/V5,TC = TC_plus, TC = TC*100) %>% dplyr::mutate(type = 'Z_2_example')
Mexample_2 = reads_Mexampleconversions_R3 %>% dplyr::mutate(TC_minus = V4/V6, TC_plus = V3/V5,TC = TC_minus+ TC_plus, TC = TC*100) %>% dplyr::mutate(type = 'Mexample')
Other_2 =  reads_other_R3 %>% dplyr::mutate(TC_minus = V4/V6, TC_plus = V3/V5,TC = TC_minus+ TC_plus, TC = TC*100) %>% dplyr::mutate(type = 'Other')


Z = data.frame(position = c(-250:-1), TC= (Z_1$TC + Z_2$TC)/2, type = 'Z')
M = data.frame(position = c(-250:-1), TC= (M_1$TC + M_2$TC)/2, type = 'M')
M_conv = data.frame(position = c(-250:-1), TC= (MConv_1$TC + MConv_2$TC)/2, type = 'M_withConv')
M_Noconv = data.frame(position = c(-250:-1), TC= (MNoConv_1$TC + MNoConv_2$TC)/2, type = 'M_WithoutConv')
Z_examlple = data.frame(position = c(-250:-1), TC= (Z_1_example$TC + Z_2_example$TC)/2, type = 'Zexample')
M_example = data.frame(position = c(-250:-1), TC= (Mexample_1$TC + Mexample_2$TC)/2, type = 'M_example')
other_example = data.frame(position = c(-250:-1), TC= (Other_1$TC + Other_2$TC)/2, type = 'Other')



total_M = rbind.data.frame(Z,M)
total_MwithConv = rbind.data.frame(Z,M_conv)
total_MnoConv = rbind.data.frame(Z,M_Noconv)
total_Zexample = rbind.data.frame(Z,Z_examlple)
total_Mexample = rbind.data.frame(Z,M_example)
total_other = rbind.data.frame(Z,other_example)


pval_add_M= t.test(Z$TC,M$TC)$p.value
pval_add_Mcov= t.test(Z$TC,M_conv$TC)$p.value
pval_add_MNocov= t.test(Z$TC,M_Noconv$TC)$p.value
pval_add_Zexample= t.test(Z$TC,Z_examlple$TC)$p.value
pval_add_Mexample= t.test(Z$TC,M_example$TC)$p.value
pval_add_other= t.test(Z$TC,other_example$TC)$p.value

################################################

Zgenes = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data/preDefinedZgenes.bed",
                    sep="\t",stringsAsFactors = F)
Mgenes = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data/preDefinedMgenes.bed",
                    sep="\t",stringsAsFactors = F)
MNoConversionContainign = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data/preDefinedM_noConversionsgenes.bed",
                                     sep="\t",stringsAsFactors = F)
MconversionContaining = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data/preDefinedM_Conversionsgenes.bed",
                                   sep="\t",stringsAsFactors = F)
Mexamples = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data/preDefinedMselectedGenes_Conversionsgenes.bed",
                       sep="\t",stringsAsFactors = F)

Othergenes = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/data/preDefinedOther_Conversionsgenes.bed",
                       sep="\t",stringsAsFactors = F)




pdf('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/plots/perPositionConversions.pdf', height=3,width=5)
    
    perUTR = ggplot(total_M,aes(x=position,y=TC,group=type,color=type)) + geom_line()+ scale_color_brewer(palette = 'Set1')+ 
        ggplot2::annotate(geom='text',label = paste0("pval=",pval_add_M),x=-150,y=9) + ylim(c(0,12))+ xlab("Distance from 3' end " ) + 
        theme_cowplot() + ylab('T>C conversion rate') + ggtitle('M-all')
      perUTR = perUTR + theme_ameres(type='barplot')
      print(perUTR)
    
      perUTR_Mconv = ggplot(total_MwithConv,aes(x=position,y=TC,group=type,color=type)) + geom_line()+ scale_color_brewer(palette = 'Set1')+ 
        ggplot2::annotate(geom='text',label = paste0("pval=",pval_add_Mcov),x=-150,y=9) + ylim(c(0,12))+ xlab("Distance from 3' end " ) + 
        theme_cowplot() + ylab('T>C conversion rate')+ ggtitle('M with conversions')
      perUTR_Mconv = perUTR_Mconv + theme_ameres(type='barplot')
      print(perUTR_Mconv)
      
      
      perUTR_MNoconv = ggplot(total_MnoConv,aes(x=position,y=TC,group=type,color=type)) + geom_line()+ scale_color_brewer(palette = 'Set1')+ 
        ggplot2::annotate(geom='text',label = paste0("pval=",pval_add_MNocov),x=-150,y=9) + ylim(c(0,12))+ xlab("Distance from 3' end " ) + 
        theme_cowplot() + ylab('T>C conversion rate')+ ggtitle('M without conversions')
      perUTR_MNoconv = perUTR_MNoconv + theme_ameres(type='barplot')
      print(perUTR_MNoconv)
      
      onlyM = rbind.data.frame(total_MnoConv %>% dplyr::filter(type=="Z") ,
                               total_MnoConv %>% dplyr::filter(type=="M_WithoutConv"), 
                               total_MwithConv %>% dplyr::filter(type=="M_withConv"),
                               total_Mexample %>% dplyr::filter(type=="M_example"), total_other%>% dplyr::filter(type=="Other") )
      
      perUTR_AllMexample = ggplot(onlyM,aes(x=position,y=TC,group=type,color=type)) + geom_line()+ scale_color_brewer(palette = 'Set1')+ 
        ggplot2::annotate(geom='text',label = paste0("pval M example=",pval_add_Mexample),x=-150,y=9) + 
        ggplot2::annotate(geom='text',label = paste0("pval M all=",pval_add_M),x=-150,y=8)+
        ggplot2::annotate(geom='text',label = paste0("pval M with conv =",pval_add_Mcov),x=-150,y=7) +
        ggplot2::annotate(geom='text',label = paste0("pval M without conv =",pval_add_MNocov),x=-150,y=6) +
        ggplot2::annotate(geom='text',label = paste0("pval other =",pval_add_other),x=-150,y=5) +
        ylim(c(0,20))+ xlab("Distance from 3' end " ) + 
        theme_cowplot() + ylab('T>C conversion rate')+ ggtitle('Z example')
      perUTR_AllMexample = perUTR_AllMexample + theme_ameres(type='barplot')
      print(perUTR_AllMexample)
      
      n <- 5
      nr <- nrow(onlyM)
     onlyM_split =  split(onlyM, rep(1:ceiling(nr/n), each=n, length.out=nr))
     onlyM_split = lapply(onlyM_split,function(x) x %>% dplyr::mutate(range = c(paste0(max(position),":",min(position)))) %>% dplyr::mutate(TC_mean = mean(TC)) ) 
     onlyM_split = lapply(onlyM_split,function(x) x[1,])
     onlyM_split = do.call(rbind.data.frame,onlyM_split)
     
     perUTR_AllMexample = ggplot(onlyM_split,aes(x=position,y=TC_mean,group=type,color=type)) + geom_line()+ geom_point()+scale_color_brewer(palette = 'Set1')+ 
       ggplot2::annotate(geom='text',label = paste0("pval M example=",pval_add_Mexample),x=-150,y=9) + 
       ggplot2::annotate(geom='text',label = paste0("pval M all=",pval_add_M),x=-150,y=8)+
       ggplot2::annotate(geom='text',label = paste0("pval M with conv =",pval_add_Mcov),x=-150,y=7) +
       ggplot2::annotate(geom='text',label = paste0("pval M without conv =",pval_add_MNocov),x=-150,y=6) +
       ggplot2::annotate(geom='text',label = paste0("pval other =",pval_add_other),x=-150,y=5) +
       ylim(c(0,11))+ xlab("Distance from 3' end " ) + 
       theme_cowplot() + ylab('T>C conversion rate')+ ggtitle('5 ntAverage')
     perUTR_AllMexample = perUTR_AllMexample + theme_ameres(type='barplot')
     print(perUTR_AllMexample)
     
     n <- 10
     nr <- nrow(onlyM)
     onlyM_split =  split(onlyM, rep(1:ceiling(nr/n), each=n, length.out=nr))
     onlyM_split = lapply(onlyM_split,function(x) x %>% dplyr::mutate(range = c(paste0(max(position),":",min(position)))) %>% dplyr::mutate(TC_mean = mean(TC)) ) 
     onlyM_split = lapply(onlyM_split,function(x) x[1,])
     onlyM_split = do.call(rbind.data.frame,onlyM_split)
     onlyM_split$rangeVal  = c(1:25)
     perUTR_AllMexample = ggplot(onlyM_split,aes(x=rangeVal,y=TC_mean,group=type,color=type)) + geom_line()+ geom_point()+scale_color_brewer(palette = 'Set1')+ 
       ggplot2::annotate(geom='text',label = paste0("pval M example=",pval_add_Mexample),x=15,y=9) + 
       ggplot2::annotate(geom='text',label = paste0("pval M all=",pval_add_M),x=15,y=8)+
       ggplot2::annotate(geom='text',label = paste0("pval M with conv =",pval_add_Mcov),x=15,y=7) +
       ggplot2::annotate(geom='text',label = paste0("pval M without conv =",pval_add_MNocov),x=15,y=6) +
       ggplot2::annotate(geom='text',label = paste0("pval other =",pval_add_other),x=15,y=5) +
       ylim(c(0,11))+ xlab("Distance from 3' end " ) + 
       theme_cowplot() + ylab('T>C conversion rate')+ ggtitle('10 ntAverage')
     perUTR_AllMexample = perUTR_AllMexample + theme_ameres(type='barplot')
     print(perUTR_AllMexample)
     
     
     my_comparisons = list(c('Z','M_WithoutConv'),c('Z','M_withConv'),c('Z','M_example'),c('Z','Other'))
       
     perUTR_Mexample = ggplot(total_Mexample,aes(x=position,y=TC,group=type,color=type)) + geom_line()+ scale_color_brewer(palette = 'Set1')+ 
        ggplot2::annotate(geom='text',label = paste0("pval=",pval_add_Mexample),x=-150,y=9) + ylim(c(0,12))+ xlab("Distance from 3' end " ) + 
        theme_cowplot() + ylab('T>C conversion rate')+ ggtitle('M examples')
      perUTR_Mexample = perUTR_Mexample + theme_ameres(type='barplot')
      print(perUTR_Mexample)
      
      

       perUTR_Zexample = ggplot(total_Zexample,aes(x=position,y=TC,group=type,color=type)) + geom_line()+ scale_color_brewer(palette = 'Set1')+ 
        ggplot2::annotate(geom='text',label = paste0("pval=",pval_add_Zexample),x=-150,y=9) + ylim(c(0,12))+ xlab("Distance from 3' end " ) + 
        theme_cowplot() + ylab('T>C conversion rate')+ ggtitle('Z example')
      perUTR_Zexample = perUTR_Zexample + theme_ameres(type='barplot')
      print(perUTR_Zexample)
      
dev.off()

### quantifying these in one jitter
pdf('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/plots/perpositionConversions_jitter.pdf', height=3, width=4)
    jitter_perPosition = ggplot(onlyM_split,aes(x=type,y=TC_mean)) + geom_quasirandom(varwidth = T) + 
      ggpubr::stat_compare_means(data=onlyM_split,comparisons = my_comparisons,label = "p.signif") + theme_cowplot() 
    jitter_perPosition =jitter_perPosition + theme_ameres(type = 'barplot')
    print(jitter_perPosition)
dev.off()


### getting allThe bed files.. 
allCategories= list(Zgenes,Mgenes,MNoConversionContainign,MconversionContaining,Mexamples,Othergenes)
names(allCategories) = c('Z','allM','M_noConv','M_conv','M_example','Othergenes')
allCategories = reshape2::melt(lapply(allCategories,nrow))
