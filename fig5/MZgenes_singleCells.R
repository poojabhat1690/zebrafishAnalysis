library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(reshape)

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
`%!in%` = Negate(`%in%`)

splitReplicates = function(dataFrameToSplit,condition,metadata_add){
  dataFrameToSplit_condition = dataFrameToSplit[,grep(condition,colnames(dataFrameToSplit ))]
  # dataFrameToSplit_condition_R1 = dataFrameToSplit_condition[,grep("R1",colnames(dataFrameToSplit_condition ))]
  dataFrameToSplit_condition_R2 = dataFrameToSplit_condition[,grep("R2",colnames(dataFrameToSplit_condition ))]
  dataFrameToSplit_condition_R3 = dataFrameToSplit_condition[,grep("R3",colnames(dataFrameToSplit_condition ))]
  mean_repl = (dataFrameToSplit_condition_R2+dataFrameToSplit_condition_R3)/2
  #dataFrameToSplit_condition_R1 = cbind.data.frame(dataFrameToSplit_condition_R1,metadata_add)
  dataFrameToSplit_condition_R2 = cbind.data.frame(dataFrameToSplit_condition_R2,metadata_add)
  dataFrameToSplit_condition_R3 = cbind.data.frame(dataFrameToSplit_condition_R3,metadata_add)
  mean_repl = cbind.data.frame(mean_repl,metadata_add)
  splitReplicates = list(dataFrameToSplit_condition_R2,dataFrameToSplit_condition_R3,mean_repl)
  names(splitReplicates) = c("R2","R3","mean")
  return(splitReplicates)
}

allFractions = list.files("//Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure5/data/singleCellData//fractionOfcells/",pattern = "fractionCells*")
allFractions_path = paste0("//Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure5/data/singleCellData//fractionOfcells/", allFractions)

# allFractions = list.files("//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/MZdynamics/singleCell/fractionOfcells/",pattern = "fractionCells*")
# allFractions_path = paste0("//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/MZdynamics/singleCell/fractionOfcells/", allFractions)
allFractions_data = lapply(allFractions_path,function(x) read.table(x, stringsAsFactors = F, header = T))
names(allFractions_data) = allFractions

classifiedCountingWindows = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/packages/data//countingWindows_classified_conversions.txt",sep="\t",stringsAsFactors = F, header = T)
classifiedCountingWindows = classifiedCountingWindows %>% dplyr::mutate(genes = toupper(gene)) %>% dplyr::select(c('genes','class','description'))
allFractions_data = lapply(allFractions_data,function(x) plyr::join(x,classifiedCountingWindows))
allFractions_data = lapply(allFractions_data, function(x) x %>% dplyr::mutate(fractionCells = numberOfCellsWithSignal/numberOfcells) %>%
                             dplyr::mutate(signalPerCell = sumSignal/numberOfcells) %>% dplyr::mutate(signalPerExpressedCell = 
                                                                                                        sumSignal/numberOfCellsWithSignal          )) 
metadata = allFractions_data$fractionCells_dome.txt[,c("genes","class","description")]

fractionOfcells = do.call(cbind.data.frame,lapply(allFractions_data,function(x) x$fractionCells)) %>% dplyr::select(c('fractionCells_high.txt','fractionCells_oblong.txt',  'fractionCells_dome.txt','fractionCells_epiboly30.txt','fractionCells_epiboly50.txt', 'fractionCells_shieldStage.txt', 
                                                                                                                      'fractionCells_Epiboly_60.txt','fractionCells_Epiboly_75.txt','fractionCells_Epiboly_90.txt','fractionCells_budStage.txt','fractionCells_somite3Stage.txt','fractionCells_somite6Stage.txt'))
sumSignal = do.call(cbind.data.frame,lapply(allFractions_data,function(x) x$sumSignal)) %>% dplyr::select(c('fractionCells_high.txt','fractionCells_oblong.txt',  'fractionCells_dome.txt','fractionCells_epiboly30.txt','fractionCells_epiboly50.txt', 'fractionCells_shieldStage.txt', 
                                                                                                            'fractionCells_Epiboly_60.txt','fractionCells_Epiboly_75.txt','fractionCells_Epiboly_90.txt','fractionCells_budStage.txt','fractionCells_somite3Stage.txt','fractionCells_somite6Stage.txt'))
numberOfcellsWithExpression =  do.call(cbind.data.frame,lapply(allFractions_data,function(x) x$numberOfCellsWithSignal)) %>% dplyr::select(c('fractionCells_high.txt','fractionCells_oblong.txt',  'fractionCells_dome.txt','fractionCells_epiboly30.txt','fractionCells_epiboly50.txt', 'fractionCells_shieldStage.txt', 
                                                                                                                                             'fractionCells_Epiboly_60.txt','fractionCells_Epiboly_75.txt','fractionCells_Epiboly_90.txt','fractionCells_budStage.txt','fractionCells_somite3Stage.txt','fractionCells_somite6Stage.txt'))

signalPerExpressedCell = do.call(cbind.data.frame,lapply(allFractions_data,function(x) x$signalPerExpressedCell)) %>% dplyr::select(c('fractionCells_high.txt','fractionCells_oblong.txt',  'fractionCells_dome.txt','fractionCells_epiboly30.txt','fractionCells_epiboly50.txt', 'fractionCells_shieldStage.txt', 
                                                                                                                                      'fractionCells_Epiboly_60.txt','fractionCells_Epiboly_75.txt','fractionCells_Epiboly_90.txt','fractionCells_budStage.txt','fractionCells_somite3Stage.txt','fractionCells_somite6Stage.txt'))
signalPerCell = do.call(cbind.data.frame,lapply(allFractions_data,function(x) x$signalPerCell)) %>% dplyr::select(c('fractionCells_high.txt','fractionCells_oblong.txt',  'fractionCells_dome.txt','fractionCells_epiboly30.txt','fractionCells_epiboly50.txt', 'fractionCells_shieldStage.txt', 
                                                                                                                    'fractionCells_Epiboly_60.txt','fractionCells_Epiboly_75.txt','fractionCells_Epiboly_90.txt','fractionCells_budStage.txt','fractionCells_somite3Stage.txt','fractionCells_somite6Stage.txt'))

colnames(fractionOfcells) = paste0('fractionOfcells_',c('high','oblong','dome','epiboly30','epiboly50','shield','epiboly60','epiboly75','epiboly90','bud','somite3','somite6'))
colnames(sumSignal) = paste0('sumSignal_',c('high','oblong','dome','epiboly30','epiboly50','shield','epiboly60','epiboly75','epiboly90','bud','somite3','somite6'))
colnames(numberOfcellsWithExpression) = paste0('numberOfcellsWithExpression_',c('high','oblong','dome','epiboly30','epiboly50','shield','epiboly60','epiboly75','epiboly90','bud','somite3','somite6'))
colnames(signalPerExpressedCell) = paste0('signalPerExpressedCell_',c('high','oblong','dome','epiboly30','epiboly50','shield','epiboly60','epiboly75','epiboly90','bud','somite3','somite6'))
colnames(signalPerCell) = paste0('ssignalPerCell_',c('high','oblong','dome','epiboly30','epiboly50','shield','epiboly60','epiboly75','epiboly90','bud','somite3','somite6'))

allData = cbind.data.frame(fractionOfcells,sumSignal,numberOfcellsWithExpression,signalPerExpressedCell,signalPerCell)
masterTable_numberOfcells = cbind.data.frame(metadata,allData)
masterTable_numberOfcells$numberOfubiquitousStages = masterTable_numberOfcells %>% dplyr::select(dplyr::contains('fraction')) %>%
  apply(.,1,function(x) length(which(x>0.75)))

genes_MZ = masterTable_numberOfcells %>% dplyr::filter(description == "MZ pure")
genes_MZ = genes_MZ[!duplicated(genes_MZ$genes),]
genes_MZ = genes_MZ %>% dplyr::mutate(ubiquitousness = ifelse(fractionOfcells_high>0.75,T,F))

##### also adding if these are fast or slow... 
fractionReads = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures//figure4/data//clusteredBasedOnTCaccumulation.bed",
                           sep="\t",stringsAsFactors = F)
fractionReads = fractionReads %>% dplyr::select(c(V10,V11)) %>% dplyr::mutate(genes= toupper(V10)) %>% dplyr::select(c(genes,V11))
genes_MZ = plyr::join(genes_MZ,fractionReads)
genes_MZ = genes_MZ[!duplicated(genes_MZ$genes),]

pdf('/Volumes/Macintosh HD/Users/pooja.bhat/Desktop/tmp.pdf', height = 5, width=5)
  plot(t(fractionReads[which(fractionReads$V10 == 'tob1a'),c(1:9)]),ylab='fraction TC reads')
  plot(t(fractionReads[which(fractionReads$V10 == 'xbp1'),c(1:9)]),ylab='fraction TC reads')
dev.off()
##### now also getting and processing quantSeq roms
RPMs_all = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures//figure2/data/RPM_allCws.txt",sep="\t",stringsAsFactors = F,header=T)


RPMs_all = RPMs_all %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Inj_R3_TP9), sum, na.rm = TRUE)
RPMs_all_split = splitReplicates(dataFrameToSplit = RPMs_all,condition = "Inj",metadata_add = RPMs_all[,1])
RPMs_all_split$mean$genes  = toupper(RPMs_all_split$mean$name)
genes_MZ = plyr::join(genes_MZ,RPMs_all_split$mean,by='genes')
genes_MZ = genes_MZ[!duplicated(genes_MZ$genes),]

ubiquitous_geneExpression = reshape::melt(genes_MZ %>% dplyr::filter(ubiquitousness == T) %>% dplyr::select(dplyr::contains("Inj"))) %>% dplyr::mutate(class = "Ubiquitous")
NONubiquitous_geneExpression = reshape::melt(genes_MZ %>% dplyr::filter(ubiquitousness == F) %>% dplyr::select(dplyr::contains("Inj"))) %>% dplyr::mutate(class = "NONUbiquitous")
totalGeneExpression = rbind.data.frame(ubiquitous_geneExpression,NONubiquitous_geneExpression)
totalGeneExpression$value = log10(totalGeneExpression$value+1)
totalGeneExpression = totalGeneExpression %>% dplyr::group_by(class) %>% dplyr::mutate(number = n()/9) %>% dplyr::mutate(class_number = paste(class,"_",number))

pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure5/plots//ubiquitous_nonUbiquitousGenes.pdf", height=5,width=8)
            ggpubr::ggviolin(totalGeneExpression,x='variable',y='value',fill='class_number',palette = 'Set1',add='boxplot',outlier.shape = NA,ylab='log10 (RPM)') + 
            ggpubr::stat_compare_means(aes(group = class),label = "p.signif") + theme_ameres(type = 'barplot')
          

          normalizedMZ = t(apply(genes_MZ %>% dplyr::select(dplyr::contains("Inj")) ,1,function(x) x/max(x)))
          colnames(normalizedMZ) = paste0("Normalized_",colnames(normalizedMZ))
          genes_MZ = cbind.data.frame(genes_MZ,normalizedMZ)
          
          
          ####### plotting the gene expression for ubiquitous and non ubiquitous genes
          genes_MZ_ubiquitous = genes_MZ %>% dplyr::filter(ubiquitousness == T) %>%dplyr::select(dplyr::contains('ssignalPerCell')) 
          colnames(genes_MZ_ubiquitous) = c(3.3,3.8,4.3,4.8,5.3,6,7,8,9,10,11,12)
          genes_MZ_ubiquitous = melt(genes_MZ_ubiquitous)  %>% dplyr::mutate(class='Ubiquitous')
          
          genes_MZ_NONubiquitous = genes_MZ %>% dplyr::filter(ubiquitousness == F) %>%dplyr::select(dplyr::contains('ssignalPerCell')) 
          colnames(genes_MZ_NONubiquitous) =  c(3.3,3.8,4.3,4.8,5.3,6,7,8,9,10,11,12)
          genes_MZ_NONubiquitous = melt(genes_MZ_NONubiquitous)  %>% dplyr::mutate(class='nonUbiquitous')
          
          totalExpressionPercell = rbind.data.frame(genes_MZ_NONubiquitous,genes_MZ_ubiquitous)
          ggpubr::ggviolin(totalExpressionPercell,x='variable',y='value',fill='class',palette = 'Set1',add='boxplot',ylab=' gene expression per cell',xlab = 'Time') +
            theme_ameres(type = 'barplot')
          
          
          #######=sifnal per expressed cell
          genes_MZ_ubiquitous = genes_MZ %>% dplyr::filter(ubiquitousness == T) %>%dplyr::select(dplyr::contains('signalPerExpressedCell')) 
          colnames(genes_MZ_ubiquitous) =  c(3.3,3.8,4.3,4.8,5.3,6,7,8,9,10,11,12)
          genes_MZ_ubiquitous = melt(genes_MZ_ubiquitous)  %>% dplyr::mutate(class='Ubiquitous')
          
          genes_MZ_NONubiquitous = genes_MZ %>% dplyr::filter(ubiquitousness == F) %>%dplyr::select(dplyr::contains('signalPerExpressedCell')) 
          colnames(genes_MZ_NONubiquitous) =  c(3.3,3.8,4.3,4.8,5.3,6,7,8,9,10,11,12)
          genes_MZ_NONubiquitous = melt(genes_MZ_NONubiquitous)  %>% dplyr::mutate(class='nonUbiquitous')
          
          totalExpressionPercell_signalPerExpressedCell = rbind.data.frame(genes_MZ_NONubiquitous,genes_MZ_ubiquitous)
          ggpubr::ggviolin(totalExpressionPercell_signalPerExpressedCell,x='variable',y='value',fill='class',palette = 'Set1',add='boxplot',ylab=' Signal per expressed cell',xlab = 'Time') +
            theme_ameres(type = 'barplot')
          


dev.off()

####### ubiquitous

pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure5/plots//fastAndSlowUbiquitousgenes.pdf", height = 6)
       
        genes_MZ_ubiquitous = genes_MZ %>% dplyr::filter(ubiquitousness == T)
        genes_MZ_split= split(genes_MZ,genes_MZ$V11,T)
        
        genes_MZ_ubiquitous_quantSeq = lapply(genes_MZ_split,function(x) x %>% dplyr::filter(ubiquitousness == T)%>% dplyr::select(dplyr::contains("Normalized")) )
        genes_MZ_ubiquitous_quantSeq = lapply(genes_MZ_ubiquitous_quantSeq, setNames, c(0.75,2,2.5,3,3.5,4,4.5,5,5.5))
        genes_MZ_ubiquitous_quantSeq = reshape::melt(genes_MZ_ubiquitous_quantSeq)
        
        ggpubr::ggviolin(genes_MZ_ubiquitous_quantSeq,x="variable",y='value',group="variable" ,ylab= "Normalized RPM expression",
                         title = nrow(genes_MZ_ubiquitous_quantSeq)/9,facet.by = 'L1',ylim=c(0,1), fill='L1', palette = 'Set1',trim = T, width=1.5)
        
        
        #### gene expression.. 
        genes_MZ_ubiquitous_quantSeq = lapply(genes_MZ_split,function(x) x %>% dplyr::filter(ubiquitousness == T)%>% dplyr::select(dplyr::contains("Inj_"))  ) 
        genes_MZ_ubiquitous_quantSeq = lapply(genes_MZ_ubiquitous_quantSeq,function(x) x[,c(1:9)])
        genes_MZ_ubiquitous_quantSeq = lapply(genes_MZ_ubiquitous_quantSeq, setNames, c(0.75,2,2.5,3,3.5,4,4.5,5,5.5))
        genes_MZ_ubiquitous_quantSeq = reshape::melt(genes_MZ_ubiquitous_quantSeq)
        genes_MZ_ubiquitous_quantSeq$value = log10(genes_MZ_ubiquitous_quantSeq$value)
      
          ggpubr::ggviolin(genes_MZ_ubiquitous_quantSeq,x="variable",y='value',group="variable" ,ylab= "Normalized RPM expression",
                         title = nrow(genes_MZ_ubiquitous_quantSeq)/9,facet.by = 'L1', 
                         fill='L1', palette = 'Set1',trim = T, width=1) + theme_ameres(type = 'barplot')
        
          ggpubr::ggviolin(genes_MZ_ubiquitous_quantSeq,x='L1',y='value',facet.by = 'variable',add='boxplot') + 
            ggpubr::stat_compare_means(ref.group = '1',label = "p.signif")
        
        
        genes_MZ_ubiquitous_fractionOfCells = lapply(genes_MZ_split,function(x) x %>% dplyr::filter(ubiquitousness == T) %>% dplyr::select(dplyr::contains("fractionOfcells")) )
        genes_MZ_ubiquitous_fractionOfCells = lapply(genes_MZ_ubiquitous_fractionOfCells, setNames,  c(3.3,3.8,4.3,4.8,5.3,6,7,8,9,10,11,12))
        genes_MZ_ubiquitous_fractionOfCells = reshape::melt(genes_MZ_ubiquitous_fractionOfCells)
  
        nGenes= data.frame(table(genes_MZ_ubiquitous_fractionOfCells$L1)/12)
        colnames(nGenes) = c('L1','freq')        
        genes_MZ_ubiquitous_fractionOfCells = plyr::join(genes_MZ_ubiquitous_fractionOfCells,nGenes)
        genes_MZ_ubiquitous_fractionOfCells$L1 = paste(genes_MZ_ubiquitous_fractionOfCells$L1,genes_MZ_ubiquitous_fractionOfCells$freq,sep="_")
        fractionOfcells_clusters =  ggpubr::ggviolin(genes_MZ_ubiquitous_fractionOfCells,x="variable",y='value',
                         group="variable",ylab = "Fraction of expressed cells",
                         title = nrow(genes_MZ_ubiquitous_quantSeq)/9,facet.by = 'L1',ylim=c(0,1),
                         palette = 'Set1', fill='L1',add = "median",,trim = T,width = 1) + theme_ameres(type = 'barplot') 
        print(fractionOfcells_clusters)
        summary_median = genes_MZ_ubiquitous_fractionOfCells %>% group_by(L1,variable) %>% dplyr::summarise(median_val = median(value))
        
        summary_mediaFractions  = genes_MZ_ubiquitous_fractionOfCells %>% dplyr::group_by(L1,variable) %>% dplyr::summarise(median_fraction = median(value))
        
        medianFraction = ggpubr::ggline(summary_mediaFractions,x='variable',y='median_fraction',color='L1',ylim=c(0,1),palette = 'Set1') + 
          theme_ameres(type = 'barplot') 
        print(medianFraction)
        
        genes_MZ_ubiquitous_singlecell = lapply(genes_MZ_split,function(x) as.data.frame(t(apply(x %>% dplyr::filter(ubiquitousness==T)%>% dplyr::select(dplyr::contains("ssignalPerCell")) ,1,function(y) y))))
        genes_MZ_ubiquitous_singlecell = lapply(genes_MZ_ubiquitous_singlecell, setNames, c(3.3,3.8,4.3,4.7,5.3,6,7,8,9,10,11,12))
        genes_MZ_ubiquitous_singlecell = reshape::melt(genes_MZ_ubiquitous_singlecell)
        ggpubr::ggviolin(genes_MZ_ubiquitous_singlecell,x="variable",y='value',group="variable",ylab = "Signal per cel",title = nrow(genes_MZ_ubiquitous_quantSeq)/9,facet.by = 'L1',
                         palette = 'Set1', fill='L1') + theme_ameres(type = 'barplot')
        
        genes_MZ_ubiquitous_Perexpressedsinglecell = lapply(genes_MZ_split,function(x) as.data.frame(t(apply(x %>% dplyr::filter(ubiquitousness==T)%>% dplyr::select(dplyr::contains("signalPerExpressedCell")) ,1,function(y) y))))
        genes_MZ_ubiquitous_Perexpressedsinglecell = lapply(genes_MZ_ubiquitous_Perexpressedsinglecell, setNames, c(3.3,3.8,4.3,4.7,5.3,6,7,8,9,10,11,12))
        genes_MZ_ubiquitous_Perexpressedsinglecell = reshape::melt(genes_MZ_ubiquitous_Perexpressedsinglecell)
        
        ggpubr::ggviolin(genes_MZ_ubiquitous_Perexpressedsinglecell,x="variable",y='value',group="variable",ylab = "Signal per expressed cell",
                         title = nrow(genes_MZ_ubiquitous_quantSeq)/9,
                         facet.by = 'L1', palette = 'Set1', fill='L1', trim = T)  + theme_ameres(type = 'barplot')
        
        ####
        
        genes_MZ_ubiquitous_sumSignal = lapply(genes_MZ_split,function(x) as.data.frame(t(apply(x %>% dplyr::filter(ubiquitousness==T)%>% dplyr::select(dplyr::contains("sumSignal")) ,1,function(y) y))))
        genes_MZ_ubiquitous_sumSignal = lapply(genes_MZ_ubiquitous_sumSignal, setNames, c(3.3,3.8,4.3,4.7,5.3,6,7,8,9,10,11,12))
        genes_MZ_ubiquitous_sumSignal = reshape::melt(genes_MZ_ubiquitous_sumSignal)
        genes_MZ_ubiquitous_sumSignal$value = log10(genes_MZ_ubiquitous_sumSignal$value)
       ggpubr::ggviolin(genes_MZ_ubiquitous_sumSignal,x="variable",y='value',group="variable",ylab = "Total signal",
                         title = nrow(genes_MZ_ubiquitous_quantSeq)/9,facet.by = 'L1', 
                         palette = 'Set1', fill='L1', trim = T) + theme_ameres(type = 'barplot')
        
        
        

dev.off()


clusters_genes = (lapply(genes_MZ_split,function(x) x[,c('genes','V11')]))
clusters_genes = do.call(rbind.data.frame,clusters_genes)
clusters_genes$genes = tolower(clusters_genes$genes)
write.table(clusters_genes,"/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure5/data/clusters_samples_ubiquitousSinglecell.txt",
            sep="\t",quote = F,row.names = F,col.names = F)

####### check if these go into a lineage


epiboly_final = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure5/data/connectCellsToLineage//finalInfo_shield.txt")
numberOfcells = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure5/data/connectCellsToLineage//possibleLineages_allCells_shield.txt",sep="\t",stringsAsFactors = F)
lineages = epiboly_final[,grep('Lineage',colnames(epiboly_final))]
rownames(lineages) = epiboly_final$gene

for(i in 1:ncol(lineages)){
  lineages[,i] = (lineages[,i])/numberOfcells[1,i]
}

combinedSegments_epiboly = lineages %>% dplyr::mutate(EVL = Lineage_EVL, Marginal = (Lineage_Adaxial_Cells+Lineage_Somites+Lineage_Hematopoeitic_ICM +
                                                                               Lineage_Hematopoeitic_RBI_Pronephros+ Lineage_Endoderm_Pharyngeal+Lineage_Endoderm_Pancreatic_Intestinal+
                                                                               Lineage_Heart_Primordium + Lineage_Cephalic_Mesoderm)/8, MarginalDorsal = Lineage_Prechordal_Plate, Dorsal=Lineage_Notochord,
                                              Ventral = Lineage_Tailbud, PGC = Lineage_Primordial_Germ_Cells, DorsalAnimal = (Lineage_Diencephalon+ Lineage_Optic_Cup+ Lineage_Midbrain_Neural_Crest+ Lineage_Hindbrain_R3+ Lineage_Hindbrain_R4_5_6+ Lineage_Telencephalon+ Lineage_Neural_Plate_Border+ Lineage_Placode_Adeno._Lens_Trigeminal+ 
                                                                                                                                Lineage_Placode_Epibranchial_Otic+Lineage_Placode_Olfactory)/10,
                                              VentralAnimal = Lineage_Epidermis) %>% dplyr::select(c(EVL,Marginal,MarginalDorsal,Ventral,PGC,DorsalAnimal,VentralAnimal, Dorsal))


rownames(combinedSegments_epiboly) = rownames(lineages)



getLineageInteractions = function(testData,numberOfcells_test,combinedSegments_reference,genes_MZ){
  
  lineages = testData[,grep('Lineage',colnames(testData))]
  rownames(lineages) = testData$gene
  
  for(i in 1:ncol(lineages)){
    lineages[,i] = (lineages[,i])/numberOfcells_test[1,i]
  }
  
  combinedSegments_test = lineages %>% dplyr::mutate(EVL = Lineage_EVL, Marginal = (Lineage_Adaxial_Cells+Lineage_Somites+Lineage_Hematopoeitic_ICM +
                                                                                      Lineage_Hematopoeitic_RBI_Pronephros+ Lineage_Endoderm_Pharyngeal+Lineage_Endoderm_Pancreatic_Intestinal+
                                                                                      Lineage_Heart_Primordium + Lineage_Cephalic_Mesoderm)/8, MarginalDorsal = Lineage_Prechordal_Plate, Dorsal=Lineage_Notochord,
                                                     Ventral = Lineage_Tailbud, PGC = Lineage_Primordial_Germ_Cells, DorsalAnimal = (Lineage_Diencephalon+ Lineage_Optic_Cup+ Lineage_Midbrain_Neural_Crest+ Lineage_Hindbrain_R3+ Lineage_Hindbrain_R4_5_6+ Lineage_Telencephalon+ Lineage_Neural_Plate_Border+ Lineage_Placode_Adeno._Lens_Trigeminal+ 
                                                                                                                                       Lineage_Placode_Epibranchial_Otic+Lineage_Placode_Olfactory)/10,
                                                     VentralAnimal = Lineage_Epidermis) %>% dplyr::select(c(EVL,Marginal,MarginalDorsal,Ventral,PGC,DorsalAnimal,VentralAnimal,Dorsal))
  
  combinedSegments_test = t(apply(combinedSegments_test,1,function(x) x/max(x)))
  
  rownames(combinedSegments_test) = rownames(lineages)
  
  
  
  
  combinedSegments_reference_MZ = combinedSegments_reference[rownames(combinedSegments_reference) %in% genes_MZ_ubiquitous$genes,]
  combinedSegments_test_MZ = combinedSegments_test[rownames(combinedSegments_test) %in% genes_MZ_ubiquitous$genes,]
  
  for(i in 1:ncol(combinedSegments_reference)){
    #combinedSegments_reference[,i] = log2((combinedSegments_reference[,i]+1)/(combinedSegments_test[,i]+1))
    combinedSegments_reference[,i] = combinedSegments_test[,i]
  }
  
  
  
  
  genes_MZ_split= split(genes_MZ,genes_MZ$V11,T)
  genes_MZ_split = lapply(genes_MZ_split,function(x) x %>% dplyr::filter(ubiquitousness ==T ))
  
  combinedSegments_reference$genes = rownames(combinedSegments_test)
  rownames(combinedSegments_reference) = rownames(combinedSegments_test)
  
  genes_MZ_split = lapply(genes_MZ_split,function(x) plyr::join(x,combinedSegments_reference))
  
  lineages_split = lapply(genes_MZ_split,function(x) x %>%dplyr::select(EVL:Dorsal,genes) %>% tibble::column_to_rownames(.,var = 'genes')) 
  genes_MZ_split = lapply(genes_MZ_split,function(x) x[complete.cases(x),])
  return(lineages_split)
  
}

### refernece stage ... 
finalTable = list.files(path = "/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure5/data//connectCellsToLineage/",pattern = "finalInfo")
finalTable_path = paste0("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure5/data/connectCellsToLineage/",finalTable)
finalTable_data  = lapply(finalTable_path,function(x) read.table(x,stringsAsFactors = F))
names(finalTable_data) = finalTable

#####number of cells
numberOfcells_files = list.files("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure5/data/connectCellsToLineage/",pattern = "possibleLineages_allCells")
numberOfcells_files_path = paste0("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure5/data/connectCellsToLineage/",numberOfcells_files)
numberOfcells = lapply(numberOfcells_files_path,function(x) read.table(x, stringsAsFactors = F))

names(numberOfcells) = numberOfcells_files
allList = vector('list',length(numberOfcells))
names(allList) = names(numberOfcells)

for(i in 1:length(finalTable_data)){
  allList[[i]] = getLineageInteractions(testData = finalTable_data[[i]],numberOfcells_test = numberOfcells[[i]],combinedSegments_reference = combinedSegments_epiboly,genes_MZ = genes_MZ)
  
}

clus1  = lapply(allList,function(x) x[[1]])
clus2  = lapply(allList,function(x) x[[2]])
clus3  = lapply(allList,function(x) x[[3]])
clus4  = lapply(allList,function(x) x[[4]])

### checking how many are specific to a lineage 
library(factoextra)
library(NbClust)
lineageSpecificities = lapply(allList$possibleLineages_allCells_shield.txt,function(x) apply(x,1,function(x) length((which(x>0.75)))))

shieldSpefic = do.call(rbind.data.frame,allList$possibleLineages_allCells_shield.txt)
Z_MZshieldSpecific = rbind.data.frame(clusterThis_Z,shieldSpefic)
clus = kmeans(Z_MZshieldSpecific,centers = 2)

pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure5/plots/Clusters_lineageFractions.pdf",
    height=4, width=4)
  cluster_dimensions = fviz_cluster(clus, data = Z_MZshieldSpecific,geom = 'point') + theme_cowplot()
  cluster_dimensions  = cluster_dimensions + theme_ameres(type = 'barplot')
  print(cluster_dimensions)
  clusters_elbow = fviz_nbclust(Z_MZshieldSpecific, kmeans, method = "wss") +
    geom_vline(xintercept = 2, linetype = 2) + theme_ameres(type = 'barplot')
  print(clusters_elbow)
dev.off()
  
  Z_MZshieldSpecific$cluster = clus$cluster
Z_MZshieldSpecific$class = substr(rownames(Z_MZshieldSpecific),start = 1,stop = 1)
Z_MZshieldSpecific$class = as.numeric(Z_MZshieldSpecific$class)
Z_MZshieldSpecific$class[is.na(Z_MZshieldSpecific$class)]<-0

cluster1 = Z_MZshieldSpecific[Z_MZshieldSpecific$cluster==1,c(1:8)]

cluster1_melt = melt(cluster1)
cluster1_melt$gene = rownames(cluster1)
ggplot(cluster1_melt,aes(x=variable,y=value)) + geom_violin()

##### dividing these into different cluster specific genes.

plots_lineageSpecific = vector('list',8)
names(plots_lineageSpecific) = colnames(cluster1)
linePlot_clus1 = plots_lineageSpecific

for(i in 1:length(cluster1)){
  cluster1_tmp = melt(cluster1[which(cluster1[,i]==1),])
  plots_lineageSpecific[[i]] = ggplot(cluster1_tmp,aes(x=variable,y=value)) + geom_violin() + 
     ggtitle(colnames(cluster1)[i]) + theme_cowplot()
  plots_lineageSpecific[[i]] = plots_lineageSpecific[[i]]  + theme_ameres(type = 'barplot') + 
    ylab(label = 'Normalized lineage fraction')
  linePlot_clus1[[i]] = cluster1_tmp %>% dplyr::group_by(variable) %>% dplyr::summarise(median(value))
  
}

linePlot_clus1 = melt(linePlot_clus1)
ggpubr::ggline(linePlot_clus1,x='variable',y='value',color='L1',facet.by = 'L1') + ylim(c(0,1))+ theme_ameres(type = 'barplot')


cluster2 = Z_MZshieldSpecific[Z_MZshieldSpecific$cluster==2,c(1:8)]

cluster2_melt = melt(cluster2)
cluster2_melt$gene = rownames(cluster2)
ggplot(cluster2_melt,aes(x=variable,y=value)) + geom_violin() + theme_ameres(type = 'barplot')

##### dividing these into different cluster specific genes.

plots_lineageSpecific_clus2 = vector('list',8)
names(plots_lineageSpecific_clus2) = colnames(cluster2)
linePlot_clus2 = plots_lineageSpecific_clus2
for(i in 1:length(cluster2)){
  cluster2_tmp = melt(cluster2[which(cluster2[,i]==1),])
  plots_lineageSpecific_clus2[[i]] = ggplot(cluster2_tmp,aes(x=variable,y=value)) + geom_violin() + 
    ggtitle(colnames(cluster2)[i]) + theme_cowplot()
  plots_lineageSpecific_clus2[[i]] = plots_lineageSpecific_clus2[[i]]  + theme_ameres(type = 'barplot') + 
    ylab(label = 'Normalized lineage fraction') + xlab('Lineage')
  linePlot_clus2[[i]] = cluster2_tmp %>% dplyr::group_by(variable) %>% dplyr::summarise(median(value))
}

linePlot_clus2 = melt(linePlot_clus2)
ggpubr::ggline(linePlot_clus2,x='variable',y='value',color='L1',facet.by = 'L1') + ylim(c(0,1)) + theme_ameres(type = 'barplot')


linePlot_clus1 = linePlot_clus1 %>% dplyr::mutate(cluster = 1)
linePlot_clus2 = linePlot_clus2 %>% dplyr::mutate(cluster = 2)
linePlots = rbind.data.frame(linePlot_clus1,linePlot_clus2)
linePlots$cluster = factor(linePlots$cluster)
pdf("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure5/plots/clusters_lineageFractions_linePlots.pdf", height = 8, width=8)
  ggpubr::ggline(linePlots,x='variable',y='value',facet.by = 'L1',group = 'cluster',
               color = 'cluster',palette = c("red","black"),ylab='median lineage fraction') +
    theme_ameres(type='barplot') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
##### fishers test

pdf('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure5/plots/fractionUbiquitous_lineageSpecific.pdf',
    height=3, width=3)
clus0_lineage = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==1 & Z_MZshieldSpecific$class==0),])/nrow(Z_MZshieldSpecific[which( Z_MZshieldSpecific$class==0),])
clus0_ubiquitous = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==2 & Z_MZshieldSpecific$class==0),])/nrow(Z_MZshieldSpecific[which( Z_MZshieldSpecific$class==0),])
pie(c(clus0_lineage,clus0_ubiquitous),col = c('red','black'),main = 'Cluster Z',labels = c('lineage','ubiquitous'))

clus1_lineage = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==1 & Z_MZshieldSpecific$class==1),])/nrow(Z_MZshieldSpecific[which( Z_MZshieldSpecific$class==1),])
clus1_ubiquitous= nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==2 & Z_MZshieldSpecific$class==1),])/nrow(Z_MZshieldSpecific[which( Z_MZshieldSpecific$class==1),])
pie(c(clus1_lineage,clus1_ubiquitous),col = c('red','black'),main = 'Cluster 1',labels = c('lineage','ubiquitous'))

clus2_lineage = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==1 & Z_MZshieldSpecific$class==2),])/nrow(Z_MZshieldSpecific[which( Z_MZshieldSpecific$class==2),])
clus2_ubiquitous = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==2 & Z_MZshieldSpecific$class==2),])/nrow(Z_MZshieldSpecific[which( Z_MZshieldSpecific$class==2),])
pie(c(clus2_lineage,clus2_ubiquitous),col = c('red','black'),main = 'Cluster 2',labels = c('lineage','ubiquitous'))

clus3_lineage = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==1 & Z_MZshieldSpecific$class==3),])/nrow(Z_MZshieldSpecific[which( Z_MZshieldSpecific$class==3),])
clus3_ubiquitous = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==2 & Z_MZshieldSpecific$class==3),])/nrow(Z_MZshieldSpecific[which( Z_MZshieldSpecific$class==3),])
pie(c(clus3_lineage,clus3_ubiquitous),col = c('red','black'),main = 'Cluster 3',labels = c('lineage','ubiquitous'))

clus4_lineage = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==1 & Z_MZshieldSpecific$class==4),])/nrow(Z_MZshieldSpecific[which( Z_MZshieldSpecific$class==4),])
clus4_ubiquitous  = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==2 & Z_MZshieldSpecific$class==4),])/nrow(Z_MZshieldSpecific[which( Z_MZshieldSpecific$class==4),])
pie(c(clus4_lineage,clus4_ubiquitous),col = c('red','black'),main = 'Cluster 4',labels = c('lineage','ubiquitous'))

dev.off()

##### now creating a fishers exact test .. 

### in cluster 0 and lineage specific\
### in cluster 0 and not lineage specific
### not in cluster 0 and lineage specific
### not in cluster 0 and not lineage specific 
#### Zygotic 
pval_cluster_enrich = c()
pval_cluster_depleted = c()
##### 1
incluster_lineageSpecific = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==1 & Z_MZshieldSpecific$class==1),])
incluster_notLineageSpecific =  nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==2 & Z_MZshieldSpecific$class==1),])
notIncluster_lineageSpecific = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==1 & Z_MZshieldSpecific$class!=1),])
notIncluster_notlineageSpecific = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster!=1 & Z_MZshieldSpecific$class!=1),])
pval_cluster_enrich = c(pval_cluster_enrich,fisher.test(matrix(c(incluster_lineageSpecific,incluster_notLineageSpecific,notIncluster_lineageSpecific,notIncluster_notlineageSpecific),byrow = T,nrow = 2),alternative = 'greater')$p.value)
pval_cluster_depleted = c(pval_cluster_depleted,fisher.test(matrix(c(incluster_lineageSpecific,incluster_notLineageSpecific,notIncluster_lineageSpecific,notIncluster_notlineageSpecific),byrow = T,nrow = 2),alternative = 'less')$p.value)

### 2

incluster_lineageSpecific = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==1 & Z_MZshieldSpecific$class==2),])
incluster_notLineageSpecific =  nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==2 & Z_MZshieldSpecific$class==2),])
notIncluster_lineageSpecific = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==1 & Z_MZshieldSpecific$class!=2),])
notIncluster_notlineageSpecific = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster!=1 & Z_MZshieldSpecific$class!=2),])
pval_cluster_enrich = c(pval_cluster_enrich,fisher.test(matrix(c(incluster_lineageSpecific,incluster_notLineageSpecific,notIncluster_lineageSpecific,notIncluster_notlineageSpecific),byrow = T,nrow = 2),alternative = 'greater')$p.value)
pval_cluster_depleted = c(pval_cluster_depleted,fisher.test(matrix(c(incluster_lineageSpecific,incluster_notLineageSpecific,notIncluster_lineageSpecific,notIncluster_notlineageSpecific),byrow = T,nrow = 2),alternative = 'less')$p.value)

### 3 

incluster_lineageSpecific = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==1 & Z_MZshieldSpecific$class==3),])
incluster_notLineageSpecific =  nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==2 & Z_MZshieldSpecific$class==3),])
notIncluster_lineageSpecific = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==1 & Z_MZshieldSpecific$class!=3),])
notIncluster_notlineageSpecific = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster!=1 & Z_MZshieldSpecific$class!=3),])
pval_cluster_enrich =  c(pval_cluster_enrich,fisher.test(matrix(c(incluster_lineageSpecific,incluster_notLineageSpecific,notIncluster_lineageSpecific,notIncluster_notlineageSpecific),byrow = T,nrow = 2),alternative = 'greater')$p.value)
pval_cluster_depleted = c(pval_cluster_depleted,fisher.test(matrix(c(incluster_lineageSpecific,incluster_notLineageSpecific,notIncluster_lineageSpecific,notIncluster_notlineageSpecific),byrow = T,nrow = 2),alternative = 'less')$p.value)

### 4
incluster_lineageSpecific = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==1 & Z_MZshieldSpecific$class==4),])
incluster_notLineageSpecific =  nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==2 & Z_MZshieldSpecific$class==4),])
notIncluster_lineageSpecific = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==1 & Z_MZshieldSpecific$class!=4),])
notIncluster_notlineageSpecific = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster!=1 & Z_MZshieldSpecific$class!=4),])
pval_cluster_enrich = c(pval_cluster_enrich,fisher.test(matrix(c(incluster_lineageSpecific,incluster_notLineageSpecific,notIncluster_lineageSpecific,notIncluster_notlineageSpecific),byrow = T,nrow = 2),alternative = 'greater')$p.value)
pval_cluster_depleted = c(pval_cluster_depleted,fisher.test(matrix(c(incluster_lineageSpecific,incluster_notLineageSpecific,notIncluster_lineageSpecific,notIncluster_notlineageSpecific),byrow = T,nrow = 2),alternative = 'less')$p.value)


### zygotic 


incluster_lineageSpecific = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==1 & Z_MZshieldSpecific$class==0),])
incluster_notLineageSpecific =  nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==2 & Z_MZshieldSpecific$class==0),])
notIncluster_lineageSpecific = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster==1 & Z_MZshieldSpecific$class!=0),])
notIncluster_notlineageSpecific = nrow(Z_MZshieldSpecific[which(Z_MZshieldSpecific$cluster!=1 & Z_MZshieldSpecific$class!=0),])

pval_cluster_enrich = c(pval_cluster_enrich,fisher.test(matrix(c(incluster_lineageSpecific,incluster_notLineageSpecific,notIncluster_lineageSpecific,notIncluster_notlineageSpecific),byrow = T,nrow = 2),alternative = 'greater')$p.value)
pval_cluster_depleted = c(pval_cluster_depleted,fisher.test(matrix(c(incluster_lineageSpecific,incluster_notLineageSpecific,notIncluster_lineageSpecific,notIncluster_notlineageSpecific),byrow = T,nrow = 2),alternative = 'less')$p.value)


order_cols = c('possibleLineages_allCells_high.txt','possibleLineages_allCells_oblong.txt','possibleLineages_allCells_dome.txt',
               'possibleLineages_allCells_epiboly30.txt','possibleLineages_allCells_epiboly50.txt','possibleLineages_allCells_shield.txt',
               'possibleLineages_allCells_epiboly60.txt','possibleLineages_allCells_epiboly75.txt','possibleLineages_allCells_epiboly90.txt',
               'possibleLineages_allCells_bud.txt','possibleLineages_allCells_somite3.txt','possibleLineages_allCells_somite6.txt')


#### CLus4 

combineTimepoints  = function(clus4){
  
  EVL_clus4 = lapply(clus4 ,function(x) x$EVL)
  EVL_clus4 = as.data.frame(t(do.call(rbind,EVL_clus4)))
  EVL_clus4 = EVL_clus4[,order_cols]
  rownames(EVL_clus4) = rownames(clus4$possibleLineages_allCells_bud.txt)
  #pheatmap(EVL_clus4,cluster_cols = F,cluster_rows = T,fontsize_row = 1.5,main = ' EVL')
  
  marginalL_clus4 = lapply(clus4 ,function(x) x$Marginal)
  marginalL_clus4 = as.data.frame(t(do.call(rbind,marginalL_clus4)))
  marginalL_clus4 = marginalL_clus4[,order_cols]
  rownames(marginalL_clus4) = rownames(clus4$possibleLineages_allCells_bud.txt)
  #pheatmap(marginalL_clus4,cluster_cols = F,fontsize_row = 1.5,main = ' Marginal')
  
  marginalLDorsal_clus4 = lapply(clus4 ,function(x) x$MarginalDorsal)
  marginalLDorsal_clus4 = as.data.frame(t(do.call(rbind,marginalLDorsal_clus4)))
  marginalLDorsal_clus4 = marginalLDorsal_clus4[,order_cols]
  rownames(marginalLDorsal_clus4) = rownames(clus4$possibleLineages_allCells_bud.txt)
  #pheatmap(marginalLDorsal_clus4,cluster_cols = F,fontsize_row = 1.5,main = ' Marginal Dorsal')
  
  Ventral_clus4 = lapply(clus4 ,function(x) x$Ventral)
  Ventral_clus4 = as.data.frame(t(do.call(rbind,Ventral_clus4)))
  Ventral_clus4 = Ventral_clus4[,order_cols]
  rownames(Ventral_clus4) = rownames(clus4$possibleLineages_allCells_bud.txt)
  #pheatmap(Ventral_clus4,cluster_cols = F,fontsize_row = 1.5,main = ' Ventral')
  
  PGC_clus4 = lapply(clus4 ,function(x) x$PGC)
  PGC_clus4 = as.data.frame(t(do.call(rbind,PGC_clus4)))
  PGC_clus4 = PGC_clus4[,order_cols]
  rownames(PGC_clus4) = rownames(clus4$possibleLineages_allCells_bud.txt)
  #pheatmap(PGC_clus4,cluster_cols = F,fontsize_row = 1.5,main = ' PGC')
  
  DorsalAnimal_clus4 = lapply(clus4 ,function(x) x$DorsalAnimal)
  DorsalAnimal_clus4 = as.data.frame(t(do.call(rbind,DorsalAnimal_clus4)))
  DorsalAnimal_clus4 = DorsalAnimal_clus4[,order_cols]
  rownames(DorsalAnimal_clus4) = rownames(clus4$possibleLineages_allCells_bud.txt)
  #pheatmap(DorsalAnimal_clus4,cluster_cols = F,fontsize_row = 1.5,main = ' Dorsal animal')
  
  VentralAnimal_clus4 = lapply(clus4 ,function(x) x$VentralAnimal)
  VentralAnimal_clus4 = as.data.frame(t(do.call(rbind,VentralAnimal_clus4)))
  VentralAnimal_clus4 = VentralAnimal_clus4[,order_cols]
  rownames(VentralAnimal_clus4) = rownames(clus4$possibleLineages_allCells_bud.txt)
  #pheatmap(VentralAnimal_clus4,cluster_cols = F,fontsize_row = 1.5,main = ' Ventral animal')
  
  Dorsal_clus4 = lapply(clus4 ,function(x) x$Dorsal)
  Dorsal_clus4 = as.data.frame(t(do.call(rbind,Dorsal_clus4)))
  Dorsal_clus4 = Dorsal_clus4[,order_cols]
  rownames(Dorsal_clus4) = rownames(clus4$possibleLineages_allCells_bud.txt)
  #pheatmap(Dorsal_clus4,cluster_cols = F,fontsize_row = 1.5,main = ' Dorsal')
  
  #### putting all the lineages together.. 
  colnames(EVL_clus4) = paste0('EVL_',c(3.3,3.8,4.3,4.7,5.3,6,7,8,9,10,11,12))
  colnames(marginalL_clus4) = paste0('Marginal_',c(3.3,3.8,4.3,4.7,5.3,6,7,8,9,10,11,12))
  colnames(marginalLDorsal_clus4) = paste0('MarginalDorsal_',c(3.3,3.8,4.3,4.7,5.3,6,7,8,9,10,11,12))
  colnames(Ventral_clus4) = paste0('Ventral_', c(3.3,3.8,4.3,4.7,5.3,6,7,8,9,10,11,12))
  colnames(PGC_clus4) = paste0('PGC_',c(3.3,3.8,4.3,4.7,5.3,6,7,8,9,10,11,12))
  colnames(DorsalAnimal_clus4) = paste0('DorsalAnimal_',c(3.3,3.8,4.3,4.7,5.3,6,7,8,9,10,11,12))
  colnames(VentralAnimal_clus4) =paste0('VentralAnimal_',c(3.3,3.8,4.3,4.7,5.3,6,7,8,9,10,11,12))
  colnames(Dorsal_clus4) = paste0('Dorsal_',c(3.3,3.8,4.3,4.7,5.3,6,7,8,9,10,11,12))
  
  total_lineages = cbind.data.frame(EVL_clus4,marginalL_clus4,marginalLDorsal_clus4,Ventral_clus4, PGC_clus4,
                                     DorsalAnimal_clus4, VentralAnimal_clus4,Dorsal_clus4)
  total_lineages = total_lineages[complete.cases(total_lineages),]
  pheatmap(total_lineages,cluster_cols = F,fontsize_row = 1.5,main = ' Dorsal',border_color = NA)
  
  total_lineages_ = melt(total_lineages)
  total_lineages_$lineage = unlist(lapply(strsplit(as.character(total_lineages_$variable),"_",T),function(x) x[1]))
  total_lineages_$time = unlist(lapply(strsplit(as.character(total_lineages_$variable),"_",T),function(x) x[2]))
  
  total_lineages_plot  = ggpubr::ggboxplot(total_lineages_,x='time',y='value',facet.by = 'lineage')
  print(total_lineages_plot)
 
  }


pdf('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure5/plots/fractionCellsInLineage.pdf',width=8, height = 8)

        combineTimepoints(clus4 = clus1)
        combineTimepoints(clus4 = clus2)
        combineTimepoints(clus4 = clus3)
        combineTimepoints(clus4 = clus4)

dev.off()

       


bud_ref = getLineageInteractions(testData = finalTable_data$finalInfo_bud.txt,numberOfcells_test = numberOfcells$possibleLineages_allCells_bud.txt,combinedSegments_reference = combinedSegments_epiboly,genes_MZ = genes_MZ)

somite6_ref = getLineageInteractions(testData = finalTable_data$finalInfo_somite6.txt,numberOfcells_test = numberOfcells$possibleLineages_allCells_somite6.txt,combinedSegments_reference = combinedSegments_epiboly,genes_MZ = genes_MZ)
somite3_ref = getLineageInteractions(testData = finalTable_data$finalInfo_somite3.txt,numberOfcells_test = numberOfcells$possibleLineages_allCells_somite3.txt,combinedSegments_reference = combinedSegments_epiboly,genes_MZ = genes_MZ)
epiboly60_ref = getLineageInteractions(testData = finalTable_data$finalInfo_epiboly60.txt,numberOfcells_test = numberOfcells$possibleLineages_allCells_epiboly60.txt,combinedSegments_reference = combinedSegments_epiboly,genes_MZ = genes_MZ)
epiboly90_ref = getLineageInteractions(testData = finalTable_data$finalInfo_epiboly90.txt,numberOfcells_test = numberOfcells$possibleLineages_allCells_epiboly90.txt,combinedSegments_reference = combinedSegments_epiboly,genes_MZ = genes_MZ)
high_ref = getLineageInteractions(testData = finalTable_data$finalInfo_high.txt,numberOfcells_test = numberOfcells$possibleLineages_allCells_high.txt,combinedSegments_reference = combinedSegments_epiboly,genes_MZ = genes_MZ)
shield_ref = getLineageInteractions(testData = finalTable_data$finalInfo_shield.txt,numberOfcells_test = numberOfcells$possibleLineages_allCells_shield.txt,combinedSegments_reference = combinedSegments_epiboly,genes_MZ = genes_MZ)


pheatmap::pheatmap(bud_ref$`1`[,c(1:8)],cluster_cols = F,fontsize_row = 1.5)
pheatmap::pheatmap(somite6_ref$`4`[,c(1:8)],cluster_cols = F,fontsize_row = 1.5,color=redblue(100))
pheatmap::pheatmap(somite3_ref$`4`[complete.cases(somite3_ref$`4`[,c(1:8)]),],cluster_cols = F,fontsize_row = 1.5,color=redblue(100))
pheatmap::pheatmap(epiboly60_ref$`4`[complete.cases(epiboly60_ref$`4`[,c(1:8)]),],cluster_cols = F,fontsize_row = 1.5,color=redblue(100))
pheatmap::pheatmap(epiboly90_ref$`4`[complete.cases(epiboly90_ref$`4`[,c(1:8)]),],cluster_cols = F,fontsize_row = 1.5,color=redblue(100))
pheatmap::pheatmap(epiboly90_ref$`3`[complete.cases(epiboly90_ref$`3`[,c(1:8)]),],cluster_cols = F,fontsize_row = 1.5,color=redblue(100))
pheatmap::pheatmap(epiboly90_ref$`1`[complete.cases(epiboly90_ref$`1`[,c(1:8)]),],cluster_cols = F,fontsize_row = 1.5,color=redblue(100))
pheatmap::pheatmap(high_ref$`1`,cluster_cols = F,fontsize_row = 1.5,color=redblue(100))
pheatmap::pheatmap(shield_ref$`1`,cluster_cols = F,fontsize_row = 1.5,color=redblue(100))

pdf('/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/MZdynamics/singlecellexpression/MZgenes/pheatmap_lineageFC.pdf', height = 10,width=5)
NMF::aheatmap(lineages_split$`1`[,c(1:8)],Rowv = NA,Colv = NA,color = c('white','grey','black'),border_color = NA)
    pheatmap::pheatmap(lineages_split$`1`[,c(1:8)],cluster_cols = F,fontsize_row = 1.5,color=redblue(100))
     pheatmap::pheatmap(lineages_split$`4`[,c(1:8)],cluster_cols = F ,fontsize_row = 1.5,color=redblue(100),Rowv = NA,Colv = NA)
dev.off()



########################
##### same for non ubiquitous genes
####### ubiquitous
genes_MZ_NONubiquitous = genes_MZ %>% dplyr::filter(ubiquitousness == F) %>% dplyr::filter(signalPerExpressedCell_high>5)

write.table(genes_MZ_NONubiquitous,'//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/MZdynamics/nonUbiquitousMZgenes.txt',sep="\t",quote = F)


genes_MZ_split= split(genes_MZ_NONubiquitous,genes_MZ_NONubiquitous$V11,T)


pdf('/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/MZdynamics/singlecellexpression/MZgenes/nonUbiquitousGenes_fractionOfcells.pdf', height = 4)
          genes_MZ_NONubiquitous_quantSeq = lapply(genes_MZ_split,function(x) x %>% dplyr::filter(ubiquitousness == F)%>% dplyr::select(dplyr::contains("Normalized")) )
          genes_MZ_NONubiquitous_quantSeq = lapply(genes_MZ_NONubiquitous_quantSeq, setNames, c(0.75,2,2.5,3,3.5,4,4.5,5,5.5))
          genes_MZ_NONubiquitous_quantSeq = reshape::melt(genes_MZ_NONubiquitous_quantSeq)
          ggpubr::ggboxplot(genes_MZ_NONubiquitous_quantSeq,x="variable",y='value',group="variable" ,ylab= "Normalized RPM expression",title = nrow(genes_MZ_NONubiquitous_quantSeq)/9,facet.by = 'L1', fill='L1',palette = 'Set1')
          
          genes_MZ_nonubiquitous_fractionOfCells = lapply(genes_MZ_split,function(x) x %>% dplyr::filter(ubiquitousness == F)%>% dplyr::select(dplyr::contains("fractionOfcells")) )
          genes_MZ_nonubiquitous_fractionOfCells = lapply(genes_MZ_nonubiquitous_fractionOfCells, setNames, c(3.3,3.8,4.3,4.7,5.3,6))
          genes_MZ_nonubiquitous_fractionOfCells = reshape::melt(genes_MZ_nonubiquitous_fractionOfCells)
          ggpubr::ggboxplot(genes_MZ_nonubiquitous_fractionOfCells,x="variable",y='value',group="variable",ylab = "Fraction of expressed cells",title = nrow(genes_MZ_nonubiquitous_fractionOfCells)/6,facet.by = 'L1', fill='L1',palette = 'Set1')
          
          
          genes_MZ_NONubiquitous_singlecell = lapply(genes_MZ_split,function(x) as.data.frame(t(apply(x %>% dplyr::filter(ubiquitousness==F)%>% dplyr::select(dplyr::contains("ssignalPerCell")) ,1,function(y) y/max(y)))))
          genes_MZ_NONubiquitous_singlecell = lapply(genes_MZ_NONubiquitous_singlecell, setNames, c(3.3,3.8,4.3,4.7,5.3,6))
          genes_MZ_NONubiquitous_singlecell = reshape::melt(genes_MZ_NONubiquitous_singlecell)
          ggpubr::ggboxplot(genes_MZ_NONubiquitous_singlecell,x="variable",y='value',group="variable",ylab = "Signal per cel",title = nrow(genes_MZ_nonubiquitous_fractionOfCells)/6,facet.by = 'L1', fill='L1',palette = 'Set1')
          
          genes_MZ_NONubiquitous_Perexpressedsinglecell = lapply(genes_MZ_split,function(x) as.data.frame(t(apply(x %>% dplyr::filter(ubiquitousness==F)%>% dplyr::select(dplyr::contains("signalPerExpressedCell")) ,1,function(y) y/max(y)))))
          genes_MZ_NONubiquitous_Perexpressedsinglecell = lapply(genes_MZ_NONubiquitous_Perexpressedsinglecell, setNames, c(3.3,3.8,4.3,4.7,5.3,6))
          genes_MZ_NONubiquitous_Perexpressedsinglecell = reshape::melt(genes_MZ_NONubiquitous_Perexpressedsinglecell)
          ggpubr::ggboxplot(genes_MZ_NONubiquitous_Perexpressedsinglecell,x="variable",y='value',group="variable",ylab = "Signal per expressed cell",title = nrow(genes_MZ_nonubiquitous_fractionOfCells)/6,facet.by = 'L1',ylim = c(0,1), fill='L1',palette = 'Set1')

dev.off()

genes_MZ_split = lapply(genes_MZ_split,function(x) plyr::join(x,combinedSegments_epiboly))


lineages_split = lapply(genes_MZ_split,function(x) x %>%dplyr::select(EVL:Dorsal) )
#lineages_split = lapply(lineages_split,function(x) x %>% select(-'PGC'))
rownames(lineages_split$fast) = genes_MZ_split$fast$genes
rownames(lineages_split$slow) = genes_MZ_split$slow$genes
pdf('/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/plots/MZdynamics/singlecellexpression/MZgenes/pheatmap_NoNUBIQUITOUSlineageFC.pdf', height = 10,width=5)
  pheatmap::pheatmap(lineages_split$`4`[,c(1:8)],cluster_cols = F,fontsize_row = 1.5)
  a =  pheatmap::pheatmap(lineages_split$slow[,c(1:8)],cluster_cols = F ,fontsize_row = 1.5)
dev.off()



length(which(apply(lineages_split$slow[,c(1:8)], 1,function(x) length(which(x <0)))==8))
length(which(apply(lineages_split$fast[,c(1:8)], 1,function(x) length(which(x <0)))==8))

lineages_split$fast$genes = rownames(lineages_split$fast) 
lineages_split$slow$genes = rownames(lineages_split$slow)
lineages_split$fast$type = 'fast'
lineages_split$slow$type = 'slow'
totalLineages = rbind.data.frame(lineages_split$fast,lineages_split$slow)
write.table(totalLineages,'/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/MZdynamics/singleCell/ZandMZanalysis/MZ_NONuniquitousGenes.txt',sep="\t",quote = F,row.names = F)


############ now checking of these are lineage specific... 



####### plotting the gene expression for ubiquitous and non ubiquitous genes
# genes_MZ_NONubiquitous_quantSeq = genes_MZ %>% dplyr::filter(ubiquitousness == F)  %>% dplyr::select(dplyr::contains("Normalized"))
# colnames(genes_MZ_NONubiquitous_quantSeq) = c(0.75,2,2.5,3,3.5,4,4.5,5,5.5)
# genes_MZ_NONubiquitous_quantSeq = reshape::melt(genes_MZ_NONubiquitous_quantSeq)
# ggpubr::ggboxplot(genes_MZ_NONubiquitous_quantSeq,x="variable",y='value',group="variable",ylab= "Normalized RPM expression",title = nrow(genes_MZ_NONubiquitous_quantSeq)/9)
# 
# genes_MZ_NONubiquitous_fractionOfCells = genes_MZ %>% dplyr::filter(ubiquitousness == F) %>% dplyr::select(dplyr::contains('fractionOfcells'))
# colnames(genes_MZ_NONubiquitous_fractionOfCells) = c(3.3,3.8,4.3,4.7,5.3,6)
# genes_MZ_NONubiquitous_fractionOfCells = reshape::melt(genes_MZ_NONubiquitous_fractionOfCells)
# ggpubr::ggboxplot(genes_MZ_NONubiquitous_fractionOfCells,x="variable",y='value',group="variable",ylab= "fraction of cells expressing",title = nrow(genes_MZ_NONubiquitous_quantSeq)/9)
# 
# 
# genes_MZ_NONubiquitous_singlecell = t(apply(genes_MZ %>% dplyr::filter(ubiquitousness == F) %>% dplyr::select(dplyr::contains("ssignalPerCell")),1,function(x) x/max(x)))
# colnames(genes_MZ_NONubiquitous_singlecell) = c(3.3,3.8,4.3,4.7,5.3,6)
# genes_MZ_NONubiquitous_singlecell = reshape::melt(genes_MZ_NONubiquitous_singlecell)
# ggpubr::ggboxplot(genes_MZ_NONubiquitous_singlecell,x="X2",y='value',group="X2",ylab= "Signal per cell",title = nrow(genes_MZ_NONubiquitous_quantSeq)/9)
# 
# genes_MZ_NONubiquitous_Perexpressedsinglecell = t(apply(genes_MZ %>% dplyr::filter(ubiquitousness == F) %>% dplyr::select(dplyr::contains("signalPerExpressedCell")),1,function(x) x/max(x)))
# colnames(genes_MZ_NONubiquitous_Perexpressedsinglecell) = c(3.3,3.8,4.3,4.7,5.3,6)
# genes_MZ_NONubiquitous_Perexpressedsinglecell = reshape::melt(genes_MZ_NONubiquitous_Perexpressedsinglecell)
# ggpubr::ggboxplot(genes_MZ_NONubiquitous_Perexpressedsinglecell,x="X2",y='value',group="X2",ylab= "Signal per expressing cell",title = nrow(genes_MZ_NONubiquitous_quantSeq)/9)
# 
##########


