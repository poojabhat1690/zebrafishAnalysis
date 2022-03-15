### creating a master table with RPMS of poly(A), rRNA depleted libraries, 3' end sequencing and poly(A) tail lengths. 
  ##### RNAseq data summary #### these data are still un-published #################################
                ###################### 3 replicates of this available for poly(A) selected RNAseq and rRNA depleted RNAseq
                    #t1 = 4-cell embryos (~1 hour post fertilization (hpf))
                    #t2 = 64-128 cell (~2.5 hpf)
                    #t3 = 256-512 cell (~3.5 hpf)
                    #t4 = oblong/sphere (~4.5 hpf)
               ################################################

  
  ######## mean tail lengths downloaded from Subtelney et al., 2014  at 2,4, 6 hpf ###################################
  
  
  
##########  defining functions used in this script #######################

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
      number_ticks <- function(n) {function(limits) pretty(limits, n)}
      
##########################################################################

################# loading libraries #######################################
library(dplyr)
library(patchwork)
library(ggplot2)
###########################################################################
  #### OUTDIR 
      dir.create('/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/plots/polyAtailLengthBias/')
      outdir = "/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/figures/figure1/plots/polyAtailLengthBias/"
      #################### Processing RNAseq datasets #################################
            TPM_RNAseq = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/packages/data/figure2/TPM_STAR_dr11.txt",sep="\t",stringsAsFactors = F, header = T)
      
                #TPM_RNAseq = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/externalDatat/fromAndi_ribo0VsPolyA/TPM_STAR_dr11.txt",sep="\t",stringsAsFactors = F, header = T)
                TPM_RNAseq$ensembl_gene_id = row.names(TPM_RNAseq)
                
                dr11=read.table("/Volumes/Macintosh HD/Users/pooja.bhat//Dropbox (VBC)/Paperoutline/figureOutline/packages/data/figure2/allGenes_dr11.txt",sep="\t", header=T)
                #dr11 = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/dr11/allGenes_dr11.txt",sep="\t", header = T)
                dr11 = dr11 %>% dplyr::select(c('ensembl_gene_id','external_gene_name'))
                
                T1_Ribocop = TPM_RNAseq %>% dplyr::select(dplyr::contains("T1.RiboCop")) %>% dplyr::mutate(T1_riboCop = rowMeans(.)) %>% dplyr::select(T1_riboCop)
                T2_Ribocop = TPM_RNAseq %>% dplyr::select(dplyr::contains("T2.RiboCop")) %>% dplyr::mutate(T2_riboCop = rowMeans(.))%>% dplyr::select(T2_riboCop)
                T3_Ribocop = TPM_RNAseq %>% dplyr::select(dplyr::contains("T3.RiboCop")) %>% dplyr::mutate(T3_riboCop = rowMeans(.))%>% dplyr::select(T3_riboCop)
                T4_Ribocop = TPM_RNAseq %>% dplyr::select(dplyr::contains("T4.RiboCop")) %>% dplyr::mutate(T4_riboCop = rowMeans(.))%>% dplyr::select(T4_riboCop)
                
                
                T1_polyA = TPM_RNAseq %>% dplyr::select(dplyr::contains("T1.polyA")) %>% dplyr::mutate(T1_polyA = rowMeans(.)) %>% dplyr::select(T1_polyA)
                T2_polyA = TPM_RNAseq %>% dplyr::select(dplyr::contains("T2.polyA")) %>% dplyr::mutate(T2_polyA = rowMeans(.))%>% dplyr::select(T2_polyA)
                T3_polyA = TPM_RNAseq %>% dplyr::select(dplyr::contains("T3.polyA")) %>% dplyr::mutate(T3_polyA = rowMeans(.))%>% dplyr::select(T3_polyA)
                T4_polyA = TPM_RNAseq %>% dplyr::select(dplyr::contains("T4.polyA")) %>% dplyr::mutate(T4_polyA = rowMeans(.))%>% dplyr::select(T4_polyA)
                
                TPM_RNAseq_means = data.frame(T1_Ribocop = T1_Ribocop,T2_Ribocop = T2_Ribocop,T3_Ribocop = T3_Ribocop,T4_Ribocop = T4_Ribocop,T1_polyA = T1_polyA,T2_polyA = T2_polyA,T3_polyA = T3_polyA,T4_polyA = T4_polyA, ensembl_gene_id = TPM_RNAseq$ensembl_gene_id)
                TPM_RNAseq = cbind.data.frame(TPM_RNAseq,TPM_RNAseq_means)
                
                TPM_RNAseq = plyr::join(TPM_RNAseq, dr11)
                TPM_RNAseq$name = TPM_RNAseq$external_gene_name
                TPM_RNAseq = TPM_RNAseq[!duplicated(TPM_RNAseq$name),]
                
        ########################################################################################
                #### 3' end sequencing RPMS ##################################
                QuantSeqReads = read.table("/Volumes/Macintosh HD//Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure2/data/RPM_allCws.txt",
                           sep="\t",stringsAsFactors = F, header = T)
                # QuantSeqReads = read.table("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/slamDunk_combinedStageSpecific_dr11_may2019/dataTables_expandedCounting/RPM_allCws.txt",
                #                            sep="\t",stringsAsFactors = F, header = T)
                QuantSeqReads = QuantSeqReads %>% dplyr::group_by(name)  %>%   dplyr::summarise_at(vars(Inj_R2_TP1:Untreated_TP9), sum, na.rm = TRUE)
                
                QuantSeqReads_split = splitReplicates(dataFrameToSplit = QuantSeqReads,condition = "Inj",metadata_add = QuantSeqReads$name)
                QuantSeqReads_means = QuantSeqReads_split[[3]]
                colnames(QuantSeqReads_means) = c("Inj_mean_TP1","Inj_mean_TP2","Inj_mean_TP3","Inj_mean_TP4",
                                                  "Inj_mean_TP5","Inj_mean_TP6","Inj_mean_TP7","Inj_mean_TP8","Inj_mean_TP9", "name")
                quantSeq_R1 = QuantSeqReads_split[[1]]
                quantSeq_R2 = QuantSeqReads_split[[2]]
                quantSeq_R1 = quantSeq_R1 %>% dplyr::select(dplyr::contains('Inj'))
                quantSeq_R2 = quantSeq_R2 %>% dplyr::select(dplyr::contains('Inj'))
                
                colnames(quantSeq_R1) =  c("Inj_R1_TP1","Inj_R1_TP2","Inj_R1_TP3","Inj_R1_TP4",
                                           "Inj_R1_TP5","Inj_R1_TP6","Inj_R1_TP7","Inj_R1_TP8","Inj_R1_TP9")
                
                colnames(quantSeq_R2) =  c("Inj_R2_TP1","Inj_R2_TP2","Inj_R2_TP3","Inj_R2_TP4",
                                           "Inj_R2_TP5","Inj_R2_TP6","Inj_R2_TP7","Inj_R2_TP8","Inj_R2_TP9")
                
                totalQuantSeq  = cbind.data.frame(quantSeq_R1,quantSeq_R2,QuantSeqReads_means)
                allGenes_quantSeq = data.frame(external_gene_name = totalQuantSeq[,'name'])
          ###### combining the RNAseq and QuantSeq datasets
                
                RNAseq_QuantSeq = plyr::join(totalQuantSeq,TPM_RNAseq)
                RNAseq_QuantSeq = RNAseq_QuantSeq[!duplicated(RNAseq_QuantSeq),]
                RNAseq_QuantSeq = RNAseq_QuantSeq[complete.cases(RNAseq_QuantSeq$external_gene_name),]
        
          ############################# also addding to this poly(A) tail lengths from Subtelney et al., 2014 and Chang et al, 2018 ##################
                
                samplesPolyAtail = list.files("/Volumes/Macintosh HD/Users/pooja.bhat/Downloads//polyAtailLengths_SutbtenlyEtal",pattern = "*.txt")
                
                #samplesPolyAtail = list.files("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/externalDatat/polyAtailLengths_SutbtenlyEtal",pattern = "*.txt")
                samplesPolyAtail = samplesPolyAtail[-grep(".gz",samplesPolyAtail)]
                samplesPolyAtail = samplesPolyAtail[grep("mock",samplesPolyAtail)]
                path_samplesPolyAtail = paste0("/Volumes/Macintosh HD/Users/pooja.bhat/Downloads//polyAtailLengths_SutbtenlyEtal/",samplesPolyAtail)
                
                #path_samplesPolyAtail = paste0("/Volumes/groups/ameres/Pooja/Projects/zebrafishAnnotation/zebrafish_analysis/importantDataframes/externalDatat/polyAtailLengths_SutbtenlyEtal/",samplesPolyAtail)
                
                
                polyAtailData = lapply(path_samplesPolyAtail,function(x) read.delim(x))
                names(polyAtailData) = c("2hpf","4hpf","6hpf")
                polyAtailData_meanTailLength = lapply(polyAtailData,function(x) cbind.data.frame(x$Mean.TL,x$Gene.name))
                polyAtailData_meanTailLength = lapply(polyAtailData_meanTailLength,function(x) x[-which(x$`x$Mean.TL`<=0),])
                colnames(polyAtailData_meanTailLength$`2hpf`) = c("meanTail","name")
                colnames(polyAtailData_meanTailLength$`4hpf`) = c("meanTail","name")
                colnames(polyAtailData_meanTailLength$`6hpf`) = c("meanTail","name")
                
                
                polyAtailData_meanTailLength = lapply(polyAtailData_meanTailLength,function(x) x %>% 
                                                        dplyr::mutate(quantile=dplyr::ntile(meanTail, n=20)))
                
                polyAtailData_meanTailLength = lapply(polyAtailData_meanTailLength,function(x) x %>% 
                                                        dplyr::group_by(quantile) %>% dplyr::mutate(range = paste0(min(meanTail),"-",max(meanTail))))
                
                for(i in 1:length(polyAtailData_meanTailLength)){
                  
                  polyAtailLength = as.data.frame(polyAtailData_meanTailLength[[i]])
                  colnames(polyAtailLength) = c(paste0('meanTail_',names(polyAtailData_meanTailLength)[i]),'name','quantile',paste0('range_',names(polyAtailData_meanTailLength)[i]))
                  polyAtailLength = polyAtailLength %>% dplyr::select(-quantile)
                  RNAseq_QuantSeq = plyr::join(RNAseq_QuantSeq,polyAtailLength,by='name')
                }
                
                
        ##### Also adding here the poly(A) tail lengths from Nary Kim's lab. 
                
                # tailLength_2hpf = read.table('/Volumes//groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/MZdynamics/uridylationData/polyAtailLengths_allhpf_2.txt',
                #                              sep="\t",stringsAsFactors = F, header = T)
                tailLength_2hpf = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Downloads/polyAtailLengths_allhpf_2.txt",sep="\t",stringsAsFactors = F, header = T)
                tailLength_2hpf = tailLength_2hpf %>% dplyr::group_by(external_gene_name) %>% 
                  dplyr::summarise(meanTail=mean(Atail,na.rm=T))
               
                 tailLength_2hpf = tailLength_2hpf[!duplicated(tailLength_2hpf),] %>% dplyr::mutate(kim_2hpf =meanTail ) %>% 
                  dplyr::select(kim_2hpf, external_gene_name) 
                tailLength_2hpf = tailLength_2hpf[complete.cases(tailLength_2hpf),]
                tailLength_2hpf = tailLength_2hpf %>% dplyr::mutate(quantile_kim = dplyr::ntile(kim_2hpf,20)) %>% dplyr::group_by(quantile_kim) %>%
                  dplyr::mutate(range_kim= paste0(round(min(kim_2hpf),2),"-",round(max(kim_2hpf),2))) %>% dplyr::ungroup()
                
                
                RNAseq_QuantSeq =plyr::join(RNAseq_QuantSeq,tailLength_2hpf)
                ###################################################### 
              
                RNAseq_QuantSeq = RNAseq_QuantSeq %>% dplyr::mutate(quantSeq_rRNA_2hpf = log2((Inj_mean_TP3)/(T2_riboCop)),polyA_rRNA_2hpf = log2((T2_polyA)/(T2_riboCop) ) )

                allGenes_quantSeq = plyr::join(allGenes_quantSeq,RNAseq_QuantSeq)
                allGenes_quantSeq = allGenes_quantSeq[!duplicated(allGenes_quantSeq),]
                allGenes_quantSeq = allGenes_quantSeq[!is.na(allGenes_quantSeq$Inj_R1_TP1),]
                
                write.table(allGenes_quantSeq,'/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)//Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures//figure1/data/masterTable_quantSeqRNAseqRPMs.txt',sep="\t",quote = F,row.names = F)
                
                quantSeq_polyA_rnaseq2hpf = RNAseq_QuantSeq %>% dplyr::select(range_2hpf,quantSeq_rRNA_2hpf,meanTail_2hpf,polyA_rRNA_2hpf,kim_2hpf,range_kim,quantile_kim)
                
                quantSeq_polyA_rnaseq2hpf$range_2hpf = as.factor(quantSeq_polyA_rnaseq2hpf$range_2hpf)
                
                ##### getting the p-values for subtelney et al., 2014
                  
                quantSeq_polyA_rnaseq2hpf_split = split(quantSeq_polyA_rnaseq2hpf,quantSeq_polyA_rnaseq2hpf$range_2hpf,T)
                
                pval_quantSeqRNAseq = data.frame(pVal_quantSeq = p.adjust(unlist(lapply(quantSeq_polyA_rnaseq2hpf_split,function(x) wilcox.test(x$quantSeq_rRNA_2hpf,quantSeq_polyA_rnaseq2hpf$quantSeq_rRNA_2hpf,paired = F)$p.val))))
                pval_quantSeqRNAseq = tibble::rownames_to_column(pval_quantSeqRNAseq, "range_2hpf")
                quantSeq_polyA_rnaseq2hpf =  plyr::join(quantSeq_polyA_rnaseq2hpf,pval_quantSeqRNAseq,by='range_2hpf') %>% dplyr::mutate(pval_sig_quantSeq = ifelse(pVal_quantSeq<0.01,'significant','non-significant'))
                
                
                pval_polyASeqRNAseq = data.frame(pVal_polyASeq = p.adjust(unlist(lapply(quantSeq_polyA_rnaseq2hpf_split,function(x) wilcox.test(x$polyA_rRNA_2hpf,quantSeq_polyA_rnaseq2hpf$polyA_rRNA_2hpf,paired = F)$p.val))))
                pval_polyASeqRNAseq = tibble::rownames_to_column(pval_polyASeqRNAseq, "range_2hpf")
                quantSeq_polyA_rnaseq2hpf =  plyr::join(quantSeq_polyA_rnaseq2hpf,pval_polyASeqRNAseq,by='range_2hpf') %>% dplyr::mutate(pval_sig_polyASeq = ifelse(pVal_polyASeq<0.01,'significant','non-significant'))
                
                #### getting p-values for ranges of Nary kim lab
                
                quantSeq_polyA_rnaseq2hpf_split = split(quantSeq_polyA_rnaseq2hpf,quantSeq_polyA_rnaseq2hpf$range_kim,T)
                pval_quantSeqRNAseq_kim = data.frame(pVal_quantSeq_kim = p.adjust(unlist(lapply(quantSeq_polyA_rnaseq2hpf_split,function(x) wilcox.test(x$quantSeq_rRNA_2hpf,quantSeq_polyA_rnaseq2hpf$quantSeq_rRNA_2hpf,paired = F)$p.val))))
                pval_quantSeqRNAseq_kim = tibble::rownames_to_column(pval_quantSeqRNAseq_kim, "range_kim")
                quantSeq_polyA_rnaseq2hpf =  plyr::join(quantSeq_polyA_rnaseq2hpf,pval_quantSeqRNAseq_kim,by='range_kim') %>% dplyr::mutate(pval_sig_quantSeq_kim = ifelse(pVal_quantSeq_kim<0.01,'significant','non-significant'))
                
                pval_polyASeqRNAseq_kim = data.frame(pVal_polyASeq_kim = p.adjust(unlist(lapply(quantSeq_polyA_rnaseq2hpf_split,function(x) wilcox.test(x$polyA_rRNA_2hpf,quantSeq_polyA_rnaseq2hpf$polyA_rRNA_2hpf,paired = F)$p.val))))
                pval_polyASeqRNAseq_kim = tibble::rownames_to_column(pval_polyASeqRNAseq_kim, "range_kim")
                quantSeq_polyA_rnaseq2hpf =  plyr::join(quantSeq_polyA_rnaseq2hpf,pval_polyASeqRNAseq_kim,by='range_kim') %>% dplyr::mutate(pval_sig_polyASeq_kim = ifelse(pVal_polyASeq_kim<0.01,'significant','non-significant'))
                
                
                quantSeq_polyA_rnaseq2hpf = quantSeq_polyA_rnaseq2hpf[is.finite(quantSeq_polyA_rnaseq2hpf$quantSeq_rRNA_2hpf) , ]
                quantSeq_polyA_rnaseq2hpf = quantSeq_polyA_rnaseq2hpf[is.finite(quantSeq_polyA_rnaseq2hpf$polyA_rRNA_2hpf) , ]
             
                   quantSeq_polyA_rnaseq2hpf = quantSeq_polyA_rnaseq2hpf[complete.cases(quantSeq_polyA_rnaseq2hpf),]
                  
                
                pdf(paste0(outdir,'/polyA_ribo0_quantSeq_subtelney.pdf'), height=5, width=7)
                
                           
                             quantSeq_RNAseq_boxplot =  ggpubr::ggboxplot(quantSeq_polyA_rnaseq2hpf,x='range_2hpf',y='quantSeq_rRNA_2hpf',fill ='pval_sig_quantSeq',
                                                                          bxp.errorbar = TRUE, outlier.shape = NA, palette = 'Set1',ylim=c(-10,5),title = paste0('n=',nrow(quantSeq_polyA_rnaseq2hpf))) + theme_ameres(type = 'barplot')
                              quantSeq_RNAseq_boxplot  = quantSeq_RNAseq_boxplot + theme(axis.line.x=element_blank(),axis.ticks.x=element_blank()) + theme(axis.text.x = element_text(angle = 90),
                                                                                                                                                           axis.title.x = element_blank(),legend.position = "none") +
                                ylab("3' end seq/rRNA depleted (log2)") + xlab('Poly(A) tail length')
                             
                              print(quantSeq_RNAseq_boxplot)
                              
                               
                              polyA_RNAseq_boxplot = ggpubr::ggboxplot(quantSeq_polyA_rnaseq2hpf,x='range_2hpf',y='polyA_rRNA_2hpf',fill='pval_sig_polyASeq',
                                                                       bxp.errorbar = TRUE, outlier.shape = NA, palette = 'Set1',ylim=c(-10,5),title = paste0('n=',nrow(quantSeq_polyA_rnaseq2hpf))) + theme_ameres(type = 'barplot') 
                              
                              polyA_RNAseq_boxplot  = polyA_RNAseq_boxplot + theme(axis.line.x=element_blank(),axis.ticks.x=element_blank()) + theme(axis.text.x = element_text(angle = 90),
                              axis.title.x = element_blank(),legend.position = "none") +
                                ylab("poly(A) RNAseq/rRNA depleted (log2)") + xlab('Poly(A) tail length')
                              print(polyA_RNAseq_boxplot)
                              
                              
                            
              

                dev.off()
                
                pdf(paste0(outdir,'//polyA_ribo0_quantSeq_subtelney_scatterPlot.pdf'), height=4, width=7)
                
                
                            quantSeq_RNAseq_scatterplot = ggpubr::ggscatter(quantSeq_polyA_rnaseq2hpf[is.finite(quantSeq_polyA_rnaseq2hpf$quantSeq_rRNA_2hpf) , ],
                                                                            x='meanTail_2hpf',y='quantSeq_rRNA_2hpf', alpha=0.1,ylim=c(-10,5),title = paste0('n=',nrow(quantSeq_polyA_rnaseq2hpf[is.finite(quantSeq_polyA_rnaseq2hpf$quantSeq_rRNA_2hpf) , ]))) + ggplot2::geom_hline(yintercept = 0, linetype="dashed",color='red')+ 
                              theme_ameres(type = 'barplot') +  ylab("3' end seq/rRNA depleted (log2)") + xlab('Poly(A) tail length')
                            
                            print(quantSeq_RNAseq_scatterplot) + scale_x_continuous(breaks=number_ticks(10))
                            
                            ggplot(quantSeq_polyA_rnaseq2hpf, aes(x=meanTail_2hpf,y=quantSeq_rRNA_2hpf)) + geom_density2d() + 
                              ggtitle( paste0('n=',nrow(quantSeq_polyA_rnaseq2hpf[is.finite(quantSeq_polyA_rnaseq2hpf$quantSeq_rRNA_2hpf) , ]))) +
                              ggplot2::geom_hline(yintercept = 0, linetype="dashed",color='red')+ theme_ameres(type = 'barplot')+
                              ylab("3' end seq/rRNA depleted (log2)") + xlab('Poly(A) tail length')+ scale_x_continuous(breaks=number_ticks(10))+ ylim(c(-3,3))
                            
                              
                          
                            polyA_RNAseq_scatterplot =  ggpubr::ggscatter(quantSeq_polyA_rnaseq2hpf[is.finite(quantSeq_polyA_rnaseq2hpf$polyA_rRNA_2hpf) , ],
                                                                          x='meanTail_2hpf',y='polyA_rRNA_2hpf', alpha=0.1,ylim=c(-10,5),title = paste0('n=',nrow(quantSeq_polyA_rnaseq2hpf[is.finite(quantSeq_polyA_rnaseq2hpf$polyA_rRNA_2hpf) , ])))+ ggplot2::geom_hline(yintercept = 0, linetype="dashed",color='red')+ 
                              theme_ameres(type = 'barplot') +
                              ylab("poly(A) RNAseq/rRNA depleted (log2)") + xlab('Poly(A) tail length') + scale_x_continuous(breaks=number_ticks(10))
                            print(polyA_RNAseq_scatterplot)
                            
                            
                            ggplot(quantSeq_polyA_rnaseq2hpf, aes(x=meanTail_2hpf,y=polyA_rRNA_2hpf)) + geom_density2d() + 
                              ggtitle( paste0('n=',nrow(quantSeq_polyA_rnaseq2hpf[is.finite(quantSeq_polyA_rnaseq2hpf$quantSeq_rRNA_2hpf) , ]))) +
                              ggplot2::geom_hline(yintercept = 0, linetype="dashed",color='red')+ theme_ameres(type = 'barplot')+
                              ylab("poly(A) RNAseq/rRNA depleted (log2)") + xlab('Poly(A) tail length')+ scale_x_continuous(breaks=number_ticks(10)) + ylim(c(-3,3))
                            
                            
                            
                dev.off()
                
                
                
                ### kim lab polyA
              
                pdf(paste0(outdir,'/polyA_ribo0_quantSeq_kim.pdf'), height=5, width=7)
                
                        quantSeq_polyA_rnaseq2hpf$quantile_kim = as.numeric(quantSeq_polyA_rnaseq2hpf$quantile_kim)
                        quantSeq_polyA_rnaseq2hpf = quantSeq_polyA_rnaseq2hpf[order(quantSeq_polyA_rnaseq2hpf$kim_2hpf,decreasing = F),]
                        quantSeq_RNAseq_boxplot =  ggpubr::ggboxplot(quantSeq_polyA_rnaseq2hpf,x='range_kim',y='quantSeq_rRNA_2hpf',fill ='pval_sig_quantSeq_kim',
                                                                     bxp.errorbar = TRUE, outlier.shape = NA, palette = 'Set1',ylim=c(-10,5),title = paste0('n=',nrow(quantSeq_polyA_rnaseq2hpf))) + theme_ameres(type = 'barplot')
                        quantSeq_RNAseq_boxplot  = quantSeq_RNAseq_boxplot + theme(axis.line.x=element_blank(),axis.ticks.x=element_blank()) + theme(axis.text.x = element_text(angle = 90),
                                                                                                                                                     axis.title.x = element_blank(),legend.position = "none") +
                          ylab("3' end seq/rRNA depleted (log2)") + xlab('Poly(A) tail length')
                        
                          print(quantSeq_RNAseq_boxplot)
                        
                        
                        polyA_RNAseq_boxplot = ggpubr::ggboxplot(quantSeq_polyA_rnaseq2hpf,x='range_kim',y='polyA_rRNA_2hpf',fill='pval_sig_polyASeq_kim',
                                                                 bxp.errorbar = TRUE, outlier.shape = NA, palette = 'Set1',ylim=c(-10,5),title = paste0('n=',nrow(quantSeq_polyA_rnaseq2hpf))) + theme_ameres(type = 'barplot') 
                        polyA_RNAseq_boxplot  = polyA_RNAseq_boxplot + theme(axis.line.x=element_blank(),axis.ticks.x=element_blank()) + theme(axis.text.x = element_text(angle = 90),
                                                                                                                                               axis.title.x = element_blank(),legend.position = "none") +
                              ylab("poly(A) RNAseq/rRNA depleted (log2)") + xlab('Poly(A) tail length')
                        print(polyA_RNAseq_boxplot)
                        
                        
                dev.off()
                
                #### scatter plot kim lab
                pdf(paste0(outdir,'/polyA_ribo0_quantSeq_kim_scatterplot.pdf'), height=4, width=7)      
                        quantSeq_RNAseq_scatterplot = ggpubr::ggscatter(quantSeq_polyA_rnaseq2hpf,
                                                                        x='kim_2hpf',y='quantSeq_rRNA_2hpf', alpha=0.1,ylim=c(-10,5),title = paste0('n=',nrow(quantSeq_polyA_rnaseq2hpf[is.finite(quantSeq_polyA_rnaseq2hpf$quantSeq_rRNA_2hpf) , ]))) + ggplot2::geom_hline(yintercept = 0, linetype="dashed",color='red')+ 
                          theme_ameres(type = 'barplot')+
                          ylab("3' end seq/rRNA depleted (log2)")  + xlab('Poly(A) tail length') + scale_x_continuous(breaks=number_ticks(10))
                        
                        print(quantSeq_RNAseq_scatterplot)
                        
                        ggplot(quantSeq_polyA_rnaseq2hpf, aes(x=kim_2hpf,y=quantSeq_rRNA_2hpf)) + geom_density2d() + 
                          ggtitle( paste0('n=',nrow(quantSeq_polyA_rnaseq2hpf[is.finite(quantSeq_polyA_rnaseq2hpf$quantSeq_rRNA_2hpf) , ]))) +
                          ggplot2::geom_hline(yintercept = 0, linetype="dashed",color='red')+ theme_ameres(type = 'barplot')+
                          ylab("3' end seq/rRNA depleted (log2)") + xlab('Poly(A) tail length')+ scale_x_continuous(breaks=number_ticks(10))+ ylim(c(-3,3))
                        
                        
                        
                        polyA_RNAseq_scatterplot =  ggpubr::ggscatter(quantSeq_polyA_rnaseq2hpf,
                                                                      x='kim_2hpf',y='polyA_rRNA_2hpf', alpha=0.1,ylim=c(-10,5),title = paste0('n=',nrow(quantSeq_polyA_rnaseq2hpf[is.finite(quantSeq_polyA_rnaseq2hpf$polyA_rRNA_2hpf) , ])))+ ggplot2::geom_hline(yintercept = 0, linetype="dashed",color='red')+ 
                          theme_ameres(type = 'barplot') +
                          ylab("poly(A) RNAseq/rRNA depleted (log2)") + xlab('Poly(A) tail length') + scale_x_continuous(breaks=number_ticks(10))
                        print(polyA_RNAseq_scatterplot)
                        
                        
                        
                        ggplot(quantSeq_polyA_rnaseq2hpf, aes(x=kim_2hpf,y=polyA_rRNA_2hpf)) + geom_density2d() + 
                          ggtitle( paste0('n=',nrow(quantSeq_polyA_rnaseq2hpf[is.finite(quantSeq_polyA_rnaseq2hpf$quantSeq_rRNA_2hpf) , ]))) +
                          ggplot2::geom_hline(yintercept = 0, linetype="dashed",color='red')+ theme_ameres(type = 'barplot')+
                          ylab("poly(A) RNAseq/rRNA depleted (log2)") + xlab('Poly(A) tail length')+ scale_x_continuous(breaks=number_ticks(10)) + ylim(c(-3,3))
                        
                        
                dev.off()
                
      #### some genes dont have expression in poly(A) tail.. let's remove those genes
                
      MasterTable = read.table("/Volumes/Macintosh HD/Users/pooja.bhat/Dropbox (VBC)/Paperoutline/figureOutline/paper_ZebrafishSLAMseq/other/rawfigures/figure1/data/masterTable_quantSeqRNAseqRPMs.txt",
                 sep="\t",stringsAsFactors = F, header = T)
      MasterTable = MasterTable[!is.na(MasterTable$Inj_R1_TP1),]
      