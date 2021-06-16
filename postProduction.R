
alz$other$genesOfInterest$clust_enrich <-openxlsx::read.xlsx(xlsxFile = "Tables.xlsx", sheet = "t3_clust", startRow = 9)

alz$other$genesOfInterest$clust_enrich$name <- NULL

purrr::walk2(
  .x = alz$other$genesOfInterest$clust_enrich, 
  .y = colnames(alz$other$genesOfInterest$clust_enrich),
  .f = function(geneList, name){
    readr::write_tsv(
      x = as.data.frame( na.exclude(geneList) ), 
      file = paste0("enrichr/inputs_clusters/cluster_", name, ".txt"),
      col_names = F)
    })

alz$other$genesOfInterest$clust_enrich_out <- purrr::map(
  .x = list.files(path = "enrichr/outputs"), 
  .f = function(file){
    file_ <- read.delim(file = paste0("enrichr/outputs/", file), dec = ",")
    
    file_ <- dplyr::select(.data = file_, dplyr::everything(), -`Old.P.value`, -`Old.Adjusted.P.value`, -`Odds.Ratio`, -`Combined.Score`)
    
    return(file_)
    })
names(alz$other$genesOfInterest$clust_enrich_out) <- stringr::str_remove(string = list.files(path = "enrichr/outputs"), pattern = "_output.tsv")
names(alz$other$genesOfInterest$clust_enrich_out) <- stringr::str_remove(string = names(alz$other$genesOfInterest$clust_enrich_out), pattern = "cluster_")

alz$other$genesOfInterest$clust_enrich_out <- rlist::list.clean(.data = alz$other$genesOfInterest$clust_enrich_out, fun = function(x) {
  length(x[[1]]) == 0
  })

alz$other$genesOfInterest$clust_enrich_out_main <- purrr::map(
  .x = alz$other$genesOfInterest$clust_enrich_out, 
  .f = function(dataset){
    
    dataset[dataset$Database %in% c("BioCarta_2016", "GO_Biological_Process_2018", "GO_Cellular_Component_2018", "GO_Molecular_Function_2018", "Human_Phenotype_Ontology", "HumanCyc_2016", "Jensen_COMPARTMENTS", "KEGG_2019_Human", "MGI_Mammalian_Phenotype_Level_4_2019", "NCI-Nature_2016", "Panther_2016", "Reactome_2016", "WikiPathways_2019_Human"),]
  })

alz$other$genesOfInterest$clust_enrich_out_more <- purrr::map(
  .x = alz$other$genesOfInterest$clust_enrich_out, 
  .f = function(dataset){
    
    dataset[dataset$Database %nin% c("BioCarta_2016", "GO_Biological_Process_2018", "GO_Cellular_Component_2018", "GO_Molecular_Function_2018", "Human_Phenotype_Ontology", "HumanCyc_2016", "Jensen_COMPARTMENTS", "KEGG_2019_Human", "MGI_Mammalian_Phenotype_Level_4_2019", "NCI-Nature_2016", "Panther_2016", "Reactome_2016", "WikiPathways_2019_Human"),]
  })

openxlsx::write.xlsx(x = alz$other$genesOfInterest$clust_enrich_out_main, file = "clust_enrich_out_main.xlsx")
openxlsx::write.xlsx(x = alz$other$genesOfInterest$clust_enrich_out_more, file = "clust_enrich_out_more.xlsx")






alz$other$compareTFs$toUseWholeNb <- read.delim(file = "enrichr/outputs_main_datasets/toUseWholeNb_output.tsv", dec = ",")

alz$other$compareTFs$toUseWholeNb <- subset(alz$other$compareTFs$toUseWholeNb, subset = alz$other$compareTFs$toUseWholeNb$Database %in% c("ARCHS4_IDG_Coexp", "ARCHS4_Kinases_Coexp", "ARCHS4_TFs_Coexp", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "Enrichr_Submissions_TF-Gene_Coocurrence", "ESCAPE", "Genome_Browser_PWMs", "KEA_2015", "miRTarBase_2017", "NURSA_Human_Endogenous_Complexome", "PPI_Hub_Proteins", "TargetScan_microRNA", "TargetScan_microRNA_2017", "TF_Perturbations_Followed_by_Expression", "TF-LOF_Expression_from_GEO", "Transcription_Factor_PPIs", "TRANSFAC_and_JASPAR_PWMs", "TRRUST_Transcription_Factors_2019"))

alz$other$compareTFs$toUseWholeNb$Term2 <- ifelse(
  test = alz$other$compareTFs$toUseWholeNb$Database %in% c("ARCHS4_IDG_Coexp", "ARCHS4_Kinases_Coexp", "ARCHS4_TFs_Coexp", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "TF_Perturbations_Followed_by_Expression", "TF-LOF_Expression_from_GEO", "TRANSFAC_and_JASPAR_PWMs", "TRRUST_Transcription_Factors_2019"), 
  yes = stringr::str_remove(string = alz$other$compareTFs$toUseWholeNb$Term, pattern = " .*"),
  no =  ifelse(
    
    test = alz$other$compareTFs$toUseWholeNb$Database %in% c("ESCAPE"), 
    yes = stringr::str_remove(
      string = stringr::str_remove(
        string = alz$other$compareTFs$toUseWholeNb$Term, 
        pattern = "-.*"),
      pattern = "^CHiP |^Protein "), 
    no = ifelse(
      
      test = alz$other$compareTFs$toUseWholeNb$Database %in% c("Genome_Browser_PWMs"), 
      yes = stringr::str_remove(
        string = stringr::str_remove(
          string = alz$other$compareTFs$toUseWholeNb$Term, 
          pattern = ".*\\$"), 
        pattern = " .*"), 
      no = ifelse(
        
        test = alz$other$compareTFs$toUseWholeNb$Database %in% c("miRTarBase_2017", "TargetScan_microRNA_2017"), 
        yes = stringr::str_remove(
          string = stringr::str_remove(
            string = alz$other$compareTFs$toUseWholeNb$Term, 
            pattern = "-3p$|-5p$"), 
          pattern = "^[a-z]{3}-"),
        no = ifelse(
          
          test = alz$other$compareTFs$toUseWholeNb$Database %in% c("NURSA_Human_Endogenous_Complexome"), 
          yes = stringr::str_remove(
            string = stringr::str_remove(
              string = alz$other$compareTFs$toUseWholeNb$Term, 
              pattern = "\\)$"),
            pattern = ".*\\("), 
          no = ifelse(
            
            test = alz$other$compareTFs$toUseWholeNb$Database %in% c("TargetScan_microRNA"), 
            yes = stringr::str_remove(
                string = stringr::str_remove(
                  string = stringr::str_remove(
                    string = alz$other$compareTFs$toUseWholeNb$Term, 
                    pattern = "^[A-Z]{7},"),
                  pattern = ",.*"),
                pattern = "-3p$|-5p$"), 
            no = alz$other$compareTFs$toUseWholeNb$Term))))))
