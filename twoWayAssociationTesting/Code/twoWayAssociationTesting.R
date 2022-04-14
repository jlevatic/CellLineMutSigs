library(readr)
library(effsize)

#The following cell lines are removed:
#cell line outliers as reported by Abbas Aghababazadeh et. al.
outliers.AbbasAghababazadehEtAl <- c("909253","949168","1330995","724874","909727","1330942","907295","930299","909711","1240125","687807","909255","910546","687985","905949")
#cell lines with abnormal number of mutations 
outliers.nMuts = c("905971", "1660036", "1659823", "907046", "906852", "905968", "905990", "909257", "1240222", "1480371", "908128", "906814", "753596")
#cell lines with abnormal drug response
outliers.Drugs = c("949162", "724839", "924241", "906844", "949093", "909907", "1240199")
#cell lines with likely misclassified cancer type as reported by Salvadores et. al.
outliers.Misclassified = c("687505","687800","688121","905935","905941","905957","906805","906848","906849","906999","907298","907299","908473","909697","909700","909732","909749","909753","909754","910568","910784","910927","930082","1240128","1240173","1240216","1290797","1290806","1298151","1503363")
outliers = c(outliers.AbbasAghababazadehEtAl, outliers.nMuts, outliers.Drugs, outliers.Misclassified)


#' The function performs two way association procedure, testing for statistical significance between a set of binarized features and response in two independent databases. 
#' For a given set of features (i.e., mutational signatures or CFEs) and a pair of responses, calculates p-values by randomization separately for every feature and every cancer type.
#'
#' @param folder.data Folder containing the gdsc, prism, and pscore data.
#' @param folder.output Folder where the files with results will be written.
#' @param path.to.features Path to features in long format (either mutational signatures or CFEs). The table has to contain the following columns: sample_id, cancer_type, feature, bina
#' @param analysis On of the {"prism", "pscore", "same_target"}. Specifies which type of the analysis will be performed wrt the two databases used: GDSC/PRISM, GDSC/PSCORE, or GDSC/GDSC (same target).
#' @param xaxis Name of the drug from GDSC that will be used in the x-dimension.
#' @param yaxis Name of the drug from GDSC or PRISM or gene from PSCORE that will be used in the y-dimension.
#' @param nRandomizations The number of randomizations for p-value calculation (default 100000).
#' @param eff.size.th The effect size threshold. For associations with effect size less than the specified threshold (along the x-dimension), p-values will not be calculated (default 0.2)
#'
#' @return Table with effect sizes and pvalues. Note that pvalue_pospos corresponds to sensitivity, while pvalue_negneg to resistance associations.
#' @export
#'
#' @examples
twoWayAssociations = function(folder.data, folder.output, path.to.features, analysis, xaxis, yaxis, nRandomizations = 100000, eff.size.th = 0.2) {

  if (!file.exists(folder.output)){
    dir.create(file.path(folder.output), recursive = TRUE)
  }
  
  path = paste0(folder.output, "/pvalues_", analysis, "_", sub("/","",sub(" ", "", xaxis)), "_", sub("/","",sub(" ", "", yaxis)), ".csv")
  
  # 1) LOAD DATA
  # A) Read features
  feat = read_csv(path.to.features)
  feat = feat[!(feat$sample_id %in% outliers), ] #filter outliers
  #use only cancer types with at least 8 samples
  freq_cancer.types = table(feat[!duplicated(feat$sample_id),]$cancer_type)
  freq_cancer.types = freq_cancer.types[freq_cancer.types >= 8]
  feat = feat[feat$cancer_type %in% names(freq_cancer.types),]
  
  # 2) LOAD DATA
  # B) Read gdsc data
  gdsc = read_csv(paste0(folder.data, "gdsc_drug_response.csv"))
  gdsc$DRUG_NAME = gsub(" ", "", gdsc$DRUG_NAME) ##### ADDED
  
  gdsc = gdsc[gdsc$DRUG_NAME == xaxis,]
  gdsc = gdsc[!is.na(gdsc$LN_IC50),]
  
  # sel drug id with more samples
  freq = data.frame(table(gdsc$DRUG_ID))
  gdsc = gdsc[gdsc$DRUG_ID == freq$Var1[freq$Freq == max(freq$Freq)],]
  
  # if drug in GDSC2 and GDSC1 use GDSC1
  if (length(unique(gdsc$DATASET)) == 2){gdsc = gdsc[gdsc$DATASET == "GDSC1",]}
  gdsc = gdsc[,c("COSMIC_ID", "DRUG_NAME", "LN_IC50")]
  colnames(gdsc)[3] = "value"
  xmat = gdsc
  
  # load y axis data
  if(analysis == "prism"){
    
    #read PRISM drug response
    pri = read_csv(paste0(folder.data, "prism_drug_response.csv"))
    pri$name = gsub(" ","", pri$name)
    
    pri = pri[pri$name == yaxis,]
    pri = pri[!is.na(pri$ic50),]
    pri$LN_IC50 = log(pri$ic50)
    pri = pri[,c("COSMICID", "name", "LN_IC50")]
    colnames(pri)[1] = "COSMIC_ID"
    colnames(pri)[3] = "value"
    ymat = pri
    
  }else if(analysis == "pscore"){
    
    # read PSCORE crispr values
    psco = read_csv(paste0(folder.data, "pscore_data.csv"))
    psco = psco[psco$gene == yaxis,]
    psco = psco[!is.na(psco$crispr_score),]
    psco = psco[,c("COSMIC_ID", "gene", "crispr_score")]
    colnames(psco)[3] = "value"
    ymat = psco
    
  }else if(analysis == "same_target"){
    ymat = gdsc
  }
  
  ymat = ymat[!duplicated(ymat$COSMIC_ID), ]
  cancer_types = unique(feat$cancer_type)
  
  #loop through cancer types
  res = lapply(cancer_types, function(ct){
    
    print(ct)
    
    #loop through features
    feats = unique(feat$feature[feat$cancer_type == ct])
    
    res2 = lapply(feats, function(fe){
      
      feat_red = feat[feat$cancer_type == ct & feat$feature == fe,]
      feat_x = feat_red[feat_red$sample_id %in% xmat$COSMIC_ID,]
      feat_y = feat_red[feat_red$sample_id %in% ymat$COSMIC_ID,]
      
      #check if there is at least two samples with '1' values and that all samples don't have the same value
      if(min(table(feat_x$bina)) >= 2 & min(table(feat_y$bina)) >= 2 &
         length(unique(feat_x$bina)) == 2 & length(unique(feat_y$bina)) == 2){
        
        df_res = calc_pvalue(feat_red = feat_red, xmat = xmat, ymat = ymat, iterations = nRandomizations, eff.size.th = eff.size.th)
        df_res$xaxis = xaxis
        df_res$yaxis = yaxis
        df_res$feat = fe
        df_res$cancer_type = ct
        df_res$nPosX = sum(feat_x$bina == 1)
        df_res$nPosY = sum(feat_y$bina == 1)
        df_res$nSamplesX = nrow(feat_x)
        df_res$nSamplesY = nrow(feat_y)
        return(df_res)
      }
      
    }) # 2nd lapply
    
    fin = dplyr::bind_rows(res2)
    return(fin)
  }) # 1st lapply
  
  total = dplyr::bind_rows(res)
  
  
  write.csv(total, path, row.names = F)
}



#calculate p-values for a feature by randomization procedure
calc_pvalue = function(feat_red, xmat, ymat, iterations = 100000, eff.size.th = 0.2){
  
  # Merge them
  redg = merge(feat_red, xmat, by.x = "sample_id", by.y = "COSMIC_ID", all = FALSE)
  redp = merge(feat_red, ymat, by.x = "sample_id", by.y = "COSMIC_ID", all = FALSE)
  
  real_es_g = cohen.d(redg$value, as.factor(redg$bina), pooled = TRUE) #effect size of the x-axis, real feture
  real_es_g = real_es_g$estimate
  
  real_es_p = cohen.d(redp$value, as.factor(redp$bina), pooled = TRUE) #effect size of the x-axis, real feture
  real_es_p = real_es_p$estimate
  
  real_score_min = min(c(real_es_g, real_es_p))
  real_score_max = max(c(real_es_g, real_es_p))
  
  if(abs(real_es_g) > eff.size.th) { #do randomization only if GDSC effect size > threshold
    
    ran_score_min = rep(NA, iterations)
    ran_score_max	= rep(NA, iterations)
    
    for(i in 1:iterations) {
      es_g = cohen.d(redg$value, as.factor(sample(redg$bina)), pooled = TRUE) #effect size of the x-axis, randomized feture
      es_g = es_g$estimate
      
      es_p = cohen.d(redp$value, as.factor(sample(redp$bina)), pooled = TRUE) #effect size of the y-axis, randomized feture
      es_p = es_p$estimate
      
      ran_score_min[i] = min(c(es_g, es_p))
      ran_score_max[i] = max(c(es_g, es_p))
    }
    
    pvalue_pospos = (sum(ran_score_min >= real_score_min) + 1) / (length(ran_score_min) + 1)
    pvalue_negneg = (sum(ran_score_max <= real_score_max) + 1) / (length(ran_score_max) + 1)
    
    df_return = data.frame(cohen_d_xaxis = real_es_g, cohen_d_yaxis = real_es_p, 
                           pvalue_negneg = pvalue_negneg, pvalue_pospos = pvalue_pospos, stringsAsFactors = F)
  } else {
    df_return = data.frame(cohen_d_xaxis = real_es_g, cohen_d_yaxis = real_es_p, 
                           pvalue_negneg = NA, pvalue_pospos = NA, stringsAsFactors = F)
  }
  
  return(df_return)
} # end function



