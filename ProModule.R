setwd("E:\\ProModules\\APclustering")


## read the raw miRNA expression data
rawmiRNAexp <- read.csv("BRCAmiRNAexpression.csv", header=TRUE, sep = ",", stringsAsFactors = F, check.names = F, row.names = 1 )

## read the clinical data 
clinicaldata <- read.csv('BRCA_clinical_data.csv',row.names = 1) 

## input the miRNA-gene from miRTarBase
miRNAtarget <- read.csv("BRCAmiRNAtarget.csv", header=TRUE, sep = ",", stringsAsFactors = F, check.names = F, row.names = 1 )




## main functions
ProModule <- function(rawmiRNAexp, clinicaldata, miRNAtarget) {
  
  #record running time
  ptm <- proc.time()
  
  ## computational identification of prognostic miRNA signature
  V1 <- PromiRNAs(rawmiRNAexp, clinicaldata)
  
  message("PromiRNAs over: ", nrow(V1))
  
  ## systematic detection of prognostic miRNA clusters
  V2 <- APclustering(V1)
  
  message("APclustering over: ", length(V2))
  
  ## filling the prognostic miRNA cluster by target genes
  V3 <- Addtarget(V2, miRNAtarget)
  
  message("AddTargene over: ", length(V3))
  
  print(sprintf("Time elapsed for ProModule: %.3f (mins)",
                (proc.time() - ptm)[3]/60))
  
  printModules(V3)
  
  return(V3) 
}


#univariate cox regression analysis
PromiRNAs <- function(rawmiRNAexp, clinicaldata){
 
  t_index <- which(substr(colnames(rawmiRNAexp),14,15) == "01")
  miRNAexpression <- rawmiRNAexp[,t_index]
  miRNAexpression <- t(miRNAexpression)
  
  rownames(clinicaldata) <- toupper(rownames(clinicaldata))
  
  rownames(miRNAexpression) <- substr(rownames(miRNAexpression), 1, 12)
  miRNA_clinical_expression <- merge(clinicaldata, miRNAexpression, 'row.names')

  
  miRNA_clinical_expression$age_at_initial_pathologic_diagnosis <- ifelse(as.numeric(miRNA_clinical_expression$age_at_initial_pathologic_diagnosis) > 50, 0, 1)
  
  miRNA_clinical_expression$stage_event.clinical_stage <- ifelse((miRNA_clinical_expression$stage_event.clinical_stage =="stage i" | miRNA_clinical_expression$stage_event.clinical_stage == "stage ii"), 0, 1)
  
  
  ## univariate cox regression 
  library(survival)
  pv <- data.frame(miRNA = colnames(miRNA_clinical_expression)[7:dim(miRNA_clinical_expression)[2]])
  
  for(i in 7:dim(miRNA_clinical_expression)[2]){ 
    cox_univ <- coxph(Surv(as.numeric(new_death), as.numeric(death_event)) ~ 
                        as.numeric(miRNA_clinical_expression[,i]), miRNA_clinical_expression)
    pv$pval_uni[i-6] <- summary(cox_univ)$coefficients[5]
    pv$hazardratio_uni[i-6] <- exp(summary(cox_univ)$coefficients[1])
    pv$CI_lower_uni[i-6] <- exp(confint(cox_univ)[1])
    pv$CI_higher_uni[i-6] <- exp(confint(cox_univ)[2])
  }
  length(which(pv$pval_uni < 0.05)) 
  pv <- pv[pv$pval_uni < 0.05,]
  
  ## extract the miRNAs with significant prognostic association
  miRNA_cli_exp_multi <- miRNA_clinical_expression[,c(1:6,which(colnames(miRNA_clinical_expression)
                                                                [7:ncol(miRNA_clinical_expression)] %in% pv$miRNA)+6)]
  
  ## multivariable Cox regression 
  for(i in 7:dim(miRNA_cli_exp_multi)[2]){
    cox_multi <- coxph(Surv(as.numeric(new_death), as.numeric(death_event)) ~ 
                         miRNA_cli_exp_multi[,i] + as.numeric(age_at_initial_pathologic_diagnosis) + stage_event.clinical_stage, miRNA_cli_exp_multi)
    pv$pval_multi[i-6] <- summary(cox_multi)$coefficients[13] 
    pv$hazardratio_multi[i-6] <- exp(summary(cox_multi)$coefficient[1])
    pv$CI_lower_multi[i-6] <- exp(confint(cox_multi)[1])
    pv$CI_higher_multi[i-6] <- exp(confint(cox_multi)[4])
  }
  length(which(pv$pval_multi < 0.05)) 
  pv <- pv[pv$pval_multi < 0.05,]
  
  ## further narrow the miRNA list that are independent of other clinical factors such as age of initial diagnosis, clinical pathologic stage
  promiRNAexp <- rawmiRNAexp[as.character(pv$miRNA),]
  
  return(promiRNAexp)
}


## initialize module 
module <- function(miRNA, gene){
  list(miRNA  = miRNA,
       gene   = gene)
}

## clustering stage for assigning miRNA clusters
APclustering <- function(miRNAexpression){
  
  ## maximal information coefficient 
  library(minerva)
  
  miRNAMINE <- mine(t(miRNAexpression))
  
  miRNAcoefficient <- miRNAMINE$MIC
  

  ## run affinity propagation algorithm
  library(apcluster)
  
  ## run affinity propagation 
  apresults <- apcluster(negDistMat(r=2), miRNAcoefficient, q=0.97)
  
  show(apresults) 
  
  heatmap(apresults)
  
  
  ## retain the modules with the number of miRNA more than 2
  clusters <- list()
  for(i in 1:length(apresults)){
    if(length(apresults[[i]])>=2){
      temp  <- module(apresults[[i]], NULL)
      clusters <- c(clusters, list(temp))
    }
  }
  return(clusters)
}


## adding shared target genes for each module
Addtarget <- function(clusters, miRNAtarget){
  
  for(i in 1:length(clusters)){ 
    commongenes <- c()
    for(j in 1:length(clusters[[i]]$miRNA)){
      commongenes <- c(commongenes, colnames(miRNAtarget)[which(miRNAtarget[as.numeric(clusters[[i]]$miRNA)[j],]>0)])
    }
    commongenes <- as.data.frame(table(commongenes))
    commongenes <- as.character(commongenes[which(commongenes$Freq>=2),1]) ## shared target gene frequence more than 2
    clusters[[i]]$gene <- commongenes
  }
  return(clusters)
}


## print prognostic miRNA modules ##
printModules <- function(V) {
  names(V) <- paste("M",1:length(V),sep="")
  print_helper <- function(m) {
    v <- V[[m]]
    miRNA <- names(v$miRNA)
    gene <- v$gene
    cat(sprintf("%s: \n", m))
    cat(miRNA, "\n")
    cat(gene, "\n")	
  }	
  tmp <- lapply(names(V), print_helper)
}

V <- ProModule(rawmiRNAexp, clinicaldata, miRNAtarget)


