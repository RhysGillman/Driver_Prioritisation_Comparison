#annotate_LOF_GOFs.R
# This code annotates CCLE mutations as LOF or GOF

# Load Required Packages
suppressPackageStartupMessages (library(optparse, quietly = T))
suppressPackageStartupMessages (library(tidyverse, quietly = T))
suppressPackageStartupMessages (library(data.table, quietly = T))
suppressPackageStartupMessages (library(jsonlite, quietly = T))
suppressPackageStartupMessages (library(httr, quietly = T))
suppressPackageStartupMessages (library(multidplyr, quietly = T))


# Handling input arguments
option_list = list(
  make_option(c("-t", "--threads"), type="integer", default=4, 
              help="Number of threads to use if running in parallel", metavar ="Threads")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

threads <- opt$threads


#cl <- new_cluster(threads)
#cl

#############################
# Get cell line information
###############################

sample_info <- read_csv("data/CCLE/Model.csv") %>%
  dplyr::rename(DepMap_ID=ModelID, lineage=OncotreeLineage, cell_ID=StrippedCellLineName)






#Paper: https://www.sciencedirect.com/science/article/pii/S0002929721003840#app3

#There is very limited availability of known annotated GOF vs LOF mutations. 
#The above paper has attempted to rectify this using a large NLP / textmining 
#initiative to automatically annotate these mutations using machine learning. 
#The full prediction of these annotations can be accessed here https://gitlab.com/itan-lab/logofunc-predictions

#####################################
# Fixing Itan Lab Variant Annotation
#####################################

# Variant annotation is not consistent between the LOF/GOF annotation data and CCLE mutation data
# Manually correcting these variant annotations proved very difficult
# Instead using Ensembl's Variant Recoder to acquire variant annotation synonyms based on variant IDs


## Aquiring synonymous variant annotations
# This step is very slow due to instability of the ENSEMBL server connection
# The code will continuously attempt to reconnect each time it fails
# Can take > 10 hours to complete

if(!dir.exists("cache/clinvar_annotation_synonyms")){
  suppressWarnings(dir.create("cache/clinvar_annotation_synonyms", recursive = T))
}

if(length(list.files("cache/clinvar_annotation_synonyms/")) == 0){
  clinvar <- fread("data/LOFGOF/goflof_ClinVar_v062021.csv")
  colnames(clinvar) <- c("allele_ID", "label", "chrom", "pos", "ref", "alt", "snp_ID", "gene_ID", "name")
  clinvar <- clinvar %>%
    # Getting a copy of RefSeq gene ID with and without version number
    mutate(RefSeq_v = gsub("\\(.*","",name)) %>%
    mutate(RefSeq = gsub("\\..*","",RefSeq_v)) %>%
    mutate(name = gsub(".*:","",name)) %>%
    # Retrieving cDNA and protein changes from the "name" column, and formatting
    tidyr::separate(name, into = c("cDNA_Change", "Protein_Change"), sep = " ", fill = "right") %>%
    mutate(Protein_Change = gsub("\\(|\\)","",Protein_Change)) %>%
    mutate(snp_ID = ifelse(!is.na(snp_ID), paste0("rs",snp_ID), snp_ID))
  
  
  # Because the ENSEMBL REST API imposes limits on the size of requests, need to break up data retrieval into batches so that there are < 200 requests per run
  batches = 50
  # Getting rs IDs
  id_list <- unique(na.omit(clinvar$snp_ID))
  # Splitting into 25 batches
  id_list <- split(id_list, rep(1:batches, length.out = length(id_list)))
  # Variables for contacting ENSEMBL REST API
  server <- "https://rest.ensembl.org"
  ext <- "/variant_recoder/homo_sapiens"
  
  # Setting up loop to retrieve variant synonyms and format them correctly
  
  
  
  for(batch in seq(1,batches)){
    
    # Checking connection to server. Gives up after 5 minutes of no connection
    ping <- GET("https://rest.ensembl.org/info/ping?", content_type("application/json"))
    ping <- head(fromJSON(toJSON(content(ping))))
    i=1
    while(ping != 1){
      print("No connection to ensembl REST API. Retrying in 5 seconds...")
      ping <- GET("https://rest.ensembl.org/info/ping?", content_type("application/json"))
      ping <- head(fromJSON(toJSON(content(ping))))
      Sys.sleep(10)
      i=i+1
      if(i > 30){stop("Couldn't connect to ensembl REST API")}
    }
    if(ping ==1){
      print("Successful connection to ensembl REST API.")
    }
    
    
    print(paste0("Beginning Batch: ", batch, "."))
    # Collecting batch of IDs and formatting as required for the API
    id_list_batch <- id_list[batch] %>% unlist()
    id_list_batch <- paste0('"',id_list_batch, '"')
    id_list_batch <- paste(id_list_batch, collapse = ", ")
    # Performing request and checking it worked
    annotation_synonyms <- RETRY("POST", paste(server, ext, sep = ""), content_type("application/json"), accept("application/json"), body = noquote(paste0('{ "ids" : [', id_list_batch, '] }')))
    stop_for_status(annotation_synonyms)
    
    #write_rds(annotation_synonyms, paste0("cache/clinvar_annotation_synonyms/result_batch_",batch,".rds"))
    
    rm(annotation_synonyms)
    print(paste0("Finished Batch: ", batch, ". Waiting 5 minutes and beggining next batch."))
    Sys.sleep(60)
  }
}
# HGMD

if(!dir.exists("cache/HGMD_annotation_synonyms")){
  suppressWarnings(dir.create("cache/HGMD_annotation_synonyms", recursive = T))
}


if(length(list.files("cache/HGMD_annotation_synonyms/")) == 0){

  HGMD <- fread("data/LOFGOF/goflof_HGMD2019_v032021_allfeat.csv") %>% 
    dplyr::select(gene_ID = GENE, chrom = CHROM, pos = POS, cDNA_Change = HGVSc, Protein_Change = HGVSp,snp_ID = ID, label = LABEL) %>%
    mutate(full_cDNA_Change_id = cDNA_Change) %>%
    tidyr::separate(cDNA_Change, into = c("ensembl_transcript_id_v", "cDNA_Change"), sep = ":", fill = "right") %>%
    # Removing version numbers from the ensembl transcript IDs
    mutate(ensembl_transcript_id = gsub("\\.[1-9]*$","",ensembl_transcript_id_v)) %>%
    # Removing a single row with a missing ensembl transcript ID
    filter(ensembl_transcript_id != "-") %>%
    tidyr::separate(Protein_Change, into = c("ensembl_protein_id", "Protein_Change"), sep = ":", fill = "right")
  
  # Because the ENSEMBL REST API imposes limits on the size of requests, need to break up data retrieval into batches so that there are < 200 requests per run
  batches = 100
  # Getting rs IDs
  id_list <- unique(na.omit(HGMD$full_cDNA_Change_id))
  # Splitting into 25 batches
  id_list <- split(id_list, rep(1:batches, length.out = length(id_list)))
  # Variables for contacting ENSEMBL REST API
  server <- "https://rest.ensembl.org"
  ext <- "/variant_recoder/homo_sapiens"
  
  # Setting up loop to retrieve variant synonyms and format them correctly
  
  
  
  for(batch in seq(1,batches)){
    
    # Checking connection to server. Gives up after 5 minutes of no connection
    ping <- GET("https://rest.ensembl.org/info/ping?", content_type("application/json"))
    ping <- head(fromJSON(toJSON(content(ping))))
    i=1
    while(ping != 1){
      print("No connection to ensembl REST API. Retrying in 5 seconds...")
      ping <- GET("https://rest.ensembl.org/info/ping?", content_type("application/json"))
      ping <- head(fromJSON(toJSON(content(ping))))
      Sys.sleep(10)
      i=i+1
      if(i > 30){stop("Couldn't connect to ensembl REST API")}
    }
    if(ping ==1){
      print("Successful connection to ensembl REST API.")
    }
    
    
    print(paste0("Beginning Batch: ", batch, "."))
    # Collecting batch of IDs and formatting as required for the API
    id_list_batch <- id_list[batch] %>% unlist()
    id_list_batch <- paste0('"',id_list_batch, '"')
    id_list_batch <- paste(id_list_batch, collapse = ", ")
    # Performing request and checking it worked
    annotation_synonyms <- RETRY("POST", paste(server, ext, sep = ""), content_type("application/json"), accept("application/json"), body = noquote(paste0('{ "ids" : [', id_list_batch, '] }')))
    stop_for_status(annotation_synonyms)
    
    write_rds(annotation_synonyms, paste0("cache/HGMD_annotation_synonyms/result_batch_",batch,".rds"))
    
    rm(annotation_synonyms)
    print(paste0("Finished Batch: ", batch, ". Waiting 5 minutes and beggining next batch."))
    Sys.sleep(60)
  }
}




# Reformatting ENSEMBL Recoder Results

if(!file.exists("cache/clinvar_annotation_synonyms.csv")){
  suppressWarnings(rm("clinvar_synonyms"))
  for(file in list.files("cache/clinvar_annotation_synonyms/")){
    message(paste0("reading batch ",which(list.files("cache/clinvar_annotation_synonyms/")==file),"/", length(list.files("cache/clinvar_annotation_synonyms/"))) )
    #file <- "result_batch_1.rds"
    # Read in the file
    tmp <- read_rds(paste0("cache/clinvar_annotation_synonyms/",file))
    # Extract info from JSON
    tmp <- content(tmp,"text", encoding = "UTF-8")
    tmp <- fromJSON(tmp) 
    # get elements excluding warnings
    elements <- names(tmp)[which(names(tmp)!="warnings")]
    
    suppressWarnings(rm("tmp3"))
    for(element in elements){
      #element <- "C"
      tmp2 <- tmp[element] %>% 
        unnest(everything()) %>%
        na.omit() %>%
        mutate(hgvsg_dummy=NA,
               hgvsc_dummy=NA,
               hgvsp_dummy=NA,
               spdi_dummy=NA) %>%
        unnest_wider(col = -input, names_sep = "_") %>%
        as.data.table()
      
      tmp2 <- tmp2 %>%
        data.table::melt(measure.vars=colnames(tmp2)[(which(str_detect(colnames(tmp2),"hgvsg")))],variable.name="hgvsg_ref",value.name="hgvsg") %>%
        data.table::melt(measure.vars=colnames(tmp2)[(which(str_detect(colnames(tmp2),"hgvsc")))],variable.name="hgvsc_ref",value.name="hgvsc") %>%
        data.table::melt(measure.vars=colnames(tmp2)[(which(str_detect(colnames(tmp2),"hgvsp")))],variable.name="hgvsp_ref",value.name="hgvsp") %>%
        data.table::melt(measure.vars=colnames(tmp2)[(which(str_detect(colnames(tmp2),"spdi")))],variable.name="spdi_ref",value.name="spdi") %>%
        dplyr::select(input, hgvsg, hgvsc, hgvsp, spdi)
      tmp2 <- tmp2 %>%
        data.table::melt(measure.vars=colnames(tmp2)[(which(colnames(tmp2)!="input"))],variable.name="annotation_type",value.name="variant") %>%
        unique() %>% filter(!is.na(variant))
      
      if(!exists("tmp3")){
        tmp3 <- tmp2
      }else{
        tmp3 <- rbind(tmp3,tmp2) %>% unique()
      }
    
      if(!exists("clinvar_synonyms")){
        clinvar_synonyms <- tmp3 %>% unique()
      }else{
        clinvar_synonyms <- rbind(clinvar_synonyms,tmp3) %>% unique()
      }
      message(paste0(round(100*which(elements==element)/length(elements),0),"%"))
    }
  }
  write_csv(clinvar_synonyms,"cache/clinvar_annotation_synonyms.csv")
}else{
  clinvar_synonyms <- read_csv("cache/clinvar_annotation_synonyms.csv")
}

if(!file.exists("cache/HGMD_annotation_synonyms.csv")){
  suppressWarnings(rm("HGMD_synonyms"))
  for(file in list.files("cache/HGMD_annotation_synonyms/")){
    message(paste0("reading batch ",which(list.files("cache/HGMD_annotation_synonyms/")==file),"/", length(list.files("cache/HGMD_annotation_synonyms/"))) )
    # Read in the file
    tmp <- read_rds(paste0("cache/HGMD_annotation_synonyms/",file))
    # Extract info from JSON
    tmp <- content(tmp,"text", encoding = "UTF-8")
    tmp <- fromJSON(tmp) 
    # get elements excluding warnings
    elements <- names(tmp)[which(names(tmp)!="warnings")]
    
    suppressWarnings(rm("tmp3"))
    for(element in elements){
      tmp2 <- tmp[element] %>% 
        unnest(everything()) %>%
        na.omit() %>%
        mutate(hgvsg_dummy=NA,
               hgvsc_dummy=NA,
               hgvsp_dummy=NA,
               spdi_dummy=NA,
               id_dummy=NA) %>%
        unnest_wider(col = -input, names_sep = "_") %>%
        as.data.table()
      
      tmp2 <- tmp2 %>%
        data.table::melt(measure.vars=colnames(tmp2)[(which(str_detect(colnames(tmp2),"id")))],variable.name="id_ref",value.name="id") %>%
        data.table::melt(measure.vars=colnames(tmp2)[(which(str_detect(colnames(tmp2),"hgvsg")))],variable.name="hgvsg_ref",value.name="hgvsg") %>%
        data.table::melt(measure.vars=colnames(tmp2)[(which(str_detect(colnames(tmp2),"hgvsc")))],variable.name="hgvsc_ref",value.name="hgvsc") %>%
        data.table::melt(measure.vars=colnames(tmp2)[(which(str_detect(colnames(tmp2),"hgvsp")))],variable.name="hgvsp_ref",value.name="hgvsp") %>%
        data.table::melt(measure.vars=colnames(tmp2)[(which(str_detect(colnames(tmp2),"spdi")))],variable.name="spdi_ref",value.name="spdi") %>%
        dplyr::select(input,id, hgvsg, hgvsc, hgvsp, spdi) %>%
        dplyr::filter(str_detect(id,"rs")|is.na(id))
      tmp2 <- tmp2 %>%
        data.table::melt(measure.vars=colnames(tmp2)[(which(colnames(tmp2)!="input"))],variable.name="annotation_type",value.name="variant") %>%
        unique() %>% filter(!is.na(variant))
      
      if(!exists("tmp3")){
        tmp3 <- tmp2
      }else{
        tmp3 <- rbind(tmp3,tmp2) %>% unique()
      }
      
      if(!exists("HGMD_synonyms")){
        HGMD_synonyms <- tmp3 %>% unique()
      }else{
        HGMD_synonyms <- rbind(HGMD_synonyms,tmp3) %>% unique()
      }
      message(paste0(round(100*which(elements==element)/length(elements),0),"%"))
    }
  }
  write_csv(HGMD_synonyms,"cache/HGMD_annotation_synonyms.csv")
}else{
  HGMD_synonyms <- read_csv("cache/HGMD_annotation_synonyms.csv")
}

# Adding synonyms

clinvar <- fread("data/LOFGOF/goflof_ClinVar_v062021.csv")
colnames(clinvar) <- c("allele_ID", "label", "chrom", "pos", "ref", "alt", "snp_ID", "gene_ID", "name")
clinvar <- clinvar %>%
  tidyr::separate(name, into = c("cDNA_Change", "Protein_Change"), sep = " ", fill = "right") %>%
  mutate(Protein_Change = gsub("\\(|\\)","",Protein_Change)) %>%
  mutate(snp_ID = ifelse(!is.na(snp_ID), paste0("rs",snp_ID), snp_ID))
  
clinvar_annotations <- clinvar_synonyms %>%
  filter(!annotation_type%in%c("spdi","hgvsg")) %>%
  mutate(annotation_type=ifelse(annotation_type=="hgvsc","cDNA_Change",ifelse(annotation_type=="hgvsp","Protein_Change",annotation_type)
  )) %>%
  left_join(clinvar %>% dplyr::select(snp_ID,gene_ID,label), by = "snp_ID", multiple="all") %>%
  mutate(variant=gsub(".*:","",variant)) %>%
  unique()

clinvar_annotations <- rbind(
  clinvar_annotations %>% dplyr::select(gene_ID,annotation_type,variant,label),
  clinvar_annotations %>% dplyr::select(gene_ID,snp_ID,label) %>% pivot_longer(snp_ID,names_to = "annotation_type",values_to = "variant")
) %>% 
  rbind(
    clinvar %>% dplyr::select(gene_ID,cDNA_Change,Protein_Change,snp_ID,label) %>% mutate(cDNA_Change=gsub(".*:","",cDNA_Change)) %>% pivot_longer(cols = c(cDNA_Change,Protein_Change,snp_ID),names_to = "annotation_type",values_to = "variant")
  ) %>% unique()



HGMD <- fread("data/LOFGOF/goflof_HGMD2019_v032021_allfeat.csv") %>% 
  dplyr::select(gene_ID = GENE, chrom = CHROM, pos = POS, cDNA_Change = HGVSc, Protein_Change = HGVSp,snp_ID = ID, label = LABEL) %>%
  mutate(full_cDNA_Change_id = cDNA_Change) %>%
  mutate(cDNA_Change=gsub(".*:","",cDNA_Change)) %>%
  mutate(Protein_Change=gsub(".*:","",Protein_Change))


HGMD_annotations <- HGMD_synonyms %>%
  filter(!annotation_type%in%c("spdi","hgvsg")) %>%
  mutate(annotation_type=ifelse(annotation_type=="hgvsc","cDNA_Change",
                                ifelse(annotation_type=="hgvsp","Protein_Change",ifelse(annotation_type=="id","snp_ID",annotation_type))
  )) %>%
  left_join(HGMD %>% dplyr::select(full_cDNA_Change_id,gene_ID,label), by = c("input"="full_cDNA_Change_id"), multiple="all") %>%
  mutate(variant=gsub(".*:","",variant)) %>%
  unique() %>% 
  relocate(gene_ID) %>%
  dplyr::select(-input) %>%
  rbind(
    HGMD %>% dplyr::select(gene_ID,cDNA_Change,Protein_Change,label) %>% pivot_longer(cols = c(cDNA_Change,Protein_Change),names_to = "annotation_type",values_to = "variant")
  ) %>% unique()

merged_annotations <- rbind(
  clinvar_annotations,
  HGMD_annotations
) %>% unique() %>%
  na.omit() %>%
  group_by(gene_ID,annotation_type,variant) %>%
  summarise(label=paste0(label, collapse = "")) %>%
  ungroup() %>%
  mutate(itan1_LOF=label %in% c("LOF","LOFGOF","GOFLOF"),itan1_GOF=label %in% c("GOF","LOFGOF","GOFLOF")) %>%
  dplyr::select(-label)


#####################################
# Cancer Reference Genes
#####################################

CGC <- read_csv("data/cancer_reference_genes/CGC/cancer_gene_census.csv") %>%
  # Fixing ambiguous roles
  mutate(CGC_oncogene=str_detect(tolower(`Role in Cancer`),"onco")) %>%
  mutate(CGC_tsg=str_detect(tolower(`Role in Cancer`),"tsg|suppress")) %>%
  dplyr::select(gene_ID = `Gene Symbol`,CGC_oncogene,CGC_tsg) %>%
  na.omit() %>%
  group_by(gene_ID) %>%
  summarise(CGC_oncogene=as.logical(max(CGC_oncogene)),CGC_tsg=as.logical(max(CGC_tsg))) %>%
  ungroup()


NCG <- read_tsv("data/cancer_reference_genes/NCG/NCG_cancerdrivers_annotation_supporting_evidence.tsv") %>%
  mutate(NCG_oncogene=as.logical(NCG_oncogene), NCG_tsg=as.logical(NCG_tsg)) %>%
  dplyr::select(gene_ID = symbol, NCG_oncogene, NCG_tsg) %>%
  na.omit() %>%
  group_by(gene_ID) %>%
  summarise(NCG_oncogene=as.logical(max(NCG_oncogene)),NCG_tsg=as.logical(max(NCG_tsg))) %>%
  ungroup()

CancerMine <- read_tsv("data/cancer_reference_genes/CancerMine/cancermine_collated.tsv") %>%
  # Filtering for >1 citation count as done in PersonaDrive
  filter(citation_count >2) %>%
  filter(role%in%c("Oncogene","Tumor_Suppressor")) %>%
  dplyr::select(gene_ID=gene_normalized, role, citation_count) %>%
  unique() %>%
  group_by(gene_ID,role) %>%
  summarise(citation_count=sum(citation_count)) %>%
  ungroup() %>% group_by(gene_ID) %>%
  arrange(desc(citation_count)) %>%
  filter(row_number()==1) %>%
  mutate(CancerMine_oncogene=role=="Oncogene") %>%
  mutate(CancerMine_tsg=role=="Tumor_Suppressor") %>%
  dplyr::select(gene_ID,CancerMine_oncogene,CancerMine_tsg)
  

###############################
# Annotating CCLE Mutations
################################

# Prepare function to convert 1-letter to 3-letter code for amino acids
AAs <- c("*", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
names(AAs) <- c("Ter", "Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr")

AtoAAA <- function(A,AAs){
  return(names(AAs[which(AAs==A)]))
}

# Read in CCLE mutation data
mutation <- fread("data/LOFGOF/annotated_CCLE_mutations.csv",
                  select = c("DepMap_ID","HugoSymbol","Chrom","Pos","Ref","Alt","GT","Transcript", "DNAChange",
                             "ProteinChange","DbsnpID","CivicID","VariantType","VariantInfo",
                             "AssociatedWith","LikelyLoF","LikelyGoF","id","prediction",
                             "score_neutral","score_GOF","score_LOF"))

# Subsample for testing only
#mutation <- mutation[sample(seq(from=1,to=nrow(mutation)),10000),]
  
mutation <- mutation %>%  
  # Separate multiple matching rsIDs to new rows
  separate_longer_delim(cols = DbsnpID, delim = " ") %>%
  dplyr::rename(hg38pos=Pos,  chrom=Chrom, pos=Pos, ref=Ref, alt=Alt, snp_ID=DbsnpID, gene_ID=HugoSymbol, cDNA_Change=DNAChange, Protein_Change=ProteinChange) %>%
  
  # Convert 1 letter amino acid code to 3 letter
  mutate(tmp_p1=str_match(Protein_Change,"p\\.([A-Z]+)[0-9]+[A-Z\\*]*")[,2],
         tmp_p2=str_match(Protein_Change,"p\\.[A-Z]+([0-9]+)[A-Z\\*]*")[,2],
         tmp_p3=str_match(Protein_Change,"p\\.[A-Z]+[0-9]+([A-Z\\*]*)")[,2],
         tmp_del=str_detect(Protein_Change,"del")) %>%
  # Can't fix insertion annotaions because the trailing amino acid after the insertion isn't given, required for correct annotation
  mutate(tmp_ins=nchar(tmp_p3)>nchar(tmp_p1),tmp_point=nchar(tmp_p3)==1&nchar(tmp_p1)==1) %>%

  # Fixing incorrect protein deletion annotations
  mutate(tmp_del_length=ifelse(tmp_del,nchar(tmp_p1),NA)) %>%
  mutate(tmp_del_pos1=ifelse(tmp_del,as.numeric(tmp_p2)-tmp_del_length+1,NA)) %>%
  rowwise() %>%
  mutate(Protein_Change=ifelse(tmp_point, paste0("p.",AtoAAA(tmp_p1,AAs),tmp_p2,AtoAAA(tmp_p3,AAs)),Protein_Change)) %>%
  mutate(Protein_Change=ifelse(tmp_del, paste0("p.",AtoAAA(substr(tmp_p1,1,1),AAs),tmp_del_pos1,"_",AtoAAA(substr(tmp_p1,nchar(tmp_p1),nchar(tmp_p1)),AAs),tmp_p2,"del"),Protein_Change)) %>%
  dplyr::select(-starts_with("tmp")) %>%
  # Getting LOF / GOF info from depmap annotations
  mutate(associated_with_loss = str_detect(AssociatedWith, "loss"),
         associated_with_gain = str_detect(AssociatedWith, "gain"),
         DepMap_LikelyLoF = LikelyLoF=="Y",
         DepMap_LikelyGoF = LikelyGoF=="Y",
         itan2_LOF=score_LOF>score_GOF,
         itan2_GOF=score_GOF>score_LOF
         )
# Annotating CCLE mutations with itan lab LOF/GOF annotations
mutation <- mutation %>%
  left_join(merged_annotations %>% filter(annotation_type=="cDNA_Change") %>% dplyr::select(gene_ID,variant,cDNA_itan1_LOF=itan1_LOF,cDNA_itan1_GOF=itan1_GOF), by = c("gene_ID","cDNA_Change"="variant"), multiple="all") %>%
  left_join(merged_annotations %>% filter(annotation_type=="Protein_Change") %>% dplyr::select(gene_ID,variant,Protein_itan1_LOF=itan1_LOF,Protein_itan1_GOF=itan1_GOF), by = c("gene_ID","Protein_Change"="variant"), multiple="all") %>%
  left_join(merged_annotations %>% filter(annotation_type=="snp_ID") %>% dplyr::select(variant,id_itan1_LOF=itan1_LOF,id_itan1_GOF=itan1_GOF), by = c("snp_ID"="variant"), multiple="all") %>%
  mutate(itan1_LOF= ifelse(any(!is.na(c(cDNA_itan1_LOF,Protein_itan1_LOF,id_itan1_LOF))),max(cDNA_itan1_LOF,Protein_itan1_LOF,id_itan1_LOF,na.rm = T),NA),
         itan1_GOF= ifelse(any(!is.na(c(cDNA_itan1_GOF,Protein_itan1_GOF,id_itan1_GOF))),max(cDNA_itan1_GOF,Protein_itan1_GOF,id_itan1_GOF,na.rm = T),NA)) %>%
  dplyr::select(-c(cDNA_itan1_LOF,Protein_itan1_LOF,id_itan1_LOF,cDNA_itan1_GOF,Protein_itan1_GOF,id_itan1_GOF))
                                 
# Annotating CCLE mutations with CGC, NCG, and CancerMine oncogene/TSG annotations

mutation <- mutation %>%
  left_join(CGC, by = "gene_ID", multiple = "all") %>%
  left_join(NCG, by = "gene_ID", multiple = "all") %>%
  left_join(CancerMine, by = "gene_ID", multiple = "all")



# Adding up different sources of information

CGC_weight=5
NCG_weight=5
CancerMine_weight=1

mutation <- mutation %>%
  # Gain of functin or Oncogene
  mutate(GOF_total_score=sum(associated_with_gain,DepMap_LikelyGoF,itan2_GOF,itan1_GOF,CancerMine_oncogene*CancerMine_weight,NCG_oncogene*NCG_weight,CGC_oncogene*CGC_weight, na.rm = T)) %>%
  # Loss of function or TSG
  mutate(LOF_total_score=sum(associated_with_loss,DepMap_LikelyLoF,itan2_LOF,itan1_LOF,CancerMine_tsg*CancerMine_weight,NCG_tsg*NCG_weight,CGC_tsg*CGC_weight, na.rm = T)) %>%
  mutate(final_prediction=ifelse(GOF_total_score>LOF_total_score, "GOF", "LOF"))


LOF_GOF_annotation_stats <- mutation %>%
  mutate(annotation_available=sum(associated_with_gain,DepMap_LikelyGoF,itan2_GOF,itan1_GOF,CancerMine_oncogene,NCG_oncogene,CGC_oncogene,
                                  associated_with_loss,DepMap_LikelyLoF,itan2_LOF,itan1_LOF,CancerMine_tsg,NCG_tsg,CGC_tsg, na.rm = T)) %>%
  mutate(annotation_available=annotation_available!=0)


write("Loss-of-function / Gain-of-function Mutation Predictions\n--------------\n", "log/LOF_GOF_annotation_stats.txt")
write.table(LOF_GOF_annotation_stats %>% group_by(annotation_available) %>% summarise(count=n()), "log/LOF_GOF_annotation_stats.txt", append = T, sep = "\t", row.names = F)
write("\n--------------\n", "log/LOF_GOF_annotation_stats.txt", append = T)

# For missing annotations, looking for other annotations in the same gene

missing_annotations <- LOF_GOF_annotation_stats %>%
  filter(!annotation_available)

replacement_annotations <- LOF_GOF_annotation_stats %>%
  filter(gene_ID %in% unique(missing_annotations$gene_ID)) %>%
  filter(annotation_available) %>%
  # Get unique list of annotations for each gene
  ungroup() %>%
  group_by(gene_ID,pos,ref) %>%
  summarise(LOF_total_score=max(LOF_total_score), GOF_total_score=max(GOF_total_score)) %>%
  ungroup() %>%
  unique() %>%
  # Get the sum of scores for GOF and LOF annotations
  group_by(gene_ID) %>%
  summarise(LOF_total_score=sum(LOF_total_score),GOF_total_score=sum(GOF_total_score)) %>%
  ungroup() %>%
  mutate(final_prediction=ifelse(GOF_total_score>LOF_total_score, "GOF", "LOF"),
         annotation_available="predicted at gene-level")

# Replace annotations
fixed_annotations <- missing_annotations %>%
  dplyr::select(-final_prediction,-annotation_available) %>%
  left_join(replacement_annotations %>% dplyr::select(gene_ID,final_prediction,annotation_available), by = "gene_ID") %>%
  ungroup() %>%
  # Remaining missing annotations are assigned LOF
  mutate(final_prediction = ifelse(is.na(final_prediction), "LOF",final_prediction),
         annotation_available = ifelse(is.na(annotation_available),FALSE,annotation_available)) %>%
  # Rejoin with available annotations
  rbind(LOF_GOF_annotation_stats %>% filter(annotation_available))


fixed_annotations %>% group_by(annotation_available) %>% summarise(n())

#annotation_available     `n()`
#<chr>                    <int>
#1 FALSE                    16638
#2 TRUE                    956471
#3 predicted at gene-level 419295

fixed_annotations %>% group_by(final_prediction) %>% summarise(n())

#final_prediction   `n()`
#<chr>              <int>
#1 GOF               278410
#2 LOF              1113994




# Finally, replace mutation file with these annotations


mutation <- fixed_annotations



#######################
# Copy Number Variants
#######################

cnv <- fread("data/CCLE/corrected_cnv.csv") %>%
  pivot_longer(-gene_ID, names_to = "cell_ID", values_to = "cnv") %>%
  filter(cnv != 0) %>%
  mutate(cnv_LOF_GOF=ifelse(cnv==1, "GOF","LOF"))




#######################
# Combining
#######################

all_annotations <- mutation %>%
  left_join(sample_info %>% dplyr::select(DepMap_ID,cell_ID), by = "DepMap_ID") %>%
  full_join(cnv, by=c("gene_ID","cell_ID")) %>%
                                      # Check if cnv_LOF_GOF isn't NA
  mutate(combined_prediction = ifelse(!is.na(cnv_LOF_GOF), 
                                      # TRUE: Check if SNP prediction is NA
                                      ifelse(is.na(final_prediction), 
                                             #TRUE: Make it the CNV annotation
                                             cnv_LOF_GOF, 
                                             # FALSE: Check if the CNV annotation matches the SNP annotation
                                             ifelse(final_prediction==cnv_LOF_GOF,
                                                    # TRUE: Make it the SNP annotation
                                                    final_prediction,
                                                    # FALSE: Make it "both"
                                                    "both")),
                                      # FALSE: Make it the SNP annotation
                                      final_prediction
                                      )
         
         )

all_annotations %>% group_by(combined_prediction) %>% summarise(n())

#combined_prediction   `n()`
#<chr>                 <int>
# 1 GOF                  610742
#2 LOF                 1457692
#3 both                   9501

write_csv(all_annotations,"data/LOF_GOF_annotations_all_evidence.csv")

#all_annotations <- fread("data/LOF_GOF_annotations_all_evidence.csv")


#################
# Final Output
################

LOF_GOF_annotations <- all_annotations %>%
  dplyr::select(cell_ID,gene_ID,annotation=combined_prediction) %>%
  unique() %>%
  ungroup() %>%
  # Need to indicate if a single gene has both GOF and LOF mutations
  group_by(cell_ID,gene_ID) %>%
  summarise(annotation=paste0(annotation, collapse = "and")) %>%
  ungroup() %>%
  mutate(annotation=ifelse(str_detect(annotation,"and"),"both",annotation))

#annotation   count
#<chr>        <int>
#1 GOF         555236
#2 LOF        1283672
#3 both         30989

summary <- LOF_GOF_annotations %>% group_by(annotation) %>% summarise(count=n())

write_csv(summary,"validation_data/LOF_GOF_total_summary.csv")

summary <- LOF_GOF_annotations %>% group_by(cell_ID,annotation) %>% summarise(count=n()) %>% pivot_wider(names_from = annotation, values_from = count)

write_csv(summary,"validation_data/LOF_GOF_cell_summary.csv")


write("Final predictions made. Note: Missing predictions assigned to 'LOF'\n--------------\n","log/LOF_GOF_annotation_stats.txt",append = T)
write.table(LOF_GOF_annotations %>% group_by(annotation) %>% summarise(count=n()), "log/LOF_GOF_annotation_stats.txt", append = T,sep = "\t", row.names = F)


write_csv(LOF_GOF_annotations,"data/LOF_GOF_annotations.csv")
