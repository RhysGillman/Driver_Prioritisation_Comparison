# run_sysSVM2.r


ssms_annotated = annotate_ssms(
  vcf_fn, sample = "my_sample", annovar, 
  genome_version = "hg19", 
  gene_aliases_entrez = "annotation_reference_files/gene_aliases_entrez.tsv", 
  hotspots = "annotation_reference_files/tcga_pancancer_hotspots_oncodriveclust.tsv"
)

annotate_ssms = function(
    vcf_fn,                    # File name of somatic VCF to annotate
    sample,                    # Name of this sample (NB can't process multi-sample VCFs)
    annovar_dir,               # Directory where ANNOVAR is installed (should contain table_annovar.pl)
    genome_version = "hg19",   # Version of the human genome to use for ANNOVAR - hg19 or hg38
    gene_aliases_entrez,       # gene_aliases_entrez.tsv from the sysSVM2 GitHub repository
    hotspots,                  # tcga_pancancer_hotspots_oncodriveclust.tsv from the sysSVM2 GitHub repository
    temp_dir = tempdir()       # Directory for temporary files to be created
){
  
  # Packages
  require(readr)
  require(tidyr)
  require(dplyr)
  
  
  # Run ANNOVAR on VCF
  annovar_output_fn = tempfile(tmpdir = temp_dir)
  annovar_cmd = paste0(
    "perl ", annovar_dir, "/table_annovar.pl --vcfinput ", vcf_fn, 
    " ", annovar_dir, "/humandb/ --buildver ", genome_version, " --outfile ", annovar_output_fn, " --remove ", 
    "--protocol refGene,dbnsfp35a,dbscsnv11 ", 
    "--operation g,f,f"
  )
  message("Running ANNOVAR...")
  system(annovar_cmd)
  message("Done")
  
  
  # Load annotated file and delete temporary files
  annovar_output = suppressWarnings(read_tsv(paste0(annovar_output_fn, ".", genome_version, "_multianno.txt"), col_types = cols()))
  cleanup_cmd = paste0("rm ", annovar_output_fn, "*")
  system(cleanup_cmd)
  
  
  # Clean up the annovar output
  if ("GERP++_RS" %in% colnames(annovar_output)){
    annovar_clean = annovar_output %>% rename(GERPpp_RS = `GERP++_RS`)
  } else if ("GERP.._RS" %in% colnames(annovar_output)){
    annovar_clean = annovar_output %>% rename(GERPpp_RS = `GERP.._RS`)
  } else {
    stop("GERP++ annotations missing")
  }
  annovar_clean = annovar_clean %>%
    select(
      Chr:Alt, 
      Func.refGene, Gene.refGene, ExonicFunc.refGene, AAChange.refGene,
      SIFT_pred, LRT_pred, FATHMM_pred, MutationTaster_pred, Polyphen2_HDIV_pred, Polyphen2_HVAR_pred, MutationAssessor_pred,
      phyloP100way_vertebrate, GERPpp_RS, SiPhy_29way_logOdds,
      dbscSNV_ADA_SCORE, dbscSNV_RF_SCORE
    ) %>%
    mutate_at(
      c("phyloP100way_vertebrate", "GERPpp_RS", "SiPhy_29way_logOdds", "dbscSNV_ADA_SCORE", "dbscSNV_RF_SCORE"),
      function(x) suppressWarnings(as.numeric(x))
    ) %>%
    mutate(Gene.refGene = strsplit(Gene.refGene, split = ";")) %>%
    tidyr::unnest(Gene.refGene) %>%
    unique
  
  
  # Get exonic/splicing variants within a consistently-named set of 19549 human genes
  annovar_exonic = annovar_clean %>% 
    subset(Func.refGene == "exonic" | Func.refGene %in% c("splicing", "exonic;splicing"))
  if (is.character(gene_aliases_entrez)) gene_aliases_entrez = read_tsv(gene_aliases_entrez, col_types = cols())
  gene_aliases_entrez = gene_aliases_entrez %>% select(symbol, alias, entrez)
  annovar_exonic = inner_join(
    annovar_exonic,
    gene_aliases_entrez,
    by = c("Gene.refGene" = "alias")
  )
  
  
  # Parse truncating and predicted damaging missense/splicing mutations
  damaging_trunc_ntdam = annovar_exonic %>%
    # Truncating
    mutate(dam_trunc = ExonicFunc.refGene %in% c("stopgain", "stoploss", "frameshift deletion", "frameshift insertion", "frameshift substitution")) %>%
    mutate(dam_trunc = replace_na(dam_trunc, F)) %>%
    # Functional 5/7
    mutate(
      dam_sift = SIFT_pred == "D",
      dam_lrt = LRT_pred == "D",
      dam_fathmm = FATHMM_pred == "D",
      dam_mutationtaster = MutationTaster_pred %in% c("D", "A"),
      dam_polyphen_hdiv = Polyphen2_HDIV_pred %in% c("D", "P"),
      dam_polyphen_hvar = Polyphen2_HVAR_pred %in% c("D", "P"),
      dam_mutationassessor = MutationAssessor_pred %in% c("H", "M")
    ) %>%
    mutate_at(c("dam_sift", "dam_lrt", "dam_fathmm", "dam_mutationtaster", "dam_polyphen_hdiv", "dam_polyphen_hvar", "dam_mutationassessor"),
              function(x) replace_na(x, F)) %>%
    mutate(dam_func = !is.na(ExonicFunc.refGene) & ExonicFunc.refGene == "nonsynonymous SNV" &
             dam_sift + dam_lrt + dam_fathmm + dam_mutationtaster + dam_polyphen_hdiv + dam_polyphen_hvar + dam_mutationassessor >= 5) %>%
    # Conservation 2/3
    mutate(
      dam_phylop = phyloP100way_vertebrate > 1.6,
      dam_gerp = GERPpp_RS > 4.4,
      dam_siphy = SiPhy_29way_logOdds > 12.17
    ) %>%
    mutate_at(c("dam_phylop", "dam_gerp", "dam_siphy"), function(x) replace_na(x, F)) %>%
    mutate(dam_cons = !is.na(ExonicFunc.refGene) & ExonicFunc.refGene == "nonsynonymous SNV" &
             dam_phylop + dam_gerp + dam_siphy >= 2) %>%
    # Splicing 1/2
    mutate(
      dam_ada = dbscSNV_ADA_SCORE > 0.6,
      dam_rf = dbscSNV_RF_SCORE > 0.6
    ) %>%
    mutate_at(c("dam_ada", "dam_rf"), function(x) replace_na(x, F)) %>%
    mutate(dam_splicing = dam_ada + dam_rf >= 1 & Func.refGene %in% c("exonic;splicing", "splicing"))
  
  
  # Map hotpots
  if (is.character(hotspots)) hotspots = read_tsv(hotspots, col_types = cols())
  damaging_trunc_ntdam_gof = map_hotspots(damaging_trunc_ntdam, hotspots)
  
  
  # Clean table
  damaging_ssms = damaging_trunc_ntdam_gof %>%
    unite(
      "variant_id", c(Chr, Start, End, Ref, Alt), sep = ";"
    ) %>%
    mutate(
      TRUNC_mut = dam_trunc,
      NTDam_mut = dam_func | dam_cons | dam_splicing,
      GOF_mut = is_hotspot,
      damaging = TRUNC_mut | NTDam_mut | GOF_mut,
      sample = sample
    ) %>%
    select(
      sample, variant_id, symbol, entrez, Func.refGene, ExonicFunc.refGene,
      TRUNC_mut, NTDam_mut, GOF_mut, damaging
    )
  
  
  # Output
  return(damaging_ssms)
}
