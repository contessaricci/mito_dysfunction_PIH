# mito_dysfunction_PIH

##############################################################################
###                     
### AUTHOR: CONTESSA A. RICCI, PhD ###
### MANUSCRIPT: DYSREGULATION OF MITOCHONDRIA-MEDIATED MATERNAL-FETAL
###             INTERACTIONS IN HYPERTENSIVE DISORDERS OF PREGNANCY
###             DOI:
### STUDY PURPOSE: Reanalysis of longitudinal maternal RNAseq data from
###                peripheral blood plasma and endpoint fetal RNAseq data
###                from placenta at delivery. Goal is to understand
###                consequences of mitochondrial gene dysregulation by
###                examining expression patterns of genes the dysregulated
###                mitochondrial genes are known to interact with
### DATA: GEO accession GSE154377 (maternal longitudinal data)
###       GEO accession GSE114691 (fetal placental)
###       Accessible via NCBI GEO datasets
###
##############################################################################

Each script explains the process it is meant to carry out. This read me
outlines the computational process carried out for this study


1) Obtain DEGs from maternal and fetal samples (DESeq2 analysis, mito_dysregulation_Ricci2021.r, R)
2) Identify mtDEGs from maternal and fetal samples (find_mtDEGs.py, Python)
3) Identify interaction genes present in maternal and fetal samples (find_mtDEG-INTXs.py, Python)
4) Write GMT files (actual and null model testing) for TcGSA (find_mtDEG-INTXs.py, Python)
5) Carry out TcGSA (TcGSA, mito_dysregulation_Ricci2021.r, R)
6) Carry out Ingenuity Pathway Analysis and Gene Ontology functional enrichments (externally)
