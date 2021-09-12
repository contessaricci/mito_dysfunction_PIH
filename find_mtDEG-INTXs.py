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
### NOTES: All scripts assumes all data has been compiled and is in correct
###        format. This script utilizes some files that have been compiled
###        either manually, in Python, or in R (or a combination thereof).
###
##############################################################################

# PURPOSE OF SCRIPT: 1) FIND mtDEGs FROM MATERNAL DEG FILES AND FETAL GENE LIST PROVIDED BY AUTHORS
#					 2) MAKE GENE MATRIX TRANSPOSED (GMT) FILE FOR TcGSA (carried out in R)
#					 3) MAKE NULL MODEL GMT FILE FOR TcGSA (carried out in R)
#					 4) CREATE FILES OF mtDEG INTERACTION GENES (mtDEG-INTXs) FOR MATERNAL AND
#						FETAL AT DELIVERY
#
# notes:
#		-	There are print statements during the mtDEG-INTXs searching in String database to mark progress
#		-	will use files generated by find_mtDEGs.py script:
#												- maternal_gest_mtDEG_IDs.txt (tab delimited)
#												- maternal_del_mtDEG_IDs.txt (tab delimited)
#												- fetal_mtDEG_IDs.txt (tab delimited)
#												- maternal_gest_DEGs.txt (tab delimited)
#												- maternal_del_DEGs.txt (tab delimited)
#												- fetal_del_DEGs.txt (tab delimited)
#		-	will use delivery DEG files for maternal and fetal samples:
#												- del_sigContr.csv (comma delimited)
#												- fetal_DEGs.csv (comma delimited)
#		-	will use an all contrasts file (does not need to be a specific one) for GMT files generation:
#												- T1_Contr.csv (comma delimited)
#		-	will use Ensembl BioMart files (generated using BioMart feature at ensemble.org):
#												- maternal_gestation_BioMart.txt (tab delimited)
#												- maternal_gestation_all_BioMart.txt (tab delimited)
#												- maternal_delivery_BioMart.txt (tab delimited)
#												- fetal_delivery_BioMart.txt (tab delimited)
#				note: need to use this because String database uses Ensembl protein IDs. Will
#					  generate three files - one for all DEGs during maternal gestation, one for all
#					  DEGs during maternal at delivery, one for all DEGs during fetal at delivery. Files
#					  for Ensemble BioMart searching were created by find_mtDEGs.py script (maternal_gest_DEGs.txt,
#					  maternal_del_DEGs.txt, fetal_del_DEGs.txt). Options selected for Ensembl BioMart file
#					  generation were as follows:
#												- Dataset: "Ensemble Genes 104"; "Human genes (GRCh38.p13)"
#												- Filters: Under "GENE" select "Input external references ID list"
#														   and select "Gene stable ID(s)" for maternal files or
#														   "UniProtKB/Swiss-Prot ID(s)" for fetal file. Then click
#														   on "Choose File" and upload respective file 
#														   (maternal_gest_DEGs.txt, maternal_del_DEGs.txt,fetal_del_DEGs.txt)
#												- Attributes: Under "GENE" select:
#																				- Gene stable ID
#																				- Protein stable ID
#																				- Gene name (optional)
#															  Under "EXTERNAL" scroll to "External References" and select:
#																				- UniProtKB/Swiss-Prot ID
#					  Click "Results" at the upper left corner and export results as a tab delimited file (TSV).
#					  Rename BioMart file as necessary (maternal_gestation_BioMart.txt, maternal_delivery_BioMart.txt,
#					  fetal_delivery_BioMart.txt).
#				note: will also need to used Ensembl BioMart file to convert mtDEG-INTX Ensembl protein IDs back to
#					  Ensembl gene IDs (maternal samples) or Uniprot IDs (fetal samples) for downstream analysis using
#					  orignial datasets
#		-	will use the Homo sapiens-specific String database text file (downloaded on 3/2/2021):
#												- 9606.protein.links.v11.0.txt (single space delimited) 
#				note: a "medium" interaction score is considered 0.4 or better. Text file multiplies interaction
#					  score by 10^3, so 0.4 == 400 in text file. Therefore, will pull out interactions that score
#					  greater than 400.
#		-	will search for maternal mtDEGs using gene identifiers in original files (Ensemble gene ID)
#				because MitoCarta3.0 has this ID listed. Will search for fetal mtDEGs using Uniprot IDs
#				because MitoCarta3.0 also has this ID listed.
#		-	GMT file will be mtDEG-INTXs, where the gene set name is a mtDEG and the gene set is comprised of its
#				corresponding interaction genes (genes that meet an interaction score of "medium" or greater with
#				the mtDEG that the gene set is named for).
#		-	Null model GMT file will be a size-matched GMT file comprised of non-mtDEG, non-DEG genes
#				that are randomly selected. The null model GMT file included on GitHub was created on 4/28/2021.
#				New null GMT file will produce slightly different results.
#
# OVERVIEW OF PROCESS:
#		1) FIND mtDEG-INTXs MEETING "MEDIUM" SCORE OR BETTER (> 400)
#		2) MATCH mtDEG-INTXs TO DEGs FOR: MATERNAL GESTATION, MATERNAL AT DELIVERY, FETAL AT DELIVERY
#		3) MAKE GMT FILES (ACTUAL GMT FILE AND NULL MODEL GMT FILE)
#		4) MAKE FILE OF CONTRAST RESULTS FOR MATERNAL AT DELIVERY mtDEG-INTXs
#		5) MAKE FILE OF CONTRAST RESULTS FOR FETAL AT DELIVERY mtDEG-INTXs



###################################################
### 		###### READ IN FILES ######			###
###################################################

## READ IN mtDEG IDs ##
maternal_gest_mtDEG_IDs = open("maternal_gest_mtDEG_IDs.txt", "r").read().split("\n")
maternal_del_mtDEG_IDs = open("maternal_del_mtDEG_IDs.txt", "r").read().split("\n")
fetal_mtDEG_IDs = open("fetal_mtDEG_IDs.txt", "r").read().split("\n")

## READ IN BIOMART FILES ##
# [0] = Gene stable ID
# [1] = Protein stable ID
# [2] = Gene name
# [3] = UniProtKB/Swiss-Prot ID
maternal_gestation_BioMart = open("maternal_gestation_BioMart.txt", "r").read().split("\n")

maternal_delivery_BioMart = open("maternal_delivery_BioMart.txt", "r").read().split("\n")
fetal_delivery_BioMart = open("fetal_delivery_BioMart.txt", "r").read().split("\n")

## READ IN STRING DB ##
string_db = open("9606.protein.links.v11.0.txt", "r").read().split("\n")


###############################################################
###			###### INITIALIZE DICTIONARIES ######			###
###############################################################

# Format dictionaries for: maternal gestation interactions, maternal delivery interactions, fetal delivery interactions
# Initialize with relative mtDEGs 

maternal_gestation_interactions = {}
maternal_delivery_interactions = {}
fetal_delivery_interactions = {}

## INITIALIZE DICTIONARIES WITH mtDEGs AND MATCHING ENSEMBL PROTEIN IDs ##
# note: need to collect homologous protein IDs because some IDs match multiple Ensembl protein IDs.
#		Will be populating the "Intxs_prots" list with all mtDEG-INTXs present in String database, and
#		will later populate the "Intxs_genes" list with Ensembl gene IDs (maternal) or Uniprot IDs (fetal)
#		that match the Ensembl protein IDs and are also found in BioMart file for DEGs present in samples
for mtDEG in maternal_gest_mtDEG_IDs:
	if len(mtDEG) > 1:
		maternal_gestation_interactions[mtDEG] = {"Prot_IDs":[], "Intxs_prots":[], "Intxs_genes":[]}
for entry in maternal_gestation_BioMart:
	for key, value in maternal_gestation_interactions.items():
		if key in entry and entry.split("\t")[1] not in value["Prot_IDs"]:
			value["Prot_IDs"].append(entry.split("\t")[1])

for mtDEG in maternal_del_mtDEG_IDs:
	if len(mtDEG) > 1:
		maternal_delivery_interactions[mtDEG] = {"Prot_IDs":[], "Intxs_prots":[], "Intxs_genes":[]}
for entry in maternal_delivery_BioMart:
	for key, value in maternal_delivery_interactions.items():
		if key in entry and entry.split("\t")[1] not in value["Prot_IDs"]:
			value["Prot_IDs"].append(entry.split("\t")[1])
	
for mtDEG in fetal_mtDEG_IDs:
	if len(mtDEG) > 1:
		fetal_delivery_interactions[mtDEG] = {"Prot_IDs":[], "Intxs_prots":[], "Intxs_genes":[]}
for entry in fetal_delivery_BioMart:
	for key, value in fetal_delivery_interactions.items():
		if key in entry and entry.split("\t")[1] not in value["Prot_IDs"]:
			value["Prot_IDs"].append(entry.split("\t")[1])
			
			
###########################################################################
###			###### SEARCH FOR mtDEG-INTXs IN STRING DB ######			###
###########################################################################
# note: this step takes quite a bit of time. Recommend carrying out on a server. However, this
#		script may need to be modified to accomodate this (depending on your server architecture).
#		Print statements are provided to keep track of progress
# note: excludes proteins that are homologous to mtDEG
# note: formatting of protein ID is in format 9606.XXX, where XXX is Ensembl protein ID
print "Begin searching for mtDEG-INTXs in String database. This may take a while..."

## SEARCH FOR mtDEG-INTXs IN STRING DATABASE ##
for line in string_db[1:]:	
	for key, value in maternal_gestation_interactions.items():
		for prot in value["Prot_IDs"]:
			if len(line) > 1 and int(line.split(" ")[-1]) > 400:
				protein1 = line.split(" ")[0].split(".")[1]
				protein2 = line.split(" ")[1].split(".")[1]
				if protein1 == prot and protein2 not in value["Prot_IDs"]:
					if protein2 not in value["Intxs_prots"]:
						value["Intxs_prots"].append(protein2)
				elif protein2 == prot and protein1 not in value["Prot_IDs"]:
					if protein1 not in value["Intxs_prots"]:
						value["Intxs_prots"].append(protein1)
print "All maternal gestation mtDEGs-INTXs found in String DB. Now working on maternal at delivery..."


for line in string_db[1:]:	
	for key, value in maternal_delivery_interactions.items():
		for prot in value["Prot_IDs"]:
			if len(line) > 1 and int(line.split(" ")[-1]) > 400:
				protein1 = line.split(" ")[0].split(".")[1]
				protein2 = line.split(" ")[1].split(".")[1]
				if protein1 == prot and protein2 not in value["Prot_IDs"]:
					if protein2 not in value["Intxs_prots"]:
						value["Intxs_prots"].append(protein2)
				elif protein2 == prot and protein1 not in value["Prot_IDs"]:
					if protein1 not in value["Intxs_prots"]:
						value["Intxs_prots"].append(protein1)
print "All maternal at delivery mtDEGs-INTXs found in String DB. Now working on fetal at delivery..."


for line in string_db[1:]:	
	for key, value in fetal_delivery_interactions.items():
		for prot in value["Prot_IDs"]:
			if len(line) > 1 and int(line.split(" ")[-1]) > 400:
				protein1 = line.split(" ")[0].split(".")[1]
				protein2 = line.split(" ")[1].split(".")[1]
				if protein1 == prot and protein2 not in value["Prot_IDs"]:
					if protein2 not in value["Intxs_prots"]:
						value["Intxs_prots"].append(protein2)
				elif protein2 == prot and protein1 not in value["Prot_IDs"]:
					if protein1 not in value["Intxs_prots"]:
						value["Intxs_prots"].append(protein1)
print "All fetal at delivery mtDEGs-INTXs found in String DB. Way to hang in there."
print "Almost done now. Now matching protein IDs to sample DEGs and writing output files..."


###########################################################################
###		  ###### SEARCH FOR MATCHING IDs IN BIOMART FILE ######			###
###					###### WRITE OUTPUT FILES ######					###
###########################################################################

## MAKE GMT FILES FOR TcGSA (carried out in R) ##
# note: will make actual GMT file and null model GMT file. Gene sets for both will
#		be comprised of genes sequenced from gestation. Not DEG-specific. Acutal GMT
#		file will be all of the mtDEGs found during gestation (length = 14 gene sets).
#		Each gene set will be comprised of respective mtDEG interaction genes identified
#		in String database and matched to all genes sequenced from maternal samples. 
#		Null model GMT file will be size-matched (length = 14 gene sets) with non-mtDEG,
#		non-DEG as geneset names, and populated with random non-DEG genes chosen from
#		maternal samples expression data.
				
maternal_gestation_all_BioMart = open("maternal_gestation_all_BioMart.txt", "r").read().split("\n")
for key, value in maternal_gestation_interactions.items():
	for prot in value["Intxs_prots"]:
		for entry in maternal_gestation_all_BioMart:
			if prot in entry and prot not in value["Prot_IDs"]:
				if entry.split("\t")[0] not in value["Intxs_genes"]:
					value["Intxs_genes"].append(entry.split("\t")[0])

## Write actual GMT file for TcGSA ##
gestation_GMT = open("gestation_GMT.txt", "w")
print "Writing actual GMT file"
for key, value in maternal_gestation_interactions.items():
	#for gene in value["Intxs_genes"]:
	gestation_GMT.write(key + "\t" + "\t" + ("\t").join(value["Intxs_genes"]) + "\n")

## Write null GMT file for TcGSA ##
# Find all non-DEGs
all_contrasts = open("T1_Contr.csv", "r").read().replace("\"","").split("\n")
all_geneIDs = []
for line in all_contrasts:
	if len(line) > 1 and line.split(",")[0] not in all_geneIDs:
		all_geneIDs.append(line.split(",")[0])
maternal_gestation_DEGs = open("maternal_gest_DEGs.txt", "r").read().split("\n")
non_DEGs = []
for ID in all_geneIDs:
	if len(ID) > 1 and ID not in maternal_gestation_DEGs:
		if ID not in non_DEGs:
			non_DEGs.append(ID)

# Make null GMT dictionary and populate with random non-DEG genes
from random import *
gestation_null_GMT = open("null_GMT.txt", "w")
gestation_null_GMT_dict = {}
dict_length = len(maternal_gestation_interactions)
while dict_length != 0:
    random_number = randrange(len(non_DEGs))
    for count, ID in enumerate(non_DEGs):
        if random_number == count and len(ID) > 1:
            if ID not in gestation_null_GMT_dict:
                gestation_null_GMT_dict[ID] = {"Gene_count" : 0, "Random_genes" : []}
                dict_length -= 1
# Size-matched gene sets
lengths = []
null_keys = []
for key, value in maternal_gestation_interactions.items():
	lengths.append(len(value["Intxs_genes"]))
for key in gestation_null_GMT_dict.keys():
	null_keys.append(key)

for index, length in enumerate(lengths):
	for count, entry in enumerate(null_keys):
		for key, value in gestation_null_GMT_dict.items():
			if index == count and entry == key:
				value["Gene_count"] = length

for key, value in gestation_null_GMT_dict.items():
	length = value["Gene_count"]
	while length != 0:
		random_number = randrange(len(non_DEGs))
		for count, ID in enumerate(non_DEGs):
			if random_number == count and ID not in value["Random_genes"]:
				value["Random_genes"].append(ID)
				length -= 1

print "Writing null GMT file"
for key, value in gestation_null_GMT_dict.items():
	gestation_null_GMT.write(key + "\t" + "\t" + "\t".join(value["Random_genes"]) + "\n")
	

## MAKE EXPRESSION FILES FOR MATERNAL AND FETAL AT DELIVERY mtDEG-INTXs ##
## MATERNAL AT DELIVERY ##
for key, value in maternal_delivery_interactions.items():
	for prot in value["Intxs_prots"]:
		for entry in maternal_delivery_BioMart:
			if prot in entry and prot not in value["Prot_IDs"]:
				if entry.split("\t")[0] not in value["Intxs_genes"]:
					value["Intxs_genes"].append(entry.split("\t")[0])

# Pool all maternal at delivery mtDEG-INTXs
maternal_delivery_INTXs = []
for value in maternal_delivery_interactions.values():
	for gene in value["Intxs_genes"]:
		if gene not in maternal_delivery_INTXs:
			maternal_delivery_INTXs.append(gene)

# Make file for maternal at delivery expression results (for Gene Ontology enrichments)
maternal_delivery_DEGs = open("del_sigContr.csv", "r").read().replace("\"","").split("\n")
maternal_delivery_INTXs_DEGs = []
for ID in maternal_delivery_INTXs:
	for line in maternal_delivery_DEGs:
		if ID in line and line not in maternal_delivery_INTXs_DEGs:
			maternal_delivery_INTXs_DEGs.append(line)

print "Writing maternal at delivery file..."
maternal_delivery_INTXs_expression = open("maternal_delivery_mtDEG-INTXs_expression.csv", "w")
maternal_delivery_INTXs_expression.write(maternal_delivery_DEGs[0] + "\n")
for line in maternal_delivery_INTXs_DEGs:
	maternal_delivery_INTXs_expression.write(line + "\n")


## FETAL AT DELIVERY ##
for key, value in fetal_delivery_interactions.items():
	for prot in value["Intxs_prots"]:
		for entry in fetal_delivery_BioMart:
			if prot in entry and prot not in value["Prot_IDs"]:
				if entry.split("\t")[3] not in value["Intxs_genes"]:
					value["Intxs_genes"].append(entry.split("\t")[3])

# Pool all fetal at delivery mtDEG-INTXs
fetal_delivery_INTXs = []
for value in fetal_delivery_interactions.values():
	for gene in value["Intxs_genes"]:
		if gene not in fetal_delivery_INTXs:
			fetal_delivery_INTXs.append(gene)

# Make file for fetal at delivery contrast results (for Gene Ontology enrichments)
print ("Writing fetal at delivery file...")
fetal_delivery_DEGs = open("fetal_DEGs.csv", "r").read().replace("\"","").split("\n")
fetal_delivery_INTXs_DEGs = []
for ID in fetal_delivery_INTXs:
	for line in fetal_delivery_DEGs:
		if ID in line and line not in fetal_delivery_INTXs_DEGs:
			fetal_delivery_INTXs_DEGs.append(line)
fetal_delivery_INTXs_expression = open("fetal_delivery_mtDEG-INTXs_expression.csv", "w")
fetal_delivery_INTXs_expression.write(fetal_delivery_DEGs[0] + "\n")
for line in fetal_delivery_INTXs_DEGs:
	fetal_delivery_INTXs_expression.write(line + "\n")

print "Script complete."






















