##############################################################################
###                     												   ###
### AUTHOR: CONTESSA A. RICCI, PhD 										   ###
### MANUSCRIPT: DYSREGULATION OF MITOCHONDRIA-MEDIATED MATERNAL-FETAL	   ###
###             INTERACTIONS IN HYPERTENSIVE DISORDERS OF PREGNANCY		   ###
###             DOI: (TBD)												   ###
### STUDY PURPOSE: Reanalysis of longitudinal maternal RNAseq data from    ###
###                peripheral blood plasma and endpoint fetal RNAseq data  ###
###                from placenta at delivery. Goal is to understand        ###
###                consequences of mitochondrial gene dysregulation by     ###
###                examining expression patterns of genes the dysregulated ###
###                mitochondrial genes are known to interact with		   ###
### DATA: GEO accession GSE154377 (maternal longitudinal data)			   ###
###       GEO accession GSE114691 (fetal placental)						   ###
###       Accessible via NCBI GEO datasets								   ###
###                     												   ###
### NOTES: All scripts assume all data has been compiled and is in correct ###
###        format. This script utilizes some files that have been compiled ###
###        either manually, in Python, or in R (or a combination thereof). ###
###                     												   ###
##############################################################################

# PURPOSE OF SCRIPT: 1) IDENTIFY mtDEGs FROM MATERNAL DEG FILES AND FETAL GENE LIST PROVIDED BY AUTHORS
#					 2) CREATE FILE FOR FINDING INTERACTION GENES IN STRING DATABASE
#					 3) CREATE FILES FOR ENSEMBLE BIOMART SEARCHING
# notes:
#		-	will use files for significant maternal contrasts:
#												- T1_sigContr.csv (comma delimited)
#												- T2_sigContr.csv (comma delimited)
#												- T3_sigContr.csv (comma delimited)
#												- del_sigContr.csv (comma delimited)
#		-	will use file of fetal DEG list provided by authors: 
#												- fetal_DEGs.csv (comma delimited)
#				note: only gene names were provided by authors. Therefore, gene names provided
#					  were manually searched for in Uniprot database to obtain Uniprot ID.
#		-	will use MitoCarta3.0 database (downloaded on 3/4/2021):
#												- Human.MitoCarta3.0.txt (tab delimited)
#		-	will search for maternal mtDEGs using gene identifiers in original files (Ensemble gene ID)
#				because MitoCarta3.0 has this ID listed. Will search for fetal mtDEGs using Uniprot IDs
#				because MitoCarta3.0 also has this ID listed.

#
# OVERVIEW OF PROCESS:
#		1) PULL OUT mtDEGs FROM EACH MATERNAL CONTRAST (SEARCH USING ENSEMBLE GENE ID) 
#		2) PULL OUT mtDEGs FROM FETAL DEG LIST (SEARCH USING UNIPROT ID)
#		3) POOL ALL MATERNAL mtDEGs AND OUTPUT LIST FOR SEARCHING STRING DATABASE
#		4) CREATE OUTPUT FILE FOR MATERNAL mtDEGs WITH CONTRAST RESULTS AND COLLECTION
#		   TIMEPOINT RECORDED
#		5) CREATE OUTPUT FILE FOR FETAL mtDEGs WITH P-VALUE AND DIRECTION INFO (ORIGINAL PUB)
#		6) OUTPUT LIST OF FETAL mtDEG LIST FOR SEARCHING STRING DATABASE




###################################################
### 		###### READ IN FILES ######			###
###################################################

## MITO CARTA DATABASE ##
# tab delimited file
# [14] = Ensemble gene ID *some entries have multiple IDs separated by "|" 
#						   -- is OK becase using "in" command
# [15] = Uniprot ID
MitoCarta = open("Human.MitoCarta3.0.txt", "r").read().split("\r\n")

## MATERNAL FILES ##
# comma delimited files
# [0] = Ensemble gene ID
T1_DEGs = open("T1_sigContr.csv","r").read().replace("\"","").split("\n")
T2_DEGs = open("T2_sigContr.csv","r").read().replace("\"","").split("\n")
T3_DEGs = open("T3_sigContr.csv","r").read().replace("\"","").split("\n")
del_DEGs = open("del_sigContr.csv","r").read().replace("\"","").split("\n")

## FETAL FILE ##
# comma delimited file
# [5] = Uniprot ID
fetal_DEGs = open("fetal_DEGs.csv", "r").read().split("\r\n")


###############################################################
###			###### INITIALIZE DICTIONARIES ######			###
###		  	 ###### POPULATE DICTIONARIES ######			###
###############################################################

MitoCarta_dict = {"Esemble_IDs":[], "Uniprot_IDs":[]}
maternal_dict = {"T1":[], "T2":[], "T3":[], "del":[]}
fetal_dict = {}

## POPULATE MITOCARTA DICTIONARY ##
for entry in MitoCarta[1:]:
	if len(entry.split("\t")[14]) > 1 and entry.split("\t")[14] not in MitoCarta_dict["Esemble_IDs"]:
		MitoCarta_dict["Esemble_IDs"].append(entry.split("\t")[14])
for entry in MitoCarta[1:]:
	if len(entry.split("\t")[15]) > 1 and entry.split("\t")[15] not in MitoCarta_dict["Uniprot_IDs"]:
		MitoCarta_dict["Uniprot_IDs"].append(entry.split("\t")[15])

print(MitoCarta_dict)
'''
'''
## POPULATE MATERNAL DICTIONARY ##
for entry in T1_DEGs[1:]:
	if entry.split(",")[0] in MitoCarta_dict["Esemble_IDs"] and entry.split(",")[0] not in maternal_dict["T1"]:
		maternal_dict["T1"].append(entry)
for entry in T2_DEGs[1:]:
	if entry.split(",")[0] in MitoCarta_dict["Esemble_IDs"] and entry.split(",")[0] not in maternal_dict["T2"]:
		maternal_dict["T2"].append(entry)
for entry in T3_DEGs[1:]:
	if entry.split(",")[0] in MitoCarta_dict["Esemble_IDs"] and entry.split(",")[0] not in maternal_dict["T3"]:
		maternal_dict["T3"].append(entry)
for entry in del_DEGs[1:]:
	if entry.split(",")[0] in MitoCarta_dict["Esemble_IDs"] and entry.split(",")[0] not in maternal_dict["del"]:
		maternal_dict["del"].append(entry)

## POPULATE FETAL DICTIONARY ##
for entry in fetal_DEGs[1:]:
	if entry.split(",")[5] in MitoCarta_dict["Uniprot_IDs"]:
		fetal_dict[entry] = ""
		

###########################################################
###			###### POOL MATERNAL mtDEGs ######			###
###			 ###### WRITE OUTPUT FILES ######			###
###########################################################

## POOL MATERNAL GESTATION mtDEG ENSEMBLE IDs ##
maternal_gest_mtDEGs = []
for entry in maternal_dict["T1"]:
	if entry.split(",")[0] not in maternal_gest_mtDEGs:
		maternal_gest_mtDEGs.append(entry.split(",")[0])
for entry in maternal_dict["T2"]:
	if entry.split(",")[0] not in maternal_gest_mtDEGs:
		maternal_gest_mtDEGs.append(entry.split(",")[0])
for entry in maternal_dict["T3"]:
	if entry.split(",")[0] not in maternal_gest_mtDEGs:
		maternal_gest_mtDEGs.append(entry.split(",")[0])

## CREATE LIST OF FETAL mtDEG UNIPROT IDs ##
fetal_mtDEGs = []
for key, value in fetal_dict.items():
	fetal_mtDEGs.append(key.split(",")[5])
	
## CREATE LIST OF MATERNAL GESTATION DEGs FOR ENSEMBLE BIOMART SEARCHING ##
maternal_gest = []
for entry in T1_DEGs:
	if entry.split(",")[0] not in maternal_gest:
		maternal_gest.append(entry.split(",")[0])
for entry in T2_DEGs:
	if entry.split(",")[0] not in maternal_gest:
		maternal_gest.append(entry.split(",")[0])
for entry in T3_DEGs:
	if entry.split(",")[0] not in maternal_gest:
		maternal_gest.append(entry.split(",")[0])

## CREATE LIST OF MATERNAL AT DELIVERY DEGs FOR ENSEMBLE BIOMART SEARCHING ##
maternal_del = []
for entry in del_DEGs:
	if entry.split(",")[0] not in maternal_del:
		maternal_del.append(entry.split(",")[0])

## CREATE LIST OF FETAL AT DELIVERY DEGs FOR ENSEMBLE BIOMART SEARCHING ##
fetal_del = []
for entry in fetal_DEGs:
	if entry.split(",")[5] not in fetal_del:
		fetal_del.append(entry.split(",")[5])



### WRITE OUTPUT FILES ###### WRITE OUTPUT FILES ###### WRITE OUTPUT FILES ###

## MATERNAL PLOTTING OUTPUT FILE ##
maternal_mtDEG_expression = open("maternal_mtDEG_expression.csv", "w")
maternal_mtDEG_expression.write(T1_DEGs[0] + "," + "Collection" + "\n")
for key, value in maternal_dict.items():
	for entry in value:
		maternal_mtDEG_expression.write(entry + "," + key + "\n")

## FETAL PLOTTING OUTPUT FILE ##
fetal_mtDEG_expression = open("fetal_mtDEG_expression.csv", "w")
fetal_mtDEG_expression.write(fetal_DEGs[0] + "\n")
for key, value in fetal_dict.items():
	fetal_mtDEG_expression.write(key + "\n")

## MATERNAL GESTATION mtDEG FILE FOR SEARCHING STRING DATABASE ##
maternal_gest_mtDEG_IDs = open("maternal_gest_mtDEG_IDs.txt", "w")
for ID in maternal_gest_mtDEGs:
	maternal_gest_mtDEG_IDs.write(ID + "\n")

## MATERNAL DELIVERY mtDEG FILE FOR SEARCHING STRING DATABASE ##
maternal_del_mtDEGs =[]
for entry in maternal_dict["del"]:
	if entry.split(",")[0] not in maternal_del_mtDEGs:
		maternal_del_mtDEGs.append(entry.split(",")[0])
maternal_del_mtDEG_IDs = open("maternal_del_mtDEG_IDs.txt", "w")
for ID in maternal_del_mtDEGs:
	maternal_del_mtDEG_IDs.write(ID + "\n")

## FETAL mtDEG FILE FOR SEARCHING STRING DATABASE ##
fetal_mtDEG_IDs = open("fetal_mtDEG_IDs.txt", "w")
for ID in fetal_mtDEGs:
	fetal_mtDEG_IDs.write(ID + "\n")

## MATERNAL GESTATION DEGs FILES FOR ENSEMBL BIOMART SEARCHING ##
maternal_gest_DEGs = open("maternal_gest_DEGs.txt", "w")
for ID in maternal_gest:
	maternal_gest_DEGs.write(ID + "\n")

## MATERNAL ALL GENES FOR ENSEMBL BIOMART SEARCHING (for making GMT files) ##
maternal_all_genes = open("T1_Contr.csv", "r").read().replace("\"","").split("\n")
maternal_all_geneIDs = []
for line in maternal_all_genes:
	if len(line) > 1 and line.split(",")[0] not in maternal_all_geneIDs:
		maternal_all_geneIDs.append(line.split(",")[0])
maternal_all_geneIDs_file = open("maternal_all_geneIDs.txt", "w")
for ID in maternal_all_geneIDs:
	maternal_all_geneIDs_file.write(ID + "\n")

## MATERNAL DELIVERY DEGs FILES FOR ENSEMBL BIOMART SEARCHING ##
maternal_del_DEGs = open("maternal_del_DEGs.txt", "w")
for ID in maternal_del:
	maternal_del_DEGs.write(ID + "\n")

## FETAL DEGs FILE FOR ENSEMBL BIOMART SEARCHING ##
fetal_del_DEGs = open("fetal_del_DEGs.txt", "w")
for ID in fetal_del:
	fetal_del_DEGs.write(ID + "\n")
