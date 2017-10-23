#----------
# libraries
#----------


pacman::p_load(jsonlite, httr, parallel, readxl)

loadNamespace('bbmle')
library(broom)

loadNamespace('mygene')
loadNamespace('biomaRt')

library(RMySQL)


#--------
# imports
#--------

source('_keys.R')


#--------------------#
#                    #
#   TISSUE BINNING   #
#                    #
#--------------------#


HPAFormatTissue <- function(tissue.exp) {

	tissue.exp %>%
	mutate(tissue = case_when(
		.$tissue %in% c('adipose tissue') ~ 'adipose tissue',
		.$tissue %in% c('adrenal gland') ~ 'adrenal',
		.$tissue %in% c('appendix') ~ 'appendix',
		.$tissue %in% c('urinary bladder') ~ 'bladder',
		FALSE ~ 'blood',
		.$tissue %in% c('bone marrow') ~ 'bone',
		.$tissue %in% c('cerebellum', 'cerebral cortex', 'hippocampus', 'lateral ventricle') ~ 'brain',
		.$tissue %in% c('breast') ~ 'breast',
		.$tissue %in% c('bronchus') ~ 'bronchus',
		FALSE ~ 'cerumen',
		.$tissue %in% c('cervix, uterine') ~ 'cervix',
		.$tissue %in% c('epididymis') ~ 'epididymis',
		FALSE ~ 'eye',
		.$tissue %in% c('fallopian tube') ~ 'fallopian tube',
		.$tissue %in% c('gallbladder') ~ 'gallbladder',
		.$tissue %in% c('colon', 'duodenum', 'small intestine') ~ 'gut',
		.$tissue %in% c('heart muscle') ~ 'heart',
		.$tissue %in% c('kidney') ~ 'kidney',
		.$tissue %in% c('esophagus') ~ 'esophagus',
		.$tissue %in% c('liver') ~ 'liver',
		.$tissue %in% c('lung') ~ 'lung',
		.$tissue %in% c('lymph node') ~ 'lymph node',
		.$tissue %in% c('nasopharynx') ~ 'nasopharynx',
		.$tissue %in% c('oral mucosa', 'salivary gland') ~ 'oropharynx',
		.$tissue %in% c('ovary') ~ 'ovary',
		.$tissue %in% c('pancreas') ~ 'pancreas',
		.$tissue %in% c('parathyroid gland') ~ 'parathyroid',
		.$tissue %in% c('prostate') ~ 'prostate',
		.$tissue %in% c('rectum') ~ 'rectum',
		.$tissue %in% c('seminal vesicle') ~ 'seminal',
		.$tissue %in% c('skeletal muscle') ~ 'skeletal muscle',
		.$tissue %in% c('skin', 'skin 1', 'skin 2') ~ 'skin',
		.$tissue %in% c('smooth muscle') ~ 'smooth muscle',
		.$tissue %in% c('soft tissue 1', 'soft tissue 2') ~ 'soft tissue',
		FALSE ~ 'spinal cord',
		.$tissue %in% c('spleen') ~ 'spleen',
		.$tissue %in% c('stomach', 'stomach 1', 'stomach 2') ~ 'stomach',
		FALSE ~ 'synovial fluid',
		.$tissue %in% c('testis') ~ 'testis',
		.$tissue %in% c('thyroid gland') ~ 'thyroid',
		.$tissue %in% c('tonsil') ~ 'tonsil',
		.$tissue %in% c('endometrium', 'endometrium 1', 'endometrium 2') ~ 'uterus',
		.$tissue %in% c('vagina') ~ 'vagina',
		TRUE ~ 'missing'))
}


HPMFormatTissue <- function(tissue.exp) {

	tissue.exp %>%
	mutate(tissue = case_when(
		FALSE ~ 'adipose tissue',
		.$tissue %in% c('Adult.Adrenal') ~ 'adrenal',
		FALSE ~ 'appendix',
		.$tissue %in% c('Adult.Urinary.Bladder') ~ 'bladder',
		.$tissue %in% c('B.Cells', 'CD4.Cells', 'CD8.Cells', 'Monocytes', 'NK.Cells', 'Platelets') ~ 'blood',
		FALSE ~ 'bone',
		.$tissue %in% c('Adult.Frontal.Cortex') ~ 'brain',
		FALSE ~ 'breast',
		FALSE ~ 'bronchus',
		FALSE ~ 'cerumen',
		FALSE ~ 'cervix',
		FALSE ~ 'epididymis',
		.$tissue %in% c('Adult.Retina') ~ 'eye',
		FALSE ~ 'fallopian tube',
		.$tissue %in% c('Adult.Gallbladder') ~ 'gallbladder',
		.$tissue %in% c('Adult.Colon') ~ 'gut',
		.$tissue %in% c('Adult.Heart') ~ 'heart',
		.$tissue %in% c('Adult.Kidney') ~ 'kidney',
		.$tissue %in% c('Adult.Esophagus') ~ 'esophagus',
		.$tissue %in% c('Adult.Liver') ~ 'liver',
		.$tissue %in% c('Adult.Lung') ~ 'lung',
		FALSE ~ 'lymph node',
		FALSE ~ 'nasopharynx',
		FALSE ~ 'oropharynx',
		.$tissue %in% c('Adult.Ovary') ~ 'ovary',
		.$tissue %in% c('Adult.Pancreas') ~ 'pancreas',
		FALSE ~ 'parathyroid',
		.$tissue %in% c('Adult.Prostate') ~ 'prostate',
		.$tissue %in% c('Adult.Rectum') ~ 'rectum',
		FALSE ~ 'seminal',
		FALSE ~ 'skeletal muscle',
		FALSE ~ 'skin',
		FALSE ~ 'smooth muscle',
		FALSE ~ 'soft tissue',
		.$tissue %in% c('Adult.Spinal.Cord') ~ 'spinal cord',
		FALSE ~ 'spleen',
		FALSE ~ 'stomach',
		FALSE ~ 'synovial fluid',
		.$tissue %in% c('Adult.Testis') ~ 'testis',
		FALSE ~ 'thyroid',
		FALSE ~ 'tonsil',
		FALSE ~ 'uterus',
		FALSE ~ 'vagina',
		TRUE ~ 'missing'))
}


PDBFormatTissue <- function(tissue.exp) {

	tissue.exp %>%
	mutate(tissue = case_when(
		.$tissue %in% c('adipocyte') ~ 'adipose tissue',
		.$tissue %in% c('adrenal gland') ~ 'adrenal',
		FALSE ~ 'appendix',
		.$tissue %in% c('urinary bladder', 'urine') ~ 'bladder',
		.$tissue %in% c('B-lymphocyte', 'blood', 'blood platelet', 'cytotoxic T-lymphocyte', 'helper T-lymphocyte', 'monocyte', 'natural killer cell', 'serum') ~ 'blood',
		.$tissue %in% c('bone', 'bone marrow stromal cell', 'mesenchymal stem cell') ~ 'bone',
		.$tissue %in% c('brain', 'cerebral cortex', 'prefrontal cortex') ~ 'brain',
		.$tissue %in% c('breast') ~ 'breast',
		FALSE ~ 'bronchus',
		.$tissue %in% c('cerumen') ~ 'cerumen',
		.$tissue %in% c('cervical epithelium', 'cervical mucosa', 'uterine cervix', 'uterus') ~ 'cervix',
		FALSE ~ 'epididymis',
		.$tissue %in% c('retina', 'vitreous humor') ~ 'eye',
		FALSE ~ 'fallopian tube',
		.$tissue %in% c('gall bladder') ~ 'gallbladder',
		.$tissue %in% c('colon', 'colon muscle', 'colonic epithelial cell', 'gut', 'ileum epithelial cell') ~ 'gut',
		.$tissue %in% c('heart', 'proximal fluid (coronary sinus)') ~ 'heart',
		.$tissue %in% c('kidney') ~ 'kidney',
		.$tissue %in% c('esophagus') ~ 'esophagus',
		.$tissue %in% c('bile', 'liver') ~ 'liver',
		.$tissue %in% c('lung') ~ 'lung',
		.$tissue %in% c('lymph node') ~ 'lymph node',
		.$tissue %in% c('nasopharynx') ~ 'nasopharynx',
		.$tissue %in% c('oral epithelium', 'saliva', 'salivary gland') ~ 'oropharynx',
		.$tissue %in% c('ovary') ~ 'ovary',
		.$tissue %in% c('pancreas', 'pancreatic islet', 'pancreatic juice') ~ 'pancreas',
		FALSE ~ 'parathyroid',
		.$tissue %in% c('prostate gland') ~ 'prostate',
		.$tissue %in% c('rectum') ~ 'rectum',
		.$tissue %in% c('seminal plasma', 'seminal vesicle', 'spermatozoon') ~ 'seminal',
		FALSE ~ 'skeletal muscle',
		.$tissue %in% c('hair follicle', 'skin') ~ 'skin',
		FALSE ~ 'smooth muscle',
		FALSE ~ 'soft tissue',
		.$tissue %in% c('cerebrospinal fluid', 'spinal cord') ~ 'spinal cord',
		.$tissue %in% c('spleen') ~ 'spleen',
		.$tissue %in% c('cardia', 'stomach') ~ 'stomach',
		.$tissue %in% c('synovial fluid') ~ 'synovial fluid',
		.$tissue %in% c('testis') ~ 'testis',
		.$tissue %in% c('thyroid gland') ~ 'thyroid',
		.$tissue %in% c('tonsil') ~ 'tonsil',
		.$tissue %in% c('myometrium') ~ 'uterus',
		FALSE ~ 'vagina',
		TRUE ~ 'missing'))
}




#---------#
#         #
#   HPA   #
#         #
#---------#




hpa.normal.cut <-
	read.delim('../data/HPA/normal_tissue.csv', stringsAsFactors=FALSE, sep=',') %>%
	tbl_df %>%
	select(
		ensembl.gene = Gene,
		tissue = Tissue,
		hpa.p = Level) %>%
	mutate(hpa.p =
		ifelse(hpa.p == 'High', 3,
		ifelse(hpa.p == 'Medium', 2,
		ifelse(hpa.p == 'Low', 1,
		0 )))) %>%
	filter(!tissue %in% 'placenta') %>%
	HPAFormatTissue %>%
	group_by(ensembl.gene, tissue) %>%
	mutate(hpa.p = max(hpa.p, na.rm = TRUE)) %>%
	ungroup %>%
	select(ensembl.gene, tissue, hpa.p.exp = hpa.p) %>%
	unique


hpa.rna.normal.cut <-
	read.delim('../data/HPA/rna_tissue.csv', stringsAsFactors=FALSE, sep=',') %>%
	tbl_df %>%
	select(
		ensembl.gene = Gene,
		tissue = Sample,
		hpa.r = Abundance) %>%
	mutate(hpa.r =
		ifelse(hpa.r == 'High', 3,
		ifelse(hpa.r == 'Medium', 2,
		ifelse(hpa.r == 'Low', 1,
		0 )))) %>%
	filter(!tissue %in% 'placenta') %>%
	HPAFormatTissue %>%
	group_by(ensembl.gene, tissue) %>%
	mutate(hpa.r = max(hpa.r, na.rm = TRUE)) %>%
	ungroup %>%
	select(ensembl.gene, tissue, hpa.r.exp = hpa.r) %>%
	unique


hpa.extracellular <-
	c(	# 'Aggresome',
		'Cell Junctions',
		# 'Centrosome',
		# 'Cytoplasm',
		# 'Cytoskeleton (Actin filaments)',
		'Cytoskeleton (Cytokinetic bridge)',
		# 'Cytoskeleton (Intermediate filaments)',
		# 'Cytoskeleton (Microtubule end)',
		# 'Cytoskeleton (Microtubules)',
		# 'Endoplasmic reticulum',
		'Focal Adhesions',
		# 'Golgi apparatus',
		# 'Microtubule organizing center',
		# 'Mitochondria',
		# 'Nuclear membrane',
		# 'Nucleoli',
		# 'Nucleus',
		# 'Nucleus but not nucleoli',
		'Plasma membrane'
		# 'Vesicles'
	)

hpa.location <-
	read.delim('../data/HPA/subcellular_location.csv', stringsAsFactors=FALSE, sep=',') %>%
	tbl_df %>%
	mutate(location = str_c(Main.location, ';', Other.location)) %>%
	select(ensembl.gene = Gene, location) %>%
	separate(location, into = letters[1:5], sep = '\\;') %>%
	gather(letters, location, a:e) %>%
	select(ensembl.gene, location) %>%
	arrange(ensembl.gene, location) %>%
	filter(!is.na(location) & location != '') %>%
	mutate(membrane = location %in% hpa.extracellular) %>%
	filter(membrane == TRUE) %>%
	select(ensembl.gene, location) %>%
	mutate(database = 'HPA')


#---------#
#         #
#   HPM   #
#         #
#---------#


hpm.normal <-
	read.delim('../data/HPM/HPM_gene_level_epxression_matrix_Kim_et_al_052914.csv', stringsAsFactors=FALSE, sep=',') %>%
	tbl_df %>%
	gather(tissue, hpm.p, Fetal.Heart:Platelets) %>%
	rename(hgnc = Gene) %>%
	filter(!tissue %in% c('Fetal.Brain', 'Fetal.Gut', 'Fetal.Heart', 'Fetal.Liver', 'Fetal.Ovary', 'Fetal.Testis', 'Placenta')) %>%
	HPMFormatTissue %>%
	group_by(hgnc, tissue) %>%
	mutate(hpm.p = max(hpm.p, na.rm = TRUE)) %>%
	ungroup %>%
	unique %>%
	mutate(hpm.p.log = log10(hpm.p)) %>%
	mutate(hpm.p.log = ifelse(hpm.p.log == -Inf, 0, hpm.p.log))


#---------#
#         #
#   PDB   #
#         #
#---------#


#----------
# functions
#----------


# get table back from PDB API
PDBQuery <- function(query) {

	response <- GET(query, authenticate(pdb.user, pdb.pass, 'basic'))

	stop_for_status(response)

	table <-
		content(response, as = 'text', encoding = 'UTF-8') %>%
		fromJSON %$% d %$% results

	table
}


# handle tissue request & format for serial processing
GetTissue <- function(tissue.id, tissue) {

	print(str_c(tissue.id, ' -- ', tissue))

	tissue.query <-
		str_c(
			"https://www.proteomicsdb.org/proteomicsdb/logic/api/proteinspertissue.xsodata/InputParams(TISSUE_ID='",
			tissue.id,
			"',CALCULATION_METHOD=0,SWISSPROT_ONLY=1,NO_ISOFORM=0)/Results?$select=ENTRY_NAME,UNIQUE_IDENTIFIER,DATABASE,SAMPLE_NAME,UNNORMALIZED_EXPRESSION,NORMALIZED_EXPRESSION&$format=json"
		)

	tissue.table <- PDBQuery(tissue.query)

	if(tissue.table %>% length > 0) {

		tissue.table %<>% set_names(c('metadata', 'entry.name', 'uniprot', 'database', 'sample.name', 'unnormalized.expression', 'normalized.expression'))

		tissue.uri <- tissue.table[,1]['uri'] %>% unlist %>% unname

		tissue <-
			tissue.table[,-1] %>%
			tbl_df %>%
			mutate(uri = tissue.uri) %>%
			mutate(
				tissue.id = tissue.id,
				tissue = tissue
				) %>%
			select(tissue.id, tissue, everything())

		tissue

	} else {
		data_frame(`metadata`='', `tissue.id`, `tissue`='', `entry.name`='', `uniprot`='', `database`='', `sample.name`='', `unnormalized.expression`='', `normalized.expression`='')[-1,]
	}
}


#----------------
# get all tissues
#----------------


body.query <- 'https://www.proteomicsdb.org/proteomicsdb/logic/api/tissuelist.xsodata/CA_AVAILABLEBIOLOGICALSOURCES_API?$select=TISSUE_ID,TISSUE_NAME,TISSUE_GROUP_NAME,TISSUE_CATEGORY,SCOPE_ID,SCOPE_NAME,QUANTIFICATION_METHOD_ID,QUANTIFICATION_METHOD_NAME,MS_LEVEL&$format=json'


body.table <-
	PDBQuery(body.query) %>%
	set_names( c('metadata', 'tissue.id', 'tissue', 'tissue.group.name', 'tissue.category', 'scope.id', 'scope.name', 'quantification.method.id', 'quantification.method.name', 'ms.level'))


body.uri <- body.table[,1]['uri'] %>% unlist %>% unname


body <-
	body.table[,-1] %>%
	tbl_df %>%
	mutate(uri = body.uri) %>%
	filter(!tissue.category %in% 'cell line') %>%
	filter(!tissue %in% c('breast cancer cell', 'colorectal cancer cell', 'Unknown', 'renal cell carcinoma cell', 'milk', 'arachnoid cyst', 'amniocyte', 'chronic lymphocytic leukemia cell', 'osteosarcoma cell', 'ascites')) %>%
	select(tissue.id, tissue) %>%
	unique


#---------------------------------
# get all proteins for each tissue
#---------------------------------


pdb.stack <-
	map2(body$tissue.id, body$tissue, ~ { GetTissue(.x, .y) } ) %>%
	bind_rows


pdb.normal <-
	pdb.stack %>%
	select(uniprot, tissue, pdb.p = normalized.expression) %>%
	filter(!tissue %in% c('placenta')) %>%
	PDBFormatTissue %>%
	mutate(pdb.p = as.numeric(pdb.p)) %>%
	group_by(uniprot, tissue) %>%
	mutate(pdb.p = max(pdb.p, na.rm = TRUE)) %>%
	ungroup %>%
	select(uniprot, tissue, pdb.p) %>%
	unique


#----------------------------------------
# calculate normal distribution & cutoffs
#----------------------------------------

# PDB

pdb.curve <-
	pdb.normal %>%
	unique %>%
	select(exp = pdb.p) %>%
	PlotExpDist(binwidth = 0.02)

SavePlot(gg = pdb.curve$gg.dens, file = 'pdb_exp_dens.png', type = 'png')
SavePlot(gg = pdb.curve$gg.quant, file = 'pdb_exp_quant.png', type = 'png')

pdb.mu    <- pdb.curve$ml %>% filter(term == 'mu') %$% estimate
pdb.sigma <- pdb.curve$ml %>% filter(term == 'sigma') %$% estimate

# HPM

hpm.curve <-
	hpm.normal %>%
	select(tissue, exp = hpm.p.log) %>%
	filter(exp != 0) %>%
	unique %>%
	select(-tissue) %>%
	PlotExpDist(binwidth = 0.01)

SavePlot(gg = hpm.curve$gg.dens, file = 'hpm_exp_dens.png', type = 'png')
SavePlot(gg = hpm.curve$gg.quant, file = 'hpm_exp_quant.png', type = 'png')

hpm.mu    <- hpm.curve$ml %>% filter(term == 'mu') %$% estimate
hpm.sigma <- hpm.curve$ml %>% filter(term == 'sigma') %$% estimate




#------------#
#            #
#   LOCATE   #
#            #
#------------#




#------------------
# read & format xml
#------------------


locate <-
	xml2::read_xml('../data/LOCATE/LOCATE_human_v6_20081121.xml') %>%
	xml2::xml_children(.) %>%
	xml2::as_list(.)


uids <-
	locate %>%
	map(~ { .x %>% xml2::xml_attrs(.) %>% .['uid'] }) %>%
	unlist %>%
	unname



ParseNodes <- function(entry) {

	# progress
	p.num <- which(uids == entry %>% xml2::xml_attrs(.) %>% .['uid'])

	cat(str_c('\r', p.num, ' / ', length(uids)))

	flush.console()


	# parsing
	type <-
		entry %>%
		xml2::xml_children(.) %>%
		xml2::xml_name(.)

	protein <-
		entry %>%
		xml2::as_list(.) %>%
		.[which(type == 'protein')] %>%
		map(~ unlist(.) %>% rbind %>% as_data_frame) %>%
		bind_rows %>%
		select(organism, class, id.type = source1, id = source2) %>%
		mutate(link = '')

	if('scl_prediction' %in% type) {

		scl <-
			entry %>%
			xml2::as_list(.) %>%
			.[which(type == 'scl_prediction')] %>%
			map(~ unlist(.) %>% rbind %>% as_data_frame) %>%
			bind_rows %>%
			mutate(type = 'scl') %>%
			select(type, source = source1, location = source2)

	} else { scl <- NULL }

	if('externalannot' %in% type) {

		annot <-
			entry %>%
			xml2::as_list(.) %>%
			.[[which(type == 'externalannot')]] %>%
			map(~ unlist(.) %>% rbind %>% as_data_frame) %>%
			bind_rows %>%
			mutate(type = 'annot') %>%
			select(type, source = source1, location = locations.location)

	} else { annot <- NULL }

	if('xrefs' %in% type) {

		xref <-
			entry %>%
			xml2::as_list(.) %>%
			.[[which(type == 'xrefs')]] %>%
			map(~ unlist(.) %>% rbind %>% as_data_frame) %>%
			bind_rows %>%
			mutate(type = 'xref', link = '')

		if('source1' %in% names(xref) & 'source2' %in% names(xref)) {
			xref %<>% select(x.id.type = source1, x.id = source2, link)
		}

	} else { xref <- NULL }

	# join results & return
	full <-
		bind_rows(scl, annot) %>%
		mutate(link = '') %>%
		right_join(protein, by = 'link') %>%
		right_join(xref, by = 'link') %>%
		select(-link)

	full
}


locations <-
	locate %>% mclapply(., function(entry) ParseNodes(entry), mc.cores = 4)


extracellular <-
	c(	'basolateral plasma membrane',
		'Basolateral Plasma Membrane',
		# 'centrosome',
		# 'Centrosome',
		# 'centrosome,cytoplasm',
		# 'cytoplasm',
		# 'Cytoplasm',
		# 'cytoplasm_Golgi apparatus',
		# 'cytoplasm_mitochondrion',
		# 'cytoplasm_nucleus',
		# 'cytoplasm_peroxisome',
		'cytoplasm_plasma membrane',
		'cytoplasm,endosomes,lysosomes',
		# 'cytoplasm,nucleus',
		'cytoplasm,plasma membrane',
		'cytoplasmic membrane-bound vesicle',
		# 'cytoskeleton',
		# 'Cytoskeleton',
		'cytoskeleton_plasma membrane',
		'cytoskeleton,plasma membrane',
		'early endosomes',
		'Early Endosomes',
		# 'endoplasmic reticulum',
		# 'Endoplasmic Reticulum',
		# 'endoplasmic reticulum_Golgi apparatus',
		# 'endoplasmic reticulum_mitochondrion',
		# 'endoplasmic reticulum,Golgi apparatus',
		'endosomes',
		'Endosomes',
		'endosomes,Golgi apparatus,lysosomes',
		# 'ER-Golgi intermediate compartment',
		'extracellular region',
		'extracellular region_plasma membrane',
		# 'Golgi apparatus',
		# 'Golgi Apparatus',
		'Golgi apparatus,plasma membrane',
		# 'Golgi cis cisterna',
		# 'Golgi medial cisterna',
		# 'Golgi trans cisterna',
		'late endosomes',
		'Late Endosomes',
		# 'lysosome',
		# 'lysosomes',
		# 'Lysosomes',
		# 'melanosome',
		# 'mitochondrial inner membrane',
		# 'mitochondrial outer membrane',
		# 'mitochondrion',
		# 'mitochondrion_nucleus',
		# 'mitochondrion_peroxisome',
		# 'nuclear envelope',
		# 'Nuclear Envelope',
		# 'nuclear envelope,nucleolus',
		# 'nuclear speck',
		# 'nucleolus',
		# 'Nucleolus',
		# 'nucleolus,nucleus',
		# 'nucleus',
		# 'Nucleus',
		# 'nucleus,cytoplasm',
		# 'nucleus,nucleolus',
		# 'peroxisome',
		# 'peroxisomes',
		'plasma membrane',
		'Plasma membrane',
		'Plasma Membrane',
		'plasma membrane,cytoskeleton',
		'plasma membrane,endosomes',
		'secretory granule',
		'Secretory Granule',
		'secretory vesicles',
		'synaptic vesicles',
		'Synaptic Vesicles',
		'tight junction',
		'Tight junction'
	)

locate.location <-
	locations %>%
	bind_rows %>%
	tbl_df %>%
	select(location, id.a = id, id.b = x.id, id.type.a = id.type, id.type.b = x.id.type) %>%
	gather(id.col, id, id.a:id.b) %>%
	mutate(id.type = ifelse(id.col == 'id.a', id.type.a, id.type.b)) %>%
	select(id.type, id, location) %>%
	filter(!is.na(id)) %>%
	arrange(id.type, id, location) %>%
	mutate(membrane = ifelse(location %in% extracellular, TRUE, FALSE)) %>%
	filter(membrane == TRUE) %>%
	select(id.type, id, location) %>%
	unique %>%
	mutate(id.type = ifelse(id.type == 'RefSeq Protein' & (grepl('^XP_', id) | grepl('^AP_', id)), 'RefSeq Protein X', id.type)) %>%
	filter(id.type %in% c('UniProtKB-SwissProt', 'Ensembl-Peptide Human', 'RefSeq Protein', 'RefSeq Protein X', 'Entrez Protein', 'UniProt/SPTrEMBL', 'UniGene', 'Entrez Gene', 'PDB', 'HGNC')) %>%
	mutate(id.type = ifelse(grepl('_HUMAN', id), 'UniProtKB', id.type)) %>%
	mutate(database = 'LOCATE')




#--------------------#
#                    #
#   GENE ID LOOKUP   #
#                    #
#--------------------#




#-------
# lookup
#-------

mart <- biomaRt::useEnsembl(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl', GRCh = 37, version = 88)

# attributes <- biomaRt::listAttributes(mart) %>% tbl_df

# filters <- biomaRt::listFilters(mart) %>% tbl_df

GetMartID <- function(query, id.type) {

	receipt <-
		biomaRt::getBM(
			attributes = unique(c('ensembl_gene_id', id.type)),
			filters = id.type,
			values = query,
			mart = mart
		) %>%
		tbl_df

	if(ncol(receipt) == 1) {
		receipt %<>% mutate(query = ensembl_gene_id)
	}

	receipt %>%
	set_names(c('ensembl.gene', 'query')) %>%
	mutate(ensembl.gene = as.character(ensembl.gene), query = as.character(query)) %>%
	right_join(data_frame(query = as.character(query)), by = 'query') %>%
	select(query, ensembl.gene)
}


GetMyGeneID <- function(query, id.type) {

	receipt <-
		suppressWarnings(try(
			mygene::queryMany(
				qterms    = query,
				fields    = 'ensembl.gene',
				scopes    = c('ensembl.gene', id.type),
				species   = 'human',
				returnall = TRUE )
		)) %$%
		response

		if('ensembl.gene' %in% names(receipt)) {
			receipt <- receipt[,c('query', 'ensembl.gene')]

			map2(as.list(receipt$query), receipt$ensembl.gene, ~ { data_frame(query = .x, ensembl.gene = unlist(.y)) }) %>%
			bind_rows

		} else {
			receipt <- receipt[,c('query', 'ensembl')]

			receipt.query <- as.list(receipt$query)
			receipt.ensembl <- map(receipt$ensembl, ~{
				ensembl <- unlist(.x)
				if(is.null(ensembl)) {
					ensembl <- NA
				}
				ensembl
			})

			map2(receipt.query, receipt.ensembl, ~ {
				data_frame(query = .x, ensembl.gene = unlist(.y))
			}) %>%
			bind_rows
		}
}


GetUniID <- function(query, id.type) {

	split(query, ceiling(seq_along(query) / 500)) %>%
	map( ~{

		sub.query <- .

		request <- suppressMessages(content(GET("http://www.uniprot.org/uploadlists/",
			query = list(
				query = sub.query %>% str_c(collapse = ' '),
				format = 'tab',
				from = id.type,
				to = 'ENSEMBL_ID'
			),
			add_headers('User-Agent', email)
		), 'parsed', encoding = 'UTF-8'))

		connect <- textConnection(request)
		output <-
			read.delim(connect, sep = '\t', stringsAsFactors = FALSE) %>%
			as_data_frame
		close(connect)
		output

	}) %>%
	bind_rows %>%
	set_names(c('query', 'ensembl.gene'))
}


GetUniUCSC <- function(query) {

	query.list <- split(query, ceiling(seq_along(query)/10))

	receipt <-
		query.list %>%
		map( ~ {
			ucsc.con <- dbConnect(MySQL(), user = 'genome', host = 'genome-mysql.cse.ucsc.edu', dbname = 'uniProt')

			ucsc.query <- dbSendQuery(ucsc.con, str_c("SELECT * FROM displayId WHERE val in('", str_c(.x, collapse = "','"), "')"))

			ucsc.fetch <-
				fetch(ucsc.query) %>%
				tbl_df %>%
				set_names(c('uniprotkb.entry', 'uniprotkb.name')) %>%
				full_join(data_frame(uniprotkb.name = .x), by = 'uniprotkb.name')

			suppressWarnings(dbDisconnect(ucsc.con))

			ucsc.fetch
		}) %>%
		bind_rows
}

legacy.uniprot.map <-
	read_excel('../data/UniprotKB_legacy_IDs.xlsx') %>%
	select(uniprotkb.entry, uniprotkb.entry.name, gene.names) %>%
	separate(uniprotkb.entry.name, into = LETTERS[1:6], sep = '\\,', fill = 'right') %>%
	gather(column, uniprotkb.entry.name, A:F) %>%
	select(-column) %>%
	mutate(uniprotkb.entry.name = gsub('_HUMAN', '', uniprotkb.entry.name)) %>%
	rowwise %>%
	mutate(gene.names = str_c(c(uniprotkb.entry.name, gene.names), collapse = ' ')) %>%
	separate(gene.names, into = LETTERS[1:10], sep = '\\ ', fill = 'right') %>%
	gather(column, gene.name, A:J) %>%
	select(uniprotkb.entry, gene.name) %>%
	unique %>%
	arrange(gene.name)


GetID <- function(query, id.type) {

	print(id.type)

	if(id.type == 'Ensembl-Gene Human') {
		response <-
			bind_rows(
				GetMartID(query, 'ensembl_gene_id'),
				GetMyGeneID(query, 'ensemblgene')
			)
	}

	if(id.type == 'Ensembl-Peptide Human') {
		response <-
			bind_rows(
				GetMartID(query, 'ensembl_peptide_id'),
				GetMyGeneID(query, 'ensemblprotein')
			)
	}

	if(id.type == 'Entrez Gene') {
		response <- GetMartID(query, 'entrezgene')
	}

	if(id.type == 'Entrez Protein') {
		response <- GetMartID(query, 'protein_id')
	}

	if(id.type == 'HGNC') {
		response <- GetMartID(query, 'hgnc_id')
	}

	if(id.type == 'PDB') {
		response <-
			bind_rows(
				GetMartID(query, 'pdb'),
				GetMyGeneID(query, 'pdb')
			)
	}

	if(id.type == 'RefSeq Protein') {
		response <-
			bind_rows(
				GetMartID(query, 'refseq_peptide'),
				GetMyGeneID(query, 'refseq')
			)
	}

	if(id.type == 'RefSeq Protein X') {
		response <-
			bind_rows(
				GetMartID(query, 'refseq_peptide_predicted'),
				GetMartID(query, 'refseq_peptide')
			)
	}

	if(id.type == 'symbol') {
		response <-
			bind_rows(
				GetMartID(query, 'hgnc_symbol'),
				GetMyGeneID(query, 'symbol'),
				GetMyGeneID(query, 'alias')
			)
	}

	if(id.type == 'UniGene') {
		response <-
			bind_rows(
				GetMartID(query, 'unigene'),
				GetMyGeneID(query, 'unigene')
			)
	}

	if(id.type == 'UniProt/SPTrEMBL') {
		response <- GetMartID(query, 'uniprotsptrembl')
	}

	if(id.type == 'UniProtKB-SwissProt') {
		response <-
			bind_rows(
				GetMartID(query, 'uniprotswissprot'),
				GetMyGeneID(query, 'uniprot')
			)
	}

	if(id.type == 'UniProtKB') {

		receipt <- GetUniUCSC(query)

		uni.ucsc <-
			receipt %>%
			filter(is.na(uniprotkb.entry)) %>%
			mutate(query.trunc = gsub('_HUMAN', '', uniprotkb.name)) %>%
			select(-uniprotkb.entry) %>%
			left_join(legacy.uniprot.map, by = c('query.trunc' = 'gene.name')) %>%
			filter(!is.na(uniprotkb.entry)) %>%
			select(-query.trunc)

		response <-
			bind_rows(
				GetMartID(uni.ucsc$uniprotkb.entry, 'uniprotswissprot') %>% left_join(uni.ucsc, by = c('query' = 'uniprotkb.entry')) %>% select(-query) %>% rename(query = uniprotkb.name),
				GetMyGeneID(uni.ucsc$uniprotkb.entry, 'uniprot') %>% left_join(uni.ucsc, by = c('query' = 'uniprotkb.entry')) %>% select(-query) %>% rename(query = uniprotkb.name),
				GetUniID(query, 'ACC+ID')
			)
	}
	response
}


GetHugoID <- function(query) {

	query %<>% list.filter(!is.na(.))

	hugo.mart <- 
		biomaRt::getBM(
			attributes = c('hgnc_symbol', 'ensembl_gene_id'),
			filters = 'ensembl_gene_id',
			values = query,
			mart = mart
		) %>%
		tbl_df %>%
		set_names(c('hgnc', 'ensembl.gene')) %>%
		mutate(hgnc = as.character(hgnc), ensembl.gene = as.character(ensembl.gene)) %>%
		right_join(data_frame(ensembl.gene = as.character(query)), by = 'ensembl.gene') %>%
		mutate(hgnc = ifelse(hgnc == '', NA, hgnc))

	receipt <-
		suppressWarnings(try(
			mygene::queryMany(
				qterms    = query,
				fields    = 'symbol',
				scopes    = c('symbol', 'ensembl.gene'),
				species   = 'human',
				returnall = TRUE )
		)) %$%
		response

	if('symbol' %in% names(receipt)) {
		receipt <- receipt[,c('query', 'symbol')]

		hugo.mygene <-
			map2(as.list(receipt$query), receipt$symbol, ~ { data_frame(query = .x, symbol = unlist(.y)) }) %>%
			bind_rows %>%
			set_names(c('ensembl.gene', 'hgnc'))

	} else {
		receipt <- receipt[,c('query', 'symbol')]

		receipt.query <- as.list(receipt$query)
		receipt.symbol <- map(receipt$symbol, ~{
			symbol <- unlist(.x)
			if(is.null(symbol)) {
				symbol <- NA
			}
			symbol
		})

		hugo.mygene <-
			map2(receipt.query, receipt.symbol, ~ {
				data_frame(query = .x, symbol = unlist(.y))
			}) %>%
			bind_rows %>%
			set_names(c('ensembl.gene', 'hgnc'))
	}

	bind_rows(hugo.mart, hugo.mygene) %>%
	unique %>%
	set_names(c('hgnc', 'ensembl.gene')) %>%
	right_join(data_frame(ensembl.gene = as.character(query)), by = 'ensembl.gene') %>%
	mutate(hgnc = ifelse(hgnc == '', NA, hgnc)) %>%
	group_by(ensembl.gene) %>%
	filter(n() == 1 | !all(is.na(hgnc)) & !is.na(hgnc)) %>%
	ungroup %>%
	unique
}


GetEnsemblID <- function(query) {

	query %<>% list.filter(!is.na(.))

	bind_rows(
		GetMartID(query, 'hgnc_symbol'),
		GetMyGeneID(query, 'symbol'),
		GetMyGeneID(query, 'alias')
	) %>%
	unique %>%
	set_names(c('hgnc', 'ensembl.gene')) %>%
	right_join(data_frame(hgnc = as.character(query)), by = 'hgnc') %>%
	mutate(ensembl.gene = ifelse(ensembl.gene == '', NA, ensembl.gene)) %>%
	group_by(hgnc) %>%
	filter(n() == 1 | !all(is.na(ensembl.gene)) & !is.na(ensembl.gene)) %>%
	ungroup %>%
	unique
}




#----------------#
#                #
#   CELL LINES   #
#                #
#----------------#




#---------------
# MSK AML single
#---------------


msk.aml.single <-
	read.xlsx('../data/cell_lines/2.11.15_MSK_AML_single.xlsx', sheet = 'MS_surface') %>%
	tbl_df %>%
	select(
		id            = Accession_Number,
		msk.09aml.0   = `09AML`,
		msk.kasum1.0  = KASUM1,
		msk.molm13.0  = MOLM13,
		msk.monomac.0 = MONOMAC,
		msk.tf.0      = TF,
		msk.thp1.0    = THP1 ) %>%
	filter(!grepl('-DECOY', id)) %>%
	mutate(id = gsub(' \\(\\+([0-9])\\)', '', id)) %>%
	unique


#-------------------
# MSK AML triplicate
#-------------------


msk.aml.triple <-
	read.xlsx('../data/cell_lines/5.3.15_MSK_AML_triplicate.xlsx', sheet = 'MS_surface') %>%
	tbl_df %>%
	select(
		id            = Accession_Number,
		msk.kasumi.1  = S04,
		msk.kasumi.2  = S05,
		msk.kasumi.3  = S06,
		msk.thp1.1    = S07,
		msk.thp1.2    = S08,
		msk.thp1.3    = S09,
		msk.monomac.1 = S10,
		msk.monomac.2 = S11,
		msk.monomac.3 = S12,
		msk.molm13.1  = S13,
		msk.molm13.2  = S14,
		msk.molm13.3  = S15) %>%
	filter(!grepl('-DECOY', id)) %>%
	mutate(id = gsub(' \\(\\+([0-9])\\)', '', id)) %>%
	unique


#-------------------------
# DNMT3a mutant triplicate
#-------------------------


msk.dnmt3a.triple <-
	read.xlsx('../data/cell_lines/5.22.16_MSK_DNMT3a_triplicate.xlsx', sheet = 'MS_surface') %>%
	tbl_df %>%
	select(
		id            = Accession_Number,
		msk.dnmt3a.1  = S16,
		msk.dnmt3a.2  = S17,
		msk.dnmt3a.3  = S18) %>%
	filter(!grepl('-DECOY', id)) %>%
	mutate(id = gsub(' \\(\\+([0-9])\\)', '', id)) %>%
	unique


#----------------------------
# J Proteomics AML triplicate
#----------------------------


jpro.aml.triple <-
	read.xlsx('../data/cell_lines/30.1.14_JPro_AML_triplicate.xlsx', sheet = 'Membrane', startRow = 2) %>%
	tbl_df %>%
	select( id              = Gene,
			jpro.thp1.1 = THP1_1,
			jpro.thp1.2 = THP1_2,
			jpro.thp1.3 = THP1_3)




#----------------------#
#                      #
#   build dictionary   #
#                      #
#----------------------#



id.query <-
	bind_rows(
		locate.location    %>% select(query = id, id.type),
		hpa.location       %>% select(query = ensembl.gene) %>% mutate(id.type = 'Ensembl-Gene Human'),
		hpa.normal.cut     %>% select(query = ensembl.gene) %>% mutate(id.type = 'Ensembl-Gene Human'),
		hpa.rna.normal.cut %>% select(query = ensembl.gene) %>% mutate(id.type = 'Ensembl-Gene Human'),
		pdb.normal         %>% select(query = uniprot)      %>% mutate(id.type = 'UniProtKB-SwissProt'),
		hpm.normal         %>% select(query = hgnc)         %>% mutate(id.type = 'symbol'),
		msk.aml.single     %>% select(query = id)           %>% mutate(id.type = 'UniProtKB'),
		msk.aml.triple     %>% select(query = id)           %>% mutate(id.type = 'UniProtKB'),
		msk.dnmt3a.triple  %>% select(query = id)           %>% mutate(id.type = 'UniProtKB'),
		jpro.aml.triple    %>% select(query = id)           %>% mutate(id.type = 'UniProtKB')
	) %>%
	mutate(query = ifelse(id.type == 'UniProtKB' & !grepl('_HUMAN', query), str_c(query, '_HUMAN'), query)) %>%
	unique %>%
	arrange(id.type, query)


ensembl.id.dictionary <-
	id.query %>%
	split(.$id.type) %>%
	map(~ {
		query <- .x %$% query %>% unlist
		id.type <- .x %$% id.type %>% unique
		GetID(query, id.type)
	}) %>%
	bind_rows %>%
	tbl_df %>%
	filter(!is.na(ensembl.gene)) %>%
	unique %>%
	right_join(data_frame(query = unique(id.query$query))) %>%
	bind_rows(data_frame(  # jpro missing lookups
		query = c('MT-ND4', 'NCR3LG1', 'MCEMP1', 'EMC1', 'PIEZO1', 'TREM1', 'CLEC5A'),
		ensembl.gene = c('ENSG00000198886', 'ENSG00000188211', 'ENSG00000183019', 'ENSG00000127463', 'ENSG00000103335', 'ENSG00000124731', 'ENSG00000258227')))


#----------#
#          #
# assembly #
#          #
#----------#



#---------------------------------------
# id lookups & duplicate tissue maximums
#---------------------------------------


hpm.normal.cut <-
	hpm.normal %>%
	left_join(ensembl.id.dictionary, by = c('hgnc' = 'query')) %>%
	filter(!is.na(ensembl.gene)) %>%
	group_by(ensembl.gene, tissue) %>%
	mutate(hpm.p.log.max = max(hpm.p.log)) %>%
	ungroup %>%
	select(ensembl.gene, tissue, hpm.p.log.max) %>%
	unique %>%
	mutate(hpm.p.exp =
		ifelse(hpm.p.log.max < (hpm.mu - hpm.sigma),                                          0,
		ifelse(hpm.p.log.max <  hpm.mu               & hpm.p.log.max >= (hpm.mu - hpm.sigma), 1,
		ifelse(hpm.p.log.max < (hpm.mu + hpm.sigma)  & hpm.p.log.max >=  hpm.mu,              2,
		ifelse(                                        hpm.p.log.max >= (hpm.mu + hpm.sigma), 3,
	NA))))) %>%
	select(ensembl.gene, tissue, hpm.p.exp)


pdb.normal.cut <-
	pdb.normal %>%
	left_join(ensembl.id.dictionary, by = c('uniprot' = 'query')) %>%
	filter(!is.na(ensembl.gene)) %>%
	group_by(ensembl.gene, tissue) %>%
	mutate(pdb.p.max = max(pdb.p)) %>%
	ungroup %>%
	select(ensembl.gene, tissue, pdb.p.max) %>%
	unique %>%
	mutate(pdb.p.exp =
		ifelse(pdb.p.max < (pdb.mu - pdb.sigma),                                     0,
		ifelse(pdb.p.max <  pdb.mu              & pdb.p.max >= (pdb.mu - pdb.sigma), 1,
		ifelse(pdb.p.max < (pdb.mu + pdb.sigma) & pdb.p.max >=  pdb.mu,              2,
		ifelse(                                   pdb.p.max >= (pdb.mu + pdb.sigma), 3,
	NA))))) %>%
	select(ensembl.gene, tissue, pdb.p.exp)

pdb.normal.cut.fill <-
	pdb.normal.cut %>%
	full_join(expand.grid(
		ensembl.gene = unique(pdb.normal.cut$ensembl.gene),
		tissue       = unique(pdb.normal.cut$tissue),
		pdb.p.exp    = 0,
		stringsAsFactors = FALSE)) %>%
	group_by(ensembl.gene, tissue) %>%
	mutate(pdb.p.exp = max(pdb.p.exp)) %>%
	ungroup %>%
	unique

# In general, the expression values provided by ProteomicsDB are (sum total) normalized (iBAQ or Top3)
# expression values of the proteins. They are a rough approximation of protein copy numbers in the corresponding
# tissue, cell line or body fluid and range from roughly 2-8 in log scale. Judging by eye I would estimate that
# values < 4.5 can be considered low and >6.5 high.


#--------------------------------------
# combine LOCATE & HPA subcellular data
#--------------------------------------


membrane.location <-
	locate.location %>%
	left_join(ensembl.id.dictionary, by = c('id' = 'query')) %>%
	select(ensembl.gene, location, database) %>%
	filter(!is.na(ensembl.gene)) %>%
	bind_rows(hpa.location) %>%
	arrange(ensembl.gene) %>%
	unique

membrane.location.clean <-
	membrane.location %>%
	select(ensembl.gene) %>%
	unique %>%
	mutate(membrane = TRUE)


#---------------------
# MSK AML lines single
#---------------------


msk.aml.single.clean <-
	msk.aml.single %>%
	left_join(ensembl.id.dictionary, c('id' = 'query')) %>%
	select(ensembl.gene, msk.09aml.0, msk.kasum1.0, msk.molm13.0, msk.monomac.0, msk.tf.0, msk.thp1.0) %>%
	group_by(ensembl.gene) %>%
	mutate(
		msk.09aml.0   = mean(msk.09aml.0),
		msk.kasum1.0  = mean(msk.kasum1.0),
		msk.molm13.0  = mean(msk.molm13.0),
		msk.monomac.0 = mean(msk.monomac.0),
		msk.tf.0      = mean(msk.tf.0),
		msk.thp1.0    = mean(msk.thp1.0)
	) %>%
	ungroup %>%
	unique


#-------------------------
# MSK AML lines triplicate
#-------------------------


msk.aml.triple.clean <-
	msk.aml.triple %>%
	left_join(ensembl.id.dictionary, c('id' = 'query')) %>%
	select(ensembl.gene, msk.kasumi.1, msk.kasumi.2, msk.kasumi.3, msk.thp1.1, msk.thp1.2, msk.thp1.3,
		   msk.monomac.1, msk.monomac.2, msk.monomac.3, msk.molm13.1, msk.molm13.2, msk.molm13.3) %>%
	group_by(ensembl.gene) %>%
	mutate(
		msk.1.kasumi  = mean(msk.kasumi.1),
		msk.2.kasumi  = mean(msk.kasumi.2),
		msk.3.kasumi  = mean(msk.kasumi.3),
		msk.1.thp1    = mean(msk.thp1.1),
		msk.2.thp1    = mean(msk.thp1.2),
		msk.3.thp1    = mean(msk.thp1.3),
		msk.1.monomac = mean(msk.monomac.1),
		msk.2.monomac = mean(msk.monomac.2),
		msk.3.monomac = mean(msk.monomac.3),
		msk.1.molm13  = mean(msk.molm13.1),
		msk.2.molm13  = mean(msk.molm13.2),
		msk.3.molm13  = mean(msk.molm13.3)) %>%
	ungroup %>%
	unique


#-------------------
# DNMT3a mutant line
#-------------------


msk.dnmt3a.triple.clean <-
	msk.dnmt3a.triple %>%
	left_join(ensembl.id.dictionary, by = c('id' = 'query')) %>%
	select(ensembl.gene, msk.dnmt3a.1, msk.dnmt3a.2, msk.dnmt3a.3) %>%
	group_by(ensembl.gene) %>%
	mutate(msk.dnmt3a.1 = mean(msk.dnmt3a.1)) %>%
	mutate(msk.dnmt3a.2 = mean(msk.dnmt3a.2)) %>%
	mutate(msk.dnmt3a.3 = mean(msk.dnmt3a.3)) %>%
	ungroup %>%
	unique


#----------------------------------
# J Proteomics AML lines triplicate
#----------------------------------


jpro.aml.triple.clean <-
	jpro.aml.triple %>%
	left_join(ensembl.id.dictionary, by = c('id' = 'query')) %>%
	select(ensembl.gene, jpro.thp1.1, jpro.thp1.2, jpro.thp1.3) %>%
	group_by(ensembl.gene) %>%
	mutate(jpro.thp1.1 = mean(jpro.thp1.1)) %>%
	mutate(jpro.thp1.2 = mean(jpro.thp1.2)) %>%
	mutate(jpro.thp1.3 = mean(jpro.thp1.3)) %>%
	ungroup %>%
	unique


#----------------------------------
# join datasets & normalize sources
#----------------------------------


master <-
	hpa.normal.cut %>%
	full_join(hpm.normal.cut,          by = c('ensembl.gene', 'tissue')) %>%
	full_join(pdb.normal.cut.fill,     by = c('ensembl.gene', 'tissue')) %>%
	full_join(hpa.rna.normal.cut,      by = c('ensembl.gene', 'tissue')) %>%
	full_join(membrane.location.clean, by = 'ensembl.gene') %>%
	left_join(msk.aml.single.clean,    by = 'ensembl.gene') %>%
	left_join(msk.aml.triple.clean,    by = 'ensembl.gene') %>%
	left_join(msk.dnmt3a.triple.clean,        by = 'ensembl.gene') %>%
	left_join(jpro.aml.triple.clean,   by = 'ensembl.gene') %>%
	filter(!is.na(tissue)) %>%
	select(
		ensembl.gene, tissue, membrane,
		hpa = hpa.p.exp,
		hpm = hpm.p.exp,
		pdb = pdb.p.exp,
		rna = hpa.r.exp,
		msk.09aml.0, msk.kasum1.0, msk.molm13.0, msk.monomac.0, msk.tf.0, msk.thp1.0,
		msk.kasumi.1, msk.thp1.1, msk.monomac.1, msk.molm13.1,
		msk.kasumi.2, msk.thp1.2, msk.monomac.2, msk.molm13.2,
		msk.kasumi.3, msk.thp1.3, msk.monomac.3, msk.molm13.3,
		msk.dnmt3a.1, msk.dnmt3a.2, msk.dnmt3a.3,
		jpro.thp1.1, jpro.thp1.2, jpro.thp1.3) %>%
	inner_join(GetHugoID(unique(master$ensembl.gene)), by = 'ensembl.gene') %>%
	select(ensembl.gene, hgnc, tissue, everything()) %>%
	arrange(hgnc) %>%
	group_by(ensembl.gene, tissue) %>%
	slice(1) %>%
	ungroup %>%
	group_by(ensembl.gene, tissue) %>%
	mutate(tissue.mean       = mean(c(hpa, hpm, pdb),     na.rm = TRUE)) %>%
	mutate(tissue.max        = max( c(hpa, hpm, pdb), -1, na.rm = TRUE)) %>%
	ungroup %>%
	mutate(tissue.mean       = ifelse(is.nan(tissue.mean), NA, tissue.mean)) %>%
	mutate(tissue.max        = ifelse(tissue.max == -1,    NA, tissue.max)) %>%
	group_by(ensembl.gene) %>%
	mutate(db.num = as.numeric(!all(is.na(hpa))) + as.numeric(!all(is.na(hpm))) + as.numeric(!all(is.na(pdb)))) %>%
	mutate(prot.mean         = mean(tissue.mean, na.rm = TRUE)) %>%
	mutate(rna.mean          = mean(rna,         na.rm = TRUE)) %>%
	mutate(prot.fill         = sum(!is.na(tissue.max)) / n()) %>%
	mutate(prot.rna.cor      = cor(tissue.max, rna, use = 'pairwise.complete.obs', method = 'spearman')) %>%
	mutate(prot.rna.var      = var(tissue.max, rna, na.rm = TRUE)) %>%
	mutate(prot.rna.fill     = sum(!is.na(tissue.max) & !is.na(rna)) / n()) %>%
	mutate(mean.cor          = mean(prot.rna.cor, na.rm = TRUE)) %>%
	ungroup %>%
	mutate(mean.cor          = ifelse(is.nan(mean.cor),  NA, mean.cor)) %>%
	mutate(prot.mean         = ifelse(is.nan(prot.mean), NA, prot.mean))


msk.aml.single.surface <-
	c(msk.aml.single.clean$ensembl.gene) %>%
	na.omit %>%
	unique %>%
	sort

msk.aml.triple.surface <-
	c(msk.aml.triple.clean$ensembl.gene) %>%
	na.omit %>%
	unique %>%
	sort

msk.dnmt3a.triple.surface <-
	msk.dnmt3a.triple.clean$ensembl.gene %>%
	unique %>%
	sort

jpro.aml.triple.surface <-
	jpro.aml.triple.clean$ensembl.gene %>%
	unique %>%
	sort

save(list = c('msk.aml.single.surface', 'msk.aml.triple.surface', 'msk.dnmt3a.triple.surface', 'jpro.aml.triple.surface', 'master'), file = '_clean.RData')


#--------
# objects
#--------

objects <- c(
'HPAFormatTissue',
'HPMFormatTissue',
'PDBFormatTissue',
'PDBQuery',
'GetTissue',
'ParseNodes',
'GetMartID',
'GetMyGeneID',
'GetUniID',
'GetUniUCSC',
'GetID',
'GetHugoID',
'GetEnsemblID',
'hpa.normal.cut',
'hpa.rna.normal.cut',
'hpa.extracellular',
'hpa.location',
'hpm.normal',
'body.query',
'body.table',
'body.uri',
'body',
'pdb.stack',
'pdb.normal',
'pdb.curve',
'pdb.mu',
'pdb.sigma',
'hpm.curve',
'hpm.mu',
'hpm.sigma',
'locate',
'uids',
'locations',
'extracellular',
'locate.location',
'mart',
'legacy.uniprot.map',
'msk.aml.single',
'msk.aml.triple',
'msk.dnmt3a.triple',
'jpro.aml.triple',
'id.query',
'ensembl.id.dictionary',
'hpm.normal.cut',
'pdb.normal.cut',
'pdb.normal.cut.fill',
'membrane.location',
'membrane.location.clean',
'msk.aml.single.clean',
'msk.aml.triple.clean',
'msk.dnmt3a.triple',
'jpro.aml.triple.clean',
'master',
'msk.aml.single.surface',
'msk.aml.triple.surface',
'msk.dnmt3a.triple.surface',
'jpro.aml.triple.surface')
