#--------
# imports
#--------


pacman::p_load(jsonlite, httr, parallel)

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
		.$tissue %in% c('esophagus') ~ 'laryngopharynx',
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
		.$tissue %in% c('Adult.Esophagus') ~ 'laryngopharynx',
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
		.$tissue %in% c('esophagus') ~ 'laryngopharynx',
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




hpa.normal <-
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
	HPAFormatTissue


hpa.rna.normal <-
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
	HPAFormatTissue


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
	HPMFormatTissue




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
			"',CALCULATION_METHOD=0,SWISSPROT_ONLY=1,NO_ISOFORM=1)/Results?$select=ENTRY_NAME,UNIQUE_IDENTIFIER,DATABASE,SAMPLE_NAME,UNNORMALIZED_EXPRESSION,NORMALIZED_EXPRESSION&$format=json"
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
	select(uniprot, tissue, pdb.p)


#rm(pdb.stack, body.table, body.uri, body.query)





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




#--------
# imports
#--------


loadNamespace('mygene')
loadNamespace('biomaRt')


#-------
# lookup
#-------

mart <- biomaRt::useEnsembl(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl', GRCh = 37)

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
		response <- GetMartID(query, 'uniprot_sptrembl')
	}

	if(id.type == 'UniProtKB-SwissProt') {
		response <-
			bind_rows(
				GetMartID(query, 'uniprot_swissprot'),
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
				GetMartID(uni.ucsc$uniprotkb.entry, 'uniprot_swissprot') %>% left_join(uni.ucsc, by = c('query' = 'uniprotkb.entry')) %>% select(-query) %>% rename(query = uniprotkb.name),
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


#---------------------
# MSK AML lines single
#---------------------


msk.lines.single <-
	read.xlsx('../data/cell_lines/2.11.15_MSK_AML_lines.xlsx', sheet = 'surface_proteins') %>%
	tbl_df %>%
	select(
		id            = Accession.Number,
		msk.0.09aml   = `09AML`,
		msk.0.kasum1  = KASUM1,
		msk.0.molm13  = MOLM13,
		msk.0.monomac = MONOMAC,
		msk.0.tf      = TF,
		msk.0.thp1    = THP1 ) %>%
	filter(!grepl('-DECOY', id)) %>%
	mutate(id = gsub(' \\(\\+([0-9])\\)', '', id)) %>%
	unique


#---------------------
# MSK AML lines triple
#---------------------


msk.lines.triple <-
	read.xlsx('../data/cell_lines/5.3.15_MSK_AML_lines_triplicate.xlsx', sheet = 'surface_proteins') %>%
	tbl_df %>%
	select(
		id = Accession.Number,
		msk.1.kasumi  = S04,
		msk.2.kasumi  = S05,
		msk.3.kasumi  = S06,
		msk.1.thp1    = S07,
		msk.2.thp1    = S08,
		msk.3.thp1    = S09,
		msk.1.monomac = S10,
		msk.2.monomac = S11,
		msk.3.monomac = S12,
		msk.1.molm13  = S13,
		msk.2.molm13  = S14,
		msk.3.molm13  = S15) %>%
	filter(!grepl('-DECOY', id)) %>%
	mutate(id = gsub(' \\(\\+([0-9])\\)', '', id)) %>%
	unique


#-----------------------
# J Proteomics AML lines
#-----------------------


jpro.lines <-
	read.xlsx('../data/cell_lines/30.1.14_JProteomics_AML_lines.xlsx', sheet = 1, startRow = 2) %>%
	tbl_df %>%
	select( id = Gene,
			thp1.1 = THP1_1, thp1.2 = THP1_2, thp1.3 = THP1_3 ) %>%
	rowwise %>%
	mutate(jpro.thp1 = mean(thp1.1, thp1.2, thp1.3)) %>%
	ungroup %>%
	select(id, jpro.thp1)



#----------------------#
#                      #
#   build dictionary   #
#                      #
#----------------------#




id.query <-
	bind_rows(
		locate.location  %>% select(query = id, id.type),
		hpa.location     %>% select(query = ensembl.gene) %>% mutate(id.type = 'Ensembl-Gene Human'),
		hpa.normal       %>% select(query = ensembl.gene) %>% mutate(id.type = 'Ensembl-Gene Human'),
		hpa.rna.normal   %>% select(query = ensembl.gene) %>% mutate(id.type = 'Ensembl-Gene Human'),
		pdb.normal       %>% select(query = uniprot)      %>% mutate(id.type = 'UniProtKB-SwissProt'),
		hpm.normal       %>% select(query = hgnc)         %>% mutate(id.type = 'symbol'),
		msk.lines.single %>% select(query = id)           %>% mutate(id.type = 'UniProtKB'),
		msk.lines.triple %>% select(query = id)           %>% mutate(id.type = 'UniProtKB'),
		jpro.lines       %>% select(query = id)           %>% mutate(id.type = 'UniProtKB')
	) %>%
	mutate(query = ifelse(id.type == 'UniProtKB' & !grepl('_HUMAN', query), str_c(query, '_HUMAN'), query)) %>%
	unique %>%
	arrange(id.type, query)


library(RMySQL)

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
	right_join(data_frame(query = unique(id.query$query)))




#----------#
#          #
# assembly #
#          #
#----------#




#---------------------------------------
# id lookups & duplicate tissue maximums
#---------------------------------------


hpa.normal.cut <-
	hpa.normal %>%
	group_by(ensembl.gene, tissue) %>%
	mutate(hpa.p.max = max(hpa.p)) %>%
	ungroup %>%
	unique %>%
	mutate(hpa.p.max.norm =
		ifelse(hpa.p.max < 0.5,                    0,
		ifelse(hpa.p.max < 1.5 & hpa.p.max >= 0.5, 1,
		ifelse(hpa.p.max < 2.5 & hpa.p.max >= 1.5, 2,
		ifelse(                  hpa.p.max >= 2.5, 3,
	NA))))) %>%
	select(ensembl.gene, tissue, hpa.p.max.norm)


hpm.normal.clean <-
	hpm.normal %>%
	left_join(ensembl.id.dictionary, by = c('hgnc' = 'query')) %>%
	select(ensembl.gene, tissue, hpm.p) %>%
	filter(!is.na(ensembl.gene)) %>%
	group_by(ensembl.gene, tissue) %>%
	mutate(hpm.p.max      = max(hpm.p)) %>%
	select(-hpm.p) %>%
	ungroup %>%
	unique %>%
	mutate(hpm.p.max      = log10(hpm.p.max)) %>%
	mutate(hpm.p.max      = ifelse(hpm.p.max == -Inf, 0, hpm.p.max))


hpm.normal.cut <-
	hpm.normal.clean %>%
	mutate(hpm.p.max.norm =
		ifelse(hpm.p.max < hpm.p.cut.low,                                0,
		ifelse(hpm.p.max < hpm.p.mean     & hpm.p.max >= hpm.p.cut.low,  1,
		ifelse(hpm.p.max < hpm.p.cut.high & hpm.p.max >= hpm.p.mean,     2,
		ifelse(                             hpm.p.max >= hpm.p.cut.high, 3,
	NA))))) %>%
	select(ensembl.gene, tissue, hpm.p.max.norm)


pdb.normal.clean <-
	pdb.normal %>%
	left_join(ensembl.id.dictionary, by = c('uniprot' = 'query')) %>%
	select(ensembl.gene, tissue, pdb.p) %>%
	filter(!is.na(ensembl.gene)) %>%
	group_by(ensembl.gene, tissue) %>%
	mutate(pdb.p.max = max(pdb.p)) %>%
	ungroup %>%
	unique

pdb.normal.cut <-
	pdb.normal.clean %>%
	mutate(pdb.p.max.sd   = sd(pdb.p.max)) %>%
	mutate(pdb.p.max.norm =
		ifelse(pdb.p.max < pdb.p.max.mean - pdb.p.max.sd,                                              0,
		ifelse(pdb.p.max < pdb.p.max.mean &                pdb.p.max >= pdb.p.max.mean - pdb.p.max.sd, 1,
		ifelse(pdb.p.max < pdb.p.max.mean + pdb.p.max.sd & pdb.p.max >= pdb.p.max.mean,                2,
		ifelse(                                            pdb.p.max >= pdb.p.max.mean + pdb.p.max.sd, 3,
	NA))))) %>%
	select(ensembl.gene, tissue, pdb.p.max.norm)




loadNamespace('bbmle')
library(broom)

SavePlot <- function(gg, file.name) {
	pdf(file = file.name)
		print(gg)
	dev.off()
}


PlotExp <- function(exp.table, binwidth = 0.02) {

	exp.table %<>%
		filter(!is.na(exp)) %>%
		arrange(exp) %>%
		mutate(num = row_number() - 1)

	NormFit <- function(mu, sigma) { -sum(dnorm(x, mu, sigma, log = TRUE)) }

	ml <-
		bbmle::mle2(
			NormFit,
			start     = list(mu = mean(exp.table$exp), sigma = sd(exp.table$exp)),
			data      = list(x = exp.table$exp),
			method    = 'BFGS') %>%
		tidy

	ml.mu <- ml %>% filter(term == 'mu') %$% estimate
	ml.sigma <- ml %>% filter(term == 'sigma') %$% estimate

	gg.norm <-
		ggplot(exp.table, aes(x = exp)) +
		geom_histogram(aes(y = ..density..), alpha = 0.5, binwidth = binwidth) +
		scale_x_continuous(breaks = 0:ceiling(max(exp.table$exp))) +
		scale_y_continuous(breaks = (0:20) / 20) +
		stat_function(
			fun   = dnorm,
			color = 'blue',
			args  = list(mean = ml.mu, sd = ml.sigma)) +
		geom_vline(xintercept = ml.mu, linetype = 2, color = 'black') +
		geom_vline(xintercept = ml.mu + ml.sigma, linetype = 2, color = 'purple') +
		geom_vline(xintercept = ml.mu - ml.sigma, linetype = 2, color = 'purple')

	gg.sort <-
		ggplot(exp.table) +
		geom_point(aes(x = num, y = exp), size = 0.5, alpha = 0.25) +
		scale_x_continuous(breaks = seq(0, nrow(exp.table), nrow(exp.table) / 10), labels = (0:10) / 10) +
		scale_y_continuous(breaks = pretty_breaks(8), limits = c(NA, ceiling(max(exp.table$exp)))) +
		scale_alpha_manual(0.2) +
		geom_hline(aes(yintercept = ml.mu), color = 'black', linetype = 2) +
		geom_text(
			data = data.frame(x = 0, y = ml.mu),
			aes(x, y),
			label = str_c('µ = ', round(ml.mu, 2)),
			hjust = 0,
			vjust = -0.5) +
		geom_text(
			data = data.frame(x = 0, y = ml.mu),
			aes(x, y),
			label = str_c('µ + σ = ', round(ml.mu + ml.sigma, 2)),
			hjust = 0,
			vjust = -6.4) +
		geom_text(
			data = data.frame(x = 0, y = ml.mu),
			aes(x, y),
			label = str_c('µ - σ = ', round(ml.mu - ml.sigma, 2)),
			hjust = 0,
			vjust = 5.4) +
		# stat_function(
		# 	fun   = qnorm,
		# 	color = 'blue',
		# 	args  = list(mean = ml.mu, sd = ml.sigma)) +
		geom_hline(aes(yintercept = ml.mu + ml.sigma), color = 'purple', linetype = 2) +
		geom_hline(aes(yintercept = ml.mu - ml.sigma), color = 'purple', linetype = 2)

	list(ml = ml, gg.norm = gg.norm, gg.sort = gg.sort)
}



pdb.curve <-
	pdb.normal.clean %>%
	mutate(exp = pdb.p.max) %>%
	select(exp) %>%
	unique %>%
	PlotExp(binwidth = 0.02)

hpm.curve <-
	hpm.normal.clean %>%
	mutate(exp = hpm.p.max) %>%
	select(exp) %>%
	unique %>%
	PlotExp(binwidth = 0.02)


	pdb.mu <- pdb.curve$ml %>% filter(term == 'mu') %$% estimate
	pdb.sigma <- pdb.curve$ml %>% filter(term == 'sigma') %$% estimate

	hpm.mu <- hpm.curve$ml %>% filter(term == 'mu') %$% estimate
	hpm.sigma <- hpm.curve$ml %>% filter(term == 'sigma') %$% estimate


# In general, the expression values provided by ProteomicsDB are (sum total) normalized (iBAQ or Top3)
# expression values of the proteins. They are a rough approximation of protein copy numbers in the corresponding
# tissue, cell line or body fluid and range from roughly 2-8 in log scale. Judging by eye I would estimate that
# values < 4.5 can be considered low and >6.5 high.


hpa.rna.normal.cut <-
	hpa.rna.normal %>%
	group_by(ensembl.gene, tissue) %>%
	mutate(hpa.r.max = max(hpa.r)) %>%
	ungroup %>%
	unique %>%
	mutate(hpa.r.max.norm =
		ifelse(hpa.r.max < 0.5,                    0,
		ifelse(hpa.r.max < 1.5 & hpa.r.max >= 0.5, 1,
		ifelse(hpa.r.max < 2.5 & hpa.r.max >= 1.5, 2,
		ifelse(                  hpa.r.max >= 2.5, 3,
	NA))))) %>%
	select(ensembl.gene, tissue, hpa.r.max.norm)


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

msk.lines.single.clean <-
	msk.lines.single %>%
	left_join(ensembl.id.dictionary, c('id' = 'query')) %>%
	select(ensembl.gene, msk.0.09aml, msk.0.kasum1, msk.0.molm13, msk.0.monomac, msk.0.tf, msk.0.thp1) %>%
	group_by(ensembl.gene) %>%
	mutate(
		msk.0.09aml   = mean(msk.0.09aml),
		msk.0.kasum1  = mean(msk.0.kasum1),
		msk.0.molm13  = mean(msk.0.molm13),
		msk.0.monomac = mean(msk.0.monomac),
		msk.0.tf      = mean(msk.0.tf),
		msk.0.thp1    = mean(msk.0.thp1)
	) %>%
	ungroup %>%
	unique


#---------------------
# MSK AML lines triple
#---------------------


msk.lines.triple.clean <-
	msk.lines.triple %>%
	left_join(ensembl.id.dictionary, c('id' = 'query')) %>%
	select(ensembl.gene, msk.1.kasumi, msk.2.kasumi, msk.3.kasumi, msk.1.thp1, msk.2.thp1, msk.3.thp1,
		   msk.1.monomac, msk.2.monomac, msk.3.monomac, msk.1.molm13, msk.2.molm13, msk.3.molm13) %>%
	group_by(ensembl.gene) %>%
	mutate(
		msk.1.kasumi  = mean(msk.1.kasumi),
		msk.2.kasumi  = mean(msk.2.kasumi),
		msk.3.kasumi  = mean(msk.3.kasumi),
		msk.1.thp1    = mean(msk.1.thp1),
		msk.2.thp1    = mean(msk.2.thp1),
		msk.3.thp1    = mean(msk.3.thp1),
		msk.1.monomac = mean(msk.1.monomac),
		msk.2.monomac = mean(msk.2.monomac),
		msk.3.monomac = mean(msk.3.monomac),
		msk.1.molm13  = mean(msk.1.molm13),
		msk.2.molm13  = mean(msk.2.molm13),
		msk.3.molm13  = mean(msk.3.molm13)) %>%
	ungroup %>%
	unique


#-----------------------
# J Proteomics AML lines
#-----------------------


jpro.lines.clean <-
	jpro.lines %>%
	left_join(jpro.lines.dictionary, by = c('id' = 'query')) %>%
	select(ensembl.gene, jpro.thp1) %>%
	group_by(ensembl.gene) %>%
	mutate(jpro.thp1 = mean(jpro.thp1)) %>%
	ungroup %>%
	unique


#----------------------------------
# join datasets & normalize sources
#----------------------------------


protein.data <-
	hpa.normal.cut %>%
	full_join(hpm.normal.cut,          by = c('ensembl.gene', 'tissue')) %>%
	full_join(pdb.normal.cut,          by = c('ensembl.gene', 'tissue')) %>%
	full_join(hpa.rna.normal.cut,      by = c('ensembl.gene', 'tissue')) %>%
	full_join(membrane.location.clean, by = 'ensembl.gene') %>%
	left_join(msk.lines.single.clean,  by = 'ensembl.gene') %>%
	left_join(msk.lines.triple.clean,  by = 'ensembl.gene') %>%
	left_join(jpro.lines.clean,        by = 'ensembl.gene') %>%
	filter(!is.na(tissue)) %>%
	select(
		ensembl.gene, tissue, membrane,
		hpa = hpa.p.max.norm,
		hpm = hpm.p.max.norm,
		pdb = pdb.p.max.norm,
		rna = hpa.r.max.norm,
		msk.0.09aml, msk.0.kasum1, msk.0.molm13, msk.0.monomac, msk.0.tf, msk.0.thp1,
		msk.1.kasumi, msk.1.thp1, msk.1.monomac, msk.1.molm13,
		msk.2.kasumi, msk.2.thp1, msk.2.monomac, msk.2.molm13,
		msk.3.kasumi, msk.3.thp1, msk.3.monomac, msk.3.molm13,
		jpro.thp1) %>%
		inner_join(GetHugoID(unique(protein.data$ensembl.gene)), by = 'ensembl.gene') %>%
		select(ensembl.gene, hgnc, tissue, everything()) %>%
		arrange(hgnc) %>%
		group_by(ensembl.gene, tissue) %>%
		slice(1) %>%
		ungroup



step.0 <-
	protein.data %>%
	group_by(ensembl.gene, tissue) %>%
	mutate(tissue.mean       = mean(c(hpa, hpm, pdb),     na.rm = TRUE)) %>%
	mutate(tissue.max        = max( c(hpa, hpm, pdb), -1, na.rm = TRUE)) %>%
	ungroup %>%
	mutate(tissue.mean       = ifelse(is.nan(tissue.mean), NA, tissue.mean)) %>%
	mutate(tissue.max        = ifelse(tissue.max == -1,    NA, tissue.max)) %>%
	group_by(ensembl.gene) %>%
	mutate(db.num            = (sum(0, hpa, na.rm = TRUE) > 0) +
							   (sum(0, hpm, na.rm = TRUE) > 0) +
							   (sum(0, pdb, na.rm = TRUE) > 0)) %>%
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

#	select(hgnc, ensembl.gene, tissue, db.num, mean.abundance, gross.mean, membrane, fill, hpa, hpm, pdb, tissue.max, tissue.mean, everything(), -prot.num.tissues)




protein.data.fill <-
	protein.data.slice %>%
	group_by(hgnc, ensembl.gene) %>%
	unique %>%
	ungroup %>%
	select(hgnc, ensembl.gene, gross.mean.abundance, fill) %>%
	unique %>%
	arrange(desc(fill)) %>%
	group_by(hgnc) %>%
	slice(1) %>%
	ungroup %>%
	group_by(ensembl.gene) %>%
	slice(1) %>%
	ungroup %>%
	select(hgnc, ensembl.gene, fill)




	protein.data %>%
	group_by(ensembl.gene, tissue) %>%

	ungroup %>%

	group_by(ensembl.gene) %>%

	ungroup %>%
	


	unique %>%
	rowwise %>%
	mutate(tissue.mean = mean(c(hpa, hpm, pdb), na.rm = TRUE)) %>%
	ungroup %>%

	left_join(hugo.id.dictionary, by = 'ensembl.gene') %>%
	mutate(num.tissues = tissue %>% unique %>% length) %>%
	select(hgnc, ensembl.gene, tissue, num.tissues, db.num, mean.abundance, gross.mean.abundance, membrane, hpa, hpm, pdb, tissue.max, tissue.mean, everything())











msk.surface <-
	c(msk.lines.triple.clean$ensembl.gene, msk.lines.single.clean$ensembl.gene) %>%
	na.omit %>%
	#msk.lines.triple.clean$ensembl.gene %>%
	#msk.lines.single.clean$ensembl.gene %>%
	unique %>%
	sort

jpro.surface <-
	jpro.lines.clean$ensembl.gene %>%
	sort


save(list = c('msk.surface', 'jpro.surface', 'step.0'), file = 'clean.RData')


