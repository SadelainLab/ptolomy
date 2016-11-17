



#-------------#
#             #
#  retreival  #
#             #
#-------------#




#--------
# imports
#--------


loadNamespace('mygene')
loadNamespace('biomaRt')
load('protein.data.RData')

source('secrets.R')


#-------------------------
# gene id lookup functions
#-------------------------


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
# AML candidate subset
#---------------------


msk.jpro.add <-
	c(msk.lines.clean$ensembl.gene, jpro.lines.clean$ensembl.gene, aml.additions) %>%
	unique %>%
	sort

pass.tissues <- c('bone', 'blood')


#--------------------
# candidate selection
#--------------------

id.dict <-
	protein.data$ensembl.gene %>%
	unique %>%
	GetHugoID %>%
	arrange(hgnc) %>%
	group_by(hgnc) %>%
	slice(1) %>%
	ungroup

step.0.0 <-
	protein.data %>%
	filter(!is.na(tissue)) %>%
	select(-mean.abundance) %>%
	unique %>%
	rowwise %>%
	mutate(tissue.max = max(c(hpa, hpm, pdb), -1, na.rm = TRUE)) %>%
	mutate(tissue.mean = mean(c(hpa, hpm, pdb), na.rm = TRUE)) %>%
	ungroup %>%
	mutate(tissue.max = ifelse(tissue.max == -1, NA, tissue.max)) %>%
	mutate(tissue.mean = ifelse(is.nan(tissue.mean), NA, tissue.mean)) %>%
	left_join(id.dict, by = 'ensembl.gene') %>%
	select(hgnc, ensembl.gene, db.num, gross.mean.abundance, tissue, hpa, hpm, pdb, tissue.max, everything())

step.0.1 <-
	step.0.0 %>%
	group_by(hgnc, ensembl.gene) %>%
	unique %>%
	mutate(fill = sum(!is.na(tissue.max)) / 42) %>%
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

step.0 <-
	step.0.0 %>%
	inner_join(step.0.1, by = c('hgnc', 'ensembl.gene')) %>%
	select(hgnc, ensembl.gene, db.num, gross.mean.abundance, membrane, fill, everything())

step.1 <-
	step.0 %>%
	filter(db.num >= 2) %>%
	filter(membrane == TRUE) %>%
	filter(ensembl.gene %in% msk.jpro.add) %>%
	arrange(gross.mean.abundance, hgnc) %>%
	select(hgnc, ensembl.gene, tissue, db.num, fill, hpa, hpm, pdb, tissue.max, tissue.mean, gross.mean.abundance, rna, everything(), -membrane)

step.2.0 <-
	step.1 %>%
	filter(gross.mean.abundance < 1) %>%
	group_by(hgnc) %>%
	mutate(mean.test = tissue.mean < 2.5 | tissue %in% pass.tissues) %>%
	mutate(max.test = tissue.max < 3 | tissue %in% pass.tissues) %>%
	filter(!any(mean.test == FALSE)) %>%
	#filter(!any(max.test == FALSE)) %>%
	mutate(aml.tot =
		any(!is.na(msk.0.09aml) & msk.0.09aml != 0) +
		any(!is.na(msk.0.kasum1) & msk.0.kasum1 != 0) +
		any(!is.na(msk.0.molm13) & msk.0.molm13 != 0) +
		any(!is.na(msk.0.monomac) & msk.0.monomac != 0) +
		any(!is.na(msk.0.tf) & msk.0.tf != 0) +
		any(!is.na(msk.0.thp1) & msk.0.thp1 != 0) +
		any(!is.na(jpro.thp1) & jpro.thp1 != 0)
	) %>%
	ungroup %>%
	mutate(tissue.mean.round =
		ifelse(tissue.mean < 0.5, 0,
		ifelse(tissue.mean < 1.5, 1,
		ifelse(tissue.mean < 2.5, 2,
		ifelse(tissue.mean <= 3, 3,
			NA))))) %>%
	arrange(gross.mean.abundance, ensembl.gene) #%>%
#	filter(ensembl.gene %in% picks)


# step.2.0 %>% select(hgnc, tissue, tissue.mean, tissue.max, mean.test, max.test) %>% filter(mean.test == FALSE | max.test == FALSE)
# step.2.0 %>% select(hgnc, tissue, tissue.mean, tissue.max, mean.test, max.test) %>% filter(mean.test == FALSE | max.test == FALSE)
# step.2.0 %>% select(hgnc, tissue, level = tissue.mean) %>% PlotTissue


step.2.1 <-
	step.2.0 %>%
	arrange(gross.mean.abundance) %>%
	select(hgnc) %>%
	unique %>%
	mutate(group.n = row_number())

step.2.2 <-
	step.2.0 %>%
	full_join(expand.grid(hgnc = unique(step.2.0$hgnc), tissue = unique(step.2.0$tissue), stringsAsFactors = FALSE), by = c('hgnc', 'tissue')) %>%
	left_join(step.2.1, by = 'hgnc') %>%
	arrange(hgnc, tissue)

step.2.2 %>% filter(ensembl.gene %in% picks) %$% hgnc %>% unique

#-----------

step.3 <-
	combn(1:length(unique(step.2.2$hgnc)), 2) %>%
	t %>%
	as_data_frame %>%
	set_names(c('a', 'b')) %>%
	mutate(pair.num = row_number()) %>%
	rowwise %>%
	do({
			x <- .
			bind_cols(
				step.2.2 %>% filter(group.n == x$a) %>% set_names(str_c(colnames(step.2.2), '_a')),
				step.2.2 %>% filter(group.n == x$b) %>% set_names(str_c(colnames(step.2.2), '_b'))
			) %>%
			mutate(pair.num = x$pair.num) %>%
			select(pair.num, everything())
	}) %>%
	ungroup %>%
	filter(!is.na(ensembl.gene_a) & !is.na(ensembl.gene_b))

vital <- c(
'adipose tissue',
'adrenal',
'bladder',
'brain',
'bronchus',
'eye',
'gut',
'heart',
'kidney',
'laryngopharynx',
'liver',
'lung',
'nasopharynx',
'oropharynx',
'pancreas',
'rectum',
'skeletal muscle',
'skin',
'smooth muscle',
'soft tissue',
'spinal cord',
'stomach')

non.vital <- c(
'appendix',
'bone',
'blood',
'breast',
'cerumen',
'cervix',
'epididymis',
'fallopian tube',
'gallbladder',
'lymph node',
'ovary',
'parathyroid',
'prostate',
'seminal',
'spleen',
'synovial fluid',
'testis',
'thyroid',
'tonsil',
'uterus',
'vagina')

step.4.0 <-
	step.3 %>%
	rowwise %>%
	mutate(level.1 =
		!(tissue.max_a >= 2 & tissue.max_b >= 2)) %>%
	mutate(level.2 =
		!((tissue_a %in% non.vital & tissue.max_a >= 2 & tissue.max_b >= 2) |
		(tissue_a %in% vital & tissue.max_a >= 1 & tissue.max_b >= 1))) %>%
	ungroup %>%
	select(-hpa_a, -hpa_b, -hpm_a, -hpm_b, -pdb_a, -pdb_b, -rna_a, -rna_b) %>%
	unique %>%
	select(tissue_a, tissue_b, tissue.max_a, tissue.max_b, everything()) %>%
	gather(tissue_x, tissue, tissue_a:tissue_b) %>%
	mutate(tissue_x = ifelse(tissue_x == 'tissue_a', 'a',
		ifelse(tissue_x == 'tissue_b', 'b', NA))) %>%
	gather(tissue.max_x, tissue.max, tissue.max_a:tissue.max_b) %>%
	filter((tissue_x == 'a' & tissue.max_x == 'tissue.max_a') | (tissue_x == 'b' & tissue.max_x == 'tissue.max_b')) %>%
	mutate(tissue = str_c(tissue, tissue_x, sep = '_')) %>%
	select(-tissue_x, -tissue.max_x) %>%
	unique %>%
	spread(tissue, tissue.max)

step.4.level.1 <-
	step.4.0 %>%
	group_by(pair.num) %>%
	filter(all(level.1 == TRUE)) %>%
	ungroup

step.4.level.2 <-
	step.4.0 %>%
	group_by(pair.num) %>%
	filter(all(level.2 == TRUE)) %>%
	ungroup



tt <-
step.4.0 %>% mutate(hgnc = str_c(hgnc_a, hgnc_b, sep = ' | ')) %>%
select(gene = hgnc, `adipose tissue_a`:vagina_b) %>%
gather(tissue, level, `adipose tissue_a`:vagina_b) %>%
filter(!is.na(level))

tt.1 <-
step.4.level.1 %>% mutate(hgnc = str_c(hgnc_a, hgnc_b, sep = ' | ')) %>%
select(gene = hgnc, `adipose tissue_a`:vagina_b) %>%
gather(tissue, level, `adipose tissue_a`:vagina_b) %>%
filter(!is.na(level))

tt.2 <-
step.4.level.2 %>% mutate(hgnc = str_c(hgnc_a, hgnc_b, sep = ' | ')) %>%
select(gene = hgnc, `adipose tissue_a`:vagina_b) %>%
gather(tissue, level, `adipose tissue_a`:vagina_b) %>%
filter(!is.na(level))


PlotTissue <- function(events, pdf = FALSE, file.name = 'plot_tissue.pdf', width = 20, height = 10) {

	palette <- c(
		'zero'   = '#3E5496',
		'low'    = '#E8E8E8',
		'medium' = '#C39696',
		'high'   = '#933E56'
		)

	if(all(na.omit(events$level) %% 1 == 0)) {	# check if integer, if so plot discrete
		events %<>%
			mutate(level =
				ifelse(level == 0, 'zero',
				ifelse(level == 1, 'low',
				ifelse(level == 2, 'medium',
				ifelse(level == 3, 'high',
				NA))))) %>%
			mutate(level = factor(level, levels = unique(level)))

		m.gg <-
			ggplot(events, aes(tissue, gene)) +
			geom_tile(aes(fill = level, drop = FALSE), colour = 'grey') +
			scale_fill_manual(
				breaks           = names(palette),
				values           = palette,
				na.value         = 'grey',
				drop             = FALSE,
				guide            = guide_legend(reverse = TRUE))
	} else {
		m.gg <-
			ggplot(events, aes(tissue, gene)) + 
			geom_tile(aes(fill = level), colour = 'grey') +
			scale_fill_gradientn(
				colours = palette,
				na.value = 'transparent',
				breaks = 0:3,
				labels = names(palette),
				limits = c(0, 3))
	}

	mt.gg <-
		m.gg +
		theme(
			text             = element_text(size = 10),
			axis.text.x      = element_text(angle = 90, hjust = 1),
			panel.background = element_rect(fill = 'grey'),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border     = element_rect(colour = 'black', fill = NA, size = 1))

	if(pdf == TRUE) {
		pdf(file.name, width, height)
			plot(mt.gg)
		dev.off()
	} else {
		plot(mt.gg)
	}
}

#step.2.2 %>% select(gene = hgnc, tissue, level = tissue.max) %>% PlotTissue


PlotTissue(tt)

PlotTissue(tt.1, file.name = 'combi_level_1.pdf', pdf = TRUE, width = 20, height = 2)

PlotTissue(tt.2, file.name = 'combi_level_2.pdf', pdf = TRUE, width = 20, height = 2)

