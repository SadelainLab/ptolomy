



#-------------#
#             #
#  retreival  #
#             #
#-------------#




#--------
# imports
#--------


source('secrets.R')


#----------
# variables
#----------


pass.tissue <- c('bone', 'blood')


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


#----------
# functions
#----------


palette <- c(
	'zero'   = '#202c99',
	'low'    = '#fbeaea',
	'medium' = '#e77e7e',
	'high'   = '#c12525'
	)


PlotTissue <- function(events, facet = FALSE, pdf = FALSE, file.name = 'plot_tissue.pdf', width = 20, height = 10) {

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

	if(facet == TRUE) {
		mt.gg <-
			m.gg +
			theme(
				text             = element_text(size = 10),
				axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1),
				panel.background = element_rect(fill = 'grey'),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				panel.border     = element_rect(colour = 'black', fill = NA, size = 1),
				strip.background = element_blank(),
				strip.text.x     = element_blank()) +
			facet_wrap(~ split, ncol = 1, scales = 'free_y')
	} else {
		mt.gg <-
			m.gg +
			theme(
				text             = element_text(size = 10),
				axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1),
				panel.background = element_rect(fill = 'grey'),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				panel.border     = element_rect(colour = 'black', fill = NA, size = 1))
	}

	if(pdf == TRUE) {
		pdf(file.name, width, height)
			plot(mt.gg)
		dev.off()
	} else {
		getOption('device')()
		plot(mt.gg)
	}
}


#---------------------
# AML candidate subset
#---------------------


msk.jpro.add <-
	c(msk.lines.clean$ensembl.gene, jpro.lines.clean$ensembl.gene, aml.additions) %>%
	unique %>%
	sort


#--------------------
# candidate selection
#--------------------

step.1 <-
	step.0 %>%
	filter(db.num >= 2) %>%
	filter(membrane == TRUE) %>%
	filter(ensembl.gene %in% msk.jpro.add) %>%
	arrange(gross.mean.abundance, hgnc) %>%
	select(hgnc, ensembl.gene, tissue, db.num, fill, hpa, hpm, pdb, tissue.max, tissue.mean, gross.mean.abundance, rna, everything(), -membrane)

# candidate comparison  -- rank in order of quality
step.1 %>%
filter(hgnc %in% compare.protein) %>%
arrange(gross.mean.abundance) %>%
mutate(hgnc = factor(hgnc, levels = unique(hgnc))) %>%
select(gene = hgnc, tissue, level = tissue.max) %>%
mutate(tissue = factor(tissue, levels = c(vital, non.vital))) %>%
PlotTissue(pdf = TRUE, width = 10, height = 2.5, file.name = 'compare_mean.pdf')

step.1 %>%
filter(hgnc %in% compare.protein) %>%
arrange(gross.mean.abundance) %>%
mutate(hgnc = factor(hgnc, levels = unique(hgnc))) %>%
select(gene = hgnc, tissue, level = tissue.mean) %>%
mutate(tissue = factor(tissue, levels = c(vital, non.vital))) %>%
PlotTissue(pdf = TRUE, width = 10, height = 2.5, file.name = 'compare_max.pdf')



step.2.0 <-
	step.1 %>%
	filter(gross.mean.abundance < 1 | hgnc %in% pass.protein) %>%
	group_by(hgnc) %>%
	mutate(tissue.mean.round =
		ifelse(tissue.mean <  0.5, 0,
		ifelse(tissue.mean <  1.5, 1,
		ifelse(tissue.mean <  2.5, 2,
		ifelse(tissue.mean <= 3,   3,
			NA))))) %>%
	mutate(max.test = tissue.max < 3 | tissue %in% pass.tissues | hgnc %in% pass.protein)


step.2.1 <-
	step.2.0 %>%
	filter(!any(na.omit(max.test) == FALSE) | hgnc %in% pass.protein) %>%
	mutate(aml.tot =
		any(!is.na(msk.0.09aml)   & msk.0.09aml   != 0) +
		any(!is.na(msk.0.kasum1)  & msk.0.kasum1  != 0) +
		any(!is.na(msk.0.molm13)  & msk.0.molm13  != 0) +
		any(!is.na(msk.0.monomac) & msk.0.monomac != 0) +
		any(!is.na(msk.0.tf)      & msk.0.tf      != 0) +
		any(!is.na(msk.0.thp1)    & msk.0.thp1    != 0) +
		any(!is.na(jpro.thp1)     & jpro.thp1     != 0)
	) %>%
	ungroup %>%
	arrange(gross.mean.abundance, ensembl.gene) %>%
	filter(!hgnc %in% exclude.protein)

step.2.1 %>%
select(-tissue.mean, -hpa, -pdb, -hpm, -rna, -tissue.mean.round, -max.test) %>%
spread(tissue, tissue.max) %>%
unique %>%
write_tsv('step_2.tsv')


step.2.2 <-
	step.2.1 %>%
	arrange(gross.mean.abundance) %>%
	select(hgnc) %>%
	unique %>%
	mutate(group.n = row_number())

step.2.3 <-
	step.2.1 %>%
	full_join(expand.grid(hgnc = unique(step.2.1$hgnc), tissue = unique(step.2.1$tissue), stringsAsFactors = FALSE), by = c('hgnc', 'tissue')) %>%
	left_join(step.2.2, by = 'hgnc') %>%
	arrange(hgnc, tissue)

step.2.4 <-
	step.2.3 %>%
	filter(hgnc %in% names(picks) | hgnc %in% pass.protein)

step.2.5 <-
	step.2.3 %>%
	filter(hgnc %in% names(picks))

step.2.3 %>%
filter(!hgnc %in% pass.protein) %>%
arrange(gross.mean.abundance) %>%
mutate(hgnc = factor(hgnc, levels = unique(hgnc))) %>%
select(gene = hgnc, tissue, level = tissue.max) %>%
filter(!tissue %in% c('cerumen', 'synovial fluid')) %>%
mutate(tissue = factor(tissue, levels = c(vital, non.vital))) %>%
PlotTissue(pdf = TRUE, width = 10, height = 17.5, file.name = 'step_2_max.pdf')

step.2.3 %>%
filter(!hgnc %in% pass.protein) %>%
arrange(gross.mean.abundance) %>%
mutate(hgnc = factor(hgnc, levels = unique(hgnc))) %>%
select(gene = hgnc, tissue, level = tissue.mean) %>%
filter(!tissue %in% c('cerumen', 'synovial fluid')) %>%
mutate(tissue = factor(tissue, levels = c(vital, non.vital))) %>%
PlotTissue(pdf = TRUE, width = 10, height = 17.5, file.name = 'step_2_mean.pdf')

step.2.5 %>%
filter(ensembl.gene %in% picks) %>%
filter(!hgnc %in% pass.protein) %>%
arrange(gross.mean.abundance) %>%
mutate(hgnc = factor(hgnc, levels = unique(hgnc))) %>%
select(gene = hgnc, tissue, level = tissue.max) %>%
filter(!tissue %in% c('cerumen', 'synovial fluid')) %>%
mutate(tissue = factor(tissue, levels = c(vital, non.vital))) %>%
PlotTissue(pdf = TRUE, width = 10, height = 6.2, file.name = 'step_2_picks_max.pdf')

step.2.5 %>%
filter(ensembl.gene %in% picks) %>%
filter(!hgnc %in% pass.protein) %>%
arrange(gross.mean.abundance) %>%
mutate(hgnc = factor(hgnc, levels = unique(hgnc))) %>%
select(gene = hgnc, tissue, level = tissue.mean) %>%
filter(!tissue %in% c('cerumen', 'synovial fluid')) %>%
mutate(tissue = factor(tissue, levels = c(vital, non.vital))) %>%
PlotTissue(pdf = TRUE, width = 10, height = 6.2, file.name = 'step_2_picks_mean.pdf')


#-----------

step.3 <-
	combn(1:length(unique(step.2.3$hgnc)), 2) %>%
	t %>%
	as_data_frame %>%
	set_names(c('a', 'b')) %>%
	mutate(pair.num = row_number()) %>%
	rowwise %>%
	do({
			x <- .
			bind_cols(
				step.2.3 %>% filter(group.n == x$a) %>% set_names(str_c(colnames(step.2.3), '_a')),
				step.2.3 %>% filter(group.n == x$b) %>% set_names(str_c(colnames(step.2.3), '_b'))
			) %>%
			mutate(pair.num = x$pair.num) %>%
			select(pair.num, everything())
	}) %>%
	ungroup %>%
	filter(!is.na(ensembl.gene_a) & !is.na(ensembl.gene_b))


#-------


step.4 <-
	step.3 %>%
	rowwise %>%
	mutate(max.level.1 =
		tissue.max_a < 2 | is.na(tissue.max_a) | tissue.max_b < 2 | is.na(tissue.max_b)) %>%
	mutate(max.level.2 =
		(tissue_a %in% non.vital & (tissue.max_a < 2 | is.na(tissue.max_a) | tissue.max_b < 2 | is.na(tissue.max_b))) |
		(tissue_a %in% vital     & (tissue.max_a < 1 | is.na(tissue.max_a) | tissue.max_b < 1 | is.na(tissue.max_b))) ) %>%
	ungroup %>%
	select(tissue_a, tissue_b, tissue.max_a, tissue.max_b, everything(), -hpa_a, -hpa_b, -hpm_a, -hpm_b, -pdb_a, -pdb_b, -rna_a, -rna_b) %>%
	unique %>%
	gather(tissue_x, tissue, tissue_a:tissue_b) %>%
	gather(tissue.max_x, tissue.max, tissue.max_a:tissue.max_b) %>%
	filter(
		(tissue_x == 'tissue_a' & tissue.max_x == 'tissue.max_a') |
		(tissue_x == 'tissue_b' & tissue.max_x == 'tissue.max_b') ) %>%
	mutate(tissue_x =
		ifelse(tissue_x == 'tissue_a', 'a',
		ifelse(tissue_x == 'tissue_b', 'b',
		NA))) %>%
	mutate(gross.mean.abundance = ifelse(tissue_x == 'a', gross.mean.abundance_a, gross.mean.abundance_b)) %>%
	mutate(hgnc = ifelse(tissue_x == 'a', hgnc_a, hgnc_b)) %>%
	select(pair.num, pair = tissue_x, hgnc, tissue, tissue.max, gross.mean.abundance, max.level.1, max.level.2) %>%
	unique


step.4.level.1 <-
	step.4 %>%
	group_by(pair.num) %>%
	filter(all(max.level.1 == TRUE)) %>%
	mutate(order = mean(gross.mean.abundance)) %>%
	ungroup %>%
	arrange(order) %>%
	filter(!tissue %in% c('cerumen', 'synovial fluid')) %>%
	mutate(pair.num = factor(pair.num, levels = unique(pair.num))) %>%
	mutate(tissue = factor(tissue, levels = c(vital, non.vital)))

step.4.level.2 <-
	step.4 %>%
	group_by(pair.num) %>%
	filter(all(max.level.2 == TRUE)) %>%
	mutate(order = mean(gross.mean.abundance)) %>%
	ungroup %>%
	arrange(order) %>%
	filter(!tissue %in% c('cerumen', 'synovial fluid')) %>%
	mutate(pair.num = factor(pair.num, levels = unique(pair.num))) %>%
	mutate(tissue = factor(tissue, levels = c(vital, non.vital)))



step.4 %>%
filter(pair.num == 571) %>%
filter(!tissue %in% c('cerumen', 'synovial fluid')) %>%
mutate(pair.num = factor(pair.num, levels = unique(pair.num))) %>%
mutate(tissue = factor(tissue, levels = c(vital, non.vital))) %>%
select(split = pair.num, gene = hgnc, tissue, level = tissue.max) %>%
filter(!is.na(level)) %>%
PlotTissue(facet = TRUE, pdf = TRUE, width = 10, height = 1.7, file.name = 'combi_LTB4R_CD70.pdf')



step.4 %>%
filter(pair.num == 1645) %>%
filter(!tissue %in% c('cerumen', 'synovial fluid')) %>%
mutate(pair.num = factor(pair.num, levels = unique(pair.num))) %>%
mutate(tissue = factor(tissue, levels = c(vital, non.vital))) %>%
select(split = pair.num, gene = hgnc, tissue, level = tissue.max) %>%
filter(!is.na(level)) %>%
PlotTissue(facet = TRUE, pdf = TRUE, width = 10, height = 1.7, file.name = 'combi_CD33_CD70.pdf')

tt.0 <-
	step.4.level.2 %>%
	group_by(pair.num) %>%
	mutate(pass.test = hgnc %in% pass.protein) %>%
	filter(any(pass.test == TRUE)) %>%
	arrange(pair.num, pair) %>%
	select(split = pair.num, gene = hgnc, tissue, level = tissue.max) %>%
	filter(!is.na(level))


tt.1 <-
	step.4.level.1 %>%
	group_by(pair.num) %>%
	mutate(pass.test = !hgnc %in% pass.protein) %>%
	filter(any(pass.test == TRUE)) %>%
	select(split = pair.num, gene = hgnc, tissue, level = tissue.max) %>%
	filter(!is.na(level))

tt.2 <-
	step.4.level.2 %>%
	group_by(pair.num) %>%
	mutate(pass.test = !hgnc %in% pass.protein) %>%
	filter(any(pass.test == TRUE)) %>%
	select(split = pair.num, gene = hgnc, tissue, level = tissue.max) %>%
	filter(!is.na(level))



PlotTissue(tt.1, facet = TRUE, pdf = TRUE, width = 10, height = 400, file.name = 'combi_level_1.pdf')


PlotTissue(tt.2, facet = TRUE, pdf = TRUE, width = 10, height = 80, file.name = 'combi_level_2.pdf')



tt.2 %>% group_by(split) %>% mutate(test = any(gene == 'CD70')) %>% filter(any(test == TRUE)) %>% arrange(split) %>%
PlotTissue(facet = TRUE, pdf = TRUE, width = 10, height = 3, file.name = 'combi_level_2_CD70.pdf')

tt.2 %>% group_by(split) %>% mutate(test = any(gene == 'LTB4R')) %>% filter(any(test == TRUE)) %>% arrange(split) %>%
PlotTissue(facet = TRUE, pdf = TRUE, width = 10, height = 6.5, file.name = 'combi_level_2_LTB4R.pdf')

tt.2 %>% group_by(split) %>% mutate(test = any(gene == 'ADGRE2')) %>% filter(any(test == TRUE)) %>% arrange(split) %>%
PlotTissue(facet = TRUE, pdf = TRUE, width = 10, height = 1.7, file.name = 'combi_level_2_ADGRE2.pdf')



combi.sub <-
	step.4.level.1 %>%
	group_by(pair.num) %>%
	mutate(pass.test = hgnc %in% 'NCAM1') %>%
	filter(any(pass.test == TRUE)) %>%
	arrange(pair.num, pair) %>%
	select(split = pair.num, gene = hgnc, tissue, level = tissue.max) %>%
	filter(!is.na(level))


combi.sub %>%
PlotTissue(facet = TRUE, pdf = TRUE, width = 10, height = 4.2, file.name = 'combi_level_1_NCAM1.pdf')









