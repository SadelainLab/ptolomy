
#-----------

step.4 <-
	combn(1:length(unique(step.3$hgnc)), 2) %>%
	t %>%
	as_data_frame %>%
	set_names(c('a', 'b')) %>%
	mutate(pair.num = row_number()) %>%
	rowwise %>%
	do({
			x <- .
			bind_cols(
				step.3 %>% filter(group.n == x$a) %>% set_names(str_c(colnames(step.3), '_a')),
				step.3 %>% filter(group.n == x$b) %>% set_names(str_c(colnames(step.3), '_b'))
			) %>%
			mutate(pair.num = x$pair.num) %>%
			select(pair.num, everything())
	}) %>%
	ungroup %>%
	filter(!is.na(ensembl.gene_a) & !is.na(ensembl.gene_b))



#-------


step.5 <-
	step.4 %>%
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




step.5.level.1 <-
	step.5 %>%
	group_by(pair.num) %>%
	filter(all(max.level.1 == TRUE)) %>%
	mutate(order = mean(gross.mean.abundance)) %>%
	ungroup %>%
	arrange(order) %>%
	#filter(!tissue %in% c('cerumen', 'synovial fluid')) %>%
	mutate(pair.num = factor(pair.num, levels = unique(pair.num))) %>%
	mutate(tissue = factor(tissue, levels = c(vital, non.vital))) %>%
	group_by(pair.num) %>%
	mutate(pass.test = !hgnc %in% pass.protein) %>%
	filter(any(pass.test == TRUE)) %>%
	select(split = pair.num, gene = hgnc, tissue, level = tissue.max) %>%
	filter(!is.na(level))

step.5.level.2 <-
	step.5 %>%
	group_by(pair.num) %>%
	filter(all(max.level.2 == TRUE)) %>%
	mutate(order = mean(gross.mean.abundance)) %>%
	ungroup %>%
	arrange(order) %>%
	#filter(!tissue %in% c('cerumen', 'synovial fluid')) %>%
	mutate(pair.num = factor(pair.num, levels = unique(pair.num))) %>%
	mutate(tissue = factor(tissue, levels = c(vital, non.vital))) %>%
	group_by(pair.num) %>%
	mutate(pass.test = !hgnc %in% pass.protein) %>%
	filter(any(pass.test == TRUE)) %>%
	select(split = pair.num, gene = hgnc, tissue, level = tissue.max) %>%
	filter(!is.na(level))



step.5.level.1 %>% PlotTissue(faceting = TRUE, pdf = TRUE, width = 10, height = 400, file.name = 'combi_level_1.pdf')

step.5.level.2 %>% PlotTissue(faceting = TRUE, pdf = TRUE, width = 10, height = 80, file.name = 'combi_level_2.pdf')

