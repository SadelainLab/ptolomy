#-------------#
#             #
#  retreival  #
#             #
#-------------#


#--------
# imports
#--------


#load('clean.RData')
source('_keys.R')


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

palette <- c(
'zero'   = '#202c99',
'low'    = '#fbeaea',
'medium' = '#e77e7e',
'high'   = '#c12525')


#----------
# functions
#----------


PlotTissue <- function(events, faceting = FALSE, pdf = FALSE, file.name = 'plot_tissue.pdf', width = 20, height = 10, order = TRUE) {

	events %<>% unique

	if(order == TRUE) {
		events %<>%
			bind_rows(data_frame(
				gene = as.character(unlist(events[1,'gene'])),
				tissue = c('adipose tissue', 'adrenal', 'appendix', 'bladder', 'blood', 'bone', 'brain', 'breast', 'bronchus',
						   'cerumen', 'cervix', 'epididymis', 'eye', 'fallopian tube', 'gallbladder', 'gut', 'heart', 'kidney',
						   'laryngopharynx', 'liver', 'lung', 'lymph node', 'nasopharynx', 'oropharynx', 'ovary', 'pancreas',
						   'parathyroid', 'prostate', 'rectum', 'seminal', 'skeletal muscle', 'skin', 'smooth muscle',
						   'soft tissue', 'spinal cord', 'spleen', 'stomach', 'synovial fluid', 'testis', 'thyroid', 'tonsil',
						   'uterus', 'vagina'), level = NA)) %>%
			mutate(gene = factor(gene, levels = unique(gene))) %>%
			mutate(tissue = factor(tissue, levels = c(vital, non.vital)))
	}

	if(all(na.omit(events$level) %% 1 == 0)) {  # check if integer, if so plot discrete
		events %<>%
			mutate(level =
				ifelse(level == 0, 'zero',
				ifelse(level == 1, 'low',
				ifelse(level == 2, 'medium',
				ifelse(level == 3, 'high',
				NA))))) %>%
			mutate(level = factor(level, levels = unique(level))) %>%
			arrange(desc(is.na(level)))

		m.gg <-
			ggplot(events, aes(tissue, gene)) +
			geom_tile(aes(fill = level), colour = 'grey') +
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

	if(faceting == TRUE) {
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
				legend.title     = element_blank(),
				axis.title.x     = element_blank(),
				axis.title.y     = element_blank(),
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
		dev.new(width = width, height = height)
		plot(mt.gg)
	}
}


GeneRename <- function(events) {
	events %>%
	mutate(gene =
		ifelse(gene == 'ENG', 'ENG (CD105)',
		ifelse(gene == 'CLEC12A', 'CLEC12A (CLL-1)',
		ifelse(gene == 'ADGRE2', 'ADGRE2 (EMR2)',
		ifelse(gene == 'LAIR1', 'LAIR1 (CD305)',
		ifelse(gene == 'TLR2', 'TLR2 (CD282)',
		ifelse(gene == 'P2RY13', 'P2RY13 (GPR86)',
		ifelse(gene == 'SLC19A1', 'SLC19A1 (FOLT)',
		gene))))))))
}


#-----------------
# AML surface pool
#-----------------


msk.jpro.add <-
	c(msk.surface, jpro.surface, aml.additions) %>%
	unique %>%
	sort


#--------------------
# candidate selection
#--------------------

step.1 <-
	step.0 %>%
	filter(ensembl.gene %in% msk.jpro.add) %>%
	filter(db.num >= 2) %>%
	filter(membrane == TRUE) %>%
	filter(!hgnc %in% exclude.protein) %>%
	arrange(prot.mean, hgnc) %>%
	select(-membrane) %>%
	mutate(prot.mean.cut = step.1 %>% select(ensembl.gene, prot.mean) %>% unique %$% prot.mean %>% sort %>% mean(na.rm = TRUE) )


# candidate comparison  -- rank in order of quality




prot.mean.cut.plot <-
	step.1 %>%
	select(ensembl.gene, prot.mean, prot.mean.cut) %>%
	unique %>%
	arrange(prot.mean) %>%
	mutate(prot.mean.num = row_number())

prot.mean.cut.gg <-
	ggplot(prot.mean.cut.plot) +
	geom_point(aes(x = prot.mean.num, y = prot.mean), size = 0.5, alpha = 0.25) +
	scale_y_continuous(breaks = pretty_breaks(9), limits = c(0, 3)) +
	scale_x_continuous(breaks = c(0, nrow(prot.mean.cut.plot))) +
	geom_rect(data = prot.mean.cut.plot[1,], aes(xmin = -Inf, xmax = Inf, ymin = unique(prot.mean.cut), ymax = Inf), alpha = 0.2, fill = 'blue') +
	scale_alpha_manual(0.2) +
	geom_hline(aes(yintercept = unique(prot.mean.cut)), colour = 'blue', linetype = 2) +
	geom_text(
		data = data.frame(x = 0, y = unique(prot.mean.cut.plot$prot.mean.cut)),
		aes(x, y),
		label = str_c('mean = ', round(unique(prot.mean.cut.plot$prot.mean.cut), 2)),
		hjust = -0.1,
		vjust = -1)

pdf('gross_mean_mean_cutoff.pdf', width = 6, height = 10)
	plot(prot.mean.cut.gg)
dev.off()




step.2 <-
	step.1 %>%
	filter(prot.mean < prot.mean.cut) %>% # | hgnc %in% pass.protein) %>%
	group_by(ensembl.gene) %>%
	mutate(max.test = tissue.max < 3 | tissue %in% pass.tissue) %>%
	filter(!any(na.omit(max.test) == FALSE)) %>% # | hgnc   %in% pass.protein) %>%
	mutate(aml.tot =
		any(!is.na(msk.0.09aml)   & msk.0.09aml   != 0) +
		any(!is.na(msk.0.kasum1)  & msk.0.kasum1  != 0) +
		any(!is.na(msk.0.molm13)  & msk.0.molm13  != 0) +
		any(!is.na(msk.0.monomac) & msk.0.monomac != 0) +
		any(!is.na(msk.0.tf)      & msk.0.tf      != 0) +
		any(!is.na(msk.0.thp1)    & msk.0.thp1    != 0) +
		any(!is.na(jpro.thp1)     & jpro.thp1     != 0)
	) %>%
	ungroup


	arrange(prot.mean, ensembl.gene) %>%
	full_join(expand.grid(hgnc = unique(.$hgnc), tissue = unique(.$tissue), stringsAsFactors = FALSE), by = c('hgnc', 'tissue')) %>%
	mutate(group.n = as.integer(factor(hgnc, levels=unique(hgnc)))) %>%
	arrange(hgnc, tissue) #%>%
	#filter(!hgnc %in% exclude.protein)

step.3 <-
	step.2 %>%
	filter(!hgnc %in% exclude.protein) %>%
	group_by(hgnc) %>%
	mutate(prot.zero.test = !all(c(
		msk.0.09aml   <= 0,
		msk.0.kasum1  <= 0,
		msk.0.molm13  <= 0,
		msk.0.monomac <= 0,
		msk.0.tf      <= 0,
		msk.0.thp1    <= 0,
		msk.1.kasumi  <= 0,
		msk.2.kasumi  <= 0,
		msk.3.kasumi  <= 0,
		msk.1.thp1    <= 0,
		msk.2.thp1    <= 0,
		msk.3.thp1    <= 0,
		msk.1.monomac <= 0,
		msk.2.monomac <= 0,
		msk.3.monomac <= 0,
		msk.1.molm13  <= 0,
		msk.2.molm13  <= 0,
		msk.3.molm13  <= 0,
		jpro.thp1     <= 0),
	na.rm = TRUE)) %>%
	ungroup %>%
	filter(prot.zero.test != FALSE)


step.2 %>%
select(-tissue.mean, -hpa, -pdb, -hpm, -rna, -tissue.mean.round, -max.test) %>%
spread(tissue, tissue.max) %>%
unique %>%
#filter(!hgnc %in% pass.protein) %>%
write.xlsx('step_2.xlsx')


# step.2.picks.pass <-
# 	step.2 %>%
# 	filter(hgnc %in% names(picks) | hgnc %in% pass.protein)

# step.2.picks <-
# 	step.2 %>%
# 	filter(hgnc %in% names(picks))

step.2 %>%
filter(!hgnc %in% pass.protein) %>%
arrange(gross.mean.abundance) %>%
select(gene = hgnc, tissue, level = tissue.max) %>%
GeneRename %>%
mutate(gene = factor(gene, levels = unique(gene))) %>%
mutate(tissue = factor(tissue, levels = c(vital, non.vital))) %>%
PlotTissue(pdf = TRUE, width = 10, height = 29, file.name = 'step_2_max.pdf')

step.2 %>%
filter(!hgnc %in% pass.protein) %>%
arrange(gross.mean.abundance) %>%
select(gene = hgnc, tissue, level = tissue.mean) %>%
GeneRename %>%
mutate(gene = factor(gene, levels = unique(gene))) %>%
mutate(tissue = factor(tissue, levels = c(vital, non.vital))) %>%
PlotTissue(pdf = TRUE, width = 10, height = 29, file.name = 'step_2_mean.pdf')

message('done')
stop()

step.2 %>%
filter(ensembl.gene %in% picks) %>%
filter(!hgnc %in% pass.protein) %>%
arrange(gross.mean.abundance) %>%
select(gene = hgnc, tissue, level = tissue.max) %>%
GeneRename %>%
mutate(gene = factor(gene, levels = unique(gene))) %>%
mutate(tissue = factor(tissue, levels = c(vital, non.vital))) %>%
PlotTissue(pdf = TRUE, width = 10, height = 5.6, file.name = 'step_2_picks_max.pdf')

step.2 %>%
filter(ensembl.gene %in% picks) %>%
filter(!hgnc %in% pass.protein) %>%
arrange(gross.mean.abundance) %>%
select(gene = hgnc, tissue, level = tissue.mean) %>%
GeneRename %>%
mutate(gene = factor(gene, levels = unique(gene))) %>%
mutate(tissue = factor(tissue, levels = c(vital, non.vital))) %>%
PlotTissue(pdf = TRUE, width = 10, height = 5.6, file.name = 'step_2_picks_mean.pdf')








