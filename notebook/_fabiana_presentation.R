

step.2 %>%
filter(ensembl.gene %in% picks) %>%
select(hgnc, gross.mean.abundance) %>%
arrange(gross.mean.abundance) %>%
mutate(hgnc = factor(hgnc, levels = unique(hgnc))) %>%
unique %>%
filter(!is.na(gross.mean.abundance)) %>%
mutate(tissue = 'gross mean') %>%
rename(gene = hgnc, level = gross.mean.abundance) %>%
PlotTissue(pdf = TRUE, file.name = '26_picks_gross_mean.pdf', width = 2.35, height = 10)


palette <- c(
'zero'   = '#f7d4d4',
'low'    = '#eb9494',
'medium' = '#de5454',
'high'   = '#c12525')


PlotTissue <- function(events, faceting = FALSE, pdf = FALSE, file.name = 'plot_tissue.pdf', width = 20, height = 10, order = TRUE) {

	if(order == TRUE) {
		events <-
			bind_rows(data_frame(
				gene = as.character(unlist(events[1,'gene'])),
				tissue = c('adipose tissue', 'adrenal', 'appendix', 'bladder', 'blood', 'bone', 'brain', 'breast', 'bronchus', 'cerumen', 'cervix', 'epididymis', 'eye', 'fallopian tube', 'gallbladder', 'gut', 'heart', 'kidney', 'laryngopharynx', 'liver', 'lung', 'lymph node', 'nasopharynx', 'oropharynx', 'ovary', 'pancreas', 'parathyroid', 'prostate', 'rectum', 'seminal', 'skeletal muscle', 'skin', 'smooth muscle', 'soft tissue', 'spinal cord', 'spleen', 'stomach', 'synovial fluid', 'testis', 'thyroid', 'tonsil', 'uterus', 'vagina'),
				level = NA), events) %>%
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
		#getOption('device')()
		plot(mt.gg)
	}
}


line.order = c(
'msk.0.09aml',
'msk.0.tf',
'jpro.thp1',
'msk.0.thp1',
'msk.1.thp1',
'msk.2.thp1',
'msk.3.thp1',
'msk.0.kasum1',
'msk.1.kasumi',
'msk.2.kasumi',
'msk.3.kasumi',
'msk.0.monomac',
'msk.1.monomac',
'msk.2.monomac',
'msk.3.monomac',
'msk.0.molm13',
'msk.1.molm13',
'msk.2.molm13',
'msk.3.molm13')


step.2 %>%
filter(ensembl.gene %in% picks) %>%
arrange(gross.mean.abundance) %>%
select(hgnc,
msk.0.09aml,
msk.0.kasum1,
msk.0.molm13,
msk.0.monomac,
msk.0.tf,
msk.0.thp1,
msk.1.kasumi,
msk.2.kasumi,
msk.3.kasumi,
msk.1.thp1,
msk.2.thp1,
msk.3.thp1,
msk.1.monomac,
msk.2.monomac,
msk.3.monomac,
msk.1.molm13,
msk.2.molm13,
msk.3.molm13,
jpro.thp1) %>%
unique %>%
mutate(hgnc = factor(hgnc, levels = unique(hgnc))) %>%
unique %>%
gather(tissue, level, msk.0.09aml:jpro.thp1) %>%
arrange(tissue) %>%
mutate(tissue = factor(tissue, line.order)) %>%
rename(gene = hgnc) %>%
filter(!is.na(tissue)) %>%
filter(!is.na(level)) %>%
select(gene, level, tissue) %>%
PlotTissue(pdf = TRUE, file.name = 'malignant_picks.pdf', width = 9, height = 10, order = FALSE)


step.2 %>%
filter(ensembl.gene %in% picks) %>%
arrange(gross.mean.abundance) %>%
select(hgnc,
msk.1.kasumi,
msk.2.kasumi,
msk.3.kasumi,
msk.1.thp1,
msk.2.thp1,
msk.3.thp1,
msk.1.monomac,
msk.2.monomac,
msk.3.monomac,
msk.1.molm13,
msk.2.molm13,
msk.3.molm13) %>%
unique %>%
mutate(hgnc = factor(hgnc, levels = unique(hgnc))) %>%
unique %>%
gather(tissue, level, msk.1.kasumi:msk.3.molm13) %>%
rename(gene = hgnc) %>%
filter(!is.na(tissue)) %>%
filter(!is.na(level)) %>%
select(gene, level, tissue) %>%
PlotTissue(pdf = TRUE, file.name = 'triplicate_malignant_picks.pdf', width = 7, height = 10, order = FALSE)


step.2 %>%
filter(ensembl.gene %in% picks) %>%
arrange(gross.mean.abundance) %>%
rowwise %>%
mutate(msk.kasumi = mean(c(msk.1.kasumi, msk.2.kasumi, msk.3.kasumi))) %>%
mutate(msk.thp1 = mean(c(msk.1.thp1, msk.2.thp1, msk.3.thp1))) %>%
mutate(msk.monomac = mean(c(msk.1.monomac, msk.2.monomac, msk.3.monomac))) %>%
mutate(msk.molm13 = mean(c(msk.1.molm13, msk.2.molm13, msk.3.molm13))) %>%
select(hgnc,
msk.kasumi,
msk.thp1,
msk.monomac,
msk.molm13) %>%
unique %>%
mutate(hgnc = factor(hgnc, levels = unique(hgnc))) %>%
unique %>%
gather(tissue, level, msk.kasumi:msk.molm13) %>%
rename(gene = hgnc) %>%
filter(!is.na(tissue)) %>%
filter(!is.na(level)) %>%
select(gene, level, tissue) %>%
PlotTissue(pdf = TRUE, file.name = 'triplicate_collapsed_malignant_picks.pdf', width = 3.55, height = 10, order = FALSE)


step.2 %>%
filter(ensembl.gene %in% picks) %>%
arrange(gross.mean.abundance) %>%
rowwise %>%
mutate(level = mean(c(msk.1.kasumi, msk.2.kasumi, msk.3.kasumi, msk.1.thp1, msk.2.thp1, msk.3.thp1, msk.1.monomac, msk.2.monomac, msk.3.monomac, msk.1.molm13, msk.2.molm13, msk.3.molm13))) %>%
select(hgnc,
level) %>%
unique %>%
mutate(hgnc = factor(hgnc, levels = unique(hgnc))) %>%
unique %>%
rename(gene = hgnc) %>%
filter(!is.na(level)) %>%
mutate(tissue = 'malignant mean') %>%
select(gene, level, tissue) %>%
PlotTissue(pdf = TRUE, file.name = 'triplicate_mean_malignant_picks.pdf', width = 2.2, height = 10, order = FALSE)


micro <-
	read_tsv('micro_patient.txt') %>%
	rename(tissue = patient) %>%
	gather(gene, level, GAGE1:MMP14) %>%
	filter(gene %in% c(names(picks), 'EMR2', 'GPR86')) %>%
	mutate(level = level - min(level)) %>%
	mutate(level = level * (3/max(level))) %>%
	group_by(gene) %>%
	mutate(rank = mean(level)) %>%
	arrange(rank) %>%
	ungroup %>%
	mutate(gene = factor(gene, unique(gene)))

micro %>%
	PlotTissue(pdf = TRUE, file.name = 'micro_array_26.pdf', width = 15, height = 5, order = FALSE)

foo <- c(
"ABCC4",
"ANK1",
"ARID2",
"ATP11A",
"CBL",
"CCDC88A",
"CCR1",
"CD209",
"CD84",
"CD96",
"DOCK10",
"DOCK11",
"DTNA",
"ENG",
"EPB41",
"FCAR",
"GYPA",
"ITGA4",
"ITGB3",
"KIT",
"LILRA6",
"LILRB2",
"LILRB4",
"MTHFR",
"NOTCH2",
"PLXNC1",
"RABGAP1L",
"SIGLEC9",
"SLC16A7",
"SLC2A9",
"SLC31A1",
"SLC4A7",
"SORT1",
"ST14",
"VCPIP1",
"ZZEF1")


rna <-
	read_tsv('rna_seq_08242015.txt') %>%
	rename(DNMT3a_mut = `DNMT3a mut`) %>%
	group_by(gene) %>%
	mutate(DNMT3a_mut = mean(DNMT3a_mut), s_DNMT3a_WT = mean(s_DNMT3a_WT), s_MIGR1 = mean(s_MIGR1)) %>%
	unique %>%
	ungroup %>%
	gather(tissue, level, DNMT3a_mut:s_MIGR1) %>%
	mutate(level = log10(level)) %>%
	mutate(level = level * 3/max(na.omit(level))) %>%
	filter(gene %in% foo) %>%
	group_by(gene) %>%
	mutate(rank = mean(level)) %>%
	arrange(rank) %>%
	ungroup %>%
	mutate(gene = factor(gene, unique(gene)))

rna %>%
	PlotTissue(pdf = TRUE, file.name = 'rna_seq_36.pdf', width = 2.65, height = 10, order = FALSE)






