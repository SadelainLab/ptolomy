
#--------
# imports
#--------


library(matrixStats)
library(scales)


#----------
# functions
#----------


RowToNames <- function(df, row.num = 1) {

	col.names <- df[row.num,]

	df %<>% set_names(col.names)

	df[-row.num,]
}

strip <- function(string) { gsub('^\\s+|\\s+$', '', string) }

FormatBlood <- function(blood) {

	blood %>%
	map( ~ {
		base <-
			.x %>%
			t %>%
			as_data_frame %>%
			RowToNames

		vals <-
			base %>%
			select(-1) %>%
			mutate_each(funs(2^as.numeric(.)))

		exp.max <-
			vals %>%
			as.matrix %>%
			rowMaxs

		hgnc <-
			names(base)[[1]] %>%
			str_split(' ') %>%
			unlist %>%
			head(1)

		master <-
			base %>%
			select(1) %>%
			set_names('cell') %>%
			mutate(cell = toupper(cell)) %>%
			mutate(hgnc = hgnc) %>%
			bind_cols(vals) %>%
			mutate(exp.max = exp.max) %>%
			arrange(hgnc, cell)

		master %>% select(cell, hgnc, exp.max, everything())
	}) %>%
	bind_rows %>%
	mutate(cell = ifelse(is.na(cell), 'nan', cell)) %>%
	select(cell, hgnc, exp = exp.max) %>%
	mutate(cell = ifelse(substr(strip(cell), 1, 3) == '7.0', 7, cell)) %>%
	mutate(cell = ifelse(substr(strip(cell), 1, 3) == '8.0', 8, cell)) %>%
	group_by(cell, hgnc) %>%
	mutate(mean.exp = mean(exp)) %>%
	ungroup %>%
	select(cell, hgnc, mean.exp) %>%
	unique %>%
	group_by(hgnc) %>%
	mutate(gene.order = mean(mean.exp)) %>%
	ungroup %>%
	arrange(desc(gene.order)) %>%
	select(-gene.order)
}

ExpPlot <- function(exp.table, plot.title, height = 4, width = 4, facet = TRUE, y.min = 0, y.max = 1, colors) {

	gg <-
		ggplot(exp.table , aes(x = hgnc, y = values, fill = color, width = 0.9)) +
		geom_bar(stat = 'identity', data = exp.table) +
		scale_y_continuous(limits = c(y.min, y.max)) +
		scale_fill_manual(values = colors)

	if(facet ==TRUE) {
		gg <- gg + facet_wrap(~ cell, ncol = 3)
	}

	ggt <-
		gg +
		theme(legend.title        = element_blank(),
			  panel.grid.major    = element_blank(),
			  panel.grid.minor    = element_blank(),
			  text                = element_text(size = 14),
			  axis.title.x        = element_blank(),
			  axis.title.y        = element_blank(),
			  axis.text.x         = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
			  axis.text.y         = element_text(size = 10),
			  legend.key          = element_rect(colour = 'white', fill = NULL, size = 0.1),
			  legend.key.size     = unit(1.4, 'lines'),
			  legend.text         = element_text(size = 10),
			  strip.text.x        = element_text(colour = 'black', size = 10),
			  strip.background    = element_rect(fill = 'white'),
			  plot.margin         = unit(c(1,1,1,1), 'cm'))

	pdf(str_c(gsub(' ', '_', plot.title), '.pdf'), height = height, width = width)
		plot(ggt)
	dev.off()
}


#----------
# variables
#----------


exclude.lines <- c('MPP', 'CMP', 'GMP', 'MEP', 'EARLY_PM', 'LATE_PM', 'MY', 'MM', 'BC', 'PMN', 'MONO', 'MDS', 'T(15;17)', 'ALL')


#----------
# read data
#----------


# bloodspot.normal.karyotype <-
# 	list.files('../data/bloodspot/normal_karyotype/', full.names = TRUE, pattern = '*csv') %>%
# 	map(~ { read.delim(.x, sep = ',', stringsAsFactors = FALSE, header = FALSE) })


bloodspot.pooled <-
	list.files('../data/bloodspot/pooled', full.names = TRUE, pattern = '*csv') %>%
	map(~ { read.delim(.x, sep = ',', stringsAsFactors = FALSE, header = FALSE) })


#----------
# read data
#----------

# normal karyotype

# blood.kf <-
# 	bloodspot.normal.karyotype %>%
# 	FormatBlood %>%
# 	mutate(
# 		hgnc = factor(hgnc, levels = unique(hgnc)),
# 		color = 'mean expression') %>%
# 	rename(values = mean.exp)

# blood.kf %>% ExpPlot(plot.title = 'Bloodspot normal karyotype RNA expression', height = 10, width = 20, y.min = 0, y.max = 15)


# pooled

blood.pf <-
	bloodspot.pooled %>%
	FormatBlood %>%
	filter(!cell %in% exclude.lines) %>%
	mutate(hgnc = case_when(
		.$hgnc %in% 'EMR2'     ~ 'ADGRE2',
		.$hgnc %in% 'ITFG3'    ~ 'FAM234A',
		.$hgnc %in% 'C16ORF54' ~ 'C16orf54',
		.$hgnc %in% 'C16ORF88' ~ 'KNOP1',
		TRUE ~ hgnc))

blood.pf.exp <-
	full_join(
		blood.pf %>% filter(!cell == 'HSC') %>% rename(cell.exp = mean.exp),
		blood.pf %>% filter(cell == 'HSC') %>% rename(hsc.exp = mean.exp) %>% select(-cell),
		by = 'hgnc') %>%
	rowwise %>%
	mutate(exp.ratio               = cell.exp / hsc.exp) %>%
	mutate(exp.ratio.log           = log10(exp.ratio)) %>%
	ungroup %>%
	group_by(hgnc) %>%
	mutate(gene.mean.exp.ratio     = mean(exp.ratio, na.rm = TRUE)) %>%
	mutate(gene.mean.exp.ratio.log = mean(exp.ratio.log, na.rm = TRUE)) %>%
	ungroup %>%
	# mutate(mean.exp.ratio          = mean(exp.ratio, na.rm = TRUE)) %>%
	# mutate(mean.exp.ratio.log      = mean(exp.ratio.log, na.rm = TRUE)) %>%
	mutate(mean.exp.ratio        = mean(exp.ratio, na.rm = TRUE)) %>%
	mutate(mean.exp.ratio.log    = mean(exp.ratio.log, na.rm = TRUE)) %>%
	mutate(sd.exp.ratio.log        = sd(exp.ratio.log, na.rm = TRUE)) %>%
	mutate(cutoff.log              = mean.exp.ratio.log + sd.exp.ratio.log) %>%
	mutate(cutoff                  = 10^cutoff.log) %>%
	ungroup %>%
	arrange(exp.ratio.log) %>%
	mutate(exp.ratio.num = row_number())

GetID((blood.pf.exp$hgnc %>% unique), 'symbol')




bloodspot.rna.cutoff <-
	ggplot(blood.pf.exp) +
	geom_point(aes(x = exp.ratio.num, y = exp.ratio), size = 0.5, alpha = 0.25) +
	scale_y_log10(breaks = pretty_breaks(9), limits = c(NA, 70)) +
	scale_x_continuous(breaks = c(0, nrow(blood.pf.exp))) +
	geom_rect(data = blood.pf.exp[1,], aes(xmin = -Inf, xmax = Inf, ymin = unique(cutoff), ymax = Inf), alpha = 0.2, fill = 'blue') +
	scale_alpha_manual(0.2) +
	geom_hline(aes(yintercept = unique(cutoff)), colour = 'blue', linetype = 2) +
	geom_hline(aes(yintercept = unique(mean.exp.ratio)), color = 'black', linetype = 2) +
	geom_text(data = data.frame(x = 0, y = unique(blood.pf.exp$cutoff)), aes(x, y), label = str_c('mean + 1 SD = ', round(unique(blood.pf.exp$cutoff), 2)), hjust = -0.1, vjust = -1) +
	geom_text(data = data.frame(x = 0, y = unique(blood.pf.exp$mean.exp.ratio)), aes(x, y), label = str_c('mean = ', round(unique(blood.pf.exp$mean.exp.ratio), 2)), hjust = -0.1, vjust = -1)

pdf('bloodspot_rna_cutoff.pdf', width = 6, height = 10)
	plot(bloodspot.rna.cutoff)
dev.off()

pf.join <-
	blood.pf.exp %>%
	left_join(ensembl.id.dict, by = 'hgnc') %>%
	select(-hgnc) %>%
	filter(exp.ratio.log >= cutoff.log) %>%
	inner_join(step.3, by = 'ensembl.gene') %>%
	arrange(mean.cor) %>%
	unique


ensembl.id.dict <- GetID(blood.pf.exp$hgnc %>% unique, 'symbol') %>% rename(hgnc = query)


message('done')
stop()



<- tryCatch(stop(), error = function(e) { message('hi') })



top.10 <- c('GAGE1', 'SLC2A9', 'C16orf54', 'CCR1', 'LILRB4', 'LTB4R', 'CD96', 'CLEC12A', 'CD84', 'P2RY13', 'SLC4A7')

color.4 <- c(`Excluded candidates` = 'grey', `Lowest expression in normal tissues` = 'blue', `Remaining positive correlations` = 'purple', `Known CAR targets` = 'orange')


pf.join %>%
select(hgnc, cell, mean.exp.ratio, exp.ratio, exp.ratio.log, gene.mean.exp.ratio, gene.mean.exp.ratio.log, prot.rna.cor, mean.cor) %>%
arrange(desc(gene.mean.exp.ratio)) %>%
mutate(hgnc = factor(hgnc, levels = unique(hgnc))) %>%
mutate(cell = ifelse(cell == 'NORMAL', 'normal karyotype', cell)) %>%
mutate(cell = ifelse(cell == 'COMPLEX', 'complex karyotype', cell)) %>%
mutate(cell = ifelse(cell == 'T(11Q23)/MLL', 'T(11q23)/MLL', cell)) %>%
mutate(cell = factor(cell, levels = unique(c('normal karyotype', 'complex karyotype', 'INV(16)', 'T(11q23)/MLL', cell)))) %>%
mutate(color =
	ifelse(hgnc %in% c('CLEC12A', 'FUT3'), 'Known CAR targets',
	ifelse(hgnc %in% top.10, 'Lowest expression in normal tissues',
	ifelse(is.na(prot.rna.cor) | prot.rna.cor <= 0, 'Excluded candidates',
	'Remaining positive correlations')))) %>%
mutate(color = factor(color, levels = unique(color))) %>%
rename(values = exp.ratio.log) %>%
select(hgnc, cell, values, color) %>%
ExpPlot(plot.title = 'Bloodspot RNA expression ratios', facet = TRUE, height = 15, width = 20, y.max = 3, colors = color.4)


pf.join %>%
select(hgnc, cell, mean.exp.ratio, exp.ratio, exp.ratio.log, gene.mean.exp.ratio, gene.mean.exp.ratio.log, prot.rna.cor, mean.cor) %>%
arrange(desc(prot.rna.cor)) %>%
mutate(hgnc = factor(hgnc, levels = unique(hgnc))) %>%
mutate(color =
	ifelse(hgnc %in% c('CLEC12A', 'FUT3'), 'Known CAR targets',
	ifelse(hgnc %in% top.10, 'Lowest expression in normal tissues',
	ifelse(prot.rna.cor <= 0, 'Excluded candidates',
	'Remaining positive correlations')))) %>%
mutate(color = factor(color, levels = unique(color))) %>%
rename(values = prot.rna.cor) %>%
select(hgnc, values, color) %>%
unique %>%
ExpPlot(plot.title = 'Bloodspot protein RNA correlation [cor sorted]', facet = FALSE, height = 6, width = 14, y.min = -0.5, y.max = 1, colors = color.4)



picks.29 <-
	pf.join %>%
	filter(!hgnc %in% c('GYPA', 'CDCA8')) %>%
	filter(mean.cor >= 0 | is.na(mean.cor))


picks.29 %>%
arrange(gross.mean.abundance) %>%
select(gene = hgnc, tissue, level = tissue.max) %>%
GeneRename %>%
PlotTissue(pdf = TRUE, width = 10, height = 5.6, file.name = 'picks_29_max.pdf')


picks.29 %>%
arrange(gross.mean.abundance) %>%
select(gene = hgnc, tissue, level = tissue.mean) %>%
GeneRename %>%
PlotTissue(pdf = TRUE, width = 10, height = 5.6, file.name = 'picks_29_mean.pdf')


exp.table <-
	bt.join %>%
	rowwise %>%
	mutate(mean.hsc.exp.2 = 2^mean.hsc.exp) %>%
	mutate(mean.cell.exp.2 = 2^mean.cell.exp) %>%
	mutate(exp.ratio.2 = mean.cell.exp.2 / mean.hsc.exp.2) %>%
	mutate(exp.ratio.2.log = log10(exp.ratio.2)) %>%
	ungroup %>%
	filter(hgnc %in% picks.29$hgnc) %>%
	select(hgnc, cell, exp.ratio, prot.rna.cor) %>%
	unique %>%
	filter(exp.ratio >= 1.2) %>%
	group_by(hgnc) %>%
	mutate(mean.exp.ratio = mean(exp.ratio)) %>%
	arrange(desc(mean.exp.ratio)) %>%
	ungroup %>%
	mutate(hgnc = factor(hgnc, levels = unique(hgnc))) %>%
	mutate(cell = ifelse(cell == 'NORMAL', 'normal karyotype', cell)) %>%
	mutate(cell = ifelse(cell == 'COMPLEX', 'complex karyotype', cell)) %>%
	mutate(cell = ifelse(cell == 'T(11Q23)/MLL', 'T(11q23)/MLL', cell)) %>%
	filter(!cell %in% c('MDS', 'ALL')) %>%
	mutate(cell = factor(cell, levels = unique(c('normal karyotype', 'complex karyotype', 'INV(16)', 'T(11q23)/MLL', cell)))) %>%
	mutate(color = 'blue') %>%
	mutate(color = factor(color, levels = unique(color))) %>%
	rename(values = exp.ratio)






exp.table %>% ExpPlot(plot.title = '29 Picks Bloodspot RNA expression ratios', height = 20, width = 20, y.max = 3)


#-----------------





exp.table %>% ExpPlot(plot.title = '23 Picks Bloodspot RNA expression ratios', height = 20, width = 20, y.max = 40)





exp.table <-
	bt.join %>%
	rowwise %>%
	mutate(mean.hsc.exp.2 = 2^mean.hsc.exp) %>%
	mutate(mean.cell.exp.2 = 2^mean.cell.exp) %>%
	mutate(exp.ratio.2 = mean.cell.exp.2 / mean.hsc.exp.2) %>%
	mutate(exp.ratio.2.log = log10(exp.ratio.2)) %>%
	ungroup %>%
	filter(hgnc %in% picks.23$hgnc) %>%
	select(hgnc, cell, exp.ratio.2, exp.ratio.2.log, prot.rna.cor) %>%
	unique %>%
	filter(exp.ratio.2 >= 1.2) %>%
	group_by(hgnc) %>%
	mutate(mean.exp.ratio.2 = mean(exp.ratio.2)) %>%
	arrange(desc(mean.exp.ratio.2)) %>%
	ungroup %>%
	mutate(hgnc = factor(hgnc, levels = unique(hgnc))) %>%
	mutate(cell = ifelse(cell == 'NORMAL', 'normal karyotype', cell)) %>%
	mutate(cell = ifelse(cell == 'COMPLEX', 'complex karyotype', cell)) %>%
	mutate(cell = ifelse(cell == 'T(11Q23)/MLL', 'T(11q23)/MLL', cell)) %>%
	filter(!cell %in% c('MDS', 'ALL')) %>%
	mutate(cell = factor(cell, levels = unique(c('normal karyotype', 'complex karyotype', 'INV(16)', 'T(11q23)/MLL', cell)))) %>%
	mutate(color = ifelse(hgnc %in% c('CCR1', 'P2RY13', 'ADGRE2', 'LTB4R'), 'black', 'red')) %>%
	mutate(color = factor(color, levels = unique(color))) %>%
	rename(values = exp.ratio.2.log)


exp.table %>% ExpPlot(plot.title = '23 Picks Bloodspot RNA expression ratios (log)', height = 20, width = 20, y.max = 2)



exp.table <-
	bt.join %>%
	rowwise %>%
	mutate(mean.hsc.exp.2 = 2^mean.hsc.exp) %>%
	mutate(mean.cell.exp.2 = 2^mean.cell.exp) %>%
	mutate(exp.ratio.2 = mean.cell.exp.2 / mean.hsc.exp.2) %>%
	mutate(exp.ratio.2.log = log10(exp.ratio.2)) %>%
	ungroup %>%
	filter(hgnc %in% picks.23$hgnc) %>%
	select(hgnc, cell, exp.ratio.2, exp.ratio.2.log, prot.rna.cor) %>%
	unique %>%
	filter(exp.ratio.2 >= 1.2) %>%
	group_by(hgnc) %>%
	mutate(mean.exp.ratio.2 = mean(exp.ratio.2)) %>%
	arrange(desc(mean.exp.ratio.2)) %>%
	ungroup %>%
	mutate(hgnc = factor(hgnc, levels = unique(hgnc))) %>%
	mutate(cell = ifelse(cell == 'NORMAL', 'normal karyotype', cell)) %>%
	mutate(cell = ifelse(cell == 'COMPLEX', 'complex karyotype', cell)) %>%
	mutate(cell = ifelse(cell == 'T(11Q23)/MLL', 'T(11q23)/MLL', cell)) %>%
	filter(!cell %in% c('MDS', 'ALL')) %>%
	mutate(cell = factor(cell, levels = unique(c('normal karyotype', 'complex karyotype', 'INV(16)', 'T(11q23)/MLL', cell)))) %>%
	mutate(color = ifelse(hgnc %in% c('CCR1', 'P2RY13', 'ADGRE2', 'LTB4R'), 'black', 'red')) %>%
	mutate(color = factor(color, levels = unique(color))) %>%
	rename(values = exp.ratio.2)


exp.table %>% ExpPlot(plot.title = '23 Picks Bloodspot RNA expression ratios', height = 20, width = 20, y.max = 40)

#-------------------------------


# palette
colors <- c(
"#38c76d",
"#e06ad6",
"#6ac35a",
"#9e7aed",
"#49a23f",
"#e957aa",
"#26a55e",
"#f54684",
"#46c693",
"#ee4d69",
"#38b094",
"#ee5c57",
"#2bc0cd",
"#f16944",
"#5e8aef",
"#73ba76",
"#c56cce",
"#519a5a",
"#ea88e2",
"#5c9766",
"#e35784",
"#7db986",
"#d76bac",
"#39adde",
"#d76548",
"#5c95df",
"#f58568",
"#a099e2",
"#d3785e",
"#b280d8",
"#ee9884",
"#a97fbd",
"#d26c69",
"#da99df",
"#c47569",
"#c576b5",
"#e66879",
"#f090bb",
"#e47e86",
"#f275a0",
"#cc7296")


zero.list <- c(
'SLC2A9',
'C16orf54',
'LILRB4',
'CLEC12A',
'SLC4A7',
'PIEZO1',
'ZNHIT2',
'CLEC11A',
'SLC19A1',
'KNOP1',
'TEX10',
'TLR2',
'SLC5A6',
'TTYH3',
'DOCK10',
'FAM234A',
'ABCC4',
'TRMT11',
'SLC7A6',
'EMB',
'ATP11A',
'LAIR1',
'MAEA',
'INTS7',
'PREX1',
'OSBPL3',
'ZZEF1',
'PLCB3',
'KIT',
'PIK3AP1',
'TOR2A',
'EPB41',
'CCDC88A',
'CBL',
'ENG',
'CD82',
'IL10RB',
'SLC35F2',
'SLC16A7',
'ABCC1',
'ARID2',
'TEX264',
'MFGE8',
'RABGAP1L',
'SLC6A6',
'ADGRE2',
'VCPIP1',
'DAGLB',
'CD70',
'ST14',
'NOTCH2',
'SLC31A1',
'SORT1',
'SLC2A3',
'NRD1',
'ITGA4',
'MTHFR',
'ATP13A1',
'FCGR1A',
'DOCK11',
'ASH2L',
'MMP14',
'DTNA')

zero.list <- c(
'GAGE1',
'SLC2A9',
'C16orf54',
'CCR1',
'LILRB4',
'LTB4R',
'CLEC12A',
'CD84',
'SLC4A7',
'SEMA4A',
'SDK2',
'PIEZO1',
'VWA3B',
'ZNHIT2',
'CLEC11A',
'PLXNC1',
'SLC19A1',
'KNOP1',
'CD209',
'GYPA',
'TEX10',
'TLR2',
'IFRD2',
'LTBR',
'SLC5A6',
'CARD11',
'ITGAX',
'TTYH3',
'DOCK10',
'FAM234A',
'GUCY2D',
'ABCC4',
'TRMT11',
'CDCA8',
'SLC7A6',
'EMB',
'ATP11A',
'LAIR1',
'MAEA',
'INTS7',
'PREX1',
'OSBPL3',
'ZZEF1',
'PLCB3',
'KIT',
'PIK3AP1',
'TOR2A',
'EPB41',
'CCDC88A',
'CBL',
'ENG',
'CD82',
'IL10RB',
'ANK1',
'LILRA6',
'SLC35F2',
'SLC16A7',
'ABCC1',
'ARID2',
'TEX264',
'PPP6R2',
'SLAMF9',
'MFGE8',
'RABGAP1L',
'SLC6A6',
'ADGRE2',
'TNFRSF1B',
'VCPIP1',
'PTPRG',
'DAGLB',
'CD70',
'ST14',
'EMC10',
'NOTCH2',
'SLC31A1',
'LILRB2',
'SORT1',
'SLC2A3',
'NRD1',
'ITGA4',
'MTHFR',
'ATP13A1',
'FCGR1A',
'DOCK11',
'ASH2L',
'MMP14',
'DIS3L',
'DTNA')





bar <-
foo %>%
#filter(hgnc %in% zero.list) %>%
filter(hgnc %in% one.list) %>%
# filter(exp.ratio.2 >= 4) %$% hgnc %>% unique %>% sort
select(hgnc, prot.rna.cor) %>%
unique


one.list <- c(
'CD96',
'P2RY13',
'FCAR',
'ICAM5',
'LILRA2',
'ITGB3',
'SIGLEC9')










