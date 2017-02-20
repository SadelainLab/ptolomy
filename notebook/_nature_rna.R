nature <-
	read.delim('nature.txt', stringsAsFactors = FALSE, header = FALSE) %>%
	unlist %>%
	unname


# in nature
 [1] "ADAM19"   "ADM"      "AHSP"     "AIF1L"    "AKR1C3"   "ALAS2"    "ANGPT1"   "AQP9"     "ARHGAP22" "ATP8B4"   "BCL6"     "BIVM"     "CCDC109B"
[14] "CD14"     "CD247"    "CD34"     "CD48"     "CDK6"     "CECR1"    "CHST15"   "CKAP4"    "COL24A1"  "CPXM1"    "CTSH"     "CXCL16"   "CXCR4"
[27] "DNMT3B"   "DPYSL3"   "E2F2"     "EMP1"     "EPDR1"    "FAM30A"   "FAM69B"   "FCGR2A"   "FCN1"     "FCRL3"    "FCRLA"    "FGR"      "FLT3"
[40] "FSCN1"    "GATA2"    "GIMAP4"   "GNLY"     "GPSM1"    "GUCY1A3"  "GZMA"     "GZMB"     "GZMH"     "H2AFY2"   "HBA1"     "HBA2"     "HBB"
[53] "HBM"      "HOXA5"    "HOXA6"    "HOXA9"    "IL18RAP"  "IL2RB"    "ITPR3"    "JAZF1"    "KCNK17"   "KLRB1"    "LAPTM4B"  "LILRA5"   "MAMDC2"
[66] "MMRN1"    "MTSS1"    "MYCN"     "NPL"      "NYNRIN"   "RBM38"    "SHANK3"   "SLC15A3"  "SLC25A37" "SLC4A1"   "SLC7A7"   "SOCS2"    "SPINK2"
[79] "SPNS2"    "TFPI"     "TNFRSF4"  "VWF"      "ZBTB46"

# not in nature
 [1] "AIMZ"      "C1Oorf140" "C1orf150"  "C3orf54"   "CTSL1"     "FLJ22662"  "GPR56"     "IL1ORA"    "ISG2O"     "KIAAO125"  "LOC284422" "LOC642113"
[13] "LOC647450" "LOC647506" "LOC652493" "LOC652694" "LOC654103" "NGFRAP1"   "PRSSL1"    "RAGE"      "SGK"

# nature lookups
 [1] "ADAMTS13" "ADGRG1"   "AGER"     "BEX3"     "CTSL"     "FAM212A"  "GCSAML"   "KIAA0125" "KLK10"    "LILRP2"   "LRP5"     "MCUB"     "MOK"
[14] "MYCNOS"   "NDUFA2"   "PRSS57"   "SGK1"     "TFPI2"

nature.all <-
	GetEnsemblID(nature) %$% ensembl.gene %>% GetHugoID %$% hgnc %>% c(nature) %>% unique %>% sort



nature.s1 <-
	step.0 %>%
	filter(db.num >= 2) %>%
	filter(membrane == TRUE) %>%
	filter(hgnc %in% nature.all) %>%
	arrange(gross.mean.abundance, hgnc) %>%
	select(
		hgnc, ensembl.gene, tissue, db.num, fill, hpa, hpm, pdb,
		tissue.max, tissue.mean, gross.mean.abundance, rna, everything(), -membrane)


nature.s1 %>% select(hgnc) %>% unique %>% unlist %>% unname %>% sort

 [1] "ADAM19"   "ADAMTS13" "ADGRG1"   "AGER"     "AIF1L"    "ANGPT1"   "AQP9"     "ATP8B4"   "BIVM"     "CCDC109B" "CD14"     "CD247"    "CD34"
[14] "CD48"     "COL24A1"  "CTSL"     "CXCL16"   "CXCR4"    "DPYSL3"   "EPDR1"    "FCGR2A"   "FCN1"     "FCRL3"    "FCRLA"    "FGR"      "FLT3"
[27] "FSCN1"    "GNLY"     "GZMA"     "GZMB"     "GZMH"     "HBA1"     "HBA2"     "HBB"      "IL18RAP"  "IL2RB"    "ITPR3"    "KLK10"    "KLRB1"
[40] "LILRA5"   "LRP5"     "MAMDC2"   "MMRN1"    "MOK"      "NPL"      "SGK1"     "SLC25A37" "SLC7A7"   "SOCS2"    "SPINK2"   "TFPI"     "TFPI2"
[53] "VWF"





nature.s2 <-
	nature.s1 %>%
	filter(gross.mean.abundance < 1) %>%
	group_by(hgnc) %>%
	mutate(tissue.mean.round =
		ifelse(tissue.mean <  0.5, 0,
		ifelse(tissue.mean <  1.5, 1,
		ifelse(tissue.mean <  2.5, 2,
		ifelse(tissue.mean <= 3,   3,
			NA))))) %>%
	mutate(max.test = tissue.max < 3 | tissue %in% pass.tissue) %>%
	filter(!any(na.omit(max.test) == FALSE)) %>%
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
	arrange(gross.mean.abundance) %>%
	full_join(expand.grid(hgnc = unique(.$hgnc), tissue = unique(.$tissue), stringsAsFactors = FALSE), by = c('hgnc', 'tissue')) %>%
	mutate(group.n = as.integer(factor(hgnc, levels=unique(hgnc)))) %>%
	arrange(hgnc, tissue)

nature.s2 %>% select(hgnc) %>% unique %>% unlist %>% unname %>% sort

 [1] "ADAM19"   "AQP9"     "BIVM"     "CCDC109B" "CD48"     "FCRL3"    "GNLY"     "GZMB"     "GZMH"     "IL18RAP"  "IL2RB"    "LILRA5"   "MMRN1"
[14] "SLC7A7"  { "TFPI2" }










step.2 %>% select(hgnc) %>% unique %>% unlist %>% unname -> msk.96

# step 0 nature

  [1] "CD34"      "SPINK2"    "LAPTM4B"   "HOXA5"     "GUCY1A3"   "SHANK3"    "ANGPT1"    "ARHGAP22"  "LOC284422" "MYCN"      "MAMDC2"    "PRSSL1"
 [13] "KIAAO125"  "GPSM1"     "HOXA9"     "MMRN1"     "FSCN1"     "DNMT3B"    "HOXA6"     "AIF1L"     "SOCS2"     "CDK6"      "FAM69B"    "NGFRAP1"
 [25] "C3orf54"   "CPXM1"     "TNFRSF4"   "ZBTB46"    "DPYSL3"    "NYNRIN"    "COL24A1"   "FAM30A"    "C1Oorf140" "SPNS2"     "GPR56"     "AKR1C3"
 [37] "FLT3"      "TFPI"      "KCNK17"    "EPDR1"     "C1orf150"  "BIVM"      "H2AFY2"    "VWF"       "EMP1"      "RAGE"      "ATP8B4"    "GATA2"
 [49] "SLC25A37"  "SGK"       "LOC652694" "ITPR3"     "LOC654103" "CXCR4"     "FCRL3"     "RBM38"     "LILRA5"    "IL18RAP"   "CCDC109B"  "ISG2O"
 [61] "MTSS1"     "CECR1"     "ADAM19"    "FCGR2A"    "AIMZ"      "NPL"       "IL1ORA"    "CTSL1"     "GNLY"      "CKAP4"     "ADM"       "KLRB1"
 [73] "SLC15A3"   "FGR"       "FCRLA"     "IL2RB"     "CXCL16"    "SLC4A1"    "GZMH"      "FLJ22662"  "LOC647506" "GIMAP4"    "JAZF1"     "CTSH"
 [85] "GZMA"      "CHST15"    "AQP9"      "CD247"     "BCL6"      "SLC7A7"    "E2F2"      "LOC647450" "GZMB"      "LOC652493" "HBM"       "CD14"
 [97] "ALAS2"     "HBB"       "LOC642113" "AHSP"      "FCN1"      "CD48"      "HBA2"      "HBA1"

# step 1 nature

nature.s1$hgnc %>% unique -> s1.nature

 [1] "ADAM19"   "AIF1L"    "ANGPT1"   "AQP9"     "ATP8B4"   "BIVM"     "CCDC109B" "CD14"     "CD247"    "CD34"     "CD48"     "COL24A1"  "CXCL16"
[14] "CXCR4"    "DPYSL3"   "EPDR1"    "FCGR2A"   "FCN1"     "FCRL3"    "FCRLA"    "FGR"      "FLT3"     "FSCN1"    "GNLY"     "GZMA"     "GZMB"
[27] "GZMH"     "HBA1"     "HBA2"     "HBB"      "IL18RAP"  "IL2RB"    "ITPR3"    "KLRB1"    "LILRA5"   "MAMDC2"   "MMRN1"    "NPL"      "SLC25A37"
[40] "SLC7A7"   "SOCS2"    "SPINK2"   "TFPI"     "VWF"

# step 2 nature

nature.s2$hgnc %>% unique -> s2.nature

[1] "ADAM19"   "AQP9"     "BIVM"     "CCDC109B" "CD48"     "FCRL3"    "GNLY"     "GZMB"     "GZMH"     "IL18RAP"  "IL2RB"    "LILRA5"   "MMRN1"
[14] "SLC7A7"

# our picks

msk.96

 [1] "ABCC1"    "ABCC4"    "ADGRE2"   "ANK1"     "ARID2"    "ASH2L"    "ATP11A"   "ATP13A1"  "C16orf54" "CARD11"   "CBL"      "CCDC88A"  "CCR1"
[14] "CD209"    "CD70"     "CD82"     "CD84"     "CD96"     "CDCA8"    "CLEC11A"  "CLEC12A"  "DAGLB"    "DIS3L"    "DOCK10"   "DOCK11"   "DTNA"
[27] "EMB"      "EMC10"    "ENG"      "EPB41"    "FAM234A"  "FCAR"     "FCGR1A"   "FUT3"     "GAGE1"    "GUCY2D"   "GYPA"     "ICAM5"    "IFRD2"
[40] "IL10RB"   "INTS7"    "ITGA4"    "ITGAX"    "ITGB3"    "KIT"      "KNOP1"    "LAIR1"    "LILRA2"   "LILRA6"   "LILRB2"   "LILRB4"   "LTB4R"
[53] "LTBR"     "MAEA"     "MFGE8"    "MMP14"    "MTHFR"    "NOTCH2"   "NRD1"     "OSBPL3"   "P2RY13"   "PIEZO1"   "PIK3AP1"  "PLCB3"    "PLXNC1"
[66] "PPP6R2"   "PREX1"    "PTPRG"    "RABGAP1L" "SDK2"     "SEMA4A"   "SIGLEC9"  "SLAMF9"   "SLC16A7"  "SLC19A1"  "SLC2A3"   "SLC2A9"   "SLC31A1"
[79] "SLC35F2"  "SLC4A7"   "SLC5A6"   "SLC6A6"   "SLC7A6"   "SORT1"    "ST14"     "TEX10"    "TEX264"   "TLR2"     "TNFRSF1B" "TOR2A"    "TRMT11"
[92] "TTYH3"    "VCPIP1"   "VWA3B"    "ZNHIT2"   "ZZEF1"

















dnt.mut <-
	read.delim('../tables/DNMT3A_mutant_proteomics.txt', stringsAsFactors = FALSE, header = FALSE) %>%
	unlist %>%
	unname %>%
	unique

dnt.mut %<>% GetHugoID %$% hgnc %>% unique %>% list.filter(!is.na(.)) %>% sort



dnt.s1 <-
	step.0 %>%
	filter(db.num >= 2) %>%
	filter(membrane == TRUE) %>%
	filter(hgnc %in% dnt.mut) %>%
	arrange(gross.mean.abundance, hgnc) %>%
	select(
		hgnc, ensembl.gene, tissue, db.num, fill, hpa, hpm, pdb,
		tissue.max, tissue.mean, gross.mean.abundance, rna, everything(), -membrane)


dnt.s1 %>% select(hgnc) %>% unique %>% unlist %>% unname %>% sort








dnt.s2 <-
	dnt.s1 %>%
	filter(gross.mean.abundance < 1) %>%
	group_by(hgnc) %>%
	mutate(tissue.mean.round =
		ifelse(tissue.mean <  0.5, 0,
		ifelse(tissue.mean <  1.5, 1,
		ifelse(tissue.mean <  2.5, 2,
		ifelse(tissue.mean <= 3,   3,
			NA))))) %>%
	mutate(max.test = tissue.max < 3 | tissue %in% pass.tissue) %>%
	filter(!any(na.omit(max.test) == FALSE)) %>%
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
	arrange(gross.mean.abundance) %>%
	full_join(expand.grid(hgnc = unique(.$hgnc), tissue = unique(.$tissue), stringsAsFactors = FALSE), by = c('hgnc', 'tissue')) %>%
	mutate(group.n = as.integer(factor(hgnc, levels=unique(hgnc)))) %>%
	arrange(hgnc, tissue)

dnt.s2 %>% select(hgnc) %>% unique %>% unlist %>% unname %>% sort

 [1] "ABCC1"    "ABCC4"    "ADGRB1"   "ARID2"    "C16orf54" "CARD9"    "CCDC88A"  "CD84"     "CDK2"     "CLEC4A"   "DDX51"    "DDX54"    "DNAAF5"
[14] "DOCK11"   "DSC1"     "ENG"      "EPB41"    "FANCD2"   "GNL2"     "GRK5"     "INTS3"    "ITGA2B"   "ITGA4"    "LAIR1"    "LRRC41"   "MMP14"
[27] "NRD1"     "PDCD11"   "PLCB3"    "PLXNC1"   "POLE"     "PRPF38B"  "PSMG1"    "SEZ6L2"   "SIGLEC9"  "SLC16A7"  "SLC35A1"  "SLC4A7"   "SLC7A6"
[40] "SLX4"     "STAB2"    "TEX10"    "TLR2"     "TRMT61A"  "TYW3"     "UTP15"    "ZC3H7B"


 dnt.s2 %>% select(hgnc) %>% unique %>% unlist %>% unname %>% sort %>% intersect(msk.96)

 [1] "ABCC1"    "ABCC4"    "ARID2"    "C16orf54" "CCDC88A"  "CD84"     "DOCK11"   "ENG"      "EPB41"    "ITGA4"    "LAIR1"    "MMP14"    "NRD1"
[14] "PLCB3"    "PLXNC1"   "SIGLEC9"  "SLC16A7"  "SLC4A7"   "SLC7A6"   "TEX10"    "TLR2"



step.2 %>%
filter(hgnc %in% dnt.picks) %>%
arrange(gross.mean.abundance) %>%
select(gene = hgnc, tissue, level = tissue.max) %>%
PlotTissue(pdf = TRUE, width = 10, height = 5, file.name = 'dnt_mutant_intersect_proteomics_max.pdf')

step.2 %>%
filter(hgnc %in% dnt.picks) %>%
arrange(gross.mean.abundance) %>%
select(gene = hgnc, tissue, level = tissue.mean) %>%
PlotTissue(pdf = TRUE, width = 10, height = 5, file.name = 'dnt_mutant_intersect_proteomics_mean.pdf')


CD14, CDK6, CKAP4, FSCN1, HBA1, HBA2, HBB, PRSS57











dnt.rna <-
	read.delim('../tables/dnt_rna_seq.txt', stringsAsFactors = FALSE) %>%
	tbl_df %>%
	rowwise %>%
	mutate(DNMT3a_MUT = mean(s_DNMT3a_MUT_10_23, s_DNMT3a_mutant)) %>%
	mutate(MIGR1_control = mean(s_MIGR1, s_MIGR1_11_6)) %>%
	ungroup %>%
	select(hugo, dnt.wt = s_DNMT3a_WT, mll = s_MLLAF9, dnt.mut = DNMT3a_MUT, mig.control = MIGR1_control)


> nature.all %>% intersect(dnt.rna.cgenes)
 [1] "ADAM19"   "ADGRG1"   "AHSP"     "AKR1C3"   "ALAS2"    "ANGPT1"   "AQP9"     "BCL6"     "CCDC109B" "CD14"     "CDK6"     "CECR1"    "CHST15"
[14] "CKAP4"    "COL24A1"  "CTSH"     "CTSL"     "CXCL16"   "CXCR4"    "DNMT3B"   "E2F2"     "EMP1"     "EPDR1"    "FAM212A"  "FAM30A"   "FCGR2A"
[27] "FCN1"     "FCRL3"    "FCRLA"    "GIMAP4"   "GNLY"     "GPR56"    "GUCY1A3"  "GZMA"     "HBA1"     "HBA2"     "HBB"      "HBM"      "IL2RB"
[40] "JAZF1"    "KCNK17"   "KIAA0125" "KLRB1"    "LILRA5"   "LRP5"     "MCUB"     "MMRN1"    "MTSS1"    "NPL"      "NYNRIN"   "RBM38"    "SLC15A3"
[53] "SLC25A37" "SLC4A1"   "SLC7A7"   "TFPI"     "VWF"








dnt.rna.fold <-
	dnt.rna %>% mutate(dnt.mut.mig = dnt.mut/mig.control) %>% filter(dnt.mut.mig > 1.5) %>% arrange(desc(dnt.mut.mig))

dnt.rna.genes <-
	dnt.rna.fold %$% hugo %>% unique

dnt.rna.cgenes <-
	dnt.rna.genes %>% GetEnsemblID %$% ensembl.gene %>% GetHugoID %$% hgnc %>% c(dnt.rna.genes) %>% unique %>% sort


nature.all %>% intersect(dnt.rna.cgenes)


# > 1x
"ADAM19"   "ADGRG1"   "AHSP"     "AKR1C3"   "ALAS2"    "ANGPT1"   "AQP9"     "BCL6"     "CCDC109B" "CD14"     "CDK6"     "CECR1"    "CHST15"
"CKAP4"    "COL24A1"  "CTSH"     "CTSL"     "CXCL16"   "CXCR4"    "DNMT3B"   "E2F2"     "EMP1"     "EPDR1"    "FAM212A"  "FAM30A"   "FCGR2A"
"FCN1"     "FCRL3"    "FCRLA"    "GIMAP4"   "GNLY"     "GPR56"    "GUCY1A3"  "GZMA"     "HBA1"     "HBA2"     "HBB"      "HBM"      "IL2RB"
"JAZF1"    "KCNK17"   "KIAA0125" "KLRB1"    "LILRA5"   "LRP5"     "MCUB"     "MMRN1"    "MTSS1"    "NPL"      "NYNRIN"   "RBM38"    "SLC15A3"
"SLC25A37" "SLC4A1"   "SLC7A7"   "TFPI"     "VWF"

# > 1.5x
"ADAM19"   "ADGRG1"   "AHSP"     "ALAS2"    "ANGPT1"   "AQP9"     "BCL6"     "CD14"     "CDK6"     "CECR1"    "CHST15"   "CTSL"     "CXCL16"
"EMP1"     "FAM30A"   "FCN1"     "FCRL3"    "FCRLA"    "GIMAP4"   "GNLY"     "GPR56"    "GUCY1A3"  "GZMA"     "HBA1"     "HBA2"     "HBB"
"HBM"      "JAZF1"    "KCNK17"   "KIAA0125" "KLRB1"    "LILRA5"   "MMRN1"    "MTSS1"    "NPL"      "NYNRIN"   "SLC25A37" "SLC4A1"   "SLC7A7"
"VWF"

# > 2x
"ADAM19" "AHSP"   "ALAS2"  "AQP9"   "CD14"   "CECR1"  "CHST15" "FCN1"   "FCRL3"  "GIMAP4" "GNLY"   "GZMA"   "HBA1"   "HBA2"   "HBB"    "HBM"
"JAZF1"  "KLRB1"  "LILRA5" "MMRN1"  "MTSS1"  "NPL"    "SLC4A1" "VWF"



dnt.rna.agenes <-
	dnt.rna$hugo %>% GetEnsemblID %$% ensembl.gene %>% GetHugoID %$% hgnc %>% c(dnt.rna$hugo) %>% unique %>% sort



> dnt.rna.cgenes %>% intersect(msk.96)
"ABCC4"    "ANK1"     "ARID2"    "ATP11A"   "CBL"      "CCDC88A"  "CCR1"     "CD209"    "CD84"     "CD96"     "DOCK10"   "DOCK11"   "DTNA"
"ENG"      "EPB41"    "FCAR"     "GYPA"     "ITGA4"    "ITGB3"    "KIT"      "LILRA6"   "LILRB2"   "LILRB4"   "MTHFR"    "NOTCH2"   "PLXNC1"
"RABGAP1L" "SIGLEC9"  "SLC16A7"  "SLC2A9"   "SLC31A1"  "SLC4A7"   "SORT1"    "ST14"     "VCPIP1"   "ZZEF1"



dnt.rna.s1 <-
	step.0 %>%
	filter(db.num >= 2) %>%
	filter(membrane == TRUE) %>%
	filter(hgnc %in% dnt.rna.agenes) %>%
	arrange(gross.mean.abundance, hgnc) %>%
	select(
		hgnc, ensembl.gene, tissue, db.num, fill, hpa, hpm, pdb,
		tissue.max, tissue.mean, gross.mean.abundance, rna, everything(), -membrane)





dnt.rna.s2 <-
	dnt.rna.s1 %>%
	filter(gross.mean.abundance < 1) %>%
	group_by(hgnc) %>%
	mutate(tissue.mean.round =
		ifelse(tissue.mean <  0.5, 0,
		ifelse(tissue.mean <  1.5, 1,
		ifelse(tissue.mean <  2.5, 2,
		ifelse(tissue.mean <= 3,   3,
			NA))))) %>%
	mutate(max.test = tissue.max < 3 | tissue %in% pass.tissue) %>%
	filter(!any(na.omit(max.test) == FALSE)) %>%
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
	arrange(gross.mean.abundance) %>%
	full_join(expand.grid(hgnc = unique(.$hgnc), tissue = unique(.$tissue), stringsAsFactors = FALSE), by = c('hgnc', 'tissue')) %>%
	mutate(group.n = as.integer(factor(hgnc, levels=unique(hgnc)))) %>%
	arrange(hgnc, tissue)

dnt.rna.s2 %>% select(hgnc) %>% unique %>% unlist %>% unname %>% sort



#----------------------------------------------------------------------




mll.rna.fold <-
	dnt.rna %>% mutate(mll.mig = mll/mig.control) %>% filter(mll.mig > 1.5) %>% arrange(desc(mll.mig))


mll.rna.genes <-
	mll.rna.fold %$% hugo %>% unique

mll.rna.cgenes <-
	mll.rna.genes %>% GetEnsemblID %$% ensembl.gene %>% GetHugoID %$% hgnc %>% c(dnt.rna.genes) %>% unique %>% sort


nature.all %>% intersect(mll.rna.cgenes)

 [1] "ADAM19"   "AHSP"     "ALAS2"    "ANGPT1"   "AQP9"     "BCL6"     "CD14"     "CDK6"     "CECR1"    "CHST15"   "CKAP4"    "CTSL"     "CXCL16"
[14] "EMP1"     "FCN1"     "FCRL3"    "FCRLA"    "GIMAP4"   "GNLY"     "GPR56"    "GUCY1A3"  "GZMA"     "GZMB"     "HBA1"     "HBA2"     "HBB"
[27] "HBM"      "ITPR3"    "JAZF1"    "KCNK17"   "KIAA0125" "KLRB1"    "LILRA5"   "MMRN1"    "MTSS1"    "NPL"      "NYNRIN"   "SLC25A37" "SLC4A1"
[40] "SLC7A7"   "VWF"


rna.agenes <-
	dnt.rna$hugo %>% GetEnsemblID %$% ensembl.gene %>% GetHugoID %$% hgnc %>% c(dnt.rna$hugo) %>% unique %>% sort



rna.cgenes %>% intersect(msk.96)



mll.rna.s1 <-
	step.0 %>%
	filter(db.num >= 2) %>%
	filter(membrane == TRUE) %>%
	filter(hgnc %in% mll.rna.cgenes) %>%
	arrange(gross.mean.abundance, hgnc) %>%
	select(
		hgnc, ensembl.gene, tissue, db.num, fill, hpa, hpm, pdb,
		tissue.max, tissue.mean, gross.mean.abundance, rna, everything(), -membrane)





mll.rna.s2 <-
	mll.rna.s1 %>%
	filter(gross.mean.abundance < 1) %>%
	group_by(hgnc) %>%
	mutate(tissue.mean.round =
		ifelse(tissue.mean <  0.5, 0,
		ifelse(tissue.mean <  1.5, 1,
		ifelse(tissue.mean <  2.5, 2,
		ifelse(tissue.mean <= 3,   3,
			NA))))) %>%
	mutate(max.test = tissue.max < 3 | tissue %in% pass.tissue) %>%
	filter(!any(na.omit(max.test) == FALSE)) %>%
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
	arrange(gross.mean.abundance) %>%
	full_join(expand.grid(hgnc = unique(.$hgnc), tissue = unique(.$tissue), stringsAsFactors = FALSE), by = c('hgnc', 'tissue')) %>%
	mutate(group.n = as.integer(factor(hgnc, levels=unique(hgnc)))) %>%
	arrange(hgnc, tissue)

dnt.rna.s2 %>% select(hgnc) %>% unique %>% unlist %>% unname %>% sort



> mll.rna.s2 %>% select(hgnc) %>% unique %>% unlist %>% unname %>% intersect(msk.96)
 [1] "ABCC4"    "ANK1"     "ARID2"    "ATP11A"   "CBL"      "CCDC88A"  "CCR1"     "CD209"    "CD84"     "CD96"     "DOCK10"   "DOCK11"   "DTNA"
[14] "ENG"      "EPB41"    "FAM234A"  "FCAR"     "GUCY2D"   "GYPA"     "ICAM5"    "ITGA4"    "ITGAX"    "ITGB3"    "KIT"      "LILRA6"   "LILRB2"
[27] "LILRB4"   "MFGE8"    "MTHFR"    "NOTCH2"   "OSBPL3"   "PIEZO1"   "PLXNC1"   "RABGAP1L" "SIGLEC9"  "SLC16A7"  "SLC2A3"   "SLC2A9"   "SLC31A1"
[40] "SLC4A7"   "SORT1"    "ST14"     "VCPIP1"   "ZZEF1"


blood.aml.t %>% intersect(msk.96)

"ADGRE2" "CCR1"   "LTB4R"  "P2RY13"


step.2 %>%
filter(hgnc %in% c('ADGRE2', 'CCR1', 'LTB4R', 'P2RY13')) %>%
arrange(gross.mean.abundance) %>%
select(gene = hgnc, tissue, level = tissue.mean) %>%
PlotTissue(pdf = TRUE, width = 10, height = 1.95, file.name = 'blood_aml_picks_mean.pdf')

step.2 %>%
filter(hgnc %in% c('ADGRE2', 'CCR1', 'LTB4R', 'P2RY13')) %>%
arrange(gross.mean.abundance) %>%
select(gene = hgnc, tissue, level = tissue.max) %>%
PlotTissue(pdf = TRUE, width = 10, height = 1.95, file.name = 'blood_aml_picks_max.pdf')


step.0 %>%
filter(hgnc %in% (blood.aml.t %>% list.filter(!. %in% c('ADGRE2', 'CCR1', 'LTB4R', 'P2RY13', 'GPR141', 'ADORA2A')))) %>%
arrange(gross.mean.abundance) %>%
select(gene = hgnc, tissue, level = tissue.mean) %>%
PlotTissue(pdf = TRUE, width = 10, height = 8.8, file.name = 'blood_aml_nonpicks_mean.pdf')

step.0 %>%
filter(hgnc %in% (blood.aml.t %>% list.filter(!. %in% c('ADGRE2', 'CCR1', 'LTB4R', 'P2RY13', 'GPR141', 'ADORA2A')))) %>%
arrange(gross.mean.abundance) %>%
select(gene = hgnc, tissue, level = tissue.max) %>%
PlotTissue(pdf = TRUE, width = 10, height = 8.8, file.name = 'blood_aml_nonpicks_max.pdf')






#------------------------------




blood.aml <-
	read.delim('../tables/blood_aml_genes.tsv', stringsAsFactors = FALSE, header = FALSE) %>%
	tbl_df %>%
	unlist %>%
	unname %>%
	sort

blood.aml.t <-
	blood.aml %>% GetEnsemblID %$% ensembl.gene %>% GetHugoID %$% hgnc %>% c(blood.aml) %>% unique %>% sort
