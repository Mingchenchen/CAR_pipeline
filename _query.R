#-------------#
#             #
#  retreival  #
#             #
#-------------#


#--------
# imports
#--------


library(broom)
library(matrixStats)
library(scales)

#load('clean.RData')
source('_convenience.R')
source('_keys.R')


#----------
# bloodspot
#----------


blood.pooled <-
	list.files('../data/bloodspot', full.names = TRUE, pattern = '*csv') %>%
	map( ~ {
		base <-
			read.delim(.x, sep = ',', stringsAsFactors = FALSE, header = FALSE) %>%
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

		blood.master <-
			base %>%
			select(1) %>%
			set_names('cell') %>%
			mutate(cell = toupper(cell)) %>%
			mutate(hgnc = hgnc) %>%
			bind_cols(vals) %>%
			mutate(exp.max = exp.max) %>%
			arrange(hgnc, cell) %>%
			mutate(cell = ifelse(is.na(cell), 'nan', cell))

		blood.master %>%
		select(cell, hgnc, exp.max) %>%
		unique %>%
		as.data.frame
	}) %>%
	data.table::rbindlist(., fill = TRUE) %>%
	unique %>%
	tbl_df %>%
	filter(!cell %in% exclude.lines) %>%
	mutate(cell = ifelse(substr(strip(cell), 1, 3) == '7.0', 7, cell)) %>%
	mutate(cell = ifelse(substr(strip(cell), 1, 3) == '8.0', 8, cell)) %>%
	group_by(cell, hgnc) %>%
	mutate(max.exp = max(exp.max)) %>%
	ungroup %>%
	select(cell, hgnc, max.exp) %>%
	unique %>%
	group_by(hgnc) %>%
	mutate(gene.order = mean(max.exp)) %>%
	ungroup %>%
	arrange(desc(gene.order)) %>%
	select(-gene.order) %>%
	mutate(hgnc = case_when(
		.$hgnc %in% 'EMR2'     ~ 'ADGRE2',
		.$hgnc %in% 'ITFG3'    ~ 'FAM234A',
		.$hgnc %in% 'C16ORF54' ~ 'C16orf54',
		.$hgnc %in% 'C16ORF88' ~ 'KNOP1',
		.$hgnc %in% 'TNFRSF1B'  ~ 'MIR7846',
		TRUE ~ hgnc))


#-----------
# rna ratios
#-----------


blood.pro.ratio <-
	full_join(
		blood.pooled %>% filter(!cell %in% progeniter.lines) %>% rename(cell.exp = max.exp),
		blood.pooled %>% filter(cell %in% progeniter.lines) %>%
			group_by(hgnc, cell) %>% mutate(pro.exp = max(max.exp, na.rm = TRUE)) %>% ungroup %>%
			group_by(hgnc) %>% mutate(pro.exp = mean(max.exp, na.rm = TRUE)) %>% ungroup %>%
			select(-max.exp, -cell) %>% unique,
		by = 'hgnc') %>%
	rowwise %>%
	mutate(exp.pro.ratio     = cell.exp / pro.exp) %>%
	mutate(exp.pro.ratio.log = log10(exp.pro.ratio)) %>%
	ungroup


#------------
# ratio curve
#------------


blood.pro.ratio.curve <-
	blood.pro.ratio %>%
	select(hgnc, exp.pro.ratio.log) %>%
	unique %>%
	select(exp = exp.pro.ratio.log) %>%
	PlotExpDist(binwidth = 0.02, sigma.scalar = 2)

SavePlot(gg = blood.pro.ratio.curve$gg.dens, file = 'blood_pro_ratio_dens.png', type = 'png')
SavePlot(gg = blood.pro.ratio.curve$gg.quant, file = 'blood_pro_ratio_quant.png', type = 'png')

blood.pro.ratio.mu    <- blood.pro.ratio.curve$ml %>% filter(term == 'mu') %$% estimate
blood.pro.ratio.sigma <- blood.pro.ratio.curve$ml %>% filter(term == 'sigma') %$% estimate

10^(blood.pro.ratio.mu + (blood.pro.ratio.sigma))
10^(blood.pro.ratio.mu + (blood.pro.ratio.sigma * 2))


blood.ensembl.id.dictionary <-
	GetID(unique(blood.pooled$hgnc), 'symbol') %>%
	bind_rows(ensembl.id.dictionary) %>%
	rename(hgnc = query) %>%
	unique %>%
	group_by(hgnc) %>%
	filter(!is.na(ensembl.gene) & n() != 1) %>%
	ungroup


blood.exp.cut <-
	blood.pro.ratio %>%
	left_join(blood.ensembl.id.dictionary, by = 'hgnc') %>%
	select(cell, hgnc, ensembl.gene, everything()) %>%
	filter(exp.pro.ratio.log >= (blood.pro.ratio.mu + (blood.pro.ratio.sigma * 2)))


#-----------------
# AML Surface Pool
#-----------------


aml.surface.proteins <-
	#msk.dnmt3a.triple.surface %>%
	c(msk.aml.single.surface, msk.aml.triple.surface, jpro.aml.triple.surface, aml.additions) %>%
	unique %>%
	sort


#-------------------
# Surface Proteomics
#-------------------


#compare.protein <- compare.protein.dnmt3a

step.0 <-
	master %>%
	mutate(aml.surface.test = ensembl.gene %in% aml.surface.proteins) %>%
	filter(aml.surface.test == TRUE | hgnc %in% compare.protein)


#----
# RNA
#----


step.1 <-
	step.0 %>%
	left_join(select(blood.exp.cut, -hgnc), by = 'ensembl.gene') %>%
	mutate(blood.test = hgnc %in% blood.exp.cut$hgnc) %>%
	filter(blood.test == TRUE | hgnc %in% compare.protein) %>%
	unique %>%
	arrange(prot.mean, hgnc)


#---
# QC
#---


step.2 <-
	step.1 %>%
	mutate(db.test = db.num >= 2) %>%
	mutate(exclude = hgnc %in% exclude.protein) %>%
	filter((db.test  == TRUE & membrane == TRUE & exclude  == FALSE) | hgnc %in% compare.protein) %>%
	arrange(prot.mean, hgnc)


#-----------------------------------------
# discard medium overall & high individual 
#-----------------------------------------


step.3 <-
	step.2 %>%
	mutate(two.test = prot.mean <= 2) %>%
	group_by(ensembl.gene) %>%
	mutate(max.test = tissue.max < 3 | tissue %in% pass.tissue) %>%
	mutate(max.test = all(na.omit(max.test) == TRUE)) %>%
	ungroup %>%
	#filter((two.test == TRUE & max.test == TRUE) | hgnc %in% compare.protein) %>%
	filter((two.test == TRUE & max.test == TRUE) | hgnc %in% c('IL3RA', 'CD33', 'CLEC12A')) %>%
	mutate(all.test =
	aml.surface.test == TRUE &
	blood.test == TRUE &
	db.test == TRUE &
	membrane == TRUE &
	exclude  == FALSE &
	two.test == TRUE &
	max.test  == TRUE) %>%
	mutate(tissue.mean = if_else(hgnc == 'ADGRE2' & tissue == 'skeletal muscle', 0, tissue.mean)) %>%
	mutate(tissue.max = if_else(hgnc == 'ADGRE2' & tissue == 'skeletal muscle', 0, tissue.max)) %>%
	filter(hgnc != 'P2RY13') %>%
	mutate(hgnc = if_else(hgnc == 'MIR7846', 'TNFRSF1B', hgnc))


#---------
# plotting
#---------


step.3 %>%
filter(all.test == TRUE) %>%
select(gene = hgnc, level = tissue.mean, tissue) %>%
unique %>%
GeneRename %>%
PlotTissue(pdf = TRUE, file.name = 'step_3_mean_protein_spleen.pdf', width = 10.15, height = 6, order = TRUE) 

step.3 %>%
filter(all.test == TRUE) %>%
select(gene = hgnc, level = tissue.max, tissue) %>%
unique %>%
GeneRename %>%
PlotTissue(pdf = TRUE, file.name = 'step_3_max_protein_spleen.pdf', width = 10, height = 6, order = TRUE) 


step.3 %>%
select(gene = hgnc, level = tissue.mean, tissue, all.test) %>%
mutate(split = ifelse(all.test == TRUE, '', ' ')) %>%
unique %>%
GeneRename %>%
PlotTissue(pdf = TRUE, faceting = TRUE, file.name = 'step_3_mean_protein_compare.pdf', width = 10.15, height = 6.6, order = TRUE)
#PlotTissue(pdf = TRUE, faceting = TRUE, file.name = 'step_3_DNMT3a_mean_protein_compare.pdf', width = 10, height = 2.6, order = TRUE) 


step.3 %>%
select(gene = hgnc, level = tissue.max, tissue, all.test) %>%
mutate(split = ifelse(all.test == TRUE, '', ' ')) %>%
unique %>%
GeneRename %>%
PlotTissue(pdf = TRUE, faceting = TRUE, file.name = 'step_3_max_protein_compare.pdf', width = 10, height = 6.8, order = TRUE) 
#PlotTissue(pdf = TRUE, faceting = TRUE, file.name = 'step_3_DNMT3a_max_protein_compare.pdf', width = 10, height = 2.6, order = TRUE) 


step.3 %>%
filter(!is.na(ensembl.gene)) %>%
select(gene = hgnc, level = tissue.max, tissue, all.test) %>%
filter(gene == 'IL3RA') %>%
unique %>%
GeneRename %>%
PlotTissue(pdf = TRUE, faceting = FALSE, file.name = 'IL3RA_protein_compare.pdf', width = 10, height = 1.2, order = TRUE) 




master %>%
filter(hgnc %in% c('ERBB2', 'CEACAM5', 'CA9')) %>%
select(gene = hgnc, level = tissue.max, tissue) %>%
unique %>%
GeneRename %>%
PlotTissue(pdf = TRUE, faceting = FALSE, file.name = 'trial_targets.pdf', width = 10, height = 1.55, order = TRUE)


master %>%
	filter(hgnc == 'CD19') %>%
	rename(gene = hgnc, level = tissue.max) %>%
	PlotTissue(pdf = TRUE, faceting = FALSE, file.name = 'CD19.pdf', width = 10, height = 1.2, order = TRUE)



#------
# stats
#------

# msk surface proteomics
c(msk.aml.single$id, msk.aml.triple$id) %>% unique %>% sort %>% length

surface.proteomics <-
	master %>%
	filter(ensembl.gene %in% c(msk.aml.single.surface, msk.aml.triple.surface)) %$%
	hgnc %>% unique %>% sort

reported.markers <-
	master %>%
	filter(ensembl.gene %in% c(jpro.aml.triple.surface, aml.additions)) %$%
	hgnc %>% unique %>% sort



# previously reported molecules
jpro.aml.triple$id %>% unique %>% length

# master
master %$% hgnc %>% unique %>% sort %>% length

#step 0
step.0 %>% filter(aml.surface.test == TRUE) %$% hgnc %>% unique %>% sort %>% length

# step 1
step.1 %>% filter(
	aml.surface.test == TRUE &
	blood.test == TRUE) %$% hgnc %>% unique %>% sort %>% length

# step 2
step.2 %>%
filter(
	aml.surface.test == TRUE &
	blood.test == TRUE &
	db.test == TRUE &
	membrane == TRUE &
	exclude  == FALSE) %$% hgnc %>% unique %>% sort %>% length

# step 3
step.3 %>%
filter(
	aml.surface.test == TRUE &
	blood.test == TRUE &
	db.test == TRUE &
	membrane == TRUE &
	exclude  == FALSE &
	two.test == TRUE &
	max.test  == TRUE) %$% hgnc %>% unique %>% sort %>% length


GetStats(list(master = master, `step 0` = step.0, `step 1` = step.1, `step 2` = step.2, `step 3` = step.3 %>% filter(!hgnc %in% compare.protein) ))

select(step.2, hgnc, ensembl.gene) %>%
unique %>%
filter(!hgnc %in% compare.protein) %>%
arrange(hgnc) %>%
as.data.frame

select(step.3, hgnc, ensembl.gene) %>%
unique %>%
filter(!hgnc %in% compare.protein) %>%
arrange(hgnc) %>%
as.data.frame


#-----------------------
# cell ratio comparisons
#-----------------------

cell.gg <-
	ggplot(cell.table, aes(cell, gname)) + 
	geom_tile(aes(fill = compare), colour = 'grey') +
	theme(
		legend.title     = element_blank(),
		axis.title.x     = element_blank(),
		axis.title.y     = element_blank(),
		text             = element_text(size = 10),
		axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1),
		panel.background = element_rect(fill = 'grey'),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.border     = element_rect(colour = 'black', fill = NA, size = 1)) +
		facet_grid(. ~ type, scales = 'free_y', space = 'free_y')

pdf('cell_table.pdf', 24, 10)
	plot(cell.gg)
dev.off()


#---------------------
# surface protein venn
#---------------------

require(venneuler)

Venn <- function(...) {

	args <- c(as.list(environment()), list(...))

	a <- unique(args[[2]])
	b <- unique(args[[1]])

	ab <- intersect(a, b)

	v <- c(A = length(a) - length(ab), B = length(b) - length(ab), 'A&B' = length(ab)) %>% venneuler

	v$labels <- c(names(args)[[2]], names(args)[[1]])

	v %>% plot

	text(0.2, 0.44, length(b))
	text(0.45, 0.44, length(ab))
	text(0.8, 0.44, length(a))
}

pdf('venn.pdf', width = 8, height = 8)
	Venn(`on-site surfaceomics` = surface.proteomics, `AML surfaceome\nliterature` = reported.markers)
dev.off()


#--------
# objects
#--------

objects <- c(
'pass.tissue',
'palette',
'aml.surface.proteins',
'exclude.lines',
'progeniter.lines',
'blood.pooled',
'blood.pro.ratio',
'blood.pro.ratio.curve',
'blood.pro.ratio.mu',
'blood.pro.ratio.sigma',
'blood.pro.ratio.curve2',
'blood.pro.ratio.mu',
'blood.pro.ratio.sigma',
'blood.ensembl.id.dictionary',
'blood.exp.cut')
