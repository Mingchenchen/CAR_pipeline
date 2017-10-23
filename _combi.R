

CombiPair <- function(combi.in.a, combi.in.b, strict = FALSE, select = FALSE, select.na = TRUE, tissue.rm = '') {

	if(strict == FALSE) {
		non.vital.cut <- 1
	} else {
		non.vital.cut <- 0
	}

	non.vital <- non.vital.combi %>% list.filter(!. %in% tissue.rm)
	vital <- vital %>% list.filter(!. %in% tissue.rm)

	ensembl <-
		expand.grid(unique(combi.in.a$ensembl.gene), unique(combi.in.b$ensembl.gene)) %>%
		t %>%
		as.data.frame(., stringsAsFactors = FALSE) %>%
		unlist %>%
		unname

	pair.picks <-
		data_frame(
			ensembl.gene = ensembl,
			pair = sort(rep(1:(length(ensembl) / 2), 2))) %>%
		group_by(pair) %>%
		filter(length(unique(ensembl.gene)) != 1) %>%
		mutate(match.pair = str_c(sort(unique(ensembl.gene)), collapse = '-')) %>%
		ungroup %>%
		group_by(match.pair) %>%
		arrange(match.pair, pair, ensembl.gene) %>%
		mutate(pair.num = row_number()) %>%
		filter(pair.num <= 2) %>%
		ungroup %>%
		left_join(unique(bind_rows(combi.in.a, combi.in.b)), by = 'ensembl.gene') %>%
		mutate(tissue.type =
			if_else(tissue %in% vital,     'vital',
			if_else(tissue %in% non.vital, 'non-vital',
										   'blood'))) %>%
		group_by(pair, tissue) %>%
		arrange(pair, tissue, ensembl.gene) %>%
		mutate(all.na = all(is.na(tissue.max))) %>%
		mutate(pass =
			if_else(tissue.type == 'vital',     min(tissue.max, na.rm = TRUE) == 0             | all.na == TRUE,
			if_else(tissue.type == 'non-vital', min(tissue.max, na.rm = TRUE) <= non.vital.cut | all.na == TRUE,
			TRUE))) %>%
		mutate(pass.na =
			if_else(tissue.type == 'vital',     na.omit(c(min(tissue.max, na.rm = TRUE) == 0,             all.na != TRUE))[[1]],
			if_else(tissue.type == 'non-vital', na.omit(c(min(tissue.max, na.rm = TRUE) <= non.vital.cut, all.na != TRUE))[[1]],
			TRUE))) %>%
		ungroup %>%
		group_by(pair) %>%
		filter(if_else(select    == TRUE, all(pass    == TRUE), TRUE)) %>%
		filter(if_else(select.na == TRUE, all(pass.na == TRUE), TRUE)) %>%
		mutate(exp.pair.mean = mean(tissue.max, na.rm = TRUE)) %>%
		arrange(exp.pair.mean, hgnc) %>%
		mutate(hgnc = factor(hgnc, levels = unique(hgnc))) %>%
		ungroup

	pair.picks
}



# 2SD
combi.picks.2 <- c('ADGRE2', 'CCR1', 'CD70', 'CD82', 'CD96', 'ITGB5', 'LILRB2', 'PTPRJ', 'TNFRSF1B', 'IL3RA', 'CLEC12A', 'CD33', 'FOLR2', 'IL2RA')

combi.picks.2 <- step.3 %$% hgnc %>% unique %>% sort

# --- non-strict 2SD ---


combi.in.a.2 <-
	step.3 %>%
	select(hgnc, ensembl.gene, tissue, tissue.max) %>%
	unique %>%
	filter(hgnc %in% combi.picks.2) %>%
	filter(hgnc == 'LILRB4')

combi.in.b.2 <-
	step.3 %>%
	select(hgnc, ensembl.gene, tissue, tissue.max) %>%
	unique %>%
	filter(hgnc %in% combi.picks.2)


pair.picks.2sd.na        <- CombiPair(combi.in.a.2, combi.in.b.2, strict = FALSE, select = FALSE,  select.na = TRUE)


pair.picks.2sd.na %>%
mutate(gene = str_c(hgnc, ensembl.gene, sep = ' : ')) %>%
rename(split = pair, level = tissue.max) %>%
PlotTissue(faceting = TRUE, pdf = TRUE, width = 9, height = 15, file.name = 'combi_2sd_na_LILRB4_spleen.pdf', order = TRUE)














