


blood.hem.pooled <-
	list.files('../data/bloodspot_hematopoiesis', full.names = TRUE, pattern = '*csv') %>%
	map( ~ {
		base <-
			read.delim(.x, sep = ',', stringsAsFactors = FALSE, header = FALSE) %>%
			t %>%
			as_data_frame %>%
			RowToNames

		vals <-
			base %>%
			select(-1) %>%
			mutate_all(funs(2^as.numeric(.)))

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
		TRUE ~ hgnc))

# violin

blood.hem.pooled %>%
filter(!cell %in% c('CMP', 'GMP', 'MEP')) %>%
rename(exp = max.exp) %>%
ggplot(aes(fill = hgnc, x = hgnc, y = exp)) +
geom_violin(color = NA) +
geom_point(aes(y = exp), color = "black", size = 0.5) + 
#stat_summary(fun.data = mean_sdl, geom = 'pointrange', color = 'black') +
#geom_boxplot(width = 0.05, fill = 'white') +
#theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
scale_fill_brewer(palette = 'Dark2') +
theme_minimal() +
theme(legend.position = 'none')

# bar

getPalette = colorRampPalette(brewer.pal(9, "Set1"))

blood.hem.pooled %>%
filter(!cell %in% c('CMP', 'GMP', 'MEP')) %>%
rename(exp = max.exp) %>%
group_by(cell) %>%
mutate(order = mean(exp)) %>%
arrange(order) %>%
ungroup %>%
mutate(cell = factor(cell, levels = unique(cell))) %>%
ggplot(aes(fill = cell, x = hgnc, y = exp, width = .85)) +
geom_bar(position = 'dodge', stat = 'identity') +
scale_fill_manual(values = getPalette(13)) +
theme_minimal() +
theme(axis.text.x = element_text(angle = 90, hjust = 1))

# combi bar

getPalette = colorRampPalette(brewer.pal(9, "Set1"))

blood.hem.pooled %>%
filter(!cell %in% c('CMP', 'GMP', 'MEP')) %>%
rename(exp = max.exp) %>%
group_by(cell) %>%
mutate(order = mean(exp)) %>%
arrange(order) %>%
ungroup %>%
mutate(cell = factor(cell, levels = unique(cell))) %>%
ggplot(aes(fill = cell, x = hgnc, y = exp, width = .85)) +
geom_bar(position = 'dodge', stat = 'identity') +
scale_fill_manual(values = getPalette(13)) +
theme_minimal() +
theme(axis.text.x = element_text(angle = 90, hjust = 1))

pairs <-
	data_frame(a = c('CLEC12A', 'CLEC12A', 'CLEC12A', 'CD70', 'CD33', 'CD33'), b = c('CCR1', 'CD70', 'ADGRE2', 'CD33', 'CCR1', 'ADGRE2')) %>%
	mutate(pair = str_c(a, ' + ', b)) %>%
	gather(ordinal, hgnc, a:b)


blood.hem.pooled %>%
filter(!cell %in% c('CMP', 'GMP', 'MEP')) %>%
rename(exp = max.exp) %>%
group_by(cell) %>%
mutate(order = mean(exp)) %>%
arrange(order) %>%
ungroup %>%
right_join(pairs, by = 'hgnc') %>%
mutate(cell = factor(cell, levels = unique(cell))) %>%
ggplot(aes(fill = hgnc, x = cell, y = exp, width = .85)) +
geom_bar(stat = 'identity') +
scale_fill_brewer(palette = 'Set1') +
theme_minimal() +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
facet_wrap( ~ pair)



