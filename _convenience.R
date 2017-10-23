
GetStats <- function(exp.table) {
	if(class(exp.table)[[1]] == 'list') {
		exp.table %>% map( ~ { 
			data.frame(
				ensembl = .x$ensembl.gene %>% unique %>% na.omit %>% length,
				hgnc    = .x$hgnc         %>% unique %>% na.omit %>% length
			)
		}) %>%
		bind_rows %>%
		bind_cols(data.frame(step = names(exp.table))) %>%
		select(step, ensembl, hgnc) %>%
		print(row.names = FALSE)
	} else {
		data.frame(
			ensembl = exp.table$ensembl.gene %>% unique %>% na.omit %>% length,
			hgnc    = exp.table$hgnc         %>% unique %>% na.omit %>% length
		) %>%
		print(row.names = FALSE)
	}
}


RowToNames <- function(df, row.num = 1) {

	col.names <- df[row.num,]

	df %<>% set_names(col.names)

	df[-row.num,]
}


strip <- function(string) { gsub('^\\s+|\\s+$', '', string) }


GeneRename <- function(events) {
	events %>%
	mutate(gene =
		ifelse(gene == 'MIR7846', 'TNFRSF1B',
		ifelse(gene == 'CLEC12A', 'CLEC12A (CLL1)',
		ifelse(gene == 'IL3RA', 'IL3RA (CD123)',
		ifelse(gene == 'FUT3', 'FUT3 (LeY)',
		ifelse(gene == 'HAVCR2', 'HAVCR2 (TIM3)',
		ifelse(gene == 'IL2RA', 'IL2RA (CD25)',
		ifelse(gene == 'FCGR2A', 'FCGR2A (CD32)',
		ifelse(gene == 'CA9', 'CA9 (CAIX)',
		ifelse(gene == 'ERBB2', 'ERBB2 (HER-2)',
		gene))))))))))
}
