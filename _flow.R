RowToNames <- function(x) {
	names <- unlist(x[1,])
	nonames <- x[-1,]
	set_names(as.data.frame(nonames, stringsAsFactors = FALSE), c(names))
}

combi <-
	read.xlsx('../data/flow_combi_sheets.xlsx', sheet = 1, colNames = FALSE) %>%
	t %>%
	RowToNames %>%
	tbl_df %>%
	mutate(patient = row_number()) %>%
	select(patient, everything()) %>%
	rename(CLEC12A_int_ADGRE2 = `CLEC12A & ADGRE2`, CLEC12A_uni_ADGRE2 = `CLEC12A + ADGRE2`)



#-------------------------------------------



flow <-
	read.xlsx('../data/patient_flow.xlsx') %>%
	tbl_df %>%
	gather(patient, exp.pct, A:O) %>%
	arrange(hgnc, patient) %>%
	filter(!is.na(exp.pct)) %>%
	mutate(cut = exp.pct > 50) %>%
	group_by(hgnc) %>%
	mutate(pct.patient = sum(cut) / n())

flow.pct.patient <-
	flow %>%
	select(hgnc, pct.patient) %>%
	unique %>%
	rename(level = pct.patient)





LevelBarPlot <- function(table, title = '% patients with greater 50% expression', width = 9, height = 5) {

	gg <-
		ggplot(table, aes(x = hgnc, y = level)) +
		ggtitle(title) +
		geom_bar(stat = "identity", fill = 'darkblue', colour = 'darkblue') +
		xlab("gene") +
		ylab(title)

	pdf('patients_with_greater_50_expression.pdf', width = width, height = height)
		plot(gg)
	dev.off()
}


flow.pct.patient %>%
LevelBarPlot



#-----------



LevelViolinPlot <- function(table, title = 'expression per patient', width = 9, height = 5) {

	gg <-
		ggplot(table, aes(hgnc, level)) +
		ggtitle(title) +
		geom_violin(fill = 'darkblue', colour = 'darkblue') +
		xlab("gene") +
		ylab(title)


	pdf(paste(c(title, 'pdf'), collapse = '.'), width = width, height = height)
		plot(gg)
	dev.off()
}



flow %>%
select(hgnc, level = exp.pct) %>%
LevelViolinPlot









flow <-
	read.xlsx('../data/patient_flow.xlsx') %>%
	tbl_df %>%
	gather(patient, exp.pct, A:O) %>%
	arrange(hgnc, patient) %>%
	filter(!is.na(exp.pct)) %>%
	mutate(cut.low = exp.pct >= 0 & exp.pct < (100/3)) %>%
	mutate(cut.med = exp.pct >= (100/3) & exp.pct < (200/3)) %>%
	mutate(cut.hi = exp.pct  >= (200/3) & exp.pct <= 100) %>%
	group_by(hgnc) %>%
	mutate(`% patient low` = sum(cut.low) / n()) %>%
	mutate(`% patient medium` = sum(cut.med) / n()) %>%
	mutate(`% patient high` = sum(cut.hi)  / n())

flow.pct.patient.stack <-
	flow %>%
	gather(bin, patient.pct, `% patient low`:`% patient high`) %>%
	select(hgnc, bin, patient.pct) %>%
	unique %>%
	mutate(bin = factor(bin, levels = c('% patient high', '% patient medium', '% patient low')))


LevelStackBarPlot <- function(table, title = '% patients with low, medium, high expression', width = 9, height = 5) {

	gg <-
		ggplot(table, aes(x = hgnc, y = patient.pct, fill = bin)) +
		ggtitle(title) +
		geom_bar(stat = 'identity') +
		xlab("gene") +
		ylab(title)


	pdf('patients_with_greater_low_medium_high_expression.pdf', width = width, height = height)
		plot(gg)
	dev.off()
}

flow.pct.patient.stack %>%
LevelStackBarPlot



















