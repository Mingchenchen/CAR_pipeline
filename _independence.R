


#-----------
# flow combi
#-----------


library(copula)
library(psych)

library(readxl)
library(ggExtra)



dependogram(all.flow[[1]][[2]])


cor(all.flow[[2]][[1]], method = 'spearman')
pairs.panels(all.flow[[2]][[1]])

# use CD70



#--------------------------------------------------






solo.file <- '../data/flow_solo.xlsx'

ExtractSingleFlow <- function(flow.file, sheet) {

	tbl <- read.xlsx(flow.file, sheet = sheet, colNames = FALSE) %>%
	tbl_df

	tbl %>%
	set_names(c('gene', str_c('sample_', 2:ncol(tbl) - 1)))
}


lscs <- ExtractSingleFlow(solo.file, 2)
activated.t.cells <- ExtractSingleFlow(solo.file, 4)


normal.hscs <-
	ExtractSingleFlow(solo.file, 3) %>%
	mutate(hscs.mean = select(., -gene) %>% rowMeans(na.rm = TRUE)) %>%
	select(gene, hscs.mean)


aml.bulk <-
	ExtractSingleFlow(solo.file, 1)

bulk.hscs.ratio <-
	aml.bulk %>%
	left_join(normal.hscs, by = 'gene') %>%
	rowwise %>%
	mutate_at(vars(starts_with('sample_')), funs(. / hscs.mean)) %>%
	ungroup %>%
	select(-hscs.mean)


bulk.hscs.ratio %>%
gather(sample, hscs.ratio, sample_1:sample_30) %>%
ggplot(aes(fill = sample, y = hscs.ratio, x = gene)) + 
geom_bar(position = 'dodge', stat = 'identity')



#------------------------------------



library(readxl)
library(ggExtra)

combi.file <- '../data/flow_combi.xlsx'
sheet.names <- excel_sheets(combi.file)


ExtractCombiFlow <- function(flow.file, sheet) {

	tbl <- read.xlsx(flow.file, sheet = sheet, colNames = FALSE) %>%
	tbl_df

	tbl %>%
	set_names(c('run', str_c('sample_', 2:ncol(tbl) - 1)))
}


combi.sheets <-
	sheet.names %>%
	map(., ~ {
		ExtractCombiFlow(combi.file, .x)
	}) %>%
	set_names(sheet.names)


combi.sheets[[1]] %>%
# mutate(run.mean = rowMeans(select(., -run), na.rm = TRUE)) %>%
rowwise %>%
mutate_at(vars(starts_with('sample_')), funs(. > 95)) %>%


Chi <- function(combis) {
	combis %>%
	#combi.sheets[[1]] %>%
	filter(!grepl(' ', run)) %>%
	gather(sample, fraction, starts_with('sample_')) %>%
	spread(run, fraction) %>%
	filter(!is.na(rowSums(select(., -sample)))) %>%
	as.data.frame %>%
	set_rownames(.$sample) %>%
	select(-sample) %>%
	chisq.test
}

combi.sheets %>%
map( ~ { Chi(.x) })


TwoWay <- function(name, combis) {

	data <-
		combis %>%
		#combi.sheets[[4]] %>%
		filter(!grepl(' ', run)) %>%
		gather(sample, fraction, starts_with('sample_')) %>%
		spread(run, fraction) %>%
		filter(!is.na(rowSums(select(., -sample))))

	p <-
		data %>%
		ggplot(aes_string(colnames(.)[2], colnames(.)[3])) +
		geom_point(alpha = 0.6) +
		xlim(0, 100) +
		ylim(0, 100) +
		geom_rug(col = rgb(0, 0, 0, alpha = 0.8)) +
		theme_minimal()

	pdf(str_c('independence/', name, '.pdf'), width = 8, height = 8)
		plot(p)
	dev.off()
}



map2(names(combi.sheets), combi.sheets, ~ { TwoWay(.x, .y) })













