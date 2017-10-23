
blood.pooled.sample <-
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
			arrange(hgnc, cell) %>%
			mutate(cell = ifelse(is.na(cell), 'nan', cell)) %>%
			gather(patient, exp, -cell, -hgnc) %>%
			as_data_frame
	})

blood.pooled.sample.mean <-
	blood.pooled.sample %>%
	map( ~ {
		.x %>%
		group_by(cell, hgnc, patient) %>%
		mutate(exp = mean(exp)) %>%
		unique
	}) %>%
	data.table::rbindlist(., fill = TRUE) %>%
	unique %>%
	tbl_df %>%
	filter(!cell %in% exclude.lines) %>%
	mutate(cell = ifelse(substr(strip(cell), 1, 3) == '7.0', 7, cell)) %>%
	mutate(cell = ifelse(substr(strip(cell), 1, 3) == '8.0', 8, cell)) %>%
	mutate(hgnc = case_when(
		.$hgnc %in% 'EMR2'     ~ 'ADGRE2',
		.$hgnc %in% 'ITFG3'    ~ 'FAM234A',
		.$hgnc %in% 'C16ORF54' ~ 'C16orf54',
		.$hgnc %in% 'C16ORF88' ~ 'KNOP1',
		TRUE ~ hgnc))

blood.pooled.sample.hsc.mean <-
	blood.pooled.sample.mean %>%
	filter(cell == 'HSC') %>%
	rename(hsc.exp = exp) %>%
	select(-cell)

blood.pooled.sample.mean <-
	blood.pooled.sample.mean %>%
	filter(cell != 'HSC') %>%
	full_join(blood.pooled.sample.hsc.mean) %>%
	filter(!is.na(hsc.exp))

blood.pooled.ratio <-
	blood.pooled.sample.mean %>%
	rowwise %>%
	mutate(ratio = exp / hsc.exp) %>%
	ungroup


mutate(gain.amp.fisher = fisher.test( data.frame( a = c(gain.lt.a, gain.gt.eq.a),
												  b = c(gain.lt.b, gain.gt.eq.b) ) )$p.value ) %>%
# p adjustments
mutate(amp.pos.fisher.adj  = p.adjust(amp.pos.fisher,  'BH')) %>%
# 0.05 cutoffs
mutate(amp.pos.fisher.adj.cut  = ifelse(amp.pos.fisher.adj  < 0.05, amp.pos.fisher.adj,  NA)) %>%

# log
mutate(amp.pos.fisher.adj.cut.log  = log10(amp.pos.fisher.adj.cut)) %>%


library(Hmisc)

cor(mtcars, use="pairwise.complete.obs", method="pearson") 

blood.pooled.ratio %>%
split(cell) %>%
map( ~ {
	x = .x$
	rcorr(x=1:10, y=11:20, type="pearson")
})


#-------------------------------------------------

list.files('~/Desktop/diff_exp')

library(limma)
targets <- readTargets(path = '~/Desktop/diff_exp')

# 17.3.3 The expression profiles

x <- read.ilmn(files="~/Desktop/diff_exp/probe profile.txt", ctrlfiles="~/Desktop/diff_exp/control probe profile.txt", other.columns="Detection")

table(x$genes$Status)


# 17.3.4 How many probes are truly expressed?

pe <- propexpr(x)
dim(pe) <- c(4,3)
dimnames(pe) <- list(CellType=c("MS","Stroma","ML","LP"),Donor=c(1,2,3))


# 17.3.5 Normalization and filtering

y <- neqc(x)
expressed <- rowSums(y$other$Detection < 0.05) >= 3

y <- y[expressed,]
plotMDS(y, labels=targets$CellType)


# 17.3.6 Within-patient correlations

ct <- factor(targets$Type)

design <- model.matrix(~0+ct)

colnames(design) <- levels(ct)

dupcor <- duplicateCorrelation(y, design, block = targets$Donor)

dupcor$consensus.correlation


# 17.3.7 Differential expression between cell types

fit <- lmFit(y, design, block = targets$Donor, correlation = dupcor$consensus.correlation)
contrasts <- makeContrasts(mL-MS, pL-MS, mL-pL, levels = design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend = TRUE)
summary(decideTests(fit2, method = "global"))

topTable(fit2, coef=1)


# 17.3.8 Signature genes for luminal progenitor cells

ct <- relevel(ct, ref = "pL")
design <- model.matrix(~ct)
fit <- lmFit(y, design, block = targets$Donor, correlation = dupcor$consensus.correlation)
fit2 <- fit[,c("ctMS", "ctmL")]
fit2 <- eBayes(fit2, trend = TRUE)

results <- decideTests(fit2, lfc = 1)
vennDiagram(results, include=c("up", "down"))

