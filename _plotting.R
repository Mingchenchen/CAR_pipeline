


PlotTissue <- function(events, faceting = FALSE, pdf = FALSE, file.name = 'plot_tissue.pdf', width = 20, height = 10, order = TRUE) {

	events %<>% unique

	tissues <- c('adipose tissue', 'adrenal', 'appendix', 'bladder', 'blood', 'bone', 'brain', 'breast', 'bronchus',
				 'cerumen', 'cervix', 'epididymis', 'eye', 'fallopian tube', 'gallbladder', 'gut', 'heart', 'kidney',
				 'esophagus', 'liver', 'lung', 'lymph node', 'nasopharynx', 'oropharynx', 'ovary', 'pancreas',
				 'parathyroid', 'prostate', 'rectum', 'seminal', 'skeletal muscle', 'skin', 'smooth muscle',
				 'soft tissue', 'spinal cord', 'spleen', 'stomach', 'synovial fluid', 'testis', 'thyroid', 'tonsil',
				 'uterus', 'vagina')

	if(order == TRUE) {
		events %<>%
			bind_rows(data_frame(
				gene = as.character(unlist(events[1,'gene'])),
				tissue = tissues, level = NA)) %>%
			mutate(gene = factor(gene, levels = unique(gene))) %>%
			mutate(tissue = factor(tissue, levels = c(vital, non.vital, sort(tissues %>% list.filter(! . %in% c(vital, non.vital))))))

		if(faceting == TRUE) {
			first.gene <- events[1, 'gene'] %>% unlist %>% unique %>% na.omit %>% as.character
			first.split <- events %>% filter(gene == first.gene) %$% split %>% na.omit %>% unique
			events %<>% mutate(split = ifelse(gene == first.gene & is.na(split), first.split, split))
		}
	}

	if(all(na.omit(events$level) %% 1 == 0)) {  # check if integer, if so plot discrete
		events %<>%
			mutate(level =
				ifelse(level == 0, 'not detected',
				ifelse(level == 1, 'low',
				ifelse(level == 2, 'medium',
				ifelse(level == 3, 'high',
				NA))))) %>%
			mutate(level = factor(level, levels = unique(level))) %>%
			arrange(desc(is.na(level)))

		m.gg <-
			ggplot(events, aes(tissue, gene)) +
			geom_tile(aes(fill = level), colour = 'grey') +
			scale_fill_manual(
				breaks   = names(palette),
				values   = palette,
				na.value = 'grey',
				drop     = FALSE,
				guide    = guide_legend(reverse = TRUE))
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
			facet_grid(split ~ ., scales = 'free_y', space = 'free_y')
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
		plot(mt.gg)
	}
}


PlotExpDist <- function(exp.table, binwidth = 0.02, method = 'BFGS', sigma.scalar = 1) {

	exp.table %<>%
		filter(!is.na(exp)) %>%
		arrange(exp) %>%
		mutate(num = row_number() - 1)

	NormFit <- function(mu, sigma) { -sum(log(dnorm(x, mu, sigma))) }

	ml <-
		bbmle::mle2(
			NormFit,
			start     = list(mu = mean(exp.table$exp), sigma = sd(exp.table$exp)),
			data      = list(x = exp.table$exp),
			method    = method) %>%
		tidy

	ml.mu <- ml %>% filter(term == 'mu') %$% estimate
	ml.sigma <- ml %>% filter(term == 'sigma') %$% estimate

	y.scaler <- ceiling(max(density(exp.table$exp)$y) * 10) * 0.1

	gg.dens <-
		ggplot(exp.table, aes(x = exp)) +
		geom_histogram(aes(y = ..density..), alpha = 0.5, binwidth = binwidth) +
		scale_x_continuous(breaks = floor(min(exp.table$exp)):ceiling(max(exp.table$exp))) +
		scale_y_continuous(breaks = (0:(20 * y.scaler)) / 20) +
		stat_function(
			fun   = dnorm,
			color = 'blue',
			args  = list(mean = ml.mu, sd = ml.sigma)) +
		geom_vline(xintercept = ml.mu, linetype = 2, color = 'black') +
		geom_vline(xintercept = ml.mu + ml.sigma, linetype = 2, color = 'purple') +
		geom_vline(xintercept = ml.mu - ml.sigma, linetype = 2, color = 'purple') +
		xlab("Expression") +
		ylab("Density")

	gg.quant <-
		ggplot(exp.table) +
		geom_point(aes(x = num, y = exp), size = 0.5, alpha = 0.25) +
		scale_x_continuous(breaks = seq(0, nrow(exp.table), nrow(exp.table) / 10), labels = (0:10) / 10) +
		scale_y_continuous(breaks = pretty_breaks(8), limits = c(NA, ceiling(max(exp.table$exp)))) +
		scale_alpha_manual(0.2) +
		geom_hline(aes(yintercept = ml.mu), color = 'black', linetype = 2) +
		geom_text(
			data = data.frame(x = 0, y = ml.mu),
			aes(x, y),
			label = str_c('µ = ', round(ml.mu, 2)),
			hjust = -1.65,
			vjust = -0.8) +
		geom_text(
			data = data.frame(x = 0, y = ml.mu),
			aes(x, y),
			label = str_c('µ + ', ifelse(sigma.scalar == 1, '', sigma.scalar), 'σ = ', round(ml.mu + ml.sigma * sigma.scalar, 2)),
			hjust = -1,
			vjust = -7.2) +
		geom_text(
			data = data.frame(x = 0, y = ml.mu),
			aes(x, y),
			label = str_c('µ - ', ifelse(sigma.scalar == 1, '', sigma.scalar) ,'σ = ', round(ml.mu - ml.sigma * sigma.scalar, 2)),
			hjust = -1,
			vjust = 5) +
		geom_hline(aes(yintercept = ml.mu + ml.sigma * sigma.scalar), color = 'purple', linetype = 2) +
		geom_hline(aes(yintercept = ml.mu - ml.sigma * sigma.scalar), color = 'purple', linetype = 2) +
		xlab("Fraction") +
		ylab("Expression")

	list(ml = ml, gg.dens = gg.dens, gg.quant = gg.quant)
}


PlotExpRel <- function(exp.table, plot.title, height = 4, width = 4, facet = TRUE, y.min = 0, y.max = 1, colors) {

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

SavePlot <- function(gg, file, width = 7, height = 7, type = 'pdf') {
	if(type == 'pdf') {
		pdf(file = file, width, height)
			plot(gg)
		dev.off()
	} else {
		png(file = file, width, height, units = 'in', res = 300)
			plot(gg)
		dev.off()
	}
}


CompPlot <- function(exp.table, title) {

	gg <-
		ggplot(exp.table, aes(db, n)) +
		ggtitle(title) +
		geom_bar(aes(fill = level), position = 'dodge', stat = 'identity') +
		scale_fill_brewer(palette = 'Set2', direction = -1)

	pdf(paste(c(title, 'pdf'), collapse = '.'))
		plot(gg)
	dev.off()

}
