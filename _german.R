german <- c('FUT3', 'HAVCR2', 'CD96', 'CD99', 'IL3RA', 'CLEC12A', 'FLT3', 'CD244', 'CD33')

master %>%
filter(hgnc %in% german) %>%
select(gene = hgnc, level = tissue.max, tissue) %>%
unique %>%
GeneRename %>%
PlotTissue(pdf = TRUE, faceting = FALSE, file.name = 'german_targets.pdf', width = 10, height = 2.5, order = TRUE) 

