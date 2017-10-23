pdb.user = 'bermans'
email = 'bermans@mskcc.org'
pdb.pass = '60010Violet'

palette <- c(
'not detected' = '#202c99',
'low'          = '#fbeaea',
'medium'       = '#e77e7e',
'high'         = '#c12525')

pass.tissue <- c('bone', 'blood', 'spleen')

vital <- c(
'adipose tissue',
'adrenal',
'bladder',
'brain',
'bronchus',
'eye',
'gut',
'heart',
'kidney',
'esophagus',
'liver',
'lung',
'nasopharynx',
'oropharynx',
'pancreas',
'rectum',
'skeletal muscle',
'skin',
'smooth muscle',
'soft tissue',
'spinal cord',
'stomach')

non.vital <- c(
'appendix',
#'bone',
#'blood',
'breast',
'cerumen',
'cervix',
'epididymis',
'fallopian tube',
'gallbladder',
'lymph node',
'ovary',
'parathyroid',
'prostate',
'seminal',
'spleen',
'synovial fluid',
'testis',
'thyroid',
'tonsil',
'uterus',
'vagina')


non.vital.combi <- c(
#'appendix',
#'bone',
#'blood',
'breast',
'cerumen',
'cervix',
'epididymis',
'fallopian tube',
'gallbladder',
'lymph node',
'ovary',
'parathyroid',
'prostate',
'seminal',
#'spleen',
'synovial fluid',
'testis',
'thyroid',
'tonsil',
'uterus',
'vagina')


exclude.protein <- c(
	'SP140L', 'CIRH1A', 'DDX51', 'PDCD11', 'GNL2', 'NDC80', 'DDX54', 'DNAAF5', 'IMPA2', 'INTS3',
	'RBM12B', 'GTF3C4', 'SMEK2', 'CIAPIN1', 'CDK2', 'MAD2L1', 'ANKRD28', 'PSMG1', 'CHML', 'KNOP1',
	'GYPA', 'BRD9', 'C9orf41', 'CCNB3', 'CDC123', 'CDK13', 'CSNK1G2', 'DHX40', 'DNAH11', 'DPH1',
	'ANK1','UTP15', 'WBSCR22', 'WDR43', 'ZNF148', 'ZNHIT2', 'UBFD1', 'TRMT61A', 'TOR2A', 'SCOC',
	'COL9A2', 'ESRRA', 'FAM20B', 'FAM60A', 'FBXL12', 'GTF2A2', 'HELZ2', 'RRP7A', 'RBM19', 'PRPF38B',
	'PPP4R3B', 'POLK', 'PMF1', 'MTF2', 'LRRC41', 'MAN2B2', 'MED17', 'MED8', 'METTL17', 'MRPL43',
	'MRPS18A', 'MTEF3', 'ITGB3', 'KIF22', 'LST1', 'PIGK', 'PSMG3', 'SHCBP1', 'SURF1', 'TMEM167A',
	'ADAP1', 'AGPAT1', 'MPP1', 'UBE2S', 'GUCY2D', 'ARHGBAP33', 'EPB41', 'CD4', 'IGSF1', 'INSR',
	'LMF2', 'LYST', 'MICAL1', 'ORMDL3', 'PPP1R9B', 'PPP2R5A', 'PTBP3', 'RNF168', 'SREK1', 'TBXAS1',
	'TMEM63A', 'TSC2', 'VCPIP1', '2NF598', 'HLA-DRB4', 'ZNF598', 'PDE4DIP', 'ARHGAP33', 'PTPRK',
	'CARD11', 'RASGRP3', 'MFGE8', 'MS4A3', 'CBL', 'ARRB2', 'GRK5', 'IGFBP7')

compare.protein <- c('CLEC12A', 'IL3RA', 'FOLR2', 'FUT3', 'CD33', 'CD38', 'CD44')

compare.protein.dnmt3a <- c('CD33', 'CD38', 'CD44', 'CD47', 'ADGRE5', 'CD99')

exclude.lines <- c('EARLY_PM', 'LATE_PM', 'MY', 'MM', 'BC', 'PMN', 'MONO', 'MDS', 'T(15;17)', 'ALL')

progeniter.lines <- c('HSC', 'MPP', 'CMP', 'GMP', 'MEP')

aml.additions <- c(
	'ENSG00000197635', # CD26  / DPP4
	'ENSG00000185291', # CD123 / IL3RA
	'ENSG00000154096', # CD90  / THY1
	'ENSG00000153283', # CD96
	'ENSG00000196776', # CD47
	'ENSG00000026508', # CD44
	'ENSG00000143226', # CD32  / FCGR2A
	'ENSG00000134460', # CD25  / IL2RA
	'ENSG00000172322', # CLL1  / CLEC12A
	'ENSG00000135077', # TIM3  / HAVCR2
	'ENSG00000174059', # CD34
	'ENSG00000185499', # MUC1
	'ENSG00000002586', # CD99
	'ENSG00000091409', # CD49F / ITGA6
	'ENSG00000004468', # CD38
	'ENSG00000123146', # CD97 / ADGRE5
	'ENSG00000081237', # CD45 / PTPRC
	'ENSG00000149294', # CD56 / NCAM1
	'ENSG00000171124', # FUT3
	'ENSG00000125810', # CD93
	'ENSG00000165457'  # FOLR2
)
