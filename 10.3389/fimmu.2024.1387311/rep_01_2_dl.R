# Figure 01

library(GEOquery)

data_path = './10.3389/fimmu.2024.1387311/data/'

# gse_datasets = c("GSE55235", "GSE77298", "GSE12021", "GSE55457")
# eSet <- getGEO(gse_num, destdir = 'data', getGPL= F)

GSE55235 = getGEO("GSE55235", destdir = data_path, getGPL= F)
GSE55235 = GSE55235[[1]]
exp_55235 = exprs(GSE55235)

GSE77298 = getGEO("GSE77298", destdir = data_path, getGPL= F)
GSE77298 = GSE77298[[1]]
exp_77298 = exprs(GSE77298)

GSE12021 = getGEO("GSE12021", destdir = data_path, getGPL= F)
GSE12021 = GSE12021[[1]]
exp_12021 = exprs(GSE12021)

GSE55457 = getGEO("GSE55457", destdir = data_path, getGPL= F)
GSE55457 = GSE55457[[1]]
exp_55457 = exprs(GSE55457)

headers = c('Datasets', 'Platform', 'Total samples', 'Normal', 'RA')

library(stringr)

dat_55235 = c(
  "GSE55235",
  GSE55235@annotation,
  ncol(exp_55235),
  sum(str_detect(GSE55235$characteristics_ch1, 'healthy')),
  sum(str_detect(GSE55235$characteristics_ch1, 'rheumatoid'))
  # ncol(exp_55235) - sum(str_detect(GSE55235$characteristics_ch1, 'healthy'))
)

dat_77298 = c(
  "GSE77298",
  GSE77298@annotation,
  ncol(exp_77298),
  sum(str_detect(GSE77298$characteristics_ch1, 'control')),
  sum(str_detect(GSE77298$characteristics_ch1, 'arthritis'))
)

dat_12021 = c(
  "GSE12021",
  GSE12021@annotation,
  # ncol(exp_12021),
  sum(!str_detect(GSE12021$characteristics_ch1.2, 'osteoarthritis')),
  sum(str_detect(GSE12021$characteristics_ch1.2, 'control')),
  sum(str_detect(GSE12021$characteristics_ch1.2, 'rheumatoid'))
)

dat_55457 = c(
  "GSE55457",
  GSE55457@annotation,
  sum(!str_detect(GSE55457$characteristics_ch1.2, 'osteoarthritis')),
  sum(str_detect(GSE55457$characteristics_ch1.2, 'control')),
  sum(str_detect(GSE55457$characteristics_ch1.2, 'rheumatoid'))
)

all_data = data.frame(rbind(dat_55235, dat_77298, dat_12021, dat_55457))
colnames(all_data) = headers
rownames(all_data) = c(1, 2, 3, 4)
View(all_data)

save(
  GSE55235, exp_55235,
  GSE77298, exp_77298,
  GSE12021, exp_12021,
  GSE55457, exp_55457,
  data_path,
  file = paste0(data_path, 'rep_01_2_dl.RData')
)
