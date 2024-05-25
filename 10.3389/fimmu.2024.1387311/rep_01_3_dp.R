data_path = './10.3389/fimmu.2024.1387311/data/'
load(paste0(data_path, 'rep_01_2_dl.RData'))


# first do log2 + 1

# any > 20?
if (any(exp_55235 > 20)) {
  exp_55235 = log2(exp_55235 + 1)
}

if (any(exp_12021 > 20)) {
  exp_12021 = log2(exp_12021 + 1)
}

if (any(exp_77298 > 20)) {
  exp_77298 = log2(exp_77298 + 1)
}

# probe ids to genes

GSE55235@annotation == 'GPL96'
GSE12021@annotation == 'GPL96'
GSE77298@annotation == 'GPL570'


# add probe names to each

library(tinyarray)

find_anno('GPL96')
library(hgu133a.db)
GPL96_ids <- toTable(hgu133aSYMBOL)

find_anno('GPL570')
library(hgu133plus2.db)
GPL570_ids <- toTable(hgu133plus2SYMBOL)

GPL96_ids = GPL96_ids[!duplicated(GPL96_ids$symbol),]
GPL570_ids = GPL570_ids[!duplicated(GPL570_ids$symbol),]


# drop rows that don't have a gene name
nrow(exp_55235)
exp_55235 = exp_55235[rownames(exp_55235) %in% GPL96_ids$probe_id,]
nrow(exp_55235)

nrow(exp_12021)
exp_12021 = exp_12021[rownames(exp_12021) %in% GPL96_ids$probe_id,]
nrow(exp_12021)

nrow(exp_77298)
exp_77298 = exp_77298[rownames(exp_77298) %in% GPL570_ids$probe_id,]
nrow(exp_77298)

rownames(exp_55235) = GPL96_ids$symbol[match(GPL96_ids$probe_id, rownames(exp_55235))]
rownames(exp_12021) = GPL96_ids$symbol[match(GPL96_ids$probe_id, rownames(exp_12021))]
rownames(exp_77298) = GPL570_ids$symbol[match(GPL570_ids$probe_id, rownames(exp_77298))]

# inner join the three datasets
exp_55235 = as.data.frame(exp_55235)
exp_55235$symbol = rownames(exp_55235)

exp_12021 = as.data.frame(exp_12021)
exp_12021$symbol = rownames(exp_12021)

exp_77298 = as.data.frame(exp_77298)
exp_77298$symbol = rownames(exp_77298)

merged_exp = merge(merge(exp_55235, exp_12021, by='symbol'), exp_77298, by='symbol')
rownames(merged_exp) = merged_exp$symbol
merged_exp$symbol = NULL
nrow(merged_exp)
ncol(merged_exp)

exp_55235$symbol = NULL
exp_12021$symbol = NULL
exp_77298$symbol = NULL
ncol(exp_55235) + ncol(exp_12021) + ncol(exp_77298)

boxplot(merged_exp)

# normalize it
merged_exp = limma::normalizeBetweenArrays(merged_exp)
boxplot(merged_exp)


# set group var:
# GSM : Normal / RA

pd_55235 = pData(GSE55235)
pd_12021 = pData(GSE12021)
pd_77298 = pData(GSE77298)


# drop osteoarthritis

library(stringr)

to_drop_55235 = ifelse(
  str_detect(
    pd_55235$`disease state:ch1`,
    'osteo'
  ),
  TRUE,
  FALSE
)

to_drop_12021 = ifelse(
  str_detect(
    pd_12021$`disease:ch1`,
    'osteo'
  ),
  TRUE,
  FALSE
)

na.omit(to_drop_12021)

to_drop_77298 = ifelse(
  str_detect(
    pd_77298$`disease state:ch1`,
    'osteo'
  ),
  TRUE,
  FALSE
)

to_drop = c(pd_55235$geo_accession[to_drop_55235], pd_12021$geo_accession[to_drop_12021], pd_77298$geo_accession[to_drop_77298])
to_drop = na.omit(to_drop)
merged_exp = merged_exp[, -which(colnames(merged_exp) %in% to_drop)]
boxplot(merged_exp)
ncol(merged_exp)

# now we need to set the group var

group_55235 = ifelse(
  str_detect(
    pd_55235$`disease state:ch1`,
    'healthy'
  ),
  'Control',
  'RA'
)
group_55235 = group_55235[!to_drop_55235]

group_12021 = ifelse(
  str_detect(
    pd_12021$`disease:ch1`,
    'rheumatoid'
  ),
  'RA',
  'Control'
)
group_12021 = group_12021[!to_drop_12021]
# replace NA to 'Control'
group_12021[is.na(group_12021)] = 'Control'

group_77298 = ifelse(
  str_detect(
    pd_77298$`disease state:ch1`,
    'control'
  ),
  'Control',
  'RA'
)

# nothing to drop

ncol(merged_exp) == length(group_55235) + length(group_12021) + length(group_77298)
group = c(group_55235, group_12021, group_77298)

save(
  merged_exp,
  group,
  file = paste0(data_path, 'rep_01_3_dp.RData')
)
