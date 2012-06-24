library(biomaRt)

listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)



ensid = 'ENSG00000213218'
data = getBM(attributes=c('ensembl_gene_id','human_paralog_ensembl_gene','paralog_hsap__dm_description_4014'),filters='ensembl_gene_id',values=ensid,mart=ensembl,uniqueRows = TRUE)

data = getBM(attributes=c('ensembl_gene_id', '5_utr_start', '5_utr_end'),filters='biotype',values="protein_coding",mart=ensembl,uniqueRows = TRUE)
