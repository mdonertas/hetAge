## code to prepare `brain_age_expression` dataset goes here

library(SummarizedExperiment)

counts <- readRDS('data-raw/exps_not_scaled.rds')
ages <- readRDS('data-raw/ages.rds')
names(ages) = gsub('Somel2010','Somel2011',names(ages))

for(dname in names(ages)){
    cnts <- counts[[dname]]
    agx <- (ages[[dname]]^4)/365
    agx <- agx[colnames(cnts)]
    colData <- DataFrame(ages = agx)
    assign(dname,SummarizedExperiment(assays = list(counts = cnts),
                                      colData = colData))
    save(list = dname, file = paste('data/',dname,'.RData', sep=''))
}
