if (!requireNamespace('BiocManager',quietly=TRUE)) install.packages('BiocManager',repos='https://cloud.r-project.org')
BiocManager::install(c('edgeR','statmod'),version='3.20',ask=FALSE,update=FALSE)
install.packages(c('data.table','jsonlite'),repos='https://cloud.r-project.org')
stopifnot(as.character(packageVersion('edgeR'))=='4.4.0')
