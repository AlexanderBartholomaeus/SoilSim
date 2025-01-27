# Script DESeq2 analysis

# extended manual is available here
#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts

# load library 
library(DESeq2)

# loop over KEGG KOs, modules and pathways
my_loop <- c('KO','module','pathway')
for(i in my_loop){

  # load data
  counts <- read.table(paste0('data/kegg_',i,'_sum.csv'),
                       sep=';', header=T, stringsAsFactors = F)

  # set rownames
  rownames(counts) <- counts[,1]
  counts <- counts[,-1]
  
  
  # build table for conditions ! same order needed
  cond <- data.frame(row.names = colnames(counts), 
                     location = rep(c('PA', 'SG'),each=4), 
                     tp = rep(c('T1','T3','T1','T3'),each=2),
                     plant = rep(c('insitu','plant'), 4),
                     stringsAsFactors = F)
  
  # make factor
  cond$location <- factor(cond$location)
  cond$tp <- factor(cond$tp)
  cond$plant <- factor(cond$plant)
  
  # perform DESeq analysis
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = cond,
                                design = ~ location) # location
                                #design = ~ tp) # timepoint
                                #design = ~ plant) # plant 
  dds <- DESeq(dds)
  res <- results(dds)
  
  #res
  summary(res)
  resOrdered <- res[order(res$pvalue),]
  
  # write results
  write.csv(as.data.frame(resOrdered), 
            file=paste0('DESeq2_results_',i,'_location.csv'))
}  
