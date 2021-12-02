#
# Compute from VCF files
# pairwise FST (median of local Fst)
#
# Ophélie Da Silva
# December, 2021
#
#
#
# input: VCF files
# output: txt containing pairwise Fst for all organisms
#------------------------------------------------------------

# options
options(stringsAsFactors = F)

# import libraries
library('tidyverse')

# cleaning
rm(list=ls())
graphics.off()

# function
extract_n = function(x,n) {
  k <- stringr::str_split_fixed(x, ":", 2)[,n]
  k[which(k=="NA")] <- NA
  k <- as.numeric(k)
}

# define variable
input.dir='data/raw/snps'
output.dir='data/out'

s = 4 # vertical coverage threshold : 4X
vcf.list <- list()
fst.list <- list()
write.list <- list()

# list raw vcf files
files <- list.files(input.dir, ".vcf") %>% 
  paste(input.dir,.,sep='/') # vcf files
n_org <- length(files) # number of files
org_name <- gsub("^(.+-.+)_\\d+-\\d+_snp\\.vcf$", "\\1", basename(files)) # organism names corresponding to each file

for (i in 1:n_org) {
  org <- org_name[i]
  
  # Extract the colnames from the vcf file
  
  readLines(files[i]) %>%
    .[grep('#CHROM',.)] %>%      # the line begining with '#CHROM' contains column headers
    strsplit(., "\t") %>%        # headers are tab separated
    unlist %>%                   # create a vector
    gsub('#CHROM','CHROM',.) %>% # remove '#' symbol
    gsub('^genpop_.+_\\d+-\\d+-\\d+/alignment_filt/TARA_(\\d+)/ERR\\d+\\.filter\\.sort\\.bam$',
         'TO-\\1', .             # simplify sample name 
         ) -> vcf.names          # vcf.names contains the vcf column names
  
  vcf.names %>% 
    grep("(TO-)", .) %>% 
    length -> n_stations # number of stations
  
  
  # Extract the vcf data
  
  read.table(files[i]) %>% 
    `names<-`(vcf.names) %>%                        # rename columns
    .[,grep("(CHROM)|(POS)|(TO-)",vcf.names)] %>%   # select columns (remove: ID,REF,ALT,QUAL,FILTER,INFO,FORMAT)
    lapply(.,                                       # extract coverage information from GT:PL:DP:AD
           function(x) {
             gsub('(\\d+|\\.):(\\d+,\\d+|\\.):(\\d+|\\.):(\\d+,(\\d+)|(\\.))',
                  '\\3:\\5\\6',x)
           }) %>% 
    as.data.frame -> vcf.list[[org]]
  
  # Replace by NA if no snps detected or more than 2 alleles
  vcf.list[[org]] %>%
    mutate_at(vars(contains("TO.")),
              .funs = function(x) {
                x <- ifelse(!grepl("^\\d+:\\d+$", x), NA, x)
              }) -> vcf.list[[org]]
  
  # Delete lines with only na values
  vcf.list[[org]] %>% 
    select(contains("TO.")) %>% 
    is.na() %>% 
    rowSums() != n_stations -> x
  which(x)
  
  vcf.list[[org]] <- vcf.list[[org]][which(x != n_stations),]
  
  # Clean
  rm(vcf.names, n_stations, x)
  
  # Preparation of AD and DP tables
  fst.list[[org]] <- list()
    
  # Extract AD values
  vcf.list[[org]] %>% 
    mutate_at(vars(contains("TO.")), .funs = function(x) extract_n(x,2)) -> fst.list[[org]][["AD"]]
  
  # Extract DP values
  vcf.list[[org]] %>% 
    mutate_at(vars(contains("TO.")), .funs = function(x) extract_n(x,1)) -> fst.list[[org]][["DP"]]

  # Filtering SNPs based on vertical coverage (filtering loci : x ∈ [µ±2σ] and x < s)
  
  fst.list[[org]][["DP"]] %>% 
    select(contains("TO.")) %>% 
    colMeans(na.rm = T) -> cov_mean              # mean coverage at each snps for each sample
  
  fst.list[[org]][["DP"]] %>% 
    select(contains("TO.")) %>% 
    apply(., 2, sd, na.rm=T) -> cov_sd           # coverage standard deviation at each snps for each sample
  
  up <- cov_mean + 2*cov_sd     # upper coverage threshold
  low <- cov_mean - 2*cov_sd    # lower coverage threshold
    
  fst.list[[org]][["DP"]] %>% 
    select(contains("TO.")) %>% 
    apply(., 1, function(x) {
      ifelse(x >= low & x <= up & x >= 4, x, NA)  # replace by na when extreme coverage
    }) %>%
    t %>%
    as.data.frame() -> fst.list[[org]][["DP.filt"]]
    
  # Clean
  rm(cov_mean, cov_sd, up, low)
  
  # Compute frequency matrix for the alternative allele
  fst.list[[org]][["AD"]] %>% 
    select(contains("TO.")) /
    (fst.list[[org]][["DP.filt"]] %>% 
    select(contains("TO."))) -> freq
  
  # Compute fst matrix
  
  nx <- nrow(freq)       # number of loci
  ny <- ncol(freq)       # number of pair of samples
  matrix(nrow = nx, ncol=ny*(ny-1)/2) %>%   # create a dataframe to store local Fst for pairwise samples
    data.frame -> fst.list[[org]]$fst
  
  freq %>% 
    names %>% 
    combn(. , 2) -> station_pair            # all pair of samples
  
  for (pair in 1:ncol(station_pair)) {
    progress <- (pair/ncol(station_pair)) %>% round(digits = 4)*100
    
    st1 <- station_pair[1, pair]
    st2 <- station_pair[2, pair]
    
    message("Progress: ", progress, "% (", st1, " and ", st2,")")
    
    var_loc <- apply(freq[,c(st1, st2)], 1, var)
    mean_loc <- apply(freq[,c(st1, st2)], 1, mean)

    fst.list[[org]]$fst[, pair] <- ifelse(is.na(var_loc), NA,
                                          ifelse(mean_loc*(1-mean_loc) != 0, 
                                                 var_loc/(mean_loc*(1-mean_loc)), 
                                                 NA))
    colnames(fst.list[[org]]$fst)[pair] <- paste(st1, st2, sep = "_")
    
    rm(st1, st2, var_loc)
  }
  
  # clean
  rm(freq, station_pair, pair, nx, ny)
  
  # compute median pairwise Fst (one value by pair of samples)
  fst.list[[org]]$fst %>% 
    apply(., 2, median, na.rm = T) -> fst.list[[org]]$median_fst
    
  
  # write datatable
  write.list[[org]] <- data.frame(pStation = fst.list[[org]]$median_fst %>% 
                                    names,
                                  Organism = org,
                                  Fst = fst.list[[org]]$median_fst)
}; rm(i)

# save file
fst <- do.call(rbind, write.list)
write.table(fst, 
            paste(output.dir, "median_pairwise_fst_meanloc.txt", sep = "/"),
            row.names = F, quote = F)
