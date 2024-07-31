library(readr)
library(tidyverse)
library(openxlsx)

dat <- read_xlsx("data/CORRECTED_40_fatty_amidases_from_Thierry_with_genus.xlsx") %>%
  select(p_notation,names) %>%
  separate(names, into = c("accession_number", "host"), sep = "_", extra = "merge", remove = FALSE) %>%
  mutate(host = gsub("_", " ", host)) %>%
  select(-names) %>%
  mutate(host = sub("([a-zA-Z]+ [a-zA-Z]+).*", "\\1", host))

dat$p_notation[1] <- "p204"
dat$p_notation[2] <- "p203"
dat$p_notation[40] <- "p205"
dat$host[1] <- "Gordonia terrae"
dat$host[2] <- "Corynebacterium aurimucosum"
dat$host[40] <- "Lacticaseibacillus rhamnosus"
dat$p_notation <- toupper(dat$p_notation)


output_file <- "output/20240718_table_SI_p_notation.xlsx"  # Change this to your desired file path
write.xlsx(dat, output_file, rowNames = FALSE)