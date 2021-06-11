## Prepare your datasets download from GISAID and Nextclade
# A combined dataset containing strain, date, region, country, division, 
# day from the initial sequencing, month, clade, total mutations, substitutions
# Input: metadata tsv file downloaded from GISAID and Nextclade results
# Output: 
# 1. a combined dataset containing strain, date, region, country, division, 
# day from the initial sequencing, month, clade, total mutations, substitutions
# 2. a dataset with separate substitutions (each row is a substitution)
# one strain could have more than one rows
all_information <- function(region, gisaid, nextclade){
  # GISAID variables
  meta <- read_tsv(gisaid) # read in metadata from GISAID
  meta <- meta %>% select(c("strain", "date", "region_exposure", 
                            "country_exposure", "division_exposure")) %>% 
    mutate(day = as.numeric(difftime(date, min(date), units = "days")),
           month = format(date, format = "%Y-%m"))
  # nextclade variables
  clade <- read_tsv(nextclade) # read in information from Nextclade
  clade <- clade %>% select(c("seqName", "clade", "totalSubstitutions",
                                 "substitutions"))
  combine <- inner_join(meta, clade, by = c("strain" = "seqName"))

  write_tsv(combine, paste0("data/", region, nrow(combine), ".tsv"))
  # separate substitutions
  mutation <- separate_rows(combine, "substitutions", convert = TRUE)
  mutation2 <- mutation %>% separate(substitutions, 
                                  into = c("reference", "pos_mut"),
                                  sep = 1, remove = FALSE, convert = TRUE) %>%
    separate(pos_mut, into = c("position", "mutation"), sep = -1, convert = TRUE)
  write_tsv(mutation2, paste0("data/", region, nrow(combine), "mut.tsv"))
}