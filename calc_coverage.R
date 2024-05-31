library("data.table")
library("stringr")
library("dplyr")
library("ggplot2")
library("gridExtra")


# to assess the library saturation per gene
# recover fraction of TA sites disrupted per locus (not intergenic regions, introduce too much noise)
# plot histogram of library saturation


##### batch TAHitfraction files ################################
myfiles <- list.files("C:/git/Ranalysis/tamapper_results", pattern= "__TAHitFraction.txt", full.names = TRUE)

list_id <- sub(".*results/(.*)__TA.*", "\\1", myfiles) # extract sample id from file name

ta_hit <- lapply(myfiles, function(x){
  sample_id <- sub(".*results/(.*)__TA.*", "\\1", x) # extract sample id from file name
  file <- fread(x)  # read TAHitFraction files
  file$id <- sample_id
  file
})

mydataset <- rbindlist(ta_hit)  #make a single dataset


####################################
## clean dataset ###################

#remove not used columns 3, 4, 6, 7
df <- mydataset[, -c(3,4,6,7)]

#remove IG from Locus
df <- filter(df, str_detect(df$Locus,"^IG", negate = TRUE))

#rename variable
names(df)[names(df) == "Fraction TAs hit"] <- "Fraction_TAs_hit"
# test <- filter(df, df$id == "A1")

df2 <- filter(df, df$id %in% c("A1", "A2", "A3"))
df2 <- filter(df, df$id %in% c("A4", "A5", "A6"))
df2 <- filter(df, df$id %in% c("B1", "B2", "B3"))
df2 <- filter(df, df$id %in% c("B4", "B5", "B6"))
df2 <- filter(df, df$id %in% c("C1", "C2", "C3"))
df2 <- filter(df, df$id %in% c("C4", "C5", "C6"))
df2 <- filter(df, df$id %in% c("D1", "D2", "D3"))
df2 <- filter(df, df$id %in% c("D4", "D5", "D6"))
df2 <- filter(df, df$id %in% c("D7", "D8", "D9"))
df2 <- filter(df, df$id %in% c("E1", "E2", "E3", "E4"))
df2 <- filter(df, df$id %in% c("F1", "F2", "F3", "F4"))
df2 <- filter(df, df$id %in% c("G1", "G2", "G3", "G4"))
df2 <- filter(df, df$id %in% c("H1", "H2", "H3", "H4", "H5"))

df2 <- filter(df, df$id %in% c("B1", "B2", "B3", "C1", "C2", "C3",
                              "D1", "D2", "D3"))

wt_michael <- fread("C:/git/Ranalysis/tamapper_results/WT__TAHitFraction.txt")
wt_michael$id <- "wt"
wt_michael <- wt_michael[, -c(3,4,6,7)]
wt_michael <- filter(wt_michael, str_detect(wt_michael$Locus,"^IG", negate = TRUE))
names(wt_michael)[names(wt_michael) == "Fraction TAs hit"] <- "Fraction_TAs_hit"

df2 <- wt_michael

####################################
#### PLOT LIBRARY SATURATION ######
####################################

p1 <- ggplot(df2, aes(x=Fraction_TAs_hit)) +
  geom_histogram(binwidth=0.05, color="black", fill="aquamarine3")

p2 <- ggplot(df2, aes(x=Fraction_TAs_hit)) +
  geom_histogram(binwidth=0.05, color="black", fill="cornflowerblue") +
  facet_wrap(facets = vars(factor(id)), ncol = 6)

grid.arrange(p2, p1, widths = c(2.5,1))
