# read Excel file
incl_dat <- gdata::read.xls("data/narchd.xlsx",
                            header = TRUE, skip = 0,
                            colClasses = "character", na.strings = c("", "–"))
# rename
names(incl_dat) <- gsub("\\.", "_", tolower(names(incl_dat)))
# remove entries marked for removal
incl_dat <- subset(incl_dat, is.na(remove))
# remove entries assigned to the same Google Scholar cluster
duped <- duplicated(incl_dat$google_scholar_cluster) &
  !is.na(incl_dat$google_scholar_cluster)
incl_dat <- incl_dat[-which(duped), ]
print(dim(incl_dat))
