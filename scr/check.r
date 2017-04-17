### check literature data for errors

## setup

# working directory
setwd("~/Documents/CQM/network analysis of clinical data/")
# packages
pkgs <- c("gdata", "zoo")
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    stopifnot(require(pkg, character.only = TRUE))
  }
}

# settings
source("scr/setup/specs.r")

## data

# read Excel file (CHANGE DEPENDING ON STEP)
incl_dat <- read.xls("data/narchd.xlsx",
                     header = TRUE, skip = 0,
                     colClasses = "character", na.strings = "")
# rename
names(incl_dat) <- gsub("\\.", "_", tolower(names(incl_dat)))

## standardize dates

# publication date, using latest possible (e.g. Dec 31 if only year is given)
sea_mo <- c(Win = "Feb", Spr = "May", Sum = "Aug", Aut = "Nov", Fal = "Nov")
# assign entries without dates (e.g. book chapters) to December
mos <- ifelse(is.na(incl_dat$date) | incl_dat$date == "Date", "Dec", ifelse(
  substr(incl_dat$date, 1, 3) %in% names(sea_mo),
  sea_mo[substr(incl_dat$date, 1, 3)],
  # use the second month if more than one are listed
  substr(gsub("^[A-Z][a-z]{2}-([A-Z][a-z]{2})", "\\1", incl_dat$date), 1, 3)
))
das <- as.numeric(gsub(".*( |-)([0-9]+)$", "\\2", incl_dat$date))
stopifnot(!any(is.na(mos) & !is.na(das)))
wh <- which(is.na(das))
das[wh] <- format(as.Date(paste(
  # year, *following* month, and day 1
  as.numeric(incl_dat$year[wh]) + as.numeric(mos[wh] == "Dec"),
  ifelse(mos[wh] == "Dec",
         1,
         as.numeric(format(as.Date(paste(
           incl_dat$year[wh],
           mos[wh],
           1
         ), format = "%Y %b %d"), "%m")) + 1),
  1
), format = "%Y %m %d") - 1, "%d")
date2 <- as.Date(paste(incl_dat$year, mos, das), format = "%Y %b %d")
# e-date, which is always a precise date
edate2 <- as.Date(incl_dat$edate)
# minimum of these two
incl_dat$debut <- pmin(date2, edate2, na.rm = TRUE)
# remove temp fields
rm(sea_mo, mos, das, wh, date2, edate2)
# write to CSV
write_dat <- data.frame(
  gsc = incl_dat$google_scholar_cluster,
  debut = incl_dat$debut
)
write.csv(write_dat, file = "data/gsc_debuts.csv")

## disambiguate author names

authors <- unlist(strsplit(incl_dat$authors, split = "\\|"))
uniq_authors <- unique(authors)
lastnames <- gsub("(^[^,]+),[^,]+$", "\\1", uniq_authors)
tab <- table(lastnames)[table(lastnames) > 1]
for (nm in names(tab)) {
  print(unique(authors[grep(paste0(nm, ","), authors)]))
}

## disambiguate journal/conference/book names

forums <- incl_dat$journal_conference_book
tab_forums <- table(forums)
data.frame(forum = substr(names(tab_forums[tab_forums == 1]), 1, 60),
           freq = unname(as.vector(tab_forums[tab_forums == 1])))
data.frame(forum = substr(names(tab_forums[tab_forums > 1]), 1, 60),
           freq = unname(as.vector(tab_forums[tab_forums > 1])))

rm(list = ls())
