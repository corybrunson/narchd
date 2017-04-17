# content analysis of the reviewed literature

# working directory
setwd("~/Documents/CQM/network analysis of clinical data/")
# packages
pkgs <- c("gdata", "xtable")
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    stopifnot(require(pkg, character.only = TRUE))
  }
}
github_pkgs <- c("corybrunson/context")
for (pkg in github_pkgs) {
  pkg_name <- gsub("^[A-Za-z0-9\\._]+/(.+)$", "\\1", pkg)
  if (!require(pkg_name, character.only = TRUE)) {
    devtools::install_github(pkg)
    stopifnot(require(pkg_name, character.only = TRUE))
  }
}

# settings
source("scr/setup/specs.r")

# read Excel file
incl_dat <- read.xls("data/nacd inclusion-new.xlsx",
                     header = TRUE, skip = 0,
                     colClasses = "character", na.strings = "")
# rename
names(incl_dat) <- gsub("\\.", "_", tolower(names(incl_dat)))

# conceptual scales
medium_scale <- make_scale_grep(objs = unique(incl_dat$medium),
                                patterns = c("Journal", "Book|Conference"))
#data_type_scale <- scale_standard(objs = unique(incl_dat$data_type),
#                                  scale = "nominal")
# WANT TO BE ABLE TO INPUT ORDER RELATIONS BETWEEN VALUES
construction_scale <- make_scale_grep(
  objs = incl_dat$construction,
  patterns = c("affiliation", "cooccurrence", "sequence", "designation")
)
# WANT FUNCTION TO ACCEPT LOGICAL RELATIONS AMONG UNNAMED VARIABLES
directed_scale <- make_scale_grep(objs = incl_dat$graph_type,
                                  patterns = c("directed"))
level_scale <- make_scale_grep(objs = incl_dat$level_of_analysis,
                               patterns = c("micro", "meso", "macro"))
subset_scale <- make_scale_standard(objs = incl_dat$subset, scale = "nominal")
# "NTN" contains both! so this should work fine
#theory_scale <- make_scale_grep(objs = incl_dat$theory,
#                                patterns = c("TN", "NT"))

# labels
incl_dat$label <- paste0(
  gsub("^([^,]+),.*$", "\\1", incl_dat$authors),
  ifelse(grepl("\\|", incl_dat$authors),
         ifelse(grepl("\\|.+\\|", incl_dat$authors),
                " & al ",
                paste0(" & ", gsub("^[^\\|]+\\|([^,]+),.*$",
                                   "\\1", incl_dat$authors), " ")),
         " "),
  "(", incl_dat$year, ")"
)
# label duplicates with suffixes a, b, ...
dupe_label <- function(x) {
  max_dupes <- max(table(x))
  stopifnot(max_dupes <= 26)
  if (max_dupes == 1) return(x)
  for (i in 1:(max_dupes - 1)) {
    firsts <- rev(setdiff(length(x) + 1 - which(duplicated(rev(x))),
                          which(duplicated(x))))
    x[duplicated(x)] <- paste0(gsub("[a-z]*)$", "", x[duplicated(x)]),
                               letters[i + 1], ")")
    if (i == 1) {
      x[firsts] <- paste0(gsub("[a-z]*)$", "", x[firsts]), letters[i], ")")
    }
  }
  x
}
incl_dat$label <- dupe_label(incl_dat$label)

# graph formal context
graph_ctx <- scale_by_variables(
  data = incl_dat,
  scales = list(construction = construction_scale,
                graph_type = directed_scale)
)
colnames(graph_ctx)[5] <- "directed" # ???
write_context(graph_ctx, file = "ctx/graph_ctx.csv", format = "concepts")

# level formal context
level_ctx <- scale_by_variables(incl_dat,
                                scale = list(level_of_analysis = level_scale))
write_context(level_ctx, file = "ctx/level_ctx.csv", format = "concepts")

# aggregate formal context
lit_ctx <- scale_by_variables(
  data = incl_dat,
  scales = list(medium = medium_scale,
                construction = construction_scale,
                graph_type = directed_scale,
                level_of_analysis = level_scale))
#rownames(lit_ctx) <- incl_dat$label # use format [ABC2012]
colnames(lit_ctx)[7] <- "directed" # ???
write_context(lit_ctx, file = "ctx/lit_ctx.csv", format = "concepts")

rm(list = ls())
