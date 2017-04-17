# content analysis of the reviewed literature

## setup

rm(list = ls())
# working directory
setwd("~/Documents/CQM/network analysis of clinical data/")
# packages
pkgs <- c("tidyr", "dplyr", "xtable", "ggplot2")
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    stopifnot(require(pkg, character.only = TRUE))
  }
}

# settings
source("scr/setup/specs.r")

# load and tidy data
source("scr/setup/tidy.r")

## step plot

# cumulative number of articles
cum_dat <- incl_dat[, c("google_scholar_cluster", "debut", "authors")]
cum_dat <- cum_dat[order(cum_dat$debut), ]
cum_dat$debut <- as.Date(cum_dat$debut)
cum_dat$tot <- 1:nrow(cum_dat)
# plot
step_plot <- ggplot(data = cum_dat,
                    aes(x = debut, y = tot)) +
  geom_step() +
  xlab("Date of publication") + ylab("Cumulative number of publications") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme_bw()
print(step_plot)
pdf(height = bodwid * .5, width = bodwid, file = "fig/step_plot.pdf")
print(step_plot)
dev.off()
# JAMIA figure
print(step_plot +
        xlab("Date") + ylab("Cumulative publication count") +
        scale_x_date(date_breaks = "4 years", date_labels = "%Y"))
ggsave(file = "fig/Fig2-old.eps", width = 80, height = 60, units = "mm")

## bar plot

# number of articles each year
yr_dat <- incl_dat %>% select(google_scholar_cluster, debut, type) %>%
  mutate(year = as.numeric(substr(debut, 1, 4)),
         type = factor(type, levels = c("Conference Paper", "Book Section",
                                        "Manuscript", "Journal Article")))
# plot
bar_plot <- ggplot(data = yr_dat, aes(x = year, group = type)) +
  geom_bar(aes(fill = type), colour = "black", size = 1/3) +
  scale_fill_brewer(type = "qual", palette = "Set3") +
  xlab("Year of publication") + ylab("Number of publications") +
  guides(fill = guide_legend(title = "Type of publication")) +
  theme_bw()
print(bar_plot)
pdf(height = bodwid * .5, width = bodwid, file = "fig/bar_plot.pdf")
print(bar_plot)
dev.off()
# JAMIA figure
bar_plot <- ggplot(data = yr_dat, aes(x = year)) +
  geom_bar(colour = "black", size = 1/3) +
  xlab("Year of publication") + ylab("Number of publications") +
  theme_bw()
print(bar_plot)
print(bar_plot + xlab("Year") + ylab("Number of publications"))
ggsave(file = "fig/Fig2.eps", width = 80, height = 60, units = "mm")

## label distributions

# publication types
print(table(incl_dat$type))

# whittle entries reporting the same study (currently at most 2)
print(table(is.na(incl_dat$same_study)))
incl_dat <- incl_dat[-which(!is.na(incl_dat$same_study)), ]

# distribution of group assigments
print(table(incl_dat$review_group))

# keyword frequencies
keyword_table <- sort(table(unlist(strsplit(tolower(incl_dat$keywords),
                                            split = "\\|"))),
                      decreasing = TRUE)
# table of keywords with frequencies > 2
keyword_leq2 <- which(keyword_table > 2)
keyword_data <- data.frame(
  Keyword = unname(names(keyword_table[keyword_leq2])),
  Frequency = as.vector(unname(keyword_table[keyword_leq2]))
)
print(keyword_data)
cat(print(xtable(x = keyword_data), include.rownames = FALSE),
    file = "tab/keywords.tex")

# conceptual frameworks
sort(table(incl_dat$conceptual_framework), decreasing = TRUE)

# methodological frameworks
sort(table(incl_dat$methodological_framework), decreasing = TRUE)

#  framework tables (counts)
tab <- incl_dat %>%
  filter(type == "Journal Article") %>%
  select(conceptual_framework, methodological_framework) %>%
  table()
write.csv(tab, file = paste0("tab/frameworks.csv"))
grps <- unique(incl_dat$review_group)
for (i in 1:length(grps)) {
  grp <- grps[i]
  tab <- incl_dat %>%
    filter(review_group == grp) %>%
    select(conceptual_framework, methodological_framework) %>%
    table()
  write.csv(tab, file = paste0("tab/frameworks-", i, ".csv"))
}

# framework tables (citations)
incl_dat <- mutate(incl_dat, label = paste0(
  gsub("^([^,]+),.*$", "\\1", authors),
  ifelse(grepl("\\|", authors),
         ifelse(grepl("\\|.+\\|", authors),
                " & al",
                paste0(" & ", gsub("^[^\\|]+\\|([^,]+),.*$",
                                   "\\1", authors))),
         ""),
  ", ", year
))
method_dat <- incl_dat %>%
  filter(type == "Journal Article") %>%
  select(cf = conceptual_framework, mf = methodological_framework, label) %>%
  mutate(cf = ifelse(is.na(cf), "-", cf),
         mf = ifelse(is.na(mf), "-", mf)) %>%
  group_by(cf, mf) %>%
  summarize(Studies = paste(label, collapse = "; ")) %>%
  rename(`Conceptual framework` = cf,
         `Methodological framework` = mf)
write.csv(method_dat, file = "tab/frameworks.csv")

# software
incl_dat$software %>% strsplit(" > |, ") %>% unlist() %>%
  table() %>% sort(decreasing = TRUE)

# frequency table stratified by literature subset
incl_dat$subset <- LETTERS[sample(1:4, size = nrow(incl_dat), replace = TRUE)]
strat_freq_table <- function(var, ...) {
  var_list <- strsplit(incl_dat[[var]], split = "\\|")
  var_counts <- sapply(var_list, length)
  var_set_dat <- data.frame(Item = unlist(var_list),
                            Subset = rep(incl_dat$subset, var_counts))
  var_set_tab <- table(var_set_dat)
  var_set_tab <- cbind(var_set_tab, Total = as.integer(rowSums(var_set_tab)))
  var_set_tab <- var_set_tab[rownames(var_set_tab) != "0", ]
  var_set_tab[order(var_set_tab[, "Total"], ...), ]
}
strat_freq_xtable <- function(...) {
  tab <- strat_freq_table(...)
  xtable(x = tab, align = c("l|", rep("r", ncol(tab) - 1), "|r"))
}
# level of analysis
cat(print(strat_freq_xtable(var = "level_of_analysis")),
    file = "tab/levels.tex")
cat(print(strat_freq_xtable(var = "structure_of_interest")),
    file = "tab/structure.tex")
cat(print(strat_freq_xtable(var = "tools")),
    file = "tab/tools.tex")
