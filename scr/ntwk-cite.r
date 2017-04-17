### graph representations of the revewed literature

## setup

rm(list = ls())
# working directory
setwd("~/Documents/CQM/network analysis of clinical data/")
# packages
pkgs <- c("zoo", "igraph")
for (pkg in pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    stopifnot(require(pkg, character.only = TRUE))
  }
}
github_pkgs <- c("hadley/ggplot2", "thomas85/ggforce", "thomas85/ggraph")
for (pkg in github_pkgs) {
  pkg_name <- strsplit(pkg, split = "/")[[1]][2]
  if (!require(pkg_name, character.only = TRUE)) {
    devtools::install_github(pkg)
    stopifnot(require(pkg_name, character.only = TRUE))
  }
}

# settings
source("scr/setup/specs.r")

# load and tidy data
source("scr/setup/tidy.r")

## citation graph

# nodes (unique identifiers)
cite_nodes <- incl_dat$google_scholar_cluster
# citations
win_citings <- strsplit(incl_dat$internal_citations, split = "\\|")
win_cites <- sapply(win_citings, length)
cites <- lapply(win_citings, function(x) {
  if (all(is.na(x))) return(NULL)
  i_vec <- c()
  for (y in x) {
    i <- grep(y, incl_dat$google_scholar_cluster)
    #i <- which(incl_dat$google_scholar_cluster == y)
    #if (length(i) == 0) i <- which(incl_dat$pmid == y)
    #if (length(i) == 0) i <- which(incl_dat$doi == y)
    if (length(i) == 0) next else {
      if (length(i) != 1) {
        print(x)
        print(y)
        print(i)
      }
      stopifnot(length(i) == 1)
      i_vec <- c(i_vec, i)
    }
  }
  i_vec
})
# directed edges (citations)
cite_edges <- do.call(c, lapply(1:length(cites), function(i) {
  if (length(cites[[i]]) == 0) NULL else as.vector(rbind(cites[[i]], i))
}))
# construct graph
cite_graph <- make_empty_graph(n = length(cite_nodes), directed = TRUE)
V(cite_graph)$name <- cite_nodes
cite_graph <- add_edges(graph = cite_graph, edges = cite_edges)

# which citations are not self-?
cite_edgelist <- as_edgelist(cite_graph)
article_authors <- incl_dat$authors
names(article_authors) <- incl_dat$google_scholar_cluster
cite_from_authors <- strsplit(article_authors[cite_edgelist[, 1, drop = TRUE]],
                              split = "\\|")
cite_to_authors <- strsplit(article_authors[cite_edgelist[, 2, drop = TRUE]],
                            split = "\\|")
E(cite_graph)$self <- sapply(1:ecount(cite_graph), function(i) {
  length(intersect(cite_from_authors[[i]], cite_to_authors[[i]])) > 0
})

# dates
V(cite_graph)$debut <- incl_dat$debut

# aesthetics
V(cite_graph)$win_cites <- win_cites
V(cite_graph)$all_cites <- as.numeric(incl_dat$citation_count)
V(cite_graph)$label <- paste0(
  gsub("^([^,]+),.*$", "\\1", incl_dat$authors),
  ifelse(grepl("\\|", incl_dat$authors),
         ifelse(grepl("\\|.+\\|", incl_dat$authors),
                " & al ",
                paste0(" & ", gsub("^[^\\|]+\\|([^,]+),.*$",
                                   "\\1", incl_dat$authors), " ")),
         " "),
  "(", incl_dat$year, ")"
)
V(cite_graph)$venue <- ifelse(
  grepl("(^| )([Pp]roc|[Cc]onf)", incl_dat$journal_conference_book) &
    !grepl("National Academy", incl_dat$journal_conference_book),
  "Conference",
  "Journal"
)
V(cite_graph)$literature <- ifelse(
  incl_dat$type == "Journal Article", "Peer-reviewed", "Grey"
)
V(cite_graph)$group <- incl_dat$review_group

# save the full citation graph
if (!file.exists("data/scientometrics")) dir.create("data/scientometrics")
save(cite_graph, file = "data/scientometrics/cite_graph.rdata")

# plot
# good seeds: 6, 16, 27, *38*, 41, 43, 47, 66, 84
# good seeds (start.temp * .1): 23, 39, *53*
set.seed(38)
cite_plot <- ggraph(cite_graph, layout = "fr", niter = 200) +
  geom_edge_fan(aes(alpha = ..index.., color = self)) +
  scale_edge_alpha(guide = "none") +
  scale_edge_color_manual(name = "Self-citation?",
                          labels = c("No", "Yes"),
                          values = c("black", "grey")) + #e41a1c
  geom_node_point(aes(size = V(cite_graph)$all_cites,
                      color = V(cite_graph)$group,
                      alpha = V(cite_graph)$literature)) +
  scale_size_continuous(name = "No. citations",
                        trans = "sqrt",
                        range = .3 * c(1, sqrt(max(V(cite_graph)$all_cites))),
                        breaks = seq(1, max(V(cite_graph)$all_cites),
                                     floor(max(V(cite_graph)$all_cites) /
                                             4))) +
  scale_color_brewer(name = "Group", type = "qual", palette = "Dark2") +
  scale_alpha_manual(name = "Literature", values = c(.5, 1)) +
  geom_node_text(aes(label = ifelse(
    V(cite_graph)$all_cites >= 30,
    V(cite_graph)$label, ""
  )), size = 3, repel = TRUE) +
  ggforce::theme_no_axes() +
  guides(alpha = guide_legend(override.aes = list(size = 4)),
         color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        legend.box = "horizontal")
print(cite_plot)
pdf(height = bodwid * 1.1, width = bodwid, file = "fig/cite_plot.pdf")
print(cite_plot)
dev.off()

# isolates published in the most recent 2 years
table(data.frame(isolate = degree(cite_graph) == 0,
                 recent = incl_dat$debut >= as.Date("2014-08-01")))
# statistical test of isolates and venue
tab <- table(data.frame(isolate = degree(cite_graph) == 0,
                        venue = incl_dat$type == "Journal Article"))
chisq.test(x = tab)

# subset to cited or citing articles
cite_subgraph <- delete_vertices(cite_graph, which(degree(cite_graph) == 0))
# plot
set.seed(9)
cite_subplot <- ggraph(cite_subgraph, layout = "fr", niter = 200) +
  geom_edge_fan(aes(alpha = ..index.., color = self)) +
  scale_edge_alpha(guide = "none") +
  scale_edge_color_manual(name = "Self-citation?",
                          labels = c("No", "Yes"),
                          values = c("black", "grey")) +
  geom_node_point(aes(size = V(cite_subgraph)$all_cites,
                      color = V(cite_subgraph)$group,
                      alpha = V(cite_subgraph)$literature)) +
  scale_size_continuous(name = "No. citations",
                        trans = "sqrt",
                        range = .3 * c(1, sqrt(max(V(cite_graph)$all_cites))),
                        breaks = seq(1, max(V(cite_graph)$all_cites),
                                     floor(max(V(cite_graph)$all_cites) /
                                             4))) +
  scale_color_brewer(name = "Group", type = "qual", palette = "Dark2") +
  scale_alpha_manual(name = "Literature", values = c(.5, 1)) +
  geom_node_text(aes(label = ifelse(V(cite_subgraph)$all_cites > 30,
                                    V(cite_subgraph)$label, "")),
                 repel = TRUE) +
  guides(alpha = guide_legend(override.aes = list(size = 4)),
         color = guide_legend(override.aes = list(size = 4))) +
  ggforce::theme_no_axes() +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        legend.box = "horizontal")
print(cite_subplot)
pdf(height = bodwid * 1.1, width = bodwid, file = "fig/cite_subplot.pdf")
print(cite_subplot)
dev.off()

# which articles are not linked by citations
wh_uncited <- which(degree(cite_graph) == 0)
# graph communities
cite_comm_mod <- cluster_fast_greedy(as.undirected(cite_subgraph),
                                     merges = TRUE, modularity = TRUE)
cite_comm_btw <- cluster_edge_betweenness(cite_subgraph,
                                          merges = TRUE, bridges = TRUE,
                                          modularity = TRUE, membership = TRUE)
cite_comm_wlk <- cluster_walktrap(cite_subgraph, steps = 4,
                                  merges = TRUE, modularity = TRUE,
                                  membership = TRUE)
# concordance between communities and thematic analysis groups
mclust::adjustedRandIndex(x = incl_dat$review_group[-wh_uncited],
                          y = cite_comm_mod$membership)
mclust::adjustedRandIndex(x = incl_dat$review_group[-wh_uncited],
                          y = cite_comm_btw$membership)
mclust::adjustedRandIndex(x = incl_dat$review_group[-wh_uncited],
                          y = cite_comm_wlk$membership)

# citations (arcs) of highest betweenness
arc_btw <- edge_betweenness(cite_subgraph)
arc_ends <- ends(cite_subgraph, E(cite_subgraph), names = F)
# undirected citation links of highest betweenness
edge_btw <- edge_betweenness(cite_subgraph, directed = FALSE)
edge_ends <- ends(cite_subgraph, E(cite_subgraph), names = F)
# data frame
btw_dat <- data.frame(
  Source = paste0(V(cite_subgraph)$label[edge_ends[, 1]], " (",
                 V(cite_subgraph)$name[edge_ends[, 1]], ")"),
  Target = paste0(V(cite_subgraph)$label[edge_ends[, 2]], " (",
                  V(cite_subgraph)$name[edge_ends[, 2]], ")"),
  DirBtw = arc_btw,
  UndirBtw = edge_btw
)
print(head(btw_dat[rev(order(edge_btw)), ], n = 18))
print(head(btw_dat[rev(order(arc_btw)), ], n = 18))

E(cite_subgraph)$hi_arc_btw <- arc_btw > 20
set.seed(9)
cite_subplot_alt <- ggraph(cite_subgraph, layout = "fr", niter = 200) +
  geom_edge_fan(aes(alpha = ..index.., color = hi_arc_btw)) +
  scale_edge_alpha(guide = "none") +
  scale_edge_color_manual(name = "Arc betweenness",
                          labels = c("< 20", "> 20"),
                          values = c("black", "red")) +
  #scale_edge_color_manual(name = "Self-citation?",
  #                        labels = c("No", "Yes"),
  #                        values = c("SkyBlue2", "black")) +
  geom_node_point(aes(size = V(cite_subgraph)$all_cites,
                      color = V(cite_subgraph)$group,
                      alpha = V(cite_subgraph)$literature)) +
  scale_size_continuous(name = "No. citations",
                        trans = "sqrt",
                        range = .3 * c(1, sqrt(max(V(cite_graph)$all_cites))),
                        breaks = seq(1, max(V(cite_graph)$all_cites),
                                     floor(max(V(cite_graph)$all_cites) /
                                             4))) +
  scale_color_brewer(name = "Group", type = "qual", palette = "Dark2") +
  scale_alpha_manual(name = "Literature", values = c(.5, 1)) +
  geom_node_text(aes(label = V(cite_subgraph)$label),
                 repel = TRUE) +
  guides(alpha = guide_legend(override.aes = list(size = 4)),
         color = guide_legend(override.aes = list(size = 4))) +
  ggforce::theme_no_axes() +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        legend.box = "horizontal")
pdf(height = bodwid * 5 * 1.1, width = bodwid * 5,
    file = "fig/cite_subplot_big.pdf")
print(cite_subplot_alt)
dev.off()
