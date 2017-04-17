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

## new collaborations

# list of authors, sorted by debut
credits <- strsplit(incl_dat$authors, "\\|")
names(credits) <- incl_dat$google_scholar_cluster
credits <- credits[order(incl_dat$debut)]
authors <- unique(unlist(credits))
# empty graph of authors
author_graph <- make_empty_graph(n = length(authors), directed = FALSE)
V(author_graph)$name <- authors
prev_credit <- c()
new_collab <- rep(FALSE, length(credits))
# for each article, test for redundancy involving non-isolates and add new edges
for (i in 1:length(credits)) {
  # which authors are isolated
  #wh_isolated <- which(degree(author_graph,
  #                            V(author_graph)[credits[[i]]]) == 0)
  wh_prev <- which(credits[[i]] %in% prev_credit)
  # check whether any pair of authors is not already connected
  any_disconnected <- if(length(wh_prev) == 0) FALSE else {
    ecount(induced_subgraph(
      author_graph,
      V(author_graph)[credits[[i]]][wh_prev]
    )) < choose(length(wh_prev), 2)
  }
  if (any_disconnected) {
    plot(induced_subgraph(author_graph, V(author_graph)[credits[[i]]][wh_prev]))
  }
  # add edges and credit
  author_graph <- add_edges(author_graph,
                            as.vector(combn(V(author_graph)[credits[[i]]], 2)))
  prev_credit <- unique(c(prev_credit, credits[[i]]))
  # flag new collaborations
  if (any_disconnected) new_collab[i] <- TRUE
}
print(credits[which(new_collab)])

## coauthorship graph

# nodes (distinct authors)
coauthor_freq <- table(unlist(strsplit(incl_dat$authors, split = "\\|")))
coauthor_nodes <- names(coauthor_freq)
coauthor_indices <- 1:length(coauthor_nodes)
names(coauthor_indices) <- coauthor_nodes
# undirected edges (within-article pairs of authors)
coauthor_edges <- do.call(
  c,
  lapply(strsplit(incl_dat$authors, split = "\\|"), function(x) {
    if (length(x) < 2) return(NULL)
    y <- as.vector(combn(x, 2))
    unname(coauthor_indices[coauthor_indices[y]])
  })
)
# construct graph
coauthor_graph <- make_empty_graph(n = length(coauthor_nodes), directed = FALSE)
V(coauthor_graph)$name <- coauthor_nodes
coauthor_graph <- add_edges(graph = coauthor_graph, edges = coauthor_edges)

# aesthetics
V(coauthor_graph)$coauthor_freq <- coauthor_freq
V(coauthor_graph)$label <- coauthor_nodes

# subset to collaborating authors
coauthor_subgraph <- delete_vertices(coauthor_graph,
                                     which(degree(coauthor_graph) == 0))
# plot
coauthor_plot <- ggraph(coauthor_subgraph, layout = "fr", niter = 200) +
  geom_edge_fan(edge_width = 1/3) +
  geom_node_point(aes(size = V(coauthor_subgraph)$coauthor_freq)) +
  scale_size_continuous(name = "No. coauthors",
                        trans = "sqrt",
                        range = 2 * c(1, sqrt(max(coauthor_freq))),
                        breaks = seq(1, max(coauthor_freq), 3)) +
  geom_node_text(aes(label = ifelse(V(coauthor_subgraph)$coauthor_freq > 2,
                                    V(coauthor_subgraph)$label, ""))) +
  ggforce::theme_no_axes()
print(coauthor_plot)

## authorship graph

# nodes (article identfiers and author names)
cite_nodes <- incl_dat$google_scholar_cluster
author_nodes <- c(cite_nodes, coauthor_nodes)
# edges (authorships)
stopifnot(nrow(incl_dat) == length(cite_nodes))
author_edges <- do.call(c, lapply(1:nrow(incl_dat), function(i) {
  rbind(i, coauthor_indices[coauthor_indices[strsplit(incl_dat$authors[i],
                                                      split = "\\|")[[1]]]] +
          length(cite_nodes))
}))
# construct graph
author_graph <- make_empty_graph(n = length(author_nodes), directed = FALSE)
V(author_graph)$name <- author_nodes
author_graph <- add_edges(graph = author_graph, edges = author_edges)

# bipartite
V(author_graph)$type <- c(rep(FALSE, length(cite_nodes)),
                          rep(TRUE, length(coauthor_nodes)))
# citation counts
citations <- strsplit(incl_dat$internal_citations, split = "\\|")
in_cites <- ifelse(is.na(citations), 0, sapply(citations, length))
citings <- table(unlist(citations))
out_cites <- as.vector(citings)[match(cite_nodes, names(citings))]
out_cites[is.na(out_cites)] <- 0
# number of citations
V(author_graph)$cites <- c(out_cites, rep(NA, length(coauthor_nodes)))
V(author_graph)$cited <- c(in_cites, rep(NA, length(coauthor_nodes)))
#V(author_graph)$other_cites <- c()

# aesthetics
V(author_graph)$freq <- c(in_cites, coauthor_freq)
V(author_graph)$label <- c(rep(NA, length(cite_nodes)), coauthor_nodes)
V(author_graph)$label <- ""

# number of citations for each author
author_cites <- sapply(
  neighborhood(author_graph, 1, coauthor_nodes),
  function(x) {
    wh_authored <- match(names(x[-1]), cite_nodes)
    sum(in_cites[wh_authored])
  }
)
# include labels for authors with multiple papers, collaborators, & citations
V(author_graph)$label <- c(
  rep(NA, length(cite_nodes)),
  ifelse(degree(author_graph)[V(author_graph)$type] > 1 &
           coauthor_freq > 1 & author_cites > 5, coauthor_nodes, NA)
)

# subset to components with multiple articles
# connected components and the number of articles in each
clust <- clusters(author_graph)
clust_tab <- table(clust$membership[1:length(cite_nodes)])
# identify components with only one article
clust_single <- as.numeric(names(clust_tab)[clust_tab <= 1])
# remove nodes in components with only one article
author_multiple <- delete_vertices(
  author_graph, v = which(clust$membership %in% clust_single)
)

# flag likely publications reporting the same study (connected components)
author_comps <- components(author_graph)
for (i in 1:author_comps$no) {
  wh <- which(author_comps$membership == i)
  wh <- intersect(wh, which(V(author_graph)$type == 0))
  gscs <- V(author_graph)$name[wh]
  wh_dat <- which(incl_dat$google_scholar_cluster %in% gscs)
  if (length(wh_dat) < 2) next
  print(cbind(incl_dat[wh_dat, c("google_scholar_cluster", "year")],
              lead_author = gsub("^([^\\|]+)\\|.*$", "\\1",
                                 incl_dat$authors[wh_dat])))
}

# subset to citing and cited articles
author_connected <- delete_vertices(
  author_multiple, v = which(V(author_multiple)$cites == 0 &
                               V(author_multiple)$cited == 0)
)
# remove authors with no remaining articles
author_connected <- delete_vertices(
  author_connected, v = which(degree(author_connected) == 0)
)

# subset (again) to components with multiple articles
# connected components and the number of articles in each
clust <- clusters(author_connected)
clust_tab <- table(clust$membership[
  !is.na(as.numeric(V(author_connected)$name))])
# identify components with only one article
clust_single <- as.numeric(names(clust_tab)[clust_tab <= 1])
# remove nodes in components with only one article
author_multi_conn <- delete_vertices(
  author_connected, v = which(clust$membership %in% clust_single)
)

# identify most-cited publications with no already-included neighbors

# ggraph
set.seed(10) # 7, 10, 
author_plot <- ggraph(author_multi_conn, layout = "fr") +
  geom_edge_link(edge_width = 1/3) +
  geom_node_point(aes(size = V(author_multi_conn)$freq,
                      color = V(author_multi_conn)$type,
                      shape = V(author_multi_conn)$type,
                      alpha = !is.na(V(author_multi_conn)$label))) +
  scale_size_continuous(name = "No. citations/coauthors",
                        trans = "sqrt",
                        range = c(1, sqrt(max(V(author_multi_conn)$freq))),
                        breaks = c(1, max(V(author_multi_conn)$freq))) +
  scale_color_manual(name = "Node type", labels = c("Article", "Author"),
                     values = c("#e41a1c", "#377eb8")) +
  scale_shape_manual(name = "Node type", labels = c("Article", "Author"),
                     values = c(15, 16)) +
  scale_alpha_discrete(guide = FALSE, range = c(.6, 1)) +
  geom_node_text(aes(label = V(author_multi_conn)$label),
                 size = 3, repel = TRUE) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  ggforce::theme_no_axes()
print(author_plot)
pdf(height = bodwid * .85, width = bodwid, file = "fig/author_plot.pdf")
print(author_plot)
dev.off()

# add citation edges

# define edge type (citation or authorship)

# restrict to citations between (not within) components

# ggraph


## components

comp_tab <- sapply(bipartite_projection(author_graph),
                   function(g) components(g)$csize)
comp_tab[comp_tab[, 1] > 1 & comp_tab[, 2] > 1, ]
