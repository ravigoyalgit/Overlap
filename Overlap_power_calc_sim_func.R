
generate_genetic_network <- function() {
  
  g_genetic = make_empty_graph(n = n_HIVpos + n_HIVneg, directed = FALSE)
  
  genetic_cluster_sizes = sort(genetic_cluster_sizes, decreasing = TRUE)
  
  for (i in c(1:length(genetic_cluster_sizes))) {
    genetic_cluster_size = genetic_cluster_sizes[i]
    if (genetic_cluster_size > 1) {
      if (i > 1) {
        start_node = genetic_cluster_sizes[c(1:(i-1))] %>% sum + 1
      } else {
        start_node = 1
      }
      end_node = genetic_cluster_sizes[c(1:i)] %>% sum
      genetic_edges = tuples(start_node:end_node, 2)
      g_genetic = igraph::add_edges(g_genetic, c(t(genetic_edges)))
    }
  }
  
  g_genetic  = igraph::simplify(g_genetic)
  #components(g_genetic)$csize
  return(g_genetic)
}

generate_overlap_network <- function(g) {
  
  g_edges <- igraph::as_data_frame(g) %>% rename(node1 = from, node2 = to)
  g_edges <- left_join(g_edges, node.df, by = c("node1" = "id")) %>%
    rename(HIVpos_node1 = HIVpos)
  g_edges <- left_join(g_edges, node.df, by = c("node2" = "id")) %>%
    rename(HIVpos_node2 = HIVpos)
  
  g_edges_pp <-  filter(g_edges, HIVpos_node1 == 1 & HIVpos_node2 == 1) 
  g_overlap_edges <- sample_n(g_edges_pp, prop_overlap_genetic * (g_edges_pp %>% nrow()), replace = FALSE)
  
  g_overlap = make_empty_graph(n = n_HIVpos + n_HIVneg, directed = FALSE) %>% 
    add_edges(c(t(g_overlap_edges[,c("node1", "node2")])))
  
  return(g_overlap)
}

generate_social_network <- function(g_genetic, g_overlap) {
  
  node_HIVpos.df <- filter(node.df, HIVpos == 1)
  n_HIVpos  <- nrow(node_HIVpos.df)
  
  node_HIVneg.df <- filter(node.df, HIVpos == 0)
  n_HIVneg  <- nrow(node_HIVneg.df)
  
  p_PP <-  0 #(prop_pos * mean_social_edges * n_HIVpos)/choose(n_HIVpos,2)  #number of edges / possible
  p_PN <-  ((mean_social_edges * (n_HIVpos+n_HIVneg) * .5) - (prop_pos * mean_social_edges * n_HIVpos * .5) - (prop_neg * mean_social_edges * n_HIVneg * .5))/(n_HIVpos * n_HIVneg)
  p_NN <-  (prop_neg * mean_social_edges * n_HIVneg * .5)/choose(n_HIVneg,2)
  pm <- cbind( c(p_PP, p_PN), c(p_PN, p_NN) )
  g_social <- igraph::sample_sbm(n_HIVpos+n_HIVneg, pref.matrix=pm,
                                 block.sizes=c(n_HIVpos,n_HIVneg),
                                 directed = FALSE, loops = FALSE) %>%
    set_vertex_attr("name", value = c(node_HIVpos.df$id, node_HIVneg.df$id))
  
  #Need to add HIV+ to HIV+ edges
  
  g_social = igraph::add_edges(g_social, c(t(as_edgelist(g_overlap))))
  
  edges_pp_needed = round((prop_pos * mean_social_edges * n_HIVpos * .5) - g_overlap %>% igraph::ecount())
  for (i in c(1:edges_pp_needed)) {
    invalid_edge = TRUE
    while(invalid_edge) {
      potential_edge = sample(c(1:n_HIVpos), 2, replace = FALSE)
      if ((get.edge.ids(g_genetic, potential_edge) == 0) & (get.edge.ids(g_social, potential_edge) == 0)) {
        #Valid edge to add
        g_social = igraph::add_edges(g_social, potential_edge)
        invalid_edge = FALSE
      }
    }
  }
  
  return(g_social)
  
}

simulate_RDS <- function(g_social) {
  
  net_social = intergraph::asNetwork(g_social)
  net_social = network::set.vertex.attribute(net_social, attrname = "HIVpos", value = c(rep(1,n_HIVpos), rep(0,n_HIVneg)))
                                    
  RDS_results = RDS::rdssampleC(
    net =  net_social,
    nnodes = network::network.size(net_social),
    nsamp0 = RDS_nsamp0,
    fixinitial = FALSE,
    nsamp = RDS_nsamp,
    replace = RDS_replace,
    coupons = RDS_coupons,
    seed.distribution = NULL,
    attrall = FALSE,
    trait.variable = "HIVpos",
    nsims = 1,
    seeds = NULL,
    prob.network.recall = 1,
    verbose = TRUE
  )
  
  g_RDS = make_empty_graph(n = n_HIVpos + n_HIVneg, directed = TRUE) 
  
  for (i in c(1:length(RDS_results$nsample))) {
    node_sample = RDS_results$nsample[i]
    node_nominator = RDS_results$nominators[i]
    
    if (node_nominator != 0) {
      g_RDS = g_RDS %>% add_edges(c(node_nominator, node_sample))
    }
  }
  
  RDS_results$g_RDS = g_RDS
  
  g_RDS_UD = as.undirected(
    g_RDS,
    mode = "collapse"
  )
  
  RDS_results$g_RDS_overlap = igraph::intersection(g_genetic, g_RDS_UD)

  return(RDS_results)
}
