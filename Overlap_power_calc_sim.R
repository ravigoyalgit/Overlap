
library('igraph')   
library('tidyverse')
library(mpoly)
library('intergraph')
library('RDS')

source('Overlap_power_calc_sim_params.R')
source('Overlap_power_calc_sim_func.R')

RDS_nsamp	 = 200 #SAMPLE SIZE - NEED TO SET

power = 0
for (i in c(1:100)) {
  g_genetic = generate_genetic_network()
  g_overlap  = generate_overlap_network(g_genetic)
  g_social = generate_social_network(g_genetic, g_overlap)
  
  RDS_results = simulate_RDS(g_social)
  
  expected_overlap = (igraph::ecount(g_genetic) / choose(n_HIVpos, 2)) #* (prop_pos * mean_social_edges * n_HIVpos * .5)
  obs_overlap = RDS_results$g_RDS_overlap %>% igraph::ecount()
  n_trials = igraph::induced_subgraph(RDS_results$g_RDS, v = c(1:n_HIVpos)) %>% igraph::ecount()
  
  stat_test = binom.test(obs_overlap, n_trials, p = expected_overlap) # / n_trials)
  power = power + as.numeric(stat_test$p.value < .05)
  
  print(i)
}

print(c("The power is: ", power))
