
library('igraph')   
library('tidyverse')
library(mpoly)
library('intergraph')
library('RDS')

source('Overlap_power_calc_sim_params.R')
source('Overlap_power_calc_sim_func.R')


g_genetic = generate_genetic_network()
g_overlap  = generate_overlap_network(g_genetic)
g_social = generate_social_network(g_genetic, g_overlap)

RDS_results = simulate_RDS(g_social)

