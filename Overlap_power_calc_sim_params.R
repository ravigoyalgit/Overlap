


n_HIVpos = 500
n_HIVneg = 500

node.df = data.frame(id = c(c(1:n_HIVpos), c((n_HIVpos+1):(n_HIVpos+n_HIVneg))),
                      HIVpos = c(rep(1,n_HIVpos), rep(0,n_HIVneg))
                      )

prop_pos = .7 #The probability that a connection from a HIV+ goes to HIV+; Prob(HIV+ recruitee | HIV+ recruiter)
prop_neg = .7 #The probability that a connection from a HIV- goes to HIV-; Prob(HIV- recruitee | HIV- recruiter)
mean_social_edges = 4.17934783  #2)	San Diego social network contact mean number

prop_overlap_genetic = .1 #The proportion of pos-pos edges in genetic that is part of the social network

genetic_cluster_sizes = c(rep(1,463),
                           rep(2,0),
                           rep(3,0),
                           rep(4,0),
                           rep(5,0),
                           rep(6,0),
                           rep(7,0),
                           rep(8,0),
                           rep(9,0),
                           rep(10,0),
                           rep(11,0),
                           rep(12,0),
                           rep(13,0),
                           rep(14,1),
                           rep(15,0),
                           rep(16,0),
                           rep(17,0),
                           rep(18,0),
                           rep(19,0),
                           rep(20,0),
                           rep(21,0),
                           rep(22,0),
                           rep(23,1)) #ensure that the sum adds to n_HIVpos, i.e., every HIV+ in cluster 

#RDS simulation

RDS_nsamp0 = 4 #the number of seeds to be drawn (i.e. the size of the 0th wave of sampling)
RDS_fixinitial	= NULL # a variable that indicates the distribution from which to draw the initial seeds, if the seeds variable is NULL and the seed.distribution variable is NULL
RDS_nsamp	 = 20 #number of individuals in each RDS sample
RDS_replace	= FALSE #sampling with replacement
RDS_coupons	= 3 #number of coupons

