


HIVprev = .83*.17 + .17*.37
n_HIVpos = round(51010 * HIVprev)
n_HIVneg = round(51010 * (1-HIVprev))

node.df = data.frame(id = c(c(1:n_HIVpos), c((n_HIVpos+1):(n_HIVpos+n_HIVneg))),
                      HIVpos = c(rep(1,n_HIVpos), rep(0,n_HIVneg))
                      )

prop_pos = .31 #The probability that a connection from a HIV+ goes to HIV+; Prob(HIV+ recruitee | HIV+ recruiter)
prop_neg = .95 #The probability that a connection from a HIV- goes to HIV-; Prob(HIV- recruitee | HIV- recruiter)
mean_social_edges = 4.17934783  #2)	San Diego social network contact mean number

prop_overlap_genetic = .1 #The proportion of pos-pos edges in genetic that is part of the social network

alpha = n_HIVpos / 1262
genetic_cluster_sizes = c(rep(1,665*alpha),
                           rep(2,102*alpha),
                           rep(3,25*alpha),
                           rep(4,15*alpha),
                           rep(5,8*alpha),
                           rep(6,4*alpha),
                           rep(7,3*alpha),
                           rep(8,3*alpha),
                           rep(9,4*alpha),
                           rep(10,0*alpha),
                           rep(11,1*alpha),
                           rep(12,2*alpha),
                           rep(13,0*alpha),
                           rep(14,1*alpha),
                           rep(15,0*alpha),
                           rep(16,0*alpha),
                           rep(17,0*alpha),
                           rep(18,0*alpha),
                           rep(19,0*alpha),
                           rep(20,1*alpha),
                           rep(21,1*alpha),
                           rep(22,0*alpha),
                           rep(23,1*alpha)) #ensure that the sum adds to n_HIVpos, i.e., every HIV+ in cluster 

#RDS simulation

RDS_nsamp0 = 25 #the number of seeds to be drawn (i.e. the size of the 0th wave of sampling)
RDS_fixinitial	= NULL # a variable that indicates the distribution from which to draw the initial seeds, if the seeds variable is NULL and the seed.distribution variable is NULL
RDS_replace	= FALSE #sampling with replacement
RDS_coupons	= 10 #number of coupons

#Hypothesis test


