library(igraph)

twitter <- read.csv('~/workspace/SNA-R/data/twitter.csv')
#colnames(friends) <- c('mutual', 'create_time', 'alter', 'ego')
#head(colnames)
#friends_data = friends[, c('ego','alter','mutual','create_time')]
#rm(friends)
twitter_full <- graph.data.frame(d = twitter)

write.graph(twitter_full, "~/workspace/Reachability/twitter.lgl", format=c("lgl"))

deg_twitter_in <- degree(twitter_full, mode="in") 
deg_twitter_out <- degree(twitter_full, mode="out") 

twitter_reachability=vector()
for (i in 1:vcount(twitter_full)){
  twitter_reachability[length(twitter_reachability)+1]=length(subcomponent(twitter_full, (i), mode = 'out'))
}

#hist(slashdot_reachability)
