
library(igraph)
setwd("~/workspace/SNA-R")

#nodes_num=1000 #1000
#edges_num=nodes_num*8 #8*(nodes_num-1)
prob=0.5

#prob=0.5
#nodes_num=2000
#edges_num=16000
edges=c()
  edges=c(edges, 1,2)
  edges=c(edges, 1,3)
  edges=c(edges, 3,2)
  edges=c(edges, 2,4)
#  edges=c(edges, 2,5)
#  edges=c(edges, 4,5)
#  edges=c(edges, 5,3)
#  edges=c(edges, 5,1)
# edges=c(edges, 1,2)
# edges=c(edges, 2,3)
# edges=c(edges, 2,4)
# edges=c(edges, 4,1)
# edges=c(edges, 1,3)
# edges=c(edges, 3,4)
# edges=c(edges, 4,2)
# edges=c(edges, 4,3)
# edges=c(edges, 3,5)


g <- graph(edges, directed=TRUE)

#power law random graph
#powerG <- static.power.law.game(nodes_num, edges_num, 2.45, 2.45, loop=FALSE)
#randomG <- erdos.renyi.game(nodes_num, edges_num, type=c("gnm"), directed=T)

#g=read.graph("~/workspace/Reachability/output/2000.lgl", format=c("lgl"))
N=vcount(g)
E=ecount(g)
rSample=Matrix(0, nrow=N, ncol=N, sparse=TRUE)

cascade_size = function(g, starter, nodes){
  infected=c()
  g=set.edge.attribute(g, "weight", value=toss)
  g_mod = delete.edges(g, which(E(g)$weight == 0))
  
  sizes=c()
  # BFS on linear
  for (start in starter){
    #potential = neighborhood(linear, 1, infected)
    potential = subcomponent(g_mod, start, mode=c("out"))
    sizes=c(sizes, length(potential)-1)
    infected = unique(c(infected, potential))
  }
  
  for (i in infected){
    nodes[[i]] = nodes[[i]] + 1
  }
  
  return (nodes)
}

IS=list()
CS=list()
reachability=vector()
cat("start: ", date(),"\n")
for (initial in 1:1){#nodes_num){
  IS_v=vector()
  CS_v=vector()
  
  nodes=list()
  for (i in 1:N){
    nodes[[i]]=0
  }
  for (i in 1:10000){
    #get intial number of starter
    #starter=sample(1:nodes_num, initial, replace=F)
    starter=initial
    infected=c()
    toss=sample(c(0,1),edges_num,replace=TRUE,prob=c(1-prob, prob))
    
    nodes=cascade_size(g, starter, nodes)
    
    #sizes = cascade_size(g, starter)
    #IS_v[length(IS_v)+1] = sizes[length(sizes)]
    #CS_v = c(CS_v, sizes[1:length(sizes)-1])    
    
  }
  
  #cat(initial, date(),"\n")
  for (i in 1:N){
    cat(i, nodes[[i]]/10000, "\n")
  }

  
#   IS[[initial]]=IS_v
#   CS[[initial]]=CS_v
#   
#   reachability[length(reachability)+1] = mean(IS_v)
}

#cat(reachability)

#hist(reachability)


