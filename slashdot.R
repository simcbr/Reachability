library(igraph)

slashdot = read.table("~/workspace/SNA-R/data/soc-sign-Slashdot081106.txt", 
                    sep="\t", 
                    col.names=c("ego", "alter","sign"), 
                    fill=FALSE, 
                    strip.white=TRUE)
slashdot_full <- graph.data.frame(d = slashdot)

write.graph(slashdot_full, "~/workspace/Reachability/slashdot.lgl", format=c("lgl"))

#node statistic
deg_slashdot_in <- degree(slashdot_full, mode="in") 
deg_slashdot_out <- degree(slashdot_full, mode="out") 


slashdot_reachability=vector()
for (i in 1:vcount(slashdot_full)){
  slashdot_reachability[length(slashdot_reachability)+1]=length(subcomponent(slashdot_full, (i), mode = 'out'))
}

hist(slashdot_reachability)

