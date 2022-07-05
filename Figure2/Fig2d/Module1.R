##laoding libraries
library(igraph)
library(networkD3)

###load Module 1 data
mod1 = read.table("Module1.csv",sep=",",header=T,stringsAsFactors = F)

##Plotting graph
p <- simpleNetwork(mod1)
g <- graph.data.frame(d = mod1, directed = FALSE)
myCols <- setNames(c("#DCE319FF","#DCE319FF","#55C667FF","#55C667FF","#55C667FF","#55C667FF","#55C667FF","#55C667FF","#55C667FF","#55C667FF","#55C667FF","#55C667FF","#55C667FF","#55C667FF","#55C667FF","#55C667FF"),
                   c("THAP1", "YY1","ADIPOR2","ADNP","BCKDHA","GALNT2","GPT2",   
                     "ISCU","MAST4","MON1A","PPIF","TBC1D22B",
                     "TTC31","TTLL7","UNC13B","WDR37"))
V(g)$color <- myCols[V(g)$name]
plot(g, vertex.label = V(g)$name,vertex.label.color="black")

