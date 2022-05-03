library(igraph)
library(networkD3)
mod2 = read.table("Module1.csv",sep=",",header=T,stringsAsFactors = F)
p <- simpleNetwork(mod2)

g <- graph.data.frame(d = mod2, directed = FALSE)

myCols <- setNames(c("#DCE319FF","#DCE319FF","#55C667FF","#55C667FF","#55C667FF","#55C667FF","#55C667FF","#55C667FF","#55C667FF","#55C667FF","#55C667FF","#55C667FF","#55C667FF","#55C667FF","#55C667FF","#55C667FF"),
                   c("THAP1", "YY1","ADIPOR2","ADNP","BCKDHA","GALNT2","GPT2",   
                     "ISCU","MAST4","MON1A","PPIF","TBC1D22B",
                     "TTC31","TTLL7","UNC13B","WDR37"))
V(g)$color <- myCols[V(g)$name]

plot(g, vertex.label = V(g)$name,vertex.label.color="black")

