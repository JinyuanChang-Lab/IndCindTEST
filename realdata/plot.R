rm(list=ls()) 
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
library(RColorBrewer)
library(ggpubr)



x <- read.csv('D:/0-xnc-du-dr/1-paper/0-code-new/new_FNN/real_red1.csv')
y <- as.matrix(x, 12)
data_1 <- y[1:11, 2:12]
data_2 <- y[1:11, 14:24]



colnames(data_2) <- colnames(data_1)
rownames(data_1) <- rownames(data_2) <- colnames(data_1)




colnames(data_1) <- colnames(data_2) <- c("Communication Services","Consumer Discretionary", "Consumer Staples","Engery","Financials",
                                          "Health Care", "Industrials", " Information Technology", "Materials","Real Estate","Utilities")

#c("Engery", "Materials", "Industry", "Consumer Discretionary", 
#                                          "Consumer Staples", "Health Care", "Financials", " Information Technolohy",
#                                          "Communication Services",  "Utilities", "Real Estate")

rownames(data_1) <- rownames(data_2)  <-  c("Communication Services","Consumer Discretionary", "Consumer Staples","Engery","Financials",
                                            "Health Care", "Industrials", " Information Technology", "Materials","Real Estate","Utilities")



###graph, V()-node and E()-edge 
g1 <- graph_from_adjacency_matrix(data_1, mode =  "undirected")
g2 <- graph_from_adjacency_matrix(data_2, mode =  "undirected")


#vcolor <- brewer.pal(9, "Spectral")
vcolor <- brewer.pal(11, "Spectral")
#vcolor <- brewer.pal(9, "Pastel1")

graph_gt_1 <- as_tbl_graph(g1) %>% 
  mutate(deg = centrality_degree(mode='in'))
deg <- NULL
degree <- NULL
deg <- degree(g1,mode="all")
degree <- factor(deg) 

###colour="gray10"/degree
p_1 <- ggraph(graph_gt_1) +
  geom_edge_link() +
  geom_node_point(aes(colour = as.factor(colnames(data_1)), alpha=1, shape = 'circle'), size = 4.5*deg,
                  show.legend= F) +
  geom_node_text(aes(label = rownames(data_1)), colour = "black", size=6,  repel = TRUE) + 
  geom_edge_link(colour = "gray50") + 
  guides(shape='none', color = guide_legend(override.aes = list(size = 10))) +
  theme_graph()
p_1 + expand_limits()

###change color
p_12 <- p_1 + scale_fill_manual(values=rev(vcolor))+
  scale_color_manual(values=rev(vcolor)) 
p_12





#####second
graph_gt_2 <- as_tbl_graph(g2) %>% 
  mutate(deg = centrality_degree(mode='in'))
deg <- NULL
degree <- NULL
deg <- degree(g2,mode="all")
degree <- factor(deg) 

p_2 <- ggraph(graph_gt_2,layout = 'kk') +
  geom_edge_link() +
  geom_node_point(aes(colour = as.factor(colnames(data_1)), alpha=1, shape = 'circle'), size = 4.5*deg,
                  show.legend= F) +
  geom_node_text(aes(label = rownames(data_2)), colour = "black", size=6,  repel = TRUE) + 
  geom_edge_link(colour = "gray50") + 
  guides(shape='none', color = guide_legend(override.aes = list(size = 10))) +
  theme_graph()
p_2

###change color
p_22 <- p_2 + scale_fill_manual(values=rev(vcolor))+
  scale_color_manual(values=rev(vcolor))
p_22

