library(admixtools)
library(tidypopgen)

setwd("../")
f2_blocks = f2_from_precomp("./data/f2_tidypopgen", verbose = FALSE)

neand_euras_wave <- qpwave(data = f2_blocks,
      left = c("French","Spanish","Tuscan"),
      right = c("AltaiNea","Mota", "Yoruba", "Denisova", "Mbuti")
)
neand_euras_wave

#p_nested is not significant so we rank 1 is not explaining significantly more
# rank +1 -> 1 wave

#QUESTION: Now, add the ‘Han’ and ‘Onge’ to the set of left populations and repeat the test. Does the result change? What does it show?

neand_euras_wave <- qpwave(data = f2_blocks,
                           left = c("French","Spanish","Tuscan","Han","Onge"),
                           right = c("AltaiNea","Mota", "Yoruba", "Denisova", "Mbuti")
)
neand_euras_wave


# 2 waves because we reject rank 0

french_adm <- qpadm(data = f2_blocks,
      left = c("Loschbour", "LBK", "Yamnaya"),
      right = c("Mbuti", "Mota", "Dinka", "Yoruba", "Han"),
      target= "French")

french_adm$popdrop

#QUESTION: Now, let’s look at the Basque. Can they be modelled in the same way? 
#How large is the Yamnaya component compared to the French?

basq_adm <- qpadm(data = f2_blocks,
                    left = c("Loschbour", "LBK", "Yamnaya"),
                    right = c("Mbuti", "Mota", "Dinka", "Yoruba", "Han"),
                    target= "Basque")

basq_adm$popdrop

# 001 and 100 are working models to explain Basque admixture 65% percent Yamnaya instead of 75%

base_edges <- matrix(
  c("R",	"Mbuti",
    "R", "eAfr",
    "eAfr",	"Dinka",
    "eAfr",	"outAfrica",
    "outAfrica",	"Han",
    "outAfrica",	"Loschbour"),
  ncol=2,byrow = TRUE,
  dimnames=list(NULL, c("from","to")))

base_edges %>% edges_to_igraph() %>% plot_graph()

base_edges <- matrix(
  c("R",	"Mbuti",
    "R", "eAfr",
    "eAfr",	"Dinka",
    "eAfr",	"outAfrica",
    "outAfrica",	"Han",
    "outAfrica",	"Loschbour"),
  ncol=2,
  byrow = TRUE,
  dimnames=list(NULL, c("from","to")))

base_edges

base_igraph <- base_edges %>% edges_to_igraph()

is_valid(base_igraph)

base_igraph %>% plot_graph()

base_igraph %>% plotly_graph()

base_qpgraph <- qpgraph(data = f2_blocks, graph = base_igraph)

base_qpgraph$f3

base_qpgraph$f3 %>% filter(abs(z)>2)

base_qpgraph$edges %>% plot_graph()

#QUESTION: Let’s see what happens if we choose the wrong outgroup. 
#Try creating a graph where we swap the role of Mbuti and Dinka 
#(i.e. placing Mbuti as the closest African population to European and Asians). 
#Use base_swapped_qpgraph as the name for the graph you create.



base_edges <- matrix(
  c("R",	"Dinka",
    "R", "eAfr",
    "eAfr",	"Mbuti",
    "eAfr",	"outAfrica",
    "outAfrica",	"Han",
    "outAfrica",	"Loschbour"),
  ncol=2,
  byrow = TRUE,
  dimnames=list(NULL, c("from","to")))

base_edges

base_igraph <- base_edges %>% edges_to_igraph()

is_valid(base_igraph)

base_igraph %>% plot_graph()


base_swapped_qpgraph <- qpgraph(data = f2_blocks, graph = base_igraph)

base_swapped_qpgraph$f3

fits = qpgraph_resample_multi(f2_blocks, 
                              graphlist = list(base_qpgraph[[1]], base_swapped_qpgraph[[1]]), 
                              nboot = 100)
compare_fits(fits[[1]]$score_test, fits[[2]]$score_test)


base_igraph %>% plot_graph(highlight_unidentifiable = TRUE)

yamnaya_edges <- matrix(
  c(
    "R",	"Mbuti",
    "R",	"eAfr",
    "eAfr", "Dinka",
    "eAfr", "outAfrica",
    "outAfrica", "Han",
    "outAfrica", "wEurasian",
    "wEurasian", "Yamnaya",
    "wEurasian", "Loschbour"),
  ncol = 2,
  byrow = TRUE,
  dimnames = list(NULL, c("from", "to")))
yamnaya_igraph <- yamnaya_edges %>% edges_to_igraph()
yamnaya_igraph %>% plot_graph()

#QUESTION: Add these edges and create an object names yamnaya_igraph, 
#and plot it to make sure that it matches our desired topology. 
#After checking that the graph is valid, fit the graph, inspect it, 
#and check whether it is compatible with the data: And ask whether our graph 
#is compatible with the data Does any population have a z value larger than 3? 
#What do you conclude (remember to plot and inspect the fitted graph)?

yamnaya_qpgraph <- qpgraph(data = f2_blocks, graph = yamnaya_igraph)
yamnaya_qpgraph$f3
yamnaya_qpgraph$edges %>% plot_graph()

#Compatible 
#QUESTION: Create a set of edges lbk_edges to match the graph above, 
#and then fit the graph to the data. Is the model compatible with the data?

lbk_extra_edges <- matrix(
  c(
    "R",	"Mbuti",
    "R",	"eAfr",
    "eAfr", "pBasalEurasian",
    "eAfr", "Dinka",
    "pBasalEurasian", "BasalEurasian",
    "pBasalEurasian","outAfrica",
    "outAfrica", "Han",
    "outAfrica","wEurasian",
    "wEurasian", "Yamnaya",
    "wEurasian", "pLoschbour",
    "pLoschbour", "Loschbour",
    "pLoschbour","WHG",
    "BasalEurasian", "pLBK",
    "WHG", "pLBK",
    "pLBK","LBK"),
  ncol = 2,
  byrow = TRUE,
  dimnames = list(NULL, c("from", "to")))
lbk_extra_igraph <- lbk_extra_edges %>% edges_to_igraph()
lbk_extra_igraph %>% plot_graph()

is_valid(lbk_extra_igraph)

lbk_extra_igraph %>% plot_graph(highlight_unidentifiable = TRUE)

lbk_extra_qpgraph <- qpgraph(data = f2_blocks, graph = lbk_extra_igraph)
lbk_extra_qpgraph$edges %>% plot_graph()





#QUESTION: If you have additional time, can you add an addition admixture edge, 
#exploring what happens if you add the Sardinians. Try modelling the Sardinian as a 
#direct descendant of the ancestor of LBK, and then try to add some admixture from the 
#Yamnaya (as we saw earlier that we could detect that signal). 
#Can you get away with modelling just a Neolithic component, 
#or do you need the Yamnaya admixture to make the model work. 
#If you have clear hypotheses, it is best to build alternative graphs and compare them. 
#On the other hand, if you are trying to place populations without a clear framework, 
#you would be better off automatically exploring many graphs. 
#As that exploration is computationally intensive, we will not do that during 
#this practical, but you have the building blocks on how to compare models, 
#so you would just compare a lot of different models rather than 2.


Sard_extra_edges <- matrix(
  c(
    "R",	"Mbuti",
    "R",	"eAfr",
    "eAfr", "pBasalEurasian",
    "eAfr", "Dinka",
    "pBasalEurasian", "BasalEurasian",
    "pBasalEurasian","outAfrica",
    "outAfrica", "Han",
    "outAfrica","wEurasian",
    "wEurasian", "Yamnaya",
    "wEurasian", "pLoschbour",
    "pLoschbour", "Loschbour",
    "pLoschbour","WHG",
    "BasalEurasian", "pLBK",
    "WHG", "pLBK",
    "pLBK","LBK",
    "LBK","Sardinian"),
  ncol = 2,
  byrow = TRUE,
  dimnames = list(NULL, c("from", "to")))
Sard_extra_igraph <- Sard_extra_edges %>% edges_to_igraph()
Sard_extra_igraph %>% plot_graph()

is_valid(Sard_extra_igraph)

Sard_extra_igraph %>% plot_graph(highlight_unidentifiable = TRUE)

Sard_extra_qpgraph <- qpgraph(data = f2_blocks, graph = Sard_extra_igraph)
Sard_extra_qpgraph$edges %>% plot_graph()





Sard_extra_edges_2 <- matrix(
  c(
    "R",	"Mbuti",
    "R",	"eAfr",
    "eAfr", "pBasalEurasian",
    "eAfr", "Dinka",
    "pBasalEurasian", "BasalEurasian",
    "pBasalEurasian","outAfrica",
    "outAfrica", "Han",
    "outAfrica","wEurasian",
    "wEurasian", "Yamnaya",
    "wEurasian", "pLoschbour",
    "pLoschbour", "Loschbour",
    "pLoschbour","WHG",
    "BasalEurasian", "pLBK",
    "WHG", "pLBK",
    "pLBK","LBK",
    "LBK","pSardinian",
    "Yamnaya","pSardinian",
    "pSardinian","Sardinian"),
  ncol = 2,
  byrow = TRUE,
  dimnames = list(NULL, c("from", "to")))
Sard_extra_igraph_2 <- Sard_extra_edges_2 %>% edges_to_igraph()
Sard_extra_igraph_2 %>% plot_graph()

is_valid(Sard_extra_igraph_2)

Sard_extra_igraph_2 %>% plot_graph(highlight_unidentifiable = TRUE)

Sard_extra_qpgraph_2 <- qpgraph(data = f2_blocks, graph = Sard_extra_igraph_2)
Sard_extra_qpgraph_2$edges %>% plot_graph()

fits = qpgraph_resample_multi(f2_blocks, 
                              graphlist = list(Sard_extra_qpgraph[[1]], Sard_extra_qpgraph_2[[1]]), 
                              nboot = 100)
compare_fits(fits[[1]]$score_test, fits[[2]]$score_test)



