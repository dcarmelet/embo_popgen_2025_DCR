
install.packages('tidypopgen', repos = c('https://evolecolgroup.r-universe.dev',
                                         'https://cloud.r-project.org'))
install.packages("admixtools", 
                 repos = c("https://evolecolgroup.r-universe.dev",
                           "https://cloud.r-project.org"))

library(tidypopgen)
library(admixtools)
library(ggplot2)

dir("./Me/Day5/data/")

modern_gt <- tidypopgen::gen_tibble("./Me/Day5/data//modern_samples.bed",
                                    valid_alleles = c("A","T","C","G"),
                                    missing_alleles = c("X"))

modern_gt

modern_gt %>% group_by(population) %>% tally()

loci_report <- modern_gt %>%
  group_by(population) %>%
  qc_report_loci()

autoplot(loci_report, type = "missing")

modern_gt <- modern_gt %>% 
  select_loci_if(loci_missingness(genotypes)<0.04)

ancient_gt <- tidypopgen::gen_tibble("./Me/Day5/data/ancient_samples.vcf",
                                          valid_alleles = c("A","T","C","G","X"), 
                                          quiet = TRUE)

ancient_gt <- gt_pseudohaploid(ancient_gt)

ancient_gt$id

ancient_gt$id[ancient_gt$id == "GB20"] <- "Mota"

ancient_gt$population <- ancient_gt$id

merged_dry <- rbind_dry_run(modern_gt, ancient_gt, 
                             flip_strand = TRUE)

merged_gt <- rbind(modern_gt, ancient_gt, 
                     flip_strand = TRUE,
                     backingfile = "./Me/Day5/data/merged_samples")

merged_gt <- merged_gt %>% group_by(population)

f2_dir <- "./Me/Day5/data/f2_tidypopgen"

# compute f2
f2_tidypopgen <- gt_extract_f2(merged_gt, 
                               outdir = "./Me/Day5/data/f2_tidypopgen", 
                               overwrite = TRUE)

f2_blocks = f2_from_precomp("./Me/Day5/data/f2_tidypopgen")

lbk_modern_panel <- c("Basque", "Bedouin2", "Druze", "Cypriot", "Tuscan",
  "Sardinian", "French", "Spanish", "Onge", "Han", "Mayan", "Mixe", "Surui")

lbk_f3out <- f3(data = f2_blocks, 
                pop1 = "Mbuti", 
                pop2 = "LBK",
                pop3 = lbk_modern_panel,blgsize=0.05)
lbk_f3out

lbk_f3out %>% arrange(desc(est))


ggplot(lbk_f3out, aes(pop3, est)) +
  geom_point() +
  geom_errorbar(aes(ymin = est - 2 * se, ymax = est + 2 * se)) +
  labs(y = "Shared drift with LBK", x = "populations") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


lbk_f3out$pop3<-factor(lbk_f3out$pop3, levels = lbk_f3out$pop3[order(lbk_f3out$est)])

ggplot(lbk_f3out, aes(pop3, est)) +
  geom_point() +
  geom_errorbar(aes(ymin = est - 2 * se, ymax = est + 2 * se)) +
  labs(y = "Shared drift with LBK", x = "populations") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


aa_f3admix <- f3(data = f2_blocks,
                 pop1 = "AA",
                 pop2 = "Yoruba",
                 pop3 = "French")
aa_f3admix

##QUESTION: Do we have evidence for admixture?

#Yes negative (-0.00472) and significante (Z = -39.4, p=0)

eurasian_sources <- c("French","Spanish","Sardinian","LBK")

somali_f3admix <- f3(data = f2_blocks,
                     pop1 = "Somali",
                     pop2 = eurasian_sources, 
                     pop3 = "Mota")
somali_f3admix

#QUESTION: What can you conclude? 
# f3(C;A,B) < 0 --> Admixture
#Those population show admixtures with Mota

#QUESTION: Now do the same for the Dinka. Do you get a similar result? If not, what can we conclude?

somali_f3admix2 <- f3(data = f2_blocks,
                     pop1 = "Somali",
                     pop2 = "Dinka", 
                     pop3 = "Mota")
somali_f3admix2

# No evidences of admixture, it could be a long drift of Somali or the Dinka and Mota being sister clades

neand_french_f4 <- f4(data = f2_blocks,
                     pop1 = "pan_troglodytes",
                     pop2 = "AltaiNea",
                     pop3 = "Mbuti",
                     pop4 = "French")
neand_french_f4

 f4(data = f2_blocks,
                      pop1 = "pan_troglodytes",
                      pop2 = "AltaiNea",
                      pop3 = "Mbuti",
                      pop4 = "Han")

#QUESTION: Now we want to investigate the extend of Yamnaya ancestry that came 
#into Europe during the Bronze age. If we assume that the LBK was the predominant 
#genetic ancestry in Europe before the Yamnaya arrived, test: 
#1. that “French”, “Basque” and “Spanish”  all harbour Yamnaya ancestry 


Yamnaya_f4 <- f4(data = f2_blocks,
                      pop1 = "pan_troglodytes",
                      pop2 = "Yamnaya",
                      pop3 = "LBK",
                      pop4 = c("French","Spanish","Basque","Mbuti"))

Yamnaya_f4
# Yes they all do !

#2. whether these populations differ in how much ancestry they share 
#Can you avoid writing down every comparison individually, and use vectors of populations to create 
#some simple code to ask each question (look at what we did with the outgroup f3).


pops <- c("Han", "pan_troglodytes", "AA","Yoruba","French")
qpf4ratio(data = f2_blocks, 
          pops = pops)

#QUESTION: What is the proportion of European contribution into African Americans in the dataset?

#  0.158 is alpha, the admixture  between French and Yoruba showing that the admixture of Yoruba lean more toward AA than French

#• SISTER: refers to a sister population, one is related more closely to your hypothesised source of admixture than to the ancestral population of our target (in this case Han) 
#• OUTGROUP: refers to an outgroup to all other populations (in this case Pan troglodytes) 
#• TARGET: the admixed population under investigation (in this case, African American) 
#• ANCESTRAL: refers to the ancestral population that we know to be related to our admixed target (in this case, Yoruba) 
#• HYPOTHESISED SOURCE OF ADMIXTURE: efers to the population who is suggested to have admixed into the target (in this case, French)

#QUESTION: Can you write some code to compute the f4 ratio by hand (i.e. estimate the numerator and denominator f4, and compute the ratio by hand)? 
#QUESTION: Repeat the analysis using other European populations (Sardinian, Spanish and Basque). Do the results change when using different European populations? 
#QUESTION: A final little challenge. qpf4ratio can take a matrix of populations, where each line gives the 5 pops for an f4 ratio. Can you write the code to run the 4 comparisons we just run as a single command?


lbk_modern_panel <- c("Basque", "Bedouin2", "Druze", "Cypriot", "Tuscan",
  "Sardinian", "French", "Spanish", "Onge", "Han", "Mayan", "Mixe", "Surui")
modern_panel_gt <- merged_gt %>% filter (population %in% lbk_modern_panel)
# remove monomorphic sites (ungrouping the tibble, as we want to use global frequencies
# rather than population frequencies)
modern_panel_gt <- modern_panel_gt %>% ungroup() %>% select_loci_if(loci_maf(genotypes)>0)

modern_panel_gt <- gt_update_backingfile(modern_panel_gt)
# reset the ploidy for this tibble, as it is now all diploids
attr(modern_panel_gt$genotypes, "ploidy") <- 2
modern_panel_gt <- modern_panel_gt %>% gt_impute_simple(method="mode")

modern_panel_gt <- modern_panel_gt %>% select_loci_if(loci_ld_clump(genotypes))

modern_pca <- modern_panel_gt %>% gt_pca_randomSVD()

autoplot(modern_pca, type = "scores") 

library(ggplot2)
autoplot(modern_pca, type = "scores") +
  aes(color = modern_panel_gt$population) +
  labs(color = "population")


lbk_gt <- merged_gt %>% filter(id == "LBK")
lbk_pca_scores <- predict(modern_pca, new_data = lbk_gt, project_method = "least_square")
lbk_pca_scores

autoplot(modern_pca, type = "scores") +
  aes(color = modern_panel_gt$population) +
  labs(color = "population")+
  geom_point(data=lbk_pca_scores, mapping=aes(x=.data$.PC1, y=.data$.PC2), col = "black")


autoplot(modern_pca, type = "scores") +
  aes(color = modern_panel_gt$population) +
  labs(color = "population") +
  geom_point(data=lbk_pca_scores, mapping=aes(x=.data$.PC1, y=.data$.PC2), col = "black") +
  lims(x=c(30, 70), y = c(-10, 15))

