# enplots hallmark significative
library(clusterProfiler)
library(tidyverse)
library(dplyr)
library(msigdbr)
library(enrichplot)
library(DOSE)
library(ggplot2)

gseaplot2(em, geneSetID = 1, title = em$Description[1])
ggsave("H_APICAL_JUNCTION_enplot.pdf")

gseaplot2(em, geneSetID = 2, title = em$Description[2])
ggsave("H_EPITHELIAL_MESENCHYMAL_TRANSITION_enplot.pdf")

gseaplot2(em, geneSetID = 3, title = em$Description[3])
ggsave("H_COAGULATION_enplot.pdf")

gseaplot2(em, geneSetID = 4, title = em$Description[4])
ggsave("H_INFLAMMATORY_RESPONSE_enplot.pdf")

gseaplot2(em, geneSetID = 5, title = em$Description[5])
ggsave("H_FATTY_ACID_METABOLISM_enplot.pdf")

gseaplot2(em, geneSetID = 6, title = em$Description[6])
ggsave("H_OXIDATIVE_PHOSPHORYLATION_enplot.pdf")

gseaplot2(em, geneSetID = 1:3)
ggsave("dependent_enrich_enplot.pdf")

gseaplot2(em, geneSetID = 4:6)
ggsave("independent_enrich_enplot.pdf", width = 28, height = 23, units = "cm")
