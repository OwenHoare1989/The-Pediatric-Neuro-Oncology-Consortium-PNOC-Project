PiePlotData <- read.csv("scRNA_celltype_fraction_per_treatment_treated.csv", header = T)
PiePlotData$label <- factor(PiePlotData$label, levels = c("Astrocytes",        "B cells",           "CD4+_T_cells",      "Endothelial cells", "Epithelial cells",  "M0_Macrophages",    "M1_Macrophages",   
                                                          "M2_Macrophages",    "Malignant",         "Microglia",         "Monocytes",         "NK_cells",          "Neurons",           "Oligodendrocytes" ))


library(dplyr)
library(ggplot2)

# calculate average across biological replicates for plotting pie chart
PiePlotData <- PiePlotData %>% group_by(label) %>% summarize(percent=mean(percent))

            

ggplot(PiePlotData, aes(x = factor(1), y=percent,fill=factor(label))) + geom_bar(width = 1,stat="identity")+coord_polar(theta = "y") +
  scale_fill_manual(labels = c("Malignant", "Neurons", "Endothelial", "Monocytes", "Microglia", "M2_Macrophages", "M1_Macrophages", "M0_Macrophages", "CD4", "B", "Astrocytes", "Oligodendrocytes", "NK_cells", "Epithelial"), values =c( "#CC313D", "#3B9AB2", "mediumseagreen", "limegreen", "blue", "violetred4",  "violetred",  "violet",  "#FF2700FF", "#FFBF47",  "#0062B4FF",  "#BF812D", "deeppink", "sienna4")) + 
  theme_void()+theme(legend.title = element_blank()) 

library(plotrix)



cellType = c(Malignant = "#CC313D", Neurons = "#3B9AB2", Endothelial =  "mediumseagreen", Monocytes = "limegreen", Microglia = "blue", M2_Macrophages = "violetred4", M1_Macrophages = "violetred",  M0_Macrophages = "violet",  CD4 = "#FF2700FF", B = "#FFBF47",  Astrocytes = "#0062B4FF", Oligodendrocytes = "#BF812D", NK_cells = "deeppink", Epithelial = "sienna4")
cellType2 = c(Malignant = "#CC313D", Neurons = "#3B9AB2", Endothelial =  "mediumseagreen", Monocytes = "limegreen", Microglia = "blue", M2_Macrophages = "violetred4", M1_Macrophages = "violetred",  M0_Macrophages = "violet",  CD4 = "#FF2700FF", B = "#FFBF47", Oligodendrocytes = "#BF812D", NK_cells = "deeppink", Epithelial = "sienna4")
# Explode the pie chart
library(plotrix)
library(ggpie)
PT2 <- subset(PT, subset=variable %in% c("Malignant", "Neurons", "Endothelial", "Monocytes", "Microglia", "M2_Macrophages", "M1_Macrophages", "M0_Macrophages", "CD4", "B", "Oligodendrocytes", "NK_cells", "Epithelial"))
pie3D(PT2$value, labels = PT2$variable, main = "PT 3D pie chart", explode=0.1, radius=2, col = cellType2)




library(rAmCharts4)
amPieChart(
  data = PT,
  category = "variable",
  value    = "value",
  threeD = TRUE,
  variableDepth = FALSE
)

library(magrittr) # needs to be run every time you start R and want to use %>%
#library(dplyr) 




library(ggplot2)
library(ggrepel)
#library(plyr)
library(forcats)
library(ggthreed)
library(threed)

suppressPackageStartupMessages({
  library(rgl)
  library(devout)
  library(devoutrgl)
  library(ggrgl)
  library(ggplot2)
})


PiePlotData %>% 
  mutate(csum = rev(cumsum(rev(percent))), 
         pos = percent/2 + lead(csum, 1),
         pos = if_else(is.na(pos), percent/2, pos),
         percentage = percent/sum(percent)) %>% 
  ggplot(aes(x = "", y = percent, fill = fct_inorder(label))) + 
  geom_col(width = 1, color = 1) +
  scale_fill_manual(labels = c("Astrocytes",        "B cells",           "CD4+_T_cells",      "Endothelial cells", "Epithelial cells",  "M0_Macrophages",    "M1_Macrophages",   
                               "M2_Macrophages",    "Malignant",         "Microglia",         "Monocytes",      "Neurons",    "NK_cells",                    "Oligodendrocytes"), values =c( "#0062B4FF", "#FFBF47", "#FF2700FF", "mediumseagreen",  "sienna4", "violet", "violetred", "violetred4", "#CC313D", "blue", "limegreen", "#3B9AB2",  "deeppink",  "#BF812D")) + 
  geom_label_repel(aes(y = pos,
                       label = glue::glue("{percent} ({scales::percent(percentage)})"), 
                       fill = label),
                   size = 3,
                   nudge_x = 0.5,
                   show.legend = FALSE) +
  labs(  fill = "label" ) +
  coord_polar(theta = "y") +
  theme_void() + 
  ggtitle("Cell type scRNA-seq treated") + theme(
  plot.title = element_text(color="black", size=14),
  axis.title.x = element_text(color="black", size=14),
  axis.title.y = element_text(color="black", size=14))




PiePlotData$Treatment <- factor(PiePlotData$Treatment, levels = c("Control", "Treated"))
# Visualize: Specify the comparisons you want
my_comparisons = list(c("Control", "Treated"))
bp = ggboxplot(PiePlotData, x = "label", y = "Percentage", add = "jitter", legend = "none", combine = TRUE,
          color = "label", palette = c("#0062B4FF", "#FFBF47", "#FF2700FF", "mediumseagreen",  "sienna4", "violet", "violetred", "violetred4", "#CC313D", "blue", "limegreen", "#3B9AB2","deeppink", "#BF812D"),ggtheme = theme_bw())+ 
  rotate_x_text(angle = 45)+
  geom_hline(yintercept = mean(PiePlotData$Percentage), linetype = 2)+
  stat_compare_means(comparisons = my_comparisons, paired = FALSE, method = "t.test")+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 100)  +
  ggtitle("Test") + theme(
    plot.title = element_text(color="black", size=14),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14)
  )

bp + facet_grid(. ~ Treatment)
