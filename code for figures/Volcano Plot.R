library(ggplot2)
aa <- read.delim('111.txt', sep = '\t', stringsAsFactors = FALSE)
aa[which(aa$padj < 0.05 & aa$log2FoldChange <= -1),'sig'] <- 'Down'
aa[which(aa$padj < 0.05 & aa$log2FoldChange >= 1),'sig'] <- 'Up'
aa[which(aa$padj >= 0.05 | abs(aa$log2FoldChange) < 1),'sig'] <- 'None'
ggplot(aa, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 1, size = 1.2) +
  scale_colour_manual(values  = c("xx","xx","xx"), limits = c('Up', 'Down', 'None')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5)) +
  theme(legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.93)) +
  geom_vline(xintercept = c(-1, 1), color = 'gray', size = 0.3) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.3) +
  labs(x = 'Log2 Fold Change', y = 'Log10 P', color = '', title = 'A vs B')

