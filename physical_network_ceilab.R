string <- read.csv("String Network - Homo sapiens default edge.csv")

physical <- as.data.frame(read.csv("physical_08.csv"))

physical_genename <- physical$display.name

colors <- read.delim("NodeColorsSD1ALL.tsv")

library(hashmap)
HexFromLabData <- hashmap(keys = colors$Keys, values = colors$Values)
a <- as.data.frame(HexFromLabData[[physical_genename]])

physic_net_colors <- data.frame(V1= physical$display.name, V2= physical$shared.name, V3= a)

physic_net_colors <- na.omit(physic_net_colors)

write.table(physic_net_colors, "physical_net_colors", row.names = F, col.names = T,
            quote = F, sep = "\t")
