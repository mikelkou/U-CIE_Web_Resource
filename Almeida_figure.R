# Almeida figure
almeida <- read.delim("~/Desktop/CBPP_22012021/Rotation_Simon/for_mikaela/Almeida/almeida.jgi.by_cluster.weighted_average.good_clusters.tsv")
clusters <- read.delim("~/Desktop/CBPP_22012021/Rotation_Simon/for_mikaela/Almeida/Supplementary Data 5 - Copy.txt")

genera_from_clusters = clusters[,c("clusters","g")]
colors = read.delim("~/Downloads/Hex_codes-Lab_coords-almeida.jgi.by_cluster.weighted_average.good_clusters.tsv-2021-10-22.tsv")

unique_clusters <- genera_from_clusters[!duplicated(genera_from_clusters[,'clusters']),]
genera = unique_clusters$g
genera = gsub(pattern = "s__", replacement = "", x = genera)

genera = gsub(pattern = c("_Z") , replacement = "", x = genera)
genera = gsub(pattern = c(" .*") , replacement = "", x = genera)
genera = gsub(pattern = c("RuminococcusRuminococcus") , replacement = "Ruminococcus", x = genera)
genera = gsub(pattern = c("LactobacillusLactobacillus") , replacement = "Lactobacillus", x = genera)
genera = gsub(pattern = c("RuminiclostridiumRuminiclostridium") , replacement = "Ruminiclostridium", x = genera)


uniqclusters = cbind(unique_clusters$clusters, genera)
uniqclusters = subset(uniqclusters, uniqclusters[,2]!="")

write.table(uniqclusters, "~/Desktop/species_almeida_clusters.tsv", quote = F, row.names = F, col.names = T, sep="\t")



tax_ids = read.delim("~/Desktop/taxonomy_result.txt", h=F)
tax_ids_results = read.delim("~/Desktop/taxonomy_result copy.txt", h=F)


#------------------------------------------------------------------------------#
colnames(uniqclusters) = c("clusters", "genera")

ncbi_taxonomy = read.delim("~/Desktop/CBPP_22012021/Lars_Lab/Almeida_figure_UCIE/tree_with_names.csv", sep = ",", header = T)
ncbi_taxonomy = ncbi_taxonomy[,2:5]

uniqclusters = as.data.frame(uniqclusters)
colnames(uniqclusters) = c("X2", "X1")
uniqclusters = cbind(uniqclusters[,2], uniqclusters[,1])

ncbi_taxonomy_net = ncbi_taxonomy[,c(2,4)]
tree = inner_join()
tree = rbind(ncbi_taxonomy_net, uniqclusters)

write.table(tree, "~/Desktop/CBPP_22012021/Lars_Lab/Almeida_figure_UCIE/final_tree_ncbi_and_clusters.tsv", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(uniqclusters, "~/Desktop/CBPP_22012021/Lars_Lab/Almeida_figure_UCIE/labels_clusters.tsv", sep = "\t", quote = F, row.names = F, col.names = F)



#---#
## recursion function
traverse <- function(a,i,innerl){
  if(i < (ncol(df))){
    alevelinner <- as.character(unique(df[which(as.character(df[,i])==a),i+1]))
    desc <- NULL
    if(length(alevelinner) == 1) (newickout <- traverse(alevelinner,i+1,innerl))
    else {
      for(b in alevelinner) desc <- c(desc,traverse(b,i+1,innerl))
      il <- NULL; if(innerl==TRUE) il <- a
      (newickout <- paste("(",paste(desc,collapse=","),")",il,sep=""))
    }
  }
  else { (newickout <- a) }
}

## data.frame to newick function
df2newick <- function(df, innerlabel=FALSE){
  alevel <- as.character(unique(df[,1]))
  newick <- NULL
  for(x in alevel) newick <- c(newick,traverse(x,1,innerlabel))
  (newick <- paste("(",paste(newick,collapse=","),");",sep=""))
}
#---#
df <- data.frame(x=c('A','A','B','B','B'), y=c('Ab','Ac','Ba', 'Ba','Bd'), z=c('Abb','Acc','Bad', 'Bae','Bdd'))
myNewick <- df2newick(df)
library(ape)
mytree <- read.tree(text=myNewick)
plot(mytree)




tree_with_clusters = aggregate(clusters ~., uniqclusters, toString )

library(stringr)
gggg = str_split_fixed(tree_with_clusters$clusters, ",", 85)
gggg = as.data.frame(gggg)

tree_with_clusters = cbind(tree_with_clusters[,1], gggg)
colnames(tree_with_clusters) = c(paste("X", 1:ncol(tree_with_clusters), sep = ""))

library(dplyr)
final_tree_ncbi_and_clusters = inner_join(tree_with_clusters, ncbi_taxonomy, "X1")


#--- Figure TARA Oceans ---#
transposed_mOTUs = read.delim("~/Desktop/CBPP_22012021/Lars_Lab/Figure_UCIE/TARA_OCEANS_figure/mOTU.linkage-groups.relab.release.tsv")
rownames(transposed_mOTUs) = transposed_mOTUs[,1]
transposed_mOTUs = transposed_mOTUs[,2:ncol(transposed_mOTUs)]
transposed_mOTUs = t(transposed_mOTUs)

write.table(transposed_mOTUs, "~/Desktop/transposed_mOTUs.tsv", sep = "\t", quote = F, row.names = T, col.names = T)






#--- 16S OTU ---#
taxonomic_profiles = read.delim("~/Desktop/CBPP_22012021/Lars_Lab/Figure_UCIE/TARA_OCEANS_figure/miTAG.taxonomic.profiles.release.tsv")
row_names = paste(taxonomic_profiles[,1],taxonomic_profiles[,2], taxonomic_profiles[,3], taxonomic_profiles[,4], taxonomic_profiles[,5], taxonomic_profiles[,6],taxonomic_profiles[,7], sep = "_")
row_names = as.data.frame(row_names)
taxonomic_profiles = taxonomic_profiles[,8:ncol(taxonomic_profiles)]
rownames(taxonomic_profiles) = row_names[,1]
taxonomic_profiles = t(taxonomic_profiles)
write.table(taxonomic_profiles, "~/Desktop/CBPP_22012021/Lars_Lab/Figure_UCIE/TARA_OCEANS_figure/taxonomic_profiles.tsv", sep = "\t", quote = F, row.names = T, col.names = T)
#---#
tax_prof_srf = subset(taxonomic_profiles, grepl( "SRF", rownames(taxonomic_profiles), fixed = TRUE)==TRUE) #SRF only

data <- uwot::umap(tax_prof_srf, ret_nn = TRUE, n_neighbors = 5, n_components = 3) # colors_tax_profs
data_umap_coord <- as.data.frame(data$embedding)
rownames(data_umap_coord) <- rownames(tax_prof_srf)
write.table(colors_tax_profs_srf, "~/Desktop/MAP_FIGURE_CIELAB_coords_of_taxonomic_profile_SRF_3D.tsv", quote = F, row.names = F, col.names = T, sep = "\t")

# colors_tax_profs = data2cielab(data_umap_coord)
colors_tax_profs_srf = data2cielab(data_umap_coord, Wa = 3, Wb = 2, S = 0.8, LAB_coordinates = T)
#-------------------------------------------------------------------------------#
data <- uwot::umap(transposed_mOTUs, ret_nn = TRUE, n_neighbors = 5, n_components = 3) # colors_mOTUs
data_umap_coord <- as.data.frame(data$embedding)
rownames(data_umap_coord) <- rownames(taxonomic_profiles)
library(ucie)
colors_tax_profs = data2cielab(data_umap_coord)


colors_mOTUs_online = read.delim("~/Desktop/CBPP_22012021/Lars_Lab/Figure_UCIE/TARA_OCEANS_figure/Hex_codes-Lab_coords-transposed_mOTUs.tsv-2021-11-23.tsv")
colors_mOTUs_online = colors_mOTUs_online[,1:2]

fig <- plot_ly(
  data = as.data.frame(data_umap_coord),
  x = colors_tax_profs[,2],
  y = ,
  # z = data_umap_coord[, 3],
  source = "A",
  type = "scatter3d",
  mode = "markers",
  text = c(rownames(data_umap_coord)),
  # hoverinfo = 'text',
  marker = list(
    color = colors_mOTUs[,2],
    size = 3,
    width = 2
  )
  ,
  width = 1000, height = 720
)

fig

library(ggmap)
library(ggplot2)
library(maptools)
library(maps)

# Excel file --> Table W1 columns with Latitude [degrees North]	Longitude [degrees East]	Sampling depth [m] 
coordinates = read.delim("~/Desktop/CBPP_22012021/Lars_Lab/Figure_UCIE/TARA_OCEANS_figure/Coordinates_TATA.tsv")
coordinates = coordinates[,c(1,3,7,9,10,11,12)]
coordinates_correct = gsub("-", ".", coordinates[,1])
coordinates = cbind(V1=as.data.frame(coordinates_correct), coordinates)

# coordinates_DCM = subset(coordinates, coordinates$Environmental.Feature=="(DCM) deep chlorophyll maximum layer (ENVO:01000326)")
# coordinates_SRF = subset(coordinates, coordinates$Environmental.Feature=="(SRF) surface water layer (ENVO:00002042)")
# coordinates_MES = subset(coordinates, coordinates$Environmental.Feature=="(MES) mesopelagic zone (ENVO:00000213)")
# #---#

colnames(colors_tax_profs_srf) = c("coordinates_correct","colors" )

#---#
stations = full_join(coordinates, colors_tax_profs_srf, "coordinates_correct") # colors_tax_profs or colors_mOTUs or colors_mOTUs_online
stations = na.omit(stations)
stations_DCM = subset(stations, stations$Environmental.Feature=="(DCM) deep chlorophyll maximum layer (ENVO:01000326)")
stations_SRF = subset(stations, stations$Environmental.Feature=="(SRF) surface water layer (ENVO:00002042)")
stations_MES = subset(stations, stations$Environmental.Feature=="(MES) mesopelagic zone (ENVO:00000213)")


x = as.numeric(gsub(",",".",stations[,5]))
y = as.numeric(gsub(",",".",stations[,4]))
#---#
mp <- NULL
mapWorld <- borders("world", colour="#e7e7e7", fill="#e7e7e7") # create a layer of borders
mp <- ggplot() + mapWorld

#Now Layer the cities on top
mp <- mp + geom_point(aes(x=x, y=y) ,color=stations$colors, size=3) +
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank()) + 
  theme(plot.background = element_rect(fill = 'white', colour = 'black'))

mp




#--- KEGG ---#
kegg_functional = read.delim("~/Desktop/CBPP_22012021/Lars_Lab/Figure_UCIE/TARA_OCEANS_figure/TARA243.KO.profile.release.tsv")
rownames(kegg_functional) = kegg_functional[,1]
kegg_functional = kegg_functional[,2:ncol(kegg_functional)]
kegg_functional = t(kegg_functional)
data <- uwot::umap(kegg_functional, ret_nn = TRUE, n_neighbors = 5, n_components = 3) # colors_mOTUs
data_umap_coord <- as.data.frame(data$embedding)
rownames(data_umap_coord) <- rownames(kegg_functional)
kegg_colors = data2cielab(data_umap_coord)
colnames(kegg_colors) = c("coordinates_correct", "colors")

stations = full_join(coordinates, kegg_colors, "coordinates_correct")
stations = na.omit(stations)

stations_DCM = subset(stations, stations$Environmental.Feature=="(DCM) deep chlorophyll maximum layer (ENVO:01000326)")
stations_SRF = subset(stations, stations$Environmental.Feature=="(SRF) surface water layer (ENVO:00002042)")
stations_MES = subset(stations, stations$Environmental.Feature=="(MES) mesopelagic zone (ENVO:00000213)")

x = as.numeric(gsub(",",".",stations[,5]))
y = as.numeric(gsub(",",".",stations[,4]))
#---#

mp <- NULL
mapWorld <- borders("world", colour="gray80", fill="gray80") # create a layer of borders
mp <- ggplot() + mapWorld

#Now Layer the cities on top
mp <- mp + geom_point(aes(x=x, y=y) ,color=stations$colors, size=3) +
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank()) + 
  theme(plot.background = element_rect(fill = 'white', colour = 'black'))

mp

#--- KEGG module profiles ---#
kegg_module = read.delim("~/Desktop/CBPP_22012021/Lars_Lab/Figure_UCIE/TARA_OCEANS_figure/TARA243.KO-module.profile.release.tsv")
rownames(kegg_module) = kegg_module[,1]
kegg_module = kegg_module[,2:ncol(kegg_module)]
kegg_module = t(kegg_module)
data <- uwot::umap(kegg_module, ret_nn = TRUE, n_neighbors = 5, n_components = 3) # colors_mOTUs
data_umap_coord <- as.data.frame(data$embedding)
rownames(data_umap_coord) <- rownames(kegg_module)
kegg_module_colors = data2cielab(data_umap_coord)
colnames(kegg_module_colors) = c("coordinates_correct", "colors")

stations = full_join(coordinates, kegg_module_colors, "coordinates_correct")
stations = na.omit(stations)

stations_DCM = subset(stations, stations$Environmental.Feature=="(DCM) deep chlorophyll maximum layer (ENVO:01000326)")
stations_SRF = subset(stations, stations$Environmental.Feature=="(SRF) surface water layer (ENVO:00002042)")
stations_MES = subset(stations, stations$Environmental.Feature=="(MES) mesopelagic zone (ENVO:00000213)")

x = as.numeric(gsub(",",".",stations_MES[,5]))
y = as.numeric(gsub(",",".",stations_MES[,4]))
#---#


mp <- NULL
mapWorld <- borders("world", colour="gray80", fill="gray80") # create a layer of borders
mp <- ggplot() + mapWorld

#Now Layer the cities on top
mp <- mp + geom_point(aes(x=x, y=y) ,color=stations_MES$colors, size=3) +
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank()) + 
  theme(plot.background = element_rect(fill = 'white', colour = 'black'))

mp


#------------------------------------------------------------------------------#
# Correlations Heatmap
coordinates = read.delim("~/Desktop/CBPP_22012021/Lars_Lab/Figure_UCIE/TARA_OCEANS_figure/Coordinates_TATA.tsv")
coordinates = coordinates[,c(1,6,7,9,10,11,12)]
coordinates_correct = gsub("-", ".", coordinates[,1])
coordinates = cbind(V1=as.data.frame(coordinates_correct), coordinates)
coordinates = subset(coordinates, grepl( "SRF", coordinates[,1], fixed = TRUE)==TRUE) #SRF only

colnames(colors_tax_profs_srf)[1] = "coordinates_correct"
coords_cielabcoords = full_join(coordinates, colors_tax_profs_srf, "coordinates_correct")
coords_cielabcoords = na.omit(coords_cielabcoords)
colnames(coords_cielabcoords)[3] ="PANGAEA.Sample.ID"
features = read.delim("~/Desktop/CBPP_22012021/Lars_Lab/Figure_UCIE/TARA_OCEANS_figure/pangaea_tempr_features.tsv")

merged_files_features_coords_and_cielabcoords = full_join(coords_cielabcoords, features, by = "PANGAEA.Sample.ID")
# write.table(merged_files_features_coords, "~/Desktop/CBPP_22012021/Lars_Lab/Figure_UCIE/TARA_OCEANS_figure/merged_files_features_coords.tsv", quote = F, row.names = F, col.names = T, sep = "\t")
merged_files_features_coords_and_cielabcoords_subset = merged_files_features_coords_and_cielabcoords[,c(1,3,5:ncol(merged_files_features_coords_and_cielabcoords))]

merged_files_features_coords_and_cielabcoords_subset = subset(merged_files_features_coords_and_cielabcoords_subset, !is.na(merged_files_features_coords_and_cielabcoords_subset$L))
# colnames(merged_files_features_coords_and_cielabcoords_subset)[21] = "SI"
heatmap(x = data.matrix(merged_files_features_coords_and_cielabcoords_subset[,c(7:9,14:21)],), labRow = merged_files_features_coords_and_cielabcoords_subset[,1], margins = c(10,10))
heatmap(x = data.matrix(merged_files_features_coords_and_cielabcoords_subset[,c(9,21)]), labRow = merged_files_features_coords_and_cielabcoords_subset[,1], margins = c(20,20), cexCol = 2)

temp_L = data.matrix(merged_files_features_coords_and_cielabcoords_subset[, 7:ncol(merged_files_features_coords_and_cielabcoords_subset)])
temp_L[is.na(temp_L)] <- 0
correlation = cor(temp_L, method = "pearson")

corr_crop = correlation[1:3,-c(1:6)]

df = c()
for(i in 1:ncol(corr_crop)){
  a = abs(corr_crop)
  a = max(a[,i])
  df = cbind(df, a)
}

colnames(df) = colnames(corr_crop)
df = t(df)
df = subset(df, df[,1]>0.4)

df = df[order(df[,1],decreasing = T),]
df = as.data.frame(df)


corr_crop_2 = corr_crop[,c("Mean_Temperature..deg.C..","Mean_Oxygen..umol.kg..",
                            "NO2NO3..umol.L...", "NO2..umol.L..." )]
colnames(corr_crop_2)=c("Mean Temperature [deg.C.]",  "Mean Oxygen [umol/kg]", 
                         "NO2NO3 [umol/L]", "NO2 [umol/L]")


# colnames(corr_crop) = c("Depth [m]", "Temperature [deg.C.]", "Salinity [PSU]", "Oxygen [umol/kg]", "Nitrates [umol/L]", 
#                          "NO2 [umol/L]", "PO4 [umol/L]", "NO2NO3 [umol/L]", "SI [umol/L]", "AMODIS.PAR8d.Einsteins", colnames(correlation)[23:ncol(correlation)])
heatmap(x = data.matrix(corr_crop_2*(-1)), margins = c(20,20), scale = "none",cexCol = 1, cexRow = 1.5, col= colorRampPalette(brewer.pal(8, "Spectral"))(25), 
        Rowv = NA, Colv = NA)



colors_tax_profs_srf_notcoords = read.delim("~/Desktop/CBPP_22012021/Lars_Lab/Figure_UCIE/TARA_OCEANS_figure/MAP_FIGURE_DATA_Colors_of_taxonomic_profile_SRF_3D.tsv")
colnames(colors_tax_profs_srf_notcoords)[1] = "coordinates_correct"
colors_everything_eslse = full_join(merged_files_features_coords_and_cielabcoords_subset, colors_tax_profs_srf_notcoords, "coordinates_correct")

# Scatterplots: L -- oxygen | a -- NO2NO3 | L,b -- temprature
# a NO2NO3
fig <- plot_ly(
  data = as.data.frame(merged_files_features_coords_and_cielabcoords_subset),
  x = temp_L[,"a"], # L
  y = temp_L[,"NO2NO3..umol.L..."], #Mean_Oxygen..umol.kg.. #NO2NO3..umol.L...#Mean_Temperature..deg.C..
  type = "scatter",
  mode = "markers",
  text = c(merged_files_features_coords_and_cielabcoords_subset[,1]), 
  # hoverinfo = 'text',
  marker = list(
    color = colors_everything_eslse$colors,
    size = 20,
    width = 2
  )
  ,
  width = 1000, height = 720
)
# fig <- fig %>% add_trace(y = ~na.omit(merged_files_features_coords_and_cielabcoords_subset[,"Mean_Oxygen..umol.kg.."]), name = 'trace 0',mode = 'lines')
fig = fig %>% layout(xaxis = list(title = 'a*',
                                  tickmode = 'array',
                                  tickvals = c(-50,0,50),
                                  ticktext = c("-50","0","50"), range = c(-50,50),
                                  zerolinecolor = "efefef" #'#efefef'
                                  # gridcolor = '#ffff',
                                  # gridcolor = 'ffff'
             ),
             # xaxis = list(
             #   zerolinecolor = '#ffff',
             #   zerolinewidth = 2,
             #   gridcolor = 'ffff'),
             # yaxis = list(
             #   zerolinecolor = '#ffff',
             #   zerolinewidth = 2,
             #   gridcolor = 'ffff'),
             yaxis = list(title = 'NO2NO3 [umol/L]',
                          tickmode = 'array',
                          tickvals = c(17.5,35),
                          ticktext = c("17.5", "35"), range = c(-1,40),
                          gridcolor = 'white',
                          gridcolor = 'white')
             
) %>%
  layout(autosize = T, margin=list( l = 50, r = 50, b = 100, t = 100,  pad = 4)) %>%
  layout( xaxis = list(titlefont = list(size = 22), tickfont = list(size = 22)),
          yaxis = list(titlefont = list(size = 22), tickfont = list(size = 22)) )
fig
#------------------------------------------------------------------------------#
#L Oxygen
fig <- plot_ly(
  data = as.data.frame(merged_files_features_coords_and_cielabcoords_subset),
  x = temp_L[,"L"], # L
  y = temp_L[,"Mean_Oxygen..umol.kg.."], #Mean_Oxygen..umol.kg.. #NO2NO3..umol.L...#Mean_Temperature..deg.C..
  type = "scatter",
  mode = "markers",
  text = c(merged_files_features_coords_and_cielabcoords_subset[,1]), 
  # hoverinfo = 'text',
  marker = list(
    color = colors_everything_eslse$colors,
    size = 20,
    width = 2
  )
  ,
  width = 1000, height = 720
)
# fig <- fig %>% add_trace(y = ~na.omit(merged_files_features_coords_and_cielabcoords_subset[,"Mean_Oxygen..umol.kg.."]), name = 'trace 0',mode = 'lines')
fig = fig %>% layout(xaxis = list(title = 'L*',
                                  tickmode = 'array',
                                  tickvals = c(0,50,100),
                                  ticktext = c("0","50","100"), range = c(-1,100),
                                  zerolinecolor = "efefef" #'#efefef'
                                  # gridcolor = '#ffff',
                                  # gridcolor = 'ffff'
             ),
             # xaxis = list(
             #   zerolinecolor = '#ffff',
             #   zerolinewidth = 2,
             #   gridcolor = 'ffff'),
             # yaxis = list(
             #   zerolinecolor = '#ffff',
             #   zerolinewidth = 2,
             #   gridcolor = 'ffff'),
             yaxis = list(title = 'Oxygen [umol/kg]',
                          tickmode = 'array',
                          tickvals = c(30,60),
                          ticktext = c("30", "60"), range = c(-2,60),
                          gridcolor = 'white',
                          gridcolor = 'white')
             
) %>%
  layout(autosize = T, margin=list( l = 50, r = 50, b = 100, t = 100,  pad = 4)) %>%
  layout( xaxis = list(titlefont = list(size = 22), tickfont = list(size = 22)),
          yaxis = list(titlefont = list(size = 22), tickfont = list(size = 22)) )
fig


#L Temperature
fig <- plot_ly(
  data = as.data.frame(merged_files_features_coords_and_cielabcoords_subset),
  x = temp_L[,"L"], # L
  y = temp_L[,"Mean_Temperature..deg.C.."], #Mean_Oxygen..umol.kg.. #NO2NO3..umol.L...#Mean_Temperature..deg.C..
  type = "scatter",
  mode = "markers",
  text = c(merged_files_features_coords_and_cielabcoords_subset[,1]), 
  # hoverinfo = 'text',
  marker = list(
    color = colors_everything_eslse$colors,
    size = 20,
    width = 2
  )
  ,
  width = 1000, height = 720
)
# fig <- fig %>% add_trace(y = ~na.omit(merged_files_features_coords_and_cielabcoords_subset[,"Mean_Oxygen..umol.kg.."]), name = 'trace 0',mode = 'lines')
fig = fig %>% layout(xaxis = list(title = 'L*',
                                  tickmode = 'array',
                                  tickvals = c(0,50,100),
                                  ticktext = c("0","50","100"), range = c(-1,100),
                                  zerolinecolor = "efefef" #'#efefef'
                                  # gridcolor = '#ffff',
                                  # gridcolor = 'ffff'
             ),
             # xaxis = list(
             #   zerolinecolor = '#ffff',
             #   zerolinewidth = 2,
             #   gridcolor = 'ffff'),
             # yaxis = list(
             #   zerolinecolor = '#ffff',
             #   zerolinewidth = 2,
             #   gridcolor = 'ffff'),
             yaxis = list(title = 'Temperature [deg.C]',
                          tickmode = 'array',
                          tickvals = c(25, 50),
                          ticktext = c("25", "50"), range = c(-1,60),
                          gridcolor = 'white',
                          gridcolor = 'white')
             
) %>%
  layout(autosize = T, margin=list( l = 50, r = 50, b = 100, t = 100,  pad = 4)) %>%
  layout( xaxis = list(titlefont = list(size = 22), tickfont = list(size = 22)),
          yaxis = list(titlefont = list(size = 22), tickfont = list(size = 22)) )
fig

#b temperature
fig <- plot_ly(
  data = as.data.frame(merged_files_features_coords_and_cielabcoords_subset),
  x = temp_L[,"b"], # L
  y = temp_L[,"Mean_Temperature..deg.C.."], #Mean_Oxygen..umol.kg.. #NO2NO3..umol.L...#Mean_Temperature..deg.C..
  type = "scatter",
  mode = "markers",
  text = c(merged_files_features_coords_and_cielabcoords_subset[,1]), 
  # hoverinfo = 'text',
  marker = list(
    color = colors_everything_eslse$colors,
    size = 20,
    width = 2
  )
  ,
  width = 1000, height = 720
)
# fig <- fig %>% add_trace(y = ~na.omit(merged_files_features_coords_and_cielabcoords_subset[,"Mean_Oxygen..umol.kg.."]), name = 'trace 0',mode = 'lines')
fig = fig %>% layout(xaxis = list(title = 'b*',
                                  tickmode = 'array',
                                  tickvals = c(-50,0,50),
                                  ticktext = c("-50","0","50"), range = c(-50,50),
                                  zerolinecolor = "efefef" #'#efefef'
                                  # gridcolor = '#ffff',
                                  # gridcolor = 'ffff'
             ),
             # xaxis = list(
             #   zerolinecolor = '#ffff',
             #   zerolinewidth = 2,
             #   gridcolor = 'ffff'),
             # yaxis = list(
             #   zerolinecolor = '#ffff',
             #   zerolinewidth = 2,
             #   gridcolor = 'ffff'),
             yaxis = list(title = 'Temperature [deg.C]',
                          tickmode = 'array',
                          tickvals = c(25,50),
                          ticktext = c("25","50"), range = c(-1,60),
                          gridcolor = 'white',
                          gridcolor = 'white')
             
) %>%
  layout(autosize = T, margin=list( l = 50, r = 50, b = 100, t = 100,  pad = 4)) %>%
  layout( xaxis = list(titlefont = list(size = 22), tickfont = list(size = 22)),
          yaxis = list(titlefont = list(size = 22), tickfont = list(size = 22)) )
fig


#L,b temperature
fig <- plot_ly(
  data = as.data.frame(merged_files_features_coords_and_cielabcoords_subset),
  x = temp_L[,"L"], # L
  y = temp_L[,"b"], #Mean_Oxygen..umol.kg.. #NO2NO3..umol.L...#Mean_Temperature..deg.C..
  type = "scatter",
  mode = "markers",
  text = c(merged_files_features_coords_and_cielabcoords_subset[,1]), 
  # hoverinfo = 'text',
  name = 'L',
  marker = list(
    color = colors_everything_eslse$colors,
    # symbol = "square",
    size = 20,
    width = 2
  )
  ,
  width = 1000, height = 720
)
# fig <- fig %>% add_trace(y = ~na.omit(merged_files_features_coords_and_cielabcoords_subset[,"Mean_Oxygen..umol.kg.."]), name = 'trace 0',mode = 'lines')
fig <- fig %>% add_trace(y = temp_L[,"Mean_Temperature..deg.C.."], name = 'b', mode = 'markers',
                         marker = list(symbol = 'square')) #hexagon2

fig = fig %>% layout(xaxis = list(title = 'L,b*',
                                  tickmode = 'array',
                                  tickvals = c(-50,0,50, 100),
                                  ticktext = c("-50","0","50", "100"), range = c(-50,100),
                                  zerolinecolor = "efefef" #'#efefef'
                                  # gridcolor = '#ffff',
                                  # gridcolor = 'ffff'
             ),
             # xaxis = list(
             #   zerolinecolor = '#ffff',
             #   zerolinewidth = 2,
             #   gridcolor = 'ffff'),
             # yaxis = list(
             #   zerolinecolor = '#ffff',
             #   zerolinewidth = 2,
             #   gridcolor = 'ffff'),
             yaxis = list(title = 'Temperature [deg.C]',
                          tickmode = 'array',
                          tickvals = c(25),
                          ticktext = c("25"), range = c(-1,60),
                          gridcolor = 'white',
                          gridcolor = 'white')
             
) %>%
  layout(autosize = T, margin=list( l = 50, r = 50, b = 100, t = 100,  pad = 4)) %>%
  layout( xaxis = list(titlefont = list(size = 22), tickfont = list(size = 22)),
          yaxis = list(titlefont = list(size = 22), tickfont = list(size = 22)) )
fig


# color gradients
colfunc <- colorRampPalette(c("black", "white"))

colfunc <- colorRampPalette(c("black", "gray"))
colfunc <- colorRampPalette(c("gray", "white"))

# colfunc <- colorRampPalette(c("green", "red"))

colfunc(6)
plot(rep(1,6),col=colfunc(6),pch=19,cex=3, type = "l")

#- color lines -#
suppressMessages(library(tidyverse))

slider <- data.frame("values"=seq(0, 1, by = 0.01))
dot_value <- data.frame("Result"= 0.8)


ggplot(slider, aes(x=values, y=0)) +
  geom_line(aes(color = values),
            size=1.5)+
  # geom_point(data=dot_value,
  #            mapping = aes(as.numeric(Result)),size=4)+
  # scale_color_gradient(low = "green", high = "red", mid = "gray",) +
  # scale_color_gradient(low = "green", high = "red") +
  scale_colour_gradient2(low = "darkgreen", mid="gray", high = "red", midpoint=0.5) +
  
  theme(axis.ticks.y=element_blank(),
        axis.text.y=element_blank()) +
  xlab("")+
  ylab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  


