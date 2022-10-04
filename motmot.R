install.packages("motmot")
install.packages("geiger")

library(geiger)
library(motmot)
library(rtrees)
library(ape)
# Would be useful if I had a quick lookup MSMC Curve cluster and a tree which could
# identify samples froma simple hover (possible reason why an app is good).
# See if MSMC Cluster labels give off an interesting Pagel's Lambda
# data is something that I have to come up with
# tree is something that rtrees can handle well, I just need a list of species names

# Other data to consider pulling:
# - Order/taxa
# - Geography/habitat
# - Lon./Lat./Elev.
# Read in MSMC cluster data
msmc_data <- read.csv("data/curve-cluster-table.tsv", sep="\t" , header=T)
ln_counts <- data.frame(table(msmc_data$Latin.name), stringsAsFactors = FALSE)
dup_row_df  <- ln_counts[ln_counts$Freq > 1,]
dup_row_names  <- as.character(dup_row_df$Var1) # vector of names of duplicated species
for (name in dup_row_names) {
  remover_vec <- strtoi(rownames(msmc_data[msmc_data$Latin.name == name, ])) # Vector of row idxs for removing rows with dups
  msmc_data <- msmc_data[-tail(remover_vec, -1), ] # Remove all but the first occurrence of a dup'd species
}
rownames(msmc_data) <- msmc_data$Latin.name # Set indices to latin name of organism in row

# Read in Bird Phylogeny
bird_sp_list <- read.table("data/sp_list.txt", header = F)
bird_sp_list <- bird_sp_list$V1
new_bird_sp_list <- gsub("Chroicocephalus_maculipennis", "Larus_maculipennis", bird_sp_list) # Use gsub to replace some species names
new_bird_sp_list <- gsub("Rhadina_sibilatrix", "Phylloscopus_sibilatrix", new_bird_sp_list)
test_bird_list <- sp_list_df(sp_list=new_bird_sp_list, taxon="bird") # Use bird taxon for sampling tree with rtrees
# Can interactively view trees with: https://shiny.rstudio.com/gallery/phylo-tree-view-subset.html
# Export file to wrkdir using following r func/command:
# > write.tree(phy=msmc_bird_tree_1, file="msmc_data_1_tree.phylo")
msmc_bird_tree = get_tree(sp_list = test_bird_list, # Generates a bunch of trees (see samples from rtrees docs)
                     taxon = "bird",
                     scenario = "at_basal_node",
                     show_grafted = TRUE,)
msmc_bird_tree_1 <- msmc_bird_tree$EricsonStage2_0001_1000.zip1 # Single tree

msmc_sortedData <- sortTraitData(phy = msmc_bird_tree_1, y = msmc_data, 
                                 data.name = "Labels", pass.ultrametric = TRUE)
msmc_phy <- msmc_sortedData$phy
msmc_phy_edit <- head(msmc_phy, -2) # Shortened edit on phy to match fields of anolis example
msmc_labels <- msmc_sortedData$trait

### LEFT OFF HERE: Might need to change plot labeling as I'm using categorial data (discrete traits, not continuous which function claims to use)
# Check out fitDiscrete function as opposed to he sortTraitData thing
traitData.plot(y = msmc_labels, phy = msmc_phy_edit, 
               col.label = "red",
               )


#### ANOLIS EXAMPLE STARTS HERE ####
data(anolis.data)
data(anolis.tree)
attach(anolis.data)
sortedData <- sortTraitData(phy = anolis.tree, y = anolis.data, 
                            data.name = "Male_SVL", pass.ultrametric = TRUE)


phy <- sortedData$phy
male.length <- sortedData$trait

traitData.plot(y = male.length, phy = phy, lwd.traits = 2, 
               col.label = "#00008050", tck = -0.01, mgp = c(0, 0.2, 0), 
               cex.axis = 0.5, show.tips = FALSE)

## uncomment to view the tree
plot(phy, show.tip.label=FALSE, no.margin=TRUE, edge.col="grey20")

CLADE_NUM <- 161

nodelabels(CLADE_NUM, CLADE_NUM, bg="black", col="white")
phy.clade <- extract.clade(phy, CLADE_NUM)
temp.mat <- male.length[match(phy.clade$tip.label, rownames(male.length)), ]
male.length.clade <- as.matrix(temp.mat)


lambda.ml <- transformPhylo.ML(phy = phy.clade, y = male.length.clade, 
                               model = "lambda")

lambda.ml
