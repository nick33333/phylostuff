library(motmot)
library(rtrees)
library(ape)
# Would be useful if I had a quick lookup MSMC Curve cluster and a tree which could
# identify samples froma simple hover (possible reason why an app is good).

# See if MSMC Cluster labels give off an interesting Pagel's Lambda
# data is something that I have to come up with
# tree is something that rtrees can handle well, I just need a list of species names

# Read in MSMC cluster data
msmc_data <- read.csv("Documents/data/curve-cluster-table.tsv", sep="\t" , header=T)
rownames(msmc_data) <- msmc_data$Latin.name
# Other data to consider pulling:
# - Order/taxa
# - Geography/habitat
# - Lon./Lat./Elev.

# Trying to make rows unique based on species name:
# > msmc_data[msmc_data$Latin.name == "Paradisaea_raggiana",]


# Read in Bird Phylogeny
bird_sp_list <- read.table("Documents/data/sp_list.txt", header = F)
bird_sp_list <- bird_sp_list$V1
new_bird_sp_list <- gsub("Chroicocephalus_maculipennis", "Larus_maculipennis", bird_sp_list) # Use gsub to replace some species names
new_bird_sp_list <- gsub("Rhadina_sibilatrix", "Phylloscopus_sibilatrix", new_bird_sp_list)

test_bird_list <- sp_list_df(sp_list=new_bird_sp_list, taxon="bird")


# Can interactively view trees with: https://shiny.rstudio.com/gallery/phylo-tree-view-subset.html
# Export file to wrkdir using following r func/command:
# > write.tree(phy=msmc_bird_tree_1, file="msmc_data_1_tree.phylo")

# Missing like 2 birds but whatever
msmc_bird_tree = get_tree(sp_list = test_bird_list,
                     taxon = "bird",
                     scenario = "at_basal_node",
                     show_grafted = TRUE,)

msmc_bird_tree_1 <- msmc_bird_tree$EricsonStage2_0001_1000.zip1


data(anolis.data)
data(anolis.tree)
attach(anolis.data)

sortedData <- sortTraitData(phy = anolis.tree, y = anolis.data, 
                            data.name = "Male_SVL", pass.ultrametric = TRUE)

msmc_sortedData <- sortTraitData(phy = msmc_bird_tree_1, y = anolis.data, 
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
