library(RColorBrewer)
library(phytools)
library(motmot)
library(rtrees)
library(ape)

source("Documents/phylostuff/delta_statistic/code.R")

# Read in MSMC cluster label data
msmc_data <- read.csv("Documents/phylostuff/data/curve-cluster-table.tsv", sep="\t" , header=T)
ln_counts <- data.frame(table(msmc_data$Latin.name), stringsAsFactors = FALSE)
dup_row_df  <- ln_counts[ln_counts$Freq > 1,]
dup_row_names  <- as.character(dup_row_df$Var1) # vector of names of duplicated species
for (name in dup_row_names) {
  remover_vec <- strtoi(rownames(msmc_data[msmc_data$Latin.name == name, ])) # Vector of row idxs for removing rows with dups
  msmc_data <- msmc_data[-tail(remover_vec, -1), ] # Remove all but the first occurrence of a dup'd species
}
rownames(msmc_data) <- msmc_data$Latin.name # Set indices to latin name of organism in row

# Read in Bird Phylogeny from rtrees 
bird_sp_list <- read.table("Documents/data/sp_list.txt", header = F)
bird_sp_list <- bird_sp_list$V1
new_bird_sp_list <- gsub("Chroicocephalus_maculipennis", "Larus_maculipennis", bird_sp_list) # Use gsub to replace some species names
new_bird_sp_list <- gsub("Rhadina_sibilatrix", "Phylloscopus_sibilatrix", new_bird_sp_list)
test_bird_list <- sp_list_df(sp_list=new_bird_sp_list, taxon="bird")

# Create tree with rtrees get_tree function which references bird phylos
msmc_bird_tree <- get_tree(sp_list = test_bird_list,
                          taxon = "bird",
                          scenario = "at_basal_node",
                          show_grafted = TRUE,)
msmc_bird_tree_1 <- msmc_bird_tree$EricsonStage2_0001_1000.zip1 # Sample the first tree in list (1/100)
# msmc_bird_tree_1 <- msmc_bird_tree[["EricsonStage2_0001_1000.zip1"]]
plot(msmc_bird_tree_1, show.tip.label=FALSE, no.margin=TRUE, edge.col="grey20") # Plot Tree


# Use delta-stat
msmc_bird_tree_1$edge.length[msmc_bird_tree_1$edge.length==0] <- quantile(msmc_bird_tree_1$edge.length,0.1)*0.1
trait_cluster_label <- msmc_data[msmc_bird_tree_1$tip.label, "Labels"]

fixed_phylo <- multi2di(msmc_bird_tree_1, random = TRUE, equiprob = TRUE) # Scrappy attempt to remove polytomies/multitomies in phylo
# Resulting tree looks kinda ok, close enough to the original. Look at before and
# after plots of the phylo to make the call. After quick fix, phylo looks more boxy than before
plot(fixed_phylo, show.tip.label=FALSE, no.margin=TRUE, edge.col="grey20")
  

deltaA <- delta(trait_cluster_label, fixed_phylo ,0.1,0.0589,10000,10,100) # Takes a lil bit. First try p=0.0.2709212


#DRAW THE TREE... pulled and adapted from mrborges23
ns <- length(fixed_phylo$tip.label)
par(mfrow=c(1,2))
fixed_phylo$tip.label <- rep("",ns)
plot(fixed_phylo,main="SCENARIO A")
ar <- ace(trait_cluster_label,fixed_phylo,type="discrete",method="ML",model="ARD")$lik.anc # Also takes a long time, finds internal node
# probs and entropies. See: https://academic.oup.com/bioinformatics/article/35/11/1862/5144670?login=true
nodelabels(pie = ar, cex = 1,frame="n") 
mtrait <- matrix(0,ncol=14,nrow=ns)
for ( i in 1:ns) {
  mtrait[i,trait_cluster_label[i]] <- 1
}
tiplabels(pie=mtrait,cex=0.5)


# Draw a tree with cluster labels hopefully 
plot.new()
frame()
# Randomly sampled colors from running: sample(colors(), 14)
# color_vec <- c("gray4", "green", "palevioletred" , "lightgreen", "purple", "navyblue","deeppink", "orangered", 
#                 "lightyellow2",   "brown", "lightslateblue", "lightskyblue","lightblue1","slateblue4", "white")  
# 
# fmode<-as.factor(setNames(msmc_data$Labels,rownames(msmc_data)))
# dotTree(fixed_phylo,fmode,labels=TRUE, colors=setNames(color_vec,
#                                        c("0", "1", "2", "3", "4", "5","6","7","8","9","10","11","12","13", "14")),
#                                       ftype="i", width=4, height=100)


# Might want to order colors and labels by centroid similarity of curve clusters
colors <- RColorBrewer::brewer.pal(n=12,"Paired")
colors <- append(colors, "black")
colors <- append(colors, "magenta")
colors <- append(colors, "white")

cols<-setNames(colors,
               c("0", "1", "2", "3", "4", "5","6","7","8","9","10","11","12","13", "14"))
traits2plot <- trait_cluster_label # Substitute NAs in OG trait vector with special label "14" for white
traits2plot[is.na(trait_cluster_label)] <- "14"
fmode<-as.factor(setNames(traits2plot, rownames(msmc_data)))
dotTree(msmc_bird_tree_1, fmode,labels=FALSE, colors=cols,
        ftype="i", width=4, height=10)
