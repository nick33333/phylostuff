if(!require("remotes")) install.packages("remotes")
remotes::install_github("daijiang/megatrees")
package_list <- c("tidytree", "castor", "furrr", "future", "fastmatch", "data.table") # rtrees prereqs
for (x in package_list) { 
  print(x)
  install.packages(x)
}
install.packages("rtrees", repos = 'https://daijiang.r-universe.dev')

library(rtrees)
library(ape)


bird_sp_list <- read.table("Documents/data/sp_list.txt", header = F)
bird_sp_list <- bird_sp_list$V1

test_bird_list <- sp_list_df(sp_list=bird_sp_list, taxon="bird")


test_tree = get_tree(sp_list = test_bird_list,
                     taxon = "bird",
                     scenario = "at_basal_node",
                     show_grafted = TRUE)

pdf("bird_tree.pdf", width = 8.5, height = 18)
plot(test_tree, no.margin = T)
dev.off()
