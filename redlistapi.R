# sudo apt install libcurl4-openssl-dev # Run in terminal
install.packages('curl')
install.packages('crul')
install.packages('rredlist')
install.packages('usethis')

library(rredlist)
library(usethis)
rl_use_iucn()

# API Manual: chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://cran.r-project.org/web/packages/rredlist/rredlist.pdf
# API online docs: https://apiv3.iucnredlist.org/api/v3/docs

key = "8bc3fdadc691427d5121a49b991e6596fa9de78244fd28387938e491b1ae2378"

name = "Xiphorhynchus elegans"

# Start getting data here
name_habitat <- rl_habitats(
  name=name,
  key=key
)

rl_search(
  name=name,
  key=key,
)

rl_history(
  name=name,
  key=key,
)

rl_narrative(
  name=name,
  key=key,
)


# rl_habitats('Fratercula arctica', key=key)
# rl_habitats('Fratercula arctica', region = 'europe', key=key)
# rl_habitats(id = 12392, key=key)
# rl_habitats(id = 22694927, region = 'europe', key=key)
# 
# rl_habitats_('Fratercula arctica', key=key)
# rl_habitats_(id = 12392, key=key)

