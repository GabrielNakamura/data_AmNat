####installing mcfly package, dependences and reading libraries####
library(here)
library(devtools)
install_github("sokole/MCSim")
install_github("GabrielNakamura/mcfly")

###PS: relauch R after instalation

####Load general data and test with Sigmodontinae####
comm.full<- as.matrix(read.table(here::here("data", "comm_sigmodontinae_final.txt"),  h= T)) #community matrix
env.full<- read.table(here::here("data", "Env.txt"), h=T) #environmental variables
esp.full<-as.matrix(env.full[,1:2])
phy<-ape::read.tree(here::here("data", "tree_sigmodontinae_final.tre")) #phylogenetic hypothesis for Sigmodontinae species
