require("mcfly")

####seting parameters to enter in Mcfly function#####
runs<- 100 #number of runs
OU.alpha<- c(0, 0.05, 0.25, 1) #alpha levels to be tested
pc.env<- prcomp(env.full, scale= TRUE) #summarizing environmental variables 
pc1<- pc.env$x[, 1] #extracting PC1
pc2<- pc.env$x[, 2] #extracting PC2
match<- picante::match.phylo.comm(phy, comm.full) #organizing phylogeny and community matrix
comm<- match$comm
phy<- match$phy
rownames(comm)==rownames(env.full)
div<- vegan::renyi(comm,scales=1) #calculating diversity 
mod<- lm(div ~ env.full$Elevation)
pred<- predict.lm(mod)
envir<- scales::rescale(as.matrix(env.full$Elevation), c(1, 100))
sigma= sd(envir)
ED<- cbind(envir, pred)
theta<- ED[which.max(pred),1]
root.value<- mean(envir)

###run mcfly with Sigmodontinae data####
test.sigmodontinae<- mcfly::Mcfly(comm = comm, subset = FALSE, occurrence = TRUE, env = envir, 
                                  site.coords = esp.full, tree = phy, OU.alpha = OU.alpha, sigma = sigma, theta = theta, 
                                  root.value = root.value, runs = runs, ncores = 4, W.r = 0, 
                                  output.dir.path = "Sigmodontinae")

save.image("DadosEmpiricos_Sigmodontinae.RData") #saving results

#### testing with tribes ####
sub.tree<- ape::subtrees(phy)
n.nodes<- length(sub.tree)
nspp.nodes<- matrix(NA, n.nodes, 1)
for(i in 1:n.nodes){
  nspp.nodes[i,]<- sub.tree[[i]]$Ntip
}

dim.comm.subset<- matrix(NA, nspp.nodes, 2)
for(j in 1:n.nodes){
  sub.phy<- sub.tree[[j]]
  comm.sub.phy<- comm.full[, sub.phy$tip.label]
  zero.row<- which(rowSums(comm.sub.phy)==0)
  comm.subset<- comm.sub.phy[-zero.row, ]
  dim.comm.subset[j,]<- dim(comm.subset)
}

#vizualizing tribes´ phylogenies
plot(sub.tree[[87]]) #Thomasomyini
plot(sub.tree[[114]])#Akodontini
plot(sub.tree[[20]])$tip.label#Phyllotini
plot(sub.tree[[178]])$tip.label#Oryzomyini

####Akodontini####
comm.sub.phy<- comm.full[,sub.tree[[114]]$tip.label]
zero.row<- which(rowSums(comm.sub.phy)==0)
comm.subset<- comm.sub.phy[-zero.row,]

sub.phy<- sub.tree[[114]]
ape::is.ultrametric(sub.phy)
phangorn::nnls.tree(cophenetic(sub.phy),
                    sub.phy,rooted=TRUE)
ape::is.ultrametric(sub.phy)

runs<- 100
OU.alpha<- c(0, 0.05, 0.25, 1)
comm<- comm.subset
my.species<- as.matrix(colnames(comm))
env.subset<- env.full[rownames(comm.subset), ]
pc.env<- prcomp(env.subset,scale=TRUE)
pc1<- pc.env$x[,1]
esp.subset<- esp.full[rownames(comm.subset),]
plot(esp.subset[,1], esp.subset[,2])

match<- picante::match.phylo.comm(sub.phy, comm.subset)
comm<- match$comm
sub.phy<- match$phy
rownames(comm)==rownames(env.subset)

div<- vegan::renyi(comm,scales=1)
mod<- lm(div~env.subset$Lat)
pred<- predict.lm(mod)
envir<- scales::rescale(as.matrix(env.subset$Lat), c(1,100))
sigma= sd(envir)
ED<- cbind(envir, pred)
theta<- ED[which.max(pred), 1]
root.value<- mean(envir)

test.akodontini<-mcfly::Mcfly(comm = comm, subset = FALSE, occurrence = TRUE, env = envir, 
                              site.coords = esp.full, tree = phy, OU.alpha = OU.alpha, sigma = sigma, theta = theta, 
                              root.value = root.value, runs = runs, ncores = 4, W.r = 0, 
                              output.dir.path = "Akodontini")

save.image("DadosEmpiricos_Akodontini.RData") #results for Akodontini


####Phyllotini####
comm.sub.phy<- comm.full[,sub.tree[[20]]$tip.label]
zero.row<- which(rowSums(comm.sub.phy)==0)
comm.subset<- comm.sub.phy[-zero.row,]
dim(comm.subset)

sub.phy<- sub.tree[[20]]
ape::is.ultrametric(sub.phy)
phangorn::nnls.tree(cophenetic(sub.phy),
                    sub.phy,rooted=TRUE)
ape::is.ultrametric(sub.phy)

#setting parameters
runs<- 100
OU.alpha<- c(0, 0.05, 0.25, 1)
comm<- comm.subset
my.species<- as.matrix(colnames(comm))
env.subset<- env.full[rownames(comm.subset), ]
pc.env<- prcomp(env.subset,scale=TRUE)
pc1<- pc.env$x[, 1]
esp.subset<- esp.full[rownames(comm.subset), ]
plot(esp.subset[, 1], esp.subset[, 2])

match<- picante::match.phylo.comm(sub.phy, comm.subset) #organizing data
comm<- match$comm
sub.phy<- match$phy
rownames(comm)==rownames(env.subset)

div<- vegan::renyi(comm,scales=1) #calculating diversity
mod<- nlm(div ~ env.subset$Temp)
pred<- predict.lm(mod)
envir<- scales::rescale(as.matrix(env.subset$Temp), c(1, 100))
sigma<- sd(envir)
ED<- cbind(envir, pred)
theta<- ED[which.max(pred), 1]
root.value<- mean(envir)


test.phyllotini<-mcfly::Mcfly(comm = comm, subset = FALSE, occurrence = TRUE, env = envir, 
                              site.coords = esp.full, tree = phy, OU.alpha = OU.alpha, sigma = sigma, theta = theta, 
                              root.value = root.value, runs = runs, ncores = 4, W.r = 0, 
                              output.dir.path = "Phyllotini")


####Oryzomyini####
comm.sub.phy<- comm.full[,sub.tree[[178]]$tip.label]
zero.row<- which(rowSums(comm.sub.phy)==0)
comm.subset<- comm.sub.phy[-zero.row, ]
dim(comm.subset)

sub.phy<- sub.tree[[178]]
ape::is.ultrametric(sub.phy)
phangorn::nnls.tree(cophenetic(sub.phy),
                    sub.phy,rooted=TRUE)
ape::is.ultrametric(sub.phy)

#setting parameters
runs<- 100
OU.alpha<- c(0, 0.05, 0.25, 1)
comm<- comm.subset
my.species<- as.matrix(colnames(comm))
env.subset<- env.full[rownames(comm.subset), ]
pc.env<- prcomp(env.subset, scale= TRUE)
pc1<- pc.env$x[, 1]
esp.subset<- esp.full[rownames(comm.subset), ]
plot(esp.subset[,1], esp.subset[,2])

match<- picante::match.phylo.comm(sub.phy,comm.subset)
comm<- match$comm
sub.phy<- match$phy
rownames(comm)==rownames(env.subset)

div<- vegan::renyi(comm,scales=1)
mod<- lm(div ~ env.subset$Temp_Seas)
pred<- predict.lm(mod)
envir<- scales::rescale(as.matrix(env.subset$Temp_Seas), c(1, 100))
sigma<- sd(envir)
ED<- cbind(envir, pred)
theta<- ED[which.max(pred), 1]
root.value<- mean(envir)

test.oryzomyini<-mcfly::Mcfly(comm = comm, subset = FALSE, occurrence = TRUE, env = envir, 
                              site.coords = esp.full, tree = phy, OU.alpha = OU.alpha, sigma = sigma, theta = theta, 
                              root.value = root.value, runs = runs, ncores = 4, W.r = 0, 
                              output.dir.path = "Oryzomyini")

plot(test.oryzomyini$Predicted.entropy.1[,3],test.oryzomyini$Entropy[,1]) #ploting predicted and observed entropy
save.image("DadosEmpiricos_Oryzomyini.RData") #results for Orizomyini

####Thomasomyini####
comm.sub.phy<- comm.full[,sub.tree[[87]]$tip.label]
zero.row<- which(rowSums(comm.sub.phy)==0)
comm.subset<- comm.sub.phy[-zero.row, ]

sub.phy<- sub.tree[[87]]
ape::is.ultrametric(sub.phy)
phangorn::nnls.tree(cophenetic(sub.phy),
                    sub.phy,rooted=TRUE)
ape::is.ultrametric(sub.phy)

#setting parameters
runs<- 100
OU.alpha<- c(0, 0.05, 0.25, 1)
comm<- comm.subset
my.species<- as.matrix(colnames(comm))
env.subset<- env.full[rownames(comm.subset), ]
pc.env<- prcomp(env.subset, scale=TRUE)
summary(pc.env)

biplot(pc.env)
pc1<-pc.env$x[,1]
pc2<-pc.env$x[,2]

esp.subset<-esp.full[rownames(comm.subset), ]
plot(esp.subset[,1], esp.subset[,2])

#organizing data
match<- picante::match.phylo.comm(sub.phy, comm.subset)
comm<- match$comm
sub.phy<- match$phy
rownames(comm)==rownames(env.subset)

div<- vegan::renyi(comm, scales=1)
mod<- lm(div ~ pc2)
pred<- predict.lm(mod)
envir<- scales::rescale(as.matrix(pc2), c(1, 100))
sigma<- sd(envir)
ED<- cbind(envir, pred)
theta<- ED[which.max(pred), 1]
root.value<- mean(envir)

test.thomasomyini<-mcfly::Mcfly(comm = comm, subset = FALSE, occurrence = TRUE, env = envir, 
                                site.coords = esp.full, tree = phy, OU.alpha = OU.alpha, sigma = sigma, theta = theta, 
                                root.value = root.value, runs = runs, ncores = 4, W.r = 0, 
                                output.dir.path = "Thomasomyini")

save.image("DadosEmpiricos_Thomasomyini.RData") #results for Thomasomyini tribe
