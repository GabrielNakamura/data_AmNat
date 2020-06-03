require(mcfly)
source(here::here("R", "functions", "sim.Ncomm_v7.6.1.R")) #reading function to perform analysis with simulated data

####analysis with simulated metacommunities####

###alpha 0####
alpha0.100meta.100runs.noise0.thresh09<- sim.Ncomm(scenario.ID= "species.sorting",
              sim.ID= "alpha0", output.dir.path= "OUTPUT_alpha0",
              Nmeta= 100, Nspp.phy= 100, Nspp.comm= 100, Ncomm= 50,subset= FALSE, occurrence= TRUE, alpha.comm= 0,
              alpha.values= c(0,0.05,0.25,1), noise= 0, W.r= 0, W.threshold= 0.9, n.timestep= 50, 
              runs= 100, ncores= 4)
save.image("SimulatedData_mcfly_Mar2020_pa_1_100meta_thresh09.RData") #results

####alpha 0.05####
alpha005.100meta.100runs.noise0.thresh09<- sim.Ncomm(scenario.ID= "species.sorting",
              sim.ID= "alpha005", output.dir.path= "OUTPUT_alpha005", 
              Nmeta= 100, Nspp.phy= 100, Nspp.comm= 100, Ncomm= 50, subset= TRUE, occurrence= TRUE, alpha.comm= 0.05, 
              alpha.values= c(0, 0.05, 0.25, 1), noise= 0, W.r= 0, W.threshold= 0.9, n.timestep= 50, 
              runs= 100, ncores= 4)
save.image("SimulatedData_mcfly_Mar2020_pa_4_100meta_thresh09.RData") #results

####alpha 0.25####
alpha025.100meta.100runs.subset20.noise0.thresh09<- sim.Ncomm(scenario.ID= "species.sorting", 
              sim.ID= "alpha025", output.dir.path= "OUTPUT_alpha025", 
              Nmeta= 100, Nspp.phy= 100, Nspp.comm= 100, Ncomm= 50, subset= TRUE, occurrence= TRUE, 
              alpha.comm= 0.25, alpha.values= c(0, 0.05, 0.25, 1), noise= 0, W.r= 0, W.threshold= 0.9, 
              n.timestep= 50, runs= 100, ncores= 4)
save.image("SimulatedData_mcfly_Mar2020_pa_7_100meta_thresh09.RData") #results

####alpha 1#####
alpha1.100meta.100runs.subset20.noise0.thresh09<-sim.Ncomm(scenario.ID="species.sorting",
              sim.ID="alpha1",output.dir.path="OUTPUT_alpha1",
              Nmeta=100,Nspp.phy=100,Nspp.comm=100,Ncomm=50,subset=TRUE,occurrence=TRUE,alpha.comm=1,
              alpha.values=c(0,0.05,0.25,1),noise=0,W.r=0,W.threshold=0.9,n.timestep=50,
              runs=100,ncores=4)
save.image("SimulatedData_mcfly_Mar2020_pa_10_100meta_thresh09.RData") #results
