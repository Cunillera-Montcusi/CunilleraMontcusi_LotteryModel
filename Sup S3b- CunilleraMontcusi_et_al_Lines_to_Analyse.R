##############################################################################################################
##############################################################################################################
###################### METACOMMUNITIES GENERATION  ###########################################################
##############################################################################################################

# Install packages required to run the functions if needed

#install.packages("vegan")
#install.packages("igraph")
#install.packages("sna")
#install.packages("dummies")

library(vegan)
library(dummies)

# Initial steps before metacommunity splitting
# Metacommunity marix to start simulation
No<-rep(20,200) # 20 species havin 200 individuals each
metacom<-matrix(rep(No),length(No),542) # 542 is the number of nodes (water bodies) of the metacomunity.

# output for tracking simulation convegence
out.t<-matrix(c(1,0,0),1,3)

# List of graphs to include in the function
#detach("package:sna", unload = TRUE)
#detach("package:igraph", unload = TRUE)
library(igraph)
# Chargin UTM coordinates of water bodies network
xy.albera <- read.csv2(file = "input/xy.albera.csv", row.names = 1)
as.matrix(dist(xy.albera[,1:2]))->D_ALBERA

m250 <- ifelse(D_ALBERA>250,0,1)
m500 <- ifelse(D_ALBERA>500,0,1)
m1000 <- ifelse(D_ALBERA>1000,0,1)
m1500 <- ifelse(D_ALBERA>1500,0,1)
m2000 <- ifelse(D_ALBERA>2000,0,1)
m2500 <- ifelse(D_ALBERA>2500,0,1)
m3000 <- ifelse(D_ALBERA>3000,0,1)
mpercol <- ifelse(D_ALBERA>3849.61,0,1)
m5000 <- ifelse(D_ALBERA>5000,0,1)

# Graph generation corresponiding to each used network
LG_albera250 <- graph.adjacency(m250, mode="undirected", diag=FALSE)
LG_albera500 <- graph.adjacency(m500, mode="undirected", diag=FALSE)
LG_albera1000 <- graph.adjacency(m1000, mode="undirected", diag=FALSE)
LG_albera1500 <- graph.adjacency(m1500, mode="undirected", diag=FALSE)
LG_albera2000 <- graph.adjacency(m2000, mode="undirected", diag=FALSE)
LG_albera2500 <- graph.adjacency(m2500, mode="undirected", diag=FALSE)
LG_albera3000 <- graph.adjacency(m3000, mode="undirected", diag=FALSE)
LG_alberapercol <- graph.adjacency(mpercol, mode="undirected", diag=FALSE)
LG_albera5000 <- graph.adjacency(m5000, mode="undirected", diag=FALSE)

# Joint all created graphs into a list 
all_graphs <- list(LG_albera250, LG_albera500, LG_albera1000, LG_albera1500, LG_albera2000,
                   LG_albera2500, LG_albera3000, LG_alberapercol, LG_albera5000)
length(all_graphs) # Check number of items

##############################################################################################################
##### MANUSCRIPT FIGURE 2 FIRST STEP #########################################################################
##############################################################################################################

#### METACOMMUNITY SPLITTING  ###############################################################################
##############################################################################################################

##################################################################################################################################
#### ALPHA of 0.1 ##########################

all_alphas0.1 <- rep(0.1, times=length(all_graphs)) # list of alphas per each dispersal ability group

# split matrix in two matrix. Each "dispersal group" have a fraction (alpha) of individuals realted to a
#  resource shared with other trophic groups. All individuals in the "shared pool" follow a neutral
#  dynamic in which individuals of any species can be replaced by individuals of any other. In addition,
#  the other fraction of individulas are related to resources exclusevely exploted by each dispersal group (1-alpha)
mm0.1<- split.metacomm.alpha(grafos = all_graphs ,metacom = metacom, alphas = all_alphas0.1) # Use the alphas that represent the set of scenarios to be considered

##################################################################################################################################
#### ALPHA of 0.3 ##########################

all_alphas0.3 <- rep(0.3, times=length(all_graphs))

# split matrix in two matrix. Each "dispersal group" have a fraction (alpha) of individuals realted to a
#  resource shared with other trophic groups. All individuals in the "shared pool" follow a neutral
#  dynamic in which individuals of any species can be replaced by individuals of any other. In addition,
#  the other fraction of individulas are related to resources exclusevely exploted by each dispersal group (1-alpha)
mm0.3<- split.metacomm.alpha(grafos = all_graphs ,metacom = metacom, alphas = all_alphas0.3) # Use the alphas that represent the set of scenarios to be considered

##################################################################################################################################
#### ALPHA of 0.5 ##########################

all_alphas0.5 <- rep(0.5, times=length(all_graphs))

# split matrix in two matrix. Each "dispersal group" have a fraction (alpha) of individuals realted to a
#  resource shared with other trophic groups. All individuals in the "shared pool" follow a neutral
#  dynamic in which individuals of any species can be replaced by individuals of any other. In addition,
#  the other fraction of individulas are related to resources exclusevely exploted by each dispersal group (1-alpha)
mm0.5<- split.metacomm.alpha(grafos = all_graphs ,metacom = metacom, alphas = all_alphas0.5) # Use the alphas that represent the set of scenarios to be considered

##################################################################################################################################
#### ALPHA of 0.7 ##########################

all_alphas0.7 <- rep(0.7, times=length(all_graphs))

# split matrix in two matrix. Each "dispersal group" have a fraction (alpha) of individuals realted to a
#  resource shared with other trophic groups. All individuals in the "shared pool" follow a neutral
#  dynamic in which individuals of any species can be replaced by individuals of any other. In addition,
#  the other fraction of individulas are related to resources exclusevely exploted by each dispersal group (1-alpha)
mm0.7<- split.metacomm.alpha(grafos = all_graphs ,metacom = metacom, alphas = all_alphas0.7) # Use the alphas that represent the set of scenarios to be considered

##################################################################################################################################
#### ALPHA of 0.9 ##########################

all_alphas0.9 <- rep(0.9, times=length(all_graphs))

# split matrix in two matrix. Each "dispersal group" have a fraction (alpha) of individuals realted to a
#  resource shared with other trophic groups. All individuals in the "shared pool" follow a neutral
#  dynamic in which individuals of any species can be replaced by individuals of any other. In addition,
#  the other fraction of individulas are related to resources exclusevely exploted by each dispersal group (1-alpha)
mm0.9<- split.metacomm.alpha(grafos = all_graphs ,metacom = metacom, alphas = all_alphas0.9) # Use the alphas that represent the set of scenarios to be considered



##############################################################################################################
##### MANUSCRIPT FIGURE 2 SECOND STEP ########################################################################
##############################################################################################################

#### METACOMMUNITY CREATION  #################################################################################
##############################################################################################################

##################################################################################################################################
#### ALPHA of 0.1 ##########################

# Generate the entrada list which will have all the needed entering items
entrada0.1<-list("control.S"=out.t,"Spp.groups.dummy"=0, metacom.s=mm0.1[[1]], metacom.e=mm0.1[[2]], metacom=metacom)

# Alpha of 0.1 (0.1% of the community will be shared) 
for (i in 1:25) {
  meta.comm_neutral_asimetrica( grafos= all_graphs ,            # list with all graphs 
                                it= 400 ,                       # number of iterations
                                metacom.e= entrada0.1[[4]],        # metacom exlusive 
                                metacom.s= entrada0.1[[3]],        # metacom shared
                                m.in= 0.01,                     # migraton in 
                                m.out= 0.01,                    # migration out 
                                n= 50,                          # number of individuals removed
                                out.t= entrada0.1[[1]] ,           # out where richness is saved - CONTROL
                                cada= 50 ,                      # every time data will be saved
                                alphas= all_alphas0.1)-> entrada0.1  # a vector with all alpha values (one by each graphs)
  save.image("output/entrada01.RData") # For each for we save the environment 
}

##################################################################################################################################
#### ALPHA of 0.3 ##########################

# Generate the entrada list which will have all the needed entering items
entrada0.3<-list("control.S"=out.t,"Spp.groups.dummy"=0, metacom.s=mm0.3[[1]], metacom.e=mm0.3[[2]], metacom=metacom)

# Alpha of 0.3 (30% of the community will be shared) 
for (i in 1:25) {
  meta.comm_neutral_asimetrica( grafos= all_graphs ,            # list with all graphs 
                                it= 400 ,                       # number of iterations
                                metacom.e= entrada0.3[[4]],        # metacom exlusive 
                                metacom.s= entrada0.3[[3]],        # metacom shared
                                m.in= 0.01,                     # migraton in 
                                m.out= 0.01,                    # migration out 
                                n= 50,                          # number of individuals removed
                                out.t= entrada0.3[[1]] ,           # out where richness is saved - CONTROL
                                cada= 50 ,                      # every time data will be saved
                                alphas= all_alphas0.3)-> entrada0.3  # a vector with all alpha values (one by each graphs)
  save.image("output/entrada03.RData") # For each for we save the environment 
}

##################################################################################################################################
#### ALPHA of 0.5 ##########################

# Generate the entrada list which will have all the needed entering items
entrada0.5<-list("control.S"=out.t,"Spp.groups.dummy"=0, metacom.s=mm0.5[[1]], metacom.e=mm0.5[[2]], metacom=metacom)

# Alpha of 0.5 (50% of the community will be shared) 
for (i in 1:25) {
  meta.comm_neutral_asimetrica( grafos= all_graphs ,            # list with all graphs 
                                it= 400 ,                       # number of iterations
                                metacom.e= entrada0.5[[4]],        # metacom exlusive 
                                metacom.s= entrada0.5[[3]],        # metacom shared
                                m.in= 0.01,                     # migraton in 
                                m.out= 0.01,                    # migration out 
                                n= 50,                          # number of individuals removed
                                out.t= entrada0.5[[1]] ,           # out where richness is saved - CONTROL
                                cada= 50 ,                      # every time data will be saved
                                alphas= all_alphas0.5)-> entrada0.5  # a vector with all alpha values (one by each graphs)
  save.image("output/entrada05.RData") # For each for we save the environment 
}

##################################################################################################################################
#### ALPHA of 0.7 ##########################

# Generate the entrada list which will have all the needed entering items
entrada0.7<-list("control.S"=out.t,"Spp.groups.dummy"=0, metacom.s=mm0.7[[1]], metacom.e=mm0.7[[2]], metacom=metacom)

# Alpha of 0.7 (70% of the community will be shared) 
for (i in 1:25) {
  meta.comm_neutral_asimetrica( grafos= all_graphs ,            # list with all graphs 
                                it= 400 ,                       # number of iterations
                                metacom.e= entrada0.7[[4]],        # metacom exlusive 
                                metacom.s= entrada0.7[[3]],        # metacom shared
                                m.in= 0.01,                     # migraton in 
                                m.out= 0.01,                    # migration out 
                                n= 50,                          # number of individuals removed
                                out.t= entrada0.7[[1]] ,           # out where richness is saved - CONTROL
                                cada= 50 ,                      # every time data will be saved
                                alphas= all_alphas0.7)-> entrada0.7  # a vector with all alpha values (one by each graphs)
  save.image("output/entrada07.RData") # For each for we save the environment 
}

##################################################################################################################################
#### ALPHA of 0.9 ##########################

# Generate the entrada list which will have all the needed entering items
entrada0.9<-list("control.S"=out.t,"Spp.groups.dummy"=0, metacom.s=mm0.9[[1]], metacom.e=mm0.9[[2]], metacom=metacom)

# Alpha of 0.9 (90% of the community will be shared) 
for (i in 1:25) {
  meta.comm_neutral_asimetrica( grafos= all_graphs ,            # list with all graphs 
                                it= 400 ,                       # number of iterations
                                metacom.e= entrada0.9[[4]],        # metacom exlusive 
                                metacom.s= entrada0.9[[3]],        # metacom shared
                                m.in= 0.01,                     # migraton in 
                                m.out= 0.01,                    # migration out 
                                n= 50,                          # number of individuals removed
                                out.t= entrada0.9[[1]] ,           # out where richness is saved - CONTROL
                                cada= 50 ,                      # every time data will be saved
                                alphas= all_alphas0.9)-> entrada0.9  # a vector with all alpha values (one by each graphs)
  save.image("output/entrada09.RData") # For each for we save the environment 
}


##############################################################################################################
##### MANUSCRIPT FIGURE 2 THIRD STEP #########################################################################
##############################################################################################################

#### FIRES SIMULATION  #######################################################################################
##############################################################################################################

# Load the real wildfire data 

load("Original_Fire.RData")
Original_FIRE

#Generate the 400 fire scenarios 
fake.wildfires<-rangos(XY = xy.albera, fire0 = 467, niveles = 20, efectividad = 20)
colnames(fake.wildfires) # Check all values

##############################################################################################################
##### MANUSCRIPT FIGURE 2 FOURTH STEP ########################################################################
##############################################################################################################

#### BURNING AND RECOLONIZATION DYNAMICS  ####################################################################
##############################################################################################################

#____________________________________________________________________________________________________________________________________________________________________
#### REAL WILDFIRE  ####################################################################

# Generate an output list with 9 "items"
lista.n <- list("primer","segon","tercer","quaert","cinque","sise","sete","buite","nove")

### Aplha 0.1 Original_FIRE____________________________________________________________________________________________________________________________####
# T1 ____________________________________________________________________________________________________________________________
Orig_rec_ENTRADA_0.1 <- burning_recolonization_function.2(M =Original_FIRE, 
                                                             metacom =   entrada0.1[[5]],
                                                             metacom.s = entrada0.1[[3]],
                                                             metacom.e = entrada0.1[[4]],
                                                             migra =0.01, grafos = all_graphs, it = 1, 
                                                             alphas =  rep(0.1, times=length(all_graphs)),
                                                             dummy.pool =entrada0.1[[2]],lista.n.grupos = lista.n)

save.image("output/Orig_rec_ENTRADA_0.1_100.RData")

# T100 ____________________________________________________________________________________________________________________________
Orig_rec_ENTRADA_0.1_100 <- burning_recolonization_function.2(M =Original_FIRE, 
                                                                 metacom =   entrada0.1[[5]],
                                                                 metacom.s = entrada0.1[[3]],
                                                                 metacom.e = entrada0.1[[4]],
                                                                 migra =0.01, grafos = all_graphs, it = 100, 
                                                                 alphas =  rep(0.1, times=length(all_graphs)),
                                                                 dummy.pool =entrada0.1[[2]],lista.n.grupos = lista.n)
save.image("output/Orig_rec_ENTRADA_0.1_100.RData")

### Aplha 0.3 Original_FIRE____________________________________________________________________________________________________________________________####
# T1 ____________________________________________________________________________________________________________________________
Orig_rec_ENTRADA_0.3 <- burning_recolonization_function.2(M =Original_FIRE, 
                                                             metacom =   entrada0.3[[5]],
                                                             metacom.s = entrada0.3[[3]],
                                                             metacom.e = entrada0.3[[4]],
                                                             migra =0.01, grafos = all_graphs, it = 1, 
                                                             alphas =  rep(0.3, times=length(all_graphs)),
                                                             dummy.pool =entrada0.3[[2]],lista.n.grupos = lista.n)
save.image("output/Orig_rec_ENTRADA_0.3.RData")

# T100 ____________________________________________________________________________________________________________________________
Orig_rec_ENTRADA_0.3_100 <- burning_recolonization_function.2(M =Original_FIRE, 
                                                              metacom =   entrada0.3[[5]],
                                                              metacom.s = entrada0.3[[3]],
                                                              metacom.e = entrada0.3[[4]],
                                                              migra =0.01, grafos = all_graphs, it = 100, 
                                                              alphas =  rep(0.3, times=length(all_graphs)),
                                                              dummy.pool =entrada0.3[[2]],lista.n.grupos = lista.n)
save.image("output/Orig_rec_ENTRADA_0.3_100.RData")

### Aplha 0.5 Original_FIRE____________________________________________________________________________________________________________________________####
# T1 ____________________________________________________________________________________________________________________________
Orig_rec_ENTRADA_0.5 <- burning_recolonization_function.2(M =Original_FIRE, 
                                                          metacom =   entrada0.5[[5]],
                                                          metacom.s = entrada0.5[[3]],
                                                          metacom.e = entrada0.5[[4]],
                                                          migra =0.01, grafos = all_graphs, it = 1, 
                                                          alphas =  rep(0.5, times=length(all_graphs)),
                                                          dummy.pool =entrada0.5[[2]],lista.n.grupos = lista.n)
save.image("output/Orig_rec_ENTRADA_0.5.RData")
# T100 ____________________________________________________________________________________________________________________________
Orig_rec_ENTRADA_0.5_100 <- burning_recolonization_function.2(M =Original_FIRE, 
                                                              metacom =   entrada0.5[[5]],
                                                              metacom.s = entrada0.5[[3]],
                                                              metacom.e = entrada0.5[[4]],
                                                              migra =0.01, grafos = all_graphs, it = 100, 
                                                              alphas =  rep(0.5, times=length(all_graphs)),
                                                              dummy.pool =entrada0.5[[2]],lista.n.grupos = lista.n)
save.image("output/Orig_rec_ENTRADA_0.5_100.RData")

### Aplha 0.7 Original_FIRE____________________________________________________________________________________________________________________________####
# T1 ____________________________________________________________________________________________________________________________
Orig_rec_ENTRADA_0.7 <- burning_recolonization_function.2(M =Original_FIRE, 
                                                          metacom =   entrada0.7[[5]],
                                                          metacom.s = entrada0.7[[3]],
                                                          metacom.e = entrada0.7[[4]],
                                                          migra =0.01, grafos = all_graphs, it = 1, 
                                                          alphas =  rep(0.7, times=length(all_graphs)),
                                                          dummy.pool =entrada0.7[[2]],lista.n.grupos = lista.n)
save.image("output/Orig_rec_ENTRADA_0.7.RData")
# T100 ____________________________________________________________________________________________________________________________
Orig_rec_ENTRADA_0.7_100 <- burning_recolonization_function.2(M =Original_FIRE, 
                                                              metacom =   entrada0.7[[5]],
                                                              metacom.s = entrada0.7[[3]],
                                                              metacom.e = entrada0.7[[4]],
                                                              migra =0.01, grafos = all_graphs, it = 100, 
                                                              alphas =  rep(0.7, times=length(all_graphs)),
                                                              dummy.pool =entrada0.7[[2]],lista.n.grupos = lista.n)
save.image("output/Orig_rec_ENTRADA_0.7_100.RData")

### Aplha 0.9 Original_FIRE____________________________________________________________________________________________________________________________####
# T1 ____________________________________________________________________________________________________________________________
Orig_rec_ENTRADA_0.9 <- burning_recolonization_function.2(M =Original_FIRE, 
                                                          metacom =   entrada0.9[[5]],
                                                          metacom.s = entrada0.9[[3]],
                                                          metacom.e = entrada0.9[[4]],
                                                          migra =0.01, grafos = all_graphs, it = 1, 
                                                          alphas =  rep(0.9, times=length(all_graphs)),
                                                          dummy.pool =entrada0.9[[2]],lista.n.grupos = lista.n)
save.image("output/Orig_rec_ENTRADA_0.9.RData")
# T100 ____________________________________________________________________________________________________________________________
Orig_rec_ENTRADA_0.9_100 <- burning_recolonization_function.2(M =Original_FIRE, 
                                                              metacom =   entrada0.9[[5]],
                                                              metacom.s = entrada0.9[[3]],
                                                              metacom.e = entrada0.9[[4]],
                                                              migra =0.01, grafos = all_graphs, it = 100, 
                                                              alphas =  rep(0.9, times=length(all_graphs)),
                                                              dummy.pool =entrada0.9[[2]],lista.n.grupos = lista.n)
save.image("output/Orig_rec_ENTRADA_0.9_100.RData")

#____________________________________________________________________________________________________________________________________________________________________
#### SIMULATED WILDFIRE SCENARIOS  ####################################################################

# Generate an output list with 9 "items"
lista.n <- list("primer","segon","tercer","quaert","cinque","sise","sete","buite","nove")

### Aplha 0.1 fake.fires____________________________________________________________________________________________________________________________####
load("input/0.1_inputing.RData")

# T1 ____________________________________________________________________________________________________________________________

Burn_rec_ENTRADA_0.1 <- burning_recolonization_function.2(M =fake.wildfires, 
                                                                 metacom =   entrada0.1[[5]],
                                                                 metacom.s = entrada0.1[[3]],
                                                                 metacom.e = entrada0.1[[4]],
                                                                 migra =0.01, grafos = all_graphs, it = 1, 
                                                                 alphas =  rep(0.1, times=length(all_graphs)),
                                                                 dummy.pool =entrada0.1[[2]],lista.n.grupos = lista.n)
save.image("output/Burn_rec_ENTRADA_0.1.Rdata")

# T100 ____________________________________________________________________________________________________________________________

Burn_rec_ENTRADA_0.1_100 <- burning_recolonization_function.2(M =fake.wildfires, 
                                                                 metacom =   entrada0.1[[5]],
                                                                 metacom.s = entrada0.1[[3]],
                                                                 metacom.e = entrada0.1[[4]],
                                                                 migra =0.01, grafos = all_graphs, it = 100, 
                                                                 alphas =  rep(0.1, times=length(all_graphs)),
                                                                 dummy.pool =entrada0.1[[2]],lista.n.grupos = lista.n)
save.image("output/Burn_rec_ENTRADA_0.1_100.RData")



### Aplha 0.3 fake.fires____________________________________________________________________________________________________________________________####
load("input/0.3_inputing.RData")

# T1 ____________________________________________________________________________________________________________________________

Burn_rec_ENTRADA_0.3 <- burning_recolonization_function.2(M =fake.wildfires, 
                                                                 metacom =   entrada0.3[[5]],
                                                                 metacom.s = entrada0.3[[3]],
                                                                 metacom.e = entrada0.3[[4]],
                                                                 migra =0.01, grafos = all_graphs, it = 1, 
                                                                 alphas =  rep(0.3, times=length(all_graphs)),
                                                                 dummy.pool =entrada0.3[[2]],lista.n.grupos = lista.n)
save.image("output/Burn_rec_ENTRADA_0.3.RData")

# T100 ____________________________________________________________________________________________________________________________

Burn_rec_ENTRADA_0.3_100 <- burning_recolonization_function.2(M =fake.wildfires, 
                                                                 metacom =   entrada0.3[[5]],
                                                                 metacom.s = entrada0.3[[3]],
                                                                 metacom.e = entrada0.3[[4]],
                                                                 migra =0.01, grafos = all_graphs, it = 100, 
                                                                 alphas =  rep(0.3, times=length(all_graphs)),
                                                                 dummy.pool =entrada0.3[[2]],lista.n.grupos = lista.n)
save.image("output/Burn_rec_ENTRADA_0.3_100_05.RData")



### Aplha 0.5 fake.fires____________________________________________________________________________________________________________________________####
load("input/0.5_inputing.RData")

# T1 ____________________________________________________________________________________________________________________________

Burn_rec_ENTRADA_0.5 <- burning_recolonization_function.2(M =fake.wildfires, 
                                                                 metacom =   entrada0.5[[5]],
                                                                 metacom.s = entrada0.5[[3]],
                                                                 metacom.e = entrada0.5[[4]],
                                                                 migra =0.01, grafos = all_graphs, it = 1, 
                                                                 alphas =  rep(0.5, times=length(all_graphs)),
                                                                 dummy.pool =entrada0.5[[2]],lista.n.grupos = lista.n)
save.image("output/Burn_rec_ENTRADA_0.5.RData")

# T100 ____________________________________________________________________________________________________________________________

Burn_rec_ENTRADA_0.5_100 <- burning_recolonization_function.2(M =fake.wildfires, 
                                                                 metacom =   entrada0.5[[5]],
                                                                 metacom.s = entrada0.5[[3]],
                                                                 metacom.e = entrada0.5[[4]],
                                                                 migra =0.01, grafos = all_graphs, it = 100, 
                                                                 alphas =  rep(0.5, times=length(all_graphs)),
                                                                 dummy.pool =entrada0.5[[2]],lista.n.grupos = lista.n)
save.image("output/Burn_rec_ENTRADA_0.5_100.RData")



### Aplha 0.7 fake.fires____________________________________________________________________________________________________________________________####
load("input/0.7_inputing.RData")

# T1 ____________________________________________________________________________________________________________________________

Burn_rec_ENTRADA_0.7 <- burning_recolonization_function.2(M =fake.wildfires, 
                                                                 metacom =   entrada0.7[[5]],
                                                                 metacom.s = entrada0.7[[3]],
                                                                 metacom.e = entrada0.7[[4]],
                                                                 migra =0.01, grafos = all_graphs, it = 1, 
                                                                 alphas =  rep(0.7, times=length(all_graphs)),
                                                                 dummy.pool =entrada0.7[[2]],lista.n.grupos = lista.n)
save.image("output/Burn_rec_ENTRADA_0.7.RData")




# T100 ____________________________________________________________________________________________________________________________

Burn_rec_ENTRADA_0.7_100 <- burning_recolonization_function.2(M =fake.wildfires[,c(1:3,4)], 
                                                                 metacom =   entrada0.7[[5]],
                                                                 metacom.s = entrada0.7[[3]],
                                                                 metacom.e = entrada0.7[[4]],
                                                                 migra =0.01, grafos = all_graphs, it = 100, 
                                                                 alphas =  rep(0.7, times=length(all_graphs)),
                                                                 dummy.pool =entrada0.7[[2]],lista.n.grupos = lista.n)
save.image("output/Burn_rec_ENTRADA_0.7_100.RData")



### Aplha 0.9 fake.fires____________________________________________________________________________________________________________________________####
load("input/0.9_inputing.RData")

# T1 ____________________________________________________________________________________________________________________________

Burn_rec_ENTRADA_0.9 <- burning_recolonization_function.2(M =fake.wildfires, 
                                                                 metacom =   entrada0.9[[5]],
                                                                 metacom.s = entrada0.9[[3]],
                                                                 metacom.e = entrada0.9[[4]],
                                                                 migra =0.01, grafos = all_graphs, it = 1, 
                                                                 alphas =  rep(0.9, times=length(all_graphs)),
                                                                 dummy.pool =entrada0.9[[2]],lista.n.grupos = lista.n)
save.image("output/Burn_rec_ENTRADA_0.9.RData")

# T100 ____________________________________________________________________________________________________________________________

Burn_rec_ENTRADA_0.9_100 <- burning_recolonization_function.2(M =fake.wildfires, 
                                                                 metacom =   entrada0.9[[5]],
                                                                 metacom.s = entrada0.9[[3]],
                                                                 metacom.e = entrada0.9[[4]],
                                                                 migra =0.01, grafos = all_graphs, it = 100, 
                                                                 alphas =  rep(0.9, times=length(all_graphs)),
                                                                 dummy.pool =entrada0.9[[2]],lista.n.grupos = lista.n)
save.image("output/Burn_rec_ENTRADA_0.9_100.RData")


#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
###################### RESULTS AND PLOTS REPRESENTATION #####################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

#############################################################################################################################
#############################################################################################################################
# Original wildfire_______________________________
#############################################################################################################################
#############################################################################################################################

#load the data related with the real wildfire
load("C:/Users/Cunilleramontcusi/Desktop/Burn/Original/Orig_rec_ENTRADA_0.1_100.RData")
load("C:/Users/Cunilleramontcusi/Desktop/Burn/Original/Orig_rec_ENTRADA_0.3_100.RData")
load("C:/Users/Cunilleramontcusi/Desktop/Burn/Original/Orig_rec_ENTRADA_0.5_100.RData")
load("C:/Users/Cunilleramontcusi/Desktop/Burn/Original/Orig_rec_ENTRADA_0.7_100.RData")
load("C:/Users/Cunilleramontcusi/Desktop/Burn/Original/Orig_rec_ENTRADA_0.9_100.RData")

# Cunillera_palette. 
source("C:/Users/Cunilleramontcusi/Dropbox/DAVID DOC/LLAM al DIA/CUNILLERA_palette.R")
colorts<- CUNILLERA_pal("grey")(9)

# PLot functions for ORIGINAL fire values
plot_Orig_RCperc <- function(Orig_rec, Pre_fire_Orig, títol){
  Dispersal_abilities <- c(250,500,1000,1500,2000,2500,3000,3842,5000)
  # Original_______________________________________________________________ 
  spp_disponibles <- NULL
  pre_spp <- NULL
  for (e in 1:9) {
    spp_disponibles[e] <-Orig_rec$Prop.spp.disponibles[[e]][,14]
    pre_spp[e] <- Pre_fire_Orig$Prop.spp.disponibles[[e]][,14]
    #pre_spp[e] <- mean(apply(ifelse(as.matrix(Pre_fire_Orig[[5]][which(Pre_fire_Orig[[2]][,e]==1),Original_FIRE[,3]])>0,1,0),2,sum))
  }
  spp_disponibles <- (pre_spp*100)/spp_disponibles
  spp_disponibles<- ifelse(spp_disponibles>100,100,spp_disponibles)
  # PLot line
  plot(Dispersal_abilities, spp_disponibles,
       ylim = c(0,105), xlim = c(0,5000), type="l",col="darkblue", lwd=2, lty=3,cex=2,cex.axis=1.5,
       ylab = "", xlab= "",main = "")
  # Plot points
  par(new=T)
  plot(Dispersal_abilities, spp_disponibles,
       ylim = c(0,105), xlim = c(0,5000), type="p", col="black",bg=colorts, 
       ylab = "DR", xlab= "Dispersal abilities",main = títol,
       pch=21, cex=2,cex.axis=1.5)
  
  #SD _____________________
  sd_spp_disponibles <- NULL
  sd_pre_spp <- NULL
  for (e in 1:9) {
    sd_spp_disponibles[e] <-Orig_rec$sd.spp.disponibles[[e]][,14]
    #sd_pre_spp[e] <- sd(apply(ifelse(as.matrix(Pre_fire_Orig[[e]][which(Pre_fire_Orig[[2]][,1]==e),Original_FIRE[,3]])>0,1,0),2,sum))
    sd_pre_spp[e] <- Pre_fire_Orig$sd.spp.disponibles[[e]][,14]
    sd_spp_disponibles[e] <- mean(sd_spp_disponibles[e],sd_pre_spp[e])
  }
  
  arrows(x0 = Dispersal_abilities, x1 =Dispersal_abilities,
         y0 =spp_disponibles, y1= ifelse(spp_disponibles-sd_spp_disponibles<0,0,spp_disponibles-sd_spp_disponibles),
         code=2, angle = 90, length = 0.1, col="grey30")
  arrows(x0 = Dispersal_abilities, x1 =Dispersal_abilities,
         y0 =spp_disponibles, y1= spp_disponibles+sd_spp_disponibles,
         code=2, angle = 90, length = 0.1, col="grey30")
  
}

plot_RR_iter <- function(Orig_rec_100, títol){
  Dispersal_abilities <- c(250,500,1000,1500,2000,2500,3000,3842,5000)
  #par(mfrow=c(2,1))
  # Iterations_______________________________________________________________ 
  iter_disponibles <- NULL
  for (e in 1:9) {
    iter_disponibles[e] <-Orig_rec_100$Iteraciones_para_fill[[e]][,14]  
  }
  iter_disponibles <- (100/iter_disponibles)/100
  barCenters <- barplot(iter_disponibles,1,
                        ylim = c(0,1) , col=colorts,border = "black",main = títol,cex.names=1.2,
                        ylab = "RR", xlab= "Dispersal abilities",cex.axis=1.5,cex=2)
  box(which = "plot", lty = "solid")
  text(barCenters,par("usr")[3]-0.015, srt = 60, adj= 1, xpd = TRUE, labels =Dispersal_abilities , cex=1.4)
  
  sd_iter_disponibles <- NULL
  for (e in 1:9) {
    sd_iter_disponibles[e] <-Orig_rec_100$sd.iteraciones.fill[[e]][,14]  
  }
  sd_iter_disponibles <- (100/sd_iter_disponibles)/100
  arrows(x0 =barCenters, x1 =barCenters,
         y0 =iter_disponibles, y1=ifelse(iter_disponibles+sd_iter_disponibles>100,95,iter_disponibles+sd_iter_disponibles),
         angle = 90, length = 0.1)
  #abline(0,0, col="black")  
}

#Plot a common png image for all alpha values (Supplementary S3)
png(filename = "original_fires_outputs.png",width = 18000,height =14000 ,units = "px",res = 1000)
par(mfrow=c(3,4),cex.main=2, las=1)
plot_RR_iter(Orig_rec =Orig_rec_ENTRADA_0.1_100 , títol = "A-0.1")
plot_Orig_RCperc(Orig_rec =Orig_rec_ENTRADA_0.1_100 , Pre_fire_Orig=Orig_rec_ENTRADA_0.1, títol = "B-0.1")
par(las=1)
plot_RR_iter(Orig_rec =Orig_rec_ENTRADA_0.3_100 , títol = "A-0.3")
plot_Orig_RCperc(Orig_rec =Orig_rec_ENTRADA_0.3_100 , Pre_fire_Orig=Orig_rec_ENTRADA_0.3, títol = "B-0.3")
par(las=1)
plot_RR_iter(Orig_rec =Orig_rec_ENTRADA_0.5_100 , títol = "A-0.5")
plot_Orig_RCperc(Orig_rec =Orig_rec_ENTRADA_0.5_100 , Pre_fire_Orig=Orig_rec_ENTRADA_0.5, títol = "B-0.5")
par(las=1)
plot_RR_iter(Orig_rec =Orig_rec_ENTRADA_0.7_100 , títol = "A-0.7")
plot_Orig_RCperc(Orig_rec =Orig_rec_ENTRADA_0.7_100 , Pre_fire_Orig=Orig_rec_ENTRADA_0.7, títol = "B-0.7")
par(las=1)
plot_RR_iter(Orig_rec =Orig_rec_ENTRADA_0.9_100 , títol = "A-0.9")
plot_Orig_RCperc(Orig_rec =Orig_rec_ENTRADA_0.9_100 , Pre_fire_Orig=Orig_rec_ENTRADA_0.9, títol = "B-0.9")
dev.off()

#Plot a common png image for alpha=0.5 values (Figure 4)
png(filename = "original_0.5_outputs.png",width = 15000,height = 7000,units = "px",res = 1000)
par(mfrow=c(1,2),cex.main=2, las=1, font=1, font.lab=2,cex.lab=1.2)
plot_RR_iter(Orig_rec =Orig_rec_ENTRADA_0.5_100, títol = "")
plot_Orig_RCperc(Orig_rec =Orig_rec_ENTRADA_0.5_100 , Pre_fire_Orig=Orig_rec_ENTRADA_0.5, títol="")
dev.off()

#############################################################################################################################
#############################################################################################################################
# Simulated wildfires results      _________________________

# T1 (one iteration after wildfire)_________________________
#############################################################################################################################
#############################################################################################################################

#Load data 
load("C:/Users/Cunilleramontcusi/Desktop/Burn/T1/Burn_rec_ENTRADA_0.1.RData")
load("C:/Users/Cunilleramontcusi/Desktop/Burn/T1/Burn_rec_ENTRADA_0.3.RData")
load("C:/Users/Cunilleramontcusi/Desktop/Burn/T1/Burn_rec_ENTRADA_0.5.RData")
load("C:/Users/Cunilleramontcusi/Desktop/Burn/T1/Burn_rec_ENTRADA_0.7.RData")
load("C:/Users/Cunilleramontcusi/Desktop/Burn/T1/Burn_rec_ENTRADA_0.9.RData")

Combined_article_plot <- function(Post_Val_Matrix, file_name, iterations ){
  plot_Post_fire_1 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[1]][[1]],title = iterations,subtitle = "250 m", max = 22,is.7 =1 )  
  plot_Post_fire_2 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[1]][[2]],iterations,"500 m", max = 22,is.7 =2 )  
  plot_Post_fire_3 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[1]][[3]],iterations,"1000 m", max = 22,is.7 =3 )  
  plot_Post_fire_4 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[1]][[4]],iterations,"1500 m", max = 22,is.7 =4 )  
  plot_Post_fire_5 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[1]][[5]],iterations,"2000 m", max = 22,is.7 =5 )  
  plot_Post_fire_6 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[1]][[6]],iterations,"2500 m", max = 22,is.7 =6 )  
  plot_Post_fire_7 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[1]][[7]],iterations,"3000 m", max = 22,is.7 =7 )  
  plot_Post_fire_8 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[1]][[8]],iterations,"Percolation distance (3842 m)", max = 22,is.7 =8 )  
  plot_Post_fire_9 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[1]][[9]],iterations,"5000 m", max = 24,is.7 =9 )  
  
  ####Els multiplots PER FUNCTIONAL TRAITS AFTER###
  library(gtable)    
  library(grid)
  library(gridExtra) 
  
  # Get the gtables
  gA <- ggplotGrob(plot_Post_fire_1)
  gB <- ggplotGrob(plot_Post_fire_2)
  gC <- ggplotGrob(plot_Post_fire_3)
  gD <- ggplotGrob(plot_Post_fire_4)
  gE <- ggplotGrob(plot_Post_fire_5)
  gF <- ggplotGrob(plot_Post_fire_6)
  gG <- ggplotGrob(plot_Post_fire_7)
  gH <- ggplotGrob(plot_Post_fire_8)
  gI <- ggplotGrob(plot_Post_fire_9)
  
  # Arrange the two charts
  # The legend boxes are centered
  grid.newpage()
  png(filename = paste("Figures/",file_name,".png", sep=""),width=14000,height=12000,units="px",res=800)
  grid.arrange(gA,gB,gC,gD,gE,gF,gG,gH,gI, ncol = 3, nrow=3)
  dev.off()
}
ARTICLE_FOR_GGPLOT <- function(Post_Val_Matrix, title, subtitle, max, is.7){
  # Packages needed
  require(ggplot2)
  require(devtools)
  require(metR) # IF not installed: install_github("eliocamp/metR")
  
  df <- as.data.frame((Post_Val_Matrix)) # Make the output a data.frame
  out <- list()                             # Generate and out in list format 
  out.bo <- list()                          # Generate and out in list format
  for(i in 1:length(df)){                   # Loop for each column of the data.frame (that corresponds to Wild. Size)
    out.temp <- df[1:20,i]    # Select all rows in each column and make them a data.frame
    out[[i]] <- out.temp      # Put all the rows in data.frame format into the out list 
  }
  for (e in 1:20) {             # For every 20 times 
    ee<- data.frame(out[[e]],rep(e,times=20),seq(from=0.05, to=1, by=0.05)) # create a new data.frame with the Richness/area/intensity
    colnames(ee) <- c("Richness","area","intensity") # Change column names
    out.bo[[e]] <- ee
  }
  final <- do.call(rbind,out.bo) # Joints the several lists objects created in one big data.frame
  if(is.7==7)#;min(final[,1])>22) # Condition that if the minimum richness is GREATER than 90 will follow the below path
    plot <-ggplot(final,aes(x=area, y=intensity))+
    # Raster with all the data, hjust and vjust hide the bakground and INTERPOLATE= T dissolve the rasters 
    geom_raster(aes(x = area,y=intensity,fill=Richness),hjust=1,vjust=1,interpolate = T)+
    # Scale y axis and set the levels
    scale_y_continuous(limits = c(0.05,1),expand = c(0,0.005), 
                       breaks = seq(from=0.05, to=1, by=0.05),
                       labels =c("0","10","","20","","30","","40","","50",
                                 "","60","","70","","80","","90","","100"))+
    # Scale x axis and set the levels--> WARNING: The maximum distance is 21km!! IS A LIE!
    scale_x_continuous(limits = c(1,20),expand = c(0,0.1),
                       breaks = seq(from=1, to=20, by=1),
                       labels =c("0","2","","4","","6","","8","","10",
                                 "","12","","14","","16","","18","","20"))+
    # Scalate thte colours, spectral means going from red to blue LIMITS force the range to be btween 0 & 100
    scale_fill_distiller(palette = "Spectral",direction =1, limits = c(0,max))+
    # Labs 
    labs(x="Size (km)",y="Intensity (%)",title = title,subtitle = subtitle)+
    geom_contour(aes(z=Richness), colour="black",size=0.5, lty="solid",stat="contour",
                 lineend="round",linejoin ="round",binwidth = 2)+
    geom_point(aes(x=10.5, y=0.3), col="red", size=9, shape=21 , fill="black",stroke=2)+
    # Theme stuff
    theme_classic()+
    theme(axis.line = element_line(size = 0.05, linetype = "solid"), 
          axis.ticks = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(0.4,1.5,0.4,1.5, "cm"),
          axis.title = element_text(size = 18,face = "bold"), 
          axis.text = element_text(size = 18,face = "bold"), 
          axis.text.x = element_text(size = 17, color = "black"), 
          axis.text.y = element_text(size = 17, color = "black"), 
          legend.position = "none",
          legend.text = element_blank(),
          legend.title = element_blank(),
          panel.background = element_rect(colour = "black", size = 2.5, linetype = "solid"),
          plot.background = element_rect(colour = "black", size = 1, linetype = "solid"), 
          legend.background = element_blank(),
          plot.subtitle = element_text(size = 18,face = "bold", colour = "gray0"),
          plot.title = element_text(size = 16, face = "bold"))
  else plot <-ggplot(final,aes(x=area, y=intensity))+ # Condition that if the minimum richness is SMALLER than 90 will follow 
    # Raster with all the data, hjust and vjust hide the bakground and INTERPOLATE= T dissolve the rasters 
    geom_raster(aes(x = area,y=intensity,fill=Richness),hjust=1,vjust=1,interpolate = T)+
    # Scale y axis and set the levels
    scale_y_continuous(limits = c(0.05,1),expand = c(0,0.005),
                       breaks = seq(from=0.05, to=1, by=0.05),
                       labels =c("0","10","","20","","30","","40","","50",
                                 "","60","","70","","80","","90","","100"))+
    # Scale x axis and set the levels--> WARNING: The maximum distance is 21km!! IS A LIE!
    scale_x_continuous(limits = c(1,20),expand = c(0,0.1),
                       breaks = seq(from=1, to=20, by=1),
                       labels =c("0","2","","4","","6","","8","","10",
                                 "","12","","14","","16","","18","","20"))+
    # Scalate thte colours, spectral means going from red to blue. LIMITS force the range to be btween 0 & 100
    scale_fill_distiller(palette = "Spectral",direction =1, limits = c(0,max))+
    # Labs 
    labs(x="Size (km)",y="Intensity (%)",title = title,subtitle = subtitle)+
    # Sets the contour by Richness. BINWIDTH: the value between each contour line 
    geom_contour(aes(z=Richness), colour="black",size=0.5, lty="solid",stat="contour",
                 lineend="round",linejoin ="round",binwidth = 2)+
    geom_point(aes(x=10.5, y=0.3), col="red", size=9, shape=21 , fill="black",stroke=2)+
    # Determines the labels with a white STROKE below them
    geom_text_contour(aes(z = Richness),stroke = 0.2, size= 10, rotate = F,check_overlap=T)+
    # Theme stuff
    theme_classic()+
    theme(axis.line = element_line(size = 0.05, linetype = "solid"), 
          axis.ticks = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(0.4,1.5,0.4,1.5, "cm"),
          axis.title = element_text(size = 18,face = "bold"), 
          axis.text = element_text(size = 18,face = "bold"), 
          axis.text.x = element_text(size = 17, color = "black"), 
          axis.text.y = element_text(size = 17, color = "black"), 
          legend.position = "none",
          legend.text = element_blank(),
          legend.title = element_blank(),
          panel.background = element_rect(colour = "black", size = 2.5, linetype = "solid"),
          plot.background = element_rect(colour = "black", size = 1, linetype = "solid"), 
          legend.background = element_blank(),
          plot.subtitle = element_text(size = 18,face = "bold", colour = "gray0"),
          plot.title = element_text(size = 18, face = "bold"))
}

# Plots for every dispersal ability group (9 plots together) - Figure 6B and Supplementary S4
Combined_article_plot(Burn_rec_ENTRADA_0.1, file_name ="Alpha 0.1", iterations = "T1")
Combined_article_plot(Burn_rec_ENTRADA_0.3, file_name ="Alpha 0.3", iterations = "T1")
Combined_article_plot(Burn_rec_ENTRADA_0.5, file_name ="Alpha 0.5", iterations = "T1")
Combined_article_plot(Burn_rec_ENTRADA_0.7, file_name ="Alpha 0.7", iterations = "T1")
Combined_article_plot(Burn_rec_ENTRADA_0.9, file_name ="Alpha 0.9", iterations = "T1")

# Plots the whole metacommunity (all dispersal abilities together) - Figure 6A and Supplementary S4
GLOBAL_ARTICLE_FOR_GGPLOT <- function(Post_Val_Matrix, title, subtitle){
  # Packages needed
  require(ggplot2)
  require(devtools)
  require(metR) # IF not installed: install_github("eliocamp/metR")
  
  df <- as.data.frame((Post_Val_Matrix[[1]][[1]]+Post_Val_Matrix[[1]][[2]]+
                         Post_Val_Matrix[[1]][[3]]+Post_Val_Matrix[[1]][[4]]+
                         Post_Val_Matrix[[1]][[5]]+Post_Val_Matrix[[1]][[6]]+
                         Post_Val_Matrix[[1]][[7]]+Post_Val_Matrix[[1]][[8]]+
                         Post_Val_Matrix[[1]][[9]])) # Make the output a data.frame
  out <- list()                             # Generate and out in list format 
  out.bo <- list()                          # Generate and out in list format
  for(i in 1:length(df)){                   # Loop for each column of the data.frame (that corresponds to Wild. Size)
    out.temp <- df[1:20,i]    # Select all rows in each column and make them a data.frame
    out[[i]] <- out.temp      # Put all the rows in data.frame format into the out list 
  }
  for (e in 1:20) {             # For every 20 times 
    ee<- data.frame(out[[e]],rep(e,times=20),seq(from=0.05, to=1, by=0.05)) # create a new data.frame with the Richness/area/intensity
    colnames(ee) <- c("Richness","area","intensity") # Change column names
    out.bo[[e]] <- ee
  }
  final <- do.call(rbind,out.bo) # Joints the several lists objects created in one big data.frame  
  if(min(final[,1])>90) # Condition that if the minimum richness is GREATER than 90 will follow the below path
    plot <-ggplot(final,aes(x=area, y=intensity))+
    # Raster with all the data, hjust and vjust hide the bakground and INTERPOLATE= T dissolve the rasters 
    geom_raster(aes(x = area,y=intensity,fill=Richness),hjust=1,vjust=1,interpolate = T)+
    # Scale y axis and set the levels
    scale_y_continuous(limits = c(0.05,1),expand = c(0,0.005), 
                       breaks = seq(from=0.05, to=1, by=0.05),
                       labels =c("0","10","","20","","30","","40","","50",
                                 "","60","","70","","80","","90","","100"))+
    # Scale x axis and set the levels--> WARNING: The maximum distance is 21km!! IS A LIE!
    scale_x_continuous(limits = c(1,20),expand = c(0,0.1),
                       breaks = seq(from=1, to=20, by=1),
                       labels =c("0","2","","4","","6","","8","","10",
                                 "","12","","14","","16","","18","","20"))+
    # Scalate thte colours, spectral means going from red to blue LIMITS force the range to be btween 0 & 100
    scale_fill_distiller(palette = "Spectral",direction =1, limits = c(0,150))+
    # Labs 
    labs(x="Size (km)",y="Intensity (%)",title = title,subtitle = subtitle)+
    geom_point(aes(x=10.5, y=0.3), col="red", size=9, shape=21 , fill="black",stroke=2)+
    # Theme stuff
    theme_classic()+
    theme(axis.line = element_line(size = 0.05, linetype = "solid"), 
          axis.ticks = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(0.4,1.5,0.4,1.5, "cm"),
          axis.title = element_text(size = 18,face = "bold"), 
          axis.text = element_text(size = 18,face = "bold"), 
          axis.text.x = element_text(size = 30, color = "black"), 
          axis.text.y = element_text(size = 30, color = "black"), 
          legend.position = "none",
          legend.text = element_blank(),
          legend.title = element_blank(),
          panel.background = element_rect(colour = "black", size = 2.5, linetype = "solid"),
          plot.background = element_rect(colour = "black", size = 1, linetype = "solid"), 
          legend.background = element_blank(),
          plot.subtitle = element_text(size = 16,face = "bold", colour = "gray0"),
          plot.title = element_text(size = 24, face = "bold"))
  else plot <-ggplot(final,aes(x=area, y=intensity))+ # Condition that if the minimum richness is SMALLER than 90 will follow 
    # Raster with all the data, hjust and vjust hide the bakground and INTERPOLATE= T dissolve the rasters 
    geom_raster(aes(x = area,y=intensity,fill=Richness),hjust=1,vjust=1,interpolate = T)+
    # Scale y axis and set the levels
    scale_y_continuous(limits = c(0.05,1),expand = c(0,0.005),
                       breaks = seq(from=0.05, to=1, by=0.05),
                       labels =c("0","10","","20","","30","","40","","50",
                                 "","60","","70","","80","","90","","100"))+
    # Scale x axis and set the levels--> WARNING: The maximum distance is 21km!! IS A LIE!
    scale_x_continuous(limits = c(1,20),expand = c(0,0.1),
                       breaks = seq(from=1, to=20, by=1),
                       labels =c("0","2","","4","","6","","8","","10",
                                 "","12","","14","","16","","18","","20"))+
    # Scalate thte colours, spectral means going from red to blue. LIMITS force the range to be btween 0 & 100
    scale_fill_distiller(palette = "Spectral",direction =1, limits = c(0,150))+
    # Labs 
    labs(x="Size (km)",y="Intensity (%)",title = title,subtitle = subtitle)+
    # Sets the contour by Richness. BINWIDTH: the value between each contour line 
    geom_contour(aes(z=Richness), colour="black",size=0.5,lty="solid",stat="contour",
                 lineend="round",linejoin ="round",binwidth = 10)+
    # Determines the labels with a white STROKE below them
    geom_point(aes(x=10.5, y=0.3), col="red", size=9, shape=21 , fill="black",stroke=2)+
    geom_text_contour(aes(z = Richness),stroke = 0.2, size=20, rotate = F,check_overlap=T)+
    # Theme stuff
    theme_classic()+
    theme(axis.line = element_line(size = 0.05, linetype = "solid"), 
          axis.ticks = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(0.4,1.5,0.4,1.5, "cm"),
          axis.title = element_text(size = 18,face = "bold"), 
          axis.text = element_text(size = 18,face = "bold"), 
          axis.text.x = element_text(size = 30, color = "black"), 
          axis.text.y = element_text(size = 30, color = "black"), 
          legend.position = "none",
          legend.text = element_blank(),
          legend.title = element_blank(),
          panel.background = element_rect(colour = "black", size = 2.5, linetype = "solid"),
          plot.background = element_rect(colour = "black", size = 1, linetype = "solid"), 
          legend.background = element_blank(),
          plot.subtitle = element_text(size = 16,face = "bold", colour = "gray0"),
          plot.title = element_text(size = 24, face = "bold"))
  
  plot
}

png(filename = "Figures/Globalplot_01.png",width=14000,height=12000,units="px",res=1000)
GLOBAL_ARTICLE_FOR_GGPLOT(Burn_rec_ENTRADA_0.1,"Entire metacommunity T1","Shared = 0.1")
dev.off()
png(filename = "Figures/Globalplot_03.png",width=14000,height=12000,units="px",res=1000)
GLOBAL_ARTICLE_FOR_GGPLOT(Burn_rec_ENTRADA_0.3,"Entire metacommunity T1","Shared = 0.3")
dev.off()
png(filename = "Figures/Globalplot_05.png",width=14000,height=12000,units="px",res=1000)
GLOBAL_ARTICLE_FOR_GGPLOT(Burn_rec_ENTRADA_0.5,"Entire metacommunity T1","Shared = 0.5")
dev.off()
png(filename = "Figures/Globalplot_07.png",width=14000,height=12000,units="px",res=1000)
GLOBAL_ARTICLE_FOR_GGPLOT(Burn_rec_ENTRADA_0.7,"Entire metacommunity T1","Shared = 0.7")
dev.off()
png(filename = "Figures/Globalplot_09.png",width=14000,height=12000,units="px",res=1000)
GLOBAL_ARTICLE_FOR_GGPLOT(Burn_rec_ENTRADA_0.9,"Entire metacommunity T1","Shared = 0.9")
dev.off()

## Mean richness SD plot (Supplementary S6)
Combined_article_plot_sd <- function(Post_Val_Matrix, file_name, iterations ){
  plot_Post_fire_1 <- ARTICLE_FOR_GGPLOT_sd(Post_Val_Matrix[[3]][[1]],iterations,"250 m")  
  plot_Post_fire_2 <- ARTICLE_FOR_GGPLOT_sd(Post_Val_Matrix[[3]][[2]],iterations,"500 m")  
  plot_Post_fire_3 <- ARTICLE_FOR_GGPLOT_sd(Post_Val_Matrix[[3]][[3]],iterations,"1000 m")  
  plot_Post_fire_4 <- ARTICLE_FOR_GGPLOT_sd(Post_Val_Matrix[[3]][[4]],iterations,"1500 m")  
  plot_Post_fire_5 <- ARTICLE_FOR_GGPLOT_sd(Post_Val_Matrix[[3]][[5]],iterations,"2000 m")  
  plot_Post_fire_6 <- ARTICLE_FOR_GGPLOT_sd(Post_Val_Matrix[[3]][[6]],iterations,"2500 m")  
  plot_Post_fire_7 <- ARTICLE_FOR_GGPLOT_sd(Post_Val_Matrix[[3]][[7]],iterations,"3000 m")  
  plot_Post_fire_8 <- ARTICLE_FOR_GGPLOT_sd(Post_Val_Matrix[[3]][[8]],iterations,"Percolation distance (3842 m)")  
  plot_Post_fire_9 <- ARTICLE_FOR_GGPLOT_sd(Post_Val_Matrix[[3]][[9]],iterations,"5000 m")  
  
  ####Els multiplots PER FUNCTIONAL TRAITS AFTER###
  library(gtable)    
  library(grid)
  library(gridExtra) 
  
  # Get the gtables
  gA <- ggplotGrob(plot_Post_fire_1)
  gB <- ggplotGrob(plot_Post_fire_2)
  gC <- ggplotGrob(plot_Post_fire_3)
  gD <- ggplotGrob(plot_Post_fire_4)
  gE <- ggplotGrob(plot_Post_fire_5)
  gF <- ggplotGrob(plot_Post_fire_6)
  gG <- ggplotGrob(plot_Post_fire_7)
  gH <- ggplotGrob(plot_Post_fire_8)
  gI <- ggplotGrob(plot_Post_fire_9)
  
  # Arrange the two charts
  # The legend boxes are centered
  grid.newpage()
  png(filename = paste("Figures/",file_name,".png", sep=""),width=14000,height=12000,units="px",res=800)
  grid.arrange(gA,gB,gC,gD,gE,gF,gG,gH,gI, ncol = 3, nrow=3)
  dev.off()
}
ARTICLE_FOR_GGPLOT_sd <- function(Post_Val_Matrix, title, subtitle){
  # Packages needed
  require(ggplot2)
  require(devtools)
  require(metR) # IF not installed: install_github("eliocamp/metR")
  
  df <- as.data.frame((Post_Val_Matrix)) # Make the output a data.frame
  out <- list()                             # Generate and out in list format 
  out.bo <- list()                          # Generate and out in list format
  for(i in 1:length(df)){                   # Loop for each column of the data.frame (that corresponds to Wild. Size)
    out.temp <- df[1:20,i]    # Select all rows in each column and make them a data.frame
    out[[i]] <- out.temp      # Put all the rows in data.frame format into the out list 
  }
  for (e in 1:20) {             # For every 20 times 
    ee<- data.frame(out[[e]],rep(e,times=20),seq(from=0.05, to=1, by=0.05)) # create a new data.frame with the Richness/area/intensity
    colnames(ee) <- c("Richness","area","intensity") # Change column names
    out.bo[[e]] <- ee
  }
  final <- do.call(rbind,out.bo) # Joints the several lists objects created in one big data.frame
  final[which(is.na(final)),1] <- 0
  
#return(final)
  if(min(final[,1])>90) # Condition that if the minimum richness is GREATER than 90 will follow the below path
    plot <-ggplot(final,aes(x=area, y=intensity))+
    # Raster with all the data, hjust and vjust hide the bakground and INTERPOLATE= T dissolve the rasters 
    geom_raster(aes(x = area,y=intensity,fill=Richness),hjust=1,vjust=1,interpolate = T)+
    # Scale y axis and set the levels
    scale_y_continuous(limits = c(0.05,1),expand = c(0,0.005), 
                       breaks = seq(from=0.05, to=1, by=0.05),
                       labels =c("0","10","","20","","30","","40","","50",
                                 "","60","","70","","80","","90","","100"))+
    # Scale x axis and set the levels--> WARNING: The maximum distance is 21km!! IS A LIE!
    scale_x_continuous(limits = c(1,20),expand = c(0,0.1),
                       breaks = seq(from=1, to=20, by=1),
                       labels =c("0","2","","4","","6","","8","","10",
                                 "","12","","14","","16","","18","","20"))+
    # Scalate thte colours, spectral means going from red to blue LIMITS force the range to be btween 0 & 100
    scale_fill_distiller(palette = "Spectral",direction =1, limits = c(0,100))+
    # Labs 
    labs(x="Size (km)",y="Intensity (%)",title = title,subtitle = subtitle)+
    # Theme stuff
    theme_classic()+
    theme(axis.line = element_line(size = 0.05, linetype = "solid"), 
          axis.ticks = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(0.4,1.5,0.4,1.5, "cm"),
          axis.title = element_text(size = 18,face = "bold"), 
          axis.text = element_text(size = 18,face = "bold"), 
          axis.text.x = element_text(size = 17, color = "black"), 
          axis.text.y = element_text(size = 17, color = "black"), 
          legend.position = "none",
          legend.text = element_blank(),
          legend.title = element_blank(),
          panel.background = element_rect(colour = "black", size = 2.5, linetype = "solid"),
          plot.background = element_rect(colour = "black", size = 1, linetype = "solid"), 
          legend.background = element_blank(),
          plot.subtitle = element_text(size = 18,face = "bold", colour = "gray0"),
          plot.title = element_text(size = 16, face = "bold"))
  else plot <-ggplot(final,aes(x=area, y=intensity))+ # Condition that if the minimum richness is SMALLER than 90 will follow 
    # Raster with all the data, hjust and vjust hide the bakground and INTERPOLATE= T dissolve the rasters 
    geom_raster(aes(x = area,y=intensity,fill=Richness),hjust=1,vjust=1,interpolate = T)+
    # Scale y axis and set the levels
    scale_y_continuous(limits = c(0.05,1),expand = c(0,0.005),
                       breaks = seq(from=0.05, to=1, by=0.05),
                       labels =c("0","10","","20","","30","","40","","50",
                                 "","60","","70","","80","","90","","100"))+
    # Scale x axis and set the levels--> WARNING: The maximum distance is 21km!! IS A LIE!
    scale_x_continuous(limits = c(1,20),expand = c(0,0.1),
                       breaks = seq(from=1, to=20, by=1),
                       labels =c("0","2","","4","","6","","8","","10",
                                 "","12","","14","","16","","18","","20"))+
    # Scalate thte colours, spectral means going from red to blue. LIMITS force the range to be btween 0 & 100
    scale_fill_distiller(palette = "BrBG",direction =-1, limits = c(0,20))+
    # Labs 
    labs(x="Size (km)",y="Intensity (%)",title = title, subtitle = subtitle)+
    # Sets the contour by Richness. BINWIDTH: the value between each contour line 
    geom_contour(aes(z=Richness), colour="black",size=0.5, lty="solid",stat="contour",
                 lineend="round",linejoin ="round",binwidth = 1)+
    # Determines the labels with a white STROKE below them
    geom_text_contour(aes(z = Richness),stroke = 0.2, size= 15, rotate = F,check_overlap=T)+
    # Theme stuff
    theme_classic()+
    theme(axis.line = element_line(size = 0.05, linetype = "solid"), 
          axis.ticks = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(0.4,1.5,0.4,1.5, "cm"),
          axis.title = element_text(size = 18,face = "bold"), 
          axis.text = element_text(size = 18,face = "bold"), 
          axis.text.x = element_text(size = 17, color = "black"), 
          axis.text.y = element_text(size = 17, color = "black"), 
          legend.position = "none",
          legend.text = element_blank(),
          legend.title = element_blank(),
          panel.background = element_rect(colour = "black", size = 2.5, linetype = "solid"),
          plot.background = element_rect(colour = "black", size = 1, linetype = "solid"), 
          legend.background = element_blank(),
          plot.subtitle = element_text(size = 18,face = "bold", colour = "gray0"),
          plot.title = element_text(size = 18, face = "bold"))
}

Combined_article_plot_sd(Burn_rec_ENTRADA_0.1, file_name ="SD_Alpha 0.1", iterations = "T1")
Combined_article_plot_sd(Burn_rec_ENTRADA_0.3, file_name ="SD_Alpha 0.3", iterations = "T1")
Combined_article_plot_sd(Burn_rec_ENTRADA_0.5, file_name ="SD_Alpha 0.5", iterations = "T1")
Combined_article_plot_sd(Burn_rec_ENTRADA_0.7, file_name ="SD_Alpha 0.7", iterations = "T1")
Combined_article_plot_sd(Burn_rec_ENTRADA_0.9, file_name ="SD_Alpha 0.9", iterations = "T1")

# T100_____________________________________________________________________________________________________

#Load data 
load("C:/Users/Cunilleramontcusi/Desktop/Burn/T100/0.1/Burn_rec_ENTRADA_0.1_100_JOINT.RData")
load("C:/Users/Cunilleramontcusi/Desktop/Burn/T100/0.3/Burn_rec_ENTRADA_0.3_100_JOINT.RData")
load("C:/Users/Cunilleramontcusi/Desktop/Burn/T100/0.5/Burn_rec_ENTRADA_0.5_100_JOINT.RData")
load("C:/Users/Cunilleramontcusi/Desktop/Burn/T100/0.7/Burn_rec_ENTRADA_0.7_100_JOINT.RData")
load("C:/Users/Cunilleramontcusi/Desktop/Burn/T100/0.9/Burn_rec_ENTRADA_0.9_100_JOINT.RData")

# Data treatment to joint several parts of the data 
buidatge_BURN_output <- function(Burn_05,Burn_10,Burn_15,Burn_20,Burn_25,
                                 Burn_30,Burn_35,Burn_40,Burn_45,Burn_50,
                                 Burn_55,Burn_60,Burn_65,Burn_70,Burn_75,
                                 Burn_80,Burn_85,Burn_90,Burn_95,Burn_100){
  Prop.spp.disponibles <- list("1","2","3","4","5","6","7","8","9")
  for (a in 1:9) {
    Prop.spp.disponibles[[a]] <- rbind(Burn_05[[1]][[a]],Burn_10[[1]][[a]],Burn_15[[1]][[a]],Burn_20[[1]][[a]],Burn_25[[1]][[a]],
                                       Burn_30[[1]][[a]],Burn_35[[1]][[a]],Burn_40[[1]][[a]],Burn_45[[1]][[a]],Burn_50[[1]][[a]],
                                       Burn_55[[1]][[a]],Burn_60[[1]][[a]],Burn_65[[1]][[a]],Burn_70[[1]][[a]],Burn_75[[1]][[a]],
                                       Burn_80[[1]][[a]],Burn_85[[1]][[a]],Burn_90[[1]][[a]],Burn_95[[1]][[a]],Burn_100[[1]][[a]])  
  }
  Iteraciones_para_fill<- list("1","2","3","4","5","6","7","8","9")
  for (b in 1:9) {
    Iteraciones_para_fill[[b]] <- rbind(Burn_05[[2]][[b]],Burn_10[[2]][[b]],Burn_15[[2]][[b]],Burn_20[[2]][[b]],Burn_25[[2]][[b]],
                                        Burn_30[[2]][[b]],Burn_35[[2]][[b]],Burn_40[[2]][[b]],Burn_45[[2]][[b]],Burn_50[[2]][[b]],
                                        Burn_55[[2]][[b]],Burn_60[[2]][[b]],Burn_65[[2]][[b]],Burn_70[[2]][[b]],Burn_75[[2]][[b]],
                                        Burn_80[[2]][[b]],Burn_85[[2]][[b]],Burn_90[[2]][[b]],Burn_95[[2]][[b]],Burn_100[[2]][[b]])  
  }
  sd.spp.disponibles<- list("1","2","3","4","5","6","7","8","9")
  for (c in 1:9) {
    sd.spp.disponibles[[c]] <- rbind(Burn_05[[3]][[c]],Burn_10[[3]][[c]],Burn_15[[3]][[c]],Burn_20[[3]][[c]],Burn_25[[3]][[c]],
                                     Burn_30[[3]][[c]],Burn_35[[3]][[c]],Burn_40[[3]][[c]],Burn_45[[3]][[c]],Burn_50[[3]][[c]],
                                     Burn_55[[3]][[c]],Burn_60[[3]][[c]],Burn_65[[3]][[c]],Burn_70[[3]][[c]],Burn_75[[3]][[c]],
                                     Burn_80[[3]][[c]],Burn_85[[3]][[c]],Burn_90[[3]][[c]],Burn_95[[3]][[c]],Burn_100[[3]][[c]])  
  }
  sd.iteraciones.fill<- list("1","2","3","4","5","6","7","8","9")
  for (d in 1:9) {
    sd.iteraciones.fill[[d]] <- rbind(Burn_05[[4]][[d]],Burn_10[[4]][[d]],Burn_15[[4]][[d]],Burn_20[[4]][[d]],Burn_25[[4]][[d]],
                                      Burn_30[[4]][[d]],Burn_35[[4]][[d]],Burn_40[[4]][[d]],Burn_45[[4]][[d]],Burn_50[[4]][[d]],
                                      Burn_55[[4]][[d]],Burn_60[[4]][[d]],Burn_65[[4]][[d]],Burn_70[[4]][[d]],Burn_75[[4]][[d]],
                                      Burn_80[[4]][[d]],Burn_85[[4]][[d]],Burn_90[[4]][[d]],Burn_95[[4]][[d]],Burn_100[[4]][[d]])  
  }
  out <- list(Prop.spp.disponibles,Iteraciones_para_fill,sd.spp.disponibles,sd.iteraciones.fill)
  out 
}

# generate an output 0.1 T100
Burn_0.1_T100 <- buidatge_BURN_output(Burn_rec_ENTRADA_0.1_100_05,Burn_rec_ENTRADA_0.1_100_10,Burn_rec_ENTRADA_0.1_100_15,
                                      Burn_rec_ENTRADA_0.1_100_20,Burn_rec_ENTRADA_0.1_100_25,Burn_rec_ENTRADA_0.1_100_30,
                                      Burn_rec_ENTRADA_0.1_100_35,Burn_rec_ENTRADA_0.1_100_40,Burn_rec_ENTRADA_0.1_100_45,
                                      Burn_rec_ENTRADA_0.1_100_50,Burn_rec_ENTRADA_0.1_100_55,Burn_rec_ENTRADA_0.1_100_60,
                                      Burn_rec_ENTRADA_0.1_100_65,Burn_rec_ENTRADA_0.1_100_70,Burn_rec_ENTRADA_0.1_100_75,
                                      Burn_rec_ENTRADA_0.1_100_80,Burn_rec_ENTRADA_0.1_100_85,Burn_rec_ENTRADA_0.1_100_90,
                                      Burn_rec_ENTRADA_0.1_100_95,Burn_rec_ENTRADA_0.1_100_1)

# generate an output 0.3 T100
Burn_0.3_T100 <- buidatge_BURN_output(Burn_rec_ENTRADA_0.3_100_05,Burn_rec_ENTRADA_0.3_100_10,Burn_rec_ENTRADA_0.3_100_15,
                                      Burn_rec_ENTRADA_0.3_100_20,Burn_rec_ENTRADA_0.3_100_25,Burn_rec_ENTRADA_0.3_100_30,
                                      Burn_rec_ENTRADA_0.3_100_35,Burn_rec_ENTRADA_0.3_100_40,Burn_rec_ENTRADA_0.3_100_45,
                                      Burn_rec_ENTRADA_0.3_100_50,Burn_rec_ENTRADA_0.3_100_55,Burn_rec_ENTRADA_0.3_100_60,
                                      Burn_rec_ENTRADA_0.3_100_65,Burn_rec_ENTRADA_0.3_100_70,Burn_rec_ENTRADA_0.3_100_75,
                                      Burn_rec_ENTRADA_0.3_100_80,Burn_rec_ENTRADA_0.3_100_85,Burn_rec_ENTRADA_0.3_100_90,
                                      Burn_rec_ENTRADA_0.3_100_95,Burn_rec_ENTRADA_0.3_100_1)

# generate an output 0.5 T100
Burn_0.5_T100 <- buidatge_BURN_output(Burn_rec_ENTRADA_0.5_100_05,Burn_rec_ENTRADA_0.5_100_10,Burn_rec_ENTRADA_0.5_100_15,
                                      Burn_rec_ENTRADA_0.5_100_20,Burn_rec_ENTRADA_0.5_100_25,Burn_rec_ENTRADA_0.5_100_30,
                                      Burn_rec_ENTRADA_0.5_100_35,Burn_rec_ENTRADA_0.5_100_40,Burn_rec_ENTRADA_0.5_100_45,
                                      Burn_rec_ENTRADA_0.5_100_50,Burn_rec_ENTRADA_0.5_100_55,Burn_rec_ENTRADA_0.5_100_60,
                                      Burn_rec_ENTRADA_0.5_100_65,Burn_rec_ENTRADA_0.5_100_70,Burn_rec_ENTRADA_0.5_100_75,
                                      Burn_rec_ENTRADA_0.5_100_80,Burn_rec_ENTRADA_0.5_100_85,Burn_rec_ENTRADA_0.5_100_90,
                                      Burn_rec_ENTRADA_0.5_100_95,Burn_rec_ENTRADA_0.5_100_1)

# generate an output 0.7 T100
Burn_0.7_T100 <- buidatge_BURN_output(Burn_rec_ENTRADA_0.7_100_05,Burn_rec_ENTRADA_0.7_100_10,Burn_rec_ENTRADA_0.7_100_15,
                                      Burn_rec_ENTRADA_0.7_100_20,Burn_rec_ENTRADA_0.7_100_25,Burn_rec_ENTRADA_0.7_100_30,
                                      Burn_rec_ENTRADA_0.7_100_35,Burn_rec_ENTRADA_0.7_100_40,Burn_rec_ENTRADA_0.7_100_45,
                                      Burn_rec_ENTRADA_0.7_100_50,Burn_rec_ENTRADA_0.7_100_55,Burn_rec_ENTRADA_0.7_100_60,
                                      Burn_rec_ENTRADA_0.7_100_65,Burn_rec_ENTRADA_0.7_100_70,Burn_rec_ENTRADA_0.7_100_75,
                                      Burn_rec_ENTRADA_0.7_100_80,Burn_rec_ENTRADA_0.7_100_85,Burn_rec_ENTRADA_0.7_100_90,
                                      Burn_rec_ENTRADA_0.7_100_95,Burn_rec_ENTRADA_0.7_100_1)
# generate an output 0.9 T100
Burn_0.9_T100 <- buidatge_BURN_output(Burn_rec_ENTRADA_0.9_100_05,Burn_rec_ENTRADA_0.9_100_10,Burn_rec_ENTRADA_0.9_100_15,
                                      Burn_rec_ENTRADA_0.9_100_20,Burn_rec_ENTRADA_0.9_100_25,Burn_rec_ENTRADA_0.9_100_30,
                                      Burn_rec_ENTRADA_0.9_100_35,Burn_rec_ENTRADA_0.9_100_40,Burn_rec_ENTRADA_0.9_100_45,
                                      Burn_rec_ENTRADA_0.9_100_50,Burn_rec_ENTRADA_0.9_100_55,Burn_rec_ENTRADA_0.9_100_60,
                                      Burn_rec_ENTRADA_0.9_100_65,Burn_rec_ENTRADA_0.9_100_70,Burn_rec_ENTRADA_0.9_100_75,
                                      Burn_rec_ENTRADA_0.9_100_80,Burn_rec_ENTRADA_0.9_100_85,Burn_rec_ENTRADA_0.9_100_90,
                                      Burn_rec_ENTRADA_0.9_100_95,Burn_rec_ENTRADA_0.9_100_1)
#____________________________________________________________________________________

# RR (one iteration after wildfire)_________________________
#############################################################################################################################
#############################################################################################################################

Combined_article_plot <- function(Post_Val_Matrix, file_name, iterations ){
  plot_Post_fire_1 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[2]][[1]],iterations,"250 m")  
  plot_Post_fire_2 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[2]][[2]],iterations,"500 m")  
  plot_Post_fire_3 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[2]][[3]],iterations,"1000 m")  
  plot_Post_fire_4 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[2]][[4]],iterations,"1500 m")  
  plot_Post_fire_5 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[2]][[5]],iterations,"2000 m")  
  plot_Post_fire_6 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[2]][[6]],iterations,"2500 m")  
  plot_Post_fire_7 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[2]][[7]],iterations,"3000 m")  
  plot_Post_fire_8 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[2]][[8]],iterations,"Percolation distance (3842 m)")  
  plot_Post_fire_9 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[2]][[9]],iterations,"5000 m")  
  
  ####Els multiplots PER FUNCTIONAL TRAITS AFTER###
  library(gtable)    
  library(grid)
  library(gridExtra) 
  
  # Get the gtables
  gA <- ggplotGrob(plot_Post_fire_1)
  gB <- ggplotGrob(plot_Post_fire_2)
  gC <- ggplotGrob(plot_Post_fire_3)
  gD <- ggplotGrob(plot_Post_fire_4)
  gE <- ggplotGrob(plot_Post_fire_5)
  gF <- ggplotGrob(plot_Post_fire_6)
  gG <- ggplotGrob(plot_Post_fire_7)
  gH <- ggplotGrob(plot_Post_fire_8)
  gI <- ggplotGrob(plot_Post_fire_9)
  
  # Arrange the two charts
  # The legend boxes are centered
  grid.newpage()
  png(filename = paste("C:/Users/Cunilleramontcusi/Desktop/",file_name,".png", sep=""),width=14000,height=12000,units="px",res=800)
  grid.arrange(gA,gB,gC,gD,gE,gF,gG,gH,gI, ncol = 3, nrow=3)
  dev.off()
}
ARTICLE_FOR_GGPLOT <- function(Post_Val_Matrix, title, subtitle, max, is.7){
  # Packages needed
  require(ggplot2)
  require(devtools)
  require(metR) # IF not installed: install_github("eliocamp/metR")
  
  df <- as.data.frame((Post_Val_Matrix)) # Make the output a data.frame
  out <- list()                             # Generate and out in list format 
  out.bo <- list()                          # Generate and out in list format
  for(i in 1:length(df)){                   # Loop for each column of the data.frame (that corresponds to Wild. Size)
    out.temp <- df[1:20,i]    # Select all rows in each column and make them a data.frame
    out[[i]] <- out.temp      # Put all the rows in data.frame format into the out list 
  }
  for (e in 1:20) {             # For every 20 times 
    ee<- data.frame(out[[e]],rep(e,times=20),seq(from=0.05, to=1, by=0.05)) # create a new data.frame with the Richness/area/intensity
    colnames(ee) <- c("Richness","area","intensity") # Change column names
    out.bo[[e]] <- ee
  }
  final <- do.call(rbind,out.bo) # Joints the several lists objects created in one big data.frame
  final[,1] <-(100/final[,1])/100
  
  ggplot(final,aes(x=area, y=intensity))+ # Condition that if the minimum richness is SMALLER than 90 will follow 
    # Raster with all the data, hjust and vjust hide the bakground and INTERPOLATE= T dissolve the rasters 
    geom_raster(aes(x = area,y=intensity,fill=Richness),hjust=1,vjust=1,interpolate = T)+
    # Scale y axis and set the levels
    scale_y_continuous(limits = c(0.05,1),expand = c(0,0.005),
                       breaks = seq(from=0.05, to=1, by=0.05),
                       labels =c("0","10","","20","","30","","40","","50",
                                 "","60","","70","","80","","90","","100"))+
    # Scale x axis and set the levels--> WARNING: The maximum distance is 21km!! IS A LIE!
    scale_x_continuous(limits = c(1,20),expand = c(0,0.1),
                       breaks = seq(from=1, to=20, by=1),
                       labels =c("0","2","","4","","6","","8","","10",
                                 "","12","","14","","16","","18","","20"))+
    # Scalate thte colours, spectral means going from red to blue. LIMITS force the range to be btween 0 & 100
    scale_fill_distiller(palette = "Spectral",direction =1, limits = c(0,1))+
    # Labs 
    labs(x="Size (km)",y="Intensity (%)",title = title,subtitle = subtitle)+
    # Sets the contour by Richness. BINWIDTH: the value between each contour line 
    geom_contour(aes(z=Richness), colour="black",size=0.5, lty="solid",stat="contour",
                 lineend="round",linejoin ="round",binwidth = 0.4)+
    geom_point(aes(x=10.5, y=0.3), col="red", size=9, shape=21 , fill="black",stroke=2)+
    # Determines the labels with a white STROKE below them
    geom_text_contour(aes(z = Richness),stroke = 0.2, size= 10, rotate = F,check_overlap=T)+
    # Theme stuff
    theme_classic()+
    theme(axis.line = element_line(size = 0.05, linetype = "solid"), 
          axis.ticks = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(0.4,1.5,0.4,1.5, "cm"),
          axis.title = element_text(size = 18,face = "bold"), 
          axis.text = element_text(size = 18,face = "bold"), 
          axis.text.x = element_text(size = 17, color = "black"), 
          axis.text.y = element_text(size = 17, color = "black"), 
          legend.position = "none",
          legend.text = element_blank(),
          legend.title = element_blank(),
          panel.background = element_rect(colour = "black", size = 2.5, linetype = "solid"),
          plot.background = element_rect(colour = "black", size = 1, linetype = "solid"), 
          legend.background = element_blank(),
          plot.subtitle = element_text(size = 18,face = "bold", colour = "gray0"),
          plot.title = element_text(size = 18, face = "bold"))
}

# Plots for every dispersal ability group (9 plots together) - Figure 6D and Supplementary S4
Combined_article_plot(Burn_0.1_T100, file_name ="Alpha 0.1_RR", iterations = "RR")
Combined_article_plot(Burn_0.3_T100, file_name ="Alpha 0.3_RR", iterations = "RR")
Combined_article_plot(Burn_0.5_T100, file_name ="Alpha 0.5_RR", iterations = "RR")
Combined_article_plot(Burn_0.7_T100, file_name ="Alpha 0.7_RR", iterations = "RR")
Combined_article_plot(Burn_0.9_T100, file_name ="Alpha 0.9_RR", iterations = "RR")

# DR __________________________________________________________________________________________________________________
#############################################################################################################################
#############################################################################################################################
# Richness change percentage (at 100 iterations). 
ARTICLE_FOR_GGPLOT_Rich_Change_Rate <- function(Post_Val_Matrix_T100, PreFire_Val_Matri, title, subtitle, max){
  # Packages needed
  require(ggplot2)
  require(devtools)
  require(metR) # IF not installed: install_github("eliocamp/metR")
  
  df <- as.data.frame((Post_Val_Matrix_T100)) # Make the output a data.frame
  out <- list()                             # Generate and out in list format 
  out.bo <- list()                          # Generate and out in list format
  for(i in 1:length(df)){                   # Loop for each column of the data.frame (that corresponds to Wild. Size)
    out.temp <- df[1:20,i]    # Select all rows in each column and make them a data.frame
    out[[i]] <- out.temp      # Put all the rows in data.frame format into the out list 
  }
  for (e in 1:20) {             # For every 20 times 
    ee<- data.frame(out[[e]],rep(e,times=20),seq(from=0.05, to=1, by=0.05)) # create a new data.frame with the Richness/area/intensity
    colnames(ee) <- c("Richness","area","intensity") # Change column names
    out.bo[[e]] <- ee
  }
  final <- do.call(rbind,out.bo)
  
  df <- as.data.frame((PreFire_Val_Matri)) # Make the output a data.frame
  out <- list()                             # Generate and out in list format 
  out.bo <- list()                          # Generate and out in list format
  for(i in 1:length(df)){                   # Loop for each column of the data.frame (that corresponds to Wild. Size)
    out.temp <- df[1:20,i]    # Select all rows in each column and make them a data.frame
    out[[i]] <- out.temp      # Put all the rows in data.frame format into the out list 
  }
  for (e in 1:20) {             # For every 20 times 
    ee<- data.frame(out[[e]],rep(e,times=20),seq(from=0.05, to=1, by=0.05)) # create a new data.frame with the Richness/area/intensity
    colnames(ee) <- c("Richness","area","intensity") # Change column names
    out.bo[[e]] <- ee
  }
  PreFire_final <- do.call(rbind,out.bo)
  
  final[,1] <- (PreFire_final[,1]*100)/final
  final[,1] <- ifelse(final[,1]>100,100,final[,1])
  
  ggplot(final,aes(x=area, y=intensity))+ # Condition that if the minimum richness is SMALLER than 90 will follow 
    # Raster with all the data, hjust and vjust hide the bakground and INTERPOLATE= T dissolve the rasters 
    geom_raster(aes(x = area,y=intensity,fill=Richness),hjust=1,vjust=1,interpolate = T)+
    # Scale y axis and set the levels
    scale_y_continuous(limits = c(0.05,1),expand = c(0,0.005),
                       breaks = seq(from=0.05, to=1, by=0.05),
                       labels =c("0","10","","20","","30","","40","","50",
                                 "","60","","70","","80","","90","","100"))+
    # Scale x axis and set the levels--> WARNING: The maximum distance is 21km!! IS A LIE!
    scale_x_continuous(limits = c(1,20),expand = c(0,0.1),
                       breaks = seq(from=1, to=20, by=1),
                       labels =c("0","2","","4","","6","","8","","10",
                                 "","12","","14","","16","","18","","20"))+
    # Scalate thte colours, spectral means going from red to blue. LIMITS force the range to be btween 0 & 100
    scale_fill_distiller(palette = "Spectral",direction =1, limits = c(0,100))+
    # Labs 
    labs(x="Size (km)",y="Intensity (%)",title = title,subtitle = subtitle)+
    # Sets the contour by Richness. BINWIDTH: the value between each contour line 
    geom_contour(aes(z=Richness), colour="black",size=0.5, lty="solid",stat="contour",
                 lineend="round",linejoin ="round",binwidth = 30)+
    geom_point(aes(x=10.5, y=0.3), col="red", size=9, shape=21 , fill="black",stroke=2)+
    # Determines the labels with a white STROKE below them
    geom_text_contour(aes(z = Richness),stroke = 0.2, size= 10, rotate = F,check_overlap=T)+
    # Theme stuff
    theme_classic()+
    theme(axis.line = element_line(size = 0.05, linetype = "solid"), 
          axis.ticks = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(0.4,1.5,0.4,1.5, "cm"),
          axis.title = element_text(size = 18,face = "bold"), 
          axis.text = element_text(size = 18,face = "bold"), 
          axis.text.x = element_text(size = 17, color = "black"), 
          axis.text.y = element_text(size = 17, color = "black"), 
          legend.position = "none",
          legend.text = element_blank(),
          legend.title = element_blank(),
          panel.background = element_rect(colour = "black", size = 2.5, linetype = "solid"),
          plot.background = element_rect(colour = "black", size = 1, linetype = "solid"), 
          legend.background = element_blank(),
          plot.subtitle = element_text(size = 18,face = "bold", colour = "gray0"),
          plot.title = element_text(size = 18, face = "bold"))
}
Combined_article_plot_Rich_Change_Rate <- function(Post_Val_Matrix_T100, PreFire_Val_Matri, file_name, iterations ){
  plot_Post_fire_1 <- ARTICLE_FOR_GGPLOT_Rich_Change_Rate(Post_Val_Matrix=Post_Val_Matrix_T100[[1]][[1]],
                                                          PreFire_Val_Matri=PreFire_Val_Matri[[1]][[1]], 
                                                          title = iterations,subtitle = "250 m", max=22)  
  plot_Post_fire_2 <- ARTICLE_FOR_GGPLOT_Rich_Change_Rate(Post_Val_Matrix=Post_Val_Matrix_T100[[1]][[2]],
                                                          PreFire_Val_Matri=PreFire_Val_Matri[[1]][[2]], 
                                                          iterations,"500 m", max=22)  
  plot_Post_fire_3 <- ARTICLE_FOR_GGPLOT_Rich_Change_Rate(Post_Val_Matrix=Post_Val_Matrix_T100[[1]][[3]],
                                                          PreFire_Val_Matri=PreFire_Val_Matri[[1]][[3]], 
                                                          iterations,"1000 m", max=22)  
  plot_Post_fire_4 <- ARTICLE_FOR_GGPLOT_Rich_Change_Rate(Post_Val_Matrix=Post_Val_Matrix_T100[[1]][[4]],
                                                          PreFire_Val_Matri=PreFire_Val_Matri[[1]][[4]], 
                                                          iterations,"1500 m", max=22)  
  plot_Post_fire_5 <- ARTICLE_FOR_GGPLOT_Rich_Change_Rate(Post_Val_Matrix=Post_Val_Matrix_T100[[1]][[5]],
                                                          PreFire_Val_Matri=PreFire_Val_Matri[[1]][[5]], 
                                                          iterations,"2000 m", max=22)  
  plot_Post_fire_6 <- ARTICLE_FOR_GGPLOT_Rich_Change_Rate(Post_Val_Matrix=Post_Val_Matrix_T100[[1]][[6]],
                                                          PreFire_Val_Matri=PreFire_Val_Matri[[1]][[6]], 
                                                          iterations,"2500 m", max=22)  
  plot_Post_fire_7 <- ARTICLE_FOR_GGPLOT_Rich_Change_Rate(Post_Val_Matrix=Post_Val_Matrix_T100[[1]][[7]],
                                                          PreFire_Val_Matri=PreFire_Val_Matri[[1]][[7]], 
                                                          iterations,"3000 m", max=22)  
  plot_Post_fire_8 <- ARTICLE_FOR_GGPLOT_Rich_Change_Rate(Post_Val_Matrix=Post_Val_Matrix_T100[[1]][[8]],
                                                          PreFire_Val_Matri=PreFire_Val_Matri[[1]][[8]], 
                                                          iterations,"Percolation distance (3842 m)", max=22)  
  plot_Post_fire_9 <- ARTICLE_FOR_GGPLOT_Rich_Change_Rate(Post_Val_Matrix=Post_Val_Matrix_T100[[1]][[9]],
                                                          PreFire_Val_Matri=PreFire_Val_Matri[[1]][[9]], 
                                                          iterations,"5000 m", max=24)  
  
  ####Els multiplots PER FUNCTIONAL TRAITS AFTER###
  library(gtable)    
  library(grid)
  library(gridExtra) 
  
  # Get the gtables
  gA <- ggplotGrob(plot_Post_fire_1)
  gB <- ggplotGrob(plot_Post_fire_2)
  gC <- ggplotGrob(plot_Post_fire_3)
  gD <- ggplotGrob(plot_Post_fire_4)
  gE <- ggplotGrob(plot_Post_fire_5)
  gF <- ggplotGrob(plot_Post_fire_6)
  gG <- ggplotGrob(plot_Post_fire_7)
  gH <- ggplotGrob(plot_Post_fire_8)
  gI <- ggplotGrob(plot_Post_fire_9)
  
  # Arrange the two charts
  # The legend boxes are centered
  grid.newpage()
  png(filename = paste("C:/Users/Cunilleramontcusi/Desktop/",file_name,".png", sep=""),width=14000,height=12000,units="px",res=800)
  grid.arrange(gA,gB,gC,gD,gE,gF,gG,gH,gI, ncol = 3, nrow=3)
  dev.off()
}

Combined_article_plot_Rich_Change_Rate(Burn_0.1_T100,
                                       Burn_rec_ENTRADA_0.1,
                                       file_name ="Alpha 0.1_DR", iterations = "DR")
Combined_article_plot_Rich_Change_Rate(Burn_0.3_T100,
                                       Burn_rec_ENTRADA_0.3,
                                       file_name ="Alpha 0.3_DR", iterations = "DR")
Combined_article_plot_Rich_Change_Rate(Burn_0.5_T100,
                                       Burn_rec_ENTRADA_0.5,
                                       file_name ="Alpha 0.5_DR", iterations = "DR")
Combined_article_plot_Rich_Change_Rate(Burn_0.7_T100,
                                       Burn_rec_ENTRADA_0.7,
                                       file_name ="Alpha 0.7_DR", iterations = "DR")
Combined_article_plot_Rich_Change_Rate(Burn_0.9_T100,
                                       Burn_rec_ENTRADA_0.9,
                                       file_name ="Alpha 0.9_DR", iterations = "DR")

# Richness at T100_____________________________________________________________________________________________________
#############################################################################################################################
#############################################################################################################################

#Warning! Rerun the functions below before running the rest!
Combined_article_plot <- function(Post_Val_Matrix, file_name, iterations ){
  plot_Post_fire_1 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[1]][[1]],title = iterations,subtitle = "250 m", max = 22,is.7 =1 )  
  plot_Post_fire_2 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[1]][[2]],iterations,"500 m", max = 22,is.7 =2 )  
  plot_Post_fire_3 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[1]][[3]],iterations,"1000 m", max = 22,is.7 =3 )  
  plot_Post_fire_4 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[1]][[4]],iterations,"1500 m", max = 22,is.7 =4 )  
  plot_Post_fire_5 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[1]][[5]],iterations,"2000 m", max = 22,is.7 =5 )  
  plot_Post_fire_6 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[1]][[6]],iterations,"2500 m", max = 22,is.7 =6 )  
  plot_Post_fire_7 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[1]][[7]],iterations,"3000 m", max = 22,is.7 =7 )  
  plot_Post_fire_8 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[1]][[8]],iterations,"Percolation distance (3842 m)", max = 22,is.7 =8 )  
  plot_Post_fire_9 <- ARTICLE_FOR_GGPLOT(Post_Val_Matrix[[1]][[9]],iterations,"5000 m", max = 24,is.7 =9 )  
  
  ####Els multiplots PER FUNCTIONAL TRAITS AFTER###
  library(gtable)    
  library(grid)
  library(gridExtra) 
  
  # Get the gtables
  gA <- ggplotGrob(plot_Post_fire_1)
  gB <- ggplotGrob(plot_Post_fire_2)
  gC <- ggplotGrob(plot_Post_fire_3)
  gD <- ggplotGrob(plot_Post_fire_4)
  gE <- ggplotGrob(plot_Post_fire_5)
  gF <- ggplotGrob(plot_Post_fire_6)
  gG <- ggplotGrob(plot_Post_fire_7)
  gH <- ggplotGrob(plot_Post_fire_8)
  gI <- ggplotGrob(plot_Post_fire_9)
  
  # Arrange the two charts
  # The legend boxes are centered
  grid.newpage()
  png(filename = paste("Figures/",file_name,".png", sep=""),width=14000,height=12000,units="px",res=800)
  grid.arrange(gA,gB,gC,gD,gE,gF,gG,gH,gI, ncol = 3, nrow=3)
  dev.off()
}
ARTICLE_FOR_GGPLOT <- function(Post_Val_Matrix, title, subtitle, max, is.7){
  # Packages needed
  require(ggplot2)
  require(devtools)
  require(metR) # IF not installed: install_github("eliocamp/metR")
  
  df <- as.data.frame((Post_Val_Matrix)) # Make the output a data.frame
  out <- list()                             # Generate and out in list format 
  out.bo <- list()                          # Generate and out in list format
  for(i in 1:length(df)){                   # Loop for each column of the data.frame (that corresponds to Wild. Size)
    out.temp <- df[1:20,i]    # Select all rows in each column and make them a data.frame
    out[[i]] <- out.temp      # Put all the rows in data.frame format into the out list 
  }
  for (e in 1:20) {             # For every 20 times 
    ee<- data.frame(out[[e]],rep(e,times=20),seq(from=0.05, to=1, by=0.05)) # create a new data.frame with the Richness/area/intensity
    colnames(ee) <- c("Richness","area","intensity") # Change column names
    out.bo[[e]] <- ee
  }
  final <- do.call(rbind,out.bo) # Joints the several lists objects created in one big data.frame
  if(is.7==7)#;min(final[,1])>22) # Condition that if the minimum richness is GREATER than 90 will follow the below path
    plot <-ggplot(final,aes(x=area, y=intensity))+
    # Raster with all the data, hjust and vjust hide the bakground and INTERPOLATE= T dissolve the rasters 
    geom_raster(aes(x = area,y=intensity,fill=Richness),hjust=1,vjust=1,interpolate = T)+
    # Scale y axis and set the levels
    scale_y_continuous(limits = c(0.05,1),expand = c(0,0.005), 
                       breaks = seq(from=0.05, to=1, by=0.05),
                       labels =c("0","10","","20","","30","","40","","50",
                                 "","60","","70","","80","","90","","100"))+
    # Scale x axis and set the levels--> WARNING: The maximum distance is 21km!! IS A LIE!
    scale_x_continuous(limits = c(1,20),expand = c(0,0.1),
                       breaks = seq(from=1, to=20, by=1),
                       labels =c("0","2","","4","","6","","8","","10",
                                 "","12","","14","","16","","18","","20"))+
    # Scalate thte colours, spectral means going from red to blue LIMITS force the range to be btween 0 & 100
    scale_fill_distiller(palette = "Spectral",direction =1, limits = c(0,max))+
    # Labs 
    labs(x="Size (km)",y="Intensity (%)",title = title,subtitle = subtitle)+
    geom_contour(aes(z=Richness), colour="black",size=0.5, lty="solid",stat="contour",
                 lineend="round",linejoin ="round",binwidth = 2)+
    geom_point(aes(x=10.5, y=0.3), col="red", size=9, shape=21 , fill="black",stroke=2)+
    # Theme stuff
    theme_classic()+
    theme(axis.line = element_line(size = 0.05, linetype = "solid"), 
          axis.ticks = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(0.4,1.5,0.4,1.5, "cm"),
          axis.title = element_text(size = 18,face = "bold"), 
          axis.text = element_text(size = 18,face = "bold"), 
          axis.text.x = element_text(size = 17, color = "black"), 
          axis.text.y = element_text(size = 17, color = "black"), 
          legend.position = "none",
          legend.text = element_blank(),
          legend.title = element_blank(),
          panel.background = element_rect(colour = "black", size = 2.5, linetype = "solid"),
          plot.background = element_rect(colour = "black", size = 1, linetype = "solid"), 
          legend.background = element_blank(),
          plot.subtitle = element_text(size = 18,face = "bold", colour = "gray0"),
          plot.title = element_text(size = 16, face = "bold"))
  else plot <-ggplot(final,aes(x=area, y=intensity))+ # Condition that if the minimum richness is SMALLER than 90 will follow 
    # Raster with all the data, hjust and vjust hide the bakground and INTERPOLATE= T dissolve the rasters 
    geom_raster(aes(x = area,y=intensity,fill=Richness),hjust=1,vjust=1,interpolate = T)+
    # Scale y axis and set the levels
    scale_y_continuous(limits = c(0.05,1),expand = c(0,0.005),
                       breaks = seq(from=0.05, to=1, by=0.05),
                       labels =c("0","10","","20","","30","","40","","50",
                                 "","60","","70","","80","","90","","100"))+
    # Scale x axis and set the levels--> WARNING: The maximum distance is 21km!! IS A LIE!
    scale_x_continuous(limits = c(1,20),expand = c(0,0.1),
                       breaks = seq(from=1, to=20, by=1),
                       labels =c("0","2","","4","","6","","8","","10",
                                 "","12","","14","","16","","18","","20"))+
    # Scalate thte colours, spectral means going from red to blue. LIMITS force the range to be btween 0 & 100
    scale_fill_distiller(palette = "Spectral",direction =1, limits = c(0,max))+
    # Labs 
    labs(x="Size (km)",y="Intensity (%)",title = title,subtitle = subtitle)+
    # Sets the contour by Richness. BINWIDTH: the value between each contour line 
    geom_contour(aes(z=Richness), colour="black",size=0.5, lty="solid",stat="contour",
                 lineend="round",linejoin ="round",binwidth = 2)+
    geom_point(aes(x=10.5, y=0.3), col="red", size=9, shape=21 , fill="black",stroke=2)+
    # Determines the labels with a white STROKE below them
    geom_text_contour(aes(z = Richness),stroke = 0.2, size= 10, rotate = F,check_overlap=T)+
    # Theme stuff
    theme_classic()+
    theme(axis.line = element_line(size = 0.05, linetype = "solid"), 
          axis.ticks = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(0.4,1.5,0.4,1.5, "cm"),
          axis.title = element_text(size = 18,face = "bold"), 
          axis.text = element_text(size = 18,face = "bold"), 
          axis.text.x = element_text(size = 17, color = "black"), 
          axis.text.y = element_text(size = 17, color = "black"), 
          legend.position = "none",
          legend.text = element_blank(),
          legend.title = element_blank(),
          panel.background = element_rect(colour = "black", size = 2.5, linetype = "solid"),
          plot.background = element_rect(colour = "black", size = 1, linetype = "solid"), 
          legend.background = element_blank(),
          plot.subtitle = element_text(size = 18,face = "bold", colour = "gray0"),
          plot.title = element_text(size = 18, face = "bold"))
}

# Plots for every dispersal ability group (9 plots together) - Figure 6D and Supplementary S4
Combined_article_plot(Burn_0.1_T100, file_name ="Alpha 0.1_T100", iterations = "T100")
Combined_article_plot(Burn_0.3_T100, file_name ="Alpha 0.3_T100", iterations = "T100")
Combined_article_plot(Burn_0.5_T100, file_name ="Alpha 0.5_T100", iterations = "T100")
Combined_article_plot(Burn_0.7_T100, file_name ="Alpha 0.7_T100", iterations = "T100")
Combined_article_plot(Burn_0.9_T100, file_name ="Alpha 0.9_T100", iterations = "T100")

# Plots the whole metacommunity (all dispersal abilities together) - Figure 6C and Supplementary S4
png(filename = "Figures/Globalplot_01_T100.png",width=14000,height=12000,units="px",res=1000)
GLOBAL_ARTICLE_FOR_GGPLOT(Burn_0.1_T100,"Entire metacommunity T100","Shared = 0.1")
dev.off()
png(filename = "Figures/Globalplot_03_T100.png",width=14000,height=12000,units="px",res=1000)
GLOBAL_ARTICLE_FOR_GGPLOT(Burn_0.3_T100,"Entire metacommunity T100","Shared = 0.3")
dev.off()
png(filename = "Figures/Globalplot_05_T100.png",width=14000,height=12000,units="px",res=1000)
GLOBAL_ARTICLE_FOR_GGPLOT(Burn_0.5_T100,"Entire metacommunity T100","Shared = 0.5")
dev.off()
png(filename = "Figures/Globalplot_07_T100.png",width=14000,height=12000,units="px",res=1000)
GLOBAL_ARTICLE_FOR_GGPLOT(Burn_0.7_T100,"Entire metacommunity T100","Shared = 0.7")
dev.off()
png(filename = "Figures/Globalplot_09_T100.png",width=14000,height=12000,units="px",res=1000)
GLOBAL_ARTICLE_FOR_GGPLOT(Burn_0.9_T100,"Entire metacommunity T100","Shared = 0.9")
dev.off()

## Mean richness SD (Supplementary S6)
Combined_article_plot_sd <- function(Post_Val_Matrix, file_name, iterations ){
  plot_Post_fire_1 <- ARTICLE_FOR_GGPLOT_sd(Post_Val_Matrix[[3]][[1]],iterations,"250 m")  
  plot_Post_fire_2 <- ARTICLE_FOR_GGPLOT_sd(Post_Val_Matrix[[3]][[2]],iterations,"500 m")  
  plot_Post_fire_3 <- ARTICLE_FOR_GGPLOT_sd(Post_Val_Matrix[[3]][[3]],iterations,"1000 m")  
  plot_Post_fire_4 <- ARTICLE_FOR_GGPLOT_sd(Post_Val_Matrix[[3]][[4]],iterations,"1500 m")  
  plot_Post_fire_5 <- ARTICLE_FOR_GGPLOT_sd(Post_Val_Matrix[[3]][[5]],iterations,"2000 m")  
  plot_Post_fire_6 <- ARTICLE_FOR_GGPLOT_sd(Post_Val_Matrix[[3]][[6]],iterations,"2500 m")  
  plot_Post_fire_7 <- ARTICLE_FOR_GGPLOT_sd(Post_Val_Matrix[[3]][[7]],iterations,"3000 m")  
  plot_Post_fire_8 <- ARTICLE_FOR_GGPLOT_sd(Post_Val_Matrix[[3]][[8]],iterations,"Percolation distance (3842 m)")  
  plot_Post_fire_9 <- ARTICLE_FOR_GGPLOT_sd(Post_Val_Matrix[[3]][[9]],iterations,"5000 m")  
  
  ####Els multiplots PER FUNCTIONAL TRAITS AFTER###
  library(gtable)    
  library(grid)
  library(gridExtra) 
  
  # Get the gtables
  gA <- ggplotGrob(plot_Post_fire_1)
  gB <- ggplotGrob(plot_Post_fire_2)
  gC <- ggplotGrob(plot_Post_fire_3)
  gD <- ggplotGrob(plot_Post_fire_4)
  gE <- ggplotGrob(plot_Post_fire_5)
  gF <- ggplotGrob(plot_Post_fire_6)
  gG <- ggplotGrob(plot_Post_fire_7)
  gH <- ggplotGrob(plot_Post_fire_8)
  gI <- ggplotGrob(plot_Post_fire_9)
  
  # Arrange the two charts
  # The legend boxes are centered
  grid.newpage()
  png(filename = paste("Figures/",file_name,".png", sep=""),width=14000,height=12000,units="px",res=800)
  grid.arrange(gA,gB,gC,gD,gE,gF,gG,gH,gI, ncol = 3, nrow=3)
  dev.off()
}
ARTICLE_FOR_GGPLOT_sd <- function(Post_Val_Matrix, title, subtitle){
  # Packages needed
  require(ggplot2)
  require(devtools)
  require(metR) # IF not installed: install_github("eliocamp/metR")
  
  df <- as.data.frame((Post_Val_Matrix)) # Make the output a data.frame
  out <- list()                             # Generate and out in list format 
  out.bo <- list()                          # Generate and out in list format
  for(i in 1:length(df)){                   # Loop for each column of the data.frame (that corresponds to Wild. Size)
    out.temp <- df[1:20,i]    # Select all rows in each column and make them a data.frame
    out[[i]] <- out.temp      # Put all the rows in data.frame format into the out list 
  }
  for (e in 1:20) {             # For every 20 times 
    ee<- data.frame(out[[e]],rep(e,times=20),seq(from=0.05, to=1, by=0.05)) # create a new data.frame with the Richness/area/intensity
    colnames(ee) <- c("Richness","area","intensity") # Change column names
    out.bo[[e]] <- ee
  }
  final <- do.call(rbind,out.bo) # Joints the several lists objects created in one big data.frame
  final[which(is.na(final)),1] <- 0
  
  #return(final)
  if(min(final[,1])>90) # Condition that if the minimum richness is GREATER than 90 will follow the below path
    plot <-ggplot(final,aes(x=area, y=intensity))+
    # Raster with all the data, hjust and vjust hide the bakground and INTERPOLATE= T dissolve the rasters 
    geom_raster(aes(x = area,y=intensity,fill=Richness),hjust=1,vjust=1,interpolate = T)+
    # Scale y axis and set the levels
    scale_y_continuous(limits = c(0.05,1),expand = c(0,0.005), 
                       breaks = seq(from=0.05, to=1, by=0.05),
                       labels =c("0","10","","20","","30","","40","","50",
                                 "","60","","70","","80","","90","","100"))+
    # Scale x axis and set the levels--> WARNING: The maximum distance is 21km!! IS A LIE!
    scale_x_continuous(limits = c(1,20),expand = c(0,0.1),
                       breaks = seq(from=1, to=20, by=1),
                       labels =c("0","2","","4","","6","","8","","10",
                                 "","12","","14","","16","","18","","20"))+
    # Scalate thte colours, spectral means going from red to blue LIMITS force the range to be btween 0 & 100
    scale_fill_distiller(palette = "Spectral",direction =1, limits = c(0,100))+
    # Labs 
    labs(x="Size (km)",y="Intensity (%)",title = title,subtitle = subtitle)+
    # Theme stuff
    theme_classic()+
    theme(axis.line = element_line(size = 0.05, linetype = "solid"), 
          axis.ticks = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(0.4,1.5,0.4,1.5, "cm"),
          axis.title = element_text(size = 18,face = "bold"), 
          axis.text = element_text(size = 18,face = "bold"), 
          axis.text.x = element_text(size = 17, color = "black"), 
          axis.text.y = element_text(size = 17, color = "black"), 
          legend.position = "none",
          legend.text = element_blank(),
          legend.title = element_blank(),
          panel.background = element_rect(colour = "black", size = 2.5, linetype = "solid"),
          plot.background = element_rect(colour = "black", size = 1, linetype = "solid"), 
          legend.background = element_blank(),
          plot.subtitle = element_text(size = 18,face = "bold", colour = "gray0"),
          plot.title = element_text(size = 16, face = "bold"))
  else plot <-ggplot(final,aes(x=area, y=intensity))+ # Condition that if the minimum richness is SMALLER than 90 will follow 
    # Raster with all the data, hjust and vjust hide the bakground and INTERPOLATE= T dissolve the rasters 
    geom_raster(aes(x = area,y=intensity,fill=Richness),hjust=1,vjust=1,interpolate = T)+
    # Scale y axis and set the levels
    scale_y_continuous(limits = c(0.05,1),expand = c(0,0.005),
                       breaks = seq(from=0.05, to=1, by=0.05),
                       labels =c("0","10","","20","","30","","40","","50",
                                 "","60","","70","","80","","90","","100"))+
    # Scale x axis and set the levels--> WARNING: The maximum distance is 21km!! IS A LIE!
    scale_x_continuous(limits = c(1,20),expand = c(0,0.1),
                       breaks = seq(from=1, to=20, by=1),
                       labels =c("0","2","","4","","6","","8","","10",
                                 "","12","","14","","16","","18","","20"))+
    # Scalate thte colours, spectral means going from red to blue. LIMITS force the range to be btween 0 & 100
    scale_fill_distiller(palette = "BrBG",direction =-1, limits = c(0,20))+
    # Labs 
    labs(x="Size (km)",y="Intensity (%)",title = title, subtitle = subtitle)+
    # Sets the contour by Richness. BINWIDTH: the value between each contour line 
    geom_contour(aes(z=Richness), colour="black",size=0.5, lty="solid",stat="contour",
                 lineend="round",linejoin ="round",binwidth = 1)+
    # Determines the labels with a white STROKE below them
    geom_text_contour(aes(z = Richness),stroke = 0.2, size= 15, rotate = F,check_overlap=T)+
    # Theme stuff
    theme_classic()+
    theme(axis.line = element_line(size = 0.05, linetype = "solid"), 
          axis.ticks = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(0.4,1.5,0.4,1.5, "cm"),
          axis.title = element_text(size = 18,face = "bold"), 
          axis.text = element_text(size = 18,face = "bold"), 
          axis.text.x = element_text(size = 17, color = "black"), 
          axis.text.y = element_text(size = 17, color = "black"), 
          legend.position = "none",
          legend.text = element_blank(),
          legend.title = element_blank(),
          panel.background = element_rect(colour = "black", size = 2.5, linetype = "solid"),
          plot.background = element_rect(colour = "black", size = 1, linetype = "solid"), 
          legend.background = element_blank(),
          plot.subtitle = element_text(size = 18,face = "bold", colour = "gray0"),
          plot.title = element_text(size = 18, face = "bold"))
}

Combined_article_plot_sd(Burn_0.1_T100, file_name ="SD_Alpha 0.1_T100", iterations = "T100")
Combined_article_plot_sd(Burn_0.3_T100, file_name ="SD_Alpha 0.3_T100", iterations = "T100")
Combined_article_plot_sd(Burn_0.5_T100, file_name ="SD_Alpha 0.5_T100", iterations = "T100")
Combined_article_plot_sd(Burn_0.7_T100, file_name ="SD_Alpha 0.7_T100", iterations = "T100")
Combined_article_plot_sd(Burn_0.9_T100, file_name ="SD_Alpha 0.9_T100", iterations = "T100")

## Mean iterations ti reach burned ponds (Supplemenetary S5)

Combined_article_plot_iter <- function(Post_Val_Matrix, file_name, iterations ){
  plot_Post_fire_1 <- ARTICLE_FOR_GGPLOT_iter(Post_Val_Matrix[[2]][[1]],iterations,"250 m")  
  plot_Post_fire_2 <- ARTICLE_FOR_GGPLOT_iter(Post_Val_Matrix[[2]][[2]],iterations,"500 m")  
  plot_Post_fire_3 <- ARTICLE_FOR_GGPLOT_iter(Post_Val_Matrix[[2]][[3]],iterations,"1000 m")  
  plot_Post_fire_4 <- ARTICLE_FOR_GGPLOT_iter(Post_Val_Matrix[[2]][[4]],iterations,"1500 m")  
  plot_Post_fire_5 <- ARTICLE_FOR_GGPLOT_iter(Post_Val_Matrix[[2]][[5]],iterations,"2000 m")  
  plot_Post_fire_6 <- ARTICLE_FOR_GGPLOT_iter(Post_Val_Matrix[[2]][[6]],iterations,"2500 m")  
  plot_Post_fire_7 <- ARTICLE_FOR_GGPLOT_iter(Post_Val_Matrix[[2]][[7]],iterations,"3000 m")  
  plot_Post_fire_8 <- ARTICLE_FOR_GGPLOT_iter(Post_Val_Matrix[[2]][[8]],iterations,"Percolation distance (3842 m)")  
  plot_Post_fire_9 <- ARTICLE_FOR_GGPLOT_iter(Post_Val_Matrix[[2]][[9]],iterations,"5000 m")  
  
  ####Els multiplots PER FUNCTIONAL TRAITS AFTER###
  library(gtable)    
  library(grid)
  library(gridExtra) 
  
  # Get the gtables
  gA <- ggplotGrob(plot_Post_fire_1)
  gB <- ggplotGrob(plot_Post_fire_2)
  gC <- ggplotGrob(plot_Post_fire_3)
  gD <- ggplotGrob(plot_Post_fire_4)
  gE <- ggplotGrob(plot_Post_fire_5)
  gF <- ggplotGrob(plot_Post_fire_6)
  gG <- ggplotGrob(plot_Post_fire_7)
  gH <- ggplotGrob(plot_Post_fire_8)
  gI <- ggplotGrob(plot_Post_fire_9)
  
  # Arrange the two charts
  # The legend boxes are centered
  grid.newpage()
  png(filename = paste("Figures/",file_name,".png", sep=""),width=14000,height=12000,units="px",res=800)
  grid.arrange(gA,gB,gC,gD,gE,gF,gG,gH,gI, ncol = 3, nrow=3)
  dev.off()
}
ARTICLE_FOR_GGPLOT_iter <- function(Post_Val_Matrix, title, subtitle){
  # Packages needed
  require(ggplot2)
  require(devtools)
  require(metR) # IF not installed: install_github("eliocamp/metR")
  
  df <- as.data.frame((Post_Val_Matrix)) # Make the output a data.frame
  out <- list()                             # Generate and out in list format 
  out.bo <- list()                          # Generate and out in list format
  for(i in 1:length(df)){                   # Loop for each column of the data.frame (that corresponds to Wild. Size)
    out.temp <- df[1:20,i]    # Select all rows in each column and make them a data.frame
    out[[i]] <- out.temp      # Put all the rows in data.frame format into the out list 
  }
  for (e in 1:20) {             # For every 20 times 
    ee<- data.frame(out[[e]],rep(e,times=20),seq(from=0.05, to=1, by=0.05)) # create a new data.frame with the Richness/area/intensity
    colnames(ee) <- c("Richness","area","intensity") # Change column names
    out.bo[[e]] <- ee
  }
  final <- do.call(rbind,out.bo) # Joints the several lists objects created in one big data.frame  
  if(mean(final[,1])<5) # Condition that if the minimum richness is GREATER than 90 will follow the below path
    plot <-ggplot(final,aes(x=area, y=intensity))+
    # Raster with all the data, hjust and vjust hide the bakground and INTERPOLATE= T dissolve the rasters 
    geom_raster(aes(x = area,y=intensity,fill=Richness),hjust=1,vjust=1,interpolate = T)+
    # Scale y axis and set the levels
    scale_y_continuous(limits = c(0.05,1),expand = c(0,0.005), 
                       breaks = seq(from=0.05, to=1, by=0.05),
                       labels =c("0","10","","20","","30","","40","","50",
                                 "","60","","70","","80","","90","","100"))+
    # Scale x axis and set the levels--> WARNING: The maximum distance is 21km!! IS A LIE!
    scale_x_continuous(limits = c(1,20),expand = c(0,0.1),
                       breaks = seq(from=1, to=20, by=1),
                       labels =c("0","2","","4","","6","","8","","10",
                                 "","12","","14","","16","","18","","20"))+
    # Scalate thte colours, spectral means going from red to blue LIMITS force the range to be btween 0 & 100
    scale_fill_distiller(palette = "PuOr",direction =-1, limits = c(0,100))+
    # Labs 
    labs(x="Size (km)",y="Intensity (%)",title = title,subtitle = subtitle)+
    geom_contour(aes(z=Richness), colour="black",size=0.5, lty="solid",stat="contour",
                 lineend="round",linejoin ="round",binwidth = 10)+
    geom_point(aes(x=10.5, y=0.3), col="red", size=9, shape=21 , fill="black",stroke=2)+
    # Theme stuff
    theme_classic()+
    theme(axis.line = element_line(size = 0.05, linetype = "solid"), 
          axis.ticks = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(0.4,1.5,0.4,1.5, "cm"),
          axis.title = element_text(size = 18,face = "bold"), 
          axis.text = element_text(size = 18,face = "bold"), 
          axis.text.x = element_text(size = 17, color = "black"), 
          axis.text.y = element_text(size = 17, color = "black"), 
          legend.position = "none",
          legend.text = element_blank(),
          legend.title = element_blank(),
          panel.background = element_rect(colour = "black", size = 2.5, linetype = "solid"),
          plot.background = element_rect(colour = "black", size = 1, linetype = "solid"), 
          legend.background = element_blank(),
          plot.subtitle = element_text(size = 18,face = "bold", colour = "gray0"),
          plot.title = element_text(size = 16, face = "bold"))
  else plot <-ggplot(final,aes(x=area, y=intensity))+ # Condition that if the minimum richness is SMALLER than 90 will follow 
    # Raster with all the data, hjust and vjust hide the bakground and INTERPOLATE= T dissolve the rasters 
    geom_raster(aes(x = area,y=intensity,fill=Richness),hjust=1,vjust=1,interpolate = T)+
    # Scale y axis and set the levels
    scale_y_continuous(limits = c(0.05,1),expand = c(0,0.005),
                       breaks = seq(from=0.05, to=1, by=0.05),
                       labels =c("0","10","","20","","30","","40","","50",
                                 "","60","","70","","80","","90","","100"))+
    # Scale x axis and set the levels--> WARNING: The maximum distance is 21km!! IS A LIE!
    scale_x_continuous(limits = c(1,20),expand = c(0,0.1),
                       breaks = seq(from=1, to=20, by=1),
                       labels =c("0","2","","4","","6","","8","","10",
                                 "","12","","14","","16","","18","","20"))+
    # Scalate thte colours, spectral means going from red to blue. LIMITS force the range to be btween 0 & 100
    scale_fill_distiller(palette = "PuOr",direction =-1, limits = c(0,100))+
    # Labs 
    labs(x="Size (km)",y="Intensity (%)",title = title,subtitle = subtitle)+
    # Sets the contour by Richness. BINWIDTH: the value between each contour line 
    geom_contour(aes(z=Richness), colour="black",size=0.5, lty="solid",stat="contour",
                 lineend="round",linejoin ="round",binwidth = 10)+
    geom_point(aes(x=10.5, y=0.3), col="red", size=9, shape=21 , fill="black",stroke=2)+
    # Determines the labels with a white STROKE below them
    geom_text_contour(aes(z = Richness),stroke = 0.2, size= 10, rotate = F,check_overlap=T)+
    # Theme stuff
    theme_classic()+
    theme(axis.line = element_line(size = 0.05, linetype = "solid"), 
          axis.ticks = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.margin = margin(0.4,1.5,0.4,1.5, "cm"),
          axis.title = element_text(size = 18,face = "bold"), 
          axis.text = element_text(size = 18,face = "bold"), 
          axis.text.x = element_text(size = 17, color = "black"), 
          axis.text.y = element_text(size = 17, color = "black"), 
          legend.position = "none",
          legend.text = element_blank(),
          legend.title = element_blank(),
          panel.background = element_rect(colour = "black", size = 2.5, linetype = "solid"),
          plot.background = element_rect(colour = "black", size = 1, linetype = "solid"), 
          legend.background = element_blank(),
          plot.subtitle = element_text(size = 18,face = "bold", colour = "gray0"),
          plot.title = element_text(size = 18, face = "bold"))
}

Combined_article_plot_iter(Burn_0.1_T100, file_name ="IT_Alpha 0.1_T100", iterations = "T100")
Combined_article_plot_iter(Burn_0.3_T100, file_name ="IT_Alpha 0.3_T100", iterations = "T100")
Combined_article_plot_iter(Burn_0.5_T100, file_name ="IT_Alpha 0.5_T100", iterations = "T100")
Combined_article_plot_iter(Burn_0.7_T100, file_name ="IT_Alpha 0.7_T100", iterations = "T100")
Combined_article_plot_iter(Burn_0.9_T100, file_name ="IT_Alpha 0.9_T100", iterations = "T100")


###################################################################################################################################
### END ###########################################################################################################################
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################
###################################################################################################################################

