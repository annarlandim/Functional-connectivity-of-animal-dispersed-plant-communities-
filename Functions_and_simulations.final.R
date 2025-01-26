library("vegan") 
library("Hmisc")
library("dplyr")

#### Functions ####

### Adapted: 
## Adapted from Donoso et al. 2017

# Web to estimate fruit preference in the source patch (rowSums = 1):
makeweb_r <- function(specpar = 1, birdtraits, planttraits, nicheshape="normal"){
  
  fun_pref <- function(traitdif){
    if (nicheshape == "normal"){
      prefs <- dnorm(traitdif, mean = 0, sd=1/specpar)
    }
    if (nicheshape == "skewed"){
      prefs <- dlnorm(traitdif * specpar + exp (-1)) 
    }
    prefs
  }
  Nplant <- length(planttraits)
  Nbird <- length(birdtraits)
  web <- fun_pref(outer(planttraits, birdtraits, "-")*(-1)) 
  web
}

# Web to estimate preference to plants in target patch (colSums = 1):
makeweb_c <- function(specpar = 1, birdtraits, planttraits, nicheshape="normal"){
  
  fun_pref <- function(traitdif){
    if (nicheshape == "normal"){
      prefs <- dnorm(traitdif, mean = 0, sd=1/specpar)
    }
    if (nicheshape == "skewed"){
      prefs <- dlnorm(traitdif * specpar + exp (-1)) 
    }
    prefs
  }
  Nplant <- length(planttraits)
  Nbird <- length(birdtraits)
  web <- fun_pref(outer(planttraits, birdtraits, "-")*(-1)) 
  web
}

# Trait distributions:
get_skewtr <- function(myN, tr_meanlog = 2, tr_sdlog = 1.5){
  tr <- qlnorm(seq(0, 1, length.out = myN+2), tr_meanlog, tr_sdlog)[-c(1, myN+2)]
}

# Interaction network to define dispersal vectors:
# (uses preference web and interaction frequencies)
make_trueweb_r <- function(web_p, plantfreq, birdfreq){
  web_relfreq <- (plantfreq %*% t(birdfreq)) / mean(birdfreq)
  web_p * web_relfreq
}

# Interaction network to estimate attractiveness to target patch:
# (uses preference web and interaction frequencies) 
make_trueweb_c <- function(web_p, plantfreq){
  # Since we use this to estimate attractiveness for each bird, 
  # bird frequencies are always = 1
  web_p * plantfreq
}

## Original from Donoso et al. (2017):

# Function used to sample the interaction networks:
sampleweb <- function(web,obs=NULL, int_freq, method='perplant'){
  if(method=='perplant') {
    if(length(obs)==1){obs <- obs*int_freq}
    if(length(obs)!=nrow(web)){stop("length of obs neither 1 nor matching number of
                                           species")}
    sampledweb <- sapply(1:nrow(web), function(i) {
      web.i <- as.numeric(web[i,])
      table(sample(factor(1:length(web.i)), obs[i], prob = web.i, replace = TRUE))
    })
  }
  if(method=='percons') {
    if (length(obs)==1){obs <- rep(obs,ncol(web))}
    if (length(obs)!=ncol(web)){stop("length of obs neither 1 nor matching number of
                                           species")}
    sampledweb <- sapply(1:ncol(web), function(j) {
      web.j <- web[,j]
      table(sample(factor(1:length(web.j)),obs[j],prob=web.j,replace=TRUE)) 
    })
  }
  if(method=='perweb'){
    if (length(obs)!=1){stop("obs must have length 1 with method 'perweb'")}
    Nobs <- round(obs * ncol(web))
    sampledweb.vect <- sample(as.factor(1:length(web)), size=Nobs, prob=as.numeric(web),
                              replace=TRUE)
    sampledweb <- matrix(as.numeric(table(sampledweb.vect)), nrow=nrow(web)) 
  }
  if(method =='rarefy'){
    if(is.null(obs)) obs <- min(colSums(web)) 
    if(obs<1){stop("cannot rarefy to less than 1 observation")}
    sampledweb <- sapply(1:ncol(web), function(j) {
      web.j <- web[,j]
      24
      table(sample(rep(factor(1:length(web.j)),web.j),obs,replace=FALSE)) 
    })
  }
  dimnames(sampledweb) <- dimnames(web)
  sampledweb
}

## Original from Sorensen et al. 2020

# Function to convert a web interaction network back to a single frame (column of plant and column of birds) 
# one row for each interaction, equivalent function for the frame2web function in the bipartite package
web2frame <- function(x){
  newplant <- NULL
  newbird <- NULL
  new20 <- NULL
  for (i in 1:nrow(x)) {
    for (ii in 1:ncol(x)){
      newplant <- rep(rownames(x[i,]), each = x[i,ii])
      newbird <- rep(colnames(x[ii]), each = x[i,ii])
      new10 <- cbind(newplant,newbird)
      new20 <- rbind(new20, new10)
    }
    newplant <- NULL
    newbird <- NULL
    new10 <- NULL
  }
  return(new20)
}

# Function to create a single data frame with plants and birds traits:
matchbirdplanttraits <- function(simulatedframes, t_obs, birdtraits, planttraits) {
  bird.trait <- rep(NA, t_obs) 
  fruit.trait <- rep(NA, t_obs) 
  bird.trait2 <- rep(NA, t_obs)
  if(length(birdtraits) == 3){
    for(i in 1:nrow(simulatedframes)) {
      ti <- match(simulatedframes[i,2], birdtraits$birdnumber)
      z <- match(simulatedframes[i,1], planttraits$plantnumber)
      bird.trait[i] <- birdtraits[ti,1]
      fruit.trait[i] <- planttraits[z,1]
      bird.trait2[i] <- birdtraits[ti,3]
    }
    simulatedframes <- cbind(simulatedframes, bird.trait, fruit.trait, bird.trait2)
  } else {
    for(i in 1:nrow(simulatedframes)) {
      ti <- match(simulatedframes[i,2], birdtraits$birdnumber)
      z <- match(simulatedframes[i,1], planttraits$plantnumber)
      bird.trait[i] <- birdtraits[ti,1]
      fruit.trait[i] <- planttraits[z,1]
    }
    simulatedframes <- cbind(simulatedframes, bird.trait, fruit.trait)
  }}

# Simulation of dispersal distances:
dispsimulation <- function (x, t_obs) {
  dispdist <- rep(NA, t_obs) 
  for(i in 1:nrow(x)) {
    meanGPThour <- 4.5*x[i,3]^0.5 
    meanGPT <- meanGPThour*3600
    scalevalue <- 75311/meanGPT
    shapevalue <- meanGPT^2/75311
    GPT <- rgamma(1, shape = shapevalue, scale = scalevalue)
    meanspeed <- 15.7*x[i,3]^0.17
    speed <- rnorm(1, meanspeed, 2.078)  
    max_distance <- speed*GPT
    distance <- 0.002 * max_distance 
    if (is.nan (distance)){
      distance<-NA
    }
    dispdist[i] <- distance
  }
  x <- cbind(x, dispdist)
}

### New:
## Estimating probability of seed dispersal success:

distancefilter <- function(patches, Nplant, obs, simulatedframe){
  
  dist_filt <- array(dim = c(length(patches), max(obs), Nplant), dimnames = list(paste0("patch_", patches), paste0("event_", 1:max(obs)), paste0("plant_", 1:Nplant)))
  
  for(p in 1:length(patches)){
    
    d <- patches[p]
    
    for(j in 1:Nplant){
      
      disp_dist <- t(matrix(simulatedframe[simulatedframe$newplant == j, "dispdist"]))
      
      for(i in 1:length(disp_dist)){
        if(disp_dist[i]-d < 0){
          disp_dist[i] <- 0
        } else {
          disp_dist[i] <- 1
        }
      }
      empty <- rep(NA, times = max(obs)-length(disp_dist))
      dist_filt[p,,j] <- c(disp_dist, empty)
    }
  }
  return(dist_filt)
}

## Estimating relative attractiveness of target patch:

# Random plant assembly

r_webs_rand <- function(specsparlevel, Nplant, Nbird, Nrep, rec_webs_p){
  alldata <- vector("list", length = length(specsparlevel))
  names(alldata) <- paste("Specialization.level", c(specsparlevel), sep="_")
  
  for(s in 1:length(specsparlevel)){
    specspar_data <- vector("list", length = Nplant+1)
    
    for(k in 1:(Nplant+1)){ 
      specspar_data[[k]] <- array(NA, dim = c((Nplant-k+1), Nbird, Nrep))
    }
    
    names(specspar_data) <- paste('Nplant.recipient',Nplant:0, sep='_')
    
    for(n in 1:Nrep){
      web.old <- rec_webs_p[[s]]
      
      seq.ran <- sample(rownames(rec_webs_p[[s]]))
      
      for(k in 1:(Nplant)){
        specspar_data[[k]][,,n] <- web.old
        web.old <- web.old[-which(rownames(web.old)==seq.ran[k]),, drop=FALSE]
      }
    }
    specspar_data[[(Nplant+1)]] <- array(0, dim = c(1, Nbird, Nrep))
    alldata[[s]] <- specspar_data
  }
  return(alldata)
}

# Relative attractiveness per plant richness in target patch: 

attractivenessfilter <- function(simulatedframes, source_webs_p, recipient_webs, specsparlevel, obs, Nplant, Nrep, Ncons, plant_freq){
  
  attrac_filt <- vector("list", length(specsparlevel))
  names(attrac_filt) <- paste("Specialization_level", c(specsparlevel), sep = "_")
  
  for(s in 1:length(specsparlevel)){
    
    attrac_filt_ww <- array(dim = c(length(recipient_webs[[s]]), max(obs), Nplant, Nrep), 
                            dimnames = list(paste0("re_p_comm_", 1:(length(recipient_webs[[s]]))), 
                                            paste0("event_", 1:max(obs)), paste0("plant_", 1:Nplant),
                                            paste0("iteration_",1:Nrep)))
    
    for(w in 1:Nplant){
      
      rec <- colSums(recipient_webs[[s]][[w]])
      sourc_rec <- (colSums(plant_freq*source_webs_p[[s]]) + colSums(recipient_webs[[s]][[w]]))
      attrac <- as.data.frame(rec/sourc_rec)
      attrac$cons_id <- c(1:Ncons)
      row.names(attrac) <- NULL
      
      for(j in 1:Nplant){
        
        simulatedframes_s <- simulatedframes[[s]]
        attrac_plantcomm <- array(t(matrix(simulatedframes_s[simulatedframes_s$newplant == j, "newbird"])), dim = c(1, obs[j], Nrep))
        
        attrac_plantcomm_exp <- array(NA, dim = c(1, max(obs), Nrep))
        attrac_plantcomm_exp[, 1:obs[j], ] <- attrac_plantcomm
        attrac_plantcomm <- attrac_plantcomm_exp
        
        for(i in 1:obs[j]){
          for(n in 1:Nrep){
            id <- match(attrac_plantcomm[,i,n], attrac$cons_id)
            attrac_plantcomm[,i,n] <- attrac[id, n]
          }
        }
        attrac_filt_ww[w,,j,] <- as.numeric(attrac_plantcomm)
      }
    }
    
    for(j in 1:Nplant){
      attrac_filt_ww[Nplant+1,,j,] <- array(0, dim = c(1,max(obs),n))
    }
    attrac_filt[[s]] <- attrac_filt_ww
  }
  return(attrac_filt)
}

## Combining distance and direction:

# For all species:
reach_prob_webs <- function(recipient_webs, specsparlevel, patches, obs, Nplant, Nrep,
                            filter_distance_webs, filter_attrac_webs){
  
  prob_webs <- vector("list", (Nplant+1))
  
  for(w in 1:(Nplant+1)){
    
    probs_web_w <- vector("list", length(specsparlevel))
    
    for(s in 1:length(specsparlevel)){
      
      probs_s <- array(dim = c(length(patches), max(obs), Nplant, Nrep), 
                       dimnames = list(paste0("patch_", patches), 
                                       paste0("event_", 1:max(obs)), 
                                       paste0("plant_", 1:Nplant),
                                       paste0("iteration_", 1:Nrep)))
      
      for(i in 1:Nplant){
        for(n in 1:Nrep){
          prob_plant_i <- sweep(filter_distance_webs[[s]][,,i], 2, filter_attrac_webs[[s]][w,,i,n], "*")
          
          probs_s[,,i,n] <- prob_plant_i
        }
      }
      probs_web_w[[s]] <- probs_s
    }
    names(probs_web_w) <- paste("Specialization_level", c(specsparlevel), sep = "_")
    prob_webs[[w]] <- probs_web_w  
  }
  names(prob_webs) <- paste("Plant_Rich", Nplant:0, sep = "_")
  return(prob_webs)
}

# Combining species result into community of seeds dispersed to target patches:
seeds_reaching_webs <- function(reach_prob_webs, specsparlevel, Nplant, Nrep, patches){
  
  n_seeds_webs <- vector("list", length(specsparlevel))
  
  for(s in 1:length(specsparlevel)){
    
    n_seeds <- array(NA, dim = c(length(patches), Nplant, Nrep),  
                     dimnames = list(paste0("Patch_", patches),
                                     paste0("Plant_",c(1:Nplant)), paste0("Iteration_", c(1:Nrep))))
    
    reach_prob <- reach_prob_webs[[s]]
    
    for(i in 1:Nplant){
      for(n in 1:Nrep){
        n_seeds_plant_i <- reach_prob[,,i,n]
        n_seeds_plant_i <- rowSums(n_seeds_plant_i, na.rm = TRUE)
        n_seeds[,i,n] <- n_seeds_plant_i 
      }
    }   
    n_seeds_webs[[s]] <- n_seeds
  }
  names(n_seeds_webs) <- paste("Specialization_level", c(specsparlevel), sep = "_")
  return(n_seeds_webs)
}

# Number and diversity of reaching seeds:
dispersal_dimensions <- function(Nplant, specsparlevel, Nrep, patches, n_seeds_reaching_webs, plant_trait, recipient_webs, plant_freq, obsperplant){
  Dispersed <- vector("list", length = (Nplant+1))
  
  for(w in 1:(Nplant+1)){
    
    Dispersed_webs <- vector("list", length = length(specsparlevel))
    
    for(s in 1:length(specsparlevel)){
      
      disp <- vector("list", length = Nrep)
      
      for (n in 1:Nrep){ 
        disp_it <- as.data.frame(c(patches))
        colnames(disp_it)[1] <- "Distance"
        disp_it$Spec.par.level <- as.factor(rep(specsparlevel[s], times = length(patches)))
        # Richness_recipient
        if(w==(Nplant+1)){
          disp_it$Rich_recipient <- rep((0), times = length(patches))
        } else {disp_it$Rich_recipient <- rep(as.numeric(c(Nplant:0)[w]), times = length(patches))}
        # Estimations
        disp_it$Abundance <- rowSums(n_seeds_reaching_webs[[w]][[s]][,,n]) 
        disp_it$Diversity <- exp(diversity(n_seeds_reaching_webs[[w]][[s]][,,n], index="shannon"))
        disp[[n]] <- disp_it  
      } 
      Dispersed_webs[[s]] <- bind_rows(disp)
    }
    Dispersed[[w]] <- bind_rows(Dispersed_webs)
  }
  return(bind_rows(Dispersed))
}

### Simulations ####

# Empirical plant community:
Plants <- read.csv("Manu_plant.community.csv")

# Empirical bird community:
Birds <- read.csv("Manu_bird.community.csv")

Nbird <- 60
Nplant <- 50

### Get mean and standard deviation from empirical communities:

## Traits related to birds' fruit preferences:

# Fruit width
p_fw_meanlog <- mean(log(Plants$Fruit.D1.mm))
p_fw_sdlog <- sd(log(Plants$Fruit.D1.mm))

# Bill width
b_bw_meanlog <- mean(log(Birds$Bill.width.mm)) 
b_bw_sdlog <- sd(log(Birds$Bill.width.mm))

## Traits related to primarily used forest strat:

# Height
p_h_meanlog <- mean(log(Plants$Plant.height.m))
p_h_sdlog <- sd(log(Plants$Plant.height.m))

# Wing Pointdness
b_wp_meanlog <- mean(log(Birds$Kipp_index))
b_wp_sdlog <- sd(log(Birds$Kipp_index))

### Estimate distributions with mean and standard deviation from the empirical communities:

## Traits related to birds' fruit preferences:

# Fruit width
fit_p_fw <- get_skewtr(Nplant, p_fw_meanlog, p_fw_sdlog) 

# Bill width
fit_b_bw <- get_skewtr(Nbird, b_bw_meanlog, b_bw_sdlog)

## Traits related to primarily used forest strat:

# Height
fit_p_h <- get_skewtr(Nplant, p_h_meanlog, p_h_sdlog) 

# Wing Pointdness
fit_b_wp <- get_skewtr(Nbird, b_wp_meanlog, b_wp_sdlog)

### Estimate interaction frequencies:

# We used fruit size and avian body mass to estimate interaction frequencies.

# Plants: We assume the relationship: y=1/x
p_v <- Plants$Fruit.length.mm * Plants$Fruit.D1.mm * Plants$Fruit.D2.mm

p_v_meanlog <- mean(log(na.omit(p_v)^(1/3)))
p_v_sdlog <- sd(log(na.omit(p_v)^(1/3)))
fit_p_v <- get_skewtr(Nplant, p_v_meanlog, p_v_sdlog) 

YES_pl_freq <- (1/fit_p_v)/mean((1/fit_p_v)) 

# Birds: We assume a negative relationship: y=(1/x)+b, where b is thee undercompensation parameter set to the 10% of the maximum value of 1/x
b_bm_meanlog <- mean(log(Birds$Bodymass.g^(1/3))) 
b_bm_sdlog <- sd(log(Birds$Bodymass.g^(1/3))) 
fit_b_bm <- (get_skewtr(Nbird, b_bm_meanlog, b_bm_sdlog))

v10 <- max(1/fit_b_bm)/10
YES_bird_freq <- ((1/fit_b_bm) + v10)/mean((1/fit_b_bm)+v10) 

# Set specialization parameter levels
networkreplicates <- 1
specsparlevel <- c(
#  rep(1.5, networkreplicates) # H2' = 0.1587484
 # , 
  rep(5, networkreplicates) # 0.308189
  ,
  rep(10, networkreplicates) # 0.4556191
  ,
  rep(20, networkreplicates) # 0.6532462
)

# Interaction probabilities according to trait matching (bill x fruit width and wing pointedness and plant height)
webs_p_2 <- vector("list", length(specsparlevel))
webs_p_t1 <- lapply(specsparlevel, makeweb_r, birdtraits = as.vector(decostand(fit_b_bw, "range")), planttraits = as.vector(decostand(fit_p_fw, "range")),
                    nicheshape = "skewed")
webs_p_t2 <- lapply(specsparlevel, makeweb_r, birdtraits = as.vector(decostand(fit_b_wp, "range")), planttraits = as.vector(decostand(fit_p_h, "range")),
                    nicheshape = "normal")
for(i in 1:length(specsparlevel)){
  t1t2 <- webs_p_t1[[i]]*webs_p_t2[[i]]
  t1t2 <- t1t2 / matrix(rowSums(t1t2), nrow = Nplant, ncol = Nbird)
  webs_p_2[[i]] <- t1t2
}

# Interaction networks built using interaction frequencies and interaction probabilites:
webs_2 <- lapply(webs_p_2, make_trueweb_r, birdfreq=YES_bird_freq, plantfreq=YES_pl_freq)
funct <- function(x){x/matrix(rowSums(x)/YES_pl_freq, nrow=Nplant, ncol=Nbird)}
webs_2 <- lapply(webs_2, funct)

obsperplant <- 100 

webs_largesamp_2 <- lapply(webs_2, sampleweb, obs=obsperplant, int_freq = YES_pl_freq, method='perplant')
webs_largesamp_2 <- lapply(webs_largesamp_2, t)
for(i in 1:length(webs_largesamp_2)){
  colnames(webs_largesamp_2[[i]]) <- c(1:Nbird)
  rownames(webs_largesamp_2[[i]]) <- c(1:Nplant)
  webs_largesamp_2[[i]] <- as.data.frame(webs_largesamp_2[[i]]) 
}

# Bird traits (for data frame)
b_bm_meanlog_dist <- mean(log(Birds$Bodymass.g)) 
b_bm_sdlog_dist <- sd(log(Birds$Bodymass.g)) 
fit_b_bm_dist <- (get_skewtr(Nbird, b_bm_meanlog_dist, b_bm_sdlog_dist))/1000

fit_b_bm_df <- as.data.frame(fit_b_bm_dist)
fit_b_bw_df <- as.data.frame(fit_b_bw)

birdtrait <- as.data.frame(cbind(fit_b_bm_df, as.character(c(1:Nbird)), fit_b_bw_df))
colnames(birdtrait)[2] <- "birdnumber"

# Plant trait (for data frame)
fit_p_fw_df <- as.data.frame(fit_p_fw)

planttrait <- as.data.frame(cbind(fit_p_fw_df, as.character(c(1:Nplant))))
colnames(planttrait)[2] <- "plantnumber"

# Creating data frames with simulated dispersal distances of all seed-dispersal events of each plant species in the source patch
simulatedframes_2 <- lapply(webs_largesamp_2, web2frame)
simulatedframes_2 <- lapply(simulatedframes_2, as.data.frame)
simulatedframes_2 <- lapply(simulatedframes_2, matchbirdplanttraits, t_obs = sum(rowSums(webs_largesamp_2[[1]])),
                            birdtraits = birdtrait, planttraits = planttrait) 
simulatedframes_2 <- lapply(simulatedframes_2, dispsimulation, t_obs = sum(rowSums(webs_largesamp_2[[1]])))

# since itneraction frequencies have decimals, the total number of interaction events is a bit below the sum of interaction frequencies
# Simulated between-patch distances
patches <- c(seq(from = 0, to = 600, by = 10))

# Calculate difference between seed-dispersal and between-patch distances to calculate reacheable patches
filter_distance_webs_2 <- lapply(simulatedframes_2, distancefilter, patches = patches,
                                 Nplant = Nplant, obs = obsperplant*YES_pl_freq)

# Setting number of iterations
Nrep <- 100 

# Interaction probabilities between seed dispersers and plants in the target patch (i.e., estimating fruit preferences)
rec_webs_p_2 <- vector("list", length = length(specsparlevel))

for(s in 1:length(specsparlevel)){
  rec_webs_p_t1 <- makeweb_c(specpar=specsparlevel[[s]], birdtraits = as.vector(decostand(fit_b_bw, "range")), 
                             planttraits = as.vector(decostand(fit_p_fw, "range")),
                             nicheshape = "skewed")
  
  rec_webs_p_t2 <- makeweb_c(specpar=specsparlevel[[s]], birdtraits = as.vector(decostand(fit_b_wp, "range")), 
                             planttraits = as.vector(decostand(fit_p_h, "range")),
                             nicheshape = "normal")
  
  t1t2 <- (rec_webs_p_t1*rec_webs_p_t2)
  t1t2 <- t1t2 / matrix(rowSums(t1t2), nrow = Nplant, ncol = Nbird)
  true_t1t2 <- make_trueweb_c(t1t2, YES_pl_freq)
  
  rec_webs_p_2[[s]] <- true_t1t2
  
}

for(i in 1:length(specsparlevel)) { 
  dimnames(rec_webs_p_2[[i]]) <- list(paste('p',1:Nplant,sep=''),paste('b',1:Nbird,sep='')) 
}

# Simulating plant assemblies at target patch

rec_webs_2 <- r_webs_rand(specsparlevel, Nplant, Nbird, Nrep, rec_webs_p_2)

# Calculate total plant abundance in the recipient communities

rec_plant_ab <- vector("list", length(specsparlevel))
names(rec_plant_ab) <- paste("Specialization_level", c(specsparlevel), sep = "_")

rec_plant_div <- vector("list", length(specsparlevel))
names(rec_plant_div) <- paste("Specialization_level", c(specsparlevel), sep = "_")

for(s in 1:length(specsparlevel)){
  
  plant_ab <- array(dim = c(length(rec_webs_2[[s]]), Nrep), dimnames = list(paste0("re_p_comm_", 1:(length(rec_webs_2[[s]]))),
                                                                            paste0("iteration_",1:Nrep)))
  plant_div <- array(dim = c(length(rec_webs_2[[s]]), Nrep), dimnames = list(paste0("re_p_comm_", 1:(length(rec_webs_2[[s]]))),
                                                                             paste0("iteration_",1:Nrep)))
  t1t2 <- rec_webs_p_2[[s]]
  
  for(r in 1:length(rec_webs_2[[s]])){
    
    for(n in 1:Nrep){
      
      web <- rec_webs_2[[s]][[r]][,,n]
      
      plants <- which(apply(t1t2, 1, function(row) row[1] %in% web))
      
      ab <- YES_pl_freq[plants]
      
      plant_ab[r,n] <- sum(ab)
      
      plant_div[r,n] <- diversity(ab)
    }
  }
  rec_plant_ab[[s]] <- plant_ab
  rec_plant_div[[s]] <- plant_div
}

filter_attrac_webs_2 <- attractivenessfilter(simulatedframes=simulatedframes_2, 
                                             source_webs_p = webs_p_2,
                                             recipient_webs=rec_webs_2, 
                                             specsparlevel=specsparlevel, 
                                             obs=obsperplant*YES_pl_freq, 
                                             Nplant=Nplant, 
                                             Nrep = Nrep,
                                             Ncons=Nbird, plant_freq = YES_pl_freq)

# Probability of seeds reaching the target patch per seed-dispersal event
seeds_reach_prob_webs_2 <- reach_prob_webs(recipient_webs = rec_webs_2, specsparlevel = specsparlevel, patches = patches, obs = obsperplant*YES_pl_freq, Nplant = Nplant,
                                           Nrep = Nrep, filter_distance_webs = filter_distance_webs_2, filter_attrac_webs = filter_attrac_webs_2)

# Total number of seeds reaching the target patch (sum over the probabilities per species)
n_seeds_reaching_webs_2 <- lapply(seeds_reach_prob_webs_2, seeds_reaching_webs,
                                  specsparlevel = specsparlevel, Nplant=Nplant, 
                                  Nrep=Nrep, patches=patches)

# Applying a threshold for the number of seeds that reach the recipient patch.
# *If the total number of seeds of a certain species reaching the target patch is below 1, we set the number of seeds to 0, indicating that no seeds from that plant species reached the target patch in that scenario of between-patch distance and plant richness of the target patch.

t_1 <- vector("list", length = length(n_seeds_reaching_webs_2))
for(r in 1:length(n_seeds_reaching_webs_2)){
  s_1 <- vector("list", length(specsparlevel))
  for(s in 1:length(specsparlevel)){
    s_1[[s]] <- ifelse(n_seeds_reaching_webs_2[[r]][[s]] < 1.0, 0, n_seeds_reaching_webs_2[[r]][[s]])
  }
  names(s_1) <- paste("Specialization_level", c(specsparlevel), sep = "_")
  t_1[[r]] <- s_1
}
names(t_1) <- paste("Plant_Rich", Nplant:0, sep = "_")

# Estimate plant diversity in the recipient patch 
  
# Estimate Number and Diversity
  
disp_dim_2 <- dispersal_dimensions(Nplant, specsparlevel, Nrep, patches, t_1, fit_p_fw, rec_webs_2, YES_pl_freq, obsperplant)

# Organize into one data set

map_rich_to_comm <- function(rich_recipient) {
  return((Nplant + 1) - rich_recipient)
} # mapping the Rich_recipient to the corresponding row in rec_plant_abundance

disp_dim_2$Iteration <- rep(rep(rep(1:Nrep, each = length(patches)), times = length(specsparlevel)), times = Nplant+1)

disp_dim_2$Ab_recipient <- mapply(function(rich, spec, iteration) {
  
  re_p_comm_index <- map_rich_to_comm(rich)
  
  spec_col <- paste0("Specialization_level_", spec)
  
  abundance <- rec_plant_ab[[spec_col]][re_p_comm_index, iteration]
  
  return(abundance)
  
}, disp_dim_2$Rich_recipient, disp_dim_2$Spec.par.level, disp_dim_2$Iteration)

disp_dim_2$Div_recipient <- mapply(function(rich, spec, iteration) {
  
  re_p_comm_index <- map_rich_to_comm(rich)
  
  spec_col <- paste0("Specialization_level_", spec)
  
  diversity <- rec_plant_div[[spec_col]][re_p_comm_index, iteration]
  
  return(diversity)
  
}, disp_dim_2$Rich_recipient, disp_dim_2$Spec.par.level, disp_dim_2$Iteration)

disp_dim_2$Div_recipient[disp_dim_2$Rich_recipient == 1] <- 0.0000000001

disp_dim_2$Iteration <- as.factor(disp_dim_2$Iteration)

# Save data

write.csv(disp_dim_2, "data_plotting.csv", row.names = FALSE)

#### Plots ####

library(ggplot2)
data <- read.csv("data_plotting.csv")
data$Spec.par.level <- as.factor(data$Spec.par.level)
data$Iteration <- as.factor(data$Iteration)

ggplot(data[data$Distance==100,]) +
  geom_line(aes(log(exp(Div_recipient)), y = Diversity, group = interaction(Spec.par.level, Iteration), color = Spec.par.level),linewidth = 0.1, alpha = 1) +
  geom_segment(aes(x = -Inf, y = 10, xend = +Inf, yend = 10), linetype = "dotted", linewidth = .9, colour = "red") +
  scale_color_manual(values = c("#FF595E", "#FFCA3A", "#8AC926", "#1982C4")) +  
  scale_fill_manual(values = c("#FF595E", "#FFCA3A", "#8AC926", "#1982C4")) + 
  ylim(c(0,50)) +
  # scale_y_continuous(limits = c(0, 15), breaks = seq(from=0, to=15, by = 5), labels = seq(from=0, to=15, by = 5)) + # 250m
  #scale_y_continuous(limits = c(0, 35), breaks = c(0, 10, 20, 30, 35), labels = c(0, 10, 20, 30, 50)) +
  scale_x_continuous(breaks = log(c(1, 2:6, 11, 8, 15, 25, 40)), labels = c(1, 2:6, 11, 8, 15, 25, 40)) +
  #labs(x = expression("Effective diversity of plants at the recipient patch \n          (logarithmic scale)"), y = 'Effective diversity of dispersed units', fill = "Specialization parameter") +
  guides(colour = "none") +
  theme_classic() +
  theme(axis.title = element_blank(), axis.text = element_text(size=12), legend.position = "none",
        axis.title.x = element_blank(), axis.title.y = element_blank())


##Abundance:

ggplot(data[data$Distance==250,], aes(log(exp(Div_recipient)), Abundance, group = interaction(Spec.par.level, Iteration), color = Spec.par.level)) +
  #ggtitle("10m") +
  geom_line(linewidth = 0.1, alpha = 1) +
  geom_segment(aes(x = -Inf, y = 500, xend = +Inf, yend = 500), linetype = "dotted", linewidth = .9, colour = "red") +
  scale_color_manual(values = c("#FF595E", "#FFCA3A", "#8AC926", "#1982C4")) +  
  scale_fill_manual(values = c("#FF595E", "#FFCA3A", "#8AC926", "#1982C4")) + 
  # scale_y_continuous(limits = c(0, 2500), breaks = c(0, 500, 1000, 1500, 2000, 2500), labels = c(0, 500, 1000, 1500, 2000, 2500)) + # 10 and 50m
  scale_y_continuous(limits = c(0, 800), breaks = c(0, seq(from=100, to=800, by=100)), labels = c(0,  seq(from=100, to=800, by=100))) + # 100 and 250m
  scale_x_continuous(breaks = log(c(1, 2:6, 11, 8, 15, 25, 40)), labels = c(1, 2:6, 11, 8, 15, 25, 40)) +
  labs(x = expression("Effective diversity of plants at the recipient patch \n          (logarithmic scale)"), y = 'Abundance of dispersed units', fill = "Specialization parameter") +
  guides(colour = "none") +
  theme_classic() +
  theme(axis.title = element_blank(), axis.text = element_text(size=12), legend.position = "none",
        axis.title.x = element_blank(), axis.title.y = element_blank(),  plot.margin = margin(t = 0.5, b = 1.5, l = 1.5, r = 0.5, unit = "lines"))

