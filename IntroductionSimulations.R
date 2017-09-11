#########################################
#########################################
### INTRODUCTION - SIMULATING EFFECTS ###
 ### WITH AND WITHOUT HETEROGENEITY ###
        ### JUNE 10 2017 ###
#########################################
#########################################
packages <- c('tidyverse', 'effects', 'clubSandwich', 'plm')

### UNCOMMENT LINE BELOW IF YOU HAVE *NOT* INSTALLED PACKAGES ALREADY ###
# install.packages(packages, repos = "http://cran.r-project.org")

### CALLING PACKAGES ###
for (i in packages){
  library(i,character.only = TRUE)
}


set.seed(870523)
sims <- 50000

########################################
### CREATING DATA-GENERATING PROCESS ###
   ### WITH HETEROGENEOUS EFFECTS ###
########################################

N <- 300000 # `POPULATION' SAMPLE

# 12% HIP (WINNING) T1=1, T2=1 -- typeA
# 7% LIP (WINNING) T1=1, T2=1 -- typeB
# 23% HINP T1=0, T2=1  -- typeC
# 31% LINP T1=0 T2=0 -- typeD
# 11% LIP (LOSING) T1=0 T2=0 -- typeE
# 16% HIP (LOSING) T1=0 T2=0 -- typeF

voter_types <- sample(x = c("typeA", "typeB", "typeC", "typeD", "typeE", "typeF"), 
                size = N, 
                replace = TRUE, 
                prob = c(.12, .07, .23, .31, .11, .16))

voter_types <- data.frame(voter_types)
voter_types$id <- row.names(voter_types)

### CREATE RANDOM ATTITUDE B/W 0 AND 1 ###
voter_types$attitude_t1 <- rnorm(N, .5, 2)
voter_types$attitude_t1 <- (voter_types$attitude_t1 - min(voter_types$attitude_t1))/(max(voter_types$attitude_t1) - min(voter_types$attitude_t1))

### CREATE HYPOTHETICAL INTERACTION COEFFICIENTS ###
alpha <- 0.3
attitude_coef <- 1.1
interactioncoef_all <- 1.1
interactioncoefA <- 1.01
interactioncoefB <- 1.02
interactioncoefC <- 1.5
interactioncoefD <- 1
interactioncoefE <- 1.02
interactioncoefF <- 1.01
pidwin <- 1.00
pidlose <- -1.2
nonpid <- 0

voter_types$interaction_het[voter_types$voter_types=="typeA"] <- interactioncoefA
voter_types$interaction_het[voter_types$voter_types=="typeB"] <- interactioncoefB
voter_types$interaction_het[voter_types$voter_types=="typeC"] <- interactioncoefC
voter_types$interaction_het[voter_types$voter_types=="typeD"] <- interactioncoefD
voter_types$interaction_het[voter_types$voter_types=="typeE"] <- interactioncoefE
voter_types$interaction_het[voter_types$voter_types=="typeF"] <- interactioncoefF

voter_types$pid[voter_types$voter_types=="typeA"] <- pidwin
voter_types$pid[voter_types$voter_types=="typeB"] <- pidwin
voter_types$pid[voter_types$voter_types=="typeC"] <- nonpid
voter_types$pid[voter_types$voter_types=="typeD"] <- nonpid
voter_types$pid[voter_types$voter_types=="typeE"] <- pidlose
voter_types$pid[voter_types$voter_types=="typeF"] <- pidlose

### CREATING PRED PROB OF VOTING FOR PARTY ###
### T1 ###
voter_types$predprob_all_t1 <- alpha + attitude_coef*voter_types$attitude_t1 + voter_types$pid + rnorm(N,0,2)
### T2 ###
voter_types$predprob_all_t2 <- alpha + attitude_coef*voter_types$interaction_het*voter_types$attitude_t1 + 
              voter_types$pid + rnorm(N,0,2)

### CREATING LONG DATA ###
estimates2 <- reshape(voter_types, 
                      varying = c("predprob_all_t1", "predprob_all_t2"), 
                      v.names = "predprob",
                      timevar = "time", 
                      times = c(0,1), 
                      direction = "long")

estimates2$time <- as.factor(estimates2$time)

############################################
### BUILDING MODEL WITH NO HETEROGENEITY ###
############################################
### PRIMING ###
### ALL OBS ###
fit1 <- lm(predprob~attitude_t1*time + pid, data=estimates2)
### CLUSTERED SE ###
fit1_cr1se <- coef_test(fit1, vcov = "CR1", cluster = estimates2$id, test = "z")
### PLOTTING INTERACTION COEF ###
plot(effect("attitude_t1:time", fit1, se=TRUE), multiline=TRUE)

####################################
    ### ALL VOTERS TOGETHER ###
    ### SIMULATING N NUMBER ###
     ### OF DRAWS TO SHOW ###
  ### IMPACT OF RANDOM SAMPLING ###
### ACROSS DIFFERENT VOTER TYPES ###
####################################

### GETTING UNIQUE IDs FROM LONG DATASET FOR SAMPLING ###
subject_ids = unique(estimates2$id)

### CREATE LIST TO STORE ESTIMATES ###
outcomes_all <- list()

### TAKE RANDOM SAMPLE OF 1,100 AND RUN REGRESSION `sims': # OF TIMES ###
for (i in 1:sims){
  sample_subject_ids = sample(subject_ids, 1100)
  types = subset(estimates2, id %in% sample_subject_ids)
  fit <- lm(predprob~attitude_t1*time + pid, data=types)
  my_estimate <- coef_test(fit, vcov = "CR1", 
                           cluster = types$id, test = "z")[5,]
  names(my_estimate) <- c("est", "se", "p")
  outcomes_all[[i]] <- my_estimate
}

### ANALYSING DISTRIBUTIONS OF SIMULATED RESULTS ###

### BIND ESTIMATES LIST TOGETHER ###
homogeneous_sims <- data.frame(do.call(rbind, outcomes_all))

### CREATE MARKER FOR SIG INTERACTION COEFFICIENTS ###
homogeneous_sims$sig <- as.factor(homogeneous_sims$p <= .05)

homogeneous_sims %>%
  summarise(effect = mean(est),
            sd_estimate = sd(est),
            wrong_sign = mean(est < 0),
            success = mean(p <= 0.05 & est > 0))

### PLOT RESULTS ###
ggplot(homogeneous_sims, aes(x=est, color = sig, fill = sig)) + 
  geom_histogram(bins=1000)


################################
### HET. APPROACH TO PRIMING ###
################################

### HINP ###

### SUBSET DATA FOR JUST HINP ###

hinp <- subset(estimates2, voter_types=="typeC")

### GETTING UNIQUE IDs FROM LONG DATASET FOR SAMPLING ###
subject_ids_hinp = unique(hinp$id)

### CREATE LIST TO STORE ESTIMATES ###
outcomes_hinp <- list()

### TAKE RANDOM SAMPLE OF 1,100 AND RUN REGRESSION `sims': # OF TIMES ###
for (i in 1:sims){
  sample_subject_ids = sample(subject_ids_hinp, 1100)
  types_hinp = subset(hinp, id %in% sample_subject_ids)
  fit <- lm(predprob~attitude_t1*time + pid, data=types_hinp)
  my_estimate <- coef_test(fit, vcov = "CR1", 
                           cluster = types_hinp$id, test = "z")[5,]
  names(my_estimate) <- c("est", "se", "p")
  outcomes_hinp[[i]] <- my_estimate
}

### BIND ESTIMATES LIST TOGETHER ###
hinp_sims <- data.frame(do.call(rbind, outcomes_hinp))

### CREATE MARKER FOR SIG INTERACTION COEFFICIENTS ###
hinp_sims$sig <- as.factor(hinp_sims$p <= .05)

hinp_sims %>%
  summarise(effect = mean(est),
            sd_estimate = sd(est),
            wrong_sign = mean(est < 0),
            success = mean(p <= 0.05 & est > 0))

### PLOT RESULTS ###
ggplot(hinp_sims, aes(x=est, color = sig, fill = sig)) + 
  geom_histogram(bins=1000)

