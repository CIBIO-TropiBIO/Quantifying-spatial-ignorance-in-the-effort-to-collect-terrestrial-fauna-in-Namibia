# Author:   Javier Martinez Arribas       #
# Date:     22 - November - 2021          #
# Contact:  javimartinezarribas@gmail.com #
###########################################

###Libraries (maybe you have to add or delete some packages from the list)
load.libraries <- c('gamlss','tidyverse','ggpubr','skimr','readxl','fitdistrplus','zoib')
install.lib <- load.libraries[!load.libraries %in% installed.packages()]
for(libs in install.lib) install.packages(libs, dependences = TRUE)
sapply(load.libraries, require, character = TRUE)


#######################################
########### ANGOLA ####################
#######################################

# Load features
path_Ang_Nam = "........." #Insert the path of "Data_Ang_Nam_All.xlsx" in your PC
df_all = data.frame(read_excel(path_Ang_Nam))

################# AVES 10 Obs #############################

# Load target
path = "......."# Insert the path of "Aves_half_ign10Ignorance.rda" for Angola in your PC
load(path)
## There are duplicates due to the cells in the border between Angola and Namibia
df_class_ign = df_class_ign[!duplicated(df_class_ign$FID),]
colnames(df_class_ign)[1] = "IGNORANCE"

df_final = merge(x = df_all, y = df_class_ign, by = "FID", all.x = TRUE)

# Ignorance distribution
g1_ang <- gamlss(df_final$IGNORANCE~1, family=BEINF1)
g1_ang # You can see the parameter values estimated
BEINF1() # You can see what parameters this distribution have in order to simulate it

# Simulation for Great Escarpments observations distribution to check the distribution
# For all Angola and Namibia Territory you will have to change them.
simulation <- rBEINF1(
  n  = 2778,
  mu = exp(-0.7059)/(1+exp(-0.7059)),
  sigma = exp(0.3398)/(1+exp(0.3398)),
  nu = exp(-0.06626)
)

ggplot() +
  geom_histogram(
    data = df_final,
    aes(x = IGNORANCE),
    alpha = 0.5, fill = "gray50") +
  geom_histogram(
    data = data.frame(simulation),
    aes(x = simulation),
    alpha = 0.3, fill = "firebrick") +
  labs(title = "Real Distribution vs Beta-One Inflated Simulation") +
  theme_bw()


################### MODELS #####################################

###Initial Model
#caveat: Variance not fitted because almost all target values are equal to one
model_gamlss_ang_aves10 <- gamlss(
  formula = IGNORANCE ~ pb(RD) + pb(UD) + pb(PAD) + 
    pb(FC) + pb(POD) + pb(STP),
  nu.formula = IGNORANCE ~ pb(RD) + pb(UD) + pb(PAD) + 
    pb(FC) + pb(POD) + pb(STP),
  family  = BEINF1,
  data    = df_final,
  trace   = FALSE
)
summary(model_gamlss_ang_aves10)

# Godness of fit
# If you got an advertency here, please increase ylim.all argument
wp(model_gamlss_ang_aves10, ylim.all = 1.2)

## Predictors selector
# Features in the table in which AIC is lower means that the model without that 
# variable get better
drop1(model_gamlss_ang_aves10, parameter = "mu", parallel = "multicore", ncpus = 2)
drop1(model_gamlss_ang_aves10, parameter = "nu", parallel = "multicore", ncpus = 2)


## Final Model
# Here your features will be another ones 
# Depending on drop1 results you will need to change the formula
model_gamlss_ang_aves10_final <- gamlss(
  formula = IGNORANCE ~ pb(UD) + pb(PAD) + pb(FC),
  nu.formula = IGNORANCE ~ pb(RD) + pb(UD) + pb(FC) + pb(STP),
  family  = BEINF1,
  data    = df_final,
  trace   = FALSE
)
summary(model_gamlss_ang_aves10_final)

# Godness of fit
# If you got a warning here, please increase ylim.all argument
wp(model_gamlss_ang_aves10_final, ylim.all = 1.5)


# To see the partial influence of all independent variables over the mean and  
# the shape of the target variable IGNORANCE
term.plot(model_gamlss_ang_aves10_final, parameter = "mu", pages = 3, 
          ask = FALSE, rug = TRUE, se=FALSE)
term.plot(model_gamlss_ang_aves10_final, parameter = "nu", pages = 4, 
          ask = FALSE, rug = TRUE, se=FALSE)


################# MAMMALS 10 Obs########################## 

# Load target
path = "......."# Insert the path of "Mammalia_half_ign10Ignorance.rda" for Angola in your PC
load(path)
## There are duplicates due to the cells in the border between Angola and Namibia
df_class_ign = df_class_ign[!duplicated(df_class_ign$FID),]
colnames(df_class_ign)[1] = "IGNORANCE"

df_final = merge(x = df_all, y = df_class_ign, by = "FID", all.x = TRUE)

# Ignorance distribution
g2_ang <- gamlss(df_final$IGNORANCE~1, family=BEINF1)
g2_ang
BEINF1()

# Simulation for Great Escarpments observations distribution to check the distribution
# For all Angola and Namibia Territory you will have to change them.
simulation <- rBEINF1(
  n  = 2775,
  mu = exp(1.054)/(1+exp(1.054)),
  sigma = exp(-0.3015)/(1+exp(-0.3015)),
  nu = exp(1.681)
)

ggplot() +
  geom_histogram(
    data = df_final,
    aes(x = IGNORANCE),
    alpha = 0.5, fill = "gray50") +
  geom_histogram(
    data = data.frame(simulation),
    aes(x = simulation),
    alpha = 0.3, fill = "firebrick") +
  labs(title = "Real Distribution vs Beta-One Inflated Simulation") +
  theme_bw()


################### MODELS #####################################

###Initial Model
#caveat: Variance not fitted because almost all target values are equal to one
model_gamlss_ang_mammals10 <- gamlss(
  formula = IGNORANCE ~ pb(RD) + pb(UD) + pb(PAD) + 
    pb(FC) + pb(POD) + pb(STP),
  nu.formula = IGNORANCE ~ pb(RD) + pb(UD) + pb(PAD) + 
    pb(FC) + pb(POD) + pb(STP),
  family  = BEINF1,
  data    = df_final,
  trace   = FALSE
)
summary(model_gamlss_ang_mammals10)

# Godness of fit
# If you got a warning here, please increase ylim.all argument
wp(model_gamlss_ang_mammals10, ylim.all = 1.2)


###Predictors selector
# Features in the table in which AIC is lower means that the model without that 
# variable get better
drop1(model_gamlss_ang_mammals10, parameter = "mu", parallel = "multicore", ncpus = 2)
drop1(model_gamlss_ang_mammals10, parameter = "nu", parallel = "multicore", ncpus = 2)


###Final Model
# Depending on drop1 results you will need to change the formula
model_gamlss_ang_mammals10_final <- gamlss(
  formula = IGNORANCE ~ pb(RD) + pb(UD) + pb(PAD) + 
    pb(FC) + pb(POD) + pb(STP),
  nu.formula = IGNORANCE ~ pb(RD) + pb(UD) + pb(PAD) + 
    pb(FC) + pb(POD) + pb(STP),
  family  = BEINF1,
  data    = df_final,
  trace   = FALSE
)
summary(model_gamlss_ang_mammals10_final)

# Godness of fit
# If you got a warning here, please increase ylim.all argument
wp(model_gamlss_ang_mammals10_final, ylim.all = 1.2)

# To see the partial influence of all independent variables over the mean and  
# the shape of the target variable IGNORANCE
term.plot(model_gamlss_ang_mammals10_final, parameter = "mu", pages = 1, 
          ask = FALSE, rug = TRUE, se=FALSE)
term.plot(model_gamlss_ang_mammals10_final, parameter = "nu", pages = 1, 
          ask = FALSE, rug = TRUE, se=FALSE)


################# MAGNOLIOPSIDA 10 Obs########################## 

# Load target
path = "......."# Insert the path of "Magnoliopsida_half_ign10Ignorance.rda" for Angola in your PC
load(path)
## There are duplicates due to the cells in the border between Angola and Namibia
df_class_ign = df_class_ign[!duplicated(df_class_ign$FID),]
colnames(df_class_ign)[1] = "IGNORANCE"

df_final = merge(x = df_all, y = df_class_ign, by = "FID", all.x = TRUE)

# Ignorance distribution
g3_ang <- gamlss(df_final$IGNORANCE~1, family=BEINF1)
g3_ang
BEINF1()

# Simulation for Great Escarpments observations distribution to check the distribution
# For all Angola and Namibia Territory you will have to change them.
simulation <- rBEINF1(
  n  = 2775,
  mu = exp(0.9491)/(1+exp(0.9491)),
  sigma = exp(-0.2837)/(1+exp(-0.2837)),
  nu = exp(1.121)
)

ggplot() +
  geom_histogram(
    data = df_final,
    aes(x = IGNORANCE),
    alpha = 0.5, fill = "gray50") +
  geom_histogram(
    data = data.frame(simulation),
    aes(x = simulation),
    alpha = 0.3, fill = "firebrick") +
  labs(title = "Real Distribution vs Beta-One Inflated Simulation") +
  theme_bw()


################### MODELS #####################################

###Initial Model
#caveat: Variance not fitted because almost all target values are equal to one
model_gamlss_ang_magnoliopsida10 <- gamlss(
  formula = IGNORANCE ~ pb(RD) + pb(UD) + pb(PAD) + 
    pb(FC) + pb(POD) + pb(STP),
  nu.formula = IGNORANCE ~ pb(RD) + pb(UD) + pb(PAD) + 
    pb(FC) + pb(POD) + pb(STP),
  family  = BEINF1,
  data    = df_final,
  trace   = FALSE
)
summary(model_gamlss_ang_magnoliopsida10)

# Godness of fit
# If you have got a warning here, please increase ylim.all argument
wp(model_gamlss_ang_magnoliopsida10, ylim.all = 1.2)


###Predictors selector
# Features in the table in which AIC is lower means that the model without that 
# variable get better
drop1(model_gamlss_ang_magnoliopsida10, parameter = "mu", parallel = "multicore", ncpus = 2)
drop1(model_gamlss_ang_magnoliopsida10, parameter = "nu", parallel = "multicore", ncpus = 2)


###Final Model
# Depending on drop1 results you will need to change the formula
model_gamlss_ang_magnoliopsida10_final <- gamlss(
  formula = IGNORANCE ~ pb(RD) + pb(UD) + pb(PAD) + 
    pb(FC) + pb(POD) + pb(STP),
  nu.formula = IGNORANCE ~ pb(RD) + pb(UD) + pb(PAD) + 
    pb(FC) + pb(POD) + pb(STP),
  family  = BEINF1,
  data    = df_final,
  trace   = FALSE
)
summary(model_gamlss_ang_magnoliopsida10_final)

# Godness of fit
# If you got a warning here, please increase ylim.all argument
wp(model_gamlss_ang_magnoliopsida10_final, ylim.all = 1.2)

# To see the partial influence of all independent variables over the mean and  
# the shape of the target variable IGNORANCE
term.plot(model_gamlss_ang_magnoliopsida10_final, parameter = "mu", pages = 1, 
          ask = FALSE, rug = TRUE, se=FALSE)
term.plot(model_gamlss_ang_magnoliopsida10_final, parameter = "nu", pages = 1, 
          ask = FALSE, rug = TRUE, se=FALSE)


#######################################
########### NAMIBIA ###################
#######################################


################# AVES 10 Obs #############################

# Load target
path = "......."# Insert the path of "Aves_half_ign10Ignorance.rda" for Namibia in your PC
load(path)
## There are duplicates due to the cells in the border between Angola and Namibia
df_class_ign = df_class_ign[!duplicated(df_class_ign$FID),]
colnames(df_class_ign)[1] = "IGNORANCE"

df_final = merge(x = df_all, y = df_class_ign, by = "FID", all.x = TRUE)

# Ignorance distribution
g1_nam <- gamlss(df_final$IGNORANCE~1, family=BEINF1)
g1_nam # You can see the parameter values estimated
BEINF1() # You can see what parameters needs this distribution

# Simulation for Great Escarpments observations distribution to check the distribution
# For all Angola and Namibia Territory you will have to change them.
simulation <- rBEINF1(
  n  = 2778,
  mu = exp(-0.7059)/(1+exp(-0.7059)),
  sigma = exp(0.3398)/(1+exp(0.3398)),
  nu = exp(-0.06626)
)

ggplot() +
  geom_histogram(
    data = df_final,
    aes(x = IGNORANCE),
    alpha = 0.5, fill = "gray50") +
  geom_histogram(
    data = data.frame(simulation),
    aes(x = simulation),
    alpha = 0.3, fill = "firebrick") +
  labs(title = "Real Distribution vs Beta-One Inflated Simulation") +
  theme_bw()


################### MODELS #####################################

###Initial Model
#caveat: Variance not fitted because almost all target values are equal to one
model_gamlss_nam_aves10 <- gamlss(
  formula = IGNORANCE ~ pb(RD) + pb(UD) + pb(PAD) + 
    pb(FC) + pb(POD) + pb(STP),
  nu.formula = IGNORANCE ~ pb(RD) + pb(UD) + pb(PAD) + 
    pb(FC) + pb(POD) + pb(STP),
  family  = BEINF1,
  data    = df_final,
  trace   = FALSE
)
summary(model_gamlss_nam_aves10)

# Godness of fit
# If you have got a warning here, please increase ylim.all argument
wp(model_gamlss_nam_aves10, ylim.all = 1.2)

###Predictors selector
# Features in the table in which AIC is lower means that the model get better
# without that variable
drop1(model_gamlss_nam_aves10, parameter = "mu", parallel = "multicore", ncpus = 2)
drop1(model_gamlss_nam_aves10, parameter = "nu", parallel = "multicore", ncpus = 2)

###Final Model
# Depending on drop1 results you will need to change the formula
model_gamlss_nam_aves10_final <- gamlss(
  formula = IGNORANCE ~ pb(UD) + pb(PAD) + pb(FC),
  nu.formula = IGNORANCE ~ pb(RD) + pb(UD) + pb(FC) + pb(STP),
  family  = BEINF1,
  data    = df_final,
  trace   = FALSE
)
summary(model_gamlss_nam_aves10_final)

# Godness of fit
# If you got an advertency here, please increase ylim.all argument
wp(model_gamlss_nam_aves10_final, ylim.all = 1.5)

# To see the partial influence of all independent variables over the mean and  
# the shape of the target variable IGNORANCE
term.plot(model_gamlss_nam_aves10_final, parameter = "mu", pages = 3, 
          ask = FALSE, rug = TRUE, se=FALSE)
term.plot(model_gamlss_nam_aves10_final, parameter = "nu", pages = 4, 
          ask = FALSE, rug = TRUE, se=FALSE)


################# MAMMALS 10 Obs########################## 

# Load target
path = "......."# Insert the path of "Mammalia_half_ign10Ignorance.rda" for Namibia in your PC
load(path)
## There are duplicates due to the cells in the border between Angola and Namibia
df_class_ign = df_class_ign[!duplicated(df_class_ign$FID),]
colnames(df_class_ign)[1] = "IGNORANCE"

df_final = merge(x = df_all, y = df_class_ign, by = "FID", all.x = TRUE)

# Ignorance distribution
g2_nam <- gamlss(df_final$IGNORANCE~1, family=BEINF1)
g2_nam
BEINF1()

# Simulation for Great Escarpments observations distribution to check the distribution
# For all Angola and Namibia Territory you will have to change them.
simulation <- rBEINF1(
  n  = 2775,
  mu = exp(1.054)/(1+exp(1.054)),
  sigma = exp(-0.3015)/(1+exp(-0.3015)),
  nu = exp(1.681)
)

ggplot() +
  geom_histogram(
    data = df_final,
    aes(x = IGNORANCE),
    alpha = 0.5, fill = "gray50") +
  geom_histogram(
    data = data.frame(simulation),
    aes(x = simulation),
    alpha = 0.3, fill = "firebrick") +
  labs(title = "Real Distribution vs Beta-One Inflated Simulation") +
  theme_bw()


################### MODELS #####################################

###Initial Model
#caveat: Variance not fitted because almost all target values are equal to one
model_gamlss_nam_mammals10 <- gamlss(
  formula = IGNORANCE ~ pb(RD) + pb(UD) + pb(PAD) + 
    pb(FC) + pb(POD) + pb(STP),
  nu.formula = IGNORANCE ~ pb(RD) + pb(UD) + pb(PAD) + 
    pb(FC) + pb(POD) + pb(STP),
  family  = BEINF1,
  data    = df_final,
  trace   = FALSE
)
summary(model_gamlss_nam_mammals10)

# Godness of fit
# If you got an advertency here, please increase ylim.all argument
wp(model_gamlss_nam_mammals10, ylim.all = 1.2)


###Predictors selector
# Features in the table in which AIC is lower means that the model get better
# without that variable
drop1(model_gamlss_nam_mammals10, parameter = "mu", parallel = "multicore", ncpus = 2)
drop1(model_gamlss_nam_mammals10, parameter = "nu", parallel = "multicore", ncpus = 2)


###Final Model
# Depending on drop1 results you will need to change the formula
model_gamlss_nam_mammals10_final <- gamlss(
  formula = IGNORANCE ~ pb(RD) + pb(UD) + pb(PAD) + 
    pb(FC) + pb(POD) + pb(STP),
  nu.formula = IGNORANCE ~ pb(RD) + pb(UD) + pb(PAD) + 
    pb(FC) + pb(POD) + pb(STP),
  family  = BEINF1,
  data    = df_final,
  trace   = FALSE
)
summary(model_gamlss_nam_mammals10_final)

# Godness of fit
# If you got an advertency here, please increase ylim.all argument
wp(model_gamlss_nam_mammals10_final, ylim.all = 1.2)

# To see the partial influence of all independent variables over the mean and  
# the shape of the target variable IGNORANCE
term.plot(model_gamlss_nam_mammals10_final, parameter = "mu", pages = 1, 
          ask = FALSE, rug = TRUE, se=FALSE)
term.plot(model_gamlss_nam_mammals10_final, parameter = "nu", pages = 1, 
          ask = FALSE, rug = TRUE, se=FALSE)



################# MAGNOLIOPSIDA 10 Obs########################## 

# Load target
path = "......."# Insert the path of "Magnoliopsida_half_ign10Ignorance.rda" for Namibia in your PC
load(path)
## There are duplicates due to the cells in the border between Angola and Namibia
df_class_ign = df_class_ign[!duplicated(df_class_ign$FID),]
colnames(df_class_ign)[1] = "IGNORANCE"

df_final = merge(x = df_all, y = df_class_ign, by = "FID", all.x = TRUE)

# Ignorance distribution
g3_nam <- gamlss(df_final$IGNORANCE~1, family=BEINF1)
g3_nam
BEINF1()

# Simulation for Great Escarpments observations distribution to check the distribution
simulation <- rBEINF1(
  n  = 2775,
  mu = exp(0.9491)/(1+exp(0.9491)),
  sigma = exp(-0.2837)/(1+exp(-0.2837)),
  nu = exp(1.121)
)

ggplot() +
  geom_histogram(
    data = df_final,
    aes(x = IGNORANCE),
    alpha = 0.5, fill = "gray50") +
  geom_histogram(
    data = data.frame(simulation),
    aes(x = simulation),
    alpha = 0.3, fill = "firebrick") +
  labs(title = "Real Distribution vs Beta-One Inflated Simulation") +
  theme_bw()

################### MODELS #####################################

###Initial Model
#caveat: Variance not fitted because almost all target values are equal to one
model_gamlss_nam_magnoliopsida10 <- gamlss(
  formula = IGNORANCE ~ pb(RD) + pb(UD) + pb(PAD) + 
    pb(FC) + pb(POD) + pb(STP),
  nu.formula = IGNORANCE ~ pb(RD) + pb(UD) + pb(PAD) + 
    pb(FC) + pb(POD) + pb(STP),
  family  = BEINF1,
  data    = df_final,
  trace   = FALSE
)
summary(model_gamlss_nam_magnoliopsida10)

# Godness of fit
# If you got an advertency here, please increase ylim.all argument
wp(model_gamlss_nam_magnoliopsida10, ylim.all = 1.2)


###Predictors selector
# Features in the table in which AIC is lower means that the model get better
# without that variable
drop1(model_gamlss_nam_magnoliopsida10, parameter = "mu", parallel = "multicore", ncpus = 2)
drop1(model_gamlss_nam_magnoliopsida10, parameter = "nu", parallel = "multicore", ncpus = 2)


###Final Model
# Here your features will be other ones. 
# Depending on drop1 results you will need to change the formula
model_gamlss_nam_magnoliopsida10_final <- gamlss(
  formula = IGNORANCE ~ pb(RD) + pb(UD) + pb(PAD) + 
    pb(FC) + pb(POD) + pb(STP),
  nu.formula = IGNORANCE ~ pb(RD) + pb(UD) + pb(PAD) + 
    pb(FC) + pb(POD) + pb(STP),
  family  = BEINF1,
  data    = df_final,
  trace   = FALSE
)
summary(model_gamlss_nam_magnoliopsida10_final)

# Godness of fit
# If you got an advertency here, please increase ylim.all argument
wp(model_gamlss_nam_magnoliopsida10_final, ylim.all = 1.2)

# To see the partial influence of all independent variables over the mean and  
# the shape of the target variable IGNORANCE
term.plot(model_gamlss_nam_magnoliopsida10_final, parameter = "mu", pages = 1, 
          ask = FALSE, rug = TRUE, se=FALSE)
term.plot(model_gamlss_nam_magnoliopsida10_final, parameter = "nu", pages = 1, 
          ask = FALSE, rug = TRUE, se=FALSE)

