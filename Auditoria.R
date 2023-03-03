# Author:   Javier Martinez Arribas       #
# Date:     9 - September - 2022          #
# Contact:  javimartinezarribas@gmail.com #
###########################################

###Libraries (maybe you have to add or delete some packages from the list)
load.libraries <- c('gamlss','tidyverse','ggpubr','skimr','readxl','fitdistrplus','zoib')
install.lib <- load.libraries[!load.libraries %in% installed.packages()]
for(libs in install.lib) install.packages(libs, dependences = TRUE)
sapply(load.libraries, require, character = TRUE)


### File for species
file1 = "NamibiaAmphibiaLogNormIgn.xls"  
file2 = "NamibiaAvesLogNormIgn.xls"   
file3 = "NamibiaInsectaLogNormIgn.xls"
file4 = "NamibiaLiliopsidaLogNormIgn.xls"
file5 = "NamibiaMagnoliopsidaLogNormIgn.xls"
file6 = "NamibiaMammaliaLogNormIgn.xls"
file7 = "NamibiaReptiliaLogNormIgn.xls"

df1 = data.frame(read_excel(file1))   #Amphibia
df2 = data.frame(read_excel(file2))   #Aves
df3 = data.frame(read_excel(file3))   #Insecta
df4 = data.frame(read_excel(file4))   #Liliopsida
df5 = data.frame(read_excel(file5))   #Magnoliopsida
df6 = data.frame(read_excel(file6))   #Mammalia
df7 = data.frame(read_excel(file7))   #Reptilia


###Distribution of Response Variable
# Ignorance distribution
g1_amph <- gamlss(df2$ignScore~1, family=BEINF1)
g1_amph # You can see the parameter values estimated
BEINF1() # You can see what parameters this distribution have in order to simulate it


#### 3. Checking Response Variable Distribution
# Simulation for Great Escarpments observations distribution to check the distribution
# For all Angola and Namibia Territory you will have to change them.
simu1 <- rBEINF1(
  n  = dim(df2)[1], # number of rows in df1
  mu = exp(0.7509)/(1+exp(0.7509)),
  sigma = exp(-0.4998)/(1+exp(-0.4998)),
  nu = exp(0.1004)
)

###Overlap both distributions real vs. simu1
ggplot() +
  geom_histogram(
    data = df1,
    aes(x = ignScore),
    alpha = 0.5, fill = "yellow") +
  geom_histogram(
    data = data.frame(simu1),
    aes(x = simu1),
    alpha = 0.3, fill = "firebrick") +
  labs(title = "Real Distribution vs Beta-One Inflated Simulation") +
  theme_bw()



################### MODELS #####################################

###Initial Model
#caveat: Variance not fitted because almost all target values are equal to one
mod_amph <- gamlss(
  formula = ignScore ~ pb(RD) + pb(UD) + pb(PAD) + pb(FC) + pb(POD),
  sigma.formula = ignScore ~ 1,
  nu.formula = ignScore ~ pb(RD) + pb(UD) + pb(PAD) + pb(FC) + pb(POD),
  family  = BEINF,
  data    = df1,
  trace   = FALSE
)
summary(mod_amph)

# Godness of fit
# If you got an advertency here, please increase ylim.all argument
wp(mod_amph, ylim.all = 1.2)


## Predictors selector
# Features in the table in which AIC is lower means that the model without that 
# variable get better
#drop1(mod_amph, parameter = "mu", parallel = "multicore", ncpus = 2)
mod_amph_mu <- stepGAIC(mod_amph, what="mu", parallel="snow",  ncpus=4)
#drop1(mod_amph, parameter = "nu", parallel = "multicore", ncpus = 2)
mod_amph_nu <- stepGAIC(mod_amph, what="nu", parallel="snow",  ncpus=4)

## Final Model
# Here your features will be another ones 
# Depending on drop1 results you will need to change the formula
mod_amph_final <- gamlss(
  formula = ignScore ~ pb(RD) + pb(UD) + pb(FC) + pb(PAD),
  nu.formula = ignScore ~ 1,
  family  = BEINF,
  data    = df1,
  trace   = FALSE,
  control = gamlss.control(mu.step=0.1)
)
summary(mod_amph_final)
AIC(mod_amph,mod_amph_final)
# Godness of fit
# If you got a warning here, please increase ylim.all argument
wpAmph <- wp(mod_amph_final, ylim.all = 1.5)


# To see the partial influence of all independent variables over the mean and  
# the shape of the target variable IGNORANCE
termplot(mod_amph_final, parameter = "mu", pages = 4, 
                      ask = FALSE, rug = TRUE, se=FALSE)
termplot(mod_amph_final, parameter = "nu", pages = 1, 
        ask = FALSE, rug = TRUE, se=FALSE, ylim='free',
        ylabs=rep("Logit of Ignorance Score",4),
        xlabs=c("Road Density","University Distance","Forest Cover","Protected Area Distance"),
        main="Partial Influence over Ignorance Scores")


## Aves

mod_aves <- gamlss(
  formula = ignScore ~ pb(RD) + pb(UD) + pb(PAD) + pb(FC) + pb(POD),
  sigma.formula = ignScore ~ 1,
  nu.formula = ignScore ~ pb(RD) + pb(UD) + pb(PAD) + pb(FC) + pb(POD),
  family  = BEINF,
  data    = df2,
  trace   = FALSE
)
summary(mod_aves)



# Godness of fit
# If you got an advertency here, please increase ylim.all argument
wp(mod_aves, ylim.all = 1.2)

#drop1(mod_aves, parameter = "mu", parallel = "multicore", ncpus = 2)
mod_aves_mu <- stepGAIC(mod_aves, what="mu", parallel="snow",  ncpus=4)
#drop1(mod_aves, parameter = "nu", parallel = "multicore", ncpus = 2)
mod_aves_nu <- stepGAIC(mod_aves, what="nu", parallel="snow",  ncpus=4)



## Final Model
# Here your features will be another ones 
# Depending on drop1 results you will need to change the formula
mod_aves_final <- gamlss(
  formula = ignScore ~ pb(RD) + pb(UD) + pb(PAD) + pb(FC),
  sigma.formula = ignScore ~ 1,
  nu.formula = ignScore ~ 1,
  family  = BEINF1,
  data    = df2,
  trace   = FALSE
)
summary(mod_aves_final)

# Godness of fit
# If you got a warning here, please increase ylim.all argument
wpAves <- wp(mod_aves_final, ylim.all = 1.5)

# To see the partial influence of all independent variables over the mean and  
# the shape of the target variable IGNORANCE
tpAves_mu <- termplot(mod_aves_final, parameter = "mu", pages = 6, 
                      ask = FALSE, rug = TRUE, se=FALSE)
tpAves_nu <- termplot(mod_aves_final, parameter = "nu", pages = 3, 
                      ask = FALSE, rug = TRUE, se=FALSE)

## Insecta

mod_insec <- gamlss(
  formula = ignScore ~ pb(RD) + pb(UD) + pb(PAD) + pb(FC) + pb(POD),
  sigma.formula = ignScore ~ 1,
  nu.formula = ignScore ~ pb(RD) + pb(UD) + pb(PAD) + pb(FC) + pb(POD),
  family  = BEINF,
  data    = df3,
  trace   = FALSE
)
summary(mod_insec)

# Godness of fit
# If you got an advertency here, please increase ylim.all argument
wp(mod_insec, ylim.all = 1.2)

#drop1(mod_insec, parameter = "mu", parallel = "multicore", ncpus = 2)
mod_insec_mu <- stepGAIC(mod_insec, what="mu", parallel="snow",  ncpus=4)

#drop1(mod_insec, parameter = "nu", parallel = "multicore", ncpus = 2)
mod_insec_nu <- stepGAIC(mod_insec, what="nu", parallel="snow",  ncpus=4)

## Final Model
# Here your features will be another ones 
# Depending on drop1 results you will need to change the formula
mod_insec_final <- gamlss(
  formula = ignScore ~ pb(RD) + pb(PAD) +  pb(POD),
  sigma.formula = ignScore ~ 1,
  nu.formula = ignScore ~ 1,
  family  = BEINF,
  data    = df3,
  trace   = FALSE
)
summary(mod_insec_final)

# Godness of fit
# If you got a warning here, please increase ylim.all argument
wpInsec <- wp(mod_insec_final, ylim.all = 1.5)


# To see the partial influence of all independent variables over the mean and  
# the shape of the target variable IGNORANCE
tpInsec_mu <- termplot(mod_insec_final, parameter = "mu", pages = 6, 
                      ask = FALSE, rug = TRUE, se=FALSE)
tpInsec_nu <- termplot(mod_insec_final, parameter = "nu", pages = 3, 
                      ask = FALSE, rug = TRUE, se=FALSE)


## Liliopsida

mod_lilio <- gamlss(
  formula = ignScore ~ pb(RD) + pb(UD) + pb(PAD) + pb(FC) + pb(POD),
  sigma.formula = ignScore ~ 1,
  nu.formula = ignScore ~ pb(RD) + pb(UD) + pb(PAD) + pb(FC) + pb(POD),
  family  = BEINF,
  data    = df4,
  trace   = FALSE
)
summary(mod_lilio)

# Godness of fit
# If you got an advertency here, please increase ylim.all argument
wp(mod_lilio, ylim.all = 1.2)

#drop1(mod_lilio, parameter = "mu", parallel = "multicore", ncpus = 2)
mod_lilio_mu <- stepGAIC(mod_lilio, what="mu", parallel="snow",  ncpus=4)
#drop1(mod_lilio, parameter = "nu", parallel = "multicore", ncpus = 2)
mod_lilio_nu <- stepGAIC(mod_lilio, what="nu", parallel="snow",  ncpus=4)


## Final Model
# Here your features will be another ones 
# Depending on drop1 results you will need to change the formula
mod_lilio_final <- gamlss(
  formula = ignScore ~ 1,
  sigma.formula = ignScore ~ 1,
  nu.formula = ignScore ~ 1,
  family  = BEINF,
  data    = df4,
  trace   = FALSE
)
summary(mod_lilio_final)

# Godness of fit
# If you got a warning here, please increase ylim.all argument
wpLilio <- wp(mod_lilio_final, ylim.all = 1.5)


# To see the partial influence of all independent variables over the mean and  
# the shape of the target variable IGNORANCE
tpLilio_mu <- termplot(mod_lilio_final, parameter = "mu", pages = 6, 
                      ask = FALSE, rug = TRUE, se=FALSE)
tpLilio_nu <- termplot(mod_lilio_final, parameter = "nu", pages = 3, 
                      ask = FALSE, rug = TRUE, se=FALSE)


## Magnoliopsida

mod_magno <- gamlss(
  formula = ignScore ~ pb(RD) + pb(UD) + pb(PAD) + pb(FC) + pb(POD),
  sigma.formula = ignScore ~ 1,
  nu.formula = ignScore ~ pb(RD) + pb(UD) + pb(PAD) + pb(FC) + pb(POD),
  family  = BEINF,
  data    = df5,
  trace   = FALSE
)
summary(mod_magno)

# Godness of fit
# If you got an advertency here, please increase ylim.all argument
wp(mod_magno, ylim.all = 1.2)

#drop1(mod_magno, parameter = "mu", parallel = "multicore", ncpus = 2)
#drop1(mod_magno, parameter = "nu", parallel = "multicore", ncpus = 2)
mod_magno_mu <- stepGAIC(mod_magno, what="mu", parallel="snow",  ncpus=4)
#drop1(mod_lilio, parameter = "nu", parallel = "multicore", ncpus = 2)
mod_magno_nu <- stepGAIC(mod_magno, what="nu", parallel="snow",  ncpus=4)


## Final Model
# Here your features will be another ones 
# Depending on drop1 results you will need to change the formula
mod_magno_final <- gamlss(
  formula = ignScore ~ pb(RD) + pb(UD) + pb(PAD),
  sigma.formula = ignScore ~ 1,
  nu.formula = ignScore ~ 1,
  family  = BEINF,
  data    = df5,
  trace   = FALSE
)
summary(mod_magno_final)

# Godness of fit
# If you got a warning here, please increase ylim.all argument
wpMagno <- wp(mod_magno_final, ylim.all = 1.5)


# To see the partial influence of all independent variables over the mean and  
# the shape of the target variable IGNORANCE
tpMagno_mu <- termplot(mod_magno_final, parameter = "mu", pages = 6, 
                      ask = FALSE, rug = TRUE, se=FALSE)
tpMagno_nu <- termplot(mod_magno_final, parameter = "nu", pages = 3, 
                      ask = FALSE, rug = TRUE, se=FALSE)


## Mammalia

mod_mamma <- gamlss(
  formula = ignScore ~ pb(RD) + pb(UD) + pb(PAD) + pb(FC) + pb(POD),
  sigma.formula = ignScore ~ 1,
  nu.formula = ignScore ~ pb(RD) + pb(UD) + pb(PAD) + pb(FC) + pb(POD),
  family  = BEINF,
  data    = df6,
  trace   = FALSE
)
summary(mod_mamma)

# Godness of fit
# If you got an advertency here, please increase ylim.all argument
wp(mod_mamma, ylim.all = 1.2)

#drop1(mod_mamma, parameter = "mu", parallel = "multicore", ncpus = 2)
#drop1(mod_mamma, parameter = "nu", parallel = "multicore", ncpus = 2)
mod_mamma_mu <- stepGAIC(mod_mamma, what="mu", parallel="snow",  ncpus=4)

mod_mamma_nu <- stepGAIC(mod_mamma, what="nu", parallel="snow",  ncpus=4)


## Final Model
# Here your features will be another ones 
# Depending on drop1 results you will need to change the formula
mod_mamma_final <- gamlss(
  formula = ignScore ~ pb(RD) + pb(UD) + pb(PAD) + pb(FC),
  sigma.formula = ignScore ~ 1,
  nu.formula = ignScore ~ 1,
  family  = BEINF,
  data    = df6,
  trace   = FALSE
)
summary(mod_mamma_final)

# Godness of fit
# If you got a warning here, please increase ylim.all argument
wpMamma <- wp(mod_mamma_final, ylim.all = 1.5)


# To see the partial influence of all independent variables over the mean and  
# the shape of the target variable IGNORANCE
tpMamma_mu <- termplot(mod_mamma_final, parameter = "mu", pages = 6, 
                      ask = FALSE, rug = TRUE, se=FALSE)
tpMamma_nu <- termplot(mod_mamma_final, parameter = "nu", pages = 3, 
                      ask = FALSE, rug = TRUE, se=FALSE)



## Reptilia

mod_repti <- gamlss(
  formula = ignScore ~ pb(RD) + pb(UD) + pb(PAD) + pb(FC) + pb(POD),
  sigma.formula = ignScore ~ 1,
  nu.formula = ignScore ~ pb(RD) + pb(UD) + pb(PAD) + pb(FC) + pb(POD),
  family  = BEINF,
  data    = df7,
  trace   = FALSE
)
summary(mod_repti)

# Godness of fit
# If you got an advertency here, please increase ylim.all argument
wp(mod_repti, ylim.all = 1.2)

#drop1(mod_repti, parameter = "mu", parallel = "multicore", ncpus = 2)
#drop1(mod_repti, parameter = "nu", parallel = "multicore", ncpus = 2)
mod_repti_mu <- stepGAIC(mod_repti, what="mu", parallel="snow",  ncpus=4)

mod_repti_nu <- stepGAIC(mod_repti, what="nu", parallel="snow",  ncpus=4)

## Final Model
# Here your features will be another ones 
# Depending on drop1 results you will need to change the formula
mod_repti_final <- gamlss(
  formula = ignScore ~ pb(RD) + pb(PAD) + pb(POD),
  sigma.formula = ignScore ~ 1,
  nu.formula = ignScore ~ 1,
  family  = BEINF,
  data    = df7,
  trace   = FALSE
)
summary(mod_insec_final)

# Godness of fit
# If you got a warning here, please increase ylim.all argument
wpRepti <- wp(mod_repti_final, ylim.all = 1.5)


# To see the partial influence of all independent variables over the mean and  
# the shape of the target variable IGNORANCE
tpRepti_mu <- termplot(mod_repti_final, parameter = "mu", pages = 6, 
                      ask = FALSE, rug = TRUE, se=FALSE)
tpRepti_nu <- termplot(mod_repti_final, parameter = "nu", pages = 3, 
                      ask = FALSE, rug = TRUE, se=FALSE)
