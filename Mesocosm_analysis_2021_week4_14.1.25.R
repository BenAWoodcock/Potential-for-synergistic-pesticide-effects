# NEC06567 NERC-RP Synergistic mixture effects and sensitisation
# Analysis of synergistic effects between cypermethrin and the fungicides
# prochloraz and azoxystrobin within complex arable communites

library(lme4) # Mixed effects models
library(pbkrtest)  #  Kenward-Rogers degrees freedom for significance testing of GLMM
library(emmeans)   #   estiamted model means

library(ggplot2) #  graphical outputs
library(DHARMa) #  goodness of fit tests
library(car)   #  model simplification
library(glmmTMB) # for neg bin models
library(dplyr)

# code for model simpification in pbkrtest
KRSumFun <- function(object, objectDrop, ...) {
  krnames <- c("ndf","ddf","Fstat","p.value","F.scaling")
  r <- if (missing(objectDrop)) {
    setNames(rep(NA,length(krnames)),krnames)
  } else {
    krtest <- KRmodcomp(object,objectDrop)
    unlist(krtest$stats[krnames])
  }
  attr(r,"method") <- c("Kenward-Roger via pbkrtest package")
  r
}
# drop1(gm1,test="F",sumFun=KRSumFun)


#####################################################################################################################################
#  Main data based on direct counts of abundances and summed biomass within 
# the mesocosms (35 mesococsm smaples on three data, weeks 1, 2, and 4),
#  data on pest control of aphids (Aphis fabe) and counts of the key ladybird predator (Coccinella septempunctata). 
# mesocosms and sample date (weeks 1-4)
setwd("P:\\NEC06567 NERC-RP Synergistic mixture effects and sensitisation\\WP4_Synergistic effects in the field\\CHEMMIX_year 1 paper FINAL\\")
dir()
Main_data<-read.csv("Main_data_R_week1to4.csv")  #   other directly monitored data
Main_data$BLOCK <-as.factor(Main_data$BLOCK )
# Reorder levels
Main_data$Treatment <- factor(Main_data$Treatment, levels = c("Cont.", "Cy.", "Az.", "Pr.", "Pr.X.Az", "Cy.X.Az", "Cy.X.Pr"))



#  Subsets of data for tests of synergisms based on the principal of response additivity assessed by testing linear interaction effect 
#  Slinker, B.K. (1998) The Statistics of Synergism. Journal of Molecular and Cellular Cardiology, 30, 723-731.

web_metrics$AZO_0.5<-as.factor(web_metrics$AZO_0.5)	#   Where AZO dicatates the presence or absence of azoxystrobin
web_metrics$PCZ_0.5<-as.factor(web_metrics$PCZ_0.5)	#   Where PCZ dicatates the presence or absence of prochloraz
web_metrics$CYP_0.5<-as.factor(web_metrics$CYP_0.5)	#   Where CYP dicatates the presence or absence of cypermethrin
# We will check for synergisms only in the final week
Main_data_W4<-Main_data %>% filter(Sampling_week %in% c("W4"))

Main_data_W4_Cy_Az<- Main_data_W4 %>% filter(Treatment %in% c("Cont.", "Cy.", "Az.","Cy.X.Az"))
Main_data_W4_Cy_Pr<- Main_data_W4 %>% filter(Treatment %in% c("Cont.", "Cy.", "Pr.","Cy.X.Pr"))
Main_data_W4_Pr_Az<- Main_data_W4 %>% filter(Treatment %in% c("Cont.", "Pr.", "Az.","Pr.X.Az"))

#-----------------------------------------------
# Summary model structure

# subject to change of error distributions etc, but these are typical fixed (Treatment with 7 levels) and random effect (BLOCK and 
# nested within this individual Mesocosm to account for the repeated measures).   Decisiosn on th error distribution
# are made on a case by case basis in response to model checks done in Dharma. These are applied in more than one 
# package, for exampel lme4 for poisson and binomial, but Car using the glmmTMB function for negative binomial.
# Models including data measures at weeks 1,2 and 4 

gm1 <-glmer(Total__Aphis_Fabae  ~ Treatment + Sampling_week + Sampling_week*Treatment  +(1|BLOCK/Mesocosm), family=poisson, data=Main_data)
# models which have  single time point estiamte, e.g. earthworm coccon production is measured once in week 8
gm1 <-glmer(Coccons ~ Treatment  +(1|BLOCK), family=poisson, data=earthworms)
# model simplification is by deletion of lease sig effects
Anova(gm1, type = "III")  # Type III analysis of fixed effects, using 'car' package, used for non-normal models
#drop1(gm1,test="user",sumFun=KRSumFun)   # using F tests  with Type III analysis - model norm
# checks of model fit are undertaken in Dharma package
sim_res <- simulateResiduals(fittedModel = gm1)  #  DHARMa tests of goodness of fit
plot(sim_res)
summary(gm1)
emmeans(gm1, ~ Treatment) 




#--------------------------------------------------------------------------------------------
#Model for abundance and mass community metrics from the mesocosms
#----------------------------------------------------------------------------------------

#Aphid abundance  (wingless adults)
gm1 <-glmmTMB(Wingless_Aphis_Fabae ~ Treatment + Sampling_week  +(1|BLOCK/Mesocosm), family=nbinom2, data=Main_data)
Anova(gm1, type = "III")  # Type III analysis of fixed effects, using 'car' package, used for non-normal models
#drop1(gm1,test="user",sumFun=KRSumFun)   # using F tests  with Type III analysis - model norm
sim_res <- simulateResiduals(fittedModel = gm1)  #  DHARMa tests of goodness of fit
plot(sim_res) # some quantile issues but accepatable
emmeans(gm1, ~ Sampling_week*Treatment) 


ggplot(data=Main_data, aes(y=log(Wingless_Aphis_Fabae+1), x=interaction(Sampling_week, Treatment), fill=Treatment)) +  # Color by treatment
  geom_boxplot() +
  theme_bw()+   #  white background
  labs(
    x = "Treatment",                                           # X-axis legend
    y = "Aphid abundance (log N+1)",                                        # Y-axis legend
    fill = "Treatment"                                    # Legend title
  )  +
  theme(
    panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1) , # Rotate x-axis text
    axis.text = element_text(size = 12),                      # Adjust axis text size
    axis.title = element_text(size = 16)                      # Adjust axis title size
  )        


# test of synergism for aphid abundance
# Azoxyxystrobin x Cypermethrin
gm1 <-glmmTMB(Wingless_Aphis_Fabae ~ AZO_0.5+ CYP_0.5 + AZO_0.5*CYP_0.5  +(1|BLOCK), family=nbinom2, data=Main_data_W4_Cy_Az)
Anova(gm1, type = "III") 
summary(gm1)

# Cypermethrin X prochloraz
gm1 <-glmmTMB(Wingless_Aphis_Fabae ~ PCZ_0.5+ CYP_0.5 + PCZ_0.5*CYP_0.5  +(1|BLOCK), family=nbinom2, data=Main_data_W4_Cy_Pr)
Anova(gm1, type = "III") 
summary(gm1)

#Prochloraz X Azoxystroin
gm1 <-glmmTMB(Wingless_Aphis_Fabae~ AZO_0.5+ PCZ_0.5 + AZO_0.5*PCZ_0.5  +(1|BLOCK), family=nbinom2, data=Main_data_W4_Pr_Az)
Anova(gm1, type = "III") 
summary(gm1)




#-----------------------------------------------------------------
#Ladybird abundance
# Very variable  -  may be problematic as lots of huge outliers
gm1 <-glmer(Ladybirds ~ Treatment + Sampling_week + Sampling_week*Treatment  +(1|BLOCK/Mesocosm), family=poisson, data=Main_data)
# try simplifying the model losing higher order interactions to assess if there is a treatment effect
gm1 <-glmer(Ladybirds ~ Treatment  +(1|BLOCK/Mesocosm), family=poisson, data=Main_data)
Anova(gm1, type = "III")  # Type III analysis of fixed effects, using 'car' package, used for non-normal models
#drop1(gm1,test="user",sumFun=KRSumFun)   # using F tests  with Type III analysis - model norm
sim_res <- simulateResiduals(fittedModel = gm1)  #  DHARMa tests of goodness of fit
# KS test suggest some deviation of residual from uniformity (p=0.04)
# will try neg.bn distribution to see if over dispersion is an issue.

gm1 <-glmmTMB(Ladybirds ~ Treatment + Sampling_week + Sampling_week*Treatment  +(1|BLOCK/Mesocosm), family=nbinom2, data=Main_data)  #  does not converge
Anova(gm1, type = "III")  # Type III analysis of fixed effects, using 'car' package, used for non-normal models
#drop1(gm1,test="user",sumFun=KRSumFun)   # using F tests  with Type III analysis - model norm
sim_res <- simulateResiduals(fittedModel = gm1)  #  DHARMa tests of goodness of fit
plot(sim_res)    # this looks fine now
summary(gm1)
emmeans(gm1, ~ Sampling_week*Treatment) 






ggplot(data=Main_data, aes(y=log(Ladybirds+1), x=interaction(Sampling_week, Treatment), fill=Treatment)) +  # Color by treatment
  geom_boxplot() +
  theme_bw()+   #  white background
  labs(
    x = "Treatment",                                           # X-axis legend
    y = "Ladybird abundance (log N+1)",                                        # Y-axis legend
    fill = "Treatment"                                    # Legend title
  )  +
  theme(
    panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1) , # Rotate x-axis text
    axis.text = element_text(size = 12),                      # Adjust axis text size
    axis.title = element_text(size = 16)                      # Adjust axis title size
  )      

# test of synergism for ladybird abundance
# cypermethrin X azoxystrobin
gm1 <-glmmTMB(Ladybirds  ~ AZO_0.5+ CYP_0.5 + AZO_0.5*CYP_0.5  +(1|BLOCK), family=nbinom2, data=Main_data_W4_Cy_Az)
Anova(gm1, type = "III") 

# Cypermethrin X prochloraz
gm1 <-glmmTMB(Ladybirds  ~ PCZ_0.5+ CYP_0.5 + PCZ_0.5*CYP_0.5  +(1|BLOCK), family=nbinom2, data=Main_data_W4_Cy_Pr)
Anova(gm1, type = "III") 

#Prochloraz X Azoxystroin
gm1 <-glmmTMB(Ladybirds  ~ AZO_0.5+ PCZ_0.5 + AZO_0.5*PCZ_0.5  +(1|BLOCK), family=nbinom2, data=Main_data_W4_Pr_Az)
Anova(gm1, type = "III") 


#--------------------------------------------------------------------------------------
# detritivore biomass
gm1 <-lmer(sum_det_mass_mg  ~   Sampling_week  +(1|BLOCK/Mesocosm), data=Main_data)
#Anova(gm1, type = "III")  # Type III analysis of fixed effects, using 'car' package, used for non-normal models
drop1(gm1,test="user",sumFun=KRSumFun)   # using F tests  with Type III analysis - model norm
sim_res <- simulateResiduals(fittedModel = gm1)  #  DHARMa tests of goodness of fit
plot(sim_res)    # KS deviation significant -  less of an issue for gausian models.
# try a transformation to improve fit.
gm1 <-lmer(log(sum_det_mass_mg+1)~ Sampling_week + (1|BLOCK/Mesocosm), data=Main_data)
drop1(gm1,test="user",sumFun=KRSumFun)   # using F tests  with Type III analysis - model norm
sim_res <- simulateResiduals(fittedModel = gm1)  #  DHARMa tests of goodness of fit
plot(sim_res)    #This has significanly improved the fit of this model.  Looks acceptabel with the log mass transformation
summary(gm1)
emmeans(gm1, ~ Sampling_week) 


ggplot(data=Main_data, aes(y=log(sum_det_mass_mg+1), x=Sampling_week, fill=Sampling_week)) +  # Color by treatment
  geom_boxplot() +
  theme_bw()+   #  white background
  labs(
    x = "Sampling week",                                           # X-axis legend
    y = "Mass (log mg)",                                        # Y-axis legend
    fill = "Treatment"                                    # Legend title
  )  +
  theme(
    panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
    axis.text = element_text(size = 12),                      # Adjust axis text size
    axis.title = element_text(size = 16)                      # Adjust axis title size
  )    

# no test for synergism as no overal treatment effect

#--------------------------------------------------------------------------------------
# herbivore biomass
gm1 <-lmer(sum_herb_mass_mg~ Treatment + Sampling_week +(1|BLOCK/Mesocosm), data=Main_data)
drop1(gm1,test="user",sumFun=KRSumFun)   # using F tests  with Type III analysis - model norm
sim_res <- simulateResiduals(fittedModel = gm1)  #  DHARMa tests of goodness of fit
plot(sim_res)    # very strong quantile deviations detected
# try a transformation on the response log(N+1)
gm1 <-lmer(log(sum_herb_mass_mg+1)~ Treatment + Sampling_week + Treatment*Sampling_week +  (1|BLOCK/Mesocosm), data=Main_data)
drop1(gm1,test="user",sumFun=KRSumFun) 
sim_res <- simulateResiduals(fittedModel = gm1)  #  DHARMa tests of goodness of fit
plot(sim_res)    # Measures of model fit all look fine now
summary(gm1)
emmeans(gm1, ~ Treatment*Sampling_week) 
emmeans(gm1, ~ Treatment) 

ggplot(data=Main_data, aes(y=log(sum_herb_mass_mg+1), x=interaction(Sampling_week, Treatment), fill=Treatment)) +  # Color by treatment
  geom_boxplot() +
  theme_bw()+   #  white background
  labs(
    x = "Treatment",                                           # X-axis legend
    y = "Mass (log mg)",                                        # Y-axis legend
    fill = "Treatment"                                    # Legend title
  )  +
  theme(
    panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1) , # Rotate x-axis text
    axis.text = element_text(size = 12),                      # Adjust axis text size
    axis.title = element_text(size = 16)                      # Adjust axis title size
  )   


# test of synergism for herbivore mass
# cypermethrin X azoxystrobin
gm1 <-lmer(log(sum_herb_mass_mg+1)  ~ AZO_0.5+ CYP_0.5 + AZO_0.5*CYP_0.5  +(1|BLOCK), data=Main_data_W4_Cy_Az)
drop1(gm1,test="user",sumFun=KRSumFun)

# Cypermethrin X prochloraz
gm1 <-lmer(log(sum_herb_mass_mg+1)  ~ PCZ_0.5+ CYP_0.5 + PCZ_0.5*CYP_0.5  +(1|BLOCK), data=Main_data_W4_Cy_Pr)
drop1(gm1,test="user",sumFun=KRSumFun)

#Prochloraz X Azoxystroin
gm1 <-lmer(log(sum_herb_mass_mg+1)  ~ AZO_0.5+ PCZ_0.5 + AZO_0.5*PCZ_0.5  +(1|BLOCK),data=Main_data_W4_Pr_Az)
drop1(gm1,test="user",sumFun=KRSumFun)



#--------------------------------------------------------------------------------------
# Predator biomass
gm1 <-lmer(sum_pred_mass_mg  ~ Treatment + Sampling_week + Sampling_week*Treatment   +(1|BLOCK/Mesocosm), data=Main_data)
#Anova(gm1, type = "III")  # Type III analysis of fixed effects, using 'car' package, used for non-normal models
drop1(gm1,test="user",sumFun=KRSumFun)   # using F tests  with Type III analysis - model norm
sim_res <- simulateResiduals(fittedModel = gm1)  #  DHARMa tests of goodness of fit
plot(sim_res)    # KS deviation significant -  quantile deviations detected
# try tranformation ot improve model fit
gm1 <-lmer(log(sum_pred_mass_mg+1)  ~ Treatment + Sampling_week + Sampling_week*Treatment   +(1|BLOCK/Mesocosm), data=Main_data)
drop1(gm1,test="user",sumFun=KRSumFun)   # using F tests  with Type III analysis - model norm
sim_res <- simulateResiduals(fittedModel = gm1)  #  DHARMa tests of goodness of fit
plot(sim_res)    # Transformation worked well.  Model assumption checks are ok
plot(gm1)
summary(gm1)
emmeans(gm1, ~ Treatment*Sampling_week) 

ggplot(data=Main_data, aes(y=log(sum_pred_mass_mg+1), x=interaction(Sampling_week, Treatment), fill=Treatment)) +  # Color by treatment
  geom_boxplot() +
  theme_bw()+   #  white background
  labs(
    x = "Treatment",                                           # X-axis legend
    y = "Mass (log mg)",                                        # Y-axis legend
    fill = "Treatment"                                    # Legend title
  )  +
  theme(
    panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1) , # Rotate x-axis text
    axis.text = element_text(size = 12),                      # Adjust axis text size
    axis.title = element_text(size = 16)                      # Adjust axis title size
  )   


# test of synergism for predator mass
# cypermethrin X azoxystrobin
gm1 <-lmer(log(sum_pred_mass_mg+1)  ~ AZO_0.5+ CYP_0.5 + AZO_0.5*CYP_0.5  +(1|BLOCK), data=Main_data_W4_Cy_Az)
drop1(gm1,test="user",sumFun=KRSumFun)

# Cypermethrin X prochloraz
gm1 <-lmer(log(sum_pred_mass_mg+1)  ~ PCZ_0.5+ CYP_0.5 + PCZ_0.5*CYP_0.5  +(1|BLOCK), data=Main_data_W4_Cy_Pr)
drop1(gm1,test="user",sumFun=KRSumFun)

#Prochloraz X Azoxystroin
gm1 <-lmer(log(sum_pred_mass_mg+1)  ~ AZO_0.5+ PCZ_0.5 + AZO_0.5*PCZ_0.5  +(1|BLOCK),  data=Main_data_W4_Pr_Az)
drop1(gm1,test="user",sumFun=KRSumFun)

##############################################################################
earthworms<-read.csv("R_earthworms_week8.csv")
earthworms$BLOCK <-as.factor(earthworms$BLOCK )
# Reorder levels
earthworms$Treatment <- factor(earthworms$Treatment, levels = c("Cont.", "Cy.", "Az.", "Pr.", "Pr.X.Az", "Cy.X.Az", "Cy.X.Pr"))

gm1 <-glmer(Coccons ~ Treatment  +(1|BLOCK), family=poisson, data=earthworms)
# failed to converge
Anova(gm1, type = "III")  # Type III analysis of fixed effects, using 'car' package, used for non-normal models
#drop1(gm1,test="user",sumFun=KRSumFun)   # using F tests  with Type III analysis - model norm
sim_res <- simulateResiduals(fittedModel = gm1)  #  DHARMa tests of goodness of fit
plot(sim_res)
# KS test suggest some deviation of residual from uniformity and dispersion test failed
# will try neg.bn distribution to see if over dispersion is an issue.
gm1 <-glmmTMB(Coccons ~ Treatment + (1|BLOCK), family=nbinom2, data=earthworms)  #  does not converge
Anova(gm1, type = "III")  # Type III analysis of fixed effects, using 'car' package, used for non-normal models
# Non significant effect of treatment
sim_res <- simulateResiduals(fittedModel = gm1)  #  DHARMa tests of goodness of fit
plot(sim_res)    # this looks fine now
summary(gm1)
emmeans(gm1, ~ Treatment) 



ggplot(data=earthworms, aes(y=Coccons, x=Treatment, fill=Treatment)) +  # Color by treatment
  geom_boxplot() +
  theme_bw()+   #  white background
  labs(
    x = "Treatment",                                           # X-axis legend
    y = "Earthworm Coccon abundance",                                        # Y-axis legend
    fill = "Treatment"                                    # Legend title
  )  +
  theme(
    panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1) , # Rotate x-axis text
    axis.text = element_text(size = 14),                      # Adjust axis text size
    axis.title = element_text(size = 16)                      # Adjust axis title size
  )      


# no test for synergism as no overal treatment effect

##############################################################################
decomp<-read.csv("R_Decomposition.csv")   #   bait lamina decomposition assessments for week 4
decomp$BLOCK <-as.factor(decomp$BLOCK )
# Reorder levels
decomp$Treatment <- factor(decomp$Treatment, levels = c("Cont.", "Cy.", "Az.", "Pr.", "Pr.X.Az", "Cy.X.Az", "Cy.X.Pr"))

gm1 <-lmer(Total_average_number_removed~ Treatment  +(1|BLOCK),  data=decomp)
Anova(gm1, type = "III")  # Type III analysis of fixed effects, using 'car' package, used for non-normal models
drop1(gm1,test="user",sumFun=KRSumFun)   # using F tests  with Type III analysis - model norm
sim_res <- simulateResiduals(fittedModel = gm1)  #  DHARMa tests of goodness of fit
plot(sim_res)

summary(gm1)


ggplot(data=decomp, aes(y=Total_average_number_removed, x=Treatment, fill=Treatment)) +  # Color by treatment
  geom_boxplot() +
  theme_bw()+   #  white background
  labs(
    x = "Treatment",                                           # X-axis legend
    y = "Decomposition rate",                                        # Y-axis legend
    fill = "Treatment"                                    # Legend title
  )  +
  theme(
    panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1) , # Rotate x-axis text
    axis.text = element_text(size = 14),                      # Adjust axis text size
    axis.title = element_text(size = 16)                      # Adjust axis title size
  )


# no test for synergism as no overal treatment effect

#####################################################################################

library(cheddar) # Food web analysis
library(dplyr) #  data processing
library(plyr) # data processing
library(tibble)  #  dtat processing - row names to a column
library(purrr) #  map2 function for food web manipulation

#  Using the WebBuilder package, observed feeding relationships, published feeding relationships
# and taxonomy based generalisations about how consumers feed and resources are used
# we create predicted food webs of the arthropods monitored within the mesocosms
# after 4 weeks of exposure to the treatments.
# See below ref for Webuilder package
# Gray et al 2015. Joining the dots: An automated method for constructing food webs from compendia of published interactions. Food Webs 5, 11-20.

setwd("P:\\NEC06567 NERC-RP Synergistic mixture effects and sensitisation\\WP4_Synergistic effects in the field\\DATA\\CHEMMIX_year 1 paper FINAL\\")
dir()
source("WebBuilder.R")    #  WebBuilder package  -  see Gray et al 2015 for details.

properties<- list("Mixtures", "mg", "m^-2")   #  specifying titles of arthropod data
names(properties) <- c("title", "M.units", "N.units")  #  Identifying units for arthropod data for WebBuilder

# Establish your registry of feeding relationships
registry<-read.csv("Registry_known_links.csv")

# Directly observed trophic links within the mesocosms
trophic.links_directly_Observed<-read.csv("trophic.links.csv")  #  based on overall list -  need to subset


# Mesocosm meta data
treatment<-as.matrix(read.csv("Treatments.csv"))

# Resource consumer generalization list
# This is the list of  trophic generalizations, i.e. 
#would a consumer of this species be expected to consume any other
# species in the same genus, same family, same order or same class.
# generalization of the consumer and resource (minimum.res.method & minimum.cons.method)
res_cons_list<-read.csv("res_cons_list_week4.csv")

# Read in overall mesocosm (35 samples) for a single mesocosms sampled week 4
# 8 treatment levels of all combinations of with and without cypermethrin, prochloraz and azoxystrobin
# 5 replicate blocks
# The subsequent analysis will subset each mesocosms -  remove zero values for species
# and implement the webBuilder package to infer trophic interactions, from this
# we will use the Cheddar package to derive metrics of trophic strucutre.  This subseted data will
# represent the Nodes of the subsequent webs (e.g. a node being a species with a feeding relationship
# with another species as predicted within WebBuilder).
# NOTE - Nomenclature is at the end of this data table at L41-49
#--------------------------------------------------------------------
#  this will be removed in subsequent processing.

nodes_all_plots_dates<-read.csv("nodes_week4.csv")   # node data frame containing all species and all N values for the 120 plots x sample dates
nodes_just_mass<-nodes_all_plots_dates[-(1:35), ]  #  subsets taxonomy and mass data for each species L121-129. 

# create data matrices to store the web metrics to be produced using the 
# Cheadar pckage from each of the individually infered food webs (120 total).
Web_metrics<-matrix(ncol =8, nrow = 35)  #  empty data frame to put in the web metrics for each combination of sites and iterations
web_metrics<-cbind(Web_metrics,treatment)
colnames(web_metrics) <- c("Convex_hull_area", "allom_degree_dist_int", "allom_degree_dist_slope", "NvM_Slope", "NvM_Int", "SpeciesRich", "Linkages", "Connectacne",
                           "Plot_code", "Sampling_week",	"Treatment",	"Mesocosm",	"BLOCK",	"AZO_0.5",	"PCZ_0.5",	"CYP_0.5")
head(web_metrics)

# The following loop: 1) subsets the overall dataset containing all 120 mesocosms (40 mesocosms X 3 sample dates)
# to a single mesocosm (mesocode 1 to 35) and prunes this data and prepares it for use in deriving the 
# predicted food webs using the WebBuilder package; 2) For each web implements the taxonomic 
# generalisation combined with known and literature based trophic links between species to predict
# trophic interactions and so produce a mesocosms and date and formats this for the Cheddar package; 3)
# Derives metrics of food we structure and enters them into the pre-prepared matrix web_metrics
# for subsequent analysis using GLMM.


mesocode<- 1   # unique number for each mesocosm x sample date
for (mesocode in 1:35) {
print(paste("mesocode count out of 35:", mesocode))  
N_for_mesocosmXdate<-nodes_all_plots_dates[mesocode, ]
nodes<-rbind(N_for_mesocosmXdate,nodes_just_mass)   # combining the plot specific N data bask in with the rest of the node data
nodes<-t(nodes[, -(1:5)])   # pruning of the columns with NA values and then transposing the data
new_colnames <- as.character(unlist(nodes[1, ]))  #need to move the top row to become column names
nodes <- nodes[-1, ]   ## Step 2: Remove the top row from the data frame
colnames(nodes ) <- new_colnames# Step 3: Assign the new column names to the data frame
nodes<-cbind(rownames(nodes ), nodes)   #  add a column with the values in the row titles
colnames(nodes)[1] <- "node" # give a name to this newly created column
rownames(nodes) <- NULL    # remove row names
nodes<-as.data.frame(nodes)   #convert to a data frame
nodes <- nodes %>%    filter(N != 0)   # remove 0 roes where there is a 0 value for abundance (row N)
nodes$M<-as.numeric(nodes$M)
nodes$N<-as.numeric(nodes$N)

#subset known trophic links by the nodes present in the mesocosm  x date sub-data set
# first this is done for any of the species that are 'resources' in the known trophic.links file
x<-nodes
x$resource<-x$node    #  so create a common list of species by which we can subset trophic links by the nodes present in the data set, e.g. resource
x2<-merge(x,trophic.links_directly_Observed,by="resource", all=FALSE) # merges the nodes list of species in a mesocosm (under DF x)
# with the species list of known trophic links (trophic.links_directly_Observed) but excludes any species not present in the nodes DF.
x2<-cbind(x2$resource, x2$consumer)   #  new DF (X2) that has sub-setted all possible trophic links (trophic.links_directly_Observed)
# by the species (nodes) that are present within a particular mesocosm.
colnames(x2) <- c("resource", "consumer") #  effectively DF x2 is the trophic links DF now.
# Now repeat the process for any of the species in the nodes field (called temporarily DF 'x') 
# for all the consumer species in the known trophic.links_directly_Observed file
x<-nodes
x$consumer<-x$node
x3<-merge(x,x2,by="consumer", all=FALSE) #  note we are sub-setting the new trophic.link fiel (called DF x2) which
# contains only those resource species found in the mesocosms (as defined by the nodes DF) by those predatory  
# consumer species also in the mesocosm. i.e. its list of resource species present in a mesocosm filtered to include only those predated by a 
# consumer also in that mesocosm.
trophic.links_directly_Observed_subset<-cbind(x3$resource, x3$consumer) #  this is now defined as the new tailored known trophic links 
# for the species of resources and consumers in that mesocosm.  Called DF trophic.links_directly_Observed_subset
colnames(trophic.links_directly_Observed_subset) <- c("resource", "consumer") #  gives correct column titles

# merges the sub-setted mesocosm x data nodes and the directly observed feeding relationships in a format for WebBuilder
Mixtures.exp<- Community(properties=properties,  nodes=nodes,   trophic.links=trophic.links_directly_Observed_subset)  # merges these in a format for WebBuilder
# If you had no trophic.links_directly_Observed data then the alternative code for this would be
# Mixtures.exp<- Community(properties=properties,  nodes=nodes)
Mixtures.exp<- OrderCommunity(Mixtures.exp, 'M') #  order by body mass
# removing nodes with no taxonomic resolution -  e.g. organic matter
#Mixtures.exp<-RemoveNodes(community=Mixtures.exp, remove=c('Organic_matter')) 
nps<-NPS(Mixtures.exp) # pulls out the nodes in the WebBuilder  format (with nodal species as rows)
nodes<-cbind(nps[,c(1,7:10)])   # subsets this to produce a new nodes file containing only node, class, order, family, genus,

#subset resource -  consumer list for the nodes present in the data set.  
nodes<-merge(nodes,res_cons_list,by="node", all=FALSE) # new nodes field which includes the level of trophic

# Now use WebBuilder to infer the links from the nodes (and levels of trophic generalization within) and the registry
# of published known feeding relationships
links.inf<-WebBuilder(nodes, registry, method=c('exact','genus','family','order','class'))   
# now combine these into a new community file for use in the Cheddar package.  Note CPS(Mixtures.exp) gives the units,
# e.g. title, mass units and abundance units as defined above.
Mixtures.infered <- Community(properties = CPS(Mixtures.exp) ,   nodes = NPS(Mixtures.exp),  trophic.links = links.inf)
Mixtures.noVicia<-RemoveNodes(community=Mixtures.infered, remove=c('Vicia_faba'))  # remove the plant species as no actual biomass

#Metrics
#Convex hull Leaper, R. and Raffaelli, D. (1999) Defining the abundance body-size constraint space: data from a real food web. Ecology Letters 2, 3, 191–199.
Convex_hull<- NvMConvexHull(Mixtures.noVicia)# Returns the points and area of the minimum convex hull (a polygon in log10-transformed numerical abundance versus log10-transformed body mass space) that bounds all the species within the
web_metrics[mesocode,1]<-Convex_hull$area
# allometric degree distributions
# how species' numbers of trophic links scale with their log-transformed body masses.
# Jonsson et al., 2005; Otto et al., 2007; Digel et al., 2011
allom_degree_dist<- lm(OutDegree(Mixtures.noVicia) ~ Log10M(Mixtures.noVicia))
web_metrics[mesocode,2]<-allom_degree_dist$coefficients[1]
web_metrics[mesocode,3]<-allom_degree_dist$coefficients[2]
web_metrics[mesocode,4]<-NvMSlope(Mixtures.noVicia) # linear regressions fitted to log10- transformed numerical abundance versus log10-transformed body mass.
web_metrics[mesocode,5]<-NvMIntercept(Mixtures.noVicia) # linear regressions fitted to log10- transformed numerical abundance versus log10-transformed body mass.
# summary statistics
summary_web<-summary(Mixtures.noVicia) 
web_metrics[mesocode,6]<-summary_web$S    # number of species
web_metrics[mesocode,7]<-if (is.null(summary_web$L)) 0 else summary_web$L     # Number of trophic links -  set to 0 if none, e.g. just detritivores  
web_metrics[mesocode,8]<-if (is.null(summary_web$C)) 0 else summary_web$C     # Connectance -  set to 0 if no trophic links  
}

write.csv(web_metrics, "web_metrics.csv")


###############################################################################
# Web metric responses to treatments
# Based on derived web structure using WebBuilder
# Focuses on web structure assessed 4 weeks after exposure to pesticides.

#read in exported web metric file
setwd("P:\\NEC06567 NERC-RP Synergistic mixture effects and sensitisation\\WP4_Synergistic effects in the field\\DATA\\CHEMMIX_year 1 paper\\")
dir()
web_metrics<-read.csv("web_metrics.csv")# derived web metric data for each mesocosm
web_metrics$BLOCK <-as.factor(web_metrics$BLOCK )
web_metrics$X <-as.factor(web_metrics$X )   #  observation level classifier -  i.e. unique code for each of the 35 mesocosms

head(web_metrics)

# Reorder levels
web_metrics$Treatment <- factor(web_metrics$Treatment, levels = c("Cont.", "Cy.", "Az.", "Pr.", 
                                                                  "Az.X.Pr", "Cy.X.Az", "Cy.X.Pr"))



#  Subsets of data for tests of synergisms based on the principal of response additivity assessed by testing linear interaction effect 
#  Slinker, B.K. (1998) The Statistics of Synergism. Journal of Molecular and Cellular Cardiology, 30, 723-731.
web_metrics$AZO_0.5<-as.factor(web_metrics$AZO_0.5)	#   Where AZO dicatates the presence or absence of azoxystrobin
web_metrics$PCZ_0.5<-as.factor(web_metrics$PCZ_0.5)		#   Where PCZ dicatates the presence or absence of prochloraz
web_metrics$CYP_0.5<-as.factor(web_metrics$CYP_0.5)	#   Where CYP dicatates the presence or absence of cypermethrin
web_metrics_Cy_Az<- web_metrics %>% filter(Treatment %in% c("Cont.", "Cy.", "Az.","Cy.X.Az"))
web_metrics_Cy_Pr<- web_metrics %>% filter(Treatment %in% c("Cont.", "Cy.", "Pr.","Cy.X.Pr"))
web_metrics_Pr_Az<- web_metrics %>% filter(Treatment %in% c("Cont.", "Pr.", "Az.","Az.X.Pr"))
  
# Check the new levels
levels(web_metrics$Treatment)

# Number of nodes (Species Richness in web)
gm1 <-glmer(SpeciesRich   ~ Treatment  +(1|BLOCK), family=poisson,  data=web_metrics)
Anova(gm1, type = "III")  # Type III analysis of fixed effects, using 'car' package, used for non-normal models
sim_res <- simulateResiduals(fittedModel = gm1)  #  DHARMa tests of goodness of fit
plot(sim_res)   # evidence of overdipersion (p<0.05)
#try neg.bin
gm1_nb <-glmmTMB(SpeciesRich   ~ Treatment  +(1|BLOCK), family=nbinom2,  data=web_metrics)
Anova(gm1_nb, type = "III")  # Type III analysis of fixed effects, using 'car' package, used for non-normal models
sim_res <- simulateResiduals(fittedModel = gm1_nb)  #  DHARMa tests of goodness of fit
plot(sim_res)   # Still evidence of overdipersion (p<0.05)
#    adding an observation-level random effect to models overdispersion as individual variation.
gm1_indivrandeffect<-gm1 <-glmer(SpeciesRich   ~ Treatment  +(1|BLOCK)+(1|Mesocosm), family=poisson,  data=web_metrics)
sim_res <- simulateResiduals(fittedModel = gm1_indivrandeffect)  #  DHARMa tests of goodness of fit
plot(sim_res) # evidence of overddipersion
# Compare AIC values
AIC(gm1, gm1_nb, gm1_indivrandeffect)   #  all similar
# all evidence of overdispersion -  based on dhamra analytics poisson is best.

gm1 <-glmer(SpeciesRich   ~ Treatment  +(1|BLOCK), family=poisson,  data=web_metrics)
Anova(gm1, type = "III")  # Type III analysis of fixed effects, using 'car' package, used for non-normal models
emmeans(gm1, ~ Treatment) 


ggplot(data=web_metrics, aes(y=SpeciesRich, x=Treatment, fill=Treatment)) +  # Color by treatment
  geom_boxplot() +
  theme_bw()+   #  white background
  labs(
    x = "Treatment",                                           # X-axis legend
    y = "Nodes (species richness)",                                        # Y-axis legend
    fill = "Treatment"                                    # Legend title
  )  +
  theme(
    panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1) , # Rotate x-axis text
    axis.text = element_text(size = 14),                      # Adjust axis text size
    axis.title = element_text(size = 16)                      # Adjust axis title size
  )                                         


# No tests for synergisms undertaken as no overal treatment effect identified
  

#-----------------------------------------------------------------------------------------------------
#Connectance
gm1 <-lmer(Connectacne~ Treatment  +(1|BLOCK),  data=web_metrics)
drop1(gm1,test="user",sumFun=KRSumFun)   # using F tests  with Type III analysis - model normal
sim_res <- simulateResiduals(fittedModel = gm1)  #  DHARMa tests of goodness of fit
plot(sim_res)
summary(gm1)
# azXcy treatment 
1 - pt(2.37, 5)     #  pt() is the Cumulative Distribution Function (CDF) of the t-distribution. It calculates the probability that a t-distributed random variable is less than or equal to a given value.
emmeans(gm1, ~ Treatment) 
# Create the box plot
ggplot(data=web_metrics, aes(y=Connectacne, x=Treatment, fill=Treatment)) +  # Color by treatment
  geom_boxplot() +
  theme_bw()+   #  white background
  labs(
    x = "Treatment",                                           # X-axis legend
    y = "Web connectance",                                        # Y-axis legend
    fill = "Treatment"                                    # Legend title
  )  +
  theme(
    panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1) , # Rotate x-axis text
    axis.text = element_text(size = 14),                      # Adjust axis text size
    axis.title = element_text(size = 16)                      # Adjust axis title size
  )                                         

# test of synergism for connectance in food web
# Azoxyxystrobin x Cypermethrin
gm1 <-lmer(Connectacne~ AZO_0.5+ CYP_0.5 + AZO_0.5*CYP_0.5 +(1|BLOCK),  data=web_metrics_Cy_Az)
drop1(gm1,test="user",sumFun=KRSumFun)   

# Cypermethrin X prochloraz
gm1 <-lmer(Connectacne~ PCZ_0.5+ CYP_0.5 + PCZ_0.5*CYP_0.5 +(1|BLOCK),  data=web_metrics_Cy_Pr)
drop1(gm1,test="user",sumFun=KRSumFun)  

#Prochloraz X Azoxystroin
gm1 <-lmer(Connectacne~ AZO_0.5+ PCZ_0.5 + AZO_0.5*PCZ_0.5 +(1|BLOCK),  data=web_metrics_Pr_Az)
drop1(gm1,test="user",sumFun=KRSumFun)   

#----------------------------------------------------------------------
# Convex hull - The points and area of the minimum convex hull (a polygon in log10-transformed numerical abundance versus log10-transformed body mass space) that bounds all the species
gm1 <-lmer(Convex_hull_area~ Treatment  +(1|BLOCK),  data=web_metrics)
drop1(gm1,test="user",sumFun=KRSumFun)   # using F tests  with Type III analysis - model norm
sim_res <- simulateResiduals(fittedModel = gm1)  #  DHARMa tests of goodness of fit
plot(sim_res)
# cy treatment 
1 - pt(2.33, 5)#  pt() is the Cumulative Distribution Function (CDF) of the t-distribution. It calculates the probability that a t-distributed random variable is less than or equal to a given value.
# azXcy treatment 
1 - pt(3.23, 5)#  pt() is the Cumulative Distribution Function (CDF) of the t-distribution. It calculates the probability that a t-distributed random variable is less than or equal to a given value.
# azXpr treatment 
1 - pt(1.98, 5)#  pt() is the Cumulative Distribution Function (CDF) of the t-distribution. It calculates the probability that a t-distributed random variable is less than or equal to a given value.

web_metrics$Treatment <- factor(web_metrics$Treatment, levels = c("Cy.", "Cont.", "Az.", "Pr.", 
                                                                   "Az.X.Pr", "Cy.X.Az", "Cy.X.Pr"))   
gm1 <-lmer(Convex_hull_area~ Treatment  +(1|BLOCK),  data=web_metrics)
summary(gm1)
1 - pt(0.90, 5)  #  pt() is the Cumulative Distribution Function (CDF) of the t-distribution. It calculates the probability that a t-distributed random variable is less than or equal to a given value.

summary(gm1)
emmeans(gm1, ~ Treatment) 
ggplot(data=web_metrics, aes(y=Convex_hull_area, x=Treatment, fill=Treatment)) +  # Color by treatment
  geom_boxplot() +
  theme_bw()+   #  white background
  labs(
    x = "Treatment",                                           # X-axis legend
    y = "Convex hull area",                                        # Y-axis legend
    fill = "Treatment"                                    # Legend title
  )  +
  theme(
    panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1) , # Rotate x-axis text
    axis.text = element_text(size = 14),                      # Adjust axis text size
    axis.title = element_text(size = 16)                      # Adjust axis title size
  )                                         

# test for synergism for convex hull area
# Azoxyxystrobin x Cypermethrin
gm1 <-lmer(Convex_hull_area~ AZO_0.5+ CYP_0.5 + AZO_0.5*CYP_0.5 +(1|BLOCK),  data=web_metrics_Cy_Az)
drop1(gm1,test="user",sumFun=KRSumFun)   

# Cypermethrin X prochloraz
gm1 <-lmer(Convex_hull_area~ PCZ_0.5+ CYP_0.5 + PCZ_0.5*CYP_0.5 +(1|BLOCK),  data=web_metrics_Cy_Pr)
drop1(gm1,test="user",sumFun=KRSumFun)  

#Prochloraz X Azoxystroin
gm1 <-lmer(Convex_hull_area~ AZO_0.5+ PCZ_0.5 + AZO_0.5*PCZ_0.5 +(1|BLOCK),  data=web_metrics_Pr_Az)
drop1(gm1,test="user",sumFun=KRSumFun)   

  
#----------------------------------------------------------------
# allom_degree_dist_slope
# describe how species numbers of trophic links scale with their log-transformed body masses.
gm1 <-lmer(allom_degree_dist_slope~ Treatment  +(1|BLOCK),  data=web_metrics)
drop1(gm1,test="user",sumFun=KRSumFun)   # using F tests  with Type III analysis - model norm
sim_res <- simulateResiduals(fittedModel = gm1)  #  DHARMa tests of goodness of fit
plot(sim_res)
summary(gm1)
emmeans(gm1, ~ Treatment) 

ggplot(data=web_metrics, aes(y=allom_degree_dist_slope, x=Treatment, fill=Treatment)) +  # Color by treatment
  geom_boxplot() +
  theme_bw()+   #  white background
  labs(
    x = "Treatment",                                           # X-axis legend
    y = "Allometric degree slope",                                        # Y-axis legend
    fill = "Treatment"                                    # Legend title
  )  +
  theme(
    panel.border = element_rect(color = "black", fill = NA),   # Solid axes lines
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1) , # Rotate x-axis text
    axis.text = element_text(size = 14),                      # Adjust axis text size
    axis.title = element_text(size = 16)                      # Adjust axis title size
  )                                         

  # note no overal significant treatment effect so no test for synergism


