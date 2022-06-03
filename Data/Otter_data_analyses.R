library(MASS)
library(psych)
library(ggplot2)
library(lme4)
library(lmerTest)
library(plyr)

##ICC analyses
latencyICC<-read.csv("Latency ICC.csv",header=T)
interactionICC<-read.csv("Interaction ICC.csv",header=T)
totaltimeICC<-read.csv("Total time ICC.csv",header=T)
ICC(latencyICC)
ICC(interactionICC)
ICC(totaltimeICC)

##Reading in all data
###Observation data
otters<-read.csv("Rock juggling frequency data.csv",header=T) 
otterdatav2<-subset(otters, total.obs>230) # removing one individual with only 200 observation bouts
summary(otterdatav2)
mean(otterdatav2$rj.total)
sd(otterdatav2$rj.total)
otterminbubble<-subset(otterdatav2, rj.total<59) # removing outlier, extreme juggling frequency 
otterminbubble
summary(otterminbubble)
mean(otterminbubble$rj.total)
sd(otterminbubble$rj.total) # comparing means+SDs to demonstrate extreme outlier

###Den camera trap data
dendata<-read.csv("Den camera trap data.csv",header=T)
den<-na.omit(dendata) #remove NAs
summary(den)

###Hunger data
otterhungrydata<-read.csv("Hunger data.csv", header=T)

###Puzzle toy data
otterpuzzledata<-read.csv("Puzzle data.csv",header=T) #outlier of 1200 seconds in latency variable was manually removed
otterpuzzledatav2<-na.omit(otterpuzzledata) #remove NAs
summary(otterpuzzledatav2)
otterpuzzledatav2$jugglerate = otterpuzzledatav2$rj.total/otterpuzzledatav2$total.obs #generating juggle rate
ottersolvedata<-subset(otterpuzzledatav2, solve=='1') #subset data to give only puzzles successfully solved
ascdata<-subset(otterpuzzledatav2,species=="ASC") #subset data to give only Asian small-clawed otter data
scodata<-subset(otterpuzzledatav2,species=="SCO") #subset data to give on smooth-coated otter data

###Rock juggling observation analyses

##Is otter rock juggling frequency correlated with age, age category(i.e. young and old), species, sex, wildlife park (site) or the interactions between age and: 
##age category, species and site?
model1<-glm.nb(rj.total~age*agecat+species+site+age:site+age:species+sex+offset(log(total.obs)), data=otterminbubble)
summary(model1)
plot(model1)

model2<-glm.nb(rj.total~age*agecat+species+site+age:site+sex+offset(log(total.obs)), data=otterminbubble)
summary(model2)
anova(model1,model2,test="Chisq") # including the age:species interaction does not significantly improve the model fit (p=1), so it is removed from the model 

model3<-glm.nb(rj.total~age*agecat+species+site+sex+offset(log(total.obs)), data=otterminbubble)
summary(model3)
anova(model2,model3,test="Chisq") # including the age:site interaction does not significantly improve the model fit (p=0.0863), so it is removed from the model

model4<-glm.nb(rj.total~age*agecat+species+sex+offset(log(total.obs)), data=otterminbubble)
summary(model4) # reported age and age cat p-values from this model 
anova(model3,model4,test="Chisq") # including site does not significantly improve the model fit (p=0.823), so it is removed from the model

model5<-glm.nb(rj.total~age:agecat+species+sex+offset(log(total.obs)), data=otterminbubble)
summary(model5)
anova(model4,model5,test="Chisq") # including age and agecat as independent predictors, as well as their interaction, does not significantly improve the fit of the model (p=0.274), so they are removed as independent fixed effects, and only their interaction is retained in the model.

model6<-glm.nb(rj.total~age:agecat+species+offset(log(total.obs)), data=otterminbubble)
summary(model6)
anova(model5,model6,test="Chisq") # including sex does not significantly improve the model fit (x1=2.416, p=0.120), so it is removed from the model

model7<-glm.nb(rj.total~age:agecat+offset(log(total.obs)), data=otterminbubble)
anova(model6,model7,test="Chisq") # including species does not significantly improve the model fit (x1=2.889, p=0.0892), so it is removed from the model
summary(model7) #reported age:agecat slope estimates, SE, z, from this model

null<-glm.nb(rj.total~1+offset(log(total.obs)), data=otterminbubble)
summary(null)
anova(null,model7,test="Chisq") #age:agecat, very significant,x2=24.044 p<0.001

###Den camera trap analyses

##Do otters rock juggle more when "hungry"?
count(den,c("hungry","juggling")) #obtaining counts to create contingency tables
hungry = matrix(c(17,366,8,165),nrow=2,ncol=2,byrow=TRUE)
dimnames(hungry) = list(c("H","F"), c("juggle","nojuggle"))
hungry
hungerchi<-chisq.test(hungry)
hungerchi #not significant, p=1
hungerchi$observed
hungerchi$expected
count(den,c("hungry")) #total no. of observations when "hungry" and "full"
#173+383

##Is rock juggling frequency affected by pup presence? - all den camera trap observations of rock juggling were on adult otters.
count(den,c("pups","juggling")) #obtaining counts to create contingency tables
pups = matrix(c(12,375,13,156),nrow=2,ncol=2,byrow=TRUE)
dimnames(pups) = list(c("pups","nopups"),c("juggle","nojuggle"))
pups
pupschi<-chisq.test(pups)
pupschi #significant p=0.0292
pupschi$observed
pupschi$expected #otters juggled less when pups present and more when pups were absent than expected
count(den,c("pups")) #total no. of observations with and without pups
#169+387

##Addressing confounding factor - are pups more likely to be present when "hungry"?
count(den,c("hungry","pups"))
puphungry = matrix(c(269,114,118,55),nrow=2,ncol=2,byrow=TRUE)
dimnames(puphungry) = list(c("H","F"),c("pups","nopups"))
puphungry
puphungrychi<-chisq.test(puphungry)
puphungrychi #not significant, p=0.703
puphungrychi$observed
puphungrychi$expected 

###Hunger data analyses - full observational data used (i.e. NOT camera trap data)
## Is otter rock juggling correlated with hunger?
otterhungrydata<-read.csv("Hunger data.csv", header=T)
library(lme4)
hungry<-glmer.nb(rj~hungry+offset(log(total.obs))+(1|ID), data=otterhungrydata) #ID as random effect
summary(hungry)
plot(hungry)
# null "hunger data" model to test if hunger significantly improves model fit 
nohungry<-glmer.nb(rj~1+offset(log(total.obs))+(1|ID),data=otterhungrydata)
summary(nohungry)
anova(hungry,nohungry,test="Chisq") #hunger significantly improves model fit (x1=4.761, p=0.0291)

###Puzzle data analyses (with "puzzle order" included)

##Does latency to first interact with puzzles correlate with juggle rate, species, the interaction between juggle rate and species, the type of puzzle...
##...the order in which puzzles were presented, the interaction between puzzle type and puzzle order, sex, age, or site?
puzzlelatency<-lmer(latency~jugglerate*species+puzzle.type*puzzle.order+sex+age+site+(1|group.ID/ID),data=otterpuzzledatav2)
summary(puzzlelatency) 

puzzlelatency2<-lmer(latency~jugglerate+species+puzzle.type*puzzle.order+sex+age+site+(1|group.ID/ID),data=otterpuzzledatav2)
anova(puzzlelatency,puzzlelatency2) #jugglerate:species interaction does not significantly improve model fit (x1=1.802, p=0.180), so it is removed from the model
summary(puzzlelatency2) 

puzzlelatency3<-lmer(latency~jugglerate+species+puzzle.type*puzzle.order+sex+age+(1|group.ID/ID),data=otterpuzzledatav2)
anova(puzzlelatency2,puzzlelatency3) #site does not significantly improve model fit (x2=1.398, p=0.497), so it is removed from the model
summary(puzzlelatency3)

puzzlelatency4<-lmer(latency~jugglerate+puzzle.type*puzzle.order+sex+age+(1|group.ID/ID),data=otterpuzzledatav2)
anova(puzzlelatency3,puzzlelatency4) #species does not significantly improve model fit (x1=0.048, p=0.827), so it is removed from the model
summary(puzzlelatency4) 

puzzlelatency5<-lmer(latency~puzzle.type*puzzle.order+sex+age+(1|group.ID/ID),data=otterpuzzledatav2)
anova(puzzlelatency4,puzzlelatency5) #jugglerate does not significantly improve model fit (x1=0.410, p=0.522), so it is removed from the model
summary(puzzlelatency5)

puzzlelatency6<-lmer(latency~puzzle.type*puzzle.order+sex+(1|group.ID/ID),data=otterpuzzledatav2)
anova(puzzlelatency5,puzzlelatency6) #age does not significantly improve model fit (x1=1.714, p=0.1905), so it is removed from the model
summary(puzzlelatency6) #MAM

puzzlelatency7<-lmer(latency~puzzle.type*puzzle.order+(1|group.ID/ID),data=otterpuzzledatav2)
anova(puzzlelatency6,puzzlelatency7) #sex significantly improves model fit (x1=7.130, p=0.008)

puzzlelatency8<-lmer(latency~puzzle.type+puzzle.order+sex+(1|group.ID/ID),data=otterpuzzledatav2)
anova(puzzlelatency6,puzzlelatency8) #puzzle.type:puzzle.order interaction significantly improves model fit (x4=57.352, p<0.001) 

##Post-hoc test to assess the interaction between puzzle type and puzzle order
puzzleordertype<-aov(latency~puzzle.type*puzzle.order,data=otterpuzzledatav2)
TukeyHSD(puzzleordertype) #Results are only significant for factors including "ball:first". There was no data for tennis balls being presented first for ASC suggesting that these results exist within SCOs only

###Does latency differ depending on the order puzzles were presented, puzzle type and species
speciespuzzleordertype<-aov(latency~puzzle.type:puzzle.order:species,data=otterpuzzledatav2)
TukeyHSD(speciespuzzleordertype) #significant difference between tennis balls presented first vs. third, and between tennis balls and bottles in SCO only

#Does latency differ between sexes depending on species
sexspecieslatency<-aov(latency~species:sex,data=otterpuzzledatav2)
TukeyHSD(sexspecieslatency) #significant difference lies between ASC M/F and SCO females only.
#There was a SCO female that was presented tennis balls first and had the longest latency of all data points (latency=90 seconds)
#The next longest time to first interact was 54 seconds, also in a female SCO presented with tennis balls first.

##Does time spent interacting with puzzles correlate with juggle rate, species, the interaction between juggle rate and species, the type of puzzle...
##...the order in which puzzles were presented, the interaction between puzzle type and puzzle order, sex, age, or site?
puzzleinteraction<-lmer(interaction~jugglerate*species+puzzle.type*puzzle.order+sex+age+site+(1|group.ID/ID), data=otterpuzzledatav2) 
summary(puzzleinteraction) 

puzzleinteraction2<-lmer(interaction~jugglerate+species+puzzle.type*puzzle.order+sex+age+site+(1|group.ID/ID),data=otterpuzzledatav2)
anova(puzzleinteraction,puzzleinteraction2) #jugglerate:species interaction does not significantly improve model fit (X1=0.018, p=0.894), so it is removed from the model
summary(puzzleinteraction2)

puzzleinteraction3<-lmer(interaction~jugglerate+species+puzzle.type*puzzle.order+age+site+(1|group.ID/ID),data=otterpuzzledatav2)
anova(puzzleinteraction2,puzzleinteraction3) #sex does not significantly improve model fit (X1=0.006, p=0.940), so it is removed from the model
summary(puzzleinteraction3) 

puzzleinteraction4<-lmer(interaction~jugglerate+puzzle.type*puzzle.order+age+site+(1|group.ID/ID),data=otterpuzzledatav2)
anova(puzzleinteraction3,puzzleinteraction4)#species does not significantly improve model fit (X1=0.035, p=0.851), so it is removed from the model
summary(puzzleinteraction4)

puzzleinteraction5<-lmer(interaction~jugglerate+puzzle.type*puzzle.order+age+(1|group.ID/ID),data=otterpuzzledatav2)
anova(puzzleinteraction4,puzzleinteraction5) #site does not significantly improve model fit (X2=2.869, p=0.238)
summary(puzzleinteraction5) #MAM

puzzleinteraction6<-lmer(interaction~jugglerate+puzzle.order+puzzle.type+age+(1|group.ID/ID),data=otterpuzzledatav2)
anova(puzzleinteraction5,puzzleinteraction6) #puzzle.type:puzzle.order interaction is highly significant (X4=25.244, p<0.001)

#Does interaction time differ with puzzle order for each puzzle type between species
interactionspeciesordertype<-aov(interaction~species:puzzle.type:puzzle.order,data=otterpuzzledatav2)
TukeyHSD(interactionspeciesordertype)
#No significant difference between the order in which puzzles were presented for each puzzle type in ASCs
#No significant difference between puzzle types that were presented first vs. second, first vs. third, or second vs. third in ASCs
#No significant difference between the order in which puzzles were presented for each puzzle type in SCOs
#No significant difference between puzzle types that were presented first vs. second, first vs. third, or second vs. third in SCOs
#Significant results tend to lie between SCO:ball:first and other factors.
#e.g. ASC:bottle:first - SCO:ball:first, p<0.008; ASC:ball:second - SCO:ball:first, p<0.001

##Does success solving puzzles correlate with juggle rate, species, the interaction between juggle rate and species, the type of puzzle...
##...the order in which puzzles were presented, the interaction between puzzle type and puzzle order, sex, age, or site?
puzzlesolve<-glmer(solve~jugglerate*species+puzzle.type*puzzle.order+sex+age+site+(1|group.ID/ID),family=binomial,data=otterpuzzledatav2) #failed to converge
summary(puzzlesolve)#remove sex

puzzlesolve1<-glmer(solve~jugglerate*species+puzzle.type*puzzle.order+age+site+(1|group.ID/ID),family=binomial,data=otterpuzzledatav2) #failed to converge
summary(puzzlesolve1)#remove puzzle type:puzzle order interaction

puzzlesolve2<-glmer(solve~jugglerate*species+puzzle.type+puzzle.order+age+site+(1|group.ID/ID),family=binomial,data=otterpuzzledatav2) #fails to converge
summary(puzzlesolve2)#remove puzzle.order
anova(puzzlesolve1,puzzlesolve2) #puzzle order does not appear to significantly improve model fit (X4=8.986, p-0.061)

puzzlesolve3<-glmer(solve~jugglerate*species+puzzle.type+age+site+(1|group.ID/ID),family=binomial,data=otterpuzzledatav2) #failed to converge
summary(puzzlesolve3)#remove site

puzzlesolve4<-glmer(solve~jugglerate*species+puzzle.type+age+(1|group.ID/ID),family=binomial,data=otterpuzzledatav2) #failed to converge
summary(puzzlesolve4) #remove age

puzzlesolve5<-glmer(solve~jugglerate*species+puzzle.type+(1|group.ID/ID),family=binomial,data=otterpuzzledatav2) #converges
summary(puzzlesolve5) #remove juggle rate

puzzlesolve6<-glmer(solve~jugglerate:species+species+puzzle.type+(1|group.ID),family=binomial,data=otterpuzzledatav2)
anova(puzzlesolve5,puzzlesolve6) #juggle rate does not significantly improve model fit (X1=0.756, p=0.384), so it is removed
summary(puzzlesolve6)

puzzlesolve7<-glmer(solve~jugglerate:species+species+(1|group.ID/ID),family=binomial,data=otterpuzzledatav2)
anova(puzzlesolve6,puzzlesolve7) #puzzle.type significantly improves model fit so it is retained, remove jugglerate:species instead

puzzlesolve8<-glmer(solve~species+puzzle.type+(1|group.ID/ID),family=binomial,data=otterpuzzledatav2) #fails to converge, remove predictors in turn and compare to null model
summary(puzzlesolve8)

puzzlesolve9<-glmer(solve~species+(1|group.ID),family=binomial,data=otterpuzzledatav2)
summary(puzzlesolve9)

puzzlenullsolve<-glmer(solve~1+(1|group.ID/ID),family=binomial,data=otterpuzzledatav2)
anova(puzzlesolve9,puzzlenullsolve) #species does not appear significant. However, this anova reports no dfs. Remove species instead

puzzlesolve10<-glmer(solve~puzzle.type+(1|group.ID/ID),family=binomial,data=otterpuzzledatav2)
anova(puzzlesolve10,puzzlenullsolve) #puzzle type appears significant (X2=6.992, p=0.030)
#This model seems to have major convergence issues, this may be due to the model being overfitted.

###Does the time spent interacting with puzzles before succesfully solving them correlate with juggle rate, species, the interaction between juggle rate and species, the type of puzzle...
##...the order in which puzzles were presented, the interaction between puzzle type and puzzle order, sex, age, or site?
solveinteraction<-lmer(interaction~jugglerate*species+puzzle.type*puzzle.order+sex+age+site+(1|group.ID/ID), data=ottersolvedata)
summary(solveinteraction)

solveinteraction2<-lmer(interaction~jugglerate*species+puzzle.type*puzzle.order+sex+site+(1|group.ID/ID),data=ottersolvedata)
anova(solveinteraction,solveinteraction2) #age does not significantly improve model fit (X1=0.001, p=0.972), so it is removed
summary(solveinteraction2) 

solveinteraction3<-lmer(interaction~jugglerate+species+puzzle.type*puzzle.order+sex+site+(1|group.ID/ID),data=ottersolvedata) #fails to converge
anova(solveinteraction2,solveinteraction3) #jugglerate:species interaction does not appear to be significant (X1=0.183, p=0.669), so it is removed
summary(solveinteraction3)

solveinteraction4<-lmer(interaction~jugglerate+species+puzzle.type*puzzle.order+site+(1|group.ID/ID),data=ottersolvedata) #model converges
anova(solveinteraction3,solveinteraction4) #sex does not significantly improve model fit (X1=0.302, p=0.583), so it is removed
summary(solveinteraction4) 

solveinteraction5<-lmer(interaction~jugglerate+species+puzzle.type*puzzle.order+(1|group.ID/ID),data=ottersolvedata)
anova(solveinteraction4,solveinteraction5) #site does not significantly improve model fit (X2=5.701, p=0.058), so it is removed
summary(solveinteraction5)

solveinteraction6<-lmer(interaction~jugglerate+puzzle.type*puzzle.order+(1|group.ID/ID),data=ottersolvedata)
anova(solveinteraction5,solveinteraction6)#species does not significantly improve model fit (X1=1.353, p=0.245), so it is removed
summary(solveinteraction6)

solveinteraction7<-lmer(interaction~puzzle.type*puzzle.order+(1|group.ID/ID),data=ottersolvedata) #fails to converge
anova(solveinteraction6,solveinteraction7) #juggle rate does not appear to be significant (X1=0.372, p=0.542)
summary(solveinteraction7) #all predictors appear significant, try removing stepwise anyway

solveinteraction8<-lmer(interaction~puzzle.type+puzzle.order+(1|group.ID/ID),data=ottersolvedata) #converges
anova(solveinteraction7,solveinteraction8) #puzzle type:puzzle order interaction appears to be significant (X4=77.097, p<0.001)

solveinteraction9<-lmer(interaction~puzzle.type:puzzle.order+puzzle.type+(1|group.ID/ID),data=ottersolvedata) #fails to converge
anova(solveinteraction7,solveinteraction9) #removing puzzle order, fails to converge and anova looks suspect

solveinteraction10<-lmer(interaction~puzzle.type:puzzle.order+puzzle.order+(1|group.ID/ID),data=ottersolvedata) #fails to converge
anova(solveinteraction10,solveinteraction7) #removing puzzle type, fails to converge and anova looks suspect.
####Again, there are convergence issues with this model

###Does time spent interacting with puzzles to successfully solve them differ between puzzle types depending on puzzle order between species?
speciessolveinteracttypeorder<-aov(interaction~puzzle.type:puzzle.order:species,data=ottersolvedata)
TukeyHSD(speciessolveinteracttypeorder) #significant results lie in SCOs only
#A significant difference between tennis balls presented first vs. third
#A significant difference between tennis balls vs. bottles that were presented first.
#This is likely due to a result within SCOs, SCO cub 1 was presented the tennis ball first and took 819 seconds to solve, 
#the next longest interaction time to solve is 378 seconds.

##Puzzle data analyses (without "puzzle order" included)

## Does latency to first interact with puzzles correlate with juggle rate, species, the interaction between juggle rate and species, age and wildlife park (site)?
latency<-lmer(latency~jugglerate*species+puzzle.type+sex+age+site+(1|group.ID/ID),data=otterpuzzledatav2)
summary(latency)

latency2<-lmer(latency~jugglerate+species+puzzle.type+sex+age+site+(1|group.ID/ID),data=otterpuzzledatav2)
summary(latency2)
anova(latency,latency2) # jugglerate:species does not significantly improve model fit (Chisq=0.370, p=0.543), so it is removed from the model

latency3<-lmer(latency~jugglerate+species+puzzle.type+sex+age+(1|group.ID/ID),data=otterpuzzledatav2)
summary(latency3)
anova(latency2,latency3) # site does not signficantly improve model fit (Chisq=1.851, p=0.396), so it is removed from the model 

latency4<-lmer(latency~species+puzzle.type+sex+age+(1|group.ID/ID),data=otterpuzzledatav2)
summary(latency4)
anova(latency3,latency4) # jugglerate does not signficantly improve model fit (Chisq=0.060, p=0.806), so it is removed from the model

latency5<-lmer(latency~species+puzzle.type+sex+(1|group.ID/ID),data=otterpuzzledatav2)
summary(latency5)
anova(latency4,latency5) # age does not significantly improve model fit (Chisq=0.744, p=0.388), so it is removed from the model

latency6<-lmer(latency~species+puzzle.type+(1|group.ID/ID),data=otterpuzzledatav2) ###MAM
summary(latency6)
anova(latency5, latency6) # sex does not significantly improve model fit (Chisq=2.624, p=0.105), so it is removed from the model

specieslatency<-update(latency6,~.-species)
puzzletypelatency<-update(latency6,~.-puzzle.type)
anova(specieslatency,latency6) #species significantly improves model fit (Chisq=8.481, p=0.003), so it is retained

anova(puzzletypelatency,latency6) # puzzle type significantly improves model fit (Chisq=6.973, p=0.031), so it is retained

summary(latency6) #SCOs had a longer latency than ASC. Overall, latency for puzzle balls was longest compared to bottles and bricks, bricks had shortest latency
#Difference in latency between puzzle types is likely due to species difference in latency as suggested by Figure 5

##Does latency differ between species for each puzzle type
speciesandpuzzle<-aov(latency~species:puzzle.type,data=otterpuzzledatav2)
TukeyHSD(speciesandpuzzle)
#Tennis balls - diff=31.213, lwr=20.425, upr=42.002, p=<0.001
#Bottles - diff=-0.872, lwr=-12.200, upr=10.454, p=1
#Bricks - diff=-0.815, lwr=-11.526, upr=9.896, p=1

##Does latency to first interact with puzzles in Asian small-clawed otters (ASCs) correlate with juggle rate, puzzle type, sex, age or site?
asclatency<-lmer(latency~jugglerate+puzzle.type+sex+age+site+(1|group.ID/ID),data=ascdata)
summary(asclatency)

asclatency2<-lmer(latency~jugglerate+sex+age+site+(1|group.ID/ID),data=ascdata)
summary(asclatency2)
anova(asclatency,asclatency2) #puzzle type does not significantly improve model fit (X2=0.405, p=0.817), so it is removed from the model

asclatency3<-lmer(latency~sex+age+site+(1|group.ID/ID),data=ascdata)
summary(asclatency3)
anova(asclatency2,asclatency3) #juggle rate does not significantly improve model fit (X1=0.243, p=0.622), so it is removed from the model
asclatency4<-lmer(latency~age+site+(1|group.ID/ID),data=ascdata)
summary(asclatency4)
anova(asclatency3,asclatency4) #sex does not significantly improve model fit (X1=0.486, p=0.486), so it is removed from the model

asclatency5<-lmer(latency~age+(1|group.ID/ID),data=ascdata)
summary(asclatency5)
anova(asclatency4,asclatency5) #site does not significantly improve model fit (X2=4.479, p=0.107), so it is removed from the model

asclatency6<-lmer(latency~1+(1|group.ID/ID),data=ascdata) 
summary(asclatency6)
anova(asclatency5,asclatency6) #age does not significantly improve model fit (X1=2.490, p=0.115), so it is removed from the model

##Does latency to first interact with puzzles in smooth-coated otters (SCOs) correlate with juggle rate, puzzle type, sex or age (all SCOs were at the same site)
scolatency<-lmer(latency~jugglerate+puzzle.type+sex+age+(1|group.ID/ID),data=scodata)
summary(scolatency)

scolatency2<-lmer(latency~puzzle.type+sex+age+(1|group.ID/ID),data=scodata)
summary(scolatency2)
anova(scolatency,scolatency2) #juggle rate does not significantly improve model fit (X1=2.326, p=0.127), so it is removed from the model

scolatency3<-lmer(latency~puzzle.type+sex+(1|group.ID/ID),data=scodata)
summary(scolatency3)
anova(scolatency2,scolatency3) #age does not significantly improve model fit (X1=0.758, p=0.384), so it is removed from the model

scolatency4<-lmer(latency~puzzle.type+(1|group.ID/ID),data=scodata) #MAM
summary(scolatency4)
anova(scolatency3,scolatency4) #sex does not significantly improve model fit (X1=2.167, p=0.141), so it is removed from the model

scolatency5<-lmer(latency~1+(1|group.ID/ID),data=scodata)
summary(scolatency5)
anova(scolatency4,scolatency5) #puzzle type significantly improves model fit (X2=9.372, p=0.009)

##Does time spent interacting with the puzzle correlate with juggle rate, species, the interaction between juggle rate and species, puzzle type, sex, age or wildlife park (site)?
interaction<-lmer(interaction~jugglerate*species+puzzle.type+sex+age+site+(1|group.ID/ID), data=otterpuzzledatav2) 
summary(interaction)

interaction2<-lmer(interaction~jugglerate*species+puzzle.type+sex+age+(1|group.ID/ID), data=otterpuzzledatav2) 
summary(interaction2)
anova(interaction, interaction2) #site does not significantly improve model fit (X2=0.350, p=0.839), so it is removed from the model

interaction3<-lmer(interaction~jugglerate*species+puzzle.type+age+(1|group.ID/ID), data=otterpuzzledatav2) 
summary(interaction3)
anova(interaction2, interaction3) #sex does not significantly improve model fit (X1=0.192, p=0.661), so it is removed from the model

interaction4<-lmer(interaction~jugglerate+species+puzzle.type+age+(1|group.ID/ID), data=otterpuzzledatav2) 
summary(interaction4)
anova(interaction3, interaction4) #jugglerate:species interaction does not significantly improve model fit (X1=1.903, p=0.168), so it is removed from the model

interaction5<-lmer(interaction~jugglerate+puzzle.type+age+(1|group.ID/ID), data=otterpuzzledatav2) ###MAM
summary(interaction5) #increased juggle rate = increased interaction time; increase in age = decreased interaction time

juggleinteraction<-update(interaction5,~.-jugglerate)
anova(interaction5,juggleinteraction) ##juggle rate significantly improves model fit (X1=4.253, p=0.039) 

puzzletypeinteraction<-update(interaction5,~.-puzzle.type)
anova(interaction5,puzzletypeinteraction) ##puzzle type significantly improves model fit (X2=14.969, p<0.001) 

ageinteraction<-update(interaction5,~.-age)
anova(interaction5,ageinteraction) ##age significantly improves model fit (X1=5.696, p=0.017)

#Does success at solving puzzle toys correlate with juggle rate, species, type of puzzle or age?
solve<-glmer(solve~jugglerate+species+puzzle.type+age+(1|group.ID/ID),family=binomial,data=otterpuzzledatav2) #converges
summary(solve)

solve2<-glmer(solve~species+puzzle.type+age+(1|group.ID/ID),family=binomial,data=otterpuzzledatav2)
anova(solve,solve2) #juggle rate does not significantly improve the model fit (X1=0.087, p=0.768), so it is removed
summary(solve2)

solve3<-glmer(solve~species+puzzle.type+(1|group.ID/ID),family=binomial,data=otterpuzzledatav2) #fails to converge when age is removed
anova(solve2,solve3) #including 'age' does not appear to significantly improve the model fit (X1=0.819, p=0.365), so it is removed
summary(solve3)

solve4<-glmer(solve~puzzle.type+(1|group.ID/ID),family=binomial,data=otterpuzzledatav2) ##Minimal adequate model
anova(solve3,solve4) #including 'species' does not significantly improve model fit (X1=3.305, p=0.069), so it is removed

nullsolve<-glmer(solve~1+(1|group.ID/ID),family=binomial,data=otterpuzzledatav2)
anova(nullsolve,solve4) #including puzzle type significantly improves model fit (X2=6.992, p=0.030)

##Does time spent interacting with the puzzles before successfully solving them correlate with juggle rate, species, the interaction between juggle rate and species, puzzle type, sex, age or wildlife park (site)?
solveinteraction<-lmer(interaction~jugglerate*species+puzzle.type+sex+age+site+(1|group.ID/ID), data=ottersolvedata)
summary(solveinteraction)

solveinteraction2<-lmer(interaction~jugglerate*species+puzzle.type+sex+age+(1|group.ID/ID),data=ottersolvedata)
anova(solveinteraction,solveinteraction2) #including 'site' does not significantly improve model fit (X2=0.209, p=0.901), so it is removed
summary(solveinteraction2)

solveinteraction3<-lmer(interaction~jugglerate*species+puzzle.type+age+(1|group.ID/ID),data=ottersolvedata)
anova(solveinteraction2,solveinteraction3) #including 'sex' does not significantly improve model fit (X1=0.016, p=0.899), so it is removed
summary(solveinteraction3)

solveinteraction4<-lmer(interaction~jugglerate*species+puzzle.type+(1|group.ID/ID),data=ottersolvedata)
anova(solveinteraction3,solveinteraction4) #including 'age' does not significantly improve model fit (X1=0.669, p=0.413), so it is removed
summary(solveinteraction4)

solveinteraction5<-lmer(interaction~jugglerate*species+(1|group.ID/ID),data=ottersolvedata)
anova(solveinteraction4,solveinteraction5) # including 'puzzle type' does not significantly improve model fit (X2=3.129, p=0.209), so it is removed
summary(solveinteraction5)

solveinteraction6<-lmer(interaction~jugglerate+species+(1|group.ID/ID),data=ottersolvedata) 
anova(solveinteraction5,solveinteraction6) #including 'juggle rate:species' interaction does not significantly improve model fit (X1=3.237, p=0.072), so it is removed
summary(solveinteraction6)

solveinteraction7<-lmer(interaction~species+(1|group.ID/ID),data=ottersolvedata)###MAM
anova(solveinteraction6,solveinteraction7) #including 'juggle rate' does not significantly improve model fit (X1=0.010, p=0.920), so it is removed
summary(solveinteraction7)

nullsolveinteraction<-lmer(interaction~1+(1|group.ID/ID),data=ottersolvedata)
anova(solveinteraction7,nullsolveinteraction) #including 'species significantly improves model fit (X1=8.852, p=0.003), so it is retained


###Plots
# Rock juggling ~ Species (boxplot)
ggplot(data=otterminbubble,
       aes(x=species, y=rj.total)) +
  geom_boxplot(fill="gray") +
  labs(x = "Species", y = "Total number of rock juggling observations") +
  theme_classic()

# Rock juggling ~ Age:Agecat interaction
library(ggplot2)
ggplot(data=otterminbubble) +
  aes(x = age, y = rj.total, color = species, shape=agecat) +
  geom_point(size=2.5) +
  geom_smooth(method="glm.nb",aes(fill=species))+
  scale_color_manual(labels=c("ASC","SCO"),values=c("grey60","grey30"))+
  scale_shape_manual(labels=c("Old (>11 years)","Young (≤11 years)"), values=c(15,16)) +
  scale_fill_manual(labels=c("ASC","SCO"),values=c("grey60","grey30"))+
  guides(color=guide_legend("Species"),shape=guide_legend("Age Category"),fill=guide_legend("Species")) +
  labs(x = "Age (Years)", y = "Total number of rock juggling observations") +
  theme_classic()

# Rock juggling ~ Hunger
hungerorder<-c("H","F")
ggplot(data=otterhungrydata,
       aes(x=hungry, y=rj)) +
  geom_boxplot(fill="gray") +
  scale_x_discrete(limits=hungerorder,labels=c("Hungry","Full"))+
  labs(x = "Hunger", y = "Rock juggling frequency") +
  theme_classic()

# Latency to first interaction with puzzle ~ puzzle type for each species plot
ggplot(data=otterpuzzledatav2,
       aes(x=species, y=latency, fill=puzzle.type)) +
  geom_boxplot(position=position_dodge(0.8)) +
  labs(x = "Species", y = "Latency to first interaction (seconds)") + 
  guides(fill=guide_legend("Puzzle Type")) + 
  scale_fill_manual(labels=c("Ball","Bottle","Brick"),values=c("grey45", "grey75","grey90")) +
  theme_classic()

###Descriptive statistics
##Medians and IQR##
#Age Cat
describe(otterminbubble$rj.total[otterminbubble$agecat=="YOUNG"])
median(otterminbubble$rj.total[otterminbubble$agecat=="YOUNG"])
quantile(otterminbubble$rj.total[otterminbubble$agecat=="YOUNG"])
describe(otterminbubble$rj.total[otterminbubble$agecat=="OLD"])
quantile(otterminbubble$rj.total[otterminbubble$agecat=="OLD"])

#Species
describe(otterminbubble$rj.total[otterminbubble$species=="SCO"])
quantile(otterminbubble$rj.total[otterminbubble$species=="SCO"])
describe(otterminbubble$rj.total[otterminbubble$species=="ASC"])
quantile(otterminbubble$rj.total[otterminbubble$species=="ASC"])

#Sex - SCO
describe(otterminbubble$rj.total[otterminbubble$sex=="M"&otterminbubble$species=="SCO"])
quantile(otterminbubble$rj.total[otterminbubble$sex=="M"&otterminbubble$species=="SCO"])
describe(otterminbubble$rj.total[otterminbubble$sex=="F"&otterminbubble$species=="SCO"])
quantile(otterminbubble$rj.total[otterminbubble$sex=="F"&otterminbubble$species=="SCO"])

#Sex - ASC
describe(otterminbubble$rj.total[otterminbubble$sex=="M"&otterminbubble$species=="ASC"])
quantile(otterminbubble$rj.total[otterminbubble$sex=="M"&otterminbubble$species=="ASC"])
describe(otterminbubble$rj.total[otterminbubble$sex=="F"&otterminbubble$species=="ASC"])
quantile(otterminbubble$rj.total[otterminbubble$sex=="F"&otterminbubble$species=="ASC"])

#Hunger
median(otterhungrydata$rj[otterhungrydata$hungry=="H"])
quantile(otterhungrydata$rj[otterhungrydata$hungry=="H"])
median(otterhungrydata$rj[otterhungrydata$hungry=="F"])
quantile(otterhungrydata$rj[otterhungrydata$hungry=="F"])

#Latency means and SE for each puzzle type in SCO and ASC
describe(scodata$latency[scodata$puzzle.type=="ball"])
describe(scodata$latency[scodata$puzzle.type=="bottle"])
describe(scodata$latency[scodata$puzzle.type=="brick"])
describe(ascdata$latency[ascdata$puzzle.type=="ball"])
describe(ascdata$latency[ascdata$puzzle.type=="bottle"])
describe(ascdata$latency[ascdata$puzzle.type=="brick"])

#Interaction time means and SE for each puzzle type
describe(otterpuzzledatav2$interaction[otterpuzzledatav2$puzzle.type=="ball"]) 
describe(otterpuzzledatav2$interaction[otterpuzzledatav2$puzzle.type=="bottle"]) 
describe(otterpuzzledatav2$interaction[otterpuzzledatav2$puzzle.type=="brick"]) 

#Interaction time to successfully solve puzzles means and SEs for each species
describe(ottersolvedata$interaction[ottersolvedata$species=="SCO"])
describe(ottersolvedata$interaction[ottersolvedata$species=="ASC"])


