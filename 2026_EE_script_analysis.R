### INTRODUCTION #####

# Welcome. For the description of the project please visit: https://github.com/diekei/2025_ENS_diekei_dicliptera_speciation
# Article is available at: 


## LIBRARY ####

library(car)
library(emmeans)
library(multcomp)
library(lme4)

## HOST PREFERENCE ####

host <- read.csv('data/2024_EE_data_host_preference.csv')
host
attach(host)
str(host)

# -determine which factor explain preferences
glm.host <- glm(cbind(mik, dic) ~ race * sex, data = host, family = quasibinomial)
summary(glm.host)
anova(glm.host, test = "Chisq")
Anova(glm.host, type="II", test = "Wald")

# -population comparison (Mrace)
glm.host.M <- glm(cbind(mik, dic) ~ pop * sex, subset = pop!="pat" & pop!="ran", data = host, family = quasibinomial)
summary(glm.host.M)
anova(glm.host.M, test = "Chisq")
Anova(glm.host.M, type="II", test = "Wald")

# -population comparison (Drace)
glm.host.D <- glm(cbind(mik, dic) ~ pop * sex, subset = pop!="pdl" & pop!="bdg" & pop!="tbr", data = host, family = quasibinomial)
summary(glm.host.D)
anova(glm.host.D, test = "Chisq")
Anova(glm.host.D, type="II", test = "Wald")

## NON CHOICE MATING ####

ncm <- read.csv('data/2024_EE_data_nonchoice_mating.csv')
ncm
attach(ncm)
str(ncm)

# -determine which factor explain mating attempts
glm.ncm.att <- glm(cbind(att, natt) ~ fem * mle, data = ncm, family = binomial)
summary(glm.ncm.att)
anova(glm.ncm.att, test = "Chisq")
Anova(glm.ncm.att, type="II", test = "Wald")

# -pairwise comparison (mating attempts)
glm.att <- glm(cbind(att, natt) ~ pair, data = ncm, family = binomial)
summary(glm.att)
anova(glm.att, test = "Chisq")
Anova(glm.att, type="II", test = "Wald")
#summary(glht(glm.att, emm(pairwise ~ pair)), test=adjusted(type="holm"))
summary(pairs(emmeans(glm.att, ~ pair)), adjust = "holm")

# -determine which factor explain mating success
glm.ncm.succ <- glm(cbind(succ, fail) ~ fem * mle, data = ncm, family = binomial)
summary(glm.ncm.succ)
anova(glm.ncm.succ, test = "Chisq")
Anova(glm.ncm.succ, type="II", test = "Wald")

# pairwise comparison (mating success)
glm.mate <- glm(cbind(succ, fail) ~ pair, data = ncm, family = binomial)
summary(glm.mate)
anova(glm.mate, test = "Chisq")
Anova(glm.mate, type="II", test = "Wald")
#summary(glht(glm.mate, emm(pairwise ~ pair)), test=adjusted(type="holm"))
summary(pairs(emmeans(glm.mate, ~ pair)), adjust = "holm")

## EGG ####

egg <- read.csv('data/2024_EE_data_egg.csv', na.strings = "na")
egg
attach(egg)
str(egg)

# -determine which factor explain number of eggs produced
glm.egp1 <- glmer.nb(prod ~ fem * mle + (1|fem:fam), data = egg, family = quasipoisson)
summary(glm.egp1)
anova(glm.egp1)
Anova(glm.egp1, type="II", test = "Chisq")

# -determine significance between pairs
glm.egp2 <- glmer.nb(prod ~ pair + (1|fem:fam), data = egg, family = quasipoisson)
summary(glm.egp2)
anova(glm.egp2)
Anova(glm.egp2, type="II", test = "Chisq")
#summary(glht(glm.egp2,emm(pairwise~pair)), test=adjusted(type="holm"))
summary(pairs(emmeans(glm.egp2, ~ pair)), adjust = "holm")

# -determine which factor explain number of eggs hatch
glm.egh1 <- glmer.nb(htc ~ fem * mle + (1|fem:fam), subset = prod>2, data = egg, family = quasipoisson)
summary(glm.egh1)
anova(glm.egh1)
Anova(glm.egh1, type="II", test = "Chisq")

# -determine significance between pairs
glm.egh2 <- glmer.nb(htc ~ pair + (1|fem:fam), subset = prod>2, data = egg, family = quasipoisson)
summary(glm.egh2)
anova(glm.egh2)
Anova(glm.egh2, type="II", test = "Chisq")
#summary(glht(glm.egh2,emm(pairwise~pair)), test=adjusted(type="holm"))
summary(pairs(emmeans(glm.egh2, ~ pair)), adjust = "holm")

# -determine significance between pairs (considering rate, using binomial)
glm.egh3 <- glm(cbind(htc, nhtc) ~ pair, subset = prod>2, data = egg, family = quasibinomial)
summary(glm.egh3)
anova(glm.egh3)
Anova(glm.egh3, type="II", test = "Wald")
#summary(glht(glm.egh3,emm(pairwise~pair)), test=adjusted(type="holm"))
summary(pairs(emmeans(glm.egh3, ~ pair)), adjust = "holm")

# -determine which factor determine hatching duration
glm.egd1 <- glmer.nb(dur ~ fem * mle + (1|fem:fam), subset = dur<25, data = egg, family = quasipoisson)
summary(glm.egd1)
anova(glm.egd1)
Anova(glm.egd1, type="II", test = "Chisq")

# -determine significance between pairs (hatching duration)
glm.egd2 <- glmer.nb(dur ~ pair + (1|fem:fam), subset = dur<25, data = egg, family = quasipoisson)
summary(glm.egd2)
anova(glm.egd2)
Anova(glm.egd2, type="II", test = "Chisq")
#summary(glht(glm.egd2,emm(pairwise~pair)), test=adjusted(type="holm"))
summary(pairs(emmeans(glm.egd2, ~ pair)), adjust = "holm")

## CHOICE MATING ####

cm <- read.csv('data/2024_EE_data_choice_mating.csv', na.strings = "na")
cm
attach(cm)
str(cm)

# -determine which factor explain mating inferences in pres/abs of host
glm.cm <- glm(cbind(par, hyb) ~ male * host, data = cm, family = binomial)
summary(glm.cm)
anova(glm.cm, test = "Chisq")
Anova(glm.cm, type="II", test = "Wald")

# -determine which factor explain mating inferences in different spatial setting
glm.field <- glm(cbind(par, hyb) ~ male + arr, data = cm, family = binomial)
summary(glm.field)
anova(glm.field, test = "Chisq")
Anova(glm.field, type="II", test = "Wald")


## SIGHTINGS, MIGRATION, AND DISPERSAL ####

fid <- read.csv('data/2024_EE_data_fidelity.csv')
fid
attach(fid)
str(fid)

# -determine which factor explain occurrences of both races
glm.sig <- glm(cbind(mik, dic) ~ race * (sex + arr), data = fid, family = binomial)
summary(glm.sig)
anova(glm.sig)
Anova(glm.sig, type="II", test = "Wald")

mig <- read.csv('data/2024_EE_data_migration.csv')
mig
attach(mig)
str(mig)

# -determine which factor explain migration frequencies of both races
glm.mig <- glm(cbind(mik, dic) ~ (race * sex) + arr, data = mig, family = binomial)
summary(glm.mig)
anova(glm.mig)
Anova(glm.mig, type="II", test = "Wald")

dis <- read.csv('data/2024_EE_data_dispersal.csv')
dis
attach(dis)
str(dis)

# -determine which factor explain dispersal of both races
glm.dis <- glm(dis ~ race * sex * arr, data = dis, family = quasipoisson)
summary(glm.dis)
anova(glm.dis, test = "Chisq")
Anova(glm.dis, type="II", test = "Wald")

## LARVAL PERFORMANCE ####

lp <- read.csv('data/2024_EE_data_larval_performance.csv')
lp
attach(lp)
str(lp)

# -determine which factor explain number of larvae survived
glm.lp1 <- glmer(cbind(succ, fail) ~ fem * mle * food + (1|fem:fam), data = lp, family = binomial)
summary(glm.lp1)
anova(glm.lp1)
Anova(glm.lp1, type="II", test = "Chisq")

# -determine significance between pairs (considering rate, using binomial)
# change the subset!
glm.lp2 <- glm(cbind(succ, fail) ~ pair, 
               subset = food!="dic" & stg!="suc" & stg!="adh", data = lp, family = quasibinomial)
summary(glm.lp2)
anova(glm.lp2, test = "Chisq")
Anova(glm.lp2, type="II", test = "Wald")
#summary(glht(glm.lp2,emm(pairwise~pair)), test=adjusted(type="holm"))
summary(pairs(emmeans(glm.lp2, ~ pair)), adjust = "holm")

# -determine significance between food (considering rate, using binomial)
# change the subset!
glm.lp3 <- glm(cbind(succ, fail) ~ food, 
               subset = pair!= "mm" & pair!="dm" & pair!="dd" & stg!="suc" & stg!="acc", data = lp, family = quasibinomial)
summary(glm.lp3)
anova(glm.lp3, test = "Chisq")
Anova(glm.lp3, type="II", test = "Wald")



