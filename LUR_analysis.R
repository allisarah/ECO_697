
hist(log(lur.dat1$concn)) 
library(ggpubr)
ggdensity(log(lur.dat1$concn))
ggqqplot(log(lur.dat1$concn))#draws the correlation between a given sample and the normal distributio
shapiro.test(log(lur.dat1$concn))

names(pm_dat)
lur.dat <- as.data.frame(pm_dat)
lur.dat1 <- lur.dat %>% select(-c("Start_date", "ID", "X", "End_date", "site_type",
           "Week", "Neighborhood", "geometry", "dis2rd_gdis", "dis2rd_rast",      
           "road.50.lnt", "ndvi")) #%>% 
              #replace(is.na(.),0)          
  
names(lur.dat1)
summary(lur.dat1)
covar <- lur.dat1[, 2:77]

mod1 <- glm(concn~., data = lur.dat1, family = gaussian)

mod1

####supervised forward linear regression
n.mod <- lm(concn ~ 1, data = lur.dat1)
n.mod
summary(n.mod)

mod1 <- lm(log(concn) ~ Temp, data = lur.dat1)
mod2 <- lm(concn ~ Temp, data = lur.dat1)
mod1
summary(mod1)
summary(mod2)
sum(resid(mod1)^2) #124.56
glance(mod1)

summary(mod1)
temp <-augment(n.mod)
sum(temp$.resid^2) #189.86
glance(mod1)

## check correlation
corrplot::corrplot(cor(lur.dat1))

#####supervised forward linear regression
#run regression of concn with each predictor
varnames <- names(covar)
model.fits <- vector(length (varnames), mode = "list")
names(model.fits) <- varnames

for(i in varnames){
    modelformula <-paste("concn ~", i)
  model.fits[[i]] <- lm(as.formula(modelformula), data = lur.dat1)
  #assign(x = paste("m", i, sep = "."), value = lm(as.formula(modelformula), data = lur.dat1))
}

## to check r2 and regression ssq
glance(model.fits[[14]])
summary(model.fits[[14]])

sum(augment(model.fits[[1]])$.resid^2)
sum(resid(model.fits[[1]])^2)

##adjusted r2 >= 0.01

#temp - 0.19
#rel..hum - 0.59
#wind.speed - 0.19
#road.50.lnt - 0.005
#ndvi.50 - 0.007
#length.50 - 0.006
#length.100 - 0.008
#length.200 - 0.003
#length.500 - 0.005
#length.1000 - 0.003
#length.2000 - 0.002
#ndvi.100 - 0.007
#ndvi.200 - 0.006
#ndvi.500 - 0.004
#ndvi.1000 - 0.004
#ndvi.2000 - 0.003
#ndvi.4000 - 0.001
#com.area.500 - 0.003
#com.area.1000 - 0.003
#com.area.2000 - 0.003
#com.area.1000 - 0.003
#pop.200 - 0.001
#pop.500 - 0.002
#pop.1000 - 0.002
#pop.2000 - 0.002
#build.count.1000 - 0.002
#build.count.2000 - 0.002
#dist_airport - 0.002
#bus.count.100 - 0.004
##bus.count.200 - 0.008
##bus.count.500 - 0.004
##bus.count.1000 - 0.003
##bus.stat.200 - 0.008
##bus.stat.500 - 0.005
##bus.stat.1000 - 0.004
##bus.stat.2000 - 0.002


#run regression of log(concn) with each predictor
varnames <- names(covar)
model.fits.log <- vector(length (varnames), mode = "list")
names(model.fits.log) <- varnames

for(i in varnames){
  modelformula <-paste("log(concn) ~", i)
  model.fits.log[[i]] <- lm(as.formula(modelformula), data = lur.dat1)
  #assign(x = paste("m", i, sep = "."), value = lm(as.formula(modelformula), data = lur.dat1))
}

## to check r2 and regression ssq
glance(model.fits[[19]])
summary(model.fits[[14]])

sum(augment(model.fits[[1]])$.resid^2)
sum(resid(model.fits[[1]])^2)

##adjusted r2 >= 0.01

#temp - 0.19
#rel..hum - 0.59
#wind.speed - 0.19
#road.50.lnt - 0.005
#ndvi.50 - 0.007
#length.50 - 0.006
#length.100 - 0.008
#length.200 - 0.003
#length.500 - 0.005
#length.1000 - 0.003
#length.2000 - 0.002
#ndvi.100 - 0.007
#ndvi.200 - 0.006
#ndvi.500 - 0.004
#ndvi.1000 - 0.004
#ndvi.2000 - 0.003
#ndvi.4000 - 0.001
#com.area.500 - 0.003
#com.area.1000 - 0.003
#com.area.2000 - 0.003
#com.area.1000 - 0.003
#pop.200 - 0.001
#pop.500 - 0.002
#pop.1000 - 0.002
#pop.2000 - 0.002
#build.count.1000 - 0.002
#build.count.2000 - 0.002
#dist_airport - 0.002
#bus.count.100 - 0.004
##bus.count.200 - 0.008
##bus.count.500 - 0.004
##bus.count.1000 - 0.003
##bus.stat.200 - 0.008
##bus.stat.500 - 0.005
##bus.stat.1000 - 0.004
##bus.stat.2000 - 0.002


#run regression of log(concn) with log of each predictor
var.name.log <- paste("log","(",names(covar),")")
model.logs <- vector(length (var.name.log), mode = "list")
names(model.logs) <- var.name.log

for(i in var.name.log){
  modelformula <-paste("log(concn) ~", i)
  model.logs[[i]] <- lm(as.formula(modelformula), data = lur.dat1, na.action = na.exclude)
}



#####works - but returns individual models = all 81############
covar.test <- names(covar)
for(i in covar.test){
  modelformula <-paste("concn ~", i)
  assign(x = paste("m", i, sep = "."), value = lm(as.formula(modelformula), data = lur.dat1))
}


########################################################################
models <- lapply(paste("log(concn)", names(lur.dat1)[-1], sep = "~"), formula)
res.models <- lapply(models, FUN = function(x) {summary(lm(formula = x, data = lur.dat1))})
names(res.models) <- paste("concn", names(lur.dat1)[-1], sep = "~")


model.test <- lapply(seq_along(res.models), function(i) {
  
  data.frame(model = names(res.models)[i],
             intercept = res.models[[i]]$coefficients[1],
             coef = res.models[[i]]$coefficients[2],
             r = res.models[[i]]$r.squared,
             stringsAsFactors = FALSE)
})

View(do.call(rbind, model.test))
purrr::map_df(res.models, broom::glance, .id = 'formula')


model.logs <- lapply(paste("log(concn)", "log(names(lur.dat1)[-1])", sep = "~"), formula)
res.log.model <- lapply(model.logs, FUN = function(x) {summary(lm(formula = x, data = lur.dat1))})
names(res.log.model) <- paste("concn", names(lur.dat1)[-1], sep = "~")

temp <- lm(log(concn)~ (1/dist.maj.rd)^2, data = lur.dat1)
glance(temp)
lur.dat1$dis

#lur.dat1$inv.dis.maj <- (1/(lur.dat1$dis.maj.rd))


lur.dat1$inv.maj.rd.sq <- (1/lur.dat1$dist.maj.rd)^2

lur.dat1$inv.maj.rd <- (1/lur.dat1$dist.maj.rd)

temp <- lm(concn~inv.maj.rd, data = lur.dat1)
glance(temp)

plot(lur.dat1$concn, lur.dat1$inv.maj.rd)

cor(log(lur.dat1$concn), log(lur.dat1$inv.maj.rd))

cor(log(lur.dat1$concn), log(1/lur.dat1$ndvi.50))
#################################################
#check correlation of concn with covariates
#cor(dat.pm1[,unlist(lapply(dat.pm1, is.numeric))])

names(dat.pm1)
lur.dat <- dat.pm1  %>% as.data.frame()%>% select(-geometry) %>%
          mutate(mon = months(Start_date, abbreviate = T)) %>% 
          mutate(season = if_else(mon %in% c("Dec", "Jan", "Feb"), "Har","Non-Har"))

class(lur.dat)
names(lur.dat)
str(lur.dat$season)
table(lur.dat$season)

lur.dat1 <- lur.dat %>% select(-c("mon", "Start_date", "ID", "X", "End_date", "site_type",
                             "Temp", "Rel..Hum.", "Wind.Speed", "Week", 
                            "Neighborhood"))
str(lur.dat1)
names(lur.dat1 )
lapply(lur.dat1 , is.numeric)

lur.dat1 $season <- factor(lur.dat1 $season, levels = c("Non-Har", "Har"))
#temp6$season <- as.factor(temp5$season)
levels(lur.dat1 $season)
names(lur.dat1 )
unique(lur.dat1 $season)

#run regression of log(concn) with each predictor - worked
varnames <- names(lur.dat1 )[-c(3, 145)]
model.fits <- vector(length (varnames), mode = "list")
names(model.fits) <- varnames

for(i in varnames){
  modelformula <-paste("log(concn) ~ season +", i)
  model.fits[[i]] <- lm(as.formula(modelformula), data = lur.dat1 )
}

#View(do.call(rbind, model.fits))

res.summary <- purrr::map_df(model.fits, broom::glance, .id = 'formula')
View(res.summary)

####regression without season#####
varnames1 <- names(lur.dat1 )[-c(3, 145)]
model.fits1 <- vector(length (varnames1), mode = "list")
names(model.fits1) <- varnames1

for(i in varnames1){
  modelformula1 <-paste("log(concn) ~", i)
  model.fits1[[i]] <- lm(as.formula(modelformula1), data = lur.dat1 )
}

res.summary1 <- purrr::map_df(model.fits1, broom::glance, .id = 'formula')
View(res.summary1)

###regression for fixed sites
unique(lur.dat$ID)
fix.site <- lur.dat %>% filter(ID %in% c("AD", "JT", "La", "Nima", "Taifa", "Ash",
                                         "EL", "N1West", "TMW", "UGH"))

unique(fix.site$ID)

fix.site1 <- fix.site %>% select(-c("mon", "Start_date", "ID", "X", "End_date", "site_type",
                                  "Temp", "Rel..Hum.", "Wind.Speed", "Week", 
                                  "Neighborhood"))
str(fix.site1)

####regression with fixed sites only#####
names(fix.site1)
varnames.fix <- names(fix.site1)[-c(3, 145)]
model.fixs <- vector(length (varnames.fix), mode = "list")
names(model.fixs) <- varnames.fix

for(i in varnames.fix){
  modelformula.fix <-paste("log(concn) ~", i)
  model.fixs[[i]] <- lm(as.formula(modelformula.fix), data = fix.site1)
}

res.summary1 <- purrr::map_df(model.fixs, broom::glance, .id = 'formula')
View(res.summary1)


####rotating sites
rot.site <- lur.dat %>% filter(!ID %in% c("AD", "JT", "La", "Nima", "Taifa", "Ash",
                                         "EL", "N1West", "TMW", "UGH"))

unique(rot.site$ID)

rot.site1 <- rot.site %>% select(-c("mon", "Start_date", "ID", "X", "End_date", "site_type",
                                    "Temp", "Rel..Hum.", "Wind.Speed", "Week", 
                                    "Neighborhood"))
str(rot.site1)

####regression with rotating sites only#####
names(rot.site1)
varnames.rot <- names(rot.site1)[-c(3, 145)]
model.rot <- vector(length (varnames.rot), mode = "list")
names(model.rot) <- varnames.rot

for(i in varnames.rot){
  modelformula.rot <-paste("log(concn) ~", i)
  model.rot[[i]] <- lm(as.formula(modelformula.rot), data = rot.site1)
}

res.summary.rot <- purrr::map_df(model.rot, broom::glance, .id = 'formula')
View(res.summary.rot)


#### regression for non-harmattan season
non.har <- lur.dat %>% filter(season %in% "Non-Har")
table(non.har$season)

non.har1 <- non.har %>% select(-c("mon", "Start_date", "ID", "X", "End_date", "site_type",
                                    "Temp", "Rel..Hum.", "Wind.Speed", "Week", 
                                    "Neighborhood", "season"))
str(non.har1)

####regression with non-harmattan season#####
names(non.har1)
varnames.non.har <- names(non.har1)[-c(3, 144)]
model.non.har <- vector(length (varnames.non.har), mode = "list")
names(model.non.har) <- varnames.non.har

for(i in varnames.non.har){
  modelformula.non.har <-paste("log(concn) ~", i)
  model.non.har[[i]] <- lm(as.formula(modelformula.non.har), data = non.har1)
}

res.summary.non.har<- purrr::map_df(model.non.har, broom::glance, .id = 'formula')
View(res.summary.non.har)

glance(lm(log(concn)~ndvi.50+maj.1000.length+dist.bus.stp+com.area.500+bus.count.1000+
            bus.count.200+roads.1000.length+build.count.1000+bus.stat.1000+other.area.500+
            dist.maj.rd+hd.area.50+pop.200+roads.50.length+com.area.100, data = non.har1))
###r^2 was 0.572 using the variables above, other variables didnt change the r^2 by 0.01

#### regression for non-harmattan season
har <- lur.dat %>% filter(season %in% "Har")
table(har$season)

har1 <- har %>% select(-c("mon", "Start_date", "ID", "X", "End_date", "site_type",
                                  "Temp", "Rel..Hum.", "Wind.Speed", "Week", 
                                  "Neighborhood", "season"))
str(har1)

####regression with harmattan season#####
names(har1)
varnames.har <- names(har1)[-c(3, 144)]
model.har <- vector(length (varnames.har), mode = "list")
names(model.har) <- varnames.har

for(i in varnames.har){
  modelformula.har <-paste("log(concn) ~", i)
  model.har[[i]] <- lm(as.formula(modelformula.har), data = har1)
}

res.summary.har<- purrr::map_df(model.har, broom::glance, .id = 'formula')
View(res.summary.har)



#######################
###regression for temporally adjusted data ---- use this for project report
names(fix.rot1)
fix.rot2 <- fix.rot1 %>% select(-c("mon", "Start_date", "ID", "X", "End_date", "site_type",
                                    "Temp", "Rel..Hum.", "Wind.Speed", "Week", 
                                    "Neighborhood", "concn", "season", "avgPM", "annual.avg", "week.cf",
                                   "week.cf2", "week.cf3", "week.cf4", "pm.cor", "pm.cor2", "pm.cor3", "pm.cor4"))
str(fix.rot2)

names(fix.rot2)
varnames.fix.rot <- names(fix.rot2)[-c(145)]
model.timeadj <- vector(length (varnames.fix.rot), mode = "list")
names(model.timeadj) <- varnames.fix.rot

for(i in varnames.fix.rot){
  modelformula.fix.rot <-paste("log(pm.cor1) ~", i)
  model.timeadj[[i]] <- lm(as.formula(modelformula.fix.rot), data = fix.rot2)
}

res.summary.time <- purrr::map_df(model.timeadj, broom::glance, .id = 'formula')
View(res.summary.time)


#################using glm
mod<- glm(log(concn)~season +roads.4000.length, data = temp6, family = "gaussian")#1
summary(mod)
augment(mod)

##################penalized regression#############
library(SIS)
p.mat<- temp6[, c(1:2, 4:145)]  #(predictors)
p.mat<- na.omit(p.mat)
Y = log(temp6[, 3])

seed = 998
Pm2.5.mod <- SIS(p.mat, Y, family='gaussian', seed=seed,tune = 'bic',
                 penalty='lasso', varISIS = 'vanilla', standardize = TRUE)
Pm2.5.mod$sis.ix0  ##SIS only
Pm2.5.mod$ix ##SIS+lasso

library(corrplot)
corrplot(temp[c(7, 14:88, 90)], method="number",type= "lower",insig = "blank", 
         number.cex = 0.6)

library(psych)
pairs.panels(temp[c(7, 13:22)])

temp3 <- corr.test(x = temp[c(7, 13:88)])
temp4 <- corr.test(x = temp[c(7, 14:88, 90)])
View(temp3)
View(temp3[["r"]])
View(temp4[["r"]])
str(temp[c(7, 13:88)])
##group by season - harmattan vs non-harmattan
dat.pm1$Start_date <- as.Date(dat.pm1$Start_date)
str(dat.pm1$Start_date)
temp <- dat.pm1 %>% mutate(mon = months(Start_date)) %>% filter(!mon %in% c("December" ,"January", "February"))
temp <- as.data.frame(temp) %>% select(-geometry) 
temp <- temp %>% mutate(logPM = log(concn))
class(temp)
glance(lm(log(temp$concn)~log(temp$maj.1000.length))) #- r^2 - 0.31
glance(lm(log(temp$concn)~temp$maj.1000.length)) #- r^2 - 0.24
glance(lm(log(temp$concn)~log(temp$maj.500.length)))  #- r^2 - 0.21
glance(lm(log(temp$concn)~log(temp$dist.maj.rd))) #- r^2 - 0.1
glance(lm(log(temp$concn)~log(temp$inv.maj.rd))) #- r^2 - 0.1
glance(lm(log(temp$concn)~temp$hd.prop.50))#0.09
glance(lm(log(temp$concn)~temp$hd.area.200))#- r^2 - 0.1
glance(lm(log(temp$concn)~temp$com.area.4000))#- r^2 - 0.15
glance(lm(log(temp$concn)~temp$ndvi.50))#- r^2 - 0.319


# 
#non-har
temp1 <- dat.pm1 %>% mutate(mon = months(Start_date)) %>% filter(mon %in% c("December" ,"January", "February"))
str(temp1$Start_date)
shapiro.test(log(temp1$concn))
glance(lm(log(temp1$concn)~log(temp1$dist.maj.rd)))
glance(lm(log(temp1$concn)~log(temp1$maj.1000.length)))
glance(lm(log(temp1$concn)~log(temp1$dist.med.rd)))





