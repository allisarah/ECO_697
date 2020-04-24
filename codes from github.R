
# codes from github -------------------------------------------------------


# test for normality ------------------------------------------------------

################################################################
# This script reads, explores, and manipulates  Euorpean_NO2_LUR
# data
################################################################

## load the NO2 data
nitrogen <- read.csv(file = "Euorpean_NO2_LUR.csv",
                     header = T,
                     stringsAsFactors = F)

## explore NO2 data
class(nitrogen)
dim(nitrogen)
names(nitrogen)
head(nitrogen)
tail(nitrogen)
summary(nitrogen)
table(nitrogen$Country)
table(nitrogen$Station.Type)
str(nitrogen)

### replace <LOD values with 5/sqrt(2)
nitrogen$NO2 <- replace(nitrogen$NO2,
                        is.na(nitrogen$NO2),
                        5/sqrt(2))

### create the Over30 variable
nitrogen$Over30 <- ifelse(nitrogen$NO2 <= 30, 0, 1)

### start modeling building by creating the data subset
high.NO2 <- subset(nitrogen, Over30 == 1)

### create a log$NO2 variable in high.NO2
high.NO2$LogNO2 <- log(high.NO2$NO2)

### reduce high.no2 to variables that will be used in
### the analyses, and complete cases.
high.NO2 <- subset(high.NO2,
                   select = c("LogNO2", "Agriculture", "Road.Med", "Temperature",
                              "Transport.Emiss","Industrial.Emiss", "Urban.Emiss"))

### get rid of incomplete rows
high.NO2 <- subset(high.NO2, complete.cases(high.NO2) == T)

# make a new object to ckeck no2 distribution
NO2 <- high.NO2$NO2

#   check summary
summary(NO2)

#   calculate Arithmetic mean, standard deviation, median,
#   95th percentile of the non-transformed values.
mean(NO2)
sd(NO2)
median(NO2)
quantile(NO2,0.95)

### plot a histogram of the non-transformed values
#   base histogram
hist(NO2)

#   base histogram with 50 bins and better labels, in blue
hist(NO2,
     breaks = 50,
     main = "NO2 Concentrations",
     ylab = "Frequency",
     xlab = "NO2 Concentration (μg/m3)",
     col = "blue")

### create a logno2 variable and test for normality
LogNO2 <- log(high.NO2$NO2)

### plot a QQ plot of the log-transformed values
#   base QQ-plot comparing log-transformed values with theoretical normal
qqnorm(LogNO2)

#   add a thick red 1:1 line to the base plot
qqline(LogNO2,
       col = "red",
       lwd = 3)

#   update labels
qqnorm(LogNO2,
       main = "Log-Transformed NO2 Concentrations",
       ylab = "Quantiles of the Nitrogen Dioxide Sample",
       xlab = "Quantiles of the Normal Distribution")
qqline(LogNO2,
       col = "red",
       lwd = 3)
### conduct a Shapiro-Wilks test on the non-transformed values and
### log-transformed values
shapiro.test(NO2)
shapiro.test(LogNO2)

### association between dependent and independent variables
### and linear regression modeling with individual variables
with(high.NO2, cor.test(LogNO2, Agriculture))
with(high.NO2, t.test(LogNO2 ~ Road.Med))
with(high.NO2, cor.test(LogNO2, Temperature))
with(high.NO2, cor.test(LogNO2, Transport.Emiss))
with(high.NO2, cor.test(LogNO2, Industrial.Emiss))
with(high.NO2, cor.test(LogNO2, Urban.Emiss))

fit1 <- lm(LogNO2 ~ Agriculture, data = high.NO2); summary(fit1); exp(fit1$coef)
fit2 <- lm(LogNO2 ~ Road.Med, data = high.NO2); summary(fit2); exp(fit2$coef)
fit3 <- lm(LogNO2 ~ Temperature, data = high.NO2); summary(fit3); exp(fit3$coef)
fit4 <- lm(LogNO2 ~ Transport.Emiss, data = high.NO2); summary(fit4); exp(fit4$coef)
fit5 <- lm(LogNO2 ~ Industrial.Emiss, data = high.NO2); summary(fit5); exp(fit5$coef)
fit6 <- lm(LogNO2 ~ Urban.Emiss, data = high.NO2); summary(fit6); exp(fit6$coef)

### association of independent variables
with(high.NO2, t.test(Agriculture ~ Road.Med))
with(high.NO2, cor.test(Agriculture , Temperature))
with(high.NO2, cor.test(Agriculture , Transport.Emiss))
with(high.NO2, cor.test(Agriculture , Industrial.Emiss))
with(high.NO2, cor.test(Agriculture , Urban.Emiss))
with(high.NO2, t.test(Temperature ~ Road.Med))
with(high.NO2, t.test(Transport.Emiss ~ Road.Med))
with(high.NO2, t.test(Industrial.Emiss ~ Road.Med))
with(high.NO2, t.test(Urban.Emiss ~ Road.Med))
with(high.NO2, cor.test(Temperature , Transport.Emiss))
with(high.NO2, cor.test(Temperature , Industrial.Emiss))
with(high.NO2, cor.test(Temperature , Urban.Emiss))
with(high.NO2, cor.test(Transport.Emiss , Industrial.Emiss))
with(high.NO2, cor.test(Transport.Emiss , Urban.Emiss))
with(high.NO2, cor.test(Industrial.Emiss , Urban.Emiss))

### pair-wise visualisation with scatterplots
pairs(high.NO2, gap = 0.2)

### model building with backward elimination
fit7 <- lm(LogNO2 ~ ., data = high.NO2)
summary(fit7); exp(fit7$coef)
fit8 <- lm(LogNO2 ~ Agriculture + Road.Med + Transport.Emiss + Industrial.Emiss + Urban.Emiss, data = high.NO2)
summary(fit8); exp(fit8$coef)
fit9 <- lm(LogNO2 ~ Agriculture + Road.Med + Transport.Emiss + Industrial.Emiss, data = high.NO2)
summary(fit9); exp(fit9$coef)
fit10 <- lm(LogNO2 ~ Agriculture + Road.Med + Transport.Emiss, data = high.NO2)
summary(fit10); exp(fit10$coef)
fit11 <- lm(LogNO2 ~ Road.Med + Transport.Emiss + Industrial.Emiss, data = high.NO2)
summary(fit11); exp(fit11$coef)
fit12 <- lm(LogNO2 ~ Road.Med + Transport.Emiss, data = high.NO2)
summary(fit12); exp(fit12$coef)

### comparison of fitted models
anova(fit7,fit8)
anova(fit7,fit9)
anova(fit7,fit10)
anova(fit7,fit11)
anova(fit8,fit9)
anova(fit8,fit10)
anova(fit8,fit11)
anova(fit9,fit10)
anova(fit9,fit11)
anova(fit10,fit11)
anova(fit11,fit12)

### stepwise regression for comparison with selected model (fit9)
fit13 <- step(fit7, direction = "both", trace = F) ; summary(fit13); exp(fit13$coef)

### comapre best fitted model from backward elimination with stepwise regression model
anova(fit11,fit13)

### use fit11 to make a low outcome prediction
exp(predict(fit11, data.frame("Road.Med" = 0,
                              "Transport.Emiss" = 5000,
                              "Industrial.Emiss" = 750),
            interval = "confidence"))

### use fit11 to make a high outcome prediction
exp(predict(fit11, data.frame("Road.Med" = 1,
                              "Transport.Emiss" = 250000,
                              "Industrial.Emiss" = 750000),
            interval = "confidence"))

### leaps regression to visualize inclusion of variables in model and change in R2
library(leaps)
x <- regsubsets(LogNO2 ~ Agriculture + Road.Med + Temperature + Transport.Emiss +
                  Industrial.Emiss + Urban.Emiss, data = high.NO2)
plot(x, scale = "adjr2")
