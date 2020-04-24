temp <- lm(concn~Rel..Hum+Temp, data = )


#build.count.100 - 0.02, 3.55
# build.count.200 - 0.02, 3.55##adjusted r2 >= 0.01, and t value
# Temp - r^2 - 0.34, t - 16.13
# Rel..Hum. - 0.64, -29.6
# Wind.speed - 0.29, -14.5
# ndvi.50 - 0.07, -6.18
# ld.area.50 - 0.009 ,-2.35
# hd.area.50 - 0.01, 2.88
# com.area.50 - 0.009, 2.42
# ld.area.100 - 0.02 ,-2.95
# hd.area.100 - 0.02 ,3.275
# pop.50 - 0.02, 3.005
# build.50 - 0.02, 2.99
# dist.maj.rd - 0.07, -6.2
# length.50 - 0.02, 3.53
# length.100 - 0.04, 4.5
# length.200 - 0.02, 3.58
# length.500 - 0.04, 4.5
# length.1000 - 0.04, 4.5
# length.2000 - 0.03, 4.2
# length.4000 - 0.01, 2.89
# ndvi.100 - 0.07, -6.09
# ndvi.200 - 0.06, -5.79
# ndvi.500 - 0.05, -5.36
# ndvi.1000 - 0.06, -6.09
# ndvi.2000 - 0.04, -4.86
# ndvi.4000 - 0.03, -3.93
# ld.area.200 - 0.01, -2.56
# hd.area.200 - 0.02, 3.38
# hd.area.500 - 0.01, 2.73
# com.area.500 - 0.04, -4.41
# other.area.500 - 0.02, -3.24
# hd.area.1000 - 0.01, 2.73
# com.area.1000 - 0.04, 4.41
# other.area.1000 - 0.02, -3.24
# hd.area.2000 - 0.01, 2.73
# com.area.2000 - 0.04, 4.41
# other.area.1000 - 0.02, -3.24
# hd.area.4000 - 0.01, 2.73
# com.area.4000 - 0.04, 4.41
# other.area.4000 - 0.02, -3.24
# pop.100 - 0.02, 2.99
# pop.200 - 0.02, 3.31
# pop.500 - 0.02, 3.61
# pop.1000 - 0.02, 3.6
# pop.2000 - 0.02, 3.44
# build.count.50 - 0.02, 3.38
# build.count.500 - 0.02, 3.54
# build.count.1000 - 0.03, 3.84
# build.count.2000 - 0.02, 3.52
# build.count.4000 - 0.02, 2.94
# dist_airport - 0.01, -2.29
# bus.count.50 - 0.01, 2.15
# bus.count.100 - 0.02, 3.41
# bus.count.200 - 0.04, 4.76
# bus.count.500 - 0.03, 4.21
# bus.count.1000 - 0.03, 4.37
# bus.count.2000 - 0.01, 2.66
# bus.count.4000 - 0.02, 2.97
# bus.stat.200 - 0.02, 3.2
# bus.stat.500 - 0.02, 3.4
# bus.stat.1000 - 0.03, 3.85
# bus.stat.2000 - 0.02, 3.55
# bus.stat.4000 - 0.02, 2.94

#############################################################
shapiro.test(log(lur.dat1$hd.area.50))
hist(log(lur.dat1$other.area.50))

shapiro.test(lur.dat1$ndvi.50)
hist(lur.dat1$Rel..Hum.)


###########################################
#retain variables with r^2 >= 0.01 increase in R^2
View(temp6)
glance(lm(log(concn) ~ season, data = temp6))#0.54
glance(lm(log(concn) ~ season + maj.50.length, data = temp6)) #0.64
glance(lm(log(concn) ~ season + maj.50.length + maj.100.length, data = temp6)) #0.65
glance(lm(log(concn) ~ season + maj.50.length + maj.100.length + ndvi.50, data = temp6)) #0.65 - drop ndvi.50
glance(lm(log(concn) ~ season + maj.50.length + maj.100.length + ndvi.100, data = temp6)) #0.65 - drop ndvi.100
glance(lm(log(concn) ~ season + maj.50.length + maj.100.length + maj.500.length, data = temp6)) #0.65 - drop maj.500.length
glance(lm(log(concn) ~ season + maj.50.length + maj.100.length + ndvi.200, data = temp6)) #0.65 - drop ndvi.200
glance(lm(log(concn) ~ season + maj.50.length + maj.100.length + ndvi.500, data = temp6)) #0.65 - drop ndvi.500
glance(lm(log(concn) ~ season + maj.50.length + maj.100.length + ndvi.1000, data = temp6))#0.65 - drop ndvi.1000
glance(lm(log(concn) ~ season + maj.50.length + maj.100.length + maj.200.length, data = temp6)) # 0.66
glance(lm(log(concn) ~ season + maj.50.length + maj.100.length + maj.200.length + maj.1000.length,
          data = temp6)) #0.66 - drop maj.1000.length

glance(lm(log(concn) ~ season + maj.50.length + maj.100.length + maj.200.length + ndvi.2000
            , data = temp6)) # 0.66 - drop ndvi.2000

glance(lm(log(concn) ~ season + maj.50.length + maj.100.length + maj.200.length + ndvi.4000
          , data = temp6)) # 0.66 - drop ndvi.4000

glance(lm(log(concn) ~ season + maj.50.length + maj.100.length + maj.200.length + maj.2000.length
          , data = temp6)) # 0.66 - drop maj.2000.length

glance(lm(log(concn) ~ season + maj.50.length + maj.100.length + maj.200.length + dist.maj.rd
          , data = temp6)) # 0.66 - drop dist.maj.rd




temp <-  lm(log(concn) ~ Rel..Hum. + Temp, data = lur.dat1)#0.639 - drop temp
temp <-  lm(log(concn) ~ Rel..Hum. + Wind.Speed, data = lur.dat1)#0.641 - drop Wind.Speed

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50, data = lur.dat1)#0.689
temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd, data = lur.dat1) #0.706
temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + ndvi.100, data = lur.dat1) #0.707 - drop ndvi.100
temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + ndvi.200, 
            data = lur.dat1) #0.71 - drop ndvi.200

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + ndvi.1000, 
            data = lur.dat1) #0.708 - drop ndvi.1000
temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + ndvi.500, 
            data = lur.dat1) #0.710 - drop ndvi.500
temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100, 
            data = lur.dat1) #0.719

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + length.500, 
            data = lur.dat1) #0.723 - drop length.500

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + length.1000, 
            data = lur.dat1) #0.724 - drop length.1000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + ndvi.2000, 
            data = lur.dat1) #0.721 - drop ndvi.2000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + com.area.500, 
            data = lur.dat1) #0.721 - drop com.area.500

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + com.area.1000, 
            data = lur.dat1) #0.721 - drop com.area.1000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + com.area.1000, 
            data = lur.dat1) #0.721 - drop com.area.1000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + com.area.2000, 
            data = lur.dat1) #0.721 - drop com.area.2000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + com.area.4000, 
            data = lur.dat1) #0.721 - drop com.area.1000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + bus.count.200, 
            data = lur.dat1) #0.720 - drop bus.count.200

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + length.2000, 
            data = lur.dat1) #0.722 - drop length.2000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + ndvi.4000, 
            data = lur.dat1) #0.721 - drop length.2000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + length.2000, 
            data = lur.dat1) #0.722 - drop length.2000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + build.count.1000, 
            data = lur.dat1) #0.719 - drop build.count.1000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + bus.count.500, 
            data = lur.dat1) #0.722 - drop bus.count.500

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + bus.count.1000, 
            data = lur.dat1) #0.721 - drop bus.count.1000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + bus.stat.1000, 
            data = lur.dat1) #0.720 - drop bus.stat.1000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + ld.area.100, 
            data = lur.dat1) #0.719 - drop ld.area.100

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + hd.area.100, 
            data = lur.dat1) #0.719 - drop hd.area.100

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + pop.50, 
            data = lur.dat1) #0.721 - drop pop.50

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + build.50, 
            data = lur.dat1) #0.721 - drop ld.area.100

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + length.50, 
            data = lur.dat1) #0.722 - drop length.50

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + length.200, 
            data = lur.dat1) #0.726 - drop length.200

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + hd.area.200, 
            data = lur.dat1) #0.719 - drop hd.area.200

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + other.area.500, 
            data = lur.dat1) #0.719 - drop other.area.500

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + other.area.1000, 
            data = lur.dat1) #0.719 - drop other.area.1000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + other.area.4000, 
            data = lur.dat1) #0.719 - drop other.area.4000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + pop.100, 
            data = lur.dat1) #0.721 - drop pop.100

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + pop.200, 
            data = lur.dat1) #0.721 - drop pop.200

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + pop.500, 
            data = lur.dat1) #0.721 - drop pop.500

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + pop.1000, 
            data = lur.dat1) #0.721 - drop pop.1000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + pop.2000, 
            data = lur.dat1) #0.721 - drop pop.500

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + build.count.50, 
            data = lur.dat1) #0.719 - drop build.count.50

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + build.count.100, 
            data = lur.dat1) #0.719 - drop build.count.100

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + build.count.200, 
            data = lur.dat1) #0.719 - drop build.count.200

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + build.count.500, 
            data = lur.dat1) #0.719 - drop build.count.500

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + build.count.1000, 
            data = lur.dat1) #0.719 - drop build.count.1000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + build.count.2000, 
            data = lur.dat1) #0.719 - drop build.count.2000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + build.count.4000, 
            data = lur.dat1) #0.719 - drop build.count.4000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + build.count.1000, 
            data = lur.dat1) #0.719 - drop build.count.1000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + bus.count.100, 
            data = lur.dat1) #0.719 - drop bus.count.100

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + bus.count.4000, 
            data = lur.dat1) #0.721 - drop bus.count.4000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + bus.stat.500, 
            data = lur.dat1) #0.720 - drop bus.stat.500

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + bus.stat.2000, 
            data = lur.dat1) #0.721 - drop bus.stat.2000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + bus.stat.4000, 
            data = lur.dat1) #0.722 - drop bus.stat.4000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + hd.area.50, 
            data = lur.dat1) #0.719 - drop hd.area.50

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + length.4000, 
            data = lur.dat1) #0.722 - drop length.4000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + ld.area.200, 
            data = lur.dat1) #0.719 - drop ld.area.200

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + hd.area.500, 
            data = lur.dat1) #0.719 - drop hd.area.500

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + hd.area.1000, 
            data = lur.dat1) #0.719 - drop hd.area.1000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + hd.area.2000, 
            data = lur.dat1) #0.719 - drop hd.area.2000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + hd.area.4000, 
            data = lur.dat1) #0.719 - drop hd.area.4000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + dist_airport, 
            data = lur.dat1) #0.719 - drop dist_airport

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + bus.count.50, 
            data = lur.dat1) #0.721 - drop bus.count.50

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + bus.count.2000, 
            data = lur.dat1) #0.721 - drop bus.count.2000

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + ld.area.50, 
            data = lur.dat1) #0.719 - drop ld.area.50

temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + dist.maj.rd + length.100 + com.area.50, 
            data = lur.dat1) #0.720 - drop com.area.50

###retain variable that makes r^2 at least 0.729 - currently at 0.719

summary(temp)
glance(temp)


temp <-  lm(log(concn) ~ Rel..Hum. + ndvi.50 + log(dist.maj.rd) + length.100 + length.200, 
            data = lur.dat1) ##r2 is 0.726


temp <-  lm(log(concn) ~  log(dist.maj.rd), 
            data = lur.dat1) ##r2 is 0.726


#- add splines if relationship is not linear,

#log transform or use squares or square-root

