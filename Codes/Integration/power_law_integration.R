#***********************************
rm(list=ls())
cat("\014")
# Power law integration from Giering et al. -  PAP cruise
#**********************************

#working directory
root.folder = "your_path_to_data"
setwd(paste0(root.folder,"Data")

#Load data
df<-read.table("pour_power_law", header=T, sep="\t", na.strings = "NA", dec = ".")
head("df")

bp <- df$mgC_m_3d_1

Z1 =  135 # Upper bound from RUBALIZ
Z2 = 726 # Lower bound from RUBALIZ

#=========================
# POWER-LAW ANT
#=========================
bp <- df$mgC_m_3d_1
uptake1<-log(bp)
z1 <- log(df$Depth_m)

lm.power <- lm(uptake1~z1)
sumary = summary(lm.power)

#=======================
# Primitive integration
#=======================

a = lm.power$coefficients[2] # Take a de ax+b
b = lm.power$coefficients[1]  # Then b
sda = sumary$coefficients[4]
sdb = sumary$coefficients[3]
sda
sdb
# As lnC = lnalpha - beta*lnz, one has:

beta = -a
alpha = exp(b)

Cintegre1 <- (alpha/(1-beta))*(Z2^(1-beta)-Z1^(1-beta))
Cintegre1
