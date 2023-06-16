######
#Code for statistical analyses associated with physiological data
######

#Sven Buerki, 30th May 2023

###
#Raw data: 
###

#Table with data on seedlings & treatments
datSeed <- read.csv("HF_GxE_experiment_full_data_Feb162020.csv")

#Fluorescence data T1 and T2
datFlu <- read.csv("Fluorescence_days_1to15.csv")

#Fluorescence T3
datFluT3 <- read.csv("Fluorescence_days_16to21.csv")

###
#Subset datasets to treatment end dates 
###

#Compare data from last day to mitigate temperature effect in greenhouse

##~~~
#For T1 and T2 - Dates
##~~~
datesExpT1T2dates <- unique(datFlu$Date)

datesExpT1T2 <- datesExpT1T2dates[length(datesExpT1T2dates)]

#redT1T2 and target 2x populations
redT1T2 <- subset(datFlu, datFlu$Subsp == "tridentata2x" & datFlu$Date == datesExpT1T2)
redT1T2$Pop <- sapply(strsplit(as.vector(redT1T2$Parent), split = "_"), "[[", 1)
#Switch treatments (T1=well-watered and T2: imposed drought)
redT1T2$Treatment <- sapply(redT1T2$Treatment, switch,  "T1"="T2", "T2"="T1")

##~~~
#For T3 - Dates
##~~~

datesExpT3dates <- unique(datFluT3$Date)
datesExpT3 <- datesExpT3dates[length(datesExpT3dates)]

#redT3 and target 2x populations
redT3 <- subset(datFluT3, datFluT3$Subsp == "tridentata2x" & datFluT3$Date == datesExpT3)

##~~~
#Exclude seedlings allocated to T3 from redT1T2 dataset
##~~~
#Identify target seedlings to remove from redT1T2
SeedToRm <- redT3$SeedlingID[which(redT3$SeedlingID %in% redT1T2$SeedlingID)]

#Exclude seedlings
redT1T2 <- redT1T2[-which(redT1T2$SeedlingID %in% SeedToRm),]
#Add PhiNO 
redT1T2$PhiNO <- 1-(as.numeric(redT1T2$Phi2)+as.numeric(redT1T2$PhiNPQ))


###
# Aggregate data at replication == Family level
###

#We have two individuals per family per treatment
# --> Infer mean of the two individuals for Phi2, PhiNPQ and leaf temp
# To aggregate we need to add a unique ID (= family + Treatment)
redT1T2$UniqueID <- paste0(redT1T2$Parent, "_", redT1T2$Treatment)

#Aggregate values at family level
FamMat <- data.frame(Group = aggregate(redT1T2$Phi2, list(redT1T2$UniqueID), mean)$Group, Phi2 = aggregate(redT1T2$Phi2, list(redT1T2$UniqueID), mean)$x, PhiNO = aggregate(redT1T2$PhiNO, list(redT1T2$UniqueID), mean)$x,
                     PhiNPQ = aggregate(redT1T2$PhiNPQ, list(redT1T2$UniqueID), mean)$x, Leaf_Temperature = aggregate(redT1T2$Leaf_Temperature, list(redT1T2$UniqueID), mean)$x)
#Add pop and treatment            
FamMat$Pop <- sapply(strsplit(FamMat$Group, "_"), "[[", 1)
FamMat$Treatment <- sapply(strsplit(FamMat$Group, "_"), "[[", 3)

#######
#GLMs - Statistical analyses
#######

###~~~
# Test for population effects within treatment (T1 and T2)
###~~~

##T1
T1 <- subset(FamMat, FamMat$Treatment == "T1")
#Phi2
# Not significant
summary(glm(as.vector(T1$Phi2) ~ T1$Pop))
#PhiNO
# Not significant
summary(glm(as.vector(T1$PhiNO) ~ T1$Pop))
#PhiNPQ
# Not significant
summary(glm(as.vector(T1$PhiNPQ) ~ T1$Pop))
#Leaf temp
# Not significant
summary(glm(as.vector(T1$Leaf_Temperature) ~ T1$Pop))

##T2
T2 <- subset(FamMat, FamMat$Treatment == "T2")
#Phi2
# P-value 0.059
summary(glm(as.vector(T2$Phi2) ~ T2$Pop))
#PhiNPQ
# P-value 0.039
summary(glm(as.vector(T2$PhiNPQ) ~ T2$Pop))
#PhiNO
# Not significant
summary(glm(as.vector(T2$PhiNO) ~ T2$Pop))
#Leaf temp
# Not significant
summary(glm(as.vector(T2$Leaf_Temperature) ~ T2$Pop))

###~~~
# Test for phenotypic plasticity at pop level
###~~~

##UTT2
UTT2 <- subset(FamMat, FamMat$Pop == "UTT2")
#Phi2
# Not significant
summary(glm(as.vector(UTT2$Phi2) ~ UTT2$Treatment))
#PhiNPQ
# Not significant
summary(glm(as.vector(UTT2$PhiNPQ) ~ UTT2$Treatment))
#PhiNO
# Not significant
summary(glm(as.vector(UTT2$PhiNO) ~ UTT2$Treatment))
#Leaf temp
# P-value: 0.002
summary(glm(as.vector(UTT2$Leaf_Temperature) ~ UTT2$Treatment))

##IDT3
IDT3 <- subset(FamMat, FamMat$Pop == "IDT3")
#Phi2
# Not significant
summary(glm(as.vector(IDT3$Phi2) ~ IDT3$Treatment))
#PhiNPQ
# Not significant
summary(glm(as.vector(IDT3$PhiNPQ) ~ IDT3$Treatment))
#PhiNO
# Not significant
summary(glm(as.vector(IDT3$PhiNO) ~ IDT3$Treatment))
#Leaf temp
# P-value: 0.022
summary(glm(as.vector(IDT3$Leaf_Temperature) ~ IDT3$Treatment))

#######
#PLOTS 
#######

#3 panels
# We are not showing PhiNO since it is correlated with PhiNPQ
pdf("2023_Figure_3_drought_GxE_revised.pdf")
par(mfrow=c(3,1), mai=c(0.45,0.45,0.1,0.2))
#Leaf temp
boxplot(FamMat$Leaf_Temperature ~ FamMat$Treatment + as.vector(FamMat$Pop), ylim = c(21, 26.5), col = c("red", "red", "blue", "blue"), names = rep(c("T1", "T2"),2), xlab = "", ylab = "", cex.lab=.9)
title(ylab=expression("Leaf temperature (°C)"), mgp=c(1.9,1,0), cex.lab=1)
#IDT3
arrows(x0 = 1, x1 = 2, y0 = 26.3, y1 = 26.3, angle = 25, length = 0.1, code = 3, col = "grey")
text(x = 1.5, y=26.55, "P-value: 0.022", cex = 0.9)
#UTT2
arrows(x0 = 3, x1 = 4, y0 = 26.3, y1 = 26.3, angle = 25, length = 0.1, code = 3, col = "grey")
text(x = 3.5, y=26.55, "P-value: 0.002", cex = 0.9)
title(xlab="Treatment", mgp=c(1.9,1,0), cex.lab=1)
plotrix::corner.label(label='A.',figcorner=T)

#Phi2
boxplot(FamMat$Phi2 ~ FamMat$Treatment + as.vector(FamMat$Pop), ylim = c(0.4, 0.8), col = c("red", "red", "blue", "blue"), names = rep(c("T1", "T2"),2), xlab = "", ylab = "", cex.lab=.9)
title(ylab=expression("Phi2"), mgp=c(1.9,1,0), cex.lab=1)
#T2
arrows(x0 = 2, x1 = 4, y0 = 0.77, y1 = 0.77, angle = 25, length = 0.1, code = 3)
text(x = 3, y=0.79, "P-value: 0.059", cex = 0.9)
title(xlab="Treatment", mgp=c(1.9,1,0), cex.lab=1)
plotrix::corner.label(label='B.',figcorner=T)

#PhiNPQ
boxplot(FamMat$PhiNPQ ~ FamMat$Treatment + as.vector(FamMat$Pop), ylim = c(0.05, 0.45), col = c("red", "red", "blue", "blue"), names = rep(c("T1", "T2"),2), xlab = "", ylab = "", cex.lab=.9)
title(ylab=expression(paste("PhiNPQ", sep="")), mgp=c(1.9,1,0), cex.lab=1)
#T2
arrows(x0 = 2, x1 = 4, y0 = 0.42, y1 = 0.42, angle = 25, length = 0.1, code = 3)
text(x = 3, y=0.44, "P-value: 0.039", cex = 0.9)
title(xlab="Treatment", mgp=c(1.9,1,0), cex.lab=1)
plotrix::corner.label(label='C.',figcorner=T)

# Add a legend
#legend("topright", 
#       legend = c("IDT3", "UTT2", "Population effect", "Phenotypic plasticity"), 
#       col = c("red", "blue", "black", "grey"),
#       lty=c(1,1,1,1),
#       lwd = rep(2, 4),
#       horiz = F , 
       #inset = c(0.1, 0.05),
#       cex = 0.6)

dev.off()
