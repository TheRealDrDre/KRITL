library(ROCR)
library(gtools)
library(infotheo)  # for mutual info
source('functions.R')
library(ggplot2)
library(retimes)

## Load the data


data <- NA  # Empty holder for data

# Loads the individual data files
for (filename in dir()[grep("*.csv", dir())]) {
	subject <- read.table(filename, header=T, sep=",")
	if ( length(dim(data)) == 0) {
		data <- subject
	} else {
		data <- merge(data, subject, all=T)
	}
}

# Remove practice trials
data <- subset(data, data$Length %in% c(2, 3, 4))

# Remove practice trials
data <- subset(data, !is.na(data$encoding.rt))


# ----------------------------------------------------------------- #
# NEW ANALYSIS
# ----------------------------------------------------------------- #

# Create additional fields 
data$LogEncoding <- log(data$encoding.rt)
data$EncodingSpeed <- 1/(data$encoding.rt)

# Chunks --- number of control chunks that separate the instructions

data$Chunks <- 0
data$Chunks[data$Structure=="BUxy"] <- "C2"
data$Chunks[data$Structure=="UBxy"] <- "C2"

data$Chunks[data$Structure=="UUBxy"] <- "C2"
data$Chunks[data$Structure=="UBUxy"] <- "C3"
data$Chunks[data$Structure=="BUUxy"] <- "C2"

data$Chunks[data$Structure=="UUUBxy"] <- "C2"
data$Chunks[data$Structure=="UUBUxy"] <- "C3"
data$Chunks[data$Structure=="UBUUxy"] <- "C3"
data$Chunks[data$Structure=="BUUUxy"] <- "C2"

# Offset --- number of operations to perform before the first binary
# operation is found.

data$Offset <- 0
data$Offset[data$Structure=="BUxy"] <- 1
data$Offset[data$Structure=="UBxy"] <- 0

data$Offset[data$Structure=="UUBxy"] <- 0
data$Offset[data$Structure=="UBUxy"] <- 1
data$Offset[data$Structure=="BUUxy"] <- 2

data$Offset[data$Structure=="UUUBxy"] <- 0
data$Offset[data$Structure=="UUBUxy"] <- 1
data$Offset[data$Structure=="UBUUxy"] <- 2
data$Offset[data$Structure=="BUUUxy"] <- 3



# General analysis and cleanup

data$participant <- factor(data$participant) 
data$Length <- factor(data$Length) 


# Cleanup the problem subjects
#data <- subset(data, !data$participant %in% c(31012, 31024, 31030))
#data <- subset(data, data$participant != 31012 & data$participant != 31030)
data <- subset(data, data$participant != 31030)
length(tapply(data$execution.corr, data$participant, mean))
mean(tapply(data$execution.corr, data$participant, mean))

acc <- aggregate(data[c('execution.corr')], list(Length=data$Length, Participant=data$participant), mean)
summary(aov(execution.corr ~ Length + Error(Participant/Length), acc))

tapply(acc$execution.corr, acc$Length, mean)


## --------------------------------------------------------------- ##
##
## ENCODING DATA ANALYSIS
## 
## --------------------------------------------------------------- ##

# Remove RTs > 10secs
cdata <- subset(data, data$execution.corr == 1)
cdata <- subset(cdata, cdata$encoding.rt < 10)

acc <- aggregate(data[c('encoding.rt')], list(Length=data$Length, Participant=data$participant), mean)
summary(aov(encoding.rt ~ Length + Error(Participant/Length), acc))

# General Effects of Length and Chunks

cdata$participant <- factor(cdata$participant)
cdata$Chunks <- factor(cdata$Chunks)
cdata$Length <- factor(cdata$Length)
cdata$Offset <- factor(cdata$Offset)


acc <-  aggregate(cdata[c('encoding.rt')], list(Length=cdata$Length, Participant=cdata$participant, Chunks=cdata$Chunks, Offset=cdata$Offset), fmu3)
summary(aov(encoding.rt ~ (Chunks * Length) + Error(Participant/(Chunks * Length)), subset(acc, acc$Length!=2)))
plot.by.2factors(acc, "encoding.rt", "Offset", "Length", rng=c(1,6,1))
plot.by.2factors(subset(acc, acc$Length!=2), "encoding.rt", "Chunks", "Length", rng=c(2.5,5.5,1))
t.test(encoding.rt ~ Chunks, paired=T, subset(acc, acc$Length!=2))

acc <- aggregate(acc[c('encoding.rt')], list(Length=acc$Length, Participant=acc$Participant, Chunks=acc$Chunks), mean)
summary(aov(encoding.rt ~ (Chunks * Length) + Error(Participant/(Chunks * Length)), subset(acc, acc$Length!=2)))
plot.by.2factors(subset(acc, acc$Length!=2), "encoding.rt", "Chunks", "Length", rng=c(2.5,5.5,1))
t.test(encoding.rt ~ Chunks, paired=T, subset(acc, acc$Length!=2))

# Divide trials by number of operations

l2 <- subset(cdata, cdata$Length == 2)
l3 <- subset(cdata, cdata$Length == 3)
l4 <- subset(cdata, cdata$Length == 4)

# Length THREE, N = 3

l3chunkagg <- aggregate(l3[c('encoding.rt', 'LogEncoding', 'EncodingSpeed')], list(Subject=l3$participant, Chunks=l3$Chunks, Offset=l3$Offset), mean)
l3chunkagg <- aggregate(l3[c('encoding.rt', 'LogEncoding', 'EncodingSpeed')], list(Subject=l3$participant, Chunks=l3$Chunks), mean)
l3chunkagg$Subject <- as.factor(l3chunkagg$Subject)
l3chunkagg$Chunks <- as.factor(l3chunkagg$Chunks)

l3chunkagg2 <- subset(acc, acc$Length==3)  # Santity check. Should be equal to l3chunkagg

plot.by.1factor(l3chunkagg, "encoding.rt", "Chunks")
plot.by.1factor(l3chunkagg, "LogEncoding", "Chunks")
plot.by.1factor(l3chunkagg, "EncodingSpeed", "Chunks")


summary(aov(encoding.rt ~ Chunks + Error(Subject/Chunks), l3chunkagg))
summary(aov(log(encoding.rt) ~ Chunks + Error(Subject/Chunks), l3chunkagg))

summary(aov(LogEncoding ~ Chunks + Error(Subject/(Chunks)), l3chunkagg))
summary(aov(EncodingSpeed ~ Chunks + Error(Subject/(Chunks)), l3chunkagg))



# Length = 4

l4chunkagg <- aggregate(l4[c('encoding.rt', 'LogEncoding', 'EncodingSpeed')], list(Subject=l4$participant, Chunks=l4$Chunks), mean)
l4chunkagg$Subject <- as.factor(l4chunkagg$Subject)
l4chunkagg$Chunks <- as.factor(l4chunkagg$Chunks)

l4chunkagg2 <- subset(acc, acc$Length==4)  # Santity check. Should be equal to l3chunkagg


plot.by.1factor(l4chunkagg, "encoding.rt", "Chunks")
plot.by.1factor(l4chunkagg, "LogEncoding", "Chunks")
plot.by.1factor(l4chunkagg, "EncodingSpeed", "Chunks")


summary(aov(encoding.rt ~ Chunks + Error(Subject/Chunks), l4chunkagg))
summary(aov(log(encoding.rt) ~ Chunks + Error(Subject/Chunks), l4chunkagg))

summary(aov(LogEncoding ~ Chunks + Error(Subject/(Chunks)), l4chunkagg))
summary(aov(EncodingSpeed ~ Chunks + Error(Subject/(Chunks)), l4chunkagg))


## --------------------------------------------------------------- ##
##
## Same thing, for EXECUTION
## 
## --------------------------------------------------------------- ##

cdata <- subset(data, data$execution.corr == 1)
cdata <- subset(cdata, cdata$execution.rt < 10)

acc <- aggregate(cdata[c('execution.rt')], list(Length=cdata$Length, Participant=cdata$participant), median)
summary(aov(execution.rt ~ Length + Error(Participant/Length), acc))

# General Effects of Length and Chunks
cdata$participant <- factor(cdata$participant)
cdata$Chunks <- factor(cdata$Chunks)
cdata$Length <- factor(cdata$Length)

acc <-  aggregate(cdata[c('execution.rt')], list(Length=cdata$Length, Participant=cdata$participant, Chunks=cdata$Chunks), median)
summary(aov(execution.rt ~ (Chunks * Length) + Error(Participant/(Chunks*Length)), subset(acc, acc$Length!=2)))

plot.by.2factors(subset(cdata, cdata$Length!=2), "execution.rt", "Chunks", "Length", rng=c(0,10,1))

# Divide trials by number of operations

l2 <- subset(cdata, cdata$Length == 2)
l3 <- subset(cdata, cdata$Length == 3)
l4 <- subset(cdata, cdata$Length == 4)

# Length THREE, N = 3

l3chunkagg <- aggregate(l3[c('execution.rt')], list(Subject=l3$participant, Chunks=l3$Chunks), mean)
l3chunkagg$Subject <- as.factor(l3chunkagg$Subject)
l3chunkagg$Chunks <- as.factor(l3chunkagg$Chunks)


plot.by.1factor(l3chunkagg, "execution.rt", "Chunks")


summary(aov(execution.rt ~ Chunks + Error(Subject/Chunks), l3chunkagg))

# Length = 4

l4chunkagg <- aggregate(l4[c('execution.rt', 'execution.rt', 'LogEncoding', 'EncodingSpeed')], list(Subject=l4$participant, Chunks=l4$Chunks), median)
l4chunkagg$Subject <- as.factor(l4chunkagg$Subject)
l4chunkagg$Chunks <- as.factor(l4chunkagg$Chunks)


plot.by.1factor(l4chunkagg, "execution.rt", "Chunks")
summary(aov(execution.rt ~ Chunks + Error(Subject/Chunks), l4chunkagg))

# Fancy stuff

library(retimes)

fmu <- function(x) {
  as.double(mexgauss(x)[1])
}

fmu2 <- function(x, iter=100) {
  k <- timefit(x, iter=iter)
  #k <- timefit(x)
  as.double(k@par[1])
}

fmu3 <- function(x, iter=0) {
  k <- timefit(x, iter=iter)
  #k <- timefit(x)
  as.double(k@par[1])
}


cdata <- subset(data, data$execution.corr == 1)
cdata <- subset(cdata, !is.na(cdata$encoding.rt))
cdata <- subset(cdata, cdata$encoding.rt < 10)

acc <- aggregate(data[c('encoding.rt')], list(Length=data$Length, Participant=data$participant), fmu)
summary(aov(encoding.rt ~ Length + Error(Participant/Length), acc))

# General Effects of Length and Chunks

cdata$participant <- factor(cdata$participant)
cdata$Chunks <- factor(cdata$Chunks)
cdata$Length <- factor(cdata$Length)
cdata$Offset <- factor(cdata$Offset)

acc <-  aggregate(cdata[c('encoding.rt')], list(Length=cdata$Length, Participant=cdata$participant, Offset=cdata$Offset, Chunks=cdata$Chunks), fmu)
summary(aov(encoding.rt ~ (Chunks * Length) + Error(Participant/(Chunks * Length)), subset(acc, acc$Length!=2)))
plot.by.2factors(acc, "encoding.rt", "Offset", "Length", rng=c(0,6,1), legpos="bottomright")

acc <-  aggregate(cdata[c('encoding.rt')], list(Length=cdata$Length, Participant=cdata$participant, Chunks=cdata$Chunks), fmu)
summary(aov(encoding.rt ~ (Chunks * Length) + Error(Participant/(Chunks * Length)), subset(acc, acc$Length!=2)))
plot.by.2factors(acc, "encoding.rt", "Chunks", "Length", rng=c(0,6,1), legpos="bottomright")


tacc <- aggregate(acc[c('encoding.rt')], list(Length=acc$Length, Participant=acc$Participant, Chunks=acc$Chunks), mean)
t.test(encoding.rt ~ Chunks, paired=T, subset(tacc, tacc$Length!=2))
t.test(encoding.rt ~ Chunks, paired=T, subset(tacc, tacc$Length!=2 & tacc$Length==3))
t.test(encoding.rt ~ Chunks, paired=T, subset(tacc, tacc$Length!=2 & tacc$Length==4))


acc <-  aggregate(cdata[c('encoding.rt')], list(Length=cdata$Length, Participant=cdata$participant, Chunks=cdata$Chunks), fmu)
summary(aov(encoding.rt ~ (Chunks * Length) + Error(Participant/(Chunks * Length)), subset(acc, acc$Length!=2)))
plot.by.2factors(subset(acc, acc$Length!=2), "encoding.rt", "Chunks", "Length", rng=c(1,6,1))

# Divide trials by number of operations

l2 <- subset(cdata, cdata$Length == 2)
l3 <- subset(cdata, cdata$Length == 3)
l4 <- subset(cdata, cdata$Length == 4)

# Length THREE, N = 3

l3chunkagg <- aggregate(l3[c('encoding.rt', 'LogEncoding', 'EncodingSpeed')], list(Subject=l3$participant, Chunks=l3$Chunks), fmu)
l3chunkagg$Subject <- as.factor(l3chunkagg$Subject)
l3chunkagg$Chunks <- as.factor(l3chunkagg$Chunks)

plot.by.1factor(l3chunkagg, "encoding.rt", "Chunks")
summary(aov(encoding.rt ~ Chunks + Error(Subject/Chunks), l3chunkagg))
t.test(encoding.rt ~ Chunks, paired=T,l3chunkagg)

# Length FOUR, N = 4

l4chunkagg <- aggregate(l4[c('encoding.rt')], list(Subject=l4$participant, Chunks=l4$Chunks), fmu)
l4chunkagg$Subject <- as.factor(l4chunkagg$Subject)
l4chunkagg$Chunks <- as.factor(l4chunkagg$Chunks)


plot.by.1factor(l4chunkagg, "encoding.rt", "Chunks")

summary(aov(encoding.rt ~ Chunks + Error(Subject/Chunks), l4chunkagg))

t.test(encoding.rt ~ Chunks, paired=T, l4chunkagg)

## --------------------------------------------------------------- ##
##
## EXECUTION, with Ex-Gauss fitting
## 
## --------------------------------------------------------------- ##

cdata <- subset(data, data$execution.corr == 1)
#cdata <- subset(cdata, cdata$execution.rt < 10)

acc <- aggregate(cdata[c('execution.rt')], list(Length=cdata$Length, Participant=cdata$participant), fmu2)
summary(aov(execution.rt ~ Length + Error(Participant/Length), acc))

# General Effects of Length and Chunks
cdata$participant <- factor(cdata$participant)
cdata$Chunks <- factor(cdata$Chunks)
cdata$Length <- factor(cdata$Length)

acc <- aggregate(cdata[c('execution.rt')], list(Length=cdata$Length, Chunks=cdata$Chunks, Participant=cdata$participant), fmu2)
summary(aov(execution.rt ~ (Chunks * Length) + Error(Participant/(Chunks * Length)), subset(acc, acc$Length != 2)))
plot.by.2factors(subset(acc, acc$Length!=2), "execution.rt", "Chunks", "Length", rng=c(1,6,1))


# Divide trials by number of operations

l2 <- subset(cdata, cdata$Length == 2)
l3 <- subset(cdata, cdata$Length == 3)
l4 <- subset(cdata, cdata$Length == 4)

# Length THREE, N = 3

l3chunkagg <- aggregate(l3[c('execution.rt')], list(Subject=l3$participant, Chunks=l3$Chunks), fmu2)
l3chunkagg$Subject <- as.factor(l3chunkagg$Subject)
l3chunkagg$Chunks <- as.factor(l3chunkagg$Chunks)

plot.by.1factor(l3chunkagg, "execution.rt", "Chunks")

summary(aov(execution.rt ~ Chunks + Error(Subject/Chunks), l3chunkagg))

# Length = 4

l4chunkagg <- aggregate(l4[c('execution.rt', 'execution.rt', 'LogEncoding', 'EncodingSpeed')], list(Subject=l4$participant, Chunks=l4$Chunks), fmu2)
l4chunkagg$Subject <- as.factor(l4chunkagg$Subject)
l4chunkagg$Chunks <- as.factor(l4chunkagg$Chunks)

plot.by.1factor(l4chunkagg, "execution.rt", "Chunks")
summary(aov(execution.rt ~ Chunks + Error(Subject/Chunks), l4chunkagg))
