library(ROCR)
library(gtools)
library(infotheo)  # for mutual info
source('functions.R')
library(ggplot2)
library(retimes)

## Load the data


data <- NA  # Empty holder for data

# Loads the individual data files
for (filename in dir()[grep("raw/*.csv", dir())]) {
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
data$EncodingSpeed <- 1 / (data$encoding.rt)
data$LogExecution <- log(data$execution.rt)
data$ExecutionSpeed <- 1 / (data$execution.rt)


# Chunks --- number of control chunks that separate the instructions

data$Chunks <- 0
data$Chunks[data$Structure=="BUxy"] <- "Chunks = 2"
data$Chunks[data$Structure=="UBxy"] <- "Chunks = 2"

data$Chunks[data$Structure=="UUBxy"] <- "Chunks = 2"
data$Chunks[data$Structure=="UBUxy"] <- "Chunks = 3"
data$Chunks[data$Structure=="BUUxy"] <- "Chunks = 2"

data$Chunks[data$Structure=="UUUBxy"] <- "Chunks = 2"
data$Chunks[data$Structure=="UUBUxy"] <- "Chunks = 3"
data$Chunks[data$Structure=="UBUUxy"] <- "Chunks = 3"
data$Chunks[data$Structure=="BUUUxy"] <- "Chunks = 2"

# Offset --- number of operations to perform before the first binary
# operation is found.

data$Offset <- 0
data$Offset[data$Structure=="BUxy"] <- "2nd"
data$Offset[data$Structure=="UBxy"] <- "1st"

data$Offset[data$Structure=="UUBxy"] <- "1st"
data$Offset[data$Structure=="UBUxy"] <- "2nd"
data$Offset[data$Structure=="BUUxy"] <- "3rd"

data$Offset[data$Structure=="UUUBxy"] <- "1st"
data$Offset[data$Structure=="UUBUxy"] <- "2nd"
data$Offset[data$Structure=="UBUUxy"] <- "3rd"
data$Offset[data$Structure=="BUUUxy"] <- "4th"

data$Length[data$Length==2] <- "Length = 2"
data$Length[data$Length==3] <- "Length = 3"
data$Length[data$Length==4] <- "Length = 4"


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

accLCO <-  aggregate(data[c('execution.corr')],
                     list(Length=data$Length,
                     Participant=data$participant,
                     Chunks=data$Chunks,
                     Offset=data$Offset),
                     mean)

accLC <-  aggregate(data[c('execution.corr')],
                     list(Length=data$Length,
                          Chunks=data$Chunks,
                          Participant=data$participant),
                     mean)

accLC34 <- subset(accLC, accLC$Length != "Length = 2")



summary(aov(asin(sqrt(execution.corr)) ~ (Length * Chunks) + Error(Participant/(Length * Chunks)), accLC34))

summary(aov(asin(sqrt(execution.corr)) ~ (Length * Chunks) + Error(Participant/(Length * Chunks)), accLC))

## --------------------------------------------------------------- ##
##
## ENCODING DATA ANALYSIS
## 
## --------------------------------------------------------------- ##


cdata <- subset(data, data$execution.corr == 1)
cdata <- subset(cdata, cdata$encoding.rt < quantile(cdata$encoding.rt, 0.975)) # 97.5th percentile.
#cdata <- subset(cdata, cdata$encoding.rt > quantile(cdata$encoding.rt, 0.025)) # 99th percentile.

cdata <- subset(cdata, !is.na(cdata$encoding.rt))


cdata$participant <- factor(cdata$participant)
cdata$Chunks <- factor(cdata$Chunks)
cdata$Length <- factor(cdata$Length)
cdata$Offset <- factor(cdata$Offset)

#cdata <- subset(cdata, cdata$participant !="31012")
#cdata <- subset(cdata, cdata$participant !="31024")
#cdata <- subset(cdata, cdata$participant !="31017")


#cdata <- subset(cdata, cdata$Length != "Length = 4" | cdata$participant !="31012")
#cdata <- subset(cdata, cdata$Length != "Length = 4" | cdata$participant !="31024")

cdata$Length <- factor(cdata$Length)

# Effects of Length

dataL <- aggregate(data[c('encoding.rt')], list(Length=data$Length, Participant=data$participant), mean)
summary(aov(encoding.rt ~ Length + Error(Participant/Length), dataL))

# General Effects of Length and Chunks

# Aggregate by Length, Chunk, Offset (LCO)

dataLCO <-  aggregate(cdata[c('encoding.rt', 'LogEncoding', 'EncodingSpeed')],
                      list(Length=cdata$Length,
                           Participant=cdata$participant,
                           Chunks=cdata$Chunks,
                           Offset=cdata$Offset),
                      median)
summary(aov(encoding.rt ~ (Chunks * Length) + Error(Participant/(Chunks * Length)), subset(dataLCO, dataLCO$Length!="Length = 2")))
summary(aov(LogEncoding ~ (Chunks * Length) + Error(Participant/(Chunks * Length)), subset(dataLCO, dataLCO$Length!="Length = 2")))

names(dataLCO)[4] <- "Binary Operator Position"
plot.by.2factors(dataLCO, "encoding.rt", "Binary Operator Position", "Length", rng=c(1,6,1), legpos="bottomright")
plot.by.2factors(dataLCO, "LogEncoding", "Binary Operator Position", "Length", rng=c(0,2,1), legpos="bottomright")

dataLCO$Offset <- dataLCO['Binary Operator Position']
dataLCO$Length <- factor(dataLCO$Length)
dataLCO$Offset <- factor(dataLCO$Offset)
plot.by.2factors(subset(dataLCO, dataLCO$Length!=2), "encoding.rt", "Chunks", "Length", rng=c(0,5.5,1))


## Better overview picture.

figure3 <- function() {
  xs <- barplot(tapply(dataLCO$encoding.rt, 
                     list(dataLCO$Offset,
                            dataLCO$Length),
                       mean),
                beside=T,
                horiz=T,
                xlim=c(0,7),
                xlab="Encoding Time (secs)",
                ylab = "Instruction Type",
                xpd=F,
                main="Encoding Times\nBy Length and Position of Binary Operator",
                legend.text=unique(dataLCO$Offset),
                args.legend=list(x="bottomright", bty="n", title="Position of Binary\nOperation"),
                col=c("grey85", "grey75", "grey65", "grey55")
  )
  mus <- tapply(dataLCO$encoding.rt, list(dataLCO$Offset, dataLCO$Length), mean)
  ses <- tapply(dataLCO$encoding.rt, list(dataLCO$Offset, dataLCO$Length), se)
  arrows(mus, xs, mus + ses, xs, angle=90, length=0.05)
  arrows(mus, xs, mus - ses, xs, angle=90, length=0.05)
  box(bty="o")
  
  text(y=xs[1,1], x=0.5, label="BU", adj=c(0.5, 0.5))
  text(y=xs[2,1], x=0.5, label="UB", adj=c(0.5, 0.5))
  text(y=xs[1,2], x=0.5, label="BUU", adj=c(0.5, 0.5))
  text(y=xs[2,2], x=0.5, label="UBU", adj=c(0.5, 0.5))
  text(y=xs[3,2], x=0.5, label="UUB", adj=c(0.5, 0.5))
  text(y=xs[1,3], x=0.5, label="BUUU", adj=c(0.5, 0.5))
  text(y=xs[2,3], x=0.5, label="UBUU", adj=c(0.5, 0.5))
  text(y=xs[3,3], x=0.5, label="UUBU", adj=c(0.5, 0.5))
  text(y=xs[4,3], x=0.5, label="UUUB", adj=c(0.5, 0.5))
}
# Aggregate by Length and Chunks (LC)

dataLC <- aggregate(dataLCO[c('encoding.rt', 'LogEncoding', 'EncodingSpeed')], list(Length=dataLCO$Length, Participant=dataLCO$Participant, Chunks=dataLCO$Chunks), mean)
summary(aov(encoding.rt ~ (Chunks * Length) + Error(Participant/(Chunks * Length)), subset(dataLC, dataLC$Length != "Length = 2")))
summary(aov(LogEncoding ~ (Chunks * Length) + Error(Participant/(Chunks * Length)), subset(dataLC, dataLC$Length != "Length = 2")))
summary(aov(EncodingSpeed ~ (Chunks * Length) + Error(Participant/(Chunks * Length)), subset(dataLC, dataLC$Length != "Length = 2")))

plot.by.2factors(subset(dataLC, dataLC$Length != "Length = 2"), 
                 "encoding.rt", "Chunks", "Length", 
                 rng=c(0,6,1),
                 xlab="Encoding Times")
t.test(encoding.rt ~ Chunks, paired=T, subset(dataLC, dataLC$Length=="Length = 3"))
t.test(encoding.rt ~ Chunks, paired=T, subset(dataLC, dataLC$Length=="Length = 4"))
dataLC34 <- subset(dataLC, dataLC$Length != "Length = 2")

tapply(dataLC34$encoding.rt, list(dataLC34$Length, dataLC34$Chunks), mean)
tapply(dataLC34$encoding.rt, list(dataLC34$Length, dataLC34$Chunks), sd)



t.test(encoding.rt ~ Chunks, paired=T, aggregate(dataLC34[c('encoding.rt')], list(Chunks=dataLC34$Chunks, Sub=dataLC34$Participant), mean))

dataLC3 <- subset(dataLC, dataLC$Length == "Length = 3")
t.test(encoding.rt ~ Chunks, paired=T, dataLC3)

dataLC4 <- subset(dataLC, dataLC$Length == "Length = 4")
t.test(encoding.rt ~ Chunks, paired=T, dataLC4)

# Directly aggregate over Chunks and Length

dataLC <-  aggregate(cdata[c('encoding.rt', 'LogEncoding', 'EncodingSpeed')],
                      list(Length=cdata$Length,
                           Participant=cdata$participant,
                           Chunks=cdata$Chunks)
                      median)
summary(aov(encoding.rt ~ (Chunks * Length) + Error(Participant/(Chunks * Length)), subset(dataLC, dataLC$Length != "Length = 2")))
summary(aov(LogEncoding ~ (Chunks * Length) + Error(Participant/(Chunks * Length)), subset(dataLC, dataLC$Length != "Length = 2")))
summary(aov(EncodingSpeed ~ (Chunks * Length) + Error(Participant/(Chunks * Length)), subset(dataLC, dataLC$Length != "Length = 2")))

plot.by.2factors(subset(dataLC, dataLC$Length != "Length = 2"), "encoding.rt", "Chunks", "Length", rng=c(0,6,1))
t.test(encoding.rt ~ Chunks, paired=T, subset(dataLC, dataLC$Length=="Length = 3"))
t.test(encoding.rt ~ Chunks, paired=T, subset(dataLC, dataLC$Length=="Length = 4"))
dataLC34 <- subset(dataLC, dataLC$Length != "Length = 2")
t.test(encoding.rt ~ Chunks, paired=T, aggregate(dataLC34[c('encoding.rt')], list(Chunks=dataLC34$Chunks, Sub=dataLC34$Participant), mean))

dataLC3 <- subset(dataLC, dataLC$Length == "Length = 3")
t.test(encoding.rt ~ Chunks, paired=T, dataLC3)

dataLC4 <- subset(dataLC, dataLC$Length == "Length = 4")
t.test(encoding.rt ~ Chunks, paired=T, dataLC4)

dataLC34$Length <- factor(dataLC34$Length)

figure04 <- function() {
  xs <- barplot(tapply(dataLC34$encoding.rt, 
                 list(dataLC34$Chunks,
                      dataLC34$Length),
                 mean),
          beside=T,
          ylim=c(1,6),
          xlab="Instruction Length",
          ylab = "Encoding Times (secs)",
          xpd=F,
          main="Encoding Times\nBy Length and Number of Chunks",
          legend.text=unique(dataLC34$Chunks),
          args.legend=list(x="topleft", bty="n", title="Number of chunks"),
          col=c("grey85","grey65")
          )
  mus <- tapply(dataLC34$encoding.rt, list(dataLC34$Chunks, dataLC34$Length), mean)
  ses <- tapply(dataLC34$encoding.rt, list(dataLC34$Chunks, dataLC34$Length), se)
  arrows(xs, mus, xs, mus + ses, angle=90, length=0.1)
  arrows(xs, mus, xs, mus - ses, angle=90, length=0.1)
  box(bty="o")
}

## --------------------------------------------------------------- ##
##
## Same thing, for EXECUTION
## 
## --------------------------------------------------------------- ##

xdata <- subset(data, data$execution.corr == 1)
xdata <- subset(xdata, !is.na(xdata$execution.rt))
xdata <- subset(xdata, xdata$execution.rt < quantile(cdata$execution.rt, 0.97))

# General Effects of Length and Chunks

xdata$participant <- factor(xdata$participant)
xdata$Chunks <- factor(xdata$Chunks)
xdata$Length <- factor(xdata$Length)
xdata$Offset <- factor(xdata$Offset)

xdataL <- aggregate(xdata[c('execution.rt')], 
                    list(Length=xdata$Length, 
                         Participant=xdata$participant), 
                    mean)
summary(aov(execution.rt ~ Length + Error(Participant/Length), xdataL))



xdataLCO <-  aggregate(xdata[c('execution.rt')], 
                       list(Length=xdata$Length,
                            Participant=xdata$participant,
                            Chunks=xdata$Chunks,
                            Offset=xdata$Offset),
                       median)
summary(aov(execution.rt ~ (Chunks * Length) + Error(Participant/(Chunks * Length)), subset(xdataLCO, xdataLCO$Length!="Length = 2")))
summary(aov(execution.rt ~ (Chunks * Length * Offset) + Error(Participant/(Chunks * Length * Offset)), subset(xdataLCO, xdataLCO$Length!="Length = 2")))



plot.by.2factors(xdataLCO, "execution.rt", "Offset", "Length", 
                 rng=c(1,7,1),
                 legpos="bottomright")


plot.by.2factors(subset(xdataLCO, xdataLCO$Length!="Length = 2"), "execution.rt", "Chunks", "Length", rng=c(2.5,5.5,1))

xdataLCO$Offset <- as.factor(xdataLCO$Offset)
xdataLCO$Length <- as.factor(xdataLCO$Length)

figure5 <- function() {
  xs <- barplot(tapply(xdataLCO$execution.rt, 
                       list(xdataLCO$Offset,
                            xdataLCO$Length),
                       mean),
                beside=T,
                horiz=T,
                xlim=c(0,7),
                xlab="Encoding Time (secs)",
                ylab = "Instruction Type",
                xpd=F,
                main="Execution Times\nBy Length and Position of Binary Operator",
                legend.text=unique(xdataLCO$Offset),
                args.legend=list(x="bottomright", bty="n", title="Position of Binary\nOperation"),
                col=c("grey85", "grey75", "grey65", "grey55")
  )
  mus <- tapply(xdataLCO$execution.rt, list(xdataLCO$Offset, xdataLCO$Length), mean)
  ses <- tapply(xdataLCO$execution.rt, list(xdataLCO$Offset, xdataLCO$Length), se)
  arrows(mus, xs, mus + ses, xs, angle=90, length=0.05)
  arrows(mus, xs, mus - ses, xs, angle=90, length=0.05)
  box(bty="o")
  
  text(y=xs[1,1], x=0.5, label="BU", adj=c(0.5, 0.5))
  text(y=xs[2,1], x=0.5, label="UB", adj=c(0.5, 0.5))
  text(y=xs[1,2], x=0.5, label="BUU", adj=c(0.5, 0.5))
  text(y=xs[2,2], x=0.5, label="UBU", adj=c(0.5, 0.5))
  text(y=xs[3,2], x=0.5, label="UUB", adj=c(0.5, 0.5))
  text(y=xs[1,3], x=0.5, label="BUUU", adj=c(0.5, 0.5))
  text(y=xs[2,3], x=0.5, label="UBUU", adj=c(0.5, 0.5))
  text(y=xs[3,3], x=0.5, label="UUBU", adj=c(0.5, 0.5))
  text(y=xs[4,3], x=0.5, label="UUUB", adj=c(0.5, 0.5))
}


xdataLC <- aggregate(xdataLCO[c('execution.rt')], list(Length=xdataLCO$Length, Participant=xdataLCO$Participant, Chunks=xdataLCO$Chunks), mean)
summary(aov(execution.rt ~ (Chunks * Length) + Error(Participant/(Chunks * Length)), subset(xdataLC, xdataLC$Length!="Length = 2")))
plot.by.2factors(subset(xdataLC, xdataLC$Length!=2), "execution.rt", "Chunks", "Length", rng=c(2.5,6.5,1))


t.test(execution.rt ~ Chunks, paired=T, subset(xdataLC, xdataLC$Length=="Length = 3"))
t.test(execution.rt ~ Chunks, paired=T, subset(xdataLC, xdataLC$Length=="Length = 4"))

xdataLC34 <- subset(xdataLC, xdataLC$Length != "Length = 2")
t.test(execution.rt ~ Chunks, paired=T, aggregate(xdataLC[c('execution.rt')], list(Chunks=xdataLC$Chunks, Sub=xdataLC$Participant), mean))

xdataLC34$Length <- factor(xdataLC34$Length)

# Direct comparison

dataLC$Phase <- "Encoding"
dataLC$Value <- dataLC$encoding.rt
xdataLC$Phase <- "Execution"
xdataLC$Value <- xdataLC$execution.rt

allLC <- merge(dataLC, xdataLC, by=c("Length", "Participant", "Chunks", "Phase", "Value"), all=T)
plot.by.3factors(subset(allLC, allLC$Length!=2), "Value", "Chunks", "Length", "Phase", rng=c(0,6.5,1))
vals <- tapply(allLC$Value, list(allLC$Chunks, allLC$Length, allLC$Phase), mean)
summary(aov(Value ~ (Phase * Chunks * Length) + Error(Participant/(Phase * Chunks * Length)), subset(allLC, allLC$Length!="Length = 2")))

allLC34 <- subset(allLC, allLC$Length!="Length = 2")
allC <- aggregate(allLC34[c("Value")], list(Participant=allLC34$Participant, Chunks=allLC34$Chunks, Phase=allLC34$Phase), mean)
summary(aov(Value ~ (Phase * Chunks) + Error(Participant/(Phase * Chunks)), allC))


figure6 <- function() {
  xs <- barplot(tapply(xdataLC34$execution.rt, 
                       list(xdataLC34$Chunks,
                            xdataLC34$Length),
                       mean),
                beside=T,
                ylim=c(1,7),
                xlab="Instruction Length",
                ylab = "Encoding Times (secs)",
                xpd=F,
                main="Execution Times\nBy Length and Number of Chunks",
                legend.text=unique(xdataLC34$Chunks),
                col=c("grey85", "grey65"),
                args.legend=list(x="topleft", bty="n", title="Number of chunks"),
  )
  mus <- tapply(xdataLC34$execution.rt, list(xdataLC34$Chunks, xdataLC34$Length), mean)
  ses <- tapply(xdataLC34$execution.rt, list(xdataLC34$Chunks, xdataLC34$Length), se)
  arrows(xs, mus, xs, mus + ses, angle=90, length=0.1)
  arrows(xs, mus, xs, mus - ses, angle=90, length=0.1)
  box(bty="o")
  
}

