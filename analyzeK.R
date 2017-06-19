library(ROCR)
library(gtools)
library(infotheo)  # for mutual info
source('functions.R')

## Load the data


data <- NA
for (filename in dir()[grep("*.csv", dir())]) {
	subject <- read.table(filename, header=T, sep=",")
	if ( length(dim(data)) == 0) {
		data <- subject
	} else {
		data <- merge(data, subject, all=T)
	}
}

data <- subset(data, data$Length %in% c(2,3,4))

for (s in unique(data$Subject)) {
	sub <- subset(data, data$participant == s)
	#val <- quantile(sub$encoding.rt, c(0.95))
	val <- mean(sub$encoding.rt) + 3 * sd(sub$encoding.rt)
	data$encoding.rt[data$Subject == s &
	                 data$encoding.rt >= val] <- NA
}
data <- subset(data, !is.na(data$encoding.rt))



sacc <- aggregate(data$execution.corr, list(Subject=data$participant), mean)

acc <- aggregate(data$execution.corr, list(Subject=data$participant, Structure=data$Structure, Length=data$Length), mean)

lacc <- aggregate(acc$x, list(Length=acc$Length, Subject=acc$Subject), min)
zacc <- aggregate(acc$x, list(Subject=acc$Subject), min)
correct <- subset(data, data$execution.corr==1)
correct$LogEncoding <- log(correct$encoding.rt)
correct$EncodingSpeed <- 1/(correct$encoding.rt)

# Additional filtering
 correct <- subset(correct, correct$participant %in% sacc$Subject[sacc$x > 0.50])

#correct <- subset(correct, correct$participant %in% unique(acc$Subject[zacc$x > 0.50]))

correct <- subset(correct, correct$encoding.rt > 0.1)
#correct <- subset(correct, correct$encoding.rt < 10)

correct$Chunks <- 0
correct$Chunks[correct$Structure=="BUxy"] <- 2
correct$Chunks[correct$Structure=="UBxy"] <- 2

correct$Chunks[correct$Structure=="UUBxy"] <- 2
correct$Chunks[correct$Structure=="UBUxy"] <- 3
correct$Chunks[correct$Structure=="BUUxy"] <- 2

correct$Chunks[correct$Structure=="UUUBxy"] <- 2
correct$Chunks[correct$Structure=="UUBUxy"] <- 3
correct$Chunks[correct$Structure=="UBUUxy"] <- 3
correct$Chunks[correct$Structure=="BUUUxy"] <- 2

agg <- aggregate(correct[c('encoding.rt', 'execution.rt', 'LogEncoding', 'EncodingSpeed')], list(Subject=correct$participant, Structure=correct$Structure, Length=correct$Length), mean)
agg$Length <- as.factor(agg$Length)
agg$Subject <- as.factor(agg$Subject)


chunk <- aggregate(correct[c('encoding.rt', 'execution.rt', 'LogEncoding', 'EncodingSpeed')], list(Subject=correct$participant, Chunks=correct$Chunks, Length=correct$Length), mean)
chunk$Length <- as.factor(chunk$Length)
chunk$Subject <- as.factor(chunk$Subject)

plot.by.2factors(agg, "LogEncoding", "Structure", "Length", rng=c(0, 2, 1))
plot.by.2factors(agg, "EncodingSpeed", "Structure", "Length", rng=c(0.1, 0.8, 0.1))

plot.by.2factors(agg, "encoding.rt", "Structure", "Length", rng=c(1, 6, 1))
plot.by.2factors(agg, "execution.rt", "Structure", "Length", rng=c(1, 7, 1))

summary(aov(EncodingSpeed~(Structure) + Error(Subject/Structure), agg))


# ANOTHER DIRECTION

data$LogEncoding <- log(data$encoding.rt)
data$EncodingSpeed <- 1/(data$encoding.rt)

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

l2 <- subset(data, data$Length == 2)
l3 <- subset(data, data$Length == 3)
l4 <- subset(data, data$Length == 4)

# Length THREE, N = 3

l3 <- subset(data, data$Length == 3)

l3acc <- aggregate(l3$execution.corr, list(Subject=l3$participant, Structure=l3$Structure), mean)
names(l3acc)[3] <- "ACC"
l3minacc <- aggregate(l3acc[c('ACC')], list(Subject=l3acc$Subject), min)

l3globalacc <- aggregate(l3$execution.corr, list(Subject=l3$participant), mean)
names(l3globalacc)[2] <- "ACC"


for (s in unique(l3$participant)) {
	sub <- subset(l3, l3$participant == s)
	val <- quantile(sub$encoding.rt, c(0.75)) + 1.5 * IQR(sub$encoding.rt)
	#val <- quantile(sub$encoding.rt, c(0.90))
	#val <- mean(sub$encoding.rt) + 3 * sd(sub$encoding.rt)  
	l3$encoding.rt[l3$participant == s &
	               l3$encoding.rt > val] <- NA
}
l3 <- subset(l3, !is.na(l3$encoding.rt))


l3correct <- subset(l3, l3$execution.corr == 1)
l3correct <- subset(l3correct, !is.na(l3correct$encoding.rt))


# Additional filtering
l3correct <- subset(l3correct, l3correct$encoding.rt > 0.1)
#l3correct <- subset(l3correct, l3correct$encoding.rt < 8)



# Remove bad subjects
#l3correct <- subset(l3correct, l3correct$participant %in% l3globalacc$Subject[l3globalacc$ACC > 8/16])
l3correct <- subset(l3correct, l3correct$participant %in% l3minacc$Subject[l3minacc$ACC > 4/16])

l3agg <- aggregate(l3correct[c('encoding.rt', 'execution.rt', 'LogEncoding', 'EncodingSpeed')], list(Subject=l3correct$participant, Structure=l3correct$Structure, Chunks=l3correct$Chunks), median)
l3agg$Subject <- as.factor(l3agg$Subject)
l3agg$Chunks <- as.factor(l3agg$Chunks)


plot.by.1factor(l3agg, "encoding.rt", "Structure")
plot.by.1factor(l3agg, "LogEncoding", "Structure")

summary(aov(encoding.rt~Structure + Error(Subject/(Structure)), l3agg))
summary(aov(log(encoding.rt)~Structure + Error(Subject/(Structure)), l3agg))

summary(aov(LogEncoding~Structure + Error(Subject/(Structure)), l3agg))
summary(aov(EncodingSpeed~Structure + Error(Subject/(Structure)), l3agg))

# By Chunks

l3chunkagg <- aggregate(l3correct[c('encoding.rt', 'execution.rt', 'LogEncoding', 'EncodingSpeed')], list(Subject=l3correct$participant, Chunks=l3correct$Chunks), median)
l3chunkagg$Subject <- as.factor(l3chunkagg$Subject)
l3chunkagg$Chunks <- as.factor(l3chunkagg$Chunks)


plot.by.1factor(l3chunkagg, "encoding.rt", "Chunks")
plot.by.1factor(l3chunkagg, "LogEncoding", "Chunks")
plot.by.1factor(l3chunkagg, "EncodingSpeed", "Chunks")


summary(aov(encoding.rt ~ Chunks + Error(Subject/Chunks), l3chunkagg))
summary(aov(log(encoding.rt) ~ Chunks + Error(Subject/Chunks), l3chunkagg))

summary(aov(LogEncoding ~ Chunks + Error(Subject/(Chunks)), l3chunkagg))
summary(aov(EncodingSpeed ~ Chunks + Error(Subject/(Chunks)), l3chunkagg))



# Length FOUR N = 4
l4 <- subset(data, data$Length == 4)

l4acc <- aggregate(l4$execution.corr, list(Subject=l4$participant, Structure=l4$Structure), mean)
names(l4acc)[3] <- "ACC"
l4minacc <- aggregate(l4acc[c('ACC')], list(Subject=l4acc$Subject), min)

l4globalacc <- aggregate(l4$execution.corr, list(Subject=l4$participant), mean)
names(l4globalacc)[2] <- "ACC"


for (s in unique(l4$participant)) {
	sub <- subset(l4, l4$participant == s)
	#val <- quantile(sub$encoding.rt, c(0.75)) + 1.5 * IQR(sub$encoding.rt)
	val <- quantile(sub$encoding.rt, c(0.90))
	#val <- mean(sub$encoding.rt) + 3 * sd(sub$encoding.rt)  
	l4$encoding.rt[l4$participant == s &
	               l4$encoding.rt > val] <- NA
}
l4 <- subset(l4, !is.na(l4$encoding.rt))


l4correct <- subset(l4, l4$execution.corr == 1)
l4correct <- subset(l4correct, !is.na(l4correct$encoding.rt))


# Additional filtering
l4correct <- subset(l4correct, l4correct$encoding.rt > 0.1)
#l4correct <- subset(l4correct, l4correct$encoding.rt < 10)



# Remove bad subjects
#l4correct <- subset(l4correct, l4correct$participant %in% l4globalacc$Subject[l4globalacc$ACC > 8/16])
l4correct <- subset(l4correct, l4correct$participant %in% l4minacc$Subject[l4minacc$ACC > 4/16])

l4agg <- aggregate(l4correct[c('encoding.rt', 'execution.rt', 'LogEncoding', 'EncodingSpeed')], list(Subject=l4correct$participant, Structure=l4correct$Structure, Chunks=l4correct$Chunks), median)
l4agg$Subject <- as.factor(l4agg$Subject)
l4agg$Chunks <- as.factor(l4agg$Chunks)


plot.by.1factor(l4agg, "encoding.rt", "Structure")
plot.by.1factor(l4agg, "LogEncoding", "Structure")

summary(aov(encoding.rt~Structure + Error(Subject/(Structure)), l4agg))
summary(aov(log(encoding.rt)~Structure + Error(Subject/(Structure)), l4agg))

summary(aov(LogEncoding~Structure + Error(Subject/(Structure)), l4agg))
summary(aov(EncodingSpeed~Structure + Error(Subject/(Structure)), l4agg))

# By Chunks

l4chunkagg <- aggregate(l4correct[c('encoding.rt', 'execution.rt', 'LogEncoding', 'EncodingSpeed')], list(Subject=l4correct$participant, Chunks=l4correct$Chunks), median)
l4chunkagg$Subject <- as.factor(l4chunkagg$Subject)
l4chunkagg$Chunks <- as.factor(l4chunkagg$Chunks)


plot.by.1factor(l4chunkagg, "encoding.rt", "Chunks")
plot.by.1factor(l4chunkagg, "LogEncoding", "Chunks")
plot.by.1factor(l4chunkagg, "EncodingSpeed", "Chunks")


summary(aov(encoding.rt ~ Chunks + Error(Subject/Chunks), l4chunkagg))
summary(aov(log(encoding.rt) ~ Chunks + Error(Subject/Chunks), l4chunkagg))

summary(aov(LogEncoding ~ Chunks + Error(Subject/(Chunks)), l4chunkagg))
summary(aov(EncodingSpeed ~ Chunks + Error(Subject/(Chunks)), l4chunkagg))

plot.all <- function(length=4, output="All_subjects.pdf") {
	pdf(output)
	layout(matrix(data=1:9, nrow=3, byrow=T))
	for (filename in dir()[grep("*.csv", dir())]) {
		print(filename)
		subject <- read.table(filename, header=T, sep=",")
		
		name <- unique(subject$participant)
		sub <- subset(subject, subject$execution.corr == 1 & !is.na(subject$encoding.rt) & subject$Length == length)
		print(name)
		plot.by.1factor(sub, "encoding.rt", "Structure", title=name)
		mtext(name, side=3, line=0)
		
		plot(encoding.rt ~ Structure, sub, main=name)
		mtext(name, side=3, line=0)
		
		plot(execution.corr ~ Structure, subset(subject, subject$Length == length), main=name)
		mtext(name, side=3, line=0)
	}
	dev.off()
}