library(ROCR)
library(gtools)
library(infotheo)  # for mutual info
source('functions.R')

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

# Divide trials by number of operations

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


l3 <- subset(l3, !is.na(l3$encoding.rt))


l3correct <- subset(l3, l3$execution.corr == 1)
l3correct <- subset(l3correct, !is.na(l3correct$encoding.rt))


# Additional filtering
l3correct <- subset(l3correct, l3correct$encoding.rt > 0.1)
l3correct <- subset(l3correct, l3correct$encoding.rt < 10)



# Remove bad subjects
#l3correct <- subset(l3correct, l3correct$participant %in% l3globalacc$Subject[l3globalacc$ACC > 8/16])
l3correct <- subset(l3correct, l3correct$participant %in% l3minacc$Subject[l3minacc$ACC > 4/16])

l3agg <- aggregate(l3correct[c('encoding.rt', 'execution.rt', 'LogEncoding', 'EncodingSpeed')], list(Subject=l3correct$participant, Structure=l3correct$Structure, Chunks=l3correct$Chunks), median)
l3agg$Subject <- as.factor(l3agg$Subject)
l3agg$Chunks <- as.factor(l3agg$Chunks)


plot.by.1factor(l3agg, "encoding.rt", "Structure")
plot.by.1factor(l3agg, "LogEncoding", "Structure")

plot.by.1factor(l3agg, "execution.rt", "Structure")

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


# ----------------------------------------------------------------- #
# LENGTH FOUR --- N = 4
# ----------------------------------------------------------------- #

l4 <- subset(data, data$Length == 4)

#l4 <- subset(l4, l4$execution.corr == 1)


l4acc <- aggregate(l4$execution.corr, list(Subject=l4$participant, Structure=l4$Structure), mean)
names(l4acc)[3] <- "ACC"
l4minacc <- aggregate(l4acc[c('ACC')], list(Subject=l4acc$Subject), min)

l4globalacc <- aggregate(l4$execution.corr, list(Subject=l4$participant), mean)
names(l4globalacc)[2] <- "ACC"

l4 <- subset(l4, !is.na(l4$encoding.rt))

l4correct <- subset(l4, l4$execution.corr == 1)
l4correct <- subset(l4correct, !is.na(l4correct$encoding.rt))


# Additional filtering
l4correct <- subset(l4correct, l4correct$encoding.rt > 0.1)
l4correct <- subset(l4correct, l4correct$encoding.rt < 10)

# Remove bad subjects
l4correct <- subset(l4correct, l4correct$participant %in% l4globalacc$Subject[l4globalacc$ACC > 4/16])

# Aggregate
l4agg <- aggregate(l4correct[c('encoding.rt', 'execution.rt', 'LogEncoding', 'EncodingSpeed')], list(Subject=l4correct$participant, Structure=l4correct$Structure), mean)
l4agg$Subject <- as.factor(l4agg$Subject)
#l4agg$Chunks <- as.factor(l4agg$Chunks)


plot.by.1factor(l4agg, "encoding.rt", "Structure")
plot.by.1factor(l4agg, "LogEncoding", "Structure")

plot.by.1factor(l4agg, "execution.rt", "Structure")

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

plot.by.1factor(l4chunkagg, "execution.rt", "Chunks")


summary(aov(encoding.rt ~ Chunks + Error(Subject/Chunks), l4chunkagg))
summary(aov(log(encoding.rt) ~ Chunks + Error(Subject/Chunks), l4chunkagg))

summary(aov(LogEncoding ~ Chunks + Error(Subject/(Chunks)), l4chunkagg))
summary(aov(EncodingSpeed ~ Chunks + Error(Subject/(Chunks)), l4chunkagg))


summary(aov(execution.rt ~ Chunks + Error(Subject/Chunks), l4chunkagg))
summary(aov(log(execution.rt) ~ Chunks + Error(Subject/Chunks), l4chunkagg))


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