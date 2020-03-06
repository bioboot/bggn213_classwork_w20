#' ---
#' title: "Class 5: Data visualization and graphics in R"
#' author: "Barry Grant"
#' date: "Jan, 2020"
#' ---

# Class 5
# Data visualization and graphics in R

plot(1:5, col="blue", typ="o")

# Section 2A.
# Read the data file "weight_chart.txt"

baby <- read.table("bimm143_05_rstats/weight_chart.txt",
                   header = TRUE)

plot(baby$Age, baby$Weight,
     type = "o", col="blue",
     pch=15, main="Some old Title",
     xlab="My cool axis label")

plot(1:5, pch=1:5, cex=1:5, col=1:5)


## Read a new file of mouse genome features

mouse <- read.table("bimm143_05_rstats/feature_counts.txt",
                    header=TRUE, sep="\t")

dotchart(mouse$Count, labels = mouse$Feature)

barplot(mouse$Count, names.arg = mouse$Feature, las=2)

# Flip this plot
old.par <- par()$mar
par(mar=c(5, 11, 2, 2))
barplot(mouse$Count, names.arg = mouse$Feature, las=2,
        horiz = TRUE)

## 
mf <- read.delim("bimm143_05_rstats/male_female_counts.txt")

barplot(mf$Count, names.arg = mf$Sample, 
        col=rainbow( nrow(mf) ))



