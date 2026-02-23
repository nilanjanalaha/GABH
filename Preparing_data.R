
######################## loading data ##############################
# Reading the data
mydata2 <-read.csv("Source/rnawholedata2.csv",
                   header=TRUE, nrows=584, na.strings=c(""," "), skipNul = TRUE)
mydata2 <- mydata2[,1:40]
# Just the data, (excluding strand, arm, id, and start-ccordinate information)
main.mydata2 <- mydata2[,-c(1:4)]

#------------------ Removing NA rows -----------------------------
#Getting rid of the genes which have no data, thus 522 genes
nullvec <- numeric()
for( i in 1:583)
{
  if (sum(is.na(main.mydata2[i,]))==36)
    nullvec <- c(i,nullvec)
}
mydata <- mydata2[-nullvec,]

# ------------------ Centromere data ---------------------------

centro_data <- read.table("Source/centromere_info.txt",
                          header=FALSE)[, c(2,3,4)]
colnames(centro_data) <- name.temp <- c("chrom", "start_coordinate", "end_coordinate")

#Shuffling the rows so X and Y chromosomes are at the end
centro_xy <- centro_data[c(10,11),]
centro_data <- centro_data[-c(10,11),]
centro_data <- rbind(centro_data, centro_xy[1,])
centro_data[23,1] <- "chrX"

#------------ Appending the centromere column--------------


infodata <- mydata[,1]
chr <- paste("chr",1:22,sep='')
chr <- c(chr, "chrX")


#appending empty column
centro <- rep(0, nrow(mydata))
mydata.main <- mydata[, -c(1,2,3,4)]
mydata.info <- mydata[, c(1,2,3,4)]
mydata <- cbind(mydata.info, centro, mydata.main)
nct <- ncol(centro_data)
emptymatp <- emptymatm <- matrix(0, nrow(centro_data), 2)

#Calculating the number of genes above or below centromere in each chromosome
# Also mydata now has a column named centro
# centro is 0 if one gene is before centromere and one otherwise
for (i in 1: 23)
{
  temp <- mydata[mydata$Chr==chr[i], 1:5]
  ind <- which(temp$start_coordinate > centro_data[i,3])
  temp$centro[ind]=1
  mydata[mydata$Chr == chr[i], 1:5] <- temp

  indp <- which(temp$start_coordinate > centro_data[i,3] & temp$strand=="+")
  nbp <- length(indp)
  nup <- nrow(temp[temp$strand == "+",]) -nbp

  indm <- which(temp$start_coordinate > centro_data[i,3] & temp$strand=="-")
  nbm <- length(indm)
  num <- nrow(temp[temp$strand == "-",]) -nbm


  emptymatp[i, ] <- c(nup, nbp)
  emptymatm[i, ] <- c(num, nbm)
}

# Preparing for display
centro_data_1 <- cbind(centro_data, rep("+", 23), emptymatp)
colnames(centro_data_1) <- c(name.temp, "strand", "above_centro", "below_centro")
centro_data_2 <- cbind(centro_data, rep("-", 23), emptymatm)
colnames(centro_data_2) <- c(name.temp, "strand", "above_centro", "below_centro")

centro_info <- rbind(centro_data_1, centro_data_2)

#The following matrix shows how many genes per strand per chromosome
# is above the centromere, or below the centromere

centro_info

#----------------- Auxilliary imputation functions -----------------------

#median imputation

med_imp <- function(x)
{
  x <- as.numeric(x)
  m <- median(x,na.rm = TRUE)
  ifelse(is.na(x),m,x)
}

mean_imp <- function(x)
{
  x <- as.numeric(x)
  m <- mean(x,na.rm = TRUE)
  ifelse(is.na(x),m,x)
}

# median impute
median_impute <- function(data)
{
  for(i in 1:nrow(data))
    data[i,] <- med_imp(as.numeric(data[i,]))
  data
}

# mean impute
mean_impute <- function(data)
{
  for(i in 1:nrow(data))
    data[i,] <- mean_imp(as.numeric(data[i,]))
  data
}

