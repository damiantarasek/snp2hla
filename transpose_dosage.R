library(sampling)

# deal with arguments
args <- commandArgs(trailingOnly = TRUE)
dosageFile <- args[1]
famFile <- args[2]
output_file <- args[3]

	
# read fam file
fam <- read.table(file = famFile, header = FALSE)
fam <- subset(fam, select = c('V1','V2'))
colnames(fam) <- c('FID','IID')

# read dosage file
dos <- read.table(file = dosageFile, header = FALSE)
dos <- dos [, -c(2:3)]
trans_dos <- t(dos)
#make first row header
colnames(trans_dos) = trans_dos[1, ]
#remove firs row
trans_dos = trans_dos [-1,]

#merge
new_dosage <- cbind (fam,trans_dos)

#write to file
write.table(new_dosage, file = output_file, quote = FALSE, col.names = TRUE, row.names = FALSE)
