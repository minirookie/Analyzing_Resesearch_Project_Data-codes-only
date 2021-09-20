rm(list = ls())

setwd("/Users/syu/OneDrive - St. Jude Children's Research Hospital/UDrive/Documents_Fischerlab_Backup/Crystallization_screens")
perform.table = read.table("LCPMosquito_PerformanceStatisticsMiteGenPlatesDec2017.csv", 
         header = T, sep = ",", as.is = T)

count.table = aggregate(data.frame(perform.table$Left.drop, perform.table$Right.drop),
          by = list(perform.table$Row, perform.table$Column), FUN = sum)
colnames(count.table) = c("Row", "Column", "Left.drop.tot", "Right.drop.tot")
count.table

left.missing.count = matrix(5 - count.table$Left.drop.tot, nrow = 8, ncol = 12)
row.names(left.missing.count) = c("A","B","C","D","E","F","G","H")
colnames(left.missing.count) = 1:12

right.missing.count = matrix(5 - count.table$Right.drop.tot, nrow = 8, ncol = 12)
row.names(right.missing.count) = c("A","B","C","D","E","F","G","H")
colnames(right.missing.count) = 1:12

left.missing.per = paste(100*left.missing.count/5, "%", sep = "")
right.missing.per = paste(100*right.missing.count/5, "%", sep = "")

left.missing.prop = matrix(left.missing.per, nrow = 8, ncol = 12)
row.names(left.missing.prop) = c("A","B","C","D","E","F","G","H")
colnames(left.missing.prop) = 1:12

right.missing.prop = matrix(right.missing.per, nrow = 8, ncol = 12)
row.names(right.missing.prop) = c("A","B","C","D","E","F","G","H")
colnames(right.missing.prop) = 1:12

left.missing.count
right.missing.count

left.missing.prop
right.missing.prop
