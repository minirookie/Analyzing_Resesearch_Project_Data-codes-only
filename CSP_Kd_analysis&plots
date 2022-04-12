rm(list = ls())
library(knitr)
knitr::opts_knit$set(progress=FALSE, verbose=FALSE)

################################# AMPPNP Titration Analysis #############################################

setwd("/Users/syu/OneDrive - St. Jude Children's Research Hospital/UDrive/Documents_Fischerlab_Backup/NMR_Data/TROSY/15NMFP013_ANPtitration_600C_Feb212022/ExportedData")
#PC path
#setwd("C:/Users/syu/OneDrive - St. Jude Children's Research Hospital/UDrive/Documents_Fischerlab_Backup/NMR_Data/TROSY/15NMFP013_ANPtitration_600C_Feb212022/ExportedData")
ANP.t0.table = read.fwf("15N-MFP013_0mMAMPPNP.txt",widths =c(18,10,11,14,10,15,15,25),
                   header = F, skip=2,sep = "\t")
ANP.t1.table = read.fwf("15N-MFP013_2mMAMPPNP.txt",widths =c(18,10,11,14,10,15,15,25),
                       header = F, skip=2,sep = "\t")
ANP.t2.table = read.fwf("15N-MFP013_4mMAMPPNP.txt",widths =c(18,10,11,14,10,15,15,25),
                       header = F, skip=2,sep = "\t")
ANP.t3.table = read.fwf("15N-MFP013_6mMAMPPNP.txt",widths =c(18,10,11,14,10,15,15,25),
                        header = F, skip=2,sep = "\t")
ANP.t4.table = read.fwf("15N-MFP013_8mMAMPPNP.txt",widths =c(18,10,11,14,10,15,15,25),
                        header = F, skip=2,sep = "\t")
ANP.t5.table = read.fwf("15N-MFP013_12mMAMPPNP.txt",widths =c(18,10,11,14,10,15,15,25),
                        header = F, skip=2,sep = "\t")
ANP.t6.table = read.fwf("15N-MFP013_16mMAMPPNP.txt",widths =c(18,10,11,14,10,15,15,25),
                        header = F, skip=2,sep = "\t")
ANP.t7.table = read.fwf("15N-MFP013_20mMAMPPNP.txt",widths =c(18,10,11,14,10,15,15,25),
                        header = F, skip=2,sep = "\t")
ANP.t8.table = read.fwf("15N-MFP013_25mMAMPPNP.txt",widths =c(18,10,11,14,10,15,15,25),
                        header = F, skip=2,sep = "\t")
ANP.t9.table = read.fwf("15N-MFP013_30mMAMPPNP.txt",widths =c(18,10,11,14,10,15,15,25),
                        header = F, skip=2,sep = "\t")
ANP.t10.table = read.fwf("15N-MFP013_40mMAMPPNP.txt",widths =c(18,10,11,14,10,15,15,25),
                        header = F, skip=2,sep = "\t")
dim(ANP.t0.table)
head(ANP.t0.table)

ANP.list = list(ANP.t0.table, ANP.t1.table, ANP.t2.table, ANP.t3.table,
                ANP.t4.table, ANP.t5.table, ANP.t6.table, ANP.t7.table,
                ANP.t8.table, ANP.t9.table, ANP.t10.table)

for (k in 1:11) {
  colnames(ANP.list[[k]]) = c("Assignment","w1","w2","Data_Height","S/N","Atom_1",
                              "Atom_2","Note")
 }

ANP.list[[11]][1:2,]
ANP.list[[1]][1:5,6:7]

for (i in 1:11) {
  ANP.list[[i]][,6] = gsub("    N","-N",ANP.list[[i]][,6]) #replace between letter space with a dash
  ANP.list[[i]][,6] = gsub("^ +", "", ANP.list[[i]][,6]) #remove space before the text string
  ANP.list[[i]][,6] = gsub("    $+", "", ANP.list[[i]][,6])  #remove spaces following the text string

  ANP.list[[i]][,7] = gsub("    H","-H",ANP.list[[i]][,7])
  ANP.list[[i]][,7] = gsub("^ +", "", ANP.list[[i]][,7])
  ANP.list[[i]][,7] = gsub("    $+", "", ANP.list[[i]][,7])
  
  ANP.list[[i]][,1] = gsub("^ +", "", ANP.list[[i]][,1])
  ANP.list[[i]][,1] = gsub("     $+", "", ANP.list[[i]][,1])
}

ANP.t0.table[1:5, 6:7]
ANP.list[[1]][1:5,6:7]

#sqrt((ANP.list[[11]][1,3] - ANP.list[[1]][1,3])^2 +
#     (ANP.list[[11]][1,2] - ANP.list[[1]][1,2])^2/5)
#CSP.t1 = vector(length = 125)
#for (i in (1:125)){
#  CSP.t1[i] = sqrt((ANP.list[[2]][i,3] - ANP.list[[1]][i,3])^2 +
#                   (ANP.list[[2]][i,2] - ANP.list[[1]][i,2])^2/5)
#  }

CSP.ANP = data.frame()
for (j in 1:10) {
  for (i in 1:125){
   CSP.ANP[i,j] = sqrt((ANP.list[[j+1]][i,3] - ANP.list[[1]][i,3])^2 +
                      (ANP.list[[j+1]][i,2] - ANP.list[[1]][i,2])^2/5)
  }
}
CSP.ANP[1:5,]

CSP.ANP = cbind(rep(0,125), CSP.ANP)
colnames(CSP.ANP)=c("0mM_ANP","2mM_ANP","4mM_ANP","6mM_ANP","8mM_ANP","12mM_ANP",
                    "16mM_ANP","20mM_ANP", "25mM_ANP","30mM_ANP","40mM_ANP")
row.names(CSP.ANP)= ANP.list[[1]][,1]

CSP.ANP[120:125,]
ANP.titr.ratio=c(0,20,40,60,80,120,160,200,250,300,400)
plot(x=ANP.titr.ratio, y=CSP.ANP[125,], type = "b", col=1, ylim = c(0,0.2),
     xlab = "AMPPNP:CTD ratio", ylab = "Residue CSP", main = "AMPNP CSP plot",
     pch =20, cex.lab=0.9)
lines(ANP.titr.ratio, CSP.ANP[124,], type = "b", col=2, pch =20)
lines(ANP.titr.ratio, CSP.ANP[123,], type = "b", col=3, pch =20)
lines(ANP.titr.ratio, CSP.ANP[122,], type = "b", col=4, pch =20)
axis(1,at=ANP.titr.ratio, labels = ANP.titr.ratio)
legend("bottomright",rownames(CSP.ANP)[122:125], lty = 1, col = 4:1, bty="n")



################################# ADP Titration Analysis #############################################

setwd("/Users/syu/OneDrive - St. Jude Children's Research Hospital/UDrive/Documents_Fischerlab_Backup/NMR_Data/TROSY/15NMFP013_ADPtitration_600B_Nov2021/ExportedData")
#PC path
#setwd("C:/Users/syu/OneDrive - St. Jude Children's Research Hospital/UDrive/Documents_Fischerlab_Backup/NMR_Data/TROSY/15NMFP013_ADPtitration_600B_Nov2021/ExportedData")
ADP.t0.table = read.fwf("15NMFP013_APO.txt",widths =c(18,10,11,14,10,15,15,25),
                        header = F, skip=7, sep = "\t")
ADP.t1.table = read.fwf("15NMFP013_2mMADP.txt",widths =c(18,10,11,14,10,15,15,25),
                        header = F, skip=8, sep = "\t")
ADP.t2.table = read.fwf("15NMFP013_4mMADP.txt",widths =c(18,10,11,14,10,15,15,25),
                        header = F, skip=9, sep = "\t")
ADP.t3.table = read.fwf("15NMFP013_6mMADP.txt",widths =c(18,10,11,14,10,15,15,25),
                        header = F, skip=9, sep = "\t")
ADP.t4.table = read.fwf("15NMFP013_8mMADP.txt",widths =c(18,10,11,14,10,15,15,25),
                        header = F, skip=9,sep = "\t")
ADP.t5.table = read.fwf("15NMFP013_12mMADP.txt",widths =c(18,10,11,14,10,15,15,25),
                        header = F, skip=9, sep = "\t")
ADP.t6.table = read.fwf("15NMFP013_16mMADP.txt",widths =c(18,10,11,14,10,15,15,25),
                        header = F, skip=9, sep = "\t")
ADP.t7.table = read.fwf("15NMFP013_20mMADP.txt",widths =c(18,10,11,14,10,15,15,25),
                        header = F, skip=9, sep = "\t")
ADP.t8.table = read.fwf("15NMFP013_25mMADP.txt",widths =c(18,10,11,14,10,15,15,25),
                        header = F, skip=9, sep = "\t")
ADP.t9.table = read.fwf("15NMFP013_30mMADP.txt",widths =c(18,10,11,14,10,15,15,25),
                        header = F, skip=9, sep = "\t")
ADP.t10.table = read.fwf("15NMFP013_40mMADP.txt",widths =c(18,10,11,14,10,15,15,25),
                         header = F, skip=9, sep = "\t")
dim(ADP.t0.table)
head(ADP.t0.table)

ADP.list = list(ADP.t0.table, ADP.t1.table, ADP.t2.table, ADP.t3.table,
                ADP.t4.table, ADP.t5.table, ADP.t6.table, ADP.t7.table,
                ADP.t8.table, ADP.t9.table, ADP.t10.table)

for (k in 1:11) {
  colnames(ADP.list[[k]]) = c("Assignment","w1","w2","Data_Height","S/N","Atom_1",
                              "Atom_2","Note")
}

ADP.list[[11]][1:2,]
ADP.list[[1]][1:5,6:7]

for (i in 1:11) {
  ADP.list[[i]][,6] = gsub("    N","-N",ADP.list[[i]][,6]) #replace with a dash
  ADP.list[[i]][,6] = gsub("^ +", "", ADP.list[[i]][,6]) #remove space before the text string
  ADP.list[[i]][,6] = gsub("    $+", "", ADP.list[[i]][,6])  #remove spaces following the text string
  
  ADP.list[[i]][,7] = gsub("    H","-H",ADP.list[[i]][,7])
  ADP.list[[i]][,7] = gsub("^ +", "", ADP.list[[i]][,7])
  ADP.list[[i]][,7] = gsub("    $+", "", ADP.list[[i]][,7])
  
  ADP.list[[i]][,1] = gsub("^ +", "", ADP.list[[i]][,1])
  ADP.list[[i]][,1] = gsub("     $+", "", ADP.list[[i]][,1])
}

ADP.t0.table[1:5, 6:7]
ADP.list[[1]][1:5,6:7]

CSP.ADP = data.frame()
for (j in 1:10) {
  for (i in 1:125){
    CSP.ADP[i,j] = sqrt((ADP.list[[j+1]][i,3] - ADP.list[[1]][i,3])^2 +
                          (ADP.list[[j+1]][i,2] - ADP.list[[1]][i,2])^2/5)
  }
}
CSP.ADP[1:5,]

CSP.ADP = cbind(rep(0,125), CSP.ADP)
colnames(CSP.ADP)=c("0mM_ADP","2mM_ADP","4mM_ADP","6mM_ADP","8mM_ADP","12mM_ADP",
                    "16mM_ADP","20mM_ADP", "25mM_ADP","30mM_ADP","40mM_ADP")
row.names(CSP.ADP)= ADP.list[[1]][,1]

ADP.titr.ratio=c(0,20,40,60,80,120,160,200,250,300,400)

plot(x=ADP.titr.ratio, y=CSP.ADP[125,], type = "b", col=1, ylim = c(0,0.45),
     xlab = "ADP:CTD ratio", ylab = "Residue CSP", main = "ADP CSP plot",
     pch =20, cex.lab=0.9)
lines(ADP.titr.ratio, CSP.ADP[124,], type = "b", col=2, pch =20)
lines(ADP.titr.ratio, CSP.ADP[123,], type = "b", col=3, pch =20)
lines(ADP.titr.ratio, CSP.ADP[122,], type = "b", col=4, pch =20)
axis(1,at=ADP.titr.ratio, labels = ADP.titr.ratio)
legend("bottomright",rownames(CSP.ADP)[122:125], lty = 1, col = 4:1, bty="n")

#library(ggplot2)
#library(reshape2)
#CSP.ADP.graph=rbind(CSP.ADP,ADP.titr.ratio)
#ggplot(CSP.ADP.graph,aes(x=CSP.ADP.graph[126,], y=CSP.ADP[125,]))+
#geom_line(aes(x=CSP.ADP.graph[126,], y=CSP.ADP[124,], colour="black")) +
#  scale_y_continuous(expand = c(0,1))

