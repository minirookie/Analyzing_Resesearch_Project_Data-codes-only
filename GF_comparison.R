# Load in raw data
setwd("/Users/syu/OneDrive - St. Jude Children's Research Hospital/UDrive/Documents_Fischerlab_Backup/SEC_data")
MFP093.first30ml= read.csv("MFP093_First30mLinGelFiltration_8-12-2021.csv", header = T, sep = "\t",
                            na.strings = "NA", fileEncoding = "UTF-16LE")
head(MFP093.first30ml,3)
  
tail(MFP093.first30ml, 3)

MFP093_097=read.csv("MFP093-097_GelFiltration_8-13-2021.csv", header = T, 
                    sep = "\t",  na.strings = "NA", fileEncoding = "UTF-16LE")

#Construct a new header of the data table
names(MFP093_097)=as.character(unlist(MFP093_097[1,]))
    
colnames(MFP093_097)[9]=c("Injection Log")
MFP093_097=MFP093_097[-1,]
head(MFP093_097, 3)

str(MFP093.first30ml[14728,1])
str(MFP093_097[-1,][,c(1,3,5,7,9)])

# normalize data for analysis and write data into a new table
Comb_093_097 = MFP093_097[-1,]

Comb_093_097$UV= as.numeric(as.character(Comb_093_097$UV))+
                  as.numeric(as.character(MFP093.first30ml[14728,1]))

Comb_093_097$Injection= as.numeric(as.character(Comb_093_097$Injection))+
                 as.numeric(as.character(MFP093.first30ml[14728,1]))

Comb_093_097$`Run Log`= as.numeric(as.character(Comb_093_097$`Run Log`))+
                as.numeric(as.character(MFP093.first30ml[14728,1]))

Comb_093_097$`UV_CUT_TEMP@100,BASEM (1)`= as.numeric(as.character(Comb_093_097$`Run Log`))+
  as.numeric(as.character(MFP093.first30ml[14728,1]))

Comb_093_097$`Injection Log`= as.numeric(as.character(Comb_093_097$`Injection Log`))+
  as.numeric(as.character(MFP093.first30ml[14728,1]))

Comb_093_097$Fraction= as.numeric(as.character(Comb_093_097$Fraction))+
  as.numeric(as.character(MFP093.first30ml[14728,1]))

head(MFP093.first30ml[-1,],3)
typeof(MFP093.first30ml$UV)
head(Comb_093_097[,c(1:8,11,12)])

nrow(MFP093.first30ml[-1,])
colnames(Comb_093_097)[c(1:8,11,12)] = colnames(MFP093.first30ml) 
colnames(Comb_093_097)[9:10]=c("Injection Log","Set Mark")
Full.093to097=data.frame(matrix(NA, nrow=61036, ncol=12))

colnames(Full.093to097) = colnames(Comb_093_097)
Comb_093_097 = data.frame(Comb_093_097)
MFP093.first30ml = data.frame(MFP093.first30ml)

typeof(MFP093.first30ml[,3])
head(Comb_093_097,5)
colnames(Full.093to097)=colnames(Comb_093_097)

Full.093to097[1,c(1:8,11,12)] = as.character(unlist(MFP093.first30ml[1,]))
Full.093to097[1,9:10] = c("ml", "")

Full.093to097[2:14728, 1] = as.numeric(as.character(MFP093.first30ml[-1,1]))
Full.093to097[14729:61036, 1] = as.numeric(as.character(Comb_093_097[,1]))

Full.093to097[2:14728, 2] = as.numeric(as.character(MFP093.first30ml[-1,2]))
Full.093to097[14729:61036, 2] = as.numeric(as.character(Comb_093_097[,2]))
Full.093to097[2:6, 3] = as.numeric(as.character(Comb_093_097[1:5,3]))
Full.093to097[2:6, 5]= as.numeric(as.character(MFP093.first30ml[2:6,5]))
Full.093to097[2:6, 6]= as.character(MFP093.first30ml[2:6,6])
Full.093to097[7:21, 5]= as.numeric(as.character(Comb_093_097[1:15,5]))
Full.093to097[7:21, 6]= as.character(Comb_093_097[1:15,6])
Full.093to097[2:21, 11]= as.numeric(as.character(MFP093.first30ml[2:21,9]))
Full.093to097[2:21, 12]= as.character(MFP093.first30ml[2:21,10])
Full.093to097[22:560, 11]= as.numeric(as.character(Comb_093_097[1:539,11]))
Full.093to097[22:560, 12]= as.character(Comb_093_097[1:539,12])
Full.093to097[2:6, 9]= as.numeric(as.character(Comb_093_097[1:5,9]))
Full.093to097[2:6, 10]= as.character(Comb_093_097[1:5,10])
Full.093to097[2:14699, 7]= as.numeric(as.character(MFP093.first30ml[2:14699,7]))
Full.093to097[2:14699, 8]= as.numeric(as.character(MFP093.first30ml[2:14699,8]))
Full.093to097[14700:14714, 7]= as.numeric(as.character(Comb_093_097[3:17,7]))
Full.093to097[14700:60973, 8]= as.numeric(as.character(Comb_093_097[3:46276,8]))
write.csv(Full.093to097, file = "MFP093-097_GF_corrected.csv", sep = ",")

# Find the peak positions for each sample
str(Full.093to097$X)
VolFlow=as.numeric(Full.093to097[-1,1])
UVreads =as.numeric(Full.093to097[-1,2])
UVProf.93to97=data.frame(VolFlow, UVreads)
head(UVProf.93to97)


peak1.UV=max(UVProf.93to97$UVreads[VolFlow<=100])
peak1.vol= UVProf.93to97$VolFlow[UVreads==peak1.UV ]

peak2.UV=max(UVProf.93to97$UVreads[VolFlow >200 & VolFlow<=300])
peak2.vol= UVProf.93to97$VolFlow[UVreads==peak2.UV ]

peak3.UV=max(UVProf.93to97$UVreads[VolFlow >350 & VolFlow<=450])
peak3.vol= UVProf.93to97$VolFlow[UVreads==peak3.UV ]

peak4.UV=max(UVProf.93to97$UVreads[VolFlow >500 & VolFlow<=600])
peak4.vol= UVProf.93to97$VolFlow[UVreads==peak4.UV ]

peak5.UV=max(UVProf.93to97$UVreads[VolFlow >650 & VolFlow<=750])
peak5.vol= UVProf.93to97$VolFlow[UVreads==peak5.UV ]

Injection.indices = as.numeric(Full.093to097$Injection.Log[2:6])

# Plot the run data with normalized peak positions
plot_93to97 = plot(VolFlow,UVreads, cex.axis=0.75, cex.lab= 0.75, xlab = "",
                   type = "l", pch=19, col="blue", las=1, ylab = "", las=1,
                   tcl=0.2, lab=c(9,11,7), lty =1, ylim = c(0,2000),
                   lwd=1, bty= "n" )

title(xlab = "mL", line = 1.75, cex.axis=0.5, font=30)
title(ylab = "mAU", line = 2.5, cex.axis=0.25)

text(x=c(peak1.vol,peak2.vol,peak3.vol, peak4.vol, peak5.vol), 
     y=c(peak1.UV, peak2.UV, peak3.UV, peak4.UV, peak5.UV)+85, cex=0.75,
     labels=round(c(peak1.vol, peak2.vol, peak3.vol, peak4.vol, peak5.vol),2))


text(x=c(peak1.vol,peak2.vol,peak3.vol, peak4.vol, peak5.vol), col = "blue",
     y=c(peak1.UV, peak2.UV, peak3.UV, peak4.UV, peak5.UV)+175, cex=0.75,
      labels=c("(MFP093)","(MFP094)","(MFP095)","(MFP096)","(MFP097)"))

# Plot chromatography profiles for each sample
MFP093_Profile= UVProf.93to97[VolFlow < Injection.indices[2],]
plot_MFP093 = plot(MFP093_Profile$VolFlow,MFP093_Profile$UVreads,
                   cex.axis=0.75, cex.lab= 0.75, xlab = "",
                   type = "l", pch=19, col="blue", las=1, ylab = "", las=1,
                   tcl=0.2, lab=c(9,11,7), lty =1, ylim = c(0,2000),
                   lwd=1, bty= "n" )

title(xlab = "mL", line = 1.75, cex.axis=0.5, font=30)
title(ylab = "mAU", line = 2.5, cex.axis=0.25)

text(x=c(peak1.vol), y=c(peak1.UV)+85, cex=0.75,
     labels=round(c(peak1.vol),2))
text(x=c(peak1.vol), y=c(peak1.UV)+200, col = "blue", cex=0.75,
     labels=c("MFP093"))


MFP094_Profile= UVProf.93to97[VolFlow > Injection.indices[2] & VolFlow < Injection.indices[3],]
dim(MFP094_Profile)
plot_MFP094 = plot(MFP094_Profile$VolFlow,MFP094_Profile$UVreads,
                   cex.axis=0.75, cex.lab= 0.75, xlab = "",
                   type = "l", pch=19, col="blue", las=1, ylab = "", las=1,
                   tcl=0.2, lab=c(9,11,7), lty =1, ylim = c(0,1800),
                   lwd=1, bty= "n" )

title(xlab = "mL", line = 1.75, cex.axis=0.5, font=30)
title(ylab = "mAU", line = 2.5, cex.axis=0.25)

text(x=c(peak2.vol), y=c(peak2.UV)+85, cex=0.75,
     labels=round(c(peak2.vol),2))
text(x=c(peak2.vol), y=c(peak2.UV)+200, col = "blue", cex=0.75,
     labels=c("MFP094"))


MFP095_Profile= UVProf.93to97[VolFlow > Injection.indices[3] & VolFlow < Injection.indices[4],]
plot_MFP095 = plot(MFP095_Profile$VolFlow,MFP095_Profile$UVreads,
                   cex.axis=0.75, cex.lab= 0.75, xlab = "",
                   type = "l", pch=19, col="blue", las=1, ylab = "", las=1,
                   tcl=0.2, lab=c(9,11,7), lty =1, ylim = c(0,1900),
                   lwd=1, bty= "n" )

title(xlab = "mL", line = 1.75, cex.axis=0.5, font=30)
title(ylab = "mAU", line = 2.5, cex.axis=0.25)

text(x=c(peak3.vol), y=c(peak3.UV)+85, cex=0.75,
     labels=round(c(peak3.vol),2))
text(x=c(peak3.vol), y=c(peak3.UV)+200, col = "blue", cex=0.75,
     labels=c("MFP095"))


MFP096_Profile= UVProf.93to97[VolFlow > Injection.indices[4] & VolFlow < Injection.indices[5],]
plot_MFP096 = plot(MFP096_Profile$VolFlow,MFP096_Profile$UVreads,
                   cex.axis=0.75, cex.lab= 0.75, xlab = "",
                   type = "l", pch=19, col="blue", las=1, ylab = "", las=1,
                   tcl=0.2, lab=c(9,11,7), lty =1, ylim = c(0,1700),
                   lwd=1, bty= "n" )

title(xlab = "mL", line = 1.75, cex.axis=0.5, font=30)
title(ylab = "mAU", line = 2.5, cex.axis=0.25)

text(x=c(peak4.vol), y=c(peak4.UV)+85, cex=0.75,
     labels=round(c(peak4.vol),2))
text(x=c(peak4.vol), y=c(peak4.UV)+200, col = "blue", cex=0.75,
     labels=c("MFP096"))


MFP097_Profile= UVProf.93to97[VolFlow > Injection.indices[5],]
plot_MFP097 = plot(MFP097_Profile$VolFlow,MFP097_Profile$UVreads,
                   cex.axis=0.75, cex.lab= 0.75, xlab = "",
                   type = "l", pch=19, col="blue", las=1, ylab = "", las=1,
                   tcl=0.2, lab=c(9,11,7), lty =1, ylim = c(0,1500),
                   lwd=1, bty= "n" )

title(xlab = "mL", line = 1.75, cex.axis=0.5, font=30)
title(ylab = "mAU", line = 2.5, cex.axis=0.25)

text(x=c(peak5.vol), y=c(peak5.UV)+85, cex=0.75,
     labels=round(c(peak5.vol),2))
text(x=c(peak5.vol), y=c(peak5.UV)+200, col = "blue", cex=0.75,
     labels=c("MFP097"))

