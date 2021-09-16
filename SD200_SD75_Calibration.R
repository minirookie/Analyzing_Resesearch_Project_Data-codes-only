
getwd()
setwd("/Users/syu/OneDrive - St. Jude Children's Research Hospital/UDrive/Documents_Fischerlab_Backup/SEC_data")

SD200_Cali= read.csv("CalibrationRun_Superdex200_8-17-2021.csv", header = T, sep = "\t",
         na.strings = "NA", fileEncoding = "UTF-16LE")
SD75_Cali= read.csv("Calibration_SD75pg_HiLoad16_600_8-18-2021.csv", header = T, sep = "\t",
                     na.strings = "NA", fileEncoding = "UTF-16LE")

# rename the header and remove redundant label
names(SD200_Cali) = as.character(unlist(SD200_Cali[1,]))
names(SD75_Cali) = as.character(unlist(SD75_Cali[1,]))
str(SD200_Cali[4,3])

equil.void.sd200 = data.frame(SD200_Cali$UV, SD200_Cali[,2])
names(equil.void.sd200) =c("Volume", "UV")
equil.void.sd200 = equil.void.sd200[-(1:2),]
head(equil.void.sd200)

equil.void.sd200$Volume = as.numeric(as.character(equil.void.sd200$Volume))
equil.void.sd200$UV = as.numeric(as.character(equil.void.sd200$UV))

# normalize data for analysis
void.sd200= equil.void.sd200[(equil.void.sd200$Volume >=0 & 
                                equil.void.sd200$Volume < as.numeric(as.character(SD200_Cali[4,3]))),]


calicur.sd200 = equil.void.sd200[equil.void.sd200$Volume >= as.numeric(as.character(SD200_Cali[4,3])),]
head(calicur.sd200)
calicur.sd200$Volume = calicur.sd200$Volume - as.numeric(as.character(SD200_Cali[4,3]))

equil.void.sd75 = data.frame(SD75_Cali$UV, SD75_Cali[,2])
names(equil.void.sd75) =c("Volume", "UV")
equil.void.sd75 = equil.void.sd75[-(1:2),]
equil.void.sd75$Volume = as.numeric(as.character(equil.void.sd75$Volume))
equil.void.sd75$UV = as.numeric(as.character(equil.void.sd75$UV))
head(equil.void.sd75)
void.sd75 = equil.void.sd75[(equil.void.sd75$Volume >=0 &
                               equil.void.sd75$Volume < as.numeric(as.character(SD75_Cali[4,3]))),]


calicur.sd75 = equil.void.sd75[equil.void.sd75$Volume >= as.numeric(as.character(SD75_Cali[4,3])),]
head(calicur.sd75)
calicur.sd75$Volume = calicur.sd75$Volume - as.numeric(as.character(SD75_Cali[4,3]))

# Plot FPLC profiles of molecular standards - void 
plot_void.sd200 = plot(void.sd200$Volume, void.sd200$UV, cex.axis=0.75, cex.lab= 0.75, xlab = "",
                   type = "l", col="blue", las=1, ylab = "", las=1, ylim = c(0,280),
                   tcl=0.2, lab=c(9,11,7), lty =1, 
                   lwd=1, bty= "n" )
title(xlab = "mL", line = 1.9, cex.axis=0.5, font=30)
title(ylab = "mAU", line = 2.5, cex.axis=0.25)

voidvolUV.sd200= round(max(void.sd200$UV),6)
voidvol.sd200 = void.sd200$Volume[void.sd200$UV == voidvolUV.sd200]

text(x=c(voidvol.sd200), y=c(voidvolUV.sd200)+20, cex=0.75,
     labels=round(c(voidvol.sd200),2))
text(x=c(voidvol.sd200), y=c(voidvolUV.sd200)+40, cex=0.75, col = "blue",
     labels= "HiLoad 16/600 SD200 Void")


# Plot FPLC profiles of molecular standards - ladder series with normalized peak position
plot_calicur.sd200 = plot(calicur.sd200$Volume, calicur.sd200$UV, cex.axis=0.75, cex.lab= 0.75, xlab = "",
                          type = "l", col="blue", las=1, ylab = "", las=1, ylim = c(0,110),
                          tcl=0.2, lab=c(18,11,7), lty =1, 
                          lwd=1, bty= "n" )
title(xlab = "mL", line = 1.9, cex.axis=0.5, font=30)
title(ylab = "mAU", line = 2.5, cex.axis=0.25)

# Locate all the peak positions
s1.sd200 = max(calicur.sd200$UV[calicur.sd200$Volume <70])
s1po.sd200 = calicur.sd200$Volume[calicur.sd200$UV == s1.sd200]
s2.sd200 = max(calicur.sd200$UV[calicur.sd200$Volume >70 & calicur.sd200$Volume < 80])
s2po.sd200 = calicur.sd200$Volume[calicur.sd200$UV == s2.sd200]
s3.sd200 = max(calicur.sd200$UV[calicur.sd200$Volume >85 & calicur.sd200$Volume < 95])
s3po.sd200 = calicur.sd200$Volume[calicur.sd200$UV == s3.sd200]
s4.sd200 = max(calicur.sd200$UV[calicur.sd200$Volume >95 & calicur.sd200$Volume < 110])
s4po.sd200 = calicur.sd200$Volume[calicur.sd200$UV == s4.sd200] 
s5.sd200 = max(calicur.sd200$UV[calicur.sd200$Volume >110 & calicur.sd200$Volume < 125])
s5po.sd200 = calicur.sd200$Volume[calicur.sd200$UV == s5.sd200] 
s6.sd200 = max(calicur.sd200$UV[calicur.sd200$Volume >130 & calicur.sd200$Volume < 140])
s6po.sd200 = calicur.sd200$Volume[calicur.sd200$UV == s6.sd200] 

text(x=c(s1po.sd200,s2po.sd200,s3po.sd200,s4po.sd200,s5po.sd200,s6po.sd200), 
     y=c(s1.sd200, s2.sd200, s3.sd200, s4.sd200, s5.sd200, s6.sd200)+15, cex=0.6,
     labels=round(c(s1po.sd200,s2po.sd200,s3po.sd200,s4po.sd200,s5po.sd200,s6po.sd200),2))

text(x=c(s1po.sd200,s2po.sd200,s3po.sd200,s4po.sd200,s5po.sd200,s6po.sd200), col = " dark green",
     y=c(s1.sd200, s2.sd200, s3.sd200, s4.sd200, s5.sd200, s6.sd200)+8, cex=0.6,
     labels=c("200kDa","150kDa","66kDa","29kDa","12.4kDa", "TCEP?"))



# Plot FPLC profiles of molecular standards - void for the second equipment
plot_void.sd75 = plot(void.sd75$Volume, void.sd75$UV, cex.axis=0.75, cex.lab= 0.75, xlab = "",
                       type = "l", col="blue", las=1, ylab = "", las=1, ylim = c(0,300),
                       tcl=0.2, lab=c(9,11,7), lty =1, 
                       lwd=1, bty= "n" )

title(xlab = "mL", line = 1.9, cex.axis=0.5, font=30)
title(ylab = "mAU", line = 2.5, cex.axis=0.25)

voidvolUV.sd75= round(max(void.sd75$UV),6)
voidvol.sd75 = void.sd75$Volume[void.sd75$UV == voidvolUV.sd75]

text(x=c(voidvol.sd75), y=c(voidvolUV.sd75)+20, cex=0.75,
     labels=round(c(voidvol.sd75),2))
text(x=c(voidvol.sd75), y=c(voidvolUV.sd75)+40, cex=0.75, col = "blue",
     labels= "HiLoad 16/600 SD75 Void")


# Plot FPLC profiles of molecular standards with normalized peak position for 2nd equipment
plot_calicur.sd75 = plot(calicur.sd75$Volume, calicur.sd75$UV, cex.axis=0.75, cex.lab= 0.75, xlab = "",
                          type = "l", col="blue", las=1, ylab = "", las=1, ylim = c(0,150),
                          tcl=0.2, lab=c(18,11,7), lty =1, 
                          lwd=1, bty= "n" )
title(xlab = "mL", line = 1.9, cex.axis=0.5, font=30)
title(ylab = "mAU", line = 2.5, cex.axis=0.25)

# Locate peak positions for the standards
s1.sd75 = max(calicur.sd75$UV[calicur.sd75$Volume <55])
s1po.sd75 = calicur.sd75$Volume[calicur.sd75$UV == s1.sd75]

s2.sd75 = max(calicur.sd75$UV[calicur.sd75$Volume >50 &
                                calicur.sd75$Volume <60])
s2po.sd75 = calicur.sd75$Volume[calicur.sd75$UV == s2.sd75]
s3.sd75 = max(calicur.sd75$UV[calicur.sd75$Volume >65 &
                                calicur.sd75$Volume <75])
s3po.sd75 = calicur.sd75$Volume[calicur.sd75$UV == s3.sd75]
s4.sd75 = max(calicur.sd75$UV[calicur.sd75$Volume >80 &
                                calicur.sd75$Volume <90])
s4po.sd75 = calicur.sd75$Volume[calicur.sd75$UV == s4.sd75]
s5.sd75 = max(calicur.sd75$UV[calicur.sd75$Volume >110 &
                                calicur.sd75$Volume <125])
s5po.sd75 = calicur.sd75$Volume[calicur.sd75$UV == s5.sd75]
s6.sd75 = max(calicur.sd75$UV[calicur.sd75$Volume >140 &
                                calicur.sd75$Volume <150])
s6po.sd75 = calicur.sd75$Volume[calicur.sd75$UV == s6.sd75]

text(x=c(s1po.sd75,s2po.sd75,s3po.sd75,s4po.sd75,s5po.sd75,s6po.sd75), 
     y=c(s1.sd75, s2.sd75, s3.sd75, s4.sd75, s5.sd75, s6.sd75)+15, cex=0.6,
     labels=round(c(s1po.sd75,s2po.sd75,s3po.sd75,s4po.sd75,s5po.sd75,s6po.sd75),2))


text(x=c(s1po.sd75,s2po.sd75,s3po.sd75,s4po.sd75,s5po.sd75,s6po.sd75), col = " dark green",
     y=c(s1.sd75, s2.sd75, s3.sd75, s4.sd75, s5.sd75, s6.sd75)+8, cex=0.6,
     labels=c("200kDa","150kDa","66kDa","29kDa","12.4kDa", "TCEP?"))

# generate a brief summary for the calibration
cali.200vs75 = data.frame(c("2000kDa","200kDa","150kDa","66kDa","29kDa","12.4kDa"),
           round(c(voidvol.sd200,s1po.sd200,s2po.sd200,s3po.sd200,s4po.sd200,s5po.sd200),2),
           round(c(voidvol.sd75,s1po.sd75,s2po.sd75,s3po.sd75,s4po.sd75,s5po.sd75),2))
names(cali.200vs75)= c("MW","SD200.Position (mL)","SD75.Position (mL)")
cali.200vs75
