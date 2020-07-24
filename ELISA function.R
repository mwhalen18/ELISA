ELISA <- function(STD1.conc, dilution.factor, data) {
  t <- as.matrix(data)
  avg.NSB <- mean(t[1,1:2])
  tmNSB <- t-avg.NSB
  avg.B0 <- mean(tmNSB[1,3:4])
  tpercBind <- tmNSB/avg.B0
  tlogit <- logit(tpercBind)
  STD1 <- STD1.conc
  dilution.factor <- dilution.factor
  output <- list()
  STD2 <- STD1 / dilution.factor
  STD3 <- STD2 / dilution.factor
  STD4 <- STD3 / dilution.factor
  STD5 <- STD4 / dilution.factor
  STD6 <- STD5 / dilution.factor
  STD7 <-  STD6 / dilution.factor
  STD.conc <- rbind(STD1, STD2, STD3, STD4, STD5, STD6, STD7)
  rm(STD1, STD2, STD3, STD4, STD5, STD6, STD7)
  #######Calculate Optical Densities For Standards #######
  STD1Density <- mean(tmNSB[2,1:2])
  STD1Densitysd <- sd(tmNSB[2,1:2])
  STD2Density <- mean(tmNSB[3,1:2])
  STD2Densitysd <- sd(tmNSB[3,1:2])
  STD3Density <- mean(tmNSB[4,1:2])
  STD3Densitysd <- sd(tmNSB[4,1:2])
  STD4Density <- mean(tmNSB[5,1:2])
  STD4Densitysd <- sd(tmNSB[5,1:2])
  STD5Density <- mean(tmNSB[6,1:2])
  STD5Densitysd <- sd(tmNSB[6,1:2])
  STD6Density <- mean(tmNSB[7,1:2])
  STD6Densitysd <- sd(tmNSB[7,1:2])
  STD7Density <- mean(tmNSB[8,1:2])
  STD7Densitysd <- sd(tmNSB[8,1:2])
  StdDensities <- rbind(STD1Density, STD2Density, STD3Density, STD4Density, STD5Density, STD6Density, STD7Density)
  STDdensitySD <- rbind(STD1Densitysd , STD2Densitysd , STD3Densitysd , STD4Densitysd , STD5Densitysd , STD6Densitysd , STD7Densitysd )
  Density.Perc.Error <- (STDdensitySD/StdDensities)*100
  #######Calculate Percent Binding For Standards #######
  STD1Binding <- mean(tpercBind[2,1:2])
  STD2Binding <- mean(tpercBind[3,1:2])
  STD3Binding <- mean(tpercBind[4,1:2])
  STD4Binding <- mean(tpercBind[5,1:2])
  STD5Binding <- mean(tpercBind[6,1:2])
  STD6Binding <- mean(tpercBind[7,1:2])
  STD7Binding <- mean(tpercBind[8,1:2])
  StdBinding <- rbind(STD1Binding, STD2Binding, STD3Binding, STD4Binding, STD5Binding, STD6Binding, STD7Binding)
  #######Calculate LOGIT for Standards #######
  STD1LOGIT <- mean(tlogit[2,1:2])
  STD2LOGIT <- mean(tlogit[3,1:2])
  STD3LOGIT <- mean(tlogit[4,1:2])
  STD4LOGIT <- mean(tlogit[5,1:2])
  STD5LOGIT <- mean(tlogit[6,1:2])
  STD6LOGIT <- mean(tlogit[7,1:2])
  STD7LOGIT <- mean(tlogit[8,1:2])
  StdLOGIT <- rbind(STD1LOGIT, STD2LOGIT, STD3LOGIT, STD4LOGIT, STD5LOGIT, STD6LOGIT, STD7LOGIT)
  ####### Create Summary Table #######
  STDsummary <- cbind(STD.conc, StdDensities, StdBinding, StdLOGIT) %>% as.data.frame() %>% rename(Density = V2, Binding = V3, LOGIT = V4, Concentration = V1)
  ####### Generate Standard Curve #######
  standard.curve <- lm(LOGIT ~ log(Concentration), STDsummary)
  y.int <- summary(standard.curve)$coefficients[1,1]
  slope <- summary(standard.curve)$coefficients[2,1]
  std.plot <- ggplot(data = STDsummary, aes(x = Concentration, y = Binding)) + 
    geom_point(shape = 2) + geom_smooth(method = lm, se = FALSE, size = 0.8) + 
    scale_x_log10() +
    ggtitle("Standard Curve") + theme(plot.title = element_text(size = 10, face = "bold")) + theme_bw()
  ####### Calculate Recovery for Standards #######
  trecovery <- exp((tlogit-y.int)/(slope))
  STD1recovery <- mean(trecovery[2,1:2])
  STD1recSD <- sd(trecovery[2,1:2])
  STD2recovery <- mean(trecovery[3,1:2])
  STD2recSD <- sd(trecovery[3,1:2])
  STD3recovery <- mean(trecovery[4,1:2])
  STD3recSD <- sd(trecovery[4,1:2])
  STD4recovery <- mean(trecovery[5,1:2])
  STD4recSD <- sd(trecovery[5,1:2])
  STD5recovery <- mean(trecovery[6,1:2])
  STD5recSD <- sd(trecovery[6,1:2])
  STD6recovery <- mean(trecovery[7,1:2])
  STD6recSD <- sd(trecovery[7,1:2])
  STD7recovery <- mean(trecovery[8,1:2])
  STD7recSD <- sd(trecovery[8,1:2])
  STDrecovery <- rbind(STD1recovery, STD2recovery, STD3recovery, STD4recovery, STD5recovery, STD6recovery, STD7recovery)
  STDrecoverySD <- rbind(STD1recSD, STD2recSD, STD3recSD, STD4recSD, STD5recSD, STD6recSD, STD7recSD)
  
  ####### Generate Standard Summary Report #######
  STDsummary <- cbind(STDsummary, STDrecovery, STDrecoverySD, Density.Perc.Error) %>% as.data.frame() %>% rename(Recovery = STDrecovery, Recovery.StDev = STDrecoverySD)
  STDsummary$Recovery.SEM <- STDsummary$Recovery.StDev / sqrt(2)
  STDsummary$Sample.Perc.CV <- (STDsummary$Recovery.StDev / STDsummary$Recovery)*100
  STDsummary$Perc.Recovery <- (STDsummary$Recovery / STDsummary$Concentration)*100
  STDsummary$Notes <- ifelse(STDsummary$Density.Perc.Error > 5, "BAD STANDARD",
                             ifelse(STDsummary$Sample.Perc.CV > 14.99, "BAD STANDARD", ""))
  #######Calculate Unknowns #######
  #######Calculate Optical Densities For Samples #######
  CtrlHIGHDensity <- mean(tmNSB[2,3:4],na.rm = TRUE)
  CtrlHIGHDensitysd <- sd(tmNSB[2,3:4],na.rm = TRUE)
  CtrlLOWDensity <- mean(tmNSB[3,3:4],na.rm = TRUE)
  CtrlLOWDensitysd <- sd(tmNSB[3,3:4],na.rm = TRUE)
  U1Density <- mean(tmNSB[4,3:4],na.rm = TRUE)
  U1Densitysd <- sd(tmNSB[4,3:4],na.rm = TRUE)
  U2Density <- mean(tmNSB[5,3:4],na.rm = TRUE)
  U2Densitysd <- sd(tmNSB[5,3:4],na.rm = TRUE)
  U3Density <- mean(tmNSB[6,3:4],na.rm = TRUE)
  U3Densitysd <- sd(tmNSB[6,3:4],na.rm = TRUE)
  U4Density <- mean(tmNSB[7,3:4],na.rm = TRUE)
  U4Densitysd <- sd(tmNSB[7,3:4],na.rm = TRUE)
  U5Density <- mean(tmNSB[8,3:4],na.rm = TRUE)
  U5Densitysd <- sd(tmNSB[8,3:4],na.rm = TRUE)
  
  U6Density <- mean(tmNSB[1,5:6],na.rm = TRUE)
  U6Densitysd <- sd(tmNSB[1,5:6],na.rm = TRUE)
  U7Density <- mean(tmNSB[2,5:6],na.rm = TRUE)
  U7Densitysd <- sd(tmNSB[2,5:6],na.rm = TRUE)
  U8Density <- mean(tmNSB[3,5:6],na.rm = TRUE)
  U8Densitysd <- sd(tmNSB[3,5:6],na.rm = TRUE)
  U9Density <- mean(tmNSB[4,5:6],na.rm = TRUE)
  U9Densitysd <- sd(tmNSB[4,5:6],na.rm = TRUE)
  U10Density <- mean(tmNSB[5,5:6],na.rm = TRUE)
  U10Densitysd <- sd(tmNSB[5,5:6],na.rm = TRUE)
  U11Density <- mean(tmNSB[6,5:6],na.rm = TRUE)
  U11Densitysd <- sd(tmNSB[6,5:6],na.rm = TRUE)
  U12Density <- mean(tmNSB[7,5:6],na.rm = TRUE)
  U12Densitysd <- sd(tmNSB[7,5:6],na.rm = TRUE)
  U13Density <- mean(tmNSB[8,5:6],na.rm = TRUE)
  U13Densitysd <- sd(tmNSB[8,5:6],na.rm = TRUE)
  
  
  U14Density <- mean(tmNSB[1,7:8],na.rm = TRUE)
  U14Densitysd <- sd(tmNSB[1,7:8],na.rm = TRUE)
  U15Density <- mean(tmNSB[2,7:8],na.rm = TRUE)
  U15Densitysd <- sd(tmNSB[2,7:8],na.rm = TRUE)
  U16Density <- mean(tmNSB[3,7:8],na.rm = TRUE)
  U16Densitysd <- sd(tmNSB[3,7:8],na.rm = TRUE)
  U17Density <- mean(tmNSB[4,7:8],na.rm = TRUE)
  U17Densitysd <- sd(tmNSB[4,7:8],na.rm = TRUE)
  U18Density <- mean(tmNSB[5,7:8],na.rm = TRUE)
  U18Densitysd <- sd(tmNSB[5,7:8],na.rm = TRUE)
  U19Density <- mean(tmNSB[6,7:8],na.rm = TRUE)
  U19Densitysd <- sd(tmNSB[6,7:8],na.rm = TRUE)
  U20Density <- mean(tmNSB[7,7:8],na.rm = TRUE)
  U20Densitysd <- sd(tmNSB[7,7:8],na.rm = TRUE)
  U21Density <- mean(tmNSB[8,7:8],na.rm = TRUE)
  U21Densitysd <- sd(tmNSB[8,7:8],na.rm = TRUE)
  
  U22Density <- mean(tmNSB[1,9:10],na.rm = TRUE)
  U22Densitysd <- sd(tmNSB[1,9:10],na.rm = TRUE)
  U23Density <- mean(tmNSB[2,9:10],na.rm = TRUE)
  U23Densitysd <- sd(tmNSB[2,9:10],na.rm = TRUE)
  U24Density <- mean(tmNSB[3,9:10],na.rm = TRUE)
  U24Densitysd <- sd(tmNSB[3,9:10],na.rm = TRUE)
  U25Density <- mean(tmNSB[4,9:10],na.rm = TRUE)
  U25Densitysd <- sd(tmNSB[4,9:10],na.rm = TRUE)
  U26Density <- mean(tmNSB[5,9:10],na.rm = TRUE)
  U26Densitysd <- sd(tmNSB[5,9:10],na.rm = TRUE)
  U27Density <- mean(tmNSB[6,9:10],na.rm = TRUE)
  U27Densitysd <- sd(tmNSB[6,9:10],na.rm = TRUE)
  U28Density <- mean(tmNSB[7,9:10],na.rm = TRUE)
  U28Densitysd <- sd(tmNSB[7,9:10],na.rm = TRUE)
  U29Density <- mean(tmNSB[8,9:10],na.rm = TRUE)
  U29Densitysd <- sd(tmNSB[8,9:10],na.rm = TRUE)
  
  U30Density <- mean(tmNSB[1,11:12],na.rm = TRUE)
  U30Densitysd <- sd(tmNSB[1,11:12],na.rm = TRUE)
  U31Density <- mean(tmNSB[2,11:12],na.rm = TRUE)
  U31Densitysd <- sd(tmNSB[2,11:12],na.rm = TRUE)
  U32Density <- mean(tmNSB[3,11:12],na.rm = TRUE)
  U32Densitysd <- sd(tmNSB[3,11:12],na.rm = TRUE)
  U33Density <- mean(tmNSB[4,11:12],na.rm = TRUE)
  U33Densitysd <- sd(tmNSB[4,11:12],na.rm = TRUE)
  U34Density <- mean(tmNSB[5,11:12],na.rm = TRUE)
  U34Densitysd <- sd(tmNSB[5,11:12],na.rm = TRUE)
  U35Density <- mean(tmNSB[6,11:12],na.rm = TRUE)
  U35Densitysd <- sd(tmNSB[6,11:12],na.rm = TRUE)
  U36Density <- mean(tmNSB[7,11:12],na.rm = TRUE)
  U36Densitysd <- sd(tmNSB[7,11:12],na.rm = TRUE)
  U37Density <- mean(tmNSB[8,11:12],na.rm = TRUE)
  U37Densitysd <- sd(tmNSB[8,11:12],na.rm = TRUE)
  
  SampleDensities <- rbind(CtrlHIGHDensity, CtrlLOWDensity,U1Density,U2Density,U3Density,U4Density,U5Density,U6Density,U7Density,U8Density,U9Density,U10Density,
                           U11Density,U12Density,U13Density,U14Density,U15Density,U16Density,U17Density,U18Density,U19Density,U20Density,U21Density,U22Density,
                           U23Density,U24Density,U25Density,U26Density,U27Density,U28Density,U29Density,U30Density,U31Density,U32Density,U33Density,U34Density,
                           U35Density,U36Density,U37Density)
  SampledensitySD <- rbind(CtrlHIGHDensitysd, CtrlLOWDensitysd,U1Densitysd,U2Densitysd,U3Densitysd,U4Densitysd,U5Densitysd,U6Densitysd,U7Densitysd,U8Densitysd,U9Densitysd,
                           U10Densitysd,U11Densitysd,U12Densitysd,U13Densitysd,U14Densitysd,U15Densitysd,U16Densitysd,U17Densitysd,U18Densitysd,U19Densitysd,U20Densitysd,U21Densitysd,
                           U22Densitysd,U23Densitysd,U24Densitysd,U25Densitysd,U26Densitysd,U27Densitysd,U28Densitysd,U29Densitysd,U30Densitysd,U31Densitysd,
                           U32Densitysd,U33Densitysd,U34Densitysd,U35Densitysd,U36Densitysd,U37Densitysd)
  Sample.Density.Perc.Error <- (SampledensitySD/SampleDensities)*100
  
  #######Calculate Percent Binding For Samples #######
  CtrlHIGHBinding <- mean(tpercBind[2,3:4],na.rm = TRUE)
  CtrlLOWBinding <- mean(tpercBind[3,3:4],na.rm = TRUE)
  U1Binding <- mean(tpercBind[4,3:4],na.rm = TRUE)
  U2Binding <- mean(tpercBind[5,3:4],na.rm = TRUE)
  U3Binding <- mean(tpercBind[6,3:4],na.rm = TRUE)
  U4Binding <- mean(tpercBind[7,3:4],na.rm = TRUE)
  U5Binding <- mean(tpercBind[8,3:4],na.rm = TRUE)
  
  U6Binding <- mean(tpercBind[1,5:6],na.rm = TRUE)
  U7Binding <- mean(tpercBind[2,5:6],na.rm = TRUE)
  U8Binding <- mean(tpercBind[3,5:6],na.rm = TRUE)
  U9Binding <- mean(tpercBind[4,5:6],na.rm = TRUE)
  U10Binding <- mean(tpercBind[5,5:6],na.rm = TRUE)
  U11Binding <- mean(tpercBind[6,5:6],na.rm = TRUE)
  U12Binding <- mean(tpercBind[7,5:6],na.rm = TRUE)
  U13Binding <- mean(tpercBind[8,5:6],na.rm = TRUE)
  
  U14Binding <- mean(tpercBind[1,7:8],na.rm = TRUE)
  U15Binding <- mean(tpercBind[2,7:8],na.rm = TRUE)
  U16Binding <- mean(tpercBind[3,7:8],na.rm = TRUE)
  U17Binding <- mean(tpercBind[4,7:8],na.rm = TRUE)
  U18Binding <- mean(tpercBind[5,7:8],na.rm = TRUE)
  U19Binding <- mean(tpercBind[6,7:8],na.rm = TRUE)
  U20Binding <- mean(tpercBind[7,7:8],na.rm = TRUE)
  U21Binding <- mean(tpercBind[8,7:8],na.rm = TRUE)
  
  U22Binding <- mean(tpercBind[1,9:10],na.rm = TRUE)
  U23Binding <- mean(tpercBind[2,9:10],na.rm = TRUE)
  U24Binding <- mean(tpercBind[3,9:10],na.rm = TRUE)
  U25Binding <- mean(tpercBind[4,9:10],na.rm = TRUE)
  U26Binding <- mean(tpercBind[5,9:10],na.rm = TRUE)
  U27Binding <- mean(tpercBind[6,9:10],na.rm = TRUE)
  U28Binding <- mean(tpercBind[7,9:10],na.rm = TRUE)
  U29Binding <- mean(tpercBind[8,9:10],na.rm = TRUE)
  
  U30Binding <- mean(tpercBind[1,11:12],na.rm = TRUE)
  U31Binding <- mean(tpercBind[2,11:12],na.rm = TRUE)
  U32Binding <- mean(tpercBind[3,11:12],na.rm = TRUE)
  U33Binding <- mean(tpercBind[4,11:12],na.rm = TRUE)
  U34Binding <- mean(tpercBind[5,11:12],na.rm = TRUE)
  U35Binding <- mean(tpercBind[6,11:12],na.rm = TRUE)
  U36Binding <- mean(tpercBind[7,11:12],na.rm = TRUE)
  U37Binding <- mean(tpercBind[8,11:12],na.rm = TRUE)
  
  SampleBindings <- rbind(CtrlHIGHBinding, CtrlLOWBinding,U1Binding,U2Binding,U3Binding,U4Binding,U5Binding,U6Binding,U7Binding,U8Binding,U9Binding,
                          U10Binding,U11Binding,U12Binding,U13Binding,U14Binding,U15Binding,U16Binding,U17Binding,U18Binding,U19Binding,U20Binding,U21Binding,
                          U22Binding,U23Binding,U24Binding,U25Binding,U26Binding,U27Binding,U28Binding,U29Binding,U30Binding,U31Binding,U32Binding,
                          U33Binding,U34Binding,U35Binding,U36Binding,U37Binding)
  #######Calculate LOGIT For Samples #######
  CtrlHIGHlogit <- mean(tlogit[2,3:4],na.rm = TRUE)
  CtrlLOWlogit <- mean(tlogit[3,3:4],na.rm = TRUE)
  U1logit <- mean(tlogit[4,3:4],na.rm = TRUE)
  U2logit <- mean(tlogit[5,3:4],na.rm = TRUE)
  U3logit <- mean(tlogit[6,3:4],na.rm = TRUE)
  U4logit <- mean(tlogit[7,3:4],na.rm = TRUE)
  U5logit <- mean(tlogit[8,3:4],na.rm = TRUE)
  
  U6logit <- mean(tlogit[1,5:6],na.rm = TRUE)
  U7logit <- mean(tlogit[2,5:6],na.rm = TRUE)
  U8logit <- mean(tlogit[3,5:6],na.rm = TRUE)
  U9logit <- mean(tlogit[4,5:6],na.rm = TRUE)
  U10logit <- mean(tlogit[5,5:6],na.rm = TRUE)
  U11logit <- mean(tlogit[6,5:6],na.rm = TRUE)
  U12logit <- mean(tlogit[7,5:6],na.rm = TRUE)
  U13logit <- mean(tlogit[8,5:6],na.rm = TRUE)
  
  U14logit <- mean(tlogit[1,7:8],na.rm = TRUE)
  U15logit <- mean(tlogit[2,7:8],na.rm = TRUE)
  U16logit <- mean(tlogit[3,7:8],na.rm = TRUE)
  U17logit <- mean(tlogit[4,7:8],na.rm = TRUE)
  U18logit <- mean(tlogit[5,7:8],na.rm = TRUE)
  U19logit <- mean(tlogit[6,7:8],na.rm = TRUE)
  U20logit <- mean(tlogit[7,7:8],na.rm = TRUE)
  U21logit <- mean(tlogit[8,7:8],na.rm = TRUE)
  
  U22logit <- mean(tlogit[1,9:10],na.rm = TRUE)
  U23logit <- mean(tlogit[2,9:10],na.rm = TRUE)
  U24logit <- mean(tlogit[3,9:10],na.rm = TRUE)
  U25logit <- mean(tlogit[4,9:10],na.rm = TRUE)
  U26logit <- mean(tlogit[5,9:10],na.rm = TRUE)
  U27logit <- mean(tlogit[6,9:10],na.rm = TRUE)
  U28logit <- mean(tlogit[7,9:10],na.rm = TRUE)
  U29logit <- mean(tlogit[8,9:10],na.rm = TRUE)
  
  U30logit <- mean(tlogit[1,11:12],na.rm = TRUE)
  U31logit <- mean(tlogit[2,11:12],na.rm = TRUE)
  U32logit <- mean(tlogit[3,11:12],na.rm = TRUE)
  U33logit <- mean(tlogit[4,11:12],na.rm = TRUE)
  U34logit <- mean(tlogit[5,11:12],na.rm = TRUE)
  U35logit <- mean(tlogit[6,11:12],na.rm = TRUE)
  U36logit <- mean(tlogit[7,11:12],na.rm = TRUE)
  U37logit <- mean(tlogit[8,11:12],na.rm = TRUE)
  
  SampleLogit <- rbind(CtrlHIGHlogit, CtrlLOWlogit,U1logit,U2logit,U3logit,U4logit,U5logit,U6logit,U7logit,U8logit,U9logit,
                       U10logit,U11logit,U12logit,U13logit,U14logit,U15logit,U16logit,U17logit,U18logit,U19logit,U20logit,U21logit,
                       U22logit,U23logit,U24logit,U25logit,U26logit,U27logit,U28logit,U29logit,U30logit,U31logit,U32logit,
                       U33logit,U34logit,U35logit,U36logit,U37logit)
  ####### Create Summary Table #######
  SAMPLESsummary <- cbind(SampleDensities, SampleBindings, SampleLogit) %>% as.data.frame() %>% rename(Density = V1, Binding = V2, LOGIT = V3)
  
  ####### Calculate Recovery for Samples #######
  CtrlHIGHrecovery <- mean(trecovery[2,3:4],na.rm = TRUE)
  CtrlHIGHrecSD <- sd(trecovery[2,3:4],na.rm = TRUE)
  CtrlLOWrecovery <- mean(trecovery[3,3:4],na.rm = TRUE)
  CtrlLOWrecSD <- sd(trecovery[3,3:4],na.rm = TRUE)
  U1recovery <- mean(trecovery[4,3:4],na.rm = TRUE)
  U1recSD <- sd(trecovery[4,3:4],na.rm = TRUE)
  U2recovery <- mean(trecovery[5,3:4],na.rm = TRUE)
  U2recSD <- sd(trecovery[5,3:4],na.rm = TRUE)
  U3recovery <- mean(trecovery[6,3:4],na.rm = TRUE)
  U3recSD <- sd(trecovery[6,3:4],na.rm = TRUE)
  U4recovery <- mean(trecovery[7,3:4],na.rm = TRUE)
  U4recSD <- sd(trecovery[7,3:4],na.rm = TRUE)
  U5recovery <- mean(trecovery[8,3:4],na.rm = TRUE)
  U5recSD <- sd(trecovery[8,3:4],na.rm = TRUE)
  
  U6recovery <- mean(trecovery[1,5:6],na.rm = TRUE)
  U6recSD <- sd(trecovery[1,5:6],na.rm = TRUE)
  U7recovery <- mean(trecovery[2,5:6],na.rm = TRUE)
  U7recSD <- sd(trecovery[2,5:6],na.rm = TRUE)
  U8recovery <- mean(trecovery[3,5:6],na.rm = TRUE)
  U8recSD <- sd(trecovery[3,5:6],na.rm = TRUE)
  U9recovery <- mean(trecovery[4,5:6],na.rm = TRUE)
  U9recSD <- sd(trecovery[4,5:6],na.rm = TRUE)
  U10recovery <- mean(trecovery[5,5:6],na.rm = TRUE)
  U10recSD <- sd(trecovery[5,5:6],na.rm = TRUE)
  U11recovery <- mean(trecovery[6,5:6],na.rm = TRUE)
  U11recSD <- sd(trecovery[6,5:6],na.rm = TRUE)
  U12recovery <- mean(trecovery[7,5:6],na.rm = TRUE)
  U12recSD <- sd(trecovery[7,5:6],na.rm = TRUE)
  U13recovery <- mean(trecovery[8,5:6],na.rm = TRUE)
  U13recSD <- sd(trecovery[8,5:6],na.rm = TRUE)
  
  
  U14recovery <- mean(trecovery[1,7:8],na.rm = TRUE)
  U14recSD <- sd(trecovery[1,7:8],na.rm = TRUE)
  U15recovery <- mean(trecovery[2,7:8],na.rm = TRUE)
  U15recSD <- sd(trecovery[2,7:8],na.rm = TRUE)
  U16recovery <- mean(trecovery[3,7:8],na.rm = TRUE)
  U16recSD <- sd(trecovery[3,7:8],na.rm = TRUE)
  U17recovery <- mean(trecovery[4,7:8],na.rm = TRUE)
  U17recSD <- sd(trecovery[4,7:8],na.rm = TRUE)
  U18recovery <- mean(trecovery[5,7:8],na.rm = TRUE)
  U18recSD <- sd(trecovery[5,7:8],na.rm = TRUE)
  U19recovery <- mean(trecovery[6,7:8],na.rm = TRUE)
  U19recSD <- sd(trecovery[6,7:8],na.rm = TRUE)
  U20recovery <- mean(trecovery[7,7:8],na.rm = TRUE)
  U20recSD <- sd(trecovery[7,7:8],na.rm = TRUE)
  U21recovery <- mean(trecovery[8,7:8],na.rm = TRUE)
  U21recSD <- sd(trecovery[8,7:8],na.rm = TRUE)
  
  U22recovery <- mean(trecovery[1,9:10],na.rm = TRUE)
  U22recSD <- sd(trecovery[1,9:10],na.rm = TRUE)
  U23recovery <- mean(trecovery[2,9:10],na.rm = TRUE)
  U23recSD <- sd(trecovery[2,9:10],na.rm = TRUE)
  U24recovery <- mean(trecovery[3,9:10],na.rm = TRUE)
  U24recSD <- sd(trecovery[3,9:10],na.rm = TRUE)
  U25recovery <- mean(trecovery[4,9:10],na.rm = TRUE)
  U25recSD <- sd(trecovery[4,9:10],na.rm = TRUE)
  U26recovery <- mean(trecovery[5,9:10],na.rm = TRUE)
  U26recSD <- sd(trecovery[5,9:10],na.rm = TRUE)
  U27recovery <- mean(trecovery[6,9:10],na.rm = TRUE)
  U27recSD <- sd(trecovery[6,9:10],na.rm = TRUE)
  U28recovery <- mean(trecovery[7,9:10],na.rm = TRUE)
  U28recSD <- sd(trecovery[7,9:10],na.rm = TRUE)
  U29recovery <- mean(trecovery[8,9:10],na.rm = TRUE)
  U29recSD <- sd(trecovery[8,9:10],na.rm = TRUE)
  
  U30recovery <- mean(trecovery[1,11:12],na.rm = TRUE)
  U30recSD <- sd(trecovery[1,11:12],na.rm = TRUE)
  U31recovery <- mean(trecovery[2,11:12],na.rm = TRUE)
  U31recSD <- sd(trecovery[2,11:12],na.rm = TRUE)
  U32recovery <- mean(trecovery[3,11:12],na.rm = TRUE)
  U32recSD <- sd(trecovery[3,11:12],na.rm = TRUE)
  U33recovery <- mean(trecovery[4,11:12],na.rm = TRUE)
  U33recSD <- sd(trecovery[4,11:12],na.rm = TRUE)
  U34recovery <- mean(trecovery[5,11:12],na.rm = TRUE)
  U34recSD <- sd(trecovery[5,11:12],na.rm = TRUE)
  U35recovery <- mean(trecovery[6,11:12],na.rm = TRUE)
  U35recSD <- sd(trecovery[6,11:12],na.rm = TRUE)
  U36recovery <- mean(trecovery[7,11:12],na.rm = TRUE)
  U36recSD <- sd(trecovery[7,11:12],na.rm = TRUE)
  U37recovery <- mean(trecovery[8,11:12],na.rm = TRUE)
  U37recSD <- sd(trecovery[8,11:12],na.rm = TRUE)
  
  sampleRecovery <- rbind(CtrlHIGHrecovery, CtrlLOWrecovery,U1recovery,U2recovery,U3recovery,U4recovery,U5recovery,U6recovery,U7recovery,U8recovery,U9recovery,U10recovery,
                          U11recovery,U12recovery,U13recovery,U14recovery,U15recovery,U16recovery,U17recovery,U18recovery,U19recovery,U20recovery,U21recovery,U22recovery,
                          U23recovery,U24recovery,U25recovery,U26recovery,U27recovery,U28recovery,U29recovery,U30recovery,U31recovery,U32recovery,U33recovery,U34recovery,
                          U35recovery,U36recovery,U37recovery)
  sampleRecoverysd <- rbind(CtrlHIGHrecSD, CtrlLOWrecSD,U1recSD,U2recSD,U3recSD,U4recSD,U5recSD,U6recSD,U7recSD,U8recSD,U9recSD,
                            U10recSD,U11recSD,U12recSD,U13recSD,U14recSD,U15recSD,U16recSD,U17recSD,U18recSD,U19recSD,U20recSD,U21recSD,
                            U22recSD,U23recSD,U24recSD,U25recSD,U26recSD,U27recSD,U28recSD,U29recSD,U30recSD,U31recSD,
                            U32recSD,U33recSD,U34recSD,U35recSD,U36recSD,U37recSD)
  Sample.Recovery.Perc.Error <- (sampleRecoverysd/SampleDensities)*100
  
  ####### Generate Sample Summary Report #######
  SAMPLESsummary <- cbind(SAMPLESsummary, sampleRecovery, sampleRecoverysd, Sample.Density.Perc.Error) %>% as.data.frame() %>% rename(Recovery = sampleRecovery, Sample.Recovery.StDev = sampleRecoverysd)
  SAMPLESsummary$Recovery.SEM <- SAMPLESsummary$Sample.Recovery.StDev / sqrt(2)
  SAMPLESsummary$Sample.Perc.CV <- (SAMPLESsummary$Sample.Recovery.StDev / SAMPLESsummary$Recovery)*100
  SAMPLESsummary$HIGH.LOW <- ifelse(SAMPLESsummary$Binding < .10, "VERY HIGH",
                                    ifelse(SAMPLESsummary$Binding < .20, "HIGH",
                                           ifelse(SAMPLESsummary$Binding > .90, "VERY LOW",
                                                  ifelse(SAMPLESsummary$Binding > .80, "LOW", ""))))
  SAMPLESsummary$Notes <- ifelse(SAMPLESsummary$Density > 0 & is.na(SAMPLESsummary$Sample.Perc.CV), "SINGLET",
                                 ifelse(SAMPLESsummary$Sample.Perc.CV > 14.99, "RE-RUN", ""))
  return(list(STDsummary, SAMPLESsummary, std.plot))
}

