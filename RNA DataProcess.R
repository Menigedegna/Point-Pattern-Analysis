library(matrixStats) # rowMins()
require( spam) #to calculate distances
library(data.table)
library(ggplot2)
library(reshape2)
#Create function to calcualte distances for a postiion table
CalculateSpotDistance<-function(PositionTableX, PositionTableY, DistanceThreshold, fileName, plotType){
  #Calculate nearest distance
  #Here Posistion Y: DATA1, because we collect @entries with @colindices
  ds <- nearest.dist(x=PositionTableX, y=PositionTableY, method = "euclidean", delta=DistanceThreshold, upper=NULL)
  #Plot histogram = optional
  ##for linux
  #tiff(paste("/home/nani/Documents/Plots/", fileName,"_",plotType ,".tiff", sep=""))
  ##for windows Z:\\Result0309\\s20_Wt"
  #tiff(paste("C:\\Documents\Plots\\", fileName,"_",plotType ,".tiff", sep=""))
  ########
  #hist(diff( ds@entries),main=paste("Distances frequency within delta - ", plotType, sep=""), sub=fileName)
  #dev.off() 
  #Get data frame containing min, max, mean distance to next point and number of spots close by (below delta threshold)
  ds@colindices<-ds@colindices[ds@entries>0]
  ds@entries<-ds@entries[ds@entries>0]
  dt <- as.data.table(ds@colindices)[, list(list(ds@entries[.I])), by = ds@colindices]
  CollectDist <- data.frame(SpotId = unlist(dt$`ds@colindices`),
                            NumbSpot = unlist(dt[,list(lapply(V1,length))]),
                            MeanDist = unlist(dt[,list(lapply(V1,mean))]),
                            MinDist = unlist(dt[,list(lapply(V1,min))]), 
                            MaxDist = unlist(dt[,list(lapply(V1,max))])
  )
  colnames(CollectDist)<-c("SpotId", "NumbSpot", "MeanDist", "MinDist", "MaxDist")
  CollectDist
}

#Calculate distance between each spot and a list of spots
CalculateMinSpotSpotDistance<-function(SpotPosition, RefPos, GetIndex, vDelta){
  ds2 <- nearest.dist(x=SpotPosition, y=RefPos, method = "euclidean", delta=vDelta, upper=NULL)
  MinDistance <- rowMins(as.matrix(ds2))
  if (GetIndex){
    MinDistanceIndex <-  as.matrix(apply( ds2, 1, which.min)) #I need the index as well to get features of the closest RefPos spot  
  }else{
    MinDistanceIndex <- "" 
  }
  list(Dist=MinDistance,Ind=MinDistanceIndex) 
}

#Calculate distance from each surface = set a high delta
CalculateSpotSurfaceDistance<-function(SpotPosition, SurfacePosition, SurfaceFeat){
  #Note: Here delta is set a high value 100: I want all the distances
  #STEP1 : Calculate the distance to center of masses of nucleus, of the nucleolus and all the chromocenters
  FirstTab<-SurfacePosition[1:2,] #this table contains postions for Nucleus and nucleolus only
  CC<-SurfacePosition[3:nrow(SurfacePosition),] #this one for the CC
  CCF<-SurfaceFeat[3:nrow(SurfaceFeat),] #Get CC features
  meanCC<-colMeans(CC) #calculate average position for all CC
  FirstTab<-rbind(FirstTab, meanCC)
  ds <- nearest.dist(x=SpotPosition, y=FirstTab, method = "euclidean", delta=100, upper=NULL)
  ds <-data.frame(as.matrix(ds))
  #STEP2: Get minimum distance spot-CC for each spot
  Res <- CalculateMinSpotSpotDistance(SpotPosition, CC, TRUE, 5)
  ds<-cbind(ds, Res$Dist)
  FeatVect<-CCF[Res$Ind,]
  rownames(FeatVect)<-seq(1, nrow(FeatVect))
  ds<-cbind(ds, FeatVect)
  colnames(ds)<-c("NucleusCMDist", "NucleolusCMDist", "MeanCCDist", "ClosestCCCMDist",paste0(names(FeatVect), "CC"))
  ds$SpotId <- seq(1, nrow(ds)) #Since I transform spam into a matrix, I get the row = SpotId
  ds
}

clamp <- function(x, xmax) { ifelse(x<=xmax, x, xmax) }
ProcessOneImage<-function(MotherFolder, ImageIdObs,ImageIdSim ,ImageIdVert, FileName){
  n <- 3 #declaring number of segemented channels, there are as many simulated spots
  delta <- 0.2 #declaring the highest distance threshold
  #I increased delta from 0.1, because with this value on my test sample, I only get 1% spots having at least one neihboring spot
  #STEP1: READ POSITION AND INTENSITY FOR EACH FOCI
  Data1 <- data.frame(matrix(,nrow=0, ncol=3)) #will contain positions of spots
  Intensity <- data.frame(matrix(,nrow=0, ncol=n+1)) #will contain intensities of spots (n+DAPI)
  ChannelType <- c()#keeps track of spot type
  CollectSepSpots <- as.list(rep("", n))
  SegmentationIntensity <- c()#keeps track of segmentation intensity
  countSpots<-1
  for (ch in seq(0,n-1))
  {
    #Read observed spots
    #FilePath <- paste(MotherFolder, "/XTCountSpotPerShell_Result/Position_SP",ch , "_Obs_", ImageIdObs,".csv", sep="") #For linux
    FilePath <- paste(MotherFolder, "\\XTCountSpotPerShell_Result\\Position_SP",ch , "_Obs_", ImageIdObs,".csv", sep="") #For windows
    Spots <- read.csv(file=FilePath, header = TRUE, sep = ",")
    #FilePath <- paste(MotherFolder, "/XTCountSpotPerShell_Result/Intensity_SP",ch , "_Obs_", ImageIdObs,".csv", sep="") #For linux
    FilePath <- paste(MotherFolder, "\\XTCountSpotPerShell_Result\\Intensity_SP",ch , "_Obs_", ImageIdObs,".csv", sep="")#For windows
    Int <- read.csv(file=FilePath, header = TRUE, sep = ",")
    removeNullIntensity <- which(Int[paste("X", ch+1, sep="")]> 0) # keep only spots with intensity >0
    Spots <- Spots[removeNullIntensity, ]
    Data1 <- rbind(Data1, Spots)
    Int <- Int[removeNullIntensity, ]
    Intensity <- rbind(Intensity, Int)
    SegmentationIntensity<-c(SegmentationIntensity, Int[[ch+1]])
    ChannelType <- c(ChannelType, rep(ch+1,nrow(Spots))) #vector to keep track of channel type
    CollectSepSpots[[countSpots]] <- Spots
    rm(Int, Spots, removeNullIntensity) #Free RAM
    # Read simulated spots
    #FilePath <- paste(MotherFolder, "/XTSimulateRandomSpots_Result/Position_SP0_ST1" , "_Sml",ch,"_Sim_10000_", ImageIdSim,".csv", sep="") #For linux
    FilePath <- paste(MotherFolder, "\\XTSimulateRandomSpots_Result\\Position_SP0_ST1" , "_Sml",ch,"_Sim_10000_", ImageIdSim,".csv", sep="") #For windows
    Spots <- read.csv(file=FilePath, header = TRUE, sep = ",")
    Spots <- as.data.frame(t(Spots)) #see Note Plugin
    names(Spots)<-names(Data1)
    #FilePath <- paste(MotherFolder, "/XTSimulateRandomSpots_Result/Intensity_SP0_ST1" , "_Sml",ch,"_Sim_10000_", ImageIdSim,".csv", sep="") #For linux
    FilePath <- paste(MotherFolder, "\\XTSimulateRandomSpots_Result\\Intensity_SP0_ST1" , "_Sml",ch,"_Sim_10000_", ImageIdSim,".csv", sep="") #For windows
    Int <- read.csv(file=FilePath, header = TRUE, sep = ",")
    Data1 <- rbind(Data1, Spots)
    names(Int)<-names(Intensity)
    Intensity <- rbind(Intensity, Int)  
    SegmentationIntensity<-c(SegmentationIntensity,rep(0, nrow(Int)))
    ChannelType <- c(ChannelType, rep(-1*(ch+1),nrow(Spots))) #vector to keep track of simulation ID
    CollectSepSpots[[countSpots+1]] <- Spots
    rm(Int, Spots) #Free RAM
    countSpots<-countSpots+2
    print(ch)
  }
  names(Intensity) = paste0("IntensityCh", seq(0,n))
  Intensity$FociType=ChannelType
  Intensity$SegmentChannel <- SegmentationIntensity
  Intensity$SpotId <- seq(1, nrow(Intensity))
  rm(ChannelType, SegmentationIntensity)
  #STEP2: CALCULATE DISTANCE BETWEEN EACH SPOT AND SURFACES
  for (surfName in c("Nucleus","Chromocenters","Nucleolus")){
    ColName <- paste0(surfName, "Dist")
    #FilePath <- paste(MotherFolder, "/XTGetSurfaceVertices_Result/", surfName, " Vertices_",ImageIdVert,".csv", sep="") #For linux  
    FilePath <- paste(MotherFolder, "\\XTGetSurfaceVertices_Result\\", surfName, " Vertices_",ImageIdVert,".csv", sep="") #For windows  
    Vertices <- read.csv(file=FilePath, header = TRUE, sep = ",")
    res<-CalculateSpotDistance(Vertices, Data1, delta, FileName, surfName)
    Features <- paste0(surfName, c("_NumbSpot","_MeanDist", "_MinDist", "_MaxDist"))
    colnames(res) <- c("SpotId", Features)
    Intensity <- merge(Intensity, res,by = "SpotId", all.x=TRUE)
    rm(Vertices) #Free RAM  
    print(surfName)
  }
  #STEP3: CALCULATE DISTANCE BETWEEN FOCI
  Intensity$SpotId=seq(1, nrow(Data1))
  s <- c(paste0("Sp", seq(1,n)), paste0("Sim", seq(1,n)))
  for (ch in seq(1, 2*n))
  {
    res<-CalculateSpotDistance(CollectSepSpots[[ch]], Data1, delta, FileName, paste("Channel", ch,sep=""))
    Features <- paste0(s[ch], c("NumbSpot","MeanDist", "MinDist", "MaxDist"))
    colnames(res) <- c("SpotId", Features)
    Intensity <- merge(Intensity, res,by = "SpotId", all.x=TRUE)
  }
  rm(CollectSepSpots, res) #free RAM space
  #STEP4 : CALCULATE DISTANCE BETWEEN EACH SPOT AND CENTER OF MASS OF EACH SURFACE
  #FilePath <- paste(MotherFolder, "/XTSimulateRandomSpots_Result/SurfaceFeatures_",ImageIdSim ,".csv", sep="") # for linux
  FilePath <- paste(MotherFolder, "\\XTSimulateRandomSpots_Result\\SurfaceFeatures_",ImageIdSim ,".csv", sep="") #for windows
  Surface <- read.csv(file=FilePath, header = TRUE, sep = ",")
  Pos <- Surface[2:4]
  Feat <- Surface[5:9]
  SurfSPot<-CalculateSpotSurfaceDistance(Data1, Pos, Feat)
  Intensity<- merge(Intensity, SurfSPot, all=TRUE)
  rm(SurfSPot, Pos,  Surface, Data1) #Free RAM
  #STEP5: ADD COLUMNS WITH NUCLEUS AND NUCLEOLUS FEATURES
  NucleusF<-Feat[rep(1, nrow(Intensity)),] #this table contains nucleus features 
  names(NucleusF)<-paste0("N1", names(NucleusF))
  Intensity<-cbind(Intensity, NucleusF)
  NucleusF<-Feat[rep(2, nrow(Intensity)),] #this table contains nucleolus features
  names(NucleusF)<-paste0("N2", names(NucleusF))
  Intensity<-cbind(Intensity, NucleusF)
  Intensity$File<-FileName
  rm(NucleusF)
  #STEP6: REMOVE OUTLIERS
  drops <- c("SpotId","FociType", 'File') 
  IntensityRf<-Intensity[ , !(names(Intensity) %in% drops)]
  Intensity.quantiles <- apply(IntensityRf, 2, function(x, prob=0.9) { quantile(x, prob, names=F, na.rm=TRUE) })
  for (j in 1:ncol(IntensityRf)) {
    IntensityRf[,j] <- clamp(IntensityRf[,j], Intensity.quantiles[j])  
  }
  Intensity<-cbind(Intensity[ , (names(Intensity) %in% drops)], IntensityRf)
  rm(Intensity.quantiles, IntensityRf)
  #STEP7: SCALE NUMERIC SPOT FEATURES
  drops <- c("SpotId","FociType",'N1Volume','N1IntensityCh0','N1IntensityCh1', 'N1IntensityCh2', 'N1IntensityCh3', 'N2Volume','N2IntensityCh0','N2IntensityCh1', 'N2IntensityCh2','N2IntensityCh3', 'File',"Sp1MeanDist", "Sp1MinDist", "Sp1MaxDist", "Sp2MeanDist", "Sp2MinDist", "Sp2MaxDist", "Sp3MeanDist", "Sp3MinDist", "Sp3MaxDist") # these features are not to be scaled per image
  NormDf<-Intensity[ , !(names(Intensity) %in% drops)]
  NormDf<-as.data.frame(scale(NormDf))
  Intensity<-cbind(Intensity[ , (names(Intensity) %in% drops)], NormDf)
  Intensity
}
ProcessOneSimImage<-function(MotherFolder, ImageId,ImageIdVert, FileName, MaskId){
  MaskId <- 1
  SimPath <- paste(MotherFolder,"XTSimulateRandomSpots_Result", sep="\\")#for windows
  n <- 3 #declaring number of segemented channels
  delta <- 0.2 #declaring the highest distance threshold
  #I increased delta from 0.1, because with this value on my test sample, I only get 1% spots having at least one neihboring spot
  #STEP1: READ POSITION AND INTENSITY FOR EACH FOCI
  Data1 <- data.frame(matrix(,nrow=0, ncol=3)) #will contain positions of spots
  Intensity <- data.frame(matrix(,nrow=0, ncol=n+1)) #will contain intensities of spots (n+DAPI)
  SimN <- c()
  CollectSepSpots <- as.list(rep("", n))
  for (ch in seq(0,n-1)) #This loop over the different simulation (total: 5)
  {
    #FilePath <- paste(SimPath, "/Position_SP0_ST",MaskId , "_Sml",ch,"_Sim_10000_", ImageId,".csv", sep="") #For linux
    FilePath <- paste(SimPath, "\\Position_SP0_ST",MaskId , "_Sml",ch,"_Sim_100000_", ImageId,".csv", sep="") #For windows
    Spots <- read.csv(file=FilePath, header = TRUE, sep = ",")
    Spots <- as.data.frame(t(Spots))
    #FilePath <- paste(SimPath, "/Intensity_SP0_ST",MaskId , "_Sml",ch,"_Sim_10000_", ImageId,".csv", sep="") #For linux
    FilePath <- paste(SimPath, "\\Intensity_SP0_ST",MaskId , "_Sml",ch,"_Sim_100000_", ImageId,".csv", sep="")#For windows
    Int <- read.csv(file=FilePath, header = TRUE, sep = ",")
    print(paste0("Sp", dim(Spots)))
    Data1 <- rbind(Data1, Spots)
    print(paste0("Sp", dim(Int)))
    Intensity <- rbind(Intensity, Int)
    SimN <- c(SimN, rep(-1*(ch+1),nrow(Spots))) #vector to keep track of simulation ID
    CollectSepSpots[[ch+1]] <- Spots
    rm(Int, Spots) #Free RAM
    print(ch)
  }
  names(Intensity) <- paste0("IntensityCh", seq(0,n))
  Intensity$FociType <- SimN
  Intensity$SegmentChannel <- 0
  Intensity$SpotId <- seq(1, nrow(Data1))
  #STEP2: CALCULATE DISTANCE BETWEEN EACH SPOT AND SURFACES
  for (surfName in c("Nucleus","Chromocenters","Nucleolus")){
    ColName <- paste0(surfName, "Dist")
    #FilePath <- paste(MotherFolder, "/XTGetSurfaceVertices_Result/", surfName, " Vertices_",ImageIdVert,".csv", sep="") #For linux  
    FilePath <- paste(MotherFolder, "\\XTGetSurfaceVertices_Result\\", surfName, " Vertices_",ImageIdVert,".csv", sep="") #For windows  
    Vertices <- read.csv(file=FilePath, header = TRUE, sep = ",")
    res<-CalculateSpotDistance(Vertices, Data1, delta, FileName, surfName)
    Features <- paste0(surfName, c("_NumbSpot","_MeanDist", "_MinDist", "_MaxDist"))
    colnames(res) <- c("SpotId", Features)
    Intensity <- merge(Intensity, res,by = "SpotId", all.x=TRUE)
    rm(Vertices) #Free RAM  
    print(surfName)
  }
  #STEP3: CALCULATE DISTANCE BETWEEN FOCI
  for (ch in seq(1, n))
  {
    res<-CalculateSpotDistance(CollectSepSpots[[ch]], Data1, delta, FileName, paste("Channel", ch,sep=""))
    Features <- paste0(paste0("Sp", ch), c("NumbSpot","MeanDist", "MinDist", "MaxDist"))
    colnames(res) <- c("SpotId", Features)
    Intensity <- merge(Intensity, res, all=TRUE)
  }
  rm(SimN, CollectSepSpots, res) #free RAM space
  #STEP4 : CALCULATE DISTANCE BETWEEN EACH SPOT AND CENTER OF MASS OF EACH SURFACE
  #FilePath <- paste(MotherFolder, "/SurfaceFeatures_",ImageId ,".csv", sep="") # for linux
  FilePath <- paste(MotherFolder, "\\SurfaceFeatures_",ImageId ,".csv", sep="") #for windows
  Surface <- read.csv(file=FilePath, header = TRUE, sep = ",")
  Pos <- Surface[2:4]
  Feat <- Surface[5:9]
  SurfSPot<-CalculateSpotSurfaceDistance(Data1, Pos, Feat)
  Intensity<- merge(Intensity, SurfSPot, all=TRUE)
  rm(SurfSPot, Pos,  Surface, Data1) #Free RAM
  #STEP5: ADD COLUMNS WITH NUCLEUS AND NUCLEOLUS FEATURES
  NucleusF<-Feat[rep(1, nrow(Intensity)),] #this table contains nucleus features 
  names(NucleusF)<-paste0("N1", names(NucleusF))
  Intensity<-cbind(Intensity, NucleusF)
  NucleusF<-Feat[rep(2, nrow(Intensity)),] #this table contains nucleolus features
  names(NucleusF)<-paste0("N2", names(NucleusF))
  Intensity<-cbind(Intensity, NucleusF)
  Intensity$File=FileName
  rm(NucleusF)
  #STEP6: REMOVE OUTLIERS
  drops <- c("SpotId","FociType", 'File') 
  IntensityRf<-Intensity[ , !(names(Intensity) %in% drops)]
  Intensity.quantiles <- apply(IntensityRf, 2, function(x, prob=0.9) { quantile(x, prob, names=F, na.rm=TRUE) })
  for (j in 1:ncol(IntensityRf)) {
    IntensityRf[,j] <- clamp(IntensityRf[,j], Intensity.quantiles[j])  
  }
  Intensity<-cbind(Intensity[ , (names(Intensity) %in% drops)], IntensityRf)
  rm(Intensity.quantiles, IntensityRf)
  #STEP7: SCALE NUMERIC SPOT FEATURES
  drops <- c("SpotId","FociType",'N1Volume','N1IntensityCh0','N1IntensityCh1', 'N1IntensityCh2', 'N1IntensityCh3', 'N2Volume','N2IntensityCh0','N2IntensityCh1', 'N2IntensityCh2','N2IntensityCh3', 'File',"Sp1MeanDist", "Sp1MinDist", "Sp1MaxDist", "Sp2MeanDist", "Sp2MinDist", "Sp2MaxDist", "Sp3MeanDist", "Sp3MinDist", "Sp3MaxDist") # these features are not to be scaled per image
  NormDf<-Intensity[ , !(names(Intensity) %in% drops)]
  NormDf<-as.data.frame(scale(NormDf))
  Intensity<-cbind(Intensity[ , (names(Intensity) %in% drops)], NormDf)
  Intensity
}
SelectFiles <- function(FolderPath)
{
  Data1 <- data.frame(matrix(,nrow=0, ncol=2)) #will contain positions of spots
  for (SubFolder in c(1,2,3,4,5))
  {
    FilePath<- paste(MotherFolder,SubFolder,  "XTCountSpotPerShell_Result","FileName.csv" , sep="/") # for linux
    # FilePath<- paste(MotherFolder,SubFolder, "XTCountSpotPerShell_Result","FileName.csv", sep="\\") # for windows
    OBsdf <- read.csv(file=FilePath, header = TRUE, sep = ",")
    FilePath<- paste(MotherFolder,SubFolder,  "XTSimulateRandomSpots_Result","FileName.csv" , sep="/") #for linux
    # FilePath<- paste(MotherFolder,SubFolder, "XTSimulateRandomSpots_Result","FileName.csv", sep="\\") # for windows
    Simdf <- read.csv(file=FilePath, header = TRUE, sep = ",")
    df<-merge(OBsdf, Simdf, all = FALSE)
    df$FolderID<-SubFolder
    Data1<-rbind(Data1, df)
  }
  names(Data1)<-c("FileName", "FolderID")
  FilePath<- paste(MotherFolder,"XTNucleiPloidy_Result", "XTNucleiPloidy_Result.csv" , sep="/")
  PLoidyDf <- read.csv(file=FilePath, header = TRUE, sep = ",")
  Data1<-merge(Data1, PLoidyDf, all = FALSE)
  Data1<-Data1[, names(Data1)%in%c("FileName","FolderID", "ClusterId")]
  Data1
}


ProcessOneFOlder<-function(FolderPath){
  #for (highFol in c(1,3,4,5))
  for (highFol in c(2))
  {
    #path in Linux **********
    #MotherFolder<- paste(FolderPath,highFol, "XTCountSpotPerShell_Result", sep="/")
    #MotherFolder<- paste(FolderPath,highFol, "XTSimulateRandomSpots_Result", sep="/") #for sim data
    #FilePath <- paste(MotherFolder, "FileName.csv", sep="/")
    #path in Windows ********
    MotherFolder<- paste(FolderPath,highFol, "XTCountSpotPerShell_Result", sep="\\")
    FilePath <- paste(MotherFolder, "FileName.csv", sep="\\")
    #************************
    FileName <- read.csv(file=FilePath, header = TRUE, sep = ",")
    NumberImage <- length(FileName$X0)
    Result <- data.frame(matrix(,nrow=0, ncol=40)) #will contain intensities of spots (n+DAPI)
    for(FileId in seq(1, NumberImage)){
      # for(MaskId in c(0,1)){ #for sim data
      for(MaskId in c(-1)){ #for sim data
        FileN <-FileName$X0[FileId]
        # res<-ProcessOneSimImage( MotherFolder, FileId, FileN, MaskId) #sim
        res<-ProcessOneImage( MotherFolder, FileId, FileN) #obs
        res$Mask=MaskId #sim
        Result<-rbind(Result, res)
      }
    }
    rm(res)
    FilePath=paste(MotherFolder, "MLResult.csv", sep="/")
    write.csv(Result, file = FilePath, quote=FALSE, row.names = FALSE) 
    rm(Result)
  }
}

MotherFolder= "Z:\\Result0309\\s20_Wt" #windows
#MotherFolder= "/home/nani/smb4k/BOTSERV4.UZH.CH/gr_ug_ext/mashen2/Result0309/s20_Wt" #linux
ProcessOneFOlder(MotherFolder)


#STEP0: TEST ONE IMAGE PROCESSING
##for linux
#MotherFolder= "/home/nani/smb4k/BOTSERV4.UZH.CH/gr_ug_ext/mashen2/Result0309/s20_Wt/1"
##for windows
MotherFolder= "Z:\\Result0309\\s20_Wt"
highFol <- "1"
MotherFolder <- paste(MotherFolder, highFol, sep="\\")#for windows
res <- ProcessOneSimImage(MotherFolder, 1, "FileName", 1)

#STEP1: CREATE FILE CONTAINING NAMES OF PROCESSED IMAGE AND PLOIDY ID
##for linux
# MotherFolder= "/home/nani/smb4k/BOTSERV4.UZH.CH/gr_ug_ext/mashen2/Result0309/s20_Wt"
##for windows
MotherFolder= "Z:\\Result0309\\s20_Wt"
FileName<-SelectFiles(FolderPath)
FilePath=paste(MotherFolder, "FileNames.csv", sep="/")
write.csv(FileName, file = FilePath, quote=FALSE, row.names = FALSE)
 
#STEP2: GET SEGEMENTATION STATUS FOR EACH PROCESSED IMAGE
#segmentation status is checked manually using the snapshots produced by the plugins
##for linux
# MotherFolder<- "/home/nani/smb4k/BOTSERV4.UZH.CH/gr_ug_ext/mashen2/Result0309/s20_Wt"
##for windows
MotherFolder= "Z:\\Result0309\\s20_Wt"
FilePath=paste(MotherFolder, "FileNames.csv", sep="/") #on linux
#FilePath=paste(MotherFolder, "FileNames.csv", sep="\\") #on windows
SelectedIMages <- read.csv(file=FilePath, header = TRUE, sep = ",")
View(SelectedIMages)
table(SelectedIMages$ClusterId)
SelectedIMages<-SelectedIMages[SelectedIMages$Nucleolusok==1, names(SelectedIMages)%in%c("FileName", "ClusterId", "FolderID")]
table(SelectedIMages$ClusterId)

