#################### FUNCTION HPLC.import ####################
# Input:
  # 1) Sample ID (e.g. C.29)
	# 2) Data type (discrete or spectrum)
	# 3) Parent directory of datafiles
# Output = read in all data into R
# Development = this will break down if >1 matching files are found... 


HPLC.import <- function (Sample,result.type,parent.dir)
{
# Navigate to data
setwd(parent.dir)
data.file<-dir()[grep(x=dir(parent.dir),pattern=result.type,ignore.case=TRUE)]
data.file<-data.file[grep(x=data.file,pattern=Sample)]

data.file.names<-unlist(strsplit(x=data.file,split="-"))[seq(1,(3*length(data.file)),3)]

# Read data into R
data<-read.csv(file=data.file,header=FALSE)

}


############################################################



#################### FUNCTION HPLC.discrete.process ####################

# Input = output of HPLC.import()
# Output = list of absorbance values for each of the three wavelengths analysed

HPLC.discrete.process<- function(data,control.data,wavelength,time,flow)

{

# Determine number of datapoints for each discrete wavelength
  points<-as.numeric(array(unlist(data[,1][16:17])))
  points<-c(
    1,          points[1],
    points[1]+1,            points[1]*2,
    points[1]*2+1,            points[1]*3)
  
    
  # Parse wavelength data
  absorbance<-as.numeric(array(unlist(data[,1][-seq(1,30,1)])))
  absorbance<-absorbance[1:points[6]]
  
  # Determine number of datapoints for each discrete wavelength
  control.points<-as.numeric(array(unlist(control.data[,1][16:17])))
  control.points<-c(
    1,          control.points[1],
    control.points[1]+1,            control.points[1]*2,
    control.points[1]*2+1,            control.points[1]*3)
  
    
  # Parse wavelength data
  control.absorbance<-as.numeric(array(unlist(control.data[,1][-seq(1,30,1)])))
  control.absorbance<-control.absorbance[1:control.points[6]]

    
  # Separate the three wavelengths and make them relative to highest peak.
  wavA<-absorbance[seq(points[1],points[2],1)]
  wavB<-absorbance[seq(points[5],points[6],1)]
  control.wavA<-control.absorbance[seq(points[1],points[2],1)]
  control.wavB<-control.absorbance[seq(points[5],points[6],1)]  
  
  int.wavA<-sum(wavA)
  int.wavB<-sum(wavB)
  int.control.wavA<-sum(control.wavA)
  int.control.wavB<-sum(control.wavB)
  
  adj.wavA<-wavA/int.wavA
  adj.wavB<-wavB/int.wavB
  adj.control.wavA<-control.wavA/int.control.wavA
  adj.control.wavB<-control.wavB/int.control.wavB
  
  max<-max(wavA,wavB,control.wavA,control.wavB)
  adj.max<-max(adj.wavA,adj.wavB,adj.control.wavA,adj.control.wavB)
  
  wav1<-wavA/max
  wav2<-wavB/max
  control.wav1<-control.wavA/max
  control.wav2<-control.wavB/max
  
	# Generate X-axis by calculating elution volume
  
  # Sampling rate = 10x per second
  
  # Flow rate = 0.3 ml / minute
  
  # Experiment time = 45 minutes
  
  # Total Volume = 45 * 0.3
  
  # Total Volume / Total Readings
  
  
  # Elution point (ml) = Time (min) * ( Volume (ml) / Time (min) )
  total.volume<-time * flow
  elution<-  total.volume / (length(wav1)-1)
  elution<- seq (0, (total.volume) , elution)
  
  
  # Order results into a matrix
  col.names<-c("Volume (ml)",paste(Sample,wavelength[1],sep=": "),paste(Sample,wavelength[2],sep=": "),
  paste(Control,wavelength[1],sep=": "),paste(Control,wavelength[2],sep=": "))
  row.names<-seq(1,length(elution),1)
  dimnames<-list(row.names,col.names)
  data<-matrix(data=c(elution,wav1,wav2,control.wav1,control.wav2),ncol=5,dimnames=dimnames)
  
  

}

############################################################



#################### FUNCTION HPLC.discrete.summary ####################

# Finds and saves elution point of greatest peak

HPLC.discrete.summary<-function()
{
# Peak
setwd(parent.dir)
wav1.peak<-data[,1][rev(order(data[,2]))[1]]
wav2.peak<-data[,1][rev(order(data[,3]))[1]]
control.wav1.peak<-data[,1][rev(order(data[,4]))[1]]
control.wav2.peak<-data[,1][rev(order(data[,4]))[1]]
row.name<-c(paste(Sample, "elution peak (ml)"),paste(Control, "elution peak(ml)"))
dimnames<-list(row.name,wavelength)
sample.summary<-matrix(c(wav1.peak,wav2.peak,control.wav1.peak,control.wav2.peak),
ncol=2,nrow=2,dimnames=dimnames)
file.name<-paste(c(Sample,"-summary.csv"),collapse="")
write.csv(x=sample.summary,file=file.name)
}
############################################################


#################### FUNCTION HPLC.discrete.plot1 ####################

# Input = output of HPLC.import()
# Output = list of absorbance values for each of the three wavelengths analysed

HPLC.discrete.plot1<- function(hplc.data,figure.dir,colours)
{
  setwd(figure.dir)
  png(file=paste(c(Sample,".png"),collapse=""), bg="transparent", width =1000, height=500,units="px",pointsize=13)
  plot(x=data[,1],y=data[,2],type='l',col=colours[1],xlab="Volume (ml)",ylab="Relative Absorbance (A.U)",main=Sample,lwd=2,ylim=c(-0.1,1),cex.main=3,cex.lab=1.2,cex.axis=1.2)
  points(x=data[,1],y=data[,3],type='l',col=colours[2],lwd=2)
  points(x=data[,1],y=data[,4],type='l',col=colours[3])
  points(x=data[,1],y=data[,5],type='l',col=colours[4])

  legend(x=c(9,11),y=c(0.8,1),legend=c(figure.leg1,figure.leg2,figure.leg3,figure.leg4),col=colours,lty=1,lwd=3,cex=1.1, bty='n')
  dev.off()
}
############################################################



# Gets values for largest 24mer and monomer peak
max.wav1<-max(data[,2]
max.wav2<-max(data[,3])
max.control.wav1<-max(data[,4])
max.control.wav2<-max(data[,5])
monomer.wav1<-max(data[12300:length(data[,2]),2])
monomer.wav2<-max(data[12300:length(data[,3]),3])
monomer.control.wav1<-max(data[12300:length(data[,4]),4])
monomer.control.wav2<-max(data[12300:length(data[,5]),5])


# EXAMPLE

# Preliminary steps
# For each of your samples, the HPLC software creates three data files. Two of these data files will be in .asc format.
# Of the two files in .asc format, the one with the shorter name contains the discrete data.
# Add the word 'Discrete' to the filename of each of the files that contain the discrete dataset.
#  e.g. 'example.asc' would become 'example - Discrete.asc'.

################# DEFINE INPUT VARIABLES

# HPLC.import
Sample.name<-"C.68 - Small - 5nm GNPs" 
Sample<-Sample.name
result.type<-"spectrum"
parent.dir<-"/Users/harryrick/Dropbox/Work/Imperial/3rd Year/Nanocage Project/Harry/Data/2013-05-29"

Control<-"C.68 - Small - Control"

# HPLC.discrete.process
wavelength<-c("280 nm","530 nm")
time<-45
flow<-0.3

################# Calculate Results
data<-HPLC.import(Sample,result.type,parent.dir)
control.data<-HPLC.import(Control,result.type,parent.dir)
data<-HPLC.discrete.process(data,control.data,wavelength,time,flow)

figure.leg1<-paste(Sample,wavelength[1],sep=": ")
figure.leg2<-paste(Sample,wavelength[2],sep=": ")

# HPLC.discrete.plot1 
figure.dir<-"/Users/harryrick/Dropbox/Work/Imperial/3rd Year/Nanocage Project/Harry/Data/"
colours<-c(3,4,5,6)



################# Calculate Results

figure.leg3<-paste(Control,wavelength[1],sep=": ")
figure.leg4<-paste(Control,wavelength[2],sep=": ")

HPLC.discrete.plot1(data,figure.dir,colours)
HPLC.discrete.summary()
#rm(list=ls())

linear.x<-data[18000:length(data[,1]),1]
linear.y<-data[18000:length(data[,1]),5]
