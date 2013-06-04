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

HPLC.discrete.process<- function(data,wavelength,time,flow)
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
  

  # Separate the three wavelengths and make them relative to highest peak.
  wavA<-absorbance[seq(points[1],points[2],1)]
  
  wav1<-wavA/max(wavA)
  
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
  col.names<-c("Volume (ml)",paste(Sample,wavelength[1],sep=": "))
  row.names<-seq(1,length(elution),1)
  dimnames<-list(row.names,col.names)
  data<-matrix(data=c(elution,wav1),ncol=2,dimnames=dimnames)
}







HPLC.discrete.get.peaks<-function(data)
{
# Fits linear model to baseline of 24mer and monomer
x.values<-unlist(array(data[,1]))
y.values<-unlist(array(data[,2]))

plot(x=x.values,y=y.values,type='l',col=2)

largemer.baseline.x.val<-x.values[grep(x=x.values<=4,pattern=TRUE)]
largemer.baseline.y.val<-y.values[grep(x=x.values<=4,pattern=TRUE)]

monomer.baseline.x.val<-x.values[grep(x=x.values>=8,pattern=TRUE)]
monomer.baseline.y.val<-y.values[grep(x=x.values>=8,pattern=TRUE)]

c.largemer<-array(summary(lm(largemer.baseline.y.val~largemer.baseline.x.val))$coefficients[,1][1])
c.monomer<-array(summary(lm(monomer.baseline.y.val~monomer.baseline.x.val))$coefficients[,1][1])

m.largemer<-array(summary(lm(largemer.baseline.y.val~largemer.baseline.x.val))$coefficients[,1][2])
m.monomer<-array(summary(lm(monomer.baseline.y.val~monomer.baseline.x.val))$coefficients[,1][2])

largemer.linear<-lm(largemer.baseline.y.val~largemer.baseline.x.val)
abline(largemer.linear,col=4)
monomer.linear<-lm(monomer.baseline.y.val~monomer.baseline.x.val)
abline(monomer.linear,col=3)

# Normalise 24mer peak relative to new baseline 
x.values.largemer<-x.values[grep(x=(x.values<=6),pattern=TRUE)]
y.values.largemer<-y.values[grep(x=(x.values<=6),pattern=TRUE)]
x.output.largemer<-x.values.largemer[rev(order(y.values.largemer))[1]]
y.output.largemer<-y.values.largemer[rev(order(y.values.largemer))[1]]
y.baseline.largemer<-m.largemer*x.output.largemer+c.largemer
normalised.ymax.largemer<-y.output.largemer-y.baseline.largemer
print(c("non-normalised peak = " , y.output.largemer))
print(c("largemer y baseline = " , y.baseline.largemer))
print(c("normalised peak height = " , normalised.ymax.largemer))



# Normalise Monomer peak relative to new baseline
x.values.monomer<-x.values[grep(x=(x.values>=6),pattern=TRUE)]
y.values.monomer<-y.values[grep(x=(x.values>=6),pattern=TRUE)]
x.output.monomer<-x.values.monomer[rev(order(y.values.monomer))[1]]
y.output.monomer<-y.values.monomer[rev(order(y.values.monomer))[1]]
y.baseline.monomer<-m.monomer*x.output.monomer+c.monomer
normalised.ymax.monomer<-y.output.monomer-y.baseline.monomer
print(c("non-normalised peak = " , y.output.monomer))
print(c("monomer y baseline = " , y.baseline.monomer))
print(c("normalised peak height = " , normalised.ymax.monomer))

colnames<-c('Normalised 24mer Peak','Normalised Monomer Peak')
rownames<-'Relative Absorbance'
dimnames<-list(rownames,colnames)
normalised.peaks<-matrix(data=c(normalised.ymax.largemer,normalised.ymax.monomer),ncol=2,dimnames=dimnames)
tot.peaks<-peaks[,2]+peaks[,1]

percent.24mer<-(normalised.ymax.largemer/(normalised.ymax.largemer+normalised.ymax.monomer))*100
percent.monomer<-(normalised.ymax.monomer/(normalised.ymax.largemer+normalised.ymax.monomer))*100
colnames<-c('24mer Peak Absorbance (%)','Monomer Peak Absorbance (%)')
rownames<-'1'
dimnames<-list(rownames,colnames)
normalised.peaks<-matrix(data=c(percent.24mer,percent.monomer),ncol=2,dimnames=dimnames)
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
  png(file=paste(c(Sample,".png"),collapse=""), bg="transparent", width =1000, height=500,units="px",pointsize=13)
  plot(x=data[,1],y=data[,2],type='l',col=colours[1],xlab="Volume (ml)",ylab="Relative Absorbance (A.U)",main=Sample,lwd=2,ylim=c(-0.1,1),cex.main=3,cex.lab=1.2,cex.axis=1.2)
  points(x=data[,1],y=data[,3],type='l',col=colours[2],lwd=2)
  points(x=data[,1],y=data[,4],type='l',col=colours[3])
  points(x=data[,1],y=data[,5],type='l',col=colours[4])
  abline(-1.380262e-06,-0.009661836)

  legend(x=c(9,11),y=c(0.8,1),legend=c(figure.leg1,figure.leg2,figure.leg3,figure.leg4),col=colours,lty=1,lwd=3,cex=1.1, bty='n')
  dev.off()
}
############################################################



# Gets values for largest 24mer and monomer peak


# EXAMPLE

# Preliminary steps
# For each of your samples, the HPLC software creates three data files. Two of these data files will be in .asc format.
# Of the two files in .asc format, the one with the shorter name contains the discrete data.
# Add the word 'Discrete' to the filename of each of the files that contain the discrete dataset.
#  e.g. 'example.asc' would become 'example - Discrete.asc'.

################# DEFINE INPUT VARIABLES

# HPLC.import
dir<-dir()
Sample.name<-"C.69" 
result.type<-"discrete"
parent.dir<-"//ic.ac.uk/homes/hfr10/2013-05-23"
figure.dir<-"//ic.ac.uk/homes/hfr10/"
wavelength<-c("280 nm")
time<-45
flow<-0.3

sample.file.names<-grep(paste(Sample.name," Run ",sep=""),dir,value=TRUE)
sample.file.names<-grep("discrete",sample.file.names,value=TRUE)
sample.numbers<-gsub(paste(Sample.name," Run ",sep=""),"",sample.file.names)
sample.numbers<-as.numeric(gsub(" - discrete.asc","",sample.numbers))
i<-1
while(i<=max(sample.numbers))
{
# HPLC.discrete.process
################# Calculate Results
Sample<-paste(Sample.name,"Run",sample.numbers[i],sep=" ")
data<-HPLC.import(Sample,result.type,parent.dir)
data<-HPLC.discrete.process(data=data,wavelength=wavelength,time=time,flow=flow)
peaks<-c(peaks,HPLC.discrete.get.peaks(data))
print(peaks)
i<-i+1
}
################# Calculate Results

HPLC.discrete.plot1(data,figure.dir,colours)
HPLC.discrete.summary()
#rm(list=ls())
