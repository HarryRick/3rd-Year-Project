require(plotrix)

#################### FUNCTION HPLC.import ####################
# Input:
  # 1) Sample ID (e.g. C.29)
	# 2) Data type (discrete or spectrum)
	# 3) Parent directory of datafiles
# Output = read in all data into R
# Development = this will break down if >1 matching files are found...

error.bar <- function(x, y, upper, lower=upper, length=0.1,...)
{
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
	stop("vectors must be same length")
	arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}




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







HPLC.discrete.get.peaks<-function(data,Sample)
{
# Fits linear model to baseline of 24mer and monomer
x.values<-unlist(array(data[,1]))
y.values<-unlist(array(data[,2]))

png(file=paste(c(Sample,".png"),collapse=""), bg="transparent", width =1000, height=500,units="px",pointsize=13)

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
abline(largemer.linear,col=4,lty=2)
monomer.linear<-lm(monomer.baseline.y.val~monomer.baseline.x.val)
abline(monomer.linear,col=3,lty=2)
dev.off()

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

percent.24mer<-(normalised.ymax.largemer/(normalised.ymax.largemer+normalised.ymax.monomer))*100
percent.monomer<-(normalised.ymax.monomer/(normalised.ymax.largemer+normalised.ymax.monomer))*100

normalised.peaks<-c(percent.24mer,percent.monomer)
}

############################################################


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
Sample.name<-"C.42" 
result.type<-"discrete"
parent.dir<-"//ic.ac.uk/homes/hfr10/2013-05-23"
figure.dir<-"//ic.ac.uk/homes/hfr10/"
wavelength<-c("280 nm")
time<-45
flow<-0.3

#Extracts the repeat numbers from the directory 
sample.file.names<-grep(paste(Sample.name," Run ",sep=""),dir,value=TRUE)
sample.file.names<-grep("discrete",sample.file.names,value=TRUE)
sample.numbers<-gsub(paste(Sample.name," Run ",sep=""),"",sample.file.names)
sample.numbers<-as.numeric(gsub(" - discrete.asc","",sample.numbers))

# Creates matrix which will contain peak data.
colnames<-c(paste(Sample.name,"24mer % absorbance"),paste(Sample.name,"Monomer % Abosrbance"))
rownames<-character()
i<-1
while(i<=length(sample.numbers))
{
	rownames[i]<-paste("Run",i)
	i<-i+1
}
dimnames<-list(rownames,colnames)
peaks.store<-matrix(ncol=2,nrow=length(sample.numbers),dimnames=dimnames)

i<-1
while(i<=max(sample.numbers))
{
	# HPLC.discrete.process
	################# Calculate Results
	Sample<-paste(Sample.name,"Run",sample.numbers[i],sep=" ")
	data<-HPLC.import(Sample,result.type,parent.dir)
	data<-HPLC.discrete.process(data=data,wavelength=wavelength,time=time,flow=flow)
	normalised.peaks<-HPLC.discrete.get.peaks(data,Sample)
	peaks.store[i,]<-normalised.peaks
	print(normalised.peaks)
	if(Sample.name=="C.42")
	{
		C.42.peaks.store<-peaks.store
	}
	if(Sample.name=="C.68")
	{
		C.68.peaks.store<-peaks.store
	}
	if(Sample.name=="C.69")
	{
		C.69.peaks.store<-peaks.store
	}
	i<-i+1
}
write.csv(x=peaks.store,file=paste(Sample.name,"peak percentages.csv"))

C.42.24mer.mean<-mean(C.42.peaks.store[,1])
C.42.monomer.mean<-mean(C.42.peaks.store[,2])
C.68.monomer.mean<-mean(C.68.peaks.store[,2])
C.68.24mer.mean<-mean(C.68.peaks.store[,1])
C.69.24mer.mean<-mean(C.69.peaks.store[,1])
C.69.monomer.mean<-mean(C.69.peaks.store[,2])

C.42.24mer.std.error<-std.error(C.42.peaks.store[,1])
C.42.monomer.std.error<-std.error(C.42.peaks.store[,2])
C.68.monomer.std.error<-std.error(C.68.peaks.store[,2])
C.68.24mer.std.error<-std.error(C.68.peaks.store[,1])
C.69.24mer.std.error<-std.error(C.69.peaks.store[,1])
C.69.monomer.std.error<-std.error(C.69.peaks.store[,2])

barplot.data<-c(C.42.24mer.mean, C.42.monomer.mean,C.68.24mer.mean,C.68.monomer.mean,
C.69.24mer.mean, C.69.monomer.mean)
error.bar.data<-c(C.42.24mer.std.error, C.42.monomer.std.error,C.68.24mer.std.error,C.68.monomer.std.error,
C.69.24mer.std.error, C.69.monomer.std.error)

colnames<-c("Wildtype GLFG","F36R GLFG","L162R GLFG")
rownames<-c("24mer percentage absorption","Monomer percentage absorption")
rownames.error<-c("24mer percentage absorption std error","Monomer percentage absorption std error")

dimnames<-list(rownames,colnames)
dimnames.error<-list(rownames.error,colnames)
barplot.data.matrix<-matrix(ncol=3,nrow=2,barplot.data,byrow=FALSE,dimnames=dimnames)
error.bar.data.matrix<-matrix(ncol=3,nrow=2,error.bar.data,byrow=FALSE,dimnames=dimnames.error)

png(file=paste(c("Barchart of Percentage aborbance values.png"),
collapse=""), bg="transparent", width =1000, height=500,units="px",pointsize=13)
barplot<-barplot(height=barplot.data.matrix,main="Effect of mutants on the ratio of 24mer against the monomer",
beside=TRUE,space=c(0,1),col=c("palegreen","darkred"),
legend.text=TRUE,ylab="Percentage absorbance",ylim=c(0,120))
error.bar(barplot,barplot.data,error.bar.data)
dev.off()

