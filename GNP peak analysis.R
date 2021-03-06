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




HPLC.import <- function (Sample.name,result.type,parent.dir,peak.size,treatment,dir)
{
# Navigate to data
	data.file<-grep(Sample.name,dir,value=TRUE)
	data.file<-grep(peak.size,data.file,value=TRUE)
	data.file<-grep(result.type,data.file,value=TRUE)
	data.file<-grep(treatment,data.file,value=TRUE)

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
  if(wavelength=="280nm")
  {
  	wavA<-absorbance[seq(points[1],points[2],1)]
  }
  if(wavelength=="497nm")
  {
  	wavA<-absorbance[seq(points[3],points[4],1)]
  }
  if(wavelength=="530nm")
  {
  	wavA<-absorbance[seq(points[5],points[6],1)]
  }
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

png(file=paste(c(Sample, peak.size,treatment,".png"),collapse=" "), bg="transparent", width =1000, height=500,units="px",pointsize=13)

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
x.values.monomer<-x.values[grep(x=(x.values>=6 & x.values<=8),pattern=TRUE)]
y.values.monomer<-y.values[grep(x=(x.values>=6 & x.values<=8),pattern=TRUE)]
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
parent.dir<-"//ic.ac.uk/homes/hfr10/2013-05-29"
dir<-dir(parent.dir)
Sample.name<-"C.69"
result.type<-"discrete"
peak.size<-"Small"
figure.dir<-"//ic.ac.uk/homes/hfr10/"
wavelength<-"497nm"
time<-45
flow<-0.3

#Extracts the repeat numbers from the directory 
data.files<-grep(Sample.name,dir,value=TRUE)
data.files<-grep(result.type,data.files,value=TRUE)
data.files<-grep(peak.size,data.files,value=TRUE)
data.files<-grep(treatment,data.files,value=TRUE)

sample.numbers<-gsub(paste(Sample.name," Run ",sep=""),"",data.files)
sample.numbers<-as.numeric(gsub(pattern=paste(" - ",peak.size," - ",treatment, " - discrete.asc",sep=""),
"",sample.numbers))

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
treatment<-"Control"
while(i<=max(sample.numbers))
{
	# HPLC.discrete.process
	################# Calculate Results
	Sample<-paste(Sample.name,"Run",sample.numbers[i],sep=" ")
	data<-HPLC.import(Sample,result.type,parent.dir,peak.size,treatment,dir)
	data<-HPLC.discrete.process(data=data,wavelength=wavelength,time=time,flow=flow)
	normalised.peaks<-HPLC.discrete.get.peaks(data,Sample)
	peaks.store[i,]<-normalised.peaks
	print(normalised.peaks)
	if(Sample.name=="C.42")
	{
		if(peak.size=="Large")
		{
			if(treatment=="Control")
			{	
				cont.L.C42.peaks<-peaks.store
			}
			else
			{
				GNP.L.C42.peaks<-peaks.store
			}
		}
		else
		{
			if(treatment=="Control")
			{	
				cont.S.C42.peaks<-peaks.store
			}
			else
			{
				GNP.S.C42.peaks<-peaks.store
			}
		}
	}
	if(Sample.name=="C.68")
	{
		if(peak.size=="Large")
		{
			if(treatment=="Control")
			{	
				cont.L.C68.peaks<-peaks.store
			}
			else
			{
				GNP.L.C68.peaks<-peaks.store
			}	
		}
		else
		{
			if(treatment=="Control")
			{	
				cont.S.C68.peaks<-peaks.store
			}
			else
			{
				GNP.S.C68.peaks<-peaks.store
			}		
		}	
	}
	if(Sample.name=="C.69")
	{	
		if(peak.size=="Large")
		{
			if(treatment=="Control")
			{	
				cont.L.C69.peaks<-peaks.store
			}
			else
			{
				GNP.L.C69.peaks<-peaks.store
			}		
		}
		else
		{
			if(treatment=="Control")
			{	
				cont.S.C69.peaks<-peaks.store
			}
			else
			{
				GNP.S.C69.peaks<-peaks.store
			}		
		}
	}
	i<-i+1
	if(treatment=="Control")
	{
		if(i>=max(sample.numbers))
		{
			treatment<-"5nm GNPs"
			i<-1
		}
	}
	p<-p+1
}
write.csv(x=peaks.store,file=paste(Sample.name,"peak percentages.csv"))

#Calculates means of all the data
cont.S.C.42.24mer.mean<-mean(cont.S.C42.peaks[,1])
cont.S.C.42.monomer.mean<-mean(cont.S.C42.peaks[,2])
cont.S.C.68.monomer.mean<-mean(cont.S.C68.peaks[,2])
cont.S.C.68.24mer.mean<-mean(cont.S.C68.peaks[,1])
cont.S.C.69.24mer.mean<-mean(cont.S.C69.peaks[,1])
cont.S.C.69.monomer.mean<-mean(cont.S.C69.peaks[,2])

GNP.S.C.42.24mer.mean<-mean(GNP.S.C42.peaks[,1])
GNP.S.C.42.monomer.mean<-mean(GNP.S.C42.peaks[,2])
GNP.S.C.68.monomer.mean<-mean(GNP.S.C68.peaks[,2])
GNP.S.C.68.24mer.mean<-mean(GNP.S.C68.peaks[,1])
GNP.S.C.69.24mer.mean<-mean(GNP.S.C69.peaks[,1])
GNP.S.C.69.monomer.mean<-mean(GNP.S.C69.peaks[,2])

cont.L.C.42.24mer.mean<-mean(cont.L.C42.peaks[,1])
cont.L.C.42.monomer.mean<-mean(cont.L.C42.peaks[,2])
cont.L.C.68.monomer.mean<-mean(cont.L.C68.peaks[,2])
cont.L.C.68.24mer.mean<-mean(cont.L.C68.peaks[,1])
cont.L.C.69.24mer.mean<-mean(cont.L.C69.peaks[,1])
cont.L.C.69.monomer.mean<-mean(cont.L.C69.peaks[,2])

GNP.L.C.42.24mer.mean<-mean(GNP.L.C42.peaks[,1])
GNP.L.C.42.monomer.mean<-mean(GNP.L.C42.peaks[,2])
GNP.L.C.68.monomer.mean<-mean(GNP.L.C68.peaks[,2])
GNP.L.C.68.24mer.mean<-mean(GNP.L.C68.peaks[,1])
GNP.L.C.69.24mer.mean<-mean(GNP.L.C69.peaks[,1])
GNP.L.C.69.monomer.mean<-mean(GNP.L.C69.peaks[,2])


# Calculates standard errors of all data
cont.S.C.42.24mer.std.error<-std.error(cont.S.C42.peaks[,1])
cont.S.C.42.monomer.std.error<-std.error(cont.S.C42.peaks[,2])
cont.S.C.68.monomer.std.error<-std.error(cont.S.C68.peaks[,2])
cont.S.C.68.24mer.std.error<-std.error(cont.S.C68.peaks[,1])
cont.S.C.69.24mer.std.error<-std.error(cont.S.C69.peaks[,1])
cont.S.C.69.monomer.std.error<-std.error(cont.S.C69.peaks[,2])

GNP.S.C.42.24mer.std.error<-std.error(GNP.S.C42.peaks[,1])
GNP.S.C.42.monomer.std.error<-std.error(GNP.S.C42.peaks[,2])
GNP.S.C.68.monomer.std.error<-std.error(GNP.S.C68.peaks[,2])
GNP.S.C.68.24mer.std.error<-std.error(GNP.S.C68.peaks[,1])
GNP.S.C.69.24mer.std.error<-std.error(GNP.S.C69.peaks[,1])
GNP.S.C.69.monomer.std.error<-std.error(GNP.S.C69.peaks[,2])

cont.L.C.42.24mer.std.error<-std.error(cont.L.C42.peaks[,1])
cont.L.C.42.monomer.std.error<-std.error(cont.L.C42.peaks[,2])
cont.L.C.68.monomer.std.error<-std.error(cont.L.C68.peaks[,2])
cont.L.C.68.24mer.std.error<-std.error(cont.L.C68.peaks[,1])
cont.L.C.69.24mer.std.error<-std.error(cont.L.C69.peaks[,1])
cont.L.C.69.monomer.std.error<-std.error(cont.L.C69.peaks[,2])

GNP.L.C.42.24mer.std.error<-std.error(GNP.L.C42.peaks[,1])
GNP.L.C.42.monomer.std.error<-std.error(GNP.L.C42.peaks[,2])
GNP.L.C.68.monomer.std.error<-std.error(GNP.L.C68.peaks[,2])
GNP.L.C.68.24mer.std.error<-std.error(GNP.L.C68.peaks[,1])
GNP.L.C.69.24mer.std.error<-std.error(GNP.L.C69.peaks[,1])
GNP.L.C.69.monomer.std.error<-std.error(GNP.L.C69.peaks[,2])


C.42.barplot.data<-c(cont.L.C.42.24mer.mean, cont.L.C.42.monomer.mean, GNP.L.C.42.24mer.mean,
GNP.L.C.42.monomer.mean, cont.S.C.42.24mer.mean, cont.S.C.42.monomer.mean, GNP.S.C.42.24mer.mean,
GNP.S.C.42.monomer.mean)

C.68.barplot.data<-c(cont.L.C.68.24mer.mean, cont.L.C.68.monomer.mean, GNP.L.C.68.24mer.mean,
GNP.L.C.68.monomer.mean, cont.S.C.68.24mer.mean, cont.S.C.68.monomer.mean, GNP.S.C.68.24mer.mean,
GNP.S.C.68.monomer.mean)

C.69.barplot.data<-c(cont.L.C.69.24mer.mean, cont.L.C.69.monomer.mean, GNP.L.C.69.24mer.mean,
GNP.L.C.69.monomer.mean, cont.S.C.69.24mer.mean, cont.S.C.69.monomer.mean, GNP.S.C.69.24mer.mean,
GNP.S.C.69.monomer.mean)

error.bar.data<-c(C.42.24mer.std.error, C.42.monomer.std.error,C.68.24mer.std.error,C.68.monomer.std.error,
C.69.24mer.std.error, C.69.monomer.std.error)

colnames<-c("Large fraction","Small Fraction")
rownames<-c("24mer without GNPs","Monomer without GNPs","24mer with GNPs","Monomer with GNPs")
rownames.error<-c("24mer percentage absorption std error","Monomer percentage absorption std error")

dimnames<-list(rownames,colnames)
dimnames.error<-list(rownames.error,colnames)

C.42.barplot.data.matrix<-matrix(ncol=2,nrow=4,C.42.barplot.data,byrow=FALSE,dimnames=dimnames)
error.bar.data.matrix<-matrix(ncol=3,nrow=2,error.bar.data,byrow=FALSE,dimnames=dimnames.error)

C.68.barplot.data.matrix<-matrix(ncol=2,nrow=4,C.68.barplot.data,byrow=FALSE,dimnames=dimnames)
error.bar.data.matrix<-matrix(ncol=3,nrow=2,error.bar.data,byrow=FALSE,dimnames=dimnames.error)

C.69.barplot.data.matrix<-matrix(ncol=2,nrow=4,C.69.barplot.data,byrow=FALSE,dimnames=dimnames)
error.bar.data.matrix<-matrix(ncol=3,nrow=2,error.bar.data,byrow=FALSE,dimnames=dimnames.error)

png(file=paste("C.42 barchart showing GNP effect.png"), bg="transparent", width =1000, height=500,units="px",pointsize=13)
C.42.barplot<-barplot(height=C.42.barplot.data.matrix,main="Effect of GNPs on
the ratios of 24mer to monomer in wildtype GLFG",
beside=TRUE,space=c(0,1),col=c("palegreen","darkgreen","salmon","darkred"),
legend.text=TRUE,ylab="Percentage absorbance",ylim=c(0,120))
error.bar(barplot,barplot.data,error.bar.data)
dev.off()

png(file="C.68 barchart showing GNP effect.png", bg="transparent", width =1000, height=500,units="px",pointsize=13)
C.68.barplot<-barplot(height=C.68.barplot.data.matrix,main="Effect of GNPs on
the ratios of 24mer to monomer in F36R GLFG",
beside=TRUE,space=c(0,1),col=c("palegreen","darkgreen","salmon","darkred"),
legend.text=TRUE,ylab="Percentage absorbance",ylim=c(0,120))
error.bar(barplot,barplot.data,error.bar.data)
dev.off()

png(file=paste("C.69 barchart showing GNP effect.png"), bg="transparent", width =1000, height=500,units="px",pointsize=13)
C.69.barplot<-barplot(height=C.69.barplot.data.matrix,main="Effect of GNPs on
the ratios of 24mer to monomer in L162R GLFG",
beside=TRUE,space=c(0,1),col=c("palegreen","darkgreen","salmon","darkred"),
legend.text=TRUE,ylab="Percentage absorbance",ylim=c(0,120))
error.bar(barplot,barplot.data,error.bar.data)
dev.off()

C.68.24.T<-t.test(C.42.peaks.store[,1],C.68.peaks.store[,1])
C.68.mono.T<-t.test(C.42.peaks.store[,2],C.68.peaks.store[,2])
C.69.24.T<-t.test(C.42.peaks.store[,1],C.69.peaks.store[,1])
C.69.mono.T<-t.test(C.42.peaks.store[,2],C.69.peaks.store[,2])
C.68and69.24.T<-t.test(C.68.peaks.store[,1],C.69.peaks.store[,1])
C.68and69.mono.T<-t.test(C.68.peaks.store[,2],C.69.peaks.store[,2])


