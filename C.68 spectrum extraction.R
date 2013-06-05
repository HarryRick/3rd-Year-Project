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

HPLC.spectrum.process<- function(data,wavelength,time,flow)

{

# Determine number of datapoints for each discrete wavelength
  points<-as.numeric(array(unlist(data[15,2])))-1
  largemer.time<-as.numeric(sub(" Min","",array(unlist(data[8055,2]))))
  smallmer.time<-as.numeric(sub(" Min","",array(unlist(data[9001,2]))))

  largemer.end<-8057+points
  smallmer.end<-9003+points 
  
  # Gets absorbance data for 24mer and monomer and wavelength
  wavelength<-as.numeric(array(unlist(data[8057:largemer.end,1])))
  largemer.absorbance<-as.numeric(array(unlist(data[8057:largemer.end,2])))
  smallmer.absorbance<-as.numeric(array(unlist(data[9003:smallmer.end,2])))
  
  max.absorbance<-max(c(largemer.absorbance,smallmer.absorbance))    
  
  adj.largemer.absorbance<-largemer.absorbance/max(largemer.absorbance)
  adj.smallmer.absorbance<-smallmer.absorbance/max(smallmer.absorbance)
  
  # Generate X-axis by calculating elution volume
  
  # Sampling rate = 10x per second
  
  # Flow rate = 0.3 ml / minute
  
  # Experiment time = 45 minutes
  
  # Total Volume = 45 * 0.3
  
  # Total Volume / Total Readings
  
  
  # Elution point (ml) = Time (min) * ( Volume (ml) / Time (min) )
  largemer.volume<-largemer.time * flow
  smallmer.volume<-smallmer.time * flow  
  
  # Order results into a matrix
  col.names<-c("Wavelength",paste("24-mer, vol=",largemer.volume,sep=""),paste("Monomer, vol=",smallmer.volume,sep=""))
  row.names<-seq(1,length(wavelength),1)
  dimnames<-list(row.names,col.names)
  data<-matrix(data=c(wavelength,adj.largemer.absorbance,adj.smallmer.absorbance),ncol=3,dimnames=dimnames)
  
}

############################################################


#################### FUNCTION HPLC.spectrum.plot1 ####################

# Input = output of HPLC.import()
# Output = list of absorbance values for each of the three wavelengths analysed

HPLC.spectrum.plot1<- function(hplc.data,figure.dir,colours)
{
  setwd(figure.dir)
  png(file="F36R GLFG Absorption Spectra.png", bg="transparent", width =1000, height=500,units="px",pointsize=13)
  plot(x=data[,1],y=data[,2],type='l',col=colours[1],xlab="Wavelength",ylab="Relative Absorbance (A.U)",main="F36R GLFG Absorption Spectra", 
  lwd=2,ylim=c(-0.1,1),cex.main=3,cex.lab=1.2,cex.axis=1.2)
  points(x=data[,1],y=data[,3],type='l',col=colours[2],lwd=2)
  leg.lab<-c(paste(Sample,"24-mer"),paste(Sample,"Monomer"))
  legend('topright',legend=leg.lab,col=colours,lty=1,lwd=3,cex=1.1, bty='n')
  dev.off()
}


############################################################


# EXAMPLE

# Preliminary steps
# For each of your samples, the HPLC software creates three data files. Two of these data files will be in .asc format.
# Of the two files in .asc format, the one with the shorter name contains the discrete data.
# Add the word 'Discrete' to the filename of each of the files that contain the discrete dataset.
#  e.g. 'example.asc' would become 'example - Discrete.asc'.

################# DEFINE INPUT VARIABLES

# HPLC.import
Sample.name<-"C.68 Run 2" 
Sample<-Sample.name
result.type<-"spectrum"
parent.dir<-"//ic.ac.uk/homes/hfr10/2013-05-23"


flow<-0.3

################# Calculate Results
data<-HPLC.import(Sample,result.type,parent.dir)
data<-HPLC.spectrum.process(data,wavelength,time,flow)

# HPLC.discrete.plot1 
figure.dir<-"//ic.ac.uk/homes/hfr10/2013-05-23"
colours<-c(3,4)

################# Calculate Results

HPLC.spectrum.plot1(data,figure.dir,colours)
#rm(list=ls())

