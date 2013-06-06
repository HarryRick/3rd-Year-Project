#################### FUNCTION HPLC.import ####################
# Input:
  # 1) Sample ID (e.g. C.29)
	# 2) Data type (discrete or spectrum)
	# 3) Parent directory of datafiles
# Output = read in all data into R
# Development = this will break down if >1 matching files are found... 


HPLC.import <- function (Sample.name,result.type,parent.dir)
{
# Navigate to data
setwd(parent.dir)
data.file<-dir()[grep(x=dir(parent.dir),pattern=result.type,ignore.case=TRUE)]
data.file<-data.file[grep(x=data.file,pattern=Sample.name)]

data.file.names<-unlist(strsplit(x=data.file,split="-"))[seq(1,(3*length(data.file)),3)]

# Read data into R
data<-read.csv(file=data.file,header=FALSE)

}


############################################################



#################### FUNCTION HPLC.discrete.process ####################

# Input = output of HPLC.import()
# Output = list of absorbance values for each of the three wavelengths analysed

HPLC.discrete.process<- function(data,wavelength,time,flow,void.vol)

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
  
  
	# Generate X-axis by calculating elution volume
  
  # Sampling rate = 10x per second
  
  # Flow rate = 0.3 ml / minute
  
  # Experiment time = 45 minutes
  
  # Total Volume = 45 * 0.3
  
  # Total Volume / Total Readings
  
  
  # Elution point (ml) = Time (min) * ( Volume (ml) / Time (min) )
  total.volume<-time * flow
  elution<-  total.volume / (length(wavA)-1)
  elution<- seq (0, (total.volume) , elution)
  
  
  # Order results into a matrix
  col.names<-c("Volume (ml)",paste(Sample,wavelength[1],sep=": "))
  row.names<-seq(1,length(elution),1)
  dimnames<-list(row.names,col.names)
  data<-matrix(data=c(elution,wavA),ncol=2,dimnames=dimnames)

}

############################################################



#################### FUNCTION HPLC.discrete.summary ####################

# Finds and saves elution point of greatest peak

HPLC.discrete.summary<-function(data,void.vol)
{
# Peak
wav1.peak<-as.numeric(data[,1][rev(order(data[,2]))[1]])/void.vol
row.name<-paste(Sample.name, "elution peak (ml)")
dimnames<-list(row.name,wavelength)
sample.summary<-matrix(c(wav1.peak),
ncol=1,nrow=1,dimnames=dimnames)
}
############################################################

HPLC.ladder.plot<-function(albumin,cytochrome.c,hsf)
{
	prot.elutions<-c(albumin,cytochrome.c,hsf)
	molecular.ws<-c(66,44,12.4)
	plot(x=prot.elutions,y=molecular.ws,xlim=c(0,13))
	abline(lm(molecular.ws~prot.elutions))
}


# EXAMPLE

# Preliminary steps
# For each of your samples, the HPLC software creates three data files. Two of these data files will be in .asc format.
# Of the two files in .asc format, the one with the shorter name contains the discrete data.
# Add the word 'Discrete' to the filename of each of the files that contain the discrete dataset.
#  e.g. 'example.asc' would become 'example - Discrete.asc'.

################# DEFINE INPUT VARIABLES

# HPLC.import
Sample.name<-"Albumun 16-08-17" 
Sample<-Sample.name
result.type<-"discrete"
parent.dir<-"/Users/harryrick/Dropbox/Work/Imperial/3rd Year/Nanocage Project/Harry/Data/Ladder Proteins"
void.vol<-0.9

# HPLC.discrete.process
wavelength<-"280 nm"
time<-45
flow<-0.3
figure.dir<-"/Users/harryrick/Dropbox/Work/Imperial/3rd Year/Nanocage Project/Harry/Data/"
colours<-c(3,4,5,6)

################# Calculate Results

data<-HPLC.import(Sample.name,result.type,parent.dir)
data<-HPLC.discrete.process(data,wavelength,time,flow,void.vol)
albumin<-HPLC.discrete.summary(data,void.vol)

Sample.name<-"Cytochrome C"
data<-HPLC.import(Sample.name,result.type,parent.dir)
data<-HPLC.discrete.process(data,wavelength,time,flow,void.vol)
cytochrome.c<-HPLC.discrete.summary(data,void.vol)

Sample.name<-"Horse Spleen Ferritin 15-16-08"
data<-HPLC.import(Sample.name,result.type,parent.dir)
data<-HPLC.discrete.process(data,wavelength,time,flow,void.vol)
hsf<-HPLC.discrete.summary(data,void.vol)

HPLC.ladder.plot(albumin,cytochrome.c,hsf)
