hydrophobes.assign<-function(val,ile,leu,met,phe,trp,cys,ala,one,x)
{
  val<-grep(x=((x$atom[,4][-one])=="VAL"),pattern=TRUE)
  ile<-grep(x=((x$atom[,4][-one])=="ILE"),pattern=TRUE)
  leu<-grep(x=((x$atom[,4][-one])=="LEU"),pattern=TRUE)
  met<-grep(x=((x$atom[,4][-one])=="MET"),pattern=TRUE)
  phe<-grep(x=((x$atom[,4][-one])=="PHE"),pattern=TRUE)
  trp<-grep(x=((x$atom[,4][-one])=="TRP"),pattern=TRUE)
  cys<-grep(x=((x$atom[,4][-one])=="CYS"),pattern=TRUE)
  ala<-grep(x=((x$atom[,4][-one])=="ALA"),pattern=TRUE)

  sort(c(ile,leu,met,phe,trp,cys,val,ala))
}

hydrophobesr.assign<-function(Hcb,Hcg,Hsd,Hce,Hcg1,Hcg2,Hcd1,Hcd2,Hce1,Hce2,Hcz,Hne1,Hce3,Hcz2,Hcz3,Hch2,Hsg,hydrophobes,x)
{
  Hcb<-na.omit(subset(hydrophobes,x$atom[,2]=="CB"))
  Hcg<-na.omit(subset(hydrophobes,x$atom[,2]=="CG"))
  Hsd<-na.omit(subset(hydrophobes,x$atom[,2]=="SD"))
  Hce<-na.omit(subset(hydrophobes,x$atom[,2]=="CE"))
  Hcg1<-na.omit(subset(hydrophobes,x$atom[,2]=="CG1"))
  Hcg2<-na.omit(subset(hydrophobes,x$atom[,2]=="CG2"))
  Hcd1<-na.omit(subset(hydrophobes,x$atom[,2]=="CD1"))
  Hcd2<-na.omit(subset(hydrophobes,x$atom[,2]=="CD2"))
  Hce1<-na.omit(subset(hydrophobes,x$atom[,2]=="CE1"))
  Hce2<-na.omit(subset(hydrophobes,x$atom[,2]=="CE2"))
  Hcz<-na.omit(subset(hydrophobes,x$atom[,2]=="CZ"))
  Hne1<-na.omit(subset(hydrophobes,x$atom[,2]=="NE1"))
  Hce3<-na.omit(subset(hydrophobes,x$atom[,2]=="CE3"))
  Hcz2<-na.omit(subset(hydrophobes,x$atom[,2]=="CZ2"))
  Hcz3<-na.omit(subset(hydrophobes,x$atom[,2]=="CZ3"))
  Hch2<-na.omit(subset(hydrophobes,x$atom[,2]=="CH2"))   
  Hsg<-na.omit(subset(hydrophobes,x$atom[,2]=="SG"))                          
  sort(c(Hcb,Hcg,Hsd,Hce,Hcg1,Hcg2,Hcd1,Hcd2,Hce1,Hce2,Hcz,Hne1,Hce3,Hcz2,Hcz3,Hch2,Hsg))
}

hydrophobes1.assign<-function(val1,ile1,leu1,met1,phe1,trp1,cys1,ala1,one,x)
{
  val1<-grep(x=((x$atom[,4][one])=="VAL"),pattern=TRUE)
  ile1<-grep(x=((x$atom[,4][one])=="ILE"),pattern=TRUE)
  leu1<-grep(x=((x$atom[,4][one])=="LEU"),pattern=TRUE)
  met1<-grep(x=((x$atom[,4][one])=="MET"),pattern=TRUE)
  phe1<-grep(x=((x$atom[,4][one])=="PHE"),pattern=TRUE)
  trp1<-grep(x=((x$atom[,4][one])=="TRP"),pattern=TRUE)
  cys1<-grep(x=((x$atom[,4][one])=="CYS"),pattern=TRUE)
  ala1<-grep(x=((x$atom[,4][one])=="ala"),pattern=TRUE)
  #All hydrophobic amino acids from the first chain are sorted
  sort(c(ile1,leu1,met1,phe1,trp1,cys1,val1,ala1))
}

hydrophobesr1.assign<-function(Hcb.1,Hcg.1,Hsd.1,Hce.1,Hcg1.1,Hcg2.1,Hcd1.1,Hcd2.1,Hce1.1,Hce2.1,
Hcz.1,Hne1.1,Hce3.1,Hcz2.1,Hcz3.1,Hch2.1,Hsg.1,hydrophobes1,x)
{
  Hcb.1<-na.omit(subset(hydrophobes1,x$atom[,2]=="CB"))
  Hcg.1<-na.omit(subset(hydrophobes1,x$atom[,2]=="CG"))
  Hsd.1<-na.omit(subset(hydrophobes1,x$atom[,2]=="SD"))
  Hce.1<-na.omit(subset(hydrophobes1,x$atom[,2]=="CE"))
  Hcg1.1<-na.omit(subset(hydrophobes1,x$atom[,2]=="CG1"))
  Hcg2.1<-na.omit(subset(hydrophobes1,x$atom[,2]=="CG2"))
  Hcd1.1<-na.omit(subset(hydrophobes1,x$atom[,2]=="CD1"))
  Hcd2.1<-na.omit(subset(hydrophobes1,x$atom[,2]=="CD2"))
  Hce1.1<-na.omit(subset(hydrophobes1,x$atom[,2]=="CE1"))
  Hce2.1<-na.omit(subset(hydrophobes1,x$atom[,2]=="CE2"))
  Hcz.1<-na.omit(subset(hydrophobes1,x$atom[,2]=="CZ"))
  Hne1.1<-na.omit(subset(hydrophobes1,x$atom[,2]=="NE1"))
  Hce3.1<-na.omit(subset(hydrophobes1,x$atom[,2]=="CE3"))
  Hcz2.1<-na.omit(subset(hydrophobes1,x$atom[,2]=="CZ2"))
  Hcz3.1<-na.omit(subset(hydrophobes1,x$atom[,2]=="CZ3"))
  Hch2.1<-na.omit(subset(hydrophobes1,x$atom[,2]=="CH2"))                    
  Hsg.1<-na.omit(subset(hydrophobes1,x$atom[,2]=="SG"))                                  
  sort(c(Hcb.1,Hcg.1,Hsd.1,Hce.1,Hcg1.1,Hcg2.1,Hcd1.1,Hcd2.1,Hce1.1,Hce2.1,Hcz.1,Hne1.1,Hce3.1,Hcz2.1,Hcz3.1,Hch2.1))
}



# New function to translate character vector of one letter amino acid sequence to three letter code (doesn't work yet)
one.to.three.translate<-function(x)
{
  i<-1
  while(i<=length(x))
	{	
		if(nchar(x[i])==1)
	 		{x[i]<-sub("A","ALA",x[i])}
	 	if(nchar(x[i])==1)
			{x[i]<-sub("R","ARG",x[i])}
  		if(nchar(x[i])==1)
			{x[i]<-sub("N","ASN",x[i])}
		if(nchar(x[i])==1)
  			{x[i]<-sub("D","ASP",x[i])}
		if(nchar(x[i])==1)
			{x[i]<-sub("C","CYS",x[i])}
 		if(nchar(x[i])==1)
	 		{x[i]<-sub("Q","GLN",x[i])}
  		if(nchar(x[i])==1)
  			{x[i]<-sub("E","GLU",x[i])}
  		if(nchar(x[i])==1)
			{x[i]<-sub("G","GLY",x[i])}
		if(nchar(x[i])==1)
			{x[i]<-sub("H","HIS",x[i])}
		if(nchar(x[i])==1)
			{x[i]<-sub("I","ILE",x[i])}
  		if(nchar(x[i])==1)
  			{x[i]<-sub("L","LEU",x[i])}
		if(nchar(x[i])==1)
			{x[i]<-sub("K","LYS",x[i])}
		if(nchar(x[i])==1)
			{x[i]<-sub("M","MET",x[i])}
  		if(nchar(x[i])==1)
  			{x[i]<-sub("F","PHE",x[i])}
  		if(nchar(x[i])==1)
  			{x[i]<-sub("P","PRO",x[i])}
  		if(nchar(x[i])==1)
  			{x[i]<-sub("S","SER",x[i])}
  		if(nchar(x[i])==1)
  			{x[i]<-sub("T","THR",x[i])}
  		if(nchar(x[i])==1)
  			{x[i]<-sub("W","TRP",x[i])}
  		if(nchar(x[i])==1)
		  	{x[i]<-sub("Y","TYR",x[i])}
		if(nchar(x[i])==1)
	 	 	{x[i]<-sub("V","VAL",x[i])}
		i<-i+1	
	}
return(x)
}

# Generates mutation codes for primerX

hydrophobic.3to1.translate<-function(x)	
{
	x<-sub("CYS-","C",x,fixed=TRUE)
	x<-sub("ALA-","A",x,fixed=TRUE)
	x<-sub("ILE-","I",x,fixed=TRUE)
	x<-sub("LEU-","L",x,fixed=TRUE)
	x<-sub("PHE-","F",x,fixed=TRUE)
	x<-sub("MET-","M",x,fixed=TRUE)
	x<-sub("TRP-","W",x,fixed=TRUE)
	x<-sub("VAL-","V",x,fixed=TRUE)

	return(x)
}


seq.split<-function(x)
	
{
	seq.as.character<-character()
	i<-1
	while(i<=nchar(x))
	{
		seq.as.character[i]<-substring(seq,i,i)	
		i<-i+1
	}

return(seq.as.character)
}
