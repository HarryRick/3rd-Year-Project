# Load R extensions
require(BoSSA)
require(rgl)
require(bio3d)
require(seqinr)
require(RCurl)
require(XML)

# Add BLAST search of pdb database and use top result from this as input to rest of script 

seq<-as.character("MSSQIRQNYSTDVEAAVNSLVNLYLQASYTYLSLGFYFDRDDVALEGVSHFFRELAEEKREGYERLLKMQNQRGGRALFQ
DIKKPAEDEWGKTPDAMKAAMALEKKLNQALLDLHALGSARTDPHLCDFLETHFLDEEVKLIKKMGDHLTNLHRLGGPEA
GLGEYLFERLTLKHD")

a<-postForm("http://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&BLAST_PROGRAMS=blastp&DATABASE=pdb&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on&BLAST_SPEC=&LINK_LOC=blasttab&LAST_PAGE=blastp&CMD=PUT",QUERY=seq,.cgifields =c("BLAST"))

# Extract the RID
RIDa<-strsplit(x=a,split="RID")[[1]][6]
RIDb<-strsplit(x=RIDa,split="value=\"")[[1]][2]
RIDc<-strsplit(x=RIDb,split="\"")[[1]]
RIDd<-RIDc[1]

# Define the URL where the results can be found
Search.output<-"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?RESULTS_FILE=on&RID=XXX&FORMAT_TYPE=XML&FORMAT_OBJECT=Alignment&ALIGNMENTS=100&CMD=Get"

# Enter the RID into the URL
NewDestination<-gsub("XXX",RIDd,Search.output)

# Read the results into R
fetched.data<-readLines(NewDestination)
not.finished<-grep(x=fetched.data,pattern="waiting",ignore.case=TRUE)

while(length(not.finished)>0)
{
	fetched.data<-readLines(NewDestination)
	not.finished<-grep(x=fetched.data,pattern="waiting",ignore.case=TRUE)
}

xml.data<-xmlTreeParse(fetched.data)

#Parsing of data (XML)


blast.hits<-xmlRoot(xml.data)[["BlastOutput_iterations"]][["Iteration"]][["Iteration_hits"]]

ids<-character()
i<-1
while(i<=xmlSize(blast.hits))
{
	ids[[i]] <- xmlValue(blast.hits[[i]][["Hit_id"]])
	i<-i+1 
}

ids1<-strsplit(ids,"|",fixed=TRUE)
pdb.ids<-grep(ids1,pattern="gi|",fixed=TRUE,value=TRUE,invert=TRUE)



#Parsing of data (TEXT)
sig.align.start<-(grep(result,pattern="Sequences producing significant alignments"))+2
sig.align.end<-(grep(result,pattern="ALIGNMENTS"))-2

all.pdb.result<-result[sig.align.start:sig.align.end]

all.pdb.result<-strsplit(all.pdb.result,"  Chain")
all.pdb.result1<-character()
i<-1
while(i<=length(all.pdb.result))
{
	all.pdb.result1<-c(all.pdb.result1,strsplit(all.pdb.result[[i]],"    ",fixed=TRUE))

	i<-i+1
}


blast.pdbs<-character()

blast.pdbs<-grep(all.pdb.result1,pattern="pdb|",fixed=TRUE,value=TRUE)
blast.pdbs<-grep(pattern=">pdb|",x=blast.pdbs,fixed=TRUE,value=TRUE,invert=TRUE)

i<-1

while(i<=length(blast.pdbs))
{
	as.list(all.pdb.result1<-grep(all.pdb.result1,pattern=blast.pdbs[i],fixed=TRUE,value=TRUE,invert=TRUE))
	i<-i+1
}

all.pdb.result1<-gsub('\"',"",all.pdb.result1)
all.pdb.result1<-gsub(', ',"  ",all.pdb.result1)
all.pdb.result1<-strsplit(all.pdb.result1,"  ")
i<-1
ch.blast.result<-character()
while(i<=length(all.pdb.result1))
{
	ch.blast.result<-c(ch.blast.result,as.character(all.pdb.result1[[i]]))
	i<-i+1
}
	
all.pdb.result1<-grep(all.pdb.result1,pattern="c( ",fixed=TRUE,value=TRUE,invert=TRUE)

# General import pdb - User enters pdb id (or is obtained from blast) - script finds relevant url. 
pdbid<-"3ajo"
pdb.url<-sub("___",pdbid,"http://www.rcsb.org/pdb/files/___.pdb1",fixed=TRUE)

x<-read.pdb(pdb.url,multi=TRUE)

###### Check that structure has correct x$atom[,5] formatting and if not correct ###### 
# Input chain number
chain.num<-24

# Where required find and store correct chain names 
if (length(unique(x$atom[,5])) < chain.num ) 
  
  { 
    new.chain.store<-numeric(0)
    chain.lib<-toupper(paste(letters[1:chain.num]))

# Calculate atoms per chain

    atoms.per.chain<-length(x$atom[,5])/chain.num
    i<-1

# Replace old chian names with new ones
    while (i<=chain.num)

    {
      new.chain<-rep(chain.lib[i],atoms.per.chain)
      new.chain.store<-c(new.chain.store,new.chain)

      i<-i+1
		}
	}
x$atom[,5]= new.chain.store

# Output of changechange without duplicates
chains<-unique(new.chain.store)

# Assigns each amino acid its specific chain.
reference.ids<-numeric(0)
i<-1
  while (i<=length(x$atom[,6]))
  {
    reference.ids<-c(reference.ids,paste(c(x$atom[,5][i],"-",x$atom[,6][i]),collapse=""))
    i<-i+1
  }


# one contains only the data from the first chain and none of the others.
one<-chains[1]
one<-grep(x=((x$atom[,5])==one),pattern=TRUE)


# Selects all atoms from hydorphobic amino acids which are not in the first chain
hydrophobes<-hydrophobes.assign(val,ile,leu,met,phe,trp,cys,ala,one,x)

#Selects only atoms of the R groups of hydrophobic amino acids
hydrophobesr<-hydrophobesr.assign(Hcb,Hcg,Hsd,Hce,Hcg1,Hcg2,Hcd1,Hcd2,Hce1,Hce2,Hcz,Hne1,
Hce3,Hcz2,Hcz3,Hch2,Hsg,hydrophobes,x)

# Groups and sorts all hydrophobic amino acids from the first chain
hydrophobes1<-hydrophobes1.assign(val1,ile1,leu1,met1,phe1,trp1,cys1,ala1,one,x)

# Groups and sorts all R group atoms from hydrophobic amino acids of the first chain
hydrophobesr1<-hydrophobesr1.assign(Hcb.1,Hcg.1,Hsd.1,Hce.1,Hcg1.1,Hcg2.1,Hcd1.1,Hcd2.1,Hce1.1,Hce2.1,
Hcz.1,Hne1.1,Hce3.1,Hcz2.1,Hcz3.1,Hch2.1,Hsg.1,hydrophobes1,x)

#Assigns all xyz data to relevant atoms
# Contains positional data for all atoms
x.val<-as.numeric(x$atom[,8])
y.val<-as.numeric(x$atom[,9])
z.val<-as.numeric(x$atom[,10])

# Contains positional data for all hydrophobes 
x.hydro<-as.numeric(x$atom[,8][-one][hydrophobesr])
y.hydro<-as.numeric(x$atom[,9][-one][hydrophobesr])
z.hydro<-as.numeric(x$atom[,10][-one][hydrophobesr])

# Contains positional data for hydrophobias from first chain
One.xhydro<-as.numeric(x$atom[,8][hydrophobesr1])
One.yhydro<-as.numeric(x$atom[,9][hydrophobesr1])
One.zhydro<-as.numeric(x$atom[,10][hydrophobesr1])

# Plots the first chain 
plot3d(x=One.xhydro,y=One.yhydro,z=One.zhydro,box=0,axes=0,col=7,type="s",radius=0.5)

# Run a loop to find the proximity of other atoms in structure 
# Future proximity residue will be stored in these
residue.match.store<-numeric(0)
residue.match.store2<-numeric(0)

j<-1

while (j<=length(One.xhydro))
  
{
  x1<-One.xhydro[j]
  y1<-One.yhydro[j]
  z1<-One.zhydro[j]
  
  
  proximity.store<-numeric(0)
  
  # Calculates distance between residues - if these residues are less than 5 A apart then their coordinates are stored in residue.match.store. 
  # While i is less than or equal to the total number of positional data points all but the first chain, the distance between the residues is calculated
  i<-1
  
  while (i<=length(x.hydro))
  {
    #Calculates distance between residues
    proximity<-(((x1-x.hydro[i])^2)+((y1-y.hydro[i])^2)+((z1-z.hydro[i])^2))^0.5
    #If distance < 5
    if (proximity < 5)
    {
      
      
      # residue contains each atom closer than 4 angstroms apart and what chain they are on.
      
      
      # Which residue has a hydrophobic interaction
      residue<-paste(c(x$atom[,5][one][hydrophobesr1][j],"-",x$atom[,6][one][hydrophobesr1][j]),collapse="")
      residue.match<-grep(pattern=TRUE,x=residue==reference.ids)
      
      # residue.match.store combines the data from residue and residue.match and eliminates any duplicate values.
      
      residue.match.store<-c(residue.match.store,residue.match)
      
      residue.match.store<-unique(residue.match.store)
      
      # this gives us all the spatial data points for each atom within 5 angstroms of one another so that they may be plotted in space
      
      interface.x<-x$atom[,8][residue.match.store]
      interface.y<-x$atom[,9][residue.match.store]
      interface.z<-x$atom[,10][residue.match.store]
      plot3d(x=interface.x,y=interface.y,z=interface.z,col=3,box=FALSE)
      
      # Find partner
      residue2<-paste(c(x$atom[,5][-one][hydrophobesr][i],"-",x$atom[,6][-one][hydrophobesr][i]),collapse="")
      residue.match2<-grep(pattern=TRUE,x=residue2==reference.ids)
      
      # residue.match.store combines the data from residue and residue.match and eliminates any duplicate values.
      
      residue.match.store2<-c(residue.match.store2,residue.match2)
      
      residue.match.store2<-unique(residue.match.store2)
      
      # this gives us all the spatial data points for each atom within 5 angstroms of one another so that they may be plotted in space
      
      interface.x2<-x$atom[,8][residue.match.store2]
      interface.y2<-x$atom[,9][residue.match.store2]
      interface.z2<-x$atom[,10][residue.match.store2]
      plot3d(x=interface.x2,y=interface.y2,z=interface.z2,col=2,box=FALSE,add=TRUE)
      print(proximity)
      
      
    }
    i<-i+1
  }
  
  j<-j+1
}




interface.x<-as.numeric(interface.x)
interface.y<-as.numeric(interface.y)
interface.z<-as.numeric(interface.z)



chain.match<-x$atom[,5][residue.match.store]

#total spatial points for first chain
onexval<-x.val[one]
oneyval<-y.val[one]
onezval<-z.val[one]

#Adds all atoms to the 3d plot - can see where hydrophobic amino acids are in relation to non hydrophobic ones.
points3d(x=onexval,y=oneyval,z=onezval,col="blue",box=FALSE,axes=FALSE)



#Gives list of hydrophobic amino acids on first chain < 5 angstroms from another hydrophobic amino acid on a different chain

atom.aminoacid.match<-numeric(0)
residue.no.store<-numeric(0)
p<-1




while(p<=length(residue.match.store))
  
{
  residue.no.store<-c(residue.no.store,x$atom[residue.match.store,6])
  atom.aminoacid.match<-c(atom.aminoacid.match,paste(c(x$atom[residue.match.store,4][p],"-",x$atom[residue.match.store,6][p]),collapse=""))
  
  p<-p+1
  
}

residue.no.store<-unique(residue.no.store)
residue.no.store<-toString(residue.no.store)
aminoacid.match<-unique(atom.aminoacid.match)
print(aminoacid.match)


atom.aminoacid.match2<-numeric(0)
residue.no.store2<-numeric(0)
p2<-1

while(p2<=length(residue.match.store2))
  
{
  residue.no.store2<-c(residue.no.store2,x$atom[residue.match.store2,6])
  atom.aminoacid.match2<-c(atom.aminoacid.match2,paste(c(x$atom[residue.match.store2,4][p2],"-",x$atom[residue.match.store2,6][p2]),collapse=""))
  
  p2<-p2+1
  
}

residue.no.store2<-unique(residue.no.store2)
toString(residue.no.store2)
aminoacid.match2<-unique(atom.aminoacid.match2)
print(aminoacid.match2)


##### Create .cmd script to open pdb file of protein in Chimera and highlight interacting hydrophobic residues #####

open<-sub("___",pdbid,"open ___")

add.sym<-"sym #0"

# Makes residue string compatible with chimera command line by deleting spaces between residue numbers

residue.no.store<-gsub(", ",",",residue.no.store,fixed=TRUE)

colour<-sub("___",residue.no.store,"colour red :___")

show<-sub("___",residue.no.store,"show :___")

compile<-c(open,add.sym,colour,show)

file.name<-sub("___",pdbid,"___.cmd")
write(x=compile,file=file.name)

# Draw interface mesh
j<-1
while (j<length(chain.match)+1)
  
{
  x1<-interface.x[j]
  y1<-interface.y[j]
  z1<-interface.z[j]
  
  proximity.store<-numeric(0)
  
  
  i<-1
  
  while (i<length(interface.x)+1)
  {
    #Calculate distance
    proximity<-(((x1-interface.x[i])^2)+((y1-interface.y[i])^2)+((z1-interface.z[i])^2))^0.5
    #Save distance
    proximity.store<-c(proximity.store,proximity)
    i<-i+1
  }
  
  sequ<-(rev(seq(0,1,(1/6)))^2)[1:6]
  
  t<-1
  while(t<6)
  {
    x3<-c(interface.x[j],interface.x[order(proximity.store)[-1][t]])
    y3<-c(interface.y[j],interface.y[order(proximity.store)[-1][t]])
    z3<-c(interface.z[j],interface.z[order(proximity.store)[-1][t]])
    plot3d(x=x3,y=y3,z=z3,add=TRUE,type="l",col=grep(pattern=chain.match[j],x=unique(chain.match)),lwd=5,alpha=sequ[t])
    t<-t+1
  }
  j<-j+1
}
