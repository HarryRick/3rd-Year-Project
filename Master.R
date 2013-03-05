# Load R extensions
require(rgl)
require(bio3d)
require(seqinr)

x<-read.pdb("http://www.rcsb.org/pdb/files/2fg8.pdb")


###### Check that structure has correct x$atom[,5] formatting and if not correct ###### 
# Input chain number
chain.num<-24

# Where required correct chain names 

if (length(unique(x$atom[,5])) < chain.num ) 
  
  { 
    new.chain.store<-numeric(0)
    chain.lib<-toupper(paste(letters[1:chain.num]))

# Calculate atoms per chain

atoms.per.chain<-length(x$atom[,5])/chain.num
i<-1
while (i<=chain.num)

{
new.chain<-rep(chain.lib[i],atoms.per.chain)
new.chain.store<-c(new.chain.store,new.chain)

i<-i+1

}
x$atom[,5]= new.chain.store
}









# Assigns each amino acid its specific chain.
reference.ids<-numeric(0)
i<-1
while (i<length(x$atom[,6])+1)
{
  reference.ids<-c(reference.ids,paste(c(x$atom[,5][i],"-",x$atom[,6][i]),collapse=""))
  i<-i+1
}

#List the unique chains in the structure
chains<-unique(x$atom[,5])


# one contains only the data from the first chain and none of the others.
one<-chains[1]
one<-grep(x=((x$atom[,5])==one),pattern=TRUE)


# Each contains all the atoms from hydrophobic amino acids (not in 1st chain)

val<-grep(x=((x$atom[,4][-one])=="VAL"),pattern=TRUE)
ile<-grep(x=((x$atom[,4][-one])=="ILE"),pattern=TRUE)
leu<-grep(x=((x$atom[,4][-one])=="LEU"),pattern=TRUE)
met<-grep(x=((x$atom[,4][-one])=="MET"),pattern=TRUE)
phe<-grep(x=((x$atom[,4][-one])=="PHE"),pattern=TRUE)
trp<-grep(x=((x$atom[,4][-one])=="TRP"),pattern=TRUE)
cys<-grep(x=((x$atom[,4][-one])=="CYS"),pattern=TRUE)
ala<-grep(x=((x$atom[,4][-one])=="ALA"),pattern=TRUE)


#All hydrophobic amino acids are grouped together

hydrophobes<-sort(c(ile,leu,met,phe,trp,cys,val,ala))

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
          
           
hydrophobesr<-sort(c(Hcb,Hcg,Hsd,Hce,Hcg1,Hcg2,Hcd1,Hcd2,Hce1,Hce2,Hcz,Hne1,Hce3,Hcz2,Hcz3,Hch2,Hsg))



#Each contains hydrophobic amino acids from the first chain
val1<-grep(x=((x$atom[,4][one])=="VAL"),pattern=TRUE)
ile1<-grep(x=((x$atom[,4][one])=="ILE"),pattern=TRUE)
leu1<-grep(x=((x$atom[,4][one])=="LEU"),pattern=TRUE)
met1<-grep(x=((x$atom[,4][one])=="MET"),pattern=TRUE)
phe1<-grep(x=((x$atom[,4][one])=="PHE"),pattern=TRUE)
trp1<-grep(x=((x$atom[,4][one])=="TRP"),pattern=TRUE)
cys1<-grep(x=((x$atom[,4][one])=="CYS"),pattern=TRUE)
ala1<-grep(x=((x$atom[,4][one])=="ala"),pattern=TRUE)


#All hydrophobic amino acids from the first chain are grouped together
hydrophobes1<-sort(c(ile1,leu1,met1,phe1,trp1,cys1,val1,ala1))

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
                
hydrophobesr1<-sort(c(Hcb.1,Hcg.1,Hsd.1,Hce.1,Hcg1.1,Hcg2.1,Hcd1.1,Hcd2.1,Hce1.1,Hce2.1,Hcz.1,Hne1.1,Hce3.1,Hcz2.1,Hcz3.1,Hch2.1))

#Assigning positional data of first chain to x,y and z axes so that atoms can be plotted in space 

# Contains positional data for hydrophobias from first chain
One.xhydro<-as.numeric(x$atom[,8][hydrophobesr1])
One.yhydro<-as.numeric(x$atom[,9][hydrophobesr1])
One.zhydro<-as.numeric(x$atom[,10][hydrophobesr1])

# Plots the first chain 
plot3d(x=One.xhydro,y=One.yhydro,z=One.zhydro,box=0,axes=0,col=7,type="s",radius=0.5)


x.val<-as.numeric(x$atom[,8])
y.val<-as.numeric(x$atom[,9])
z.val<-as.numeric(x$atom[,10])


# These contain positional data for all hydrophobes 

x.hydro<-as.numeric(x$atom[,8][-one][hydrophobesr])
y.hydro<-as.numeric(x$atom[,9][-one][hydrophobesr])
z.hydro<-as.numeric(x$atom[,10][-one][hydrophobesr])

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
  
  #i is the total number of atoms in the first chain
  
  # While i is less than or equal to the total number of positional data points in the first chain, the distance between the residues is calculated
  
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

p<-1

while(p<length(residue.match.store)+1)
  
{
  
  atom.aminoacid.match<-c(atom.aminoacid.match,paste(c(x$atom[residue.match.store,4][p],"-",x$atom[residue.match.store,6][p]),collapse=""))
  
  p<-p+1
  
}


aminoacid.match<-unique(atom.aminoacid.match)
print(aminoacid.match)

atom.aminoacid.match2<-numeric(0)

p2<-1

while(p2<length(residue.match.store2)+1)
  
{
  
  atom.aminoacid.match2<-c(atom.aminoacid.match2,paste(c(x$atom[residue.match.store2,4][p2],"-",x$atom[residue.match.store2,6][p2]),collapse=""))
  
  p2<-p2+1
  
}


aminoacid.match2<-unique(atom.aminoacid.match2)
print(aminoacid.match2)

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
