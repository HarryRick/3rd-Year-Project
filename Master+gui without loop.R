# Load R extensions
require(BoSSA)
require(rgl)
require(bio3d)
require(RCurl)
require(XML)
require(Biostrings)
require(gWidgets)
options("guiToolkit"="RGtk2")

source(file="Functions.r")

sig.pdb.ids<-""

dialog<-gbasicdialog("Select DNA Sequence input",do.buttons=FALSE)
win1<-gframe(Horizontal=TRUE,container=dialog,spacing=20)
space<-glabel("",container=win1,spacing=0)
textframe<-ggroup(horizontal=TRUE,container=win1,spacing=0,pos=0)
button<-ggroup(horizontal=TRUE,container=win1,spacing=25,pos=1)
win2<-gframe(horizontal=TRUE,container=dialog,spacing=20)
text<-gtext("Enter the DNA sequence of your protein here",
container=textframe,wrap=TRUE,width=425,height=50)
but<-gbutton("Proceed",container=button,handler=function(h,...)
{	
	if(svalue(text)=="")
	{
		gmessage("Please enter DNA sequence!",title="Error!")
	}
	else
	{	svalue(text)<<-toupper(svalue(text))
		if(grepl(" ",svalue(text))==FALSE,fixed=TRUE)
	{
		gmessage("Please ensure there are no spaces in your DNA sequence",title="Error!")
	}
	else if(grepl("^[ATCG]+$",svalue(text))==FALSE)
	{
		gmessage("Please enter a valid DNA sequence!",title="Error!")
	}
	else 
	{
		dna.seq<<-svalue(text)
		dispose(dialog)
	}
})
lab2<-glabel("Or select a plain text file containing the DNA sequence for your protein",
container=win2)
but2<-gbutton("Select file…",container=win2,handler=function(h,...)
{
	dna.seq<<-scan(gfile(),what=character())
	dispose(dialog)
}
)
visible(dialog,TRUE)

species.list<<-c("Aplysia californica","Arabidopsis thaliana","Bos taurus","Caenorhabditis elegans",
"Candida albicans","Canis familiaris","Chlamydomonas reinhardtii","Danio rerio","Dictyostelium
discoideum","Drosophila melanogaster","Equus caballus","Escherichia coli","Felis catus",
"Gallus gallus","Homo sapiens","Mus musculus","Mycoplasma pneumoniae","Oryza sativa",
"Oryzias latipes","Ovis aries","Plasmodium falciparum","Pneumocystis carinii","Rattus norvegicus",
"Saccharomyces cerevisiae","Schizosaccharomyces pombe","Sus scrofa domestica","Takifugu rubripes",
"Xenopus laevis","Zea mays")
primer.type<<-c("Complementary","Overlapping")

default.primer.question<-gbasicdialog("Select primer conditions",do.buttons=FALSE,guiToolkit="RGtk2")
but.frame<-gframe(horizontal=FALSE,container=default.primer.question,spacing=20)
label<-glabel("Please select whether you would like to use default settings for 
primer creation or whether you would like to specify your own",container=but.frame)
but1<-gbutton("Default",container=but.frame,handle=function(h,...)
{
	species<<-"Escherichia coli"
	min_Tm<<-"70"
	max_Tm<<-"90"
	min_GC<<-"40"
	max_GC<<-"60"
	min_length<<-"25"
	max_length<<-"45"
	min_5p_flank<<-"11"
	max_5p_flank<<-"21"
	min_3p_flank<<-"11"
	max_3p_flank<<-"21"
	ends_in_GC<<-"on"
	mut_at_center<<-"on"
	primer_type<<-"complementary"
	dispose(default.primer.question)
})
but2<-gbutton("Custom",container=but.frame,handle=function(h,...)
{	
	dispose(default.primer.question)
	custom<-gbasicdialog("Select primer conditions",do.buttons=FALSE)
	main<-ggroup(horizontal=FALSE,container=custom,spacing=15)
	species.frame<-gframe(horizontal=TRUE,container=main,spacing=10)
	lab1<-glabel("Expression system:",container=species.frame)
	species.drop<-gdroplist(species.list,selected=13,container=species.frame,handler=function(h,...)
	{
		species<<-svalue(species.drop)
	})
	
	Tm.frame<-gframe(horizontal=TRUE,container=main,spacing=10)
	lab2<-glabel("Minimum melting temperature (ºC):",container=Tm.frame)
	min.temp<-gtext("70",width=50,height=15,container=Tm.frame,handler=function(h,...)
	{
		min_Tm<<-svalue(min.temp)
	})
	max.temp<-gtext("90",width=50,height=15,container=Tm.frame,handler=function(h,...)
	{
		max_Tm<<-svalue(max.temp)
	})
	lab3<-glabel("Maximum melting temperature (ºC)",container=Tm.frame)
	
	gc.frame<-gframe(container=main,spacing=10)
	lab4<-glabel("              Minimum GC content (%):",container=gc.frame)
	min.gc<-gtext("40",width=50,height=15,container=gc.frame,handler=function(h,...)
	{
		min_GC<<-svalue(min.gc)
	})
	max.gc<-gtext("60",width=50,height=15,container=gc.frame,handler=function(h,...)
	{
		max_GC<<-svalue(max.gc)
	})	
	lab5<-glabel("Maximum GC content (%)",container=gc.frame)
	
	length.frame<-gframe(container=main,spacing=10)
	lab4<-glabel("                      Minimum length (bp):",container=length.frame)
	min.length<-"25",gtext(width=50,height=15,container=length.frame,handler=function(h,...)
	{
		min_length<<-svalue(min.length)
	})
	max.length<-gtext("45",width=50,height=15,container=length.frame,handler=function(h,...)
	{
		max_length<<-svalue(max.length)
	})
	lab5<-glabel("Maximum length (bp)",container=length.frame)
	
	five.flank.frame<-gframe(container=main,spacing=10)
	lab4<-glabel("     Minimum 5\'\ flanking region (bp):",container=five.flank.frame)
	min.five.flank<-gtext("11",width=50,height=15,container=five.flank.frame,handler=function(h,...)
	{
		min_5p_flank<<-svalue(min.five.flank)
	})
	max.five.flank<-gtext("21",width=50,height=15,container=five.flank.frame,handler=function(h,...)
	{
		max_5p_flank<<-svalue(max.five.flank)
	})
	lab5<-glabel("Maximum 5' flanking region (bp)",container=five.flank.frame)

	three.flank.frame<-gframe(container=main,spacing=10)
	lab5<-glabel("     Minimum 3\'\ flanking region (bp):",container=three.flank.frame)
	min.three.flank<-gtext("11",width=50,height=15,container=three.flank.frame,handler=function(h,...)
	{
		min_3p_flank<<-svalue(min.three.flank)
	})
	max.three.flank<-gtext("21",width=50,height=15,container=three.flank.frame,handler=function(h,...)
	{
		max_3p_flank<<-svalue(max.three.flank)
	})
	lab5<-glabel("Maximum 3' flanking region (bp)",container=three.flank.frame)

	GC.terminate<-gcheckbox("Primers terminate in G or C",checked=TRUE,container=main,
	handler=function(h,...)
	{
		if(svalue(GC.terminate)==TRUE)
		{
			ends_in_GC<<-"on"
		}
		else 
		{
			ends_in_GC<<-""
		}
	})
	
	mut.at.center<-gcheckbox("Mutation at centre of primer",checked=TRUE,
	container=main,
	handler=function(h,...)
	{
		if(svalue(mut.at.center)==TRUE)
		{
			mut_at_center<<-"on"
		}
		else 
		{
			mut_at_center<<-""
		}
	})
	type.and.button<-ggroup(container=main,spacing=10)
	primer.type.select<-gradio(primer.type,container=type.and.button,
	handler=function(h,...)
	{
		if(svalue(primer.type.select)=="Complementary")
		{	
			primer_type<<-"complementary"
		}
		else {primer_type<<-"overlapping"}
	})
	addSpace(type.and.button,255)
	done.button<-gbutton("ok",container=type.and.button,handler=function(h,...) 
	{
		if(svalue(min.temp)==""|svalue(max.temp)==""|svalue(min.gc)==""|svalue(max.gc)==""|
		svalue(min.length)==""|svalue(max.length)==""|svalue(min.five.flank)==""|
		svalue(max.five.flank)==""|svalue(min.three.flank)==""|svalue(max.three.flank)=="")
		{	
			gmessage("Please complete all fields in this form",title="Error!")
		}
		else
		{
			dispose(custom)
		}
	})
		
	visible(custom,TRUE)
})
	
	visible(default.primer.question,TRUE)




# Checks if pdb.id has been entered, if it hasn't, use sequence blast data.
if(nchar(sig.pdb.ids)==0)

{
	# Add BLAST search of pdb database and use top result from this as input to rest of script
	dna.seq<-"ATGTCTAGCCAAATTCGCCAGAATTACAGCACCGACGTTGAAGCGGCAGTCAACAGCCTGGTTAATCTGTACTTGCAGGCCAGCTATACGTATCTGAGCCTGGGCTTTTACTTTGACCGCGACGATGTGGCCTTGGAAGGCGTGAGCCACTTTTTCCGTGAGCTGGCGGAAGAGAAACGCGAAGGCTATGAGCGCCTGCTGAAAATGCAGAACCAACGTGGCGGTCGTGCTCTGTTCCAAGACATCAAGAAACCGGCGGAAGATGAGTGGGGTAAAACCCCGGATGCGATGAAGGCCGCAATGGCTTTGGAGAAGAAACTGAATCAGGCACTGCTGGATCTGCACGCGCTGGGTTCCGCACGTACCGACCCGCACCTGTGCGATTTCTTGGAAACGCATTTTCTGGACGAAGAGGTCAAGCTGATCAAGAAAATGGGCGACCACCTGACGAACTTGCATCGTCTGGGTGGTCCAGAGGCGGGTCTGGGTGAGTACCTGTTCGAGCGTCTGACTCTGAAGCATGATCCCGGG"
	dna.string<-DNAString(dna.seq)
	prot.seq<-translate(dna.string)
	seq<-as.character(prot.seq)

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
	blast.scores<-numeric()
	blast.align.length<-numeric()
	blast.identity<-numeric()
	blast.positive<-numeric()
	blast.hit.length<-numeric()
	ids<-character()
	i<-1
	while(i<=xmlSize(blast.hits))
	{
		ids[[i]] <- xmlValue(blast.hits[[i]][["Hit_id"]])
		blast.scores[[i]]<-xmlValue(blast.hits[[i]][["Hit_hsps"]][["Hsp"]][["Hsp_score"]])
		blast.align.length[[i]]<-xmlValue(blast.hits[[i]][["Hit_hsps"]][["Hsp"]][["Hsp_align-len"]])
		blast.positive[[i]]<-xmlValue(blast.hits[[i]][["Hit_hsps"]][["Hsp"]][["Hsp_positive"]])
		blast.hit.length[[i]]<-xmlValue(blast.hits[[i]][["Hit_len"]])
		i<-i+1 
	}
	blast.scores<-as.numeric(blast.scores)
	blast.align.length<-as.numeric(blast.align.length)
	blast.positive<-as.numeric(blast.positive)
	percentage.hits<-(blast.positive/blast.align.length)*100
	blast.hit.length<-as.numeric(blast.hit.length)

	ids1<-strsplit(ids,"|",fixed=TRUE)

	i<-1
	pdb.ids<-character()
	while(i<=xmlSize(blast.hits))
	{
		pdb.ids[[i]]<-ids1[[i]][[4]]
		i<-i+1
	}

	blast.results<-data.frame(pdb.ids,blast.scores)


	#Selecting statistically significant proteins for analysis

	sig.pdb.ids<-character()
	i<-1
	while(percentage.hits[i] >= 95)
	{
		if(((nchar(seq)-1)-blast.hit.length[i])<=5&blast.hit.length[i]-(nchar(seq)-1)<=5)
		{
			sig.pdb.ids[i]<-pdb.ids[i]
		}
	
	i<-i+1
	}
}
	

# General import pdb - User enters pdb id (or is obtained from blast) - script finds relevant url. 
aminoacid.match2.store<-character()
master.i<-1

	pdb.url<-sub("___",sig.pdb.ids[master.i],"http://www.rcsb.org/pdb/files/___.pdb1",fixed=TRUE)

	x<-read.pdb(pdb.url,multi=TRUE,rm.alt=FALSE)


	### Fixes residue numbers 
	
	# Converts seq to a character vector with each letter seperated
	seq.as.character<-seq.split(seq)
	
	# Converts seq to a three letter sequence 
	three.letter.seq<-one.to.three.translate(seq.as.character)
	
	#Goes through first 100 amino acid identifiers in x$atom and compares them with those in three.letter.seq. If any
	#differ then it will change the fix the residue number identifiers so both residue numbers are the same as seq.
	#Also adds one to residue numbers if the first entry in x$atom is not a methionine.  
	
	residue.numbers<-as.numeric(x$atom[,6])
	i<-1
	p<-as.numeric(x$atom[1,6])-1
	while(i<=150)
	{
		if(as.numeric(x$atom[i,6])>p)
		{
			if(x$atom[i,4]!=three.letter.seq[residue.numbers[p]])
			{	
				residue.numbers<-residue.numbers-residue.numbers[1]+1
				i<-150
				if(x$atom[1,4]!="MET")
				{
					residue.numbers<-residue.numbers+1
				}
			}
			else p<-p+1
		}
		i<-i+1	
	}


	###### Check that structure has correct x$atom[,5] formatting and if not correct ###### 
	# Find and store correct chain names. j is an index with x$atom[i,1]. Since x$atom[i,1] resets when there is
	# a new chain in the pdb file, j will no longer equal x$atom[i,1]. This triggers p to increase and j to be set
	# as equal to x$atom[i,1]. p acts as an index for the chain library which contains all the letters of the 
	# alphabet. Since p only increases when x$atom[i,1] resets, p can be used to create a new library containing 
	# the correct chain identifiers.
	chain.lib<-toupper(paste(letters[1:26]))

	new.chain.store<-character()
	n<-length(x$atom[,1])
	i<-1
	j<-1	
	p<-1
	while(i<=n)
	{
		if(j!=x$atom[i,1])
		{
			
			p<-p+1
			j=x$atom[i,1]
			j<-as.numeric(j)
			if(x$atom[i,1]==1)
			{
				j<-1
			}
		}

		new.chain.store[i]<-chain.lib[p]
		j<-j+1
		i<-i+1
	}
	
	

	# Output of change without duplicates
	chains<-unique(new.chain.store)

	# Assigns each amino acid its specific chain.
	reference.ids<-numeric(0)
	i<-1
	n<-length(x$atom[,6])
	while (i<=n) 
	{
		reference.ids[i]<-paste(c(new.chain.store[i],"-",x$atom[i,6]),collapse="")
		i<-i+1
	}

	primer.reference.ids<-numeric(0)
	i<-1
	n<-length(x$atom[,6])
	while (i<=n) 
	{
		primer.reference.ids[i]<-paste(c(new.chain.store[i],"-",residue.numbers[i]),collapse="")
		i<-i+1
	}
	# one contains only the data from the first chain and none of the others.
	one<-chains[1]
	one<-grep(x=((new.chain.store)==one),pattern=TRUE)

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
			#If distance < 4.5
			if (proximity < 4.5)
			{
				# residue contains each atom closer than 4 angstroms apart and what chain they are on.
				# Which residue has a hydrophobic interaction
				residue<-paste(c(new.chain.store[one][hydrophobesr1][j],"-",x$atom[,6][one][hydrophobesr1][j]),collapse="")
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
				residue2<-paste(c(new.chain.store[-one][hydrophobesr][i],"-",x$atom[,6][-one][hydrophobesr][i]),collapse="")
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

	chain.match<-new.chain.store[residue.match.store]
	
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
	
	aminoacid.match2.store<-c(aminoacid.match2.store,aminoacid.match2)
	

	#Does the same but for primer mutation codes
	primer.atom.aminoacid.match<-numeric(0)
	primer.residue.no.store<-numeric(0)
	p<-1
	while(p<=length(residue.match.store))
	{
		primer.residue.no.store<-c(primer.residue.no.store,residue.numbers[residue.match.store])
		primer.atom.aminoacid.match<-c(primer.atom.aminoacid.match,paste(c(x$atom[residue.match.store,4][p],"-",residue.numbers[residue.match.store][p]),collapse=""))
		p<-p+1
	}
	
	primer.residue.no.store<-unique(primer.residue.no.store)
	primer.aminoacid.match<-unique(primer.atom.aminoacid.match)
	print(aminoacid.match)
	

	
	##### Create .cmd script to open pdb file of protein in Chimera and highlight interacting hydrophobic residues #####
	
	open<-sub("___",sig.pdb.ids[master.i],"open ___")
	add.sym<-"sym #0"
	
	# Makes residue string compatible with chimera command line by deleting spaces between residue numbers
	
	residue.no.store<-gsub(", ",",",residue.no.store,fixed=TRUE)
	
	colour<-sub("___",residue.no.store,"colour red :___")
	
	show<-sub("___",residue.no.store,"show :___")
	
	compile<-c(open,add.sym,colour,show)
	
	file.name<-sub("___",sig.pdb.ids[master.i],"___.cmd")
	write(x=compile,file=file.name)

	##### Uses PrimerX to find mutagenesis primers for breaking apart protein.
	
	# Generates mutation codes for primerX
	one.letter.aminoacid.match<-hydrophobic.3to1.translate(primer.aminoacid.match)

	mutation.codes<-character()
	mutation.codes<-paste(one.letter.aminoacid.match,"R",sep="")


	# Generates new_AA_sequence (because parsing is difficult in this instance)
	i<-1
	new.seq.store<-rep(seq,length(primer.residue.no.store))
	while(i<=length(primer.residue.no.store))
	{
		substr(new.seq.store[i],primer.residue.no.store[i],primer.residue.no.store[i])<-"R"
		i<-i+1
	}
	

	if(master.i==1)
	{
		i<-1
		while(i<=length(mutation.codes))
		{
			primer.finder<-postForm("http://www.bioinformatics.org/primerx/cgi-bin/protein_4.cgi",
			orig_DNA_sequence=dna.seq,
			chopped_DNA_sequence="",
			new_AA_sequence=new.seq.store[i],
			species=species,
			"min_Tm"=min_Tm,
			"max_Tm"=max_Tm,
			min_GC=min_GC,
			max_GC=max_GC,
			min_length=min_length,
			max_length=max_length,
			min_5p_flank=min_5p_flank,
			max_5p_flank=max_5p_flank,
			min_3p_flank=min_5p_flank,
			max_3p_flank=max_3p_flank,
			ends_in_GC=ends_in_GC,
			mut_at_center=mut_at_center,
			primer_type=primer_type,
			"Generate primers"="Generate primers",
			protocol="Basic")
	
			primer.html.file.name<-sub("___",mutation.codes[i],"___ mutagenesis primers.html")
			write(primer.finder,file=primer.html.file.name)
			i<-i+1
		}
	}
