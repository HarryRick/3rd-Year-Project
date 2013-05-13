require("gWidgets")
options("guiToolkit"="RGtk2")

dialog<-gbasicdialog("Select DNA Sequence input",do.buttons=FALSE,guiToolkit="RGtk2")
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

default.primer.question<-gbasicdialog("Select primer conditions",do.buttons=FALSE)
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
	species.drop<-gdroplist(species.list,container=species.frame,handler=function(h,...)
	{
		species<<-svalue(species.drop)
	})
	
	Tm.frame<-gframe(horizontal=TRUE,container=main,spacing=10)
	lab2<-glabel("Minimum melting temperature (ºC):",container=Tm.frame)
	min.temp<-gtext(width=50,height=15,container=Tm.frame,handler=function(h,...)
	{
		min_Tm<<-svalue(min.temp)
	})
	max.temp<-gtext(width=50,height=15,container=Tm.frame,handler=function(h,...)
	{
		max_Tm<<-svalue(max.temp)
	})
	lab3<-glabel("Maximum melting temperature (ºC)",container=Tm.frame)
	
	gc.frame<-gframe(container=main,spacing=10)
	lab4<-glabel("              Minimum GC content (%):",container=gc.frame)
	min.gc<-gtext(width=50,height=15,container=gc.frame,handler=function(h,...)
	{
		min_GC<<-svalue(min.gc)
	})
	max.gc<-gtext(width=50,height=15,container=gc.frame,handler=function(h,...)
	{
		max_GC<<-svalue(max.gc)
	})	
	lab5<-glabel("Maximum GC content (%)",container=gc.frame)
	
	length.frame<-gframe(container=main,spacing=10)
	lab4<-glabel("                      Minimum length (bp):",container=length.frame)
	min.length<-gtext(width=50,height=15,container=length.frame,handler=function(h,...)
	{
		min_length<<-svalue(min.length)
	})
	max.length<-gtext(width=50,height=15,container=length.frame,handler=function(h,...)
	{
		max_length<<-svalue(max.length)
	})
	lab5<-glabel("Maximum length (bp)",container=length.frame)
	
	five.flank.frame<-gframe(container=main,spacing=10)
	lab4<-glabel("     Minimum 5\'\ flanking region (bp):",container=five.flank.frame)
	min.five.flank<-gtext(width=50,height=15,container=five.flank.frame,handler=function(h,...)
	{
		min_5p_flank<<-svalue(min.five.flank)
	})
	max.five.flank<-gtext(width=50,height=15,container=five.flank.frame,handler=function(h,...)
	{
		max_5p_flank<<-svalue(max.five.flank)
	})
	lab5<-glabel("Maximum 5' flanking region (bp)",container=five.flank.frame)

	three.flank.frame<-gframe(container=main,spacing=10)
	lab5<-glabel("     Minimum 3\'\ flanking region (bp):",container=three.flank.frame)
	min.three.flank<-gtext(width=50,height=15,container=three.flank.frame,handler=function(h,...)
	{
		min_3p_flank<<-svalue(min.three.flank)
	})
	max.three.flank<-gtext(width=50,height=15,container=three.flank.frame,handler=function(h,...)
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
