#============================================================================================
# Title 	 	: Analysis of malaria microscopy Systematic Review dataset (AJTMH)
# Version	 	: R version 3.6.3 (2020-02-29)
# Data version		: Recieved from Debashish Das on 19-06-2020
# Script Date		: 19-06-2020
# Analysis script	: Prabin Dahal (prabin.dahal@wwarn.org)
# Data source 		: Supplemental file 3 on the manuscript
#============================================================================================
#rm(list=ls())
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(ggalt)
library(RColorBrewer)
library(cowplot) 	
library(maps)
library(rworldmap)
library(rworldxtra)
library(tableone)

# Set working directory
setwd("C:/Users/")

#===============================================================================
# Section 1: Data management and prepare analysis dataset
#===============================================================================
#----------------------
# Read raw dataset
#----------------------
all_dat <- read.csv("Supplemental file 3.csv")

#-------------------------------------------------
# Parasite species studied in each of the studies
#-------------------------------------------------
all_dat$pf 	<- ifelse(all_dat$Plasmodium.falciparum.included.in.the.study.=="Yes", "pf", "")
all_dat$pv 	<- ifelse(all_dat$Plasmodium.vivax.included.in.the.study.=="Yes", "pv", "")
all_dat$po 	<- ifelse(all_dat$Plasmodium.ovale.included.in.the.study.=="Yes", "po", "")
all_dat$pm 	<- ifelse(all_dat$Plasmodium.malariae.included.in.the.study.=="Yes", "pm", "")
all_dat$pk 	<- ifelse(all_dat$Plasmodium.knowlesi.included.in.the.study.=="Yes", "pk", "")
all_dat$pmix 	<- ifelse(all_dat$Mixed.Infection..Pf.Pv..included.in.the.study.=="Yes", "p_mix", "")

# A single column with parasites species studied
all_dat$species <- paste(all_dat$pf, all_dat$pv,all_dat$po,all_dat$pm,all_dat$pk,all_dat$pmix, sep="+")

#------------------------------------------------------
# Recoding the parasites species as discussed with DD
#------------------------------------------------------
all_dat$species1[all_dat$species=="pf+++++"] <- "pf"
all_dat$species1[all_dat$species=="+pv++++"] <- "pv"
all_dat$species1[all_dat$species=="pf+pv++++"] <- "Both"
all_dat$species1[all_dat$species=="pf+pv++++p_mix"] <- "Both"
all_dat$species1[all_dat$species=="pf+++++p_mix"] <- "Both"
all_dat$species1[all_dat$species=="+pv++++p_mix"] <- "Both"
all_dat$species1[all_dat$species=="pf++po+pm++p_mix"] <- "Other"
all_dat$species1[all_dat$species=="pf+pv+po+pm++"] <- "Other"
all_dat$species1[all_dat$species=="pf+pv+po+pm+pk+"] <- "Other"
all_dat$species1[all_dat$species=="++++pk+"] <- "Other"
all_dat$species1[all_dat$species=="+++++"] <- ""

#-------------------
# Study meta-data
#-------------------
meta_data  <- all_dat %>% 
	  filter(is.na(Repeat.Instance))

meta_data <- subset(meta_data, 
				select=c(
						"Accession.Number.of.Publication",
						"Journal.of.Publication",
						"Region.of.Study",
						"Country.of.Study",
						"Language.of.the.full.text.PDF",
						"First.Year.of.Patient.Recruitment",
						"Last.Year.of.Patient.Recruitment"
						)
				)
#-----------------------------------
# Declaring a slide as negative
#-----------------------------------
slide <- all_dat %>% 
	  filter(is.na(Repeat.Instance))

slide$slide_neg_method <- ifelse(slide$Method.to.Declare.the.Slide.Negative=="High Power Fields (HPF)",
				paste(slide$Method.to.Declare.the.Slide.Negative,slide$Number.of.HPFs.Examined,slide$Number.of.HPFs.Counted...Specify.Other, sep="."),
					ifelse(slide$Method.to.Declare.the.Slide.Negative=="White Blood Cells (WBC)",
				paste(slide$Method.to.Declare.the.Slide.Negative,slide$Number.of.WBCs.Counted,slide$Number.of.WBCs.Counted...Specify.Other, sep="."),
			"Not Stated"
		))

# Method for estimating parasite density
slide$parasite_density_method <- ifelse(slide$Parasite.Density.Estimation=="White Blood Cells (WBC)",
							paste(slide$Parasite.Density.Estimation, slide$WBC.Method,slide$WBC.Method...Specify.Other,sep="."),
							as.character(slide$Parasite.Density.Estimation)
							)		
slide <- subset(slide, select=c("Accession.Number.of.Publication","slide_neg_method","parasite_density_method"))

#----------------------------------------------------------------------------
# Calculate sample size for each study by summing for each treatment arms
#----------------------------------------------------------------------------
samplesize <- all_dat	%>% 
		dplyr::group_by(Accession.Number.of.Publication) %>%
		dplyr::summarise(
			samplesize = sum(Number.of.Patients.Recruited.in.this.Treatment.Arm, na.rm=T)
		)

#-------------------------
# Create analysis dataset
#-------------------------
dat <- all_dat %>% 
	  filter(is.na(Repeat.Instance))

dat$Journal.of.Publication 		<- NULL
dat$Region.of.Study 			<- NULL
dat$Language.of.the.full.text.PDF 	<- NULL
dat$First.Year.of.Patient.Recruitment 	<- NULL
dat$Last.Year.of.Patient.Recruitment 	<- NULL
dat$Repeat.Instrument 			<- NULL
dat$Repeat.Instance 			<- NULL
dat$Country.of.Study 			<- NULL

dat <- merge(dat ,meta_data, by="Accession.Number.of.Publication")
dat <-merge(dat, samplesize , by="Accession.Number.of.Publication")
dat <-merge(dat , slide , by="Accession.Number.of.Publication")

dat$Region.of.Study <- ifelse(dat$Region.of.Study=="Central America", "The Americas",
					ifelse(dat$Region.of.Study=="South America", "The Americas", as.character(dat$Region.of.Study))
					)	
dat <-droplevels(dat)

#-----------------------------------------------------
# Tidy  up country names using regular expression
#-----------------------------------------------------
dat$Country<-sapply(strsplit(as.character(dat$Country.of.Study), " \\("), function(x){x[1]})

#-----------------------------------------
# Plot study locations by parasite species
#-----------------------------------------
# Lat and Lon for the study locations
all_dat_species<- all_dat[!is.na(all_dat$Latitude.of.this.location |all_dat$Longitude.of.this.location),]
all_dat_species$species1 <-NULL
 
species <- all_dat[!duplicated(all_dat$Accession.Number.of.Publication),]
species <- subset(species , select=c("Accession.Number.of.Publication", "species1"))

# Merge parasite species and Lat and Lon for the study locations
all_dat_species <- merge (all_dat_species, species, by="Accession.Number.of.Publication")

#--------------------------------
# Generate Supplemental Figure 1
#--------------------------------

worldmap <- getMap(resolution = "high")
w_map<- worldmap
w_map<- worldmap[which(worldmap$REGION %in% c("Africa","Asia","Australia","North America","South America and the Caribbean","Europe")),]               
plot(w_map, lwd=1,cex.main=1, main= "")
points(all_dat_species[which(all_dat_species$species1=="pf"),]$Longitude.of.this.location,all_dat_species[which(all_dat_species$species1=="pf"),]$Latitude.of.this.location,pch=21,bg="#2c7fb1", cex=0.8)
points(all_dat_species[which(all_dat_species$species1=="pv"),]$Longitude.of.this.location,all_dat_species[which(all_dat_species$species1=="pv"),]$Latitude.of.this.location,pch=21,bg="red", cex=0.8, col="red")
points(all_dat_species[which(all_dat_species$species1=="Both"),]$Longitude.of.this.location,all_dat_species[which(all_dat_species$species1=="Both"),]$Latitude.of.this.location,pch=21,bg="seagreen4", col="seagreen4",cex=0.8)
points(all_dat_species[which(all_dat_species$species1=="Other"),]$Longitude.of.this.location,all_dat_species[which(all_dat_species$species1=="Other"),]$Latitude.of.this.location,pch=21,bg="purple", col="purple",cex=0.8)
legend( "bottomright", c("P.f","P.v","Both P.f and P.v","Other"), col=c("#2c7fb1","red","seagreen4","purple"), pch=19)

#################################################################
# Section 2: Generate summary tables 
#################################################################

#--------------------------
# Unique number of studies 
#--------------------------
dat   %>% 
	dplyr::summarise(
	number_of_trials = length(unique(Accession.Number.of.Publication))
	)

#-------------------------------
# Regional breakdown of studies
#-------------------------------
region<- dat   %>% 
		dplyr::group_by(Region.of.Study) %>%
		dplyr::summarise(
			number_of_trials = length(unique(Accession.Number.of.Publication))
		)

#--------------------------------
# Generate Supplemental Figure 2
#--------------------------------
ggbarplot(region, 
		x 	= "Region.of.Study",
		y 	= "number_of_trials",
		fill 	= "#00AFBB",
		sort.val= "desc",
      	      	add 	= "segments",                        
	        rotate	= TRUE,                              
	        dot.size= 6,                                 
	        font.label = list(color = "white", size = 9, 
                             vjust = 0.5),               
           ggtheme = theme_pubr()                        
           )+
	ylab("Number of studies") +
	ylim(0,110)+
	xlab("")+
	ggtitle("Regional breakdown of number of articles")+
		theme(axis.text.x = element_text(size=7, angle=0),
		          axis.text.y = element_text(size=7),
			 axis.title=element_text(size=8,face="bold"))

#-----------------------------
# Study design parameters
#-----------------------------
# randomisation
dat   %>% 
		dplyr::group_by(Was.the.Study.Randomized.) %>%
		dplyr::summarise(
			number_of_trials = length(unique(Accession.Number.of.Publication))
		)

# blinded
dat   %>% 
		dplyr::group_by(Was.the.Study.Blinded.) %>%
		dplyr::summarise(
			number_of_trials = length(unique(Accession.Number.of.Publication))
		)
# species
dat   %>% 
		dplyr::group_by(species1) %>%
		dplyr::summarise(
			number_of_trials = length(unique(Accession.Number.of.Publication))
		)

# grouped parasite species
dat   %>% 
		dplyr::group_by(species1) %>%
		dplyr::summarise(
			number_of_trials = length(unique(Accession.Number.of.Publication))
		)

#-----------------------------------------
# Method.to.Declare.the.Slide.Negative
#-----------------------------------------
dat %>% 
		dplyr::group_by(Method.to.Declare.the.Slide.Negative) %>%
		dplyr::summarise(
			number_of_trials = length(unique(Accession.Number.of.Publication))
		)

# method to decalre slide negative by region 
CreateTableOne(
		vars 		= c("Method.to.Declare.the.Slide.Negative" ), 
		factorVars 	= c("Method.to.Declare.the.Slide.Negative"), 
		#strata 	= "Region.of.Study", 
		data 		= dat, 
		test		= FALSE
		)

# method to decalre slide negative by species 
CreateTableOne(
		vars 		= c("Method.to.Declare.the.Slide.Negative" ), 
		factorVars 	= c("Method.to.Declare.the.Slide.Negative"), 
		strata 	= "species1", 
		data 		= dat, 
		test		= FALSE
		)
#-----------------------------------------
# Method.to.estimate parasite density
#-----------------------------------------
dat %>% 
		dplyr::group_by(Parasite.Density.Estimation) %>%
		dplyr::summarise(
			number_of_trials = length(unique(Accession.Number.of.Publication))
		)

# method to parasite density estimation by region 
CreateTableOne(
		vars 		= c("Parasite.Density.Estimation" ), 
		factorVars 	= c("Parasite.Density.Estimation"), 
		data 		= dat, 
		test		= FALSE
		)

# method to parasite density estimation by species 
CreateTableOne(
		vars 		= c("Parasite.Density.Estimation" ), 
		factorVars 	= c("Parasite.Density.Estimation"), 
		data 		= dat, 
		test		= FALSE
		)

#-----------------------------------------
# QC parameters
#-----------------------------------------
dat %>% 
		dplyr::group_by(QA.QC.Done) %>%
		dplyr::summarise(
			number_of_trials = length(unique(Accession.Number.of.Publication))
		)

#QA/QC by region 
CreateTableOne(
		vars 		= c("QA.QC.Done" ), 
		factorVars 	= c("QA.QC.Done"), 
		#strata 	= "Region.of.Study", 
		data 		= dat, 
		test		= FALSE
		)

#QA/QC by species 
CreateTableOne(
		vars 		= c("QA.QC.Done" ), 
		factorVars 	= c("QA.QC.Done"), 
		strata 		= "species1", 
		data 		= dat, 
		test		= FALSE
		)

# Where was QA/QC Done among studies which reported QAQC being carried out
dat %>% 
		dplyr::group_by(QA.QC.Done...Where) %>%
			  filter(QA.QC.Done=="Yes")%>%
				dplyr::summarise(
			number_of_trials = length(unique(Accession.Number.of.Publication))
		)

#-----------------------------------------
# Staining Method
#-----------------------------------------
dat %>% 
		dplyr::group_by(Staining.Method) %>%
		dplyr::summarise(
			number_of_trials = length(unique(Accession.Number.of.Publication))
		)

# Staining method by region
CreateTableOne(
		vars 		= c("Staining.Method" ), 
		factorVars 	= c("Staining.Method"), 
		strata 	= "Region.of.Study", 
		data 		= dat, 
		test		= FALSE
		)

# Staining method by species
CreateTableOne(
		vars 		= c("Staining.Method" ), 
		factorVars 	= c("Staining.Method"), 
		strata 	= "species1", 
		data 		= dat, 
		test		= FALSE
		)

#===========================================================
# Section 3: Generation of Figure 2 on the manuscript
#===========================================================
dat <- dput(structure(list(group = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L
), .Label = c("density", "negative"), class = "factor"), Method = structure(c(1L, 
3L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 6L, 2L, 4L, 4L, 4L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L), .Label = c("HPF", "Not stated", "RBC", 
"WBC", "WBC or RBC", "WBC and RBC"), class = "factor"), subgroup = structure(c(11L, 
13L, 18L, 23L, 21L, 20L, 19L, 22L, 17L, 15L, 16L, 24L, 12L, 3L, 
10L, 2L, 1L, 4L, 9L, 7L, 6L, 8L, 5L), .Label = c("After 100 HPFs examined", 
"After 100 WBCs counted", "After 1000 WBCs counted", "After 200 HPFs examined", 
"After 30 HPFs examined", "After 300 HPFs examined", "After 400 HPFs examined", 
"After 50 HPFs examined", "After 500 HPFs examined", "After 500 WBCs counted", 
"High Power Fields (HPF) Method", "Not stated", "RBC: per 1000", 
"WBC or RBC", "WBC: actual counts", "WBC: not clear", "WBC: per 100", 
"WBC: per 200", "WBC: per 200 or per 1000", "WBC: per 200 or per 300", 
"WBC: per 200 or per 500", "WBC: per 300", "WBC: per 500", "WBC and RBC"
), class = "factor"), Number = c(4, 1, 76, 8, 25, 4, 1, 3, 1, 
2, 5, 10, 66, 20, 4, 1, 38, 23, 3, 1, 1, 1, 1)), row.names = c(NA, 
23L), class = "data.frame")
)
density <- dat[which(dat$group=="density"),]
negative <-  dat[which(dat$group=="negative"),]

#---------------------------------------------------------
# Table 1 data: Methods for estimating parasite density
#---------------------------------------------------------

a<-ggdotchart(density, 
		x = "subgroup",
		y = "Number",
	        color = "Method",                             # Color by groups
      	        sorting = "descending",                       # Sort value in descending order
      	     	add = "segments",                             # Add segments from y = 0 to dots
	        rotate = TRUE,                                # Rotate vertically
      	     	group = "Method",                             # Order by groups
	        dot.size = 6,                                 # Large dot size
      	     	label = density$Number,                       # Add mpg values as dot labels
	           font.label = list(color = "white", size = 9, 
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
           )+
	ylab("Number of studies") +
	xlab("")+
	ylim(0,80)+
	ggtitle("A: Parasite density estimation")+
		theme(axis.text.x = element_text(size=8, angle=0),
		          axis.text.y = element_text(size=8),
			 axis.title=element_text(size=14,face="bold"))

#---------------------------------------------------------
# Table 3 data: Methods for estimating parasite density
#---------------------------------------------------------
b<-ggdotchart(negative , 
		x = "subgroup",
		y = "Number",
	           color = "Method",                       # Color by groups
	           sorting = "descending",                 # Sort value in descending order
      	     add = "segments",                             # Add segments from y = 0 to dots
	           rotate = TRUE,                          # Rotate vertically
      	     group = "Method",                             # Order by groups
	           dot.size = 6,                           # Large dot size
      	     label = negative $Number,                     # Add mpg values as dot labels
	           font.label = list(color = "white", size = 9, 
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
           )+
	ylab("Number of studies") +
	xlab("")+
	ylim(0,80)+
	ggtitle("B: Method to declare a slide negative")+
		theme(axis.text.x = element_text(size=8, angle=0),
		          axis.text.y = element_text(size=8),
			 axis.title=element_text(size=14,face="bold"))


# Use cowplot libray to create a 2 by 2 panel & export as high resolution graph (600 dpi)

tiff(file="Figure_2_revised.tiff", 
            width=36, 
		height=16, 
		units="cm", 
            pointsize="15", 
		compression = "lzw+p", 
            bg="white",
		res=600, 
		antialias = "none"
	)
plot_grid(a,b)
dev.off() # End export

#===========================================================
# Section 4: Generation of Figure 3 on the manuscript
#===========================================================

dat <- dput(structure(list(group = structure(c(4L, 4L, 4L, 4L, 4L, 2L, 2L, 
2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 3L, 3L, 3L, 3L), .Label = c("Giemsa_duration", 
"Giemsa_stength", "Quality Control", "Staining"), class = "factor"), 
    subgroup = structure(c(12L, 11L, 16L, 13L, 17L, 3L, 6L, 8L, 
    5L, 4L, 17L, 2L, 1L, 7L, 9L, 17L, 14L, 10L, 15L, 17L), .Label = c(" 11 to < 30", 
    " 5 -10", "< 3%", "1:50, v/v", "10%", "3%", "30 to < 45", 
    "4 to < 10%", "45 - 60", "External ", "Field ", "Giemsa ", 
    "Giemsa/Field", "Internal", "Internal and external ", "Leishman", 
    "Not stated"), class = "factor"), Number = c(138L, 2L, 2L, 
    1L, 63L, 12L, 10L, 4L, 26L, 1L, 85L, 11L, 8L, 20L, 5L, 94L, 
    93L, 9L, 26L, 78L)), class = "data.frame", row.names = c(NA, 
-20L)))

stain 		<- dat[which(dat$group=="Staining"),]
giemsa_strenth 	<-  dat[which(dat$group=="Giemsa_stength"),]
giemsa_duration	<-  dat[which(dat$group=="Giemsa_duration"),]
Q_C 		<-  dat[which(dat$group=="Quality Control"),]

#------------------------------------
# Staining Method Used
#------------------------------------

a <- ggbarplot(stain , 
		x = "subgroup",
		y = "Number",
		fill= "#00AFBB",
		sort.val = "desc",
      	        add= "segments",                        
	            rotate= TRUE,                              
	            dot.size	= 6,                                 
	           font.label = list(color = "white", size = 9, 
                             vjust = 0.5),               
           ggtheme = theme_pubr()                        
           )+
	ylab("Number of studies") +
	ylim(0,150)+
	xlab("")+
	ggtitle("A: Staining method used")+
		theme(axis.text.x = element_text(size=7, angle=0),
		          axis.text.y = element_text(size=7),
			 axis.title=element_text(size=8,face="bold"))

#------------------------------------
# Giemsa duration Used
#------------------------------------
b <- ggbarplot(giemsa_duration, 
			x = "subgroup",
			y = "Number",
			fill= "#00AFBB",
			 sort.val = "desc",
      	     		add	= "segments",                        
	           rotate	= TRUE,                              
	           dot.size	= 6,                                 
	           font.label = list(color = "white", size = 9, 
                             vjust = 0.5),               
           ggtheme = theme_pubr()                        
           )+
	ylab("Number of studies") +
	ylim(0,150)+
	xlab("")+
	ggtitle("B: Duration of Giemsa stain/mins")+
		theme(axis.text.x = element_text(size=7, angle=0),
		          axis.text.y = element_text(size=7),
			 axis.title=element_text(size=8,face="bold"))

#------------------------------------
# Giemsa strenth Used
#------------------------------------
c <- ggbarplot(giemsa_strenth , 
			x = "subgroup",
			y = "Number",
			fill	= "#00AFBB",
			sort.val = "desc",
      	     		add= "segments",                        
	           rotate= TRUE,                              
	           dot.size	= 6,                                 
	           font.label = list(color = "white", size = 7, 
                             vjust = 0.5),               
           ggtheme = theme_pubr()                        
           )+
	ylab("Number of studies") +
	ylim(0,150)+
	xlab("")+
	ggtitle("C: Strength of Giemsa stain")+
		theme(axis.text.x = element_text(size=7, angle=0),
		          axis.text.y = element_text(size=7),
			 axis.title=element_text(size=8,face="bold"))
###########################################################################################
# Use cowplot libray to create a 2 by 2 panel & export as high resolution graph (600 dpi)
###########################################################################################

tiff(file="Figure_3_revised.tiff", 
            width=30, 
		height=8, 
		units="cm", 
            pointsize="15", 
		compression = "lzw+p", 
            bg="white",
		res=600, 
		antialias = "none"
	)
plot_grid(a,b,c, nrow = 1)
dev.off() # End export

# End (Not Run)
