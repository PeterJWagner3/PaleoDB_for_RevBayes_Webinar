# These code used here rely on functions in the following packages.
#	If you do not have them installed on your computer, simply highlight the "install.packages("I_need_this",dependencies=T)
#	and hit return; this should install the package on your computer.
library(combinat);		#install.packages("combinat", dependencies=TRUE);
library(dplyr);			#install.packages("dplyr", dependencies=TRUE);
library(gdata);			#install.packages("gdata", dependencies=TRUE);
library(gtools);		#install.packages("gtools", dependencies=TRUE);
library(lettercase);	#install.packages("lettercase", dependencies=TRUE);
library(paleobioDB);	#install.packages("paleobioDB", dependencies=TRUE);
library(pracma);		#install.packages("pracma", dependencies=TRUE);
library(prodlim);		#install.packages("prodlim", dependencies=TRUE);
library(raster);		#install.packages("raster", dependencies=TRUE);
library(Rcpp);			#install.packages("Rcpp", dependencies=TRUE);
library(rvest);			#install.packages("rvest", dependencies=TRUE);
library(sads);			#install.packages("sads", dependencies=TRUE);
library(stringr);		#install.packages("stringr", dependencies=TRUE);
library(subplex);		#install.packages("subplex", dependencies=TRUE);

# CHANGE DIRECTORIES TO WHERE YOU KEEP COMMON SOURCE FILES!!!!
# This will be different on everyone's computers. Within my R_Projects folder, I separate two "common area" folders that
# 	I have other progams use: Common_R_Source_Files (with functions that many routines might use) and Data_for_R (with data
#	sets that many programs might use). However, I also dress like a hobbit....
common_source_folder <- "~/Documents/R_Projects/Common_R_Source_Files/";	# directory to folder where you keep common source
data_for_R_folder <- "~/Documents/R_Projects/Data_for_R/";					# directory to folder where you keep R-Data files
setwd("~/Documents/R_Projects/Rev_Bayes_Setup");
local_directory <- "~/Documents/R_Projects/Rev_Bayes_Setup/";				# this is the directory from which you are running this R-program

# load the source code needed for the analyses
source(paste(common_source_folder,"Wagner_Kluges.r",sep=""));
source(paste(common_source_folder,"General_Plot_Templates.r",sep=""));
source(paste(common_source_folder,"paleophylogeny_routines.r",sep=""));
source(paste(common_source_folder,"Data_Downloading_v4.r",sep=""));
source(paste(common_source_folder,"Stratigraphy.r",sep=""));
source(paste(common_source_folder,"RevBayes_Setup.r",sep="")); # source(file.choose())

# these are the external databases that the program uses to fine-tune PaleoDB data
load(paste(data_for_R_folder,"Rock_Unit_Database.RData",sep=""));  		# data for rock-units including biozanation & superposition
load(paste(data_for_R_folder,"Gradstein_2012_Augmented.RData",sep="")); # refined Gradstein 2012 timescale & biozonations
load(paste(data_for_R_folder,"PaleoDB_Edits.RData",sep=""));			# edits to collections I cannot edit. # refined Gradstein 2012 timescale & biozonations

# load up the relevant databases!
rock_unit_databases <- rock_unit_data;						# information about rock ages, biozonations and sequence stratigraphy
chronostratigraphic_databases <- gradstein_2012_emended;	# information about time scales, zone ages, etc.
paleodb_fixes <- paleodb_fixes;								# edits to PaleoDB data that cannot be entered currrently
# now begins the magic.....
#proetids <- accio_PaleoDB_data_from_chosen_nexus_file(analysis_name="Proetidae",taxon_subset_file=F,rate_partition="",trend_partition="",write_data_directory="",write_scripts_directory="",local_directory="",set_wdir="",UNKNOWN=-11,INAP=-22);
cincta <- accio_PaleoDB_data_from_chosen_nexus_file(onset="Cambrian",end="Cambrian",rock_unit_databases=rock_unit_data,chronostratigraphic_databases=chronostratigraphic_databases,paleodb_fixes=paleodb_fixes,control_taxon="Echinodermata",zone_taxa=c("Trilobita","Archaeocyatha"),taxon_level="species",basic_environments=c("marine","unknown"),time_scale_stratigraphic_scale="Stage Slice",analysis_name="Cincta",bogarted=F);
write.csv(cincta$occurrences,paste(local_directory,"Cincta_Occurrences_Final.csv",sep=""),row.names=F);
write.csv(cincta$collections,paste(local_directory,"Cincta_Collections_Final.csv",sep=""),row.names=F);
write.csv(cincta$fossil_summaries,paste(local_directory,"Cincta_Stratigraphic_Range_Summaries.csv",sep=""),row.names=F);
write.csv(cincta$relv_time_scale,paste(local_directory,"Cincta_Time_Scale.csv",sep=""),row.names=F);

proetids <- accio_PaleoDB_data_from_chosen_nexus_file(onset="Cambrian",end="Mississippian",rock_unit_databases=rock_unit_data,chronostratigraphic_databases=chronostratigraphic_databases,paleodb_fixes=paleodb_fixes,control_taxon="Trilobita",zone_taxa=c("Brachiopoda","Conodonta","Graptolithina"),basic_environments=c("marine","unknown"),time_scale_stratigraphic_scale="Stage Slice",analysis_name="Proetids",bogarted=T);
write.csv(proetids$occurrences,paste(local_directory,"Proetidae_Occurrences_Final.csv",sep=""),row.names=F);
write.csv(proetids$collections,paste(local_directory,"Proetidae_Collections_Final.csv",sep=""),row.names=F);
write.csv(proetids$fossil_summaries,paste(local_directory,"Proetidae_Stratigraphic_Range_Summaries.csv",sep=""),row.names=F);
write.csv(proetids$relv_time_scale,paste(local_directory,"Proetidae_Time_Scale.csv",sep=""),row.names=F);

rangeomorphs <- accio_PaleoDB_data_from_chosen_nexus_file(onset="Proterozoic",end="Ediacaran",rock_unit_databases=rock_unit_data,chronostratigraphic_databases=chronostratigraphic_databases,paleodb_fixes=paleodb_fixes,control_taxon="Eukaryota",zone_taxa=c(""),basic_environments=c("marine","unknown"),time_scale_stratigraphic_scale="Stage Slice",analysis_name="Rangeomorpha",bogarted=F);
write.csv(rangeomorphs$occurrences,paste(local_directory,"Rangeomorpha_Occurrences_Final.csv",sep=""),row.names=F);
write.csv(rangeomorphs$collections,paste(local_directory,"Rangeomorpha_Collections_Final.csv",sep=""),row.names=F);
write.csv(rangeomorphs$fossil_summaries,paste(local_directory,"Rangeomorpha_Stratigraphic_Range_Summaries.csv",sep=""),row.names=F);
write.csv(rangeomorphs$relv_time_scale,paste(local_directory,"Rangeomorpha_Time_Scale.csv",sep=""),row.names=F);

ammonites <- accio_PaleoDB_data_from_chosen_nexus_file(onset="Oxfordian",end="Coniacian",rock_unit_databases="",chronostratigraphic_databases=chronostratigraphic_databases,paleodb_fixes=paleodb_fixes,control_taxon="Ammonoidea",zone_taxa=c("Inoceramoidea","Foraminifera"),basic_environments=c("marine","unknown"),time_scale_stratigraphic_scale="International",analysis_name="Acanthocertoids",bogarted=F);
write.csv(ammonites$occurrences,paste(local_directory,"Acanthoceratid_Occurrences_Final.csv",sep=""),row.names=F);
write.csv(ammonites$collections,paste(local_directory,"Acanthoceratid_Collections_Final.csv",sep=""),row.names=F);
write.csv(ammonites$fossil_summaries,paste(local_directory,"Acanthoceratid_Stratigraphic_Range_Summaries.csv",sep=""),row.names=F);
write.csv(ammonites$relv_time_scale,paste(local_directory,"Acanthoceratid_Time_Scale.csv",sep=""),row.names=F);

nimrods <- accio_PaleoDB_data_from_chosen_nexus_file(onset="Jurassic",end="Oligocene",rock_unit_databases="",chronostratigraphic_databases=chronostratigraphic_databases,paleodb_fixes=paleodb_fixes,control_taxon=c("Mammalia","Lepidosauromorpha"),zone_taxa=c(""),basic_environments=c("terrestrial","unknown"),time_scale_stratigraphic_scale="NALMA",analysis_name="Nimrods",bogarted=T);
write.csv(nimrods$occurrences,paste(local_directory,"Nimravine_Occurrences_Final.csv",sep=""),row.names=F);
write.csv(nimrods$collections,paste(local_directory,"Nimravine_Collections_Final.csv",sep=""),row.names=F);
write.csv(nimrods$fossil_summaries,paste(local_directory,"Nimravine_Stratigraphic_Range_Summaries.csv",sep=""),row.names=F);
write.csv(nimrods$relv_time_scale,paste(local_directory,"Nimravine_Time_Scale.csv",sep=""),row.names=F);
