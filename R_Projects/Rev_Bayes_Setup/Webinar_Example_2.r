library(RCurl);			#install.packages("RCurl", dependencies=TRUE);
library(tibble);		#install.packages("tibble", dependencies=TRUE);
hell_no <- F;
# CHANGE DIRECTORIES TO WHERE YOU KEEP COMMON SOURCE FILES!!!!
# This will be different on everyone's computers. Within my R_Projects folder, I separate two "common area" folders that
# 	I have other progams use: Common_R_Source_Files (with functions that many routines might use) and Data_for_R (with data
#	sets that many programs might use). However, I also dress like a hobbit....
# Need to know the name of the directory you want as your computer sees it? Use:
#		x <- file.choose()
# & grab a file in the directory of choice.  Excise the file name & that's the pathway name.

local_directory <- "~/Documents/R_Projects/Rev_Bayes_Setup/";				# this is the directory from which you are running this R-program
setwd(local_directory);														# make sure that we are set up to the appropriate directory
common_source_folder <- "~/Documents/R_Projects/Common_R_Source_Files/";	# directory to folder where you keep common source
data_for_R_folder <- "~/Documents/R_Projects/Data_for_R/";					# directory to folder where you keep R-Data files

# load the source code needed for the analyses
source(paste(common_source_folder,"RevBayes_Setup.r",sep="")); 					# Stuck? Use source(file.choose()) & grab it!
source(paste(common_source_folder,"Historical_Diversity_Metrics.r",sep=""));	# Stuck? Use source(file.choose()) & grab it!

# these are the external databases that the program uses to fine-tune PaleoDB data
load(paste(data_for_R_folder,"Rock_Unit_Database.RData",sep=""));  		# data for rock-units including biozanation & superposition
load(paste(data_for_R_folder,"Gradstein_2012_Augmented.RData",sep="")); # refined Gradstein 2012 timescale & biozonations
load(paste(data_for_R_folder,"PaleoDB_Edits.RData",sep=""));			# edits to collections I cannot edit. # refined Gradstein 2012 timescale & biozonations

# load up the relevant databases!
rock_unit_databases <- rock_unit_data;						# information about rock ages, biozonations and sequence stratigraphy
chronostratigraphic_databases <- gradstein_2012_emended;	# information about time scales, zone ages, etc.
time_scale <- chronostratigraphic_databases$time_scale;		# time, time, time.....
paleodb_fixes <- paleodb_fixes;								# edits to PaleoDB data that I cannot enter currrently

analysis_name <- "Flexicalymene";
onset <- "Ordovician";										# the earliest collections to sample: note that if any of your taxa last appear before this, then it will be adjusted
end <- "Ordovician";										# the latest collections to sample: note that if any of your taxa first appear after this, then it will be adjusted
control_taxon <- "Trilobita";								# Control groups for sampling & diversification
zone_taxa <- c("Conodonta","Graptolithina");				# Groups with species used for major biostratigraphic schemes: these will restrict collection ages
basic_environments <- c("marine","unknown");				# "marine", "terr" (terrestrial) & "unknown" are the basic options
time_scale_stratigraphic_scale <- "Stage Slice";			# this will set the stratigraphic scale that we'll use
temporal_precision <- 0.1;									# level of precision for trying to estimate ages;
save_files <- F;											# if true, a lot of files will be output
flexicalymene <- accio_PaleoDB_data_from_chosen_nexus_file(onset=onset,end=end,rock_unit_databases=rock_unit_databases,chronostratigraphic_databases=chronostratigraphic_databases,paleodb_fixes=paleodb_fixes,control_taxon=control_taxon,zone_taxa=zone_taxa,basic_environments=basic_environments,time_scale_stratigraphic_scale=time_scale_stratigraphic_scale,analysis_name=analysis_name,temporal_precision=temporal_precision,save_files=save_files);

bogarted <- T;
flexicalymene <- accio_PaleoDB_data_from_chosen_nexus_file(onset=onset,end=end,rock_unit_databases=rock_unit_databases,chronostratigraphic_databases=chronostratigraphic_databases,paleodb_fixes=paleodb_fixes,control_taxon=control_taxon,zone_taxa=zone_taxa,basic_environments=basic_environments,time_scale_stratigraphic_scale=time_scale_stratigraphic_scale,analysis_name=analysis_name,temporal_precision=temporal_precision,save_files=save_files,bogarted=bogarted);


#bogarted <- T;												# if true, then you wll be prompted for a file of additional occurrences
analysis_name <- "Lycopoda";
onset <- "Silurian";										# the earliest collections to sample: note that if any of your taxa last appear before this, then it will be adjusted
end <- "Carboniferous";										# the latest collections to sample: note that if any of your taxa first appear after this, then it will be adjusted
control_taxon <- "Plantae";									# Control groups for sampling & diversification
zone_taxa <- c("");											# Groups with species used for major biostratigraphic schemes: these will restrict collection ages
basic_environments <- c("terr","marine","unknown");			# "marine", "terr" (terrestrial) & "unknown" are the basic options
time_scale_stratigraphic_scale <- "Stage Slice";			# this will set the stratigraphic scale that we'll use
temporal_precision <- 0.1;									# level of precision for trying to estimate ages;
bogarted <- F;												# if true, then you wll be prompted for a file of additional occurrences
taxon_subset_file <- T;										# if true, then you wll be prompted for a file stating which OTUs to include
salad <- accio_PaleoDB_data_from_chosen_nexus_file(onset=onset,end=end,rock_unit_databases=rock_unit_databases,chronostratigraphic_databases=chronostratigraphic_databases,paleodb_fixes=paleodb_fixes,control_taxon=control_taxon,zone_taxa=zone_taxa,basic_environments=basic_environments,time_scale_stratigraphic_scale=time_scale_stratigraphic_scale,analysis_name=analysis_name,taxon_subset_file=taxon_subset_file,temporal_precision=temporal_precision);

