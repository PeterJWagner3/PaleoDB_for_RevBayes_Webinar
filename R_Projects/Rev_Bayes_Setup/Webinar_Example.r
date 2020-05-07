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

#### Specify Parameters & Load Up the Relevant Databases!  ####
rock_unit_databases <- rock_unit_data;						# information about rock ages, biozonations and superposition for Paleozoic rock units
chronostratigraphic_databases <- gradstein_2012_emended;	# information about time scales, zone ages, etc.
time_scale <- chronostratigraphic_databases$time_scale;		# time, time, time.....
paleodb_fixes <- paleodb_fixes;								# edits to PaleoDB data that I cannot enter currrently
analysis_name <- "Cincta";
onset <- "Cambrian";										# the earliest collections to sample: note that if any of your taxa last appear before this, then it will be adjusted
end <- "Cambrian";											# the latest collections to sample: note that if any of your taxa first appear after this, then it will be adjusted
time_scale_stratigraphic_scale <- "Stage Slice";			# The stratigraphic scale that you want to use for the analyses. Default is "International" (= "Standard)
control_taxon <- "Echinodermata";							# Control groups for sampling & diversification
zone_taxa <- c("Trilobita","Agnostida","Archaeocyatha");	# Groups with species used for major biostratigraphic schemes: these will restrict collection ages
basic_environments <- c("marine","unknown");				# "marine", "terr" (terrestrial) & "unknown" are the basic options
temporal_precision <- 0.1;									# level of precision for trying to estimate ages;
save_files <- F;											# if T, then a whole bunch of files will be saved!

#### Download & Refine PaleoDB Data  ####
cinctans <- accio_PaleoDB_data_from_chosen_nexus_file(onset=onset,end=end,rock_unit_databases=rock_unit_databases,chronostratigraphic_databases=chronostratigraphic_databases,paleodb_fixes=paleodb_fixes,control_taxon=control_taxon,zone_taxa=zone_taxa,basic_environments=basic_environments,time_scale_stratigraphic_scale=time_scale_stratigraphic_scale,analysis_name=analysis_name,temporal_precision=temporal_precision,save_files=save_files);
write.csv(cinctans$occurrences,paste(local_directory,analysis_name,"_Occurrences_Final.csv",sep=""),row.names=F,fileEncoding = "UTF-8");
write.csv(cinctans$collections,paste(local_directory,analysis_name,"_Collections_Final.csv",sep=""),row.names=F,fileEncoding = "UTF-8");
write.csv(cinctans$fossil_summaries,paste(local_directory,analysis_name,"_Stratigraphic_Range_Summaries.csv",sep=""),row.names=F,fileEncoding = "UTF-8");
write.csv(cinctans$relv_time_scale,paste(local_directory,analysis_name,"_Time_Scale.csv",sep=""),row.names=F,fileEncoding = "UTF-8");

#### Output Taxon File with Chronostratigraphic Information for OTUs ####
output_precision <- ceiling(abs(log10(temporal_precision)));
strat_range_summaries <- cinctans$fossil_summaries;
taxon_name <- strat_range_summaries$taxon;
taxon_names <- sapply(taxon_name,nexusify_taxon_name);
# We'll generate 5 files
# 1st: just the upper & lower bounds of the stratigraphic range using millions of years ago
fossil_intervals_for_RevBayes_raw <- data.frame(taxon=as.character(taxon_names),min=round(as.numeric(strat_range_summaries$latest_poss_la),output_precision),max=round(as.numeric(strat_range_summaries$earliest_poss_fa),output_precision),stringsAsFactors = hell_no);
# 2nd: as 1st, but with the dates reset to the latest possible date
faux_present <- round(min(abs(fossil_intervals_for_RevBayes_raw$min)),output_precision);	# set end of study to the "present"
fossil_intervals_for_RevBayes <- data.frame(taxon=as.character(taxon_names),min=round(as.numeric(strat_range_summaries$latest_poss_la-faux_present),output_precision),max=round(as.numeric(strat_range_summaries$earliest_poss_fa-faux_present),output_precision),stringsAsFactors = hell_no);
# 3rd: as 2nd, but using only first appearances and latest first appearance now the "Recent"
faux_present2 <- min(abs(fossil_intervals_for_RevBayes_raw$max));	# set last first appearance to "present
fossil_intervals_for_RevBayes_FA_only <- data.frame(taxon=as.character(taxon_names),min=round(as.numeric(strat_range_summaries$earliest_poss_fa-faux_present2),output_precision),max=round(as.numeric(strat_range_summaries$earliest_poss_fa-faux_present2),output_precision),stringsAsFactors = hell_no);
# 4th: upper and lower bounds for both first and last appearances, again set to the latest possible LA is 0. 
fossil_intervals_for_RevBayes_fuzzy <- data.frame(taxon=as.character(taxon_names),min_mn=round(as.numeric(strat_range_summaries$latest_poss_la-faux_present),output_precision),min_mx=round(as.numeric(strat_range_summaries$earliest_poss_la-faux_present),output_precision),max_mn=round(as.numeric(strat_range_summaries$latest_poss_fa-faux_present),output_precision),max_mx=round(as.numeric(strat_range_summaries$earliest_poss_fa-faux_present),output_precision),stringsAsFactors = hell_no);
# 5th: upper & lower bounds for different rock units in which taxa occur
fossil_intervals_for_RevBayes_occurrences <- cinctans$otu_rock_occurrences;
taxon_name <- cinctans$otu_rock_occurrences$taxon;
fossil_intervals_for_RevBayes_occurrences$taxon <- sapply(taxon_name,nexusify_taxon_name);
fossil_intervals_for_RevBayes_occurrences$min <- round(fossil_intervals_for_RevBayes_occurrences$min-faux_present,output_precision);
fossil_intervals_for_RevBayes_occurrences$max <- round(fossil_intervals_for_RevBayes_occurrences$max-faux_present,output_precision);

fossil_interval_file <- paste(tolower(analysis_name),"_fossil_intervals.tsv",sep="");
fossil_interval_file_FA <- paste(tolower(analysis_name),"_fossil_intervals_FA.tsv",sep="");
fossil_interval_file_raw <- paste(tolower(analysis_name),"_fossil_intervals_orig_dates.tsv",sep="");
fossil_interval_file_fuzzy <- paste(tolower(analysis_name),"_fossil_intervals_fuzzy.tsv",sep="");
fossil_interval_file_finds <- paste(tolower(analysis_name),"_fossil_intervals_finds.tsv",sep="");
	
write.table(fossil_intervals_for_RevBayes_raw,file=paste(local_directory,fossil_interval_file_raw,sep=""),row.names = F,sep="\t",quote = F);
write.table(fossil_intervals_for_RevBayes,file=paste(local_directory,fossil_interval_file,sep=""),row.names = F,sep="\t",quote = F);
write.table(fossil_intervals_for_RevBayes_FA_only,file=paste(local_directory,fossil_interval_file_FA,sep=""),row.names = F,sep="\t",quote = F);
write.table(fossil_intervals_for_RevBayes_fuzzy,file=paste(local_directory,fossil_interval_file_fuzzy,sep=""),row.names = F,sep="\t",quote = F);
write.table(fossil_intervals_for_RevBayes_occurrences,file=paste(local_directory,fossil_interval_file_finds,sep=""),row.names = F,sep="\t",quote = F);

#### Get Sampling Metrics ####
finest_chronostrat <- cinctans$relv_time_scale;
finest_chronostrat$span <- finest_chronostrat$ma_lb-finest_chronostrat$ma_ub;
taxon_sites <- cinctans$collections;
interval_sites <- tally_collections_occupied_by_subinterval(taxon_sites,hierarchical_chronostrat = finest_chronostrat);
interval_rocks <- tally_rock_units_occupied_by_subinterval(taxon_sites,finest_chronostrat);

taxon_finds <- cinctans$occurrences;
taxon_species <- sort(unique(taxon_finds$accepted_name));
sites_per_bin_spc <- rocks_per_bin_spc <- c();
for (ts in 1:length(taxon_species))	{
	species_finds <- subset(taxon_finds,taxon_finds$accepted_name==taxon_species[ts]);
	species_sites <- unique(subset(taxon_sites,taxon_sites$collection_no %in% species_finds$collection_no));
	sites_per_bin_spc <- rbind(sites_per_bin_spc,tally_collections_occupied_by_subinterval(taxon_collections=species_sites,hierarchical_chronostrat = finest_chronostrat));
	rocks_per_bin_spc <- rbind(rocks_per_bin_spc,tally_rock_units_occupied_by_subinterval(taxon_collections=species_sites,hierarchical_chronostrat = finest_chronostrat));
	}
rownames(sites_per_bin_spc) <- rownames(rocks_per_bin_spc) <- taxon_species;
sites_per_bin_spc_rnd <- round_fuzzy_finds_per_interval(finds_per_bin = sites_per_bin_spc);
rocks_per_bin_spc_rnd <- round_fuzzy_finds_per_interval(finds_per_bin = rocks_per_bin_spc);

# output this for your own records #
write.csv(interval_sites,paste(local_directory,analysis_name,"_Sites_per_Bin.csv",sep=""));
write.csv(interval_rocks,paste(local_directory,analysis_name,"_Rock_Units_per_Bin.csv",sep=""));

distributed_sampling_over_time_sites <- accio_sampling_distributions_for_RevBayes(finds_per_bin = sites_per_bin_spc_rnd,sample_units_per_bin = ceiling(interval_sites));
distributed_sampling_over_time_rocks <- accio_sampling_distributions_for_RevBayes(finds_per_bin = rocks_per_bin_spc_rnd,sample_units_per_bin = ceiling(interval_rocks));
distributed_sampling_over_time_sites$uniform <- distributed_sampling_over_time_rocks$uniform <- NULL;
all_intervals <- names(interval_rocks);

### get distributions of sampling rates a la Marcot & Wagner 2013 ###
# get per-site probabilities of finding taxa;
quantiles_for_sites <- accio_sampling_quantiles_for_all_intervals(sampling_over_time = distributed_sampling_over_time_sites,all_intervals = all_intervals,criterion="AICc",ttl_quantiles=21);
# per-rock_unit probabilities of finding taxa;
quantiles_for_rocks <- accio_sampling_quantiles_for_all_intervals(sampling_over_time = distributed_sampling_over_time_rocks,all_intervals = all_intervals,criterion="AICc",ttl_quantiles=21);

pfind_sites_quantiles <- pfind_rocks_quantiles <- exp_finds_sites_quantiles <- exp_finds_rocks_quantiles <- quantiles_for_rocks;
nbins <- nrow(quantiles_for_sites);
prob_missing_site <- prob_missing_rock <- vector(length=nbins);
for (i in 1:nbins)	{
	pfind_sites_quantiles[i,] <- simplify2array(1-((1-quantiles_for_sites[i,])^interval_sites[i]));
	prob_missing_site[i] <- mean(as.numeric(1-pfind_sites_quantiles[i,])^interval_rocks[i]);
	pfind_rocks_quantiles[i,] <- simplify2array(1-((1-quantiles_for_rocks[i,])^interval_rocks[i]));
	prob_missing_rock[i] <- mean(as.numeric(1-pfind_rocks_quantiles[i,])^interval_rocks[i]);
	exp_finds_sites_quantiles[i,] <- simplify2array(probability_to_Poisson_rate(pfind_sites_quantiles[i,]));
	exp_finds_rocks_quantiles[i,] <- simplify2array(probability_to_Poisson_rate(pfind_rocks_quantiles[i,]));
	}

# output for your own purposes	
quantiles_for_sites$mean_finds <- rowSums(exp_finds_sites_quantiles)/ncol(exp_finds_sites_quantiles);
quantiles_for_sitesprob_missing <- prob_missing_site;
filename <- paste(local_directory,analysis_name,"_per_Interval_Sampling_given_Sites.csv",sep="");
write.csv(quantiles_for_sites,filename,row.names=F,fileEncoding = "UTF-8");
quantiles_for_rocks$mean_finds <- rowSums(exp_finds_rocks_quantiles)/ncol(exp_finds_rocks_quantiles);
quantiles_for_rocks$prob_missing <- prob_missing_rock;
filename <- paste(local_directory,analysis_name,"_per_Interval_Sampling_given_Rocks.csv",sep="");
write.csv(quantiles_for_rocks,filename,row.names=F,fileEncoding = "UTF-8");

### Get initial sampling rates and final interval sampling rate
study_intervals <- finest_chronostrat$interval[min(strat_range_summaries$latest_poss_fa)<=finest_chronostrat$ma_ub];
study_intervals <- study_intervals[rowMeans(exp_finds_sites_quantiles[study_intervals,])>0]
last_interval <- study_intervals[length(study_intervals)];

# estimate probability of sampling individual from the interval in which the youngest taxon appears!
rho <- (mean(as.numeric(pfind_rocks_quantiles[match(last_interval,rownames(pfind_rocks_quantiles)),]))+mean(as.numeric(pfind_sites_quantiles[match(last_interval,rownames(pfind_rocks_quantiles)),])))/2;
# get the average sampling rate per million years
psi <- mean(((rowMeans(exp_finds_sites_quantiles[study_intervals,])+rowMeans(exp_finds_rocks_quantiles[study_intervals,]))/2)/finest_chronostrat$span[match(study_intervals,finest_chronostrat$interval)]);

#### Get Initial Estimates of Diversification ####
three_timer_setup <- setup_three_timer_analysis(samples_per_interval = sites_per_bin_spc_rnd);
# three_timer_setup$taxon_ranges			# synoptic ranges in bin numbers, with numbers matching the finest_chronostrat timescale
# three_timer_setup$sampled_in_bin;			# number of taxa *sampled* in each interval (≤ number inferred to be present)
# three_timer_setup$sepkoski_richness;		# synoptic richness a la Jack
# three_timer_setup$gappers;				# number of taxa appearing before and after an interval but not in the interval
three_timer_info <- data.frame(sampled_in_bin=as.numeric(three_timer_setup$sampled_in_bin),gappers=as.numeric(three_timer_setup$gappers),synoptic_richness=as.numeric(three_timer_setup$sepkoski_richness),stringsAsFactors=hell_no);
rownames(three_timer_info) <- names(three_timer_setup$sepkoski_richness);
filename <- paste(local_directory,analysis_name,"_Three_Timer_Ingredients.csv",sep="");
write.csv(three_timer_info,file = filename,row.names=T,fileEncoding = "UTF-8");

prob_miss_interval <- prob_missing_site;
first_bin <- min((1:nbins)[prob_miss_interval<1]);
last_bin <- max((1:nbins)[prob_miss_interval<1]);
diversification_ma_rock <- accio_diversification_per_ma_likelihoods(occr_per_bin=rocks_per_bin_spc,prob_miss_interval=prob_miss_interval,time_scale=finest_chronostrat,first_bin,last_bin);#,first_bin=(first_bin_backdrop-1));

prob_miss_interval <- prob_missing_rock;
first_bin <- min((1:nbins)[prob_miss_interval<1]);
last_bin <- max((1:nbins)[prob_miss_interval<1]);
diversification_ma_site <- accio_diversification_per_ma_likelihoods(occr_per_bin=sites_per_bin_spc,prob_miss_interval=prob_miss_interval,time_scale=finest_chronostrat,first_bin,last_bin);#,first_bin=(first_bin_backdrop-1));

filename <- paste(local_directory,analysis_name,"_Diversification_Numbers_given_Rocks.csv",sep="");
write.csv(diversification_ma_rock,filename,row.names = F);
filename <- paste(local_directory,analysis_name,"_Diversification_Numbers_given_Sites.csv",sep="");
write.csv(diversification_ma_site,filename,row.names = F);

#### Get information for RevBayes (This just reads files and can be skipped if you've done the stuff above in this run ####
# reload information 
filename <- paste(local_directory,analysis_name,"_Diversification_Numbers_given_Rocks.csv",sep="");
diversification_ma_rock <- read.csv(filename,header=T,stringsAsFactors = F);
filename <- paste(local_directory,analysis_name,"_Diversification_Numbers_given_Sites.csv",sep="");
diversification_ma_site <- read.csv(filename,header=T,stringsAsFactors = F);
filename <- paste(local_directory,analysis_name,"_per_Interval_Sampling_given_Sites.csv",sep="");
s_s_d <- read.csv(filename,header = T,stringsAsFactors = F);
site_sampling_distribution <- read.csv(filename,header = T,stringsAsFactors = F);
filename <- paste(local_directory,analysis_name,"_per_Interval_Sampling_given_Rocks.csv",sep="");
r_s_d <- read.csv(filename,header = T,stringsAsFactors = F);
rock_sampling_distribution <- read.csv(filename,header = T,stringsAsFactors = F);
filename <- paste(local_directory,analysis_name,"_Three_Timer_Ingredients.csv",sep="");
three_timer_info <- read.csv(file=filename,header = T,stringsAsFactors = F);
rownames(three_timer_info) <- rownames(diversification_ma_rock) <- rownames(diversification_ma_site) <- rownames(rock_sampling_distribution) <- rownames(site_sampling_distribution) <- rock_sampling_distribution$X;
three_timer_info$X <- diversification_ma_rock$X <- diversification_ma_site$X <- site_sampling_distribution$X <- rock_sampling_distribution$X <- NULL;
nbins <- nrow(three_timer_info);

## Recalculate sampling from uploaded files
# sampling in the latest relevant interval
rho <- (mean(as.numeric(pfind_rocks_quantiles[match(last_interval,rownames(pfind_rocks_quantiles)),]))+mean(as.numeric(pfind_sites_quantiles[match(last_interval,rownames(pfind_sites_quantiles)),])))/2;
# average sampling
study_bins <- match(study_intervals,finest_chronostrat$interval);
psi <- mean(((rowMeans(exp_finds_sites_quantiles[study_bins,])+rowMeans(exp_finds_rocks_quantiles[study_bins,]))/2)/finest_chronostrat$span[study_bins]);

#### Now, write the FBD Script ####
## First, let's get initial diversification parameters
print("First, we'll generate a script to parameterize an MCMC analysis using Fossilized-Birth-Death");
Sys.sleep(0.1);
basic_data <- accio_data_from_chosen_nexus_file();
otu_names_used <- basic_data$OTUs;
taxon_names <- otu_names_used[!tolower(otu_names_used) %in% c("outgroup","out")];
otu_names <- sapply(taxon_names,scourgify_taxon_names);
ingroup_taxa <- otu_names[!(1:length(otu_names)) %in% basic_data$Outgroup];
notu <- length(otu_names);
uncoded_taxa <- basic_data$Unscored_Taxa;
if (is.null(uncoded_taxa))
	uncoded_taxa <- "";

obins <- (1:nbins)[diversification_ma_rock$ML_Orig>0][(1:nbins)[diversification_ma_rock$ML_Orig>0] %in% (1:nbins)[three_timer_info$sampled_in_bin>0]]
ebins <- (1:nbins)[diversification_ma_rock$ML_Extn>0][(1:nbins)[diversification_ma_rock$ML_Extn>0] %in% (1:nbins)[three_timer_info$sampled_in_bin>0]]
origination <- median((diversification_ma_rock$ML_Orig[obins]+diversification_ma_site$ML_Orig[obins])/2);
extinction <- median((diversification_ma_rock$ML_Extn[ebins]+diversification_ma_site$ML_Extn[ebins])/2);
turnover <- min(1.05,extinction/origination);							# NOTE: Extinction & Origination are not independent variables; we'll vary turnover instead

## get sampling parameters; if you've run this straight through, then this will be redundant
# Cousin Matthew Parameter (calculated a la Bapst 2013)
phi <- prob_sampling_clade_bapst(p=origination,q=extinction,r=psi);		# probability of sampling a clade of unknown size 
control_taxon <- "Echinodermata";
sampling_unit <- "Rocks & Sites Averaged";
divergence_bounds <- ceiling(max(fossil_intervals_for_RevBayes_FA_only$max))+c(0,ceiling(10/(psi+(origination*phi))));
extant_taxa <- fossil_intervals_for_RevBayes_FA_only$taxon[fossil_intervals_for_RevBayes_FA_only$max==0];

### Output FBD parameterization file ###
fbd_script <- scribio_fbd_portion_of_Rev_Bayes_script(analysis_name,write_scripts_directory="",origination,extinction,psi,rho,divergence_bounds,control_taxon,extant_taxa,otu_names,uncoded_taxa,script_file_lead="scripts/");
fbd_file_name <- fbd_script$filename;

#### Finally, Output MCMC & Stepping Stone Scripts ####
scribio_RevBayes_scripts_from_chosen_nexus_file_and_existing_FBD_script_and_data(analysis_name,local_directory=local_directory,fa_info=fossil_intervals_for_RevBayes_FA_only);


#### the end.... ####
