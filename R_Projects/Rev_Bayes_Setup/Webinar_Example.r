library(RCurl);			#install.packages("RCurl", dependencies=TRUE);
library(tibble);		#install.packages("tibble", dependencies=TRUE);

# CHANGE DIRECTORIES TO WHERE YOU KEEP COMMON SOURCE FILES!!!!
# This will be different on everyone's computers. Within my R_Projects folder, I separate two "common area" folders that
# 	I have other progams use: Common_R_Source_Files (with functions that many routines might use) and Data_for_R (with data
#	sets that many programs might use). However, I also dress like a hobbit....
# Need to know the name of the directory you want as your computer sees it? Use:
#		x <- file.choose()
# & grab a file in the directory of choice.  Excise the file name & that's the pathway name.

common_source_folder <- "~/Documents/R_Projects/Common_R_Source_Files/";	# directory to folder where you keep common source
data_for_R_folder <- "~/Documents/R_Projects/Data_for_R/";					# directory to folder where you keep R-Data files
setwd("~/Documents/R_Projects/Rev_Bayes_Setup");
local_directory <- "~/Documents/R_Projects/Rev_Bayes_Setup/";				# this is the directory from which you are running this R-program

# load the source code needed for the analyses
source(paste(common_source_folder,"RevBayes_Setup.r",sep="")); # Stuck? Use source(file.choose()) & grab it!

# these are the external databases that the program uses to fine-tune PaleoDB data
load(paste(data_for_R_folder,"Rock_Unit_Database.RData",sep=""));  		# data for rock-units including biozanation & superposition
load(paste(data_for_R_folder,"Gradstein_2012_Augmented.RData",sep="")); # refined Gradstein 2012 timescale & biozonations
load(paste(data_for_R_folder,"PaleoDB_Edits.RData",sep=""));			# edits to collections I cannot edit. # refined Gradstein 2012 timescale & biozonations

# load up the relevant databases!
rock_unit_databases <- rock_unit_data;						# information about rock ages, biozonations and sequence stratigraphy
chronostratigraphic_databases <- gradstein_2012_emended;	# information about time scales, zone ages, etc.
time_scale <- chronostratigraphic_databases$time_scale;		# time, time, time.....
paleodb_fixes <- paleodb_fixes;								# edits to PaleoDB data that I cannot enter currrently
analysis_name <- "Cincta";
onset <- "Cambrian";
end <- "Cambrian";
control_taxon <- "Echinodermata";							# this will be used for initial sampling estiates
zone_taxa <- c("Trilobita","Archaeocyatha");				# this will be used to tie down collection ages
basic_environments <- c("marine","unknown");
time_scale_stratigraphic_scale <- "Stage Slice";			# this will set the stratigraphic scale that we'll use
bogarted <- F;												# if true, then you wll be prompted for a file of additional occurrences
temporal_precision <- 0.1;									# level of precision for trying to estimate ages;
cinctans <- accio_PaleoDB_data_from_chosen_nexus_file(onset=onset,end=end,rock_unit_databases=rock_unit_databases,chronostratigraphic_databases=chronostratigraphic_databases,paleodb_fixes=paleodb_fixes,control_taxon=control_taxon,zone_taxa=zone_taxa,basic_environments=basic_environments,time_scale_stratigraphic_scale=time_scale_stratigraphic_scale,analysis_name=analysis_name,bogarted=bogarted,temporal_precision = temporal_precision);
write.csv(cinctans$occurrences,paste(local_directory,analysis_name,"_Occurrences_Final.csv",sep=""),row.names=F);
write.csv(cinctans$collections,paste(local_directory,analysis_name,"_Collections_Final.csv",sep=""),row.names=F);
write.csv(cinctans$fossil_summaries,paste(local_directory,analysis_name,"_Stratigraphic_Range_Summaries.csv",sep=""),row.names=F);
write.csv(cinctans$relv_time_scale,paste(local_directory,analysis_name,"_Time_Scale.csv",sep=""),row.names=F);
strat_ranges_fuzzy <- strat_range_summaries <- cinctans$fossil_summaries;
strat_ranges_fuzzy$total_finds <- strat_ranges_fuzzy[,ncol(strat_ranges_fuzzy)] <- NULL;

#### Get Sampling Metrics ####
finest_chronostrat <- cinctans$relv_time_scale;
finest_chronostrat$span <- finest_chronostrat$ma_lb-finest_chronostrat$ma_ub;
taxon_sites <- cinctans$collections;
slice_sites <- tally_collections_occupied_by_subinterval(taxon_sites,hierarchical_chronostrat = finest_chronostrat);
slice_rocks <- tally_rock_units_occupied_by_subinterval(taxon_sites,finest_chronostrat);

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

write.csv(slice_sites,paste(local_directory,analysis_name,"_Sites_per_Bin.csv",sep=""));
write.csv(slice_rocks,paste(local_directory,analysis_name,"_Rock_Units_per_Bin.csv",sep=""));

distributed_sampling_over_time_sites <- sampling_over_time_sites <- accio_sampling_distributions_for_RevBayes(finds_per_bin = sites_per_bin_spc_rnd,sample_units_per_bin = ceiling(slice_sites));
distributed_sampling_over_time_rocks <- sampling_over_time_rocks <- accio_sampling_distributions_for_RevBayes(finds_per_bin = rocks_per_bin_spc_rnd,sample_units_per_bin = ceiling(slice_rocks));
distributed_sampling_over_time_sites$uniform <- NULL;
distributed_sampling_over_time_rocks$uniform <- NULL;
all_intervals <- names(slice_rocks);
# per-site probabilities of finding taxa;
quantiles_for_sites <- accio_sampling_quantiles_for_all_intervals(sampling_over_time = distributed_sampling_over_time_sites,all_intervals = all_intervals,criterion="AICc",ttl_quantiles=21);
# per-rock_unit probabilities of finding taxa;
quantiles_for_rocks <- accio_sampling_quantiles_for_all_intervals(sampling_over_time = distributed_sampling_over_time_rocks,all_intervals = all_intervals,criterion="AICc",ttl_quantiles=21);

pfind_sites_quantiles <- pfind_rocks_quantiles <- exp_finds_sites_quantiles <- exp_finds_rocks_quantiles <- quantiles_for_rocks;
nbins <- nrow(quantiles_for_sites);
prob_missing_site <- prob_missing_rock <- vector(length=nbins);
for (i in 1:nbins)	{
	pfind_sites_quantiles[i,] <- simplify2array(1-((1-quantiles_for_sites[i,])^slice_sites[i]));
	prob_missing_site[i] <- mean(as.numeric(1-pfind_sites_quantiles[i,])^slice_rocks[i]);
	pfind_rocks_quantiles[i,] <- simplify2array(1-((1-quantiles_for_rocks[i,])^slice_rocks[i]));
	prob_missing_rock[i] <- mean(as.numeric(1-pfind_rocks_quantiles[i,])^slice_rocks[i]);
	exp_finds_sites_quantiles[i,] <- simplify2array(probability_to_Poisson_rate(pfind_sites_quantiles[i,]));
	exp_finds_rocks_quantiles[i,] <- simplify2array(probability_to_Poisson_rate(pfind_rocks_quantiles[i,]));
	}

study_bins <- finest_chronostrat$interval[min(strat_range_summaries$latest_poss_fa)<=finest_chronostrat$ma_ub];
study_bins <- study_bins[rowMeans(exp_finds_sites_quantiles[study_bins,])>0]
last_bin <- study_bins[length(study_bins)];

# estimate probability of sampling individual from the interval in which the youngest taxon appears!
rho <- (mean(as.numeric(pfind_rocks_quantiles[match(last_bin,rownames(pfind_rocks_quantiles)),]))+mean(as.numeric(pfind_sites_quantiles[match(last_bin,rownames(pfind_rocks_quantiles)),])))/2;
# get the average sampling rate per million years
phi <- mean(((rowMeans(exp_finds_sites_quantiles[study_bins,])+rowMeans(exp_finds_rocks_quantiles[study_bins,]))/2)/finest_chronostrat$span[match(study_bins,finest_chronostrat$interval)]);

#### Get Initial Estimates of Diversification ####
three_timer_setup <- setup_three_timer_analysis(samples_per_interval = sites_per_bin_spc_rnd);
taxon_ranges <- three_timer_setup$taxon_ranges;
sampled_in_bin <- three_timer_setup$sampled_in_bin;
sepkoski_richness <- three_timer_setup$sepkoski_richness;
gappers <- three_timer_setup$gappers;

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


#### Get information for RevBayes ####
basic_data <- accio_data_from_chosen_nexus_file();
otu_names_used <- basic_data$OTUs;
taxon_names <- otu_names_used[!tolower(otu_names_used) %in% c("outgroup","out")];
otu_names <- sapply(taxon_names,scourgify_taxon_names);
ingroup_taxa <- otu_names[!(1:length(otu_names)) %in% basic_data$Outgroup];
notu <- length(otu_names);
uncoded_taxa <- basic_data$Unscored_Taxa;
if (is.null(uncoded_taxa))
	uncoded_taxa <- "";

filename <- paste(local_directory,analysis_name,"_Diversification_Numbers_given_Rocks.csv",sep="");
diversification_ma_rock <- read.csv(filename,header=T,stringsAsFactors = F);
filename <- paste(local_directory,analysis_name,"_Diversification_Numbers_given_Sites.csv",sep="");
diversification_ma_site <- read.csv(filename,header=T,stringsAsFactors = F);
filename <- paste(local_directory,analysis_name,"_per_Interval_Sampling_given_Sites.csv",sep="");
site_sampling_distribution <- read.csv(filename,header = T,stringsAsFactors = F);
filename <- paste(local_directory,analysis_name,"_per_Interval_Sampling_given_Rocks.csv",sep="");
rock_sampling_distribution <- read.csv(filename,header = T,stringsAsFactors = F);
filename <- paste(local_directory,analysis_name,"_Three_Timer_Ingredients.csv",sep="");
three_timer_info <- read.csv(file=filename,header = T,stringsAsFactors = F);
rownames(three_timer_info) <- rownames(diversification_ma_rock) <- rownames(diversification_ma_site) <- rownames(rock_sampling_distribution) <- rownames(site_sampling_distribution) <- rock_sampling_distribution$X;
three_timer_info$X <- diversification_ma_rock$X <- diversification_ma_site$X <- site_sampling_distribution$X <- rock_sampling_distribution$X <- NULL;
nbins <- nrow(three_timer_info);

otu_fossil_summary <- read.csv(paste(local_directory,analysis_name,"_Stratigraphic_Range_Summaries.csv",sep=""),header=T,stringsAsFactors = F);
#otu_fossil_summary$earliest_poss_fa;
#otu_fossil_summary$latest_poss_la;
#otu_fossil_summary <- read.csv(file.choose(),sep=""),header=T,stringsAsFactors = F);
clade_time_scale <- read.csv(paste(local_directory,analysis_name,"_Time_Scale.csv",sep=""),header=T,stringsAsFactors = F);

fossil_intervals_for_RevBayes_raw <- data.frame(taxon=as.character(otu_fossil_summary$taxon),
												min=as.numeric(otu_fossil_summary$earliest_poss_la+otu_fossil_summary$latest_poss_la)/2,
												max=as.numeric(otu_fossil_summary$earliest_poss_fa+otu_fossil_summary$latest_poss_fa)/2,
												stringsAsFactors = F);
fossil_intervals_for_RevBayes_raw$taxon <- gsub(" ","_",fossil_intervals_for_RevBayes_raw$taxon);

fossil_intervals_for_RevBayes <- fossil_intervals_for_RevBayes_raw;
fossil_intervals_for_RevBayes$max <- fossil_intervals_for_RevBayes$max-min(fossil_intervals_for_RevBayes$min);
fossil_intervals_for_RevBayes$min <- fossil_intervals_for_RevBayes$min-min(fossil_intervals_for_RevBayes$min);

fossil_intervals_for_RevBayes_FA_only <- fossil_intervals_for_RevBayes;
fossil_intervals_for_RevBayes_FA_only$max <- fossil_intervals_for_RevBayes_FA_only$max-min(fossil_intervals_for_RevBayes_FA_only$max);
fossil_intervals_for_RevBayes_FA_only$min <- fossil_intervals_for_RevBayes_FA_only$max;

fossil_interval_file <- paste(tolower(analysis_name),"_fossil_intervals.tsv",sep="");
fossil_interval_file_FA <- paste(tolower(analysis_name),"_fossil_intervals_FA.tsv",sep="");
fossil_interval_file_raw <- paste(tolower(analysis_name),"_fossil_intervals_orig_dates.tsv",sep="");

write.table(fossil_intervals_for_RevBayes,file=paste(local_directory,fossil_interval_file,sep=""),row.names = F,sep="\t",quote = F);
write.table(fossil_intervals_for_RevBayes_FA_only,file=paste(local_directory,fossil_interval_file_FA,sep=""),row.names = F,sep="\t",quote = F);
write.table(fossil_intervals_for_RevBayes_raw,file=paste(local_directory,fossil_interval_file_raw,sep=""),row.names = F,sep="\t",quote = F);

#clade_time_scale <- read.csv(file.choose(),sep=""),header=T,stringsAsFactors = F);

#### Now, write the FBD Script ####
obins <- (1:nbins)[diversification_ma_rock$ML_Orig>0][(1:nbins)[diversification_ma_rock$ML_Orig>0] %in% (1:nbins)[three_timer_info$sampled_in_bin>0]]
#origination <- sum(three_timer_info$sampled_in_bin[obins]*diversification_ma_rock$ML_Orig[obins])/sum(three_timer_info$sampled_in_bin[obins]);
ebins <- (1:nbins)[diversification_ma_rock$ML_Extn>0][(1:nbins)[diversification_ma_rock$ML_Extn>0] %in% (1:nbins)[three_timer_info$sampled_in_bin>0]]
#extinction <- sum(three_timer_info$sampled_in_bin[ebins]*diversification_ma_rock$ML_Extn[ebins])/sum(three_timer_info$sampled_in_bin[ebins]);
origination <- median((diversification_ma_rock$ML_Orig[obins]+diversification_ma_site$ML_Orig[obins])/2);
extinction <- median((diversification_ma_rock$ML_Extn[ebins]+diversification_ma_site$ML_Extn[ebins])/2);
turnover <- extinction/origination;
psi <- sampling <- (sum((1-rock_sampling_distribution$prob_missing)*three_timer_info$sampled_in_bin)/sum(three_timer_info$sampled_in_bin)+sum((1-site_sampling_distribution$prob_missing)*three_timer_info$sampled_in_bin)/sum(three_timer_info$sampled_in_bin))/2;
phi <- prob_sampling_clade_bapst(p=origination,q=extinction,r=psi);
rho <- 1-(rock_sampling_distribution$prob_missing[sum(clade_time_scale$ma_lb>=min(fossil_intervals_for_RevBayes_raw$max))]+site_sampling_distribution$prob_missing[sum(clade_time_scale$ma_lb>=min(fossil_intervals_for_RevBayes_raw$max))])/2;
control_taxon <- "Trilobita";
sampling_unit <- "rocks and sites averaged";
#divergence_bounds <- c(max(fossil_intervals_for_RevBayes_FA_only$max),clade_time_scale$ma_lb[1]-min(fossil_intervals_for_RevBayes_raw$max));
divergence_bounds <- ceiling(max(fossil_intervals_for_RevBayes_FA_only$max))+c(0,ceiling(20/(psi+(origination*phi))));
extant_taxa <- fossil_intervals_for_RevBayes_FA_only$taxon[fossil_intervals_for_RevBayes_FA_only$max==0];
fbd_script <- scribio_fbd_portion_of_Rev_Bayes_script(analysis_name,write_scripts_directory="",origination,extinction,psi,rho,divergence_bounds,control_taxon,extant_taxa,otu_names,uncoded_taxa,script_file_lead="scripts/");
#   analysis_name,write_scripts_directory,origination,extinction,psi,rho,divergence_bounds,control_taxon,extant_taxa,otu_names,uncoded_taxa="",script_file_lead="scripts/"
fbd_file_name <- fbd_script$filename;

fa_info <- fossil_intervals_for_RevBayes_FA_only;
fa_info$min <- NULL;
colnames(fa_info)[colnames(fa_info)=="max"] <- "fa";
scribio_RevBayes_scripts_from_chosen_nexus_file_and_existing_FBD_script_and_data(analysis_name,local_directory=local_directory,fa_info=fa_info);
