#### SETUP ####
# accersi: fetch/summon
# divido: divide!
# expello: banish
# mundus: clean
# percursant: scour
# revelare: reveal

### udpated with routines written for first phylogenetics course at UNL
MAXNO <- 1.797693e+308;
newick_verbotten <- c(".","?","\"","\'");
letter_states <- LETTERS[!LETTERS %in% c("I","O")];
zzzz <- 0.25;

#### HOUSE-CLEANING ####
clear_matrix_na_with_another_cell_value <- function(data,j, k)	{
size <- dim(data)
for (i in 1:size[1])	{
	if (is.na(data[i,j]))	data[i,j] <- data[i,k]
	}
return(data)
}

count_characters_in_string <- function(string_to_count)	{
j <- strsplit(string_to_count,split="",fixed=TRUE)[[1]];
return(length(j));
}

mundify_nexus_text <- function(nexus_line)	{
nexus_line <- gsub("\xd4","",nexus_line);
nexus_line <- gsub("\xd5","",nexus_line);
nexus_line <- gsub("\x87","a",nexus_line);
nexus_line <- gsub("\xfc\xbe\x8d\xa3\xa4\xbc","n",nexus_line);
nexus_line <- gsub("\xfc\xbe\x98\x93\xa0\xbc","ae",nexus_line);
nexus_line <- gsub("\xfc\xbe\x99\x83\xa0\xbc","c",nexus_line);
nexus_line <- gsub("\xfc\xbe\x98\x96\x8c\xbc","",nexus_line);
nexus_line <- gsub("\xfc\xbe\x8c\x93\xa4\xbc","\'",nexus_line);
nexus_line <- gsub("\'\'\'","\'\"",nexus_line);
nexus_line <- gsub("\'\' ","\" ",nexus_line);
nexus_line <- gsub("\xac","",nexus_line);
nexus_line <- gsub("\xa0","†",nexus_line);
nexus_line <- gsub("\x80","",nexus_line);
nexus_line <- gsub("\xd1","",nexus_line);
nexus_line <- gsub("\xc9","?",nexus_line);
nexus_line <- gsub("\xe1","a",nexus_line);
nexus_line <- gsub("\xe9","e",nexus_line);
nexus_line <- gsub("\x8e","e",nexus_line);
nexus_line <- gsub("\x8f","e",nexus_line);
nexus_line <- gsub("\x92","i",nexus_line);
nexus_line <- gsub("\xbf","o",nexus_line);
nexus_line <- gsub("\x9a","o",nexus_line);
nexus_line <- gsub("\x97","o",nexus_line);
nexus_line <- gsub("\xf6","ö",nexus_line);
nexus_line <- gsub("\xfc","ue",nexus_line);
nexus_line <- gsub("\xb0","˚",nexus_line);
nexus_line <- gsub("\xba","˚",nexus_line);
nexus_line <- gsub("\xfc\xbe\x8e\x93\xa4\xbc","o",nexus_line);
nexus_line <- gsub("\x9f","ue",nexus_line);
nexus_line <- gsub("\xd0","-",nexus_line);
nexus_line <- gsub("\xd2","\"",nexus_line);
nexus_line <- gsub("\xd3","\"",nexus_line);
nexus_line <- gsub("\xfc\xbe\x8d\x86\x90\xbc","\'",nexus_line);
nexus_line <- gsub("\xfc\xbe\x8d\x86\x8c\xbc","ƒ",nexus_line);
nexus_line <- gsub("\xdf","ß",nexus_line);
nexus_line <- gsub("\xa7","ß",nexus_line);
nexus_line <- gsub("\xfc\xbe\x8c\xa6\x88\xbc","≤",nexus_line);
nexus_line <- gsub("\xb3","≥",nexus_line);
nexus_line <- gsub("\xfc\xbe\x8d\x96\x8c\xbc","≈",nexus_line);
nexus_line <- gsub("\xfc\xbe\x98\xa6\x98\xbc","˚",nexus_line);
nexus_line <- gsub("\xb6","∂",nexus_line);
nexus_line <- gsub("\xc6","∆",nexus_line);
nexus_line <- gsub("\xfc\xbe\x8d\xb6\x88\xbc","∑",nexus_line);
nexus_line <- gsub("\xfc\xbe\x99\x86\x88\xbc","Ω",nexus_line);
nexus_line <- gsub("\xa5"," ",nexus_line);
nexus_line <- gsub("Á","A",nexus_line);
nexus_line <- gsub("Ä","A",nexus_line);
nexus_line <- gsub("ä","a",nexus_line);
nexus_line <- gsub("á","a",nexus_line);
nexus_line <- gsub("å","a",nexus_line);
nexus_line <- gsub("Ç","C",nexus_line);
nexus_line <- gsub("ç","c",nexus_line);
nexus_line <- gsub("č","c",nexus_line);
nexus_line <- gsub("é","e",nexus_line);
nexus_line <- gsub("è","e",nexus_line);
nexus_line <- gsub("ê","e",nexus_line);
nexus_line <- gsub("ė","e",nexus_line);
nexus_line <- gsub("î","i",nexus_line);
nexus_line <- gsub("Î","I",nexus_line);
nexus_line <- gsub("ñ","n",nexus_line);
nexus_line <- gsub("Ö","O",nexus_line);
nexus_line <- gsub("Ø","O",nexus_line);
nexus_line <- gsub("ø","o",nexus_line);
nexus_line <- gsub("ó","o",nexus_line);
nexus_line <- gsub("ö","o",nexus_line);
nexus_line <- gsub("õ","o",nexus_line);
nexus_line <- gsub("Š","S",nexus_line);
nexus_line <- gsub("š","s",nexus_line);
nexus_line <- gsub("ů","u",nexus_line);
nexus_line <- gsub("ü","u",nexus_line);
nexus_line <- gsub("’","’",nexus_line);
nexus_line <- gsub("\x88","a",nexus_line);
nexus_line <- gsub("Ã„","A",nexus_line);
nexus_line <- gsub("Á","A",nexus_line);
nexus_line <- gsub("Ã¡","a",nexus_line);
nexus_line <- gsub("Ã¤","a",nexus_line);
nexus_line <- gsub("Ã¥","a",nexus_line);
nexus_line <- gsub("Ã§","c",nexus_line);
nexus_line <- gsub("Ã©","e",nexus_line);
nexus_line <- gsub("Ã¨","e",nexus_line);
nexus_line <- gsub("Ã±","n",nexus_line);
nexus_line <- gsub("Ã–","O",nexus_line);
nexus_line <- gsub("Ã¸","o",nexus_line);
nexus_line <- gsub("Ã¶","o",nexus_line);
nexus_line <- gsub("Ãµ","o",nexus_line);
nexus_line <- gsub("Ã´","o",nexus_line);
nexus_line <- gsub("Ã¼","u",nexus_line);
nexus_line <- gsub("√Æ","i",nexus_line);
nexus_line <- gsub("≈†","S",nexus_line);
nexus_line <- gsub("≈°","s",nexus_line);
nexus_line <- gsub("√•","a",nexus_line);
nexus_line <- gsub("&#367;","u",nexus_line);
nexus_line <- gsub("&#945;","α",nexus_line);
return(nexus_line);
}

# Turn Lophospira sp. or Lophospira sp. A to Lophospira
reduce_genus_sp_to_genus <- function(taxon_name)	{
taxon_names <- strsplit(taxon_name," ")[[1]];
indet_species <- c("sp.");
indet_species <- c(indet_species,paste("sp.",LETTERS));
indet_species <- c(indet_species,paste("sp.",letters));
indet_species <- c(indet_species,"spp.");
indet_species <- c(indet_species,paste("spp.",LETTERS));
indet_species <- c(indet_species,paste("spp.",letters));
for (i in 1:100)	indet_species <- c(indet_species,paste("sp.",i));
for (i in 1:100)	indet_species <- c(indet_species,paste("nov.",i));
for (i in 1:100)	indet_species <- c(indet_species,paste("sp. nov.",i));
for (i in 1:100)	indet_species <- c(indet_species,paste("indet.",i));

if (sum(taxon_names %in% indet_species)>0)
	taxon_names <- taxon_names[(1:length(taxon_names))[!taxon_names %in% indet_species]];
taxon_name <- paste(taxon_names,collapse=" ");
return(taxon_name);
}

accersi_study_name <- function(nexus_file_name)	{
n_f_n <- simplify2array(strsplit(nexus_file_name,"/"));
n_f_n <- n_f_n[length(n_f_n)];
n_f_n <- gsub("\\.nex","",n_f_n);
filename_parts <- simplify2array(strsplit(n_f_n,"_"));

if (sum(filename_parts %in% "&")==1)	{
	author_start <- match("&",filename_parts)-1
	} else if (sum(filename_parts %in% c("et","al"))==2)	{
	author_start <- match("et",filename_parts)-1;
	} else	{
	author_start <- length(filename_parts)-2;
	}
return(paste(filename_parts[1:(author_start-1)]," (",paste(filename_parts[author_start:length(filename_parts)],collapse=" "),")",sep=""));
}


#### ROUTINES TO READ CHARACTER MATRIX INFORMATION IN NEXUS FILE ####
# routine to read nexus file of Mesquite or Maclade format & return important infromation
accersi_data_from_nexus_file <- function(nexus_file_name, polymorphs=T, UNKNOWN=-11, INAP=-22, rate_partitions="", trend_partitions="")	{
# nexus_file_name: name of nexus file (e.g., "Phacopidae_Snow_2000.nex")
# polymorphs: boolean, if TRUE, then recode "1,2" as "-21"; otherwise, treat as unknown
# UNKNOWN: value substituting for "?"
# INAP: value substituting for gap ("-")
# rate_partitions: nameof CHARPARTITION that you want to use for dividing characters into general rate classes.
nexus <- scan(file=nexus_file_name,what=character(),sep="\n");
output <- accersi_data_from_nexus_vector(nexus=nexus,polymorphs=polymorphs,UNKNOWN=UNKNOWN,INAP=INAP,rate_partitions=rate_partitions,trend_partitions=trend_partitions);
return(output)
}

# routine to prompt user for a nexus file of Mesquite or Maclade format & return important infromation
accersi_data_from_chosen_nexus_file <- function(polymorphs=T, UNKNOWN=-11, INAP=-22, rate_partitions="", trend_partitions="")	{
# nexus_file_name: name of nexus file (e.g., "Phacopidae_Snow_2000.nex")
# polymorphs: boolean, if TRUE, then recode "1,2" as "-21"; otherwise, treat as unknown
# UNKNOWN: value substituting for "?"
# INAP: value substituting for gap ("-")
# rate_partitions: nameof CHARPARTITION that you want to use for dividing characters into general rate classes.
#print("Choose the nexus file you with to analyze: ");
print("Choose the nexus file that you wish to analyze: ");
flush.console();
Sys.sleep(zzzz);
nexus_file_name <- file.choose();
nexus <- scan(file=nexus_file_name,what=character(),sep="\n");
output <- accersi_data_from_nexus_vector(nexus=nexus,polymorphs=polymorphs,UNKNOWN=UNKNOWN,INAP=INAP,rate_partitions=rate_partitions,trend_partitions=trend_partitions);
return(output)
}

# routine to read nexus information Mesquite or Maclade format & return important infromation
accersi_data_from_nexus_vector <- function(nexus, polymorphs=T, UNKNOWN=-11, INAP=-22, rate_partitions="", trend_partitions="")	{
# nexus_file_name: name of nexus file (e.g., "Phacopidae_Snow_2000.nex")
# polymorphs: boolean, if TRUE, then recode "1,2" as "-21"; otherwise, treat as unknown
# UNKNOWN: value substituting for "?"
# INAP: value substituting for gap ("-")
# rate_partitions: nameof CHARPARTITION that you want to use for dividing characters into general rate classes.
ml <- 0;
#i <- 1

for (i in 1:length(nexus))  {
	nexus[i] <- mundify_nexus_text(nexus_line = nexus[i]);
#	j <- simplify2array(strsplit(nexus[i],split="",fixed=TRUE))[,1];
	j <- strsplit(nexus[i],split="",fixed=TRUE)[[1]];
#	str_split(string=nexus[i],pattern="")
	if (length(j)>ml) ml <- length(j);
	}
ml <- ml+1;	# LENGTH OF LONGEST LINE
	
# file is now a vector of characters.  Turn it into a matrix with one char per cell
nexusfile <- matrix("\n",length(nexus),ml)
for (i in 1:length(nexus))  {
	j <- strsplit(nexus[i],split="",fixed=TRUE)[[1]];
	for (k in 1:length(j))      nexusfile[i,k] <- j[k]
	if ((length(j)+2)<ml)
		for (k in (length(j)+2):ml) nexusfile[i,k] <- ""
	}

top <- match("matrix",tolower(nexus));
if (is.na(top)) top <- match("\tmatrix",tolower(nexus));
if (is.na(top))	{
	top <- 0;
	ln <- 1;		# this is the row with the word "Matrix": character data starts next.
	while (top==0)	{
		em_nexus <- gsub("\t","",nexus[ln]);
		nexus_words <- simplify2array(strsplit(em_nexus," ")[[1]]);
		if (!is.na(match("matrix",tolower(nexus_words))))	{
			top <- ln;
			} else	ln <- ln+1;
		}
	}
top <- top+1;	# this will give the first row of data
# skip the comment text denoting character numbers (if present)
while(nexusfile[top,1]=="[" || nexusfile[top,1]==" ") top <- top+1;

all_states <- c();
missing <- "?";
gap <- "-";
notu <- nchars <- strat <- range <- geog <- 0;
for (i in 2:top)  {
	while ((nexusfile[i,1]=="[" || nexusfile[i,1]=="\n") && i<top)	i <- i+1;
	em_nexus <- gsub("\t","",nexus[i]);
	em_nexus <- gsub("="," = ",em_nexus);
	em_nexus <- gsub(";"," ; ",em_nexus);
	nexus_words <- simplify2array(strsplit(em_nexus," ")[[1]]);
	nexus_words <- nexus_words[nexus_words!=""];
	if (!is.na(match("ntax",tolower(nexus_words))) || !is.na(match("ntaxa",tolower(nexus_words))))	{
		j <- 1+match("ntax",tolower(nexus_words));
		if (is.na(j))	j <- 1+match("ntaxa",tolower(nexus_words));
		while(nexus_words[j]=="=")	j <- j+1;
		notu <- as.numeric(nexus_words[j]);
		}
	if (!is.na(match("nchar",tolower(nexus_words))) || !is.na(match("nchars",tolower(nexus_words))))	{
		j <- 1+match("nchar",tolower(nexus_words));
		if (is.na(j))	j <- 1+match("nchars",tolower(nexus_words));
		while(nexus_words[j]=="=")	j <- j+1;
		nchars <- as.numeric(nexus_words[j]);
		}
	if (!is.na(match("gap",tolower(nexus_words))))	{
		if (nexus_words[match("gap",tolower(nexus_words))+1]=="=")	{
			j <- 1+match("gap",tolower(nexus_words));
			while(nexus_words[j]=="=")	j <- j+1;
			gap <- nexus_words[j];
			}
		}
	if (!is.na(match("missing",tolower(nexus_words))))	{
		if (nexus_words[match("missing",tolower(nexus_words))+1]=="=")	{
			j <- 1+match("missing",tolower(nexus_words));
			while(nexus_words[j]=="=")	j <- j+1;
			missing <- nexus_words[j];
			}
		}
	if (!is.na(match("symbols",tolower(nexus_words))))	{
		j <- match("symbols",tolower(nexus_words))+1;
		while (nexus_words[j] %in% c("=","\""))	j <- j+1;
		jj <- min(((j+1):length(nexus_words))[tolower(nexus_words)[((j+1):length(nexus_words))] %in% c("missing","gap",";")]-1);
#		jj <- j+match(";",nexus_words[(j+1):length(nexus_words)])-1;
		all_states <- gsub("\"","",nexus_words[j:jj]);
#		if (tolower(all_states) %in% "")
		}
	if (!is.na(match("fa",tolower(nexus_words))) || !is.na(match("fka",tolower(nexus_words))))	{
		nexus_words[tolower(nexus_words)=="fka"] <- "fa";
		strat <- as.numeric(nexus_words[match("fa",tolower(nexus_words))-1]);
		}
	if (!is.na(match("la",tolower(nexus_words))) || !is.na(match("lka",tolower(nexus_words))))	{
		nexus_words[tolower(nexus_words)=="lka"] <- "la";
		range <- as.numeric(nexus_words[match("la",tolower(nexus_words))-1]);
		}
	if (!is.na(match("geog",tolower(nexus_words))))	{
		geog <- as.numeric(nexus_words[match("geog",tolower(nexus_words))-1]);
		} else if (!is.na(match("geography",tolower(nexus_words))))	{
		geog <- as.numeric(nexus_words[match("geography",tolower(nexus_words))-1]);
		}
	}
if (is.null(all_states))
	all_states <- 0:9;

extra <- 0;
if (strat>0)	{
	if (range>0)	{
		nchars <- nchars-2
		extra <- 2
		} else {
		nchars <- nchars-1
		extra <- 1
		} 
	strat_ranges <- matrix(0,notu,2)	
	}
if (geog>0)	{
	nchars <- nchars-1
	geography <- vector(length=notu)
	extra <- extra+1
	}
	
taxa <- vector(length=notu);
nstates <- array(0,dim=nchars);
chmatrix <- matrix(0,notu,nchars);
tx <- 1;

# look for outgroup designation
exclude <- outgroup <- -1;
if (!is.na(match("BEGIN SETS;",nexus)))	{
	tx_pt <- match("BEGIN SETS;",nexus);	# look at taxon partitions
	look_for_outgroup <- TRUE;
	while (look_for_outgroup)	{
		tx_pt <- 1+tx_pt;
		yyy <- paste(nexusfile[tx_pt,], collapse = "");
		yyy <- gsub("-"," - ",yyy);
		yyy <- gsub("- "," - ",yyy);
		yyy <- gsub("  -  "," - ",yyy);
		yyy <- gsub(";","",yyy);
		yyy <- gsub(","," ,",yyy);
		yyy <- gsub("\n","",yyy);
		yyy <- gsub("\r","",yyy);
		yyy <- gsub("\t","",yyy);
		xxx <- tolower(strsplit(yyy," ")[[1]]);
		xxx <- xxx[xxx!=""];
		if (!is.na(match("outgroup",tolower(xxx))))	{
			ttl_ln <- length(xxx);
			jj <- 1+match("outgroup",tolower(xxx));
			while (xxx[jj]==":" || xxx[jj]=="=")	jj <- jj+1;
			outgroup <- c();
			while (xxx[jj]!="," && jj<=ttl_ln)	{
				if (xxx[jj]=="-")	{
					jj <- jj+1;
					outgroup <- c(outgroup,((as.numeric(outgroup[length(outgroup)])+1):as.numeric(xxx[jj])));
					} else	{
					outgroup <- c(outgroup,xxx[jj]);
					}
				jj <- jj+1;
				}
			look_for_outgroup <- FALSE;
			} else	{
			if (tolower(nexus[tx_pt])=="end;" || tolower(nexus[tx_pt])=="\tend;")
				look_for_outgroup <- FALSE;
			}
		}

	# look for characters to exclude
	tx_pt <- match("BEGIN SETS;",nexus);
	xxx <- strsplit(paste(nexusfile[tx_pt-1,],collapse = "")," ");
	while(tolower(xxx[1])!="end")	{
		tx_pt <- tx_pt+1;
		yyy <- paste(nexusfile[tx_pt,], collapse = "");
		yyy <- gsub("- "," - ",yyy);
		yyy <- gsub(";","",yyy);
		yyy <- gsub(","," ,",yyy);
		yyy <- gsub("\n","",yyy);
		yyy <- gsub("\r","",yyy);
		yyy <- gsub("\t","",yyy);
		xxx <- tolower(strsplit(yyy," ")[[1]]);
		xxx <- xxx[xxx!=""];
		if (length(xxx)==0 || is.na(xxx))
			xxx <- "";
#		if (!is.na(xxx) && !is.null(xxx) && xxx!="")	{
		if (xxx[1]=="charpartition")	{
			if (xxx[1]=="charpartition" && !is.na(match("exclude",tolower(xxx))))	{
				ttl_ln <- length(xxx);
				jj <- 1+match("exclude",tolower(xxx));
				while (xxx[jj]==":")	jj <- jj+1;
				exclude <- c();
				while (xxx[jj]!="," && jj<ttl_ln)	{
					if (xxx[jj]=="-")	{
						jj <- jj+1;
						exclude <- c(exclude,((as.numeric(exclude[length(exclude)])+1):as.numeric(xxx[jj])));
						} else	{
						exclude <- c(exclude,as.numeric(xxx[jj]));
						}
					jj <- jj+1;
					}
				}
			}
#		xxx[1];
#		tx_pt;
		}
	}

if (rate_partitions!="")	{
	ln <- match("BEGIN SETS;",nexus);
	got_splits <- F;
	while (!got_splits)	{
		ln <- ln+1;
		breakup_this_line <- strsplit(nexus[ln],split=" ")[[1]];
		if (!is.na(match(rate_partitions,breakup_this_line)))	{
			nexus[ln] <- gsub("-"," - ",nexus[ln]);	# Mesquite often puts dashes immediately after character or taxon numbers.....
			nexus[ln] <- gsub("  -"," -",nexus[ln]);
			nexus[ln] <- gsub("-  ","- ",nexus[ln]);
			breakup_this_line <- strsplit(nexus[ln],split=" ")[[1]];
			breakup_this_line <- gsub(",","",breakup_this_line);
			breakup_this_line <- gsub(";","",breakup_this_line);
			breakup_this_line <- breakup_this_line[breakup_this_line!=""];
			breakup_this_line <- breakup_this_line[match(rate_partitions,breakup_this_line):length(breakup_this_line)];
			kk <- (1:length(breakup_this_line))[breakup_this_line %in% ":"];
			partition_names <- breakup_this_line[kk-1];
			kk <- c(kk,length(breakup_this_line)+1);	# add last numberso that we can end the partion search easily below
			character_rate_partitions <- rep("",nchars);
			for (pn in 1:length(partition_names))	{
				ll <- kk[pn]+1;
				this_part <- as.numeric(breakup_this_line[ll]);
				ll <- ll+1;
#				while (ll<(kk[pn+1]-1))	{
				if (pn < length(partition_names))	{
					break_cell <- kk[pn+1]-1;
					} else	{
					break_cell <- kk[pn+1];
					}
				while (ll<break_cell)	{
					if (breakup_this_line[ll]=="-")	{
						ll <- ll+1;
						this_part <- c(this_part,as.numeric(breakup_this_line[ll-2]:as.numeric(breakup_this_line[ll])));
						} else	{
						this_part <- c(this_part,as.numeric(breakup_this_line[ll]));
						}
					ll <- ll+1;
					}
				character_rate_partitions[this_part] <- partition_names[pn];
				}
			got_splits<- T;
			}
		}
	} else	character_rate_partitions <- rep("imagine",nchars);

if (trend_partitions!="")	{
	ln <- match("BEGIN SETS;",nexus);
	got_splits <- F;
	while (!got_splits)	{
		ln <- ln+1;
		breakup_this_line <- strsplit(nexus[ln],split=" ")[[1]];
		if (!is.na(match(trend_partitions,breakup_this_line)))	{
			nexus[ln] <- gsub("-"," - ",nexus[ln]);	# Mesquite often puts dashes immediately after character or taxon numbers.....
			nexus[ln] <- gsub("  -"," -",nexus[ln]);
			nexus[ln] <- gsub("-  ","- ",nexus[ln]);
			breakup_this_line <- strsplit(nexus[ln],split=" ")[[1]];
			breakup_this_line <- gsub(",","",breakup_this_line);
			breakup_this_line <- gsub(";","",breakup_this_line);
			breakup_this_line <- breakup_this_line[breakup_this_line!=""];
			breakup_this_line <- breakup_this_line[match(trend_partitions,breakup_this_line):length(breakup_this_line)];
			kk <- (1:length(breakup_this_line))[breakup_this_line %in% ":"];
			partition_names <- breakup_this_line[kk-1];
			kk <- c(kk,length(breakup_this_line)+1);	# add last numberso that we can end the partion search easily below
			character_trend_partitions <- rep("",nchars);
			for (pn in 1:length(partition_names))	{
				ll <- kk[pn]+1;
				this_part <- as.numeric(breakup_this_line[ll]);
				ll <- ll+1;
#				while (ll<(kk[pn+1]-1))	{
				if (pn < length(partition_names))	{
					break_cell <- kk[pn+1]-1;
					} else	{
					break_cell <- kk[pn+1];
					}
				while (ll<break_cell)	{
					if (breakup_this_line[ll]=="-")	{
						ll <- ll+1;
						this_part <- c(this_part,as.numeric(breakup_this_line[ll-2]:as.numeric(breakup_this_line[ll])));
						} else	{
						this_part <- c(this_part,as.numeric(breakup_this_line[ll]));
						}
					ll <- ll+1;
					}
				character_trend_partitions[this_part] <- partition_names[pn];
				}
			got_splits<- T;
			}
		}
	} else	character_trend_partitions <- rep("square",nchars);

state_orders <- rep("unordered",nchars);

if (!is.na(match("BEGIN ASSUMPTIONS;",nexus)))	{
	tx_pt <- 1+match("BEGIN ASSUMPTIONS;",nexus);	# look at taxon partitions
	while (tolower(nexus[tx_pt])!="end;")	{
#		yyy <- paste(nexusfile[tx_pt,], collapse = "");
		yyy <- gsub("- "," - ",nexus[tx_pt]);
		yyy <- gsub(";","",yyy);
		yyy <- gsub(","," ,",yyy);
		yyy <- gsub("\n","",yyy);
		yyy <- gsub("\r","",yyy);
		yyy <- gsub("\t","",yyy);
		xxx <- tolower(strsplit(yyy," ")[[1]]);
		xxx <- xxx[xxx!=""];
		if (!is.na(match("ord:",tolower(xxx))) && !is.na(match("revbayes",tolower(xxx))))	{
			ttl_ln <- length(xxx);
			jj <- 1+match("ord:",xxx);
			while (xxx[jj]==":")	jj <- jj+1;
			ordered <- c();
			while (xxx[jj]!="," && jj<=ttl_ln)	{
				if (xxx[jj]=="-")	{
					jj <- jj+1;
					ordered <- c(ordered,((as.numeric(ordered[length(ordered)])+1):as.numeric(xxx[jj])));
					} else	{
					ordered <- c(ordered,as.numeric(xxx[jj]));
					}
				jj <- jj+1;
				}
			state_orders[ordered] <- "ordered";
			}
		tx_pt <- 1+tx_pt;
		}
	}
mxln <- length(nexusfile[top,]);
s <- top;
# te all of the taxon names
for (tx in 1:notu)	{
	# first, read taxon name
	#### look for quotations###
	s <- top+tx-1;
	endline <- match("\n",nexusfile[s,]);
	if (is.na(endline))	endline <- length(nexusfile[s,]);
	if (nexusfile[s,1]=="'" || nexusfile[s,2]=="'")	{
		jj <- ((1:length(nexusfile[s,]))[nexusfile[s,] %in% "'"]);
		i <- max((1:length(nexusfile[s,]))[nexusfile[s,] %in% "'"])
		taxa[tx] <- pracma::strcat(nexusfile[s,(jj[1]+1):(jj[2]-1)])
		i <- i+1
		while (nexusfile[s,i]==" " && i<ncol(nexusfile))	i <- i+1
		}	else	{
		i <- 1;
		if (nexusfile[s,1]!="\"")  {
			while (nexusfile[s,i]=="\t")	i <- i+1;
			taxa[tx] <- nexusfile[s,i]
			i <- i+1
			while (nexusfile[s,i]!=" " && nexusfile[s,i]!='\t' && i<ncol(nexusfile))	{
				if (nexusfile[s,i]!="_")	{
					taxa[tx] <- paste0(taxa[tx],as.character(nexusfile[s,i]))
					} else {
					taxa[tx] <- paste0(taxa[tx]," ")
					}
				i <- i+1
				}
			}	else {
			taxa[tx] <- nexusfile[s,2];
			i <- 3;
#			while (nexusfile[s,i]!=" " && nexusfile[s,i+1]!=" " && i<ncol(nexusfile))	{
			while (nexusfile[s,i]!=" " && i<ncol(nexusfile))	{
				if (as.character(nexusfile[s,i])!="\"")
					taxa[tx] <- paste0(taxa[tx],as.character(nexusfile[s,i]))
				i <- i+1;
				#print(taxa[tx]);
				}
			}
		# now, get to characters
		i <- (i:endline)[!nexusfile[s,i:endline] %in% c(" ","\t")][1];
#		while ((nexusfile[s,i]==" " || nexusfile[s,i]=="\t") && i<ncol(nexusfile))
#			i <- i+1
		}
	k <- i;
	if ((endline-k)==(nchars+extra))	{
		# true if there are no polymorphic characters for the taxon
		dummy <- nexusfile[s,k:(endline-1)];
		dummy[dummy==missing] <- UNKNOWN;
		dummy[dummy==gap] <- INAP;
		letterstate <- dummy[!dummy %in% c(UNKNOWN,INAP)];
		dummy[!dummy %in% c(UNKNOWN,INAP)] <- sapply(letterstate,switch_letter_state_to_numeric,all_states);
		chmatrix[tx,] <- as.numeric(dummy[1:nchars]);
		if (strat>0)	{
			strat_ranges[tx,1] <- strat_ranges[tx,2] <- as.numeric(dummy[strat])
			if (range>0)	strat_ranges[tx,2] <- as.numeric(dummy[range])
			}
		if (geog>0)	geography[tx]=as.numeric(nexusfile[geog,i])
		for (c in 1:nchars)	{
			if ((chmatrix[tx,c]+1)>nstates[c]) nstates[c] <- chmatrix[tx,c]+1
			}
		} else	{
#		for (c in 1:(nchars+extra))	{
		c <- 0;
		while (c < (nchars+extra))	{
			c <- c+1;
			#print(c);
			if (c<=nchars)	{
				if (nexusfile[s,i]=="(" || nexusfile[s,i]=="{")	{
					if (polymorphs==TRUE || polymorphs==1)	{
						# added 2020-11-28: sometimes polymorphics come in out-of-order
						riteparens <- (i:endline)[nexusfile[s,i:endline] %in% c(")","}")];
						ddd <- (i+1):(riteparens[1]-1);
						polysites <- ddd[!nexusfile[s,ddd] %in% c(",","&"," ")]
						polystates <- nexusfile[s,polysites];
						for (ps in 1:length(polystates))
							if (!polystates[ps] %in% 0:9)
								polystates[ps] <- switch_letter_state_to_numeric(polystates[ps],all_states=all_states);
						nexusfile[s,polysites] <- sort(as.numeric(polystates));
#						nexusfile[s,(i+1):(riteparens[1]-1)] <- sort(nexusfile[s,(i+1):(riteparens[1]-1)]);
						i <- i+1;
						w <- as.numeric(nexusfile[s,i])
						chmatrix[tx,c] <- -1*as.numeric(nexusfile[s,i])
						if ((1+w)>nstates[c])  nstates[c] <- 1+w;
						i <- i+1
						j <- 1
						while (nexusfile[s,i]!=")" && nexusfile[s,i]!="}" && i<ncol(nexusfile))	{
							if (nexusfile[s,i]!="," && nexusfile[s,i]!=" ")	{
								w <- as.numeric(nexusfile[s,i])
								if ((w+1)>nstates[c])	nstates[c] <- w+1
								chmatrix[tx,c] <- chmatrix[tx,c]-((10^j)*w)
								i <- i+1
								j <- j+1
								} else {
								i <- i+1
								}
							}
						}	else {
						chmatrix[tx,c] <- UNKNOWN;
						while (nexusfile[s,i]!=')' && nexusfile[s,i]!="}")	i <- i+1;
						}
					} else if (nexusfile[s,i]==missing)	{
					chmatrix[tx,c] <- UNKNOWN;
					}	else if (nexusfile[s,i]==gap)	{
					chmatrix[tx,c] <- INAP;
					} else if (nexusfile[s,i]>="A" && nexusfile[s,i]<="Z")  {
					chmatrix[tx,c] <- switch_letter_state_to_numeric(nexusfile[s,i],all_states=all_states);
					}	else if (nexusfile[s,i]>="0" && nexusfile[s,i]<="9") {
					chmatrix[tx,c] <- as.numeric(nexusfile[s,i]);
					}
				if ((chmatrix[tx,c]+1)>nstates[c]) nstates[c] <- chmatrix[tx,c]+1;
				if (i < (endline-1)) i <- i+1;
				}	else {
				if (c==strat)	{
					if (nexusfile[s,i]>="0" && nexusfile[s,i]<='9')	{
						strat_ranges[tx,1]=as.numeric(nexusfile[s,i])
						} else if (nexusfile[s,i]>="A" && nexusfile[s,i]<="Z")	{
						strat_ranges[tx,1]=switch_letter_state_to_numeric(nexusfile[s,i],all_states = all_states);
						}
					if (range==0)	strat_ranges[tx,2] <- strat_ranges[tx,1]
					i <- i+1
					} else if (c==range)	{
					if (nexusfile[s,i]>="0" && nexusfile[s,i]<='9')	{
						strat_ranges[tx,2]=as.numeric(nexusfile[s,i])
						} else if (nexusfile[s,i]>="A" && nexusfile[s,i]<="Z")	{
						strat_ranges[tx,2]=switch_letter_state_to_numeric(nexusfile[s,i],all_states = all_states);
						}
					i <- i+1
					} else if (c==geog)	{
						if (nexusfile[s,i]>="0" && nexusfile[s,i]<='9')	{
							geography[tx]=as.numeric(nexusfile[s,i])
						} else if (nexusfile[s,i]>="A" && nexusfile[s,i]<="Z")	{
							geography[tx]=switch_letter_state_to_numeric(nexusfile[s,i],all_states = all_states);
						}
					}
				} # end non-morphological data
#			print(nexusfile[s,k:83]);
#			print(chmatrix[tx,])
			if (nexusfile[s,i+1]=="\n" || i==(mxln-1)) c <- nchars+extra;
			}
		}
#	chmatrix[tx,];
#	tx <- tx+1;
#	s <- s+1
	}
	#x <- list(taxa,chmatrix,strat_ranges,geography)
	#return (list(taxa,chmatrix,strat_ranges,geography))

chmatrix <- mundify_character_matrix(chmatrix,minst=0,UNKNOWN,INAP);	# clean up coding
nstates <- count_states(chmatrix,UNKNOWN,INAP);

tree_found <- 0;
while (s<length(nexus) && tree_found==0)	{
	while (nexus[s]!= "BEGIN TREES; " && s<length(nexus))	s <- s+1;
	if (s<length(nexus))	{
		while (tree_found==0 && s<length(nexus))	{
			s <- s+1
			jj <- strsplit(nexus[s],split=c("\t"," "),fixed=TRUE)[[1]];
			jj <- paste(jj,collapse="")
			jj <- strsplit(jj,split=" ",fixed=TRUE)[[1]];
			if (sum(jj=="TREE")>0 || sum(jj=="tree")>0)	tree_found <- 1;
			}
#		s <- s+notu;
#		while (jj[i]=="")	jj[i] <- NULL;
#		while (j[1]=="\t")	j <- j[2:length(j)];
#		if (j[1]=="T" && j[2]=="R" && j[3]=="E")	{
#			while (j!="(")	j <- j[2:length(j)];
#			}
		newick_string <- jj[length(jj)];
		newick_string <- fix_newick_ancestors(jj[length(jj)])
		tree <- read_newick_string(newick_string);
		tree_found <- 1
		s <- length(nexus);
		}
	}

row.names(chmatrix) <- taxa;

unscored_taxa <- c();
for (n in 1:notu)	{
	if (sum(chmatrix[n,]==UNKNOWN)==nchars)
		unscored_taxa <- c(unscored_taxa,n);
	}

if (nchars<10)	{
	colnames(chmatrix) <- 1:nchars;
	} else if (nchars<100)	{
	colnames(chmatrix) <- c(paste(0,(1:9),sep=""),10:nchars);
	} else if (nchars<1000)	{
	colnames(chmatrix) <- c(paste(00,(1:9),sep=""),paste(0,(10:99),sep=""),100:nchars);
	}
if (exclude[1]!=-1)	{
	keepers <- (1:nchars)[!(1:nchars) %in% exclude];
	chmatrix <- chmatrix[,keepers];
	nstates <- nstates[keepers];
	state_orders <- state_orders[keepers];
	character_rate_partitions <- character_rate_partitions[keepers];
	}

if (strat!=0 && geog!=0 && tree_found==1)  {
	output <- list(taxa,chmatrix,as.numeric(nstates),state_orders,strat_ranges,geography,tree,as.numeric(outgroup),unscored_taxa,character_rate_partitions,character_trend_partitions);
	names(output) <-  c("OTUs","Matrix","States","State_Types","Stratigraphic_Ranges","Geography","Tree","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
	} else if (strat!=0)  {
	if (geog!=0)  {
		output <- list(taxa,chmatrix,as.numeric(nstates),state_orders,strat_ranges,geography,as.numeric(outgroup),unscored_taxa,character_rate_partitions,character_trend_partitions);
		names(output) <-  c("OTUs","Matrix","States","State_Types","Stratigraphic_Ranges","Geography","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
		} else if (tree_found!=0)	{
		output <- list(taxa,chmatrix,as.numeric(nstates),state_orders,strat_ranges,tree,as.numeric(outgroup),unscored_taxa,character_rate_partitions,character_trend_partitions);
		names(output) <-  c("OTUs","Matrix","States","State_Types","Stratigraphic_Ranges","Tree","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
		} else	{
		output <- list(taxa,chmatrix,as.numeric(nstates),state_orders,strat_ranges,as.numeric(outgroup),unscored_taxa,character_rate_partitions,character_trend_partitions);
		names(output) <-  c("OTUs","Matrix","States","State_Types","Stratigraphic_Ranges","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
		}
	} else if (geog!=0)  {
	if (tree_found!=0)	{
		output <- list(taxa,chmatrix,as.numeric(nstates),state_orders,geography,tree,as.numeric(outgroup),unscored_taxa,character_rate_partitions,character_trend_partitions);
		names(output) <-  c("OTUs","Matrix","States","State_Types","Geography","Tree","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
		} else	{
		output <- list(taxa,chmatrix,as.numeric(nstates),state_orders,geography,as.numeric(outgroup),unscored_taxa,character_rate_partitions,character_trend_partitions);
		names(output) <-  c("OTUs","Matrix","States","State_Types","Geography","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
		}
	} else if (tree_found!=0) {
	output <- list(taxa,chmatrix,as.numeric(nstates),state_orders,tree,as.numeric(outgroup),unscored_taxa,character_rate_partitions,character_trend_partitions);
	names(output) <-  c("OTUs","Matrix","States","State_Types","Tree","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
	} else	{
	output <- list(taxa,chmatrix,as.numeric(nstates),state_orders,as.numeric(outgroup),unscored_taxa,character_rate_partitions,character_trend_partitions);
	names(output) <-  c("OTUs","Matrix","States","State_Types","Outgroup","Unscored_Taxa","Rate_Partitions","Trend_Partitions");
	}

return(output)
}

#### DEAL WITH TRICKY CHARACTERS ####
switch_letter_state_to_numeric <- function(state,all_states=c(0:9,LETTERS[!LETTERS %in% c("I","O")]))  {
# 2017-10-09: now will pass numeric characters through unchanged
# 2019-01-25: simplified greatly!
# 2020-12-01: allows for i's and o's, but assumes that they are not there!
# -1 is for 0 to be zero
return(match(state,all_states)-1);
}

switch_letter_state_to_numeric_old <- function(state)  {
# 2017-10-09: now will pass numeric characters through unchanged
# 2019-01-25: simplified greatly!
if (state > 9)	{
	state <- toupper(state)
	poss_letter_states <- toupper(letters[!letters %in% c("i","o")]);
	return(9+match(state,poss_letter_states));
	} else	{
	return(state);
	}
}

switch_numeric_state_to_letter <- function(state)  {
# 2017-10-09: now will pass numeric characters through unchanged
# 2019-01-25: simplified greatly!
if (state > 9)	{
#	state <- toupper(state)
	poss_letter_states <- toupper(letters[!letters %in% c("i","o")]);
	return(poss_letter_states[state-9]);
	} else	{
	return(state);
	}
}

unravel_polymorph_badass <- function(poly,minst=0)	{
combo <- -1*poly;
state_test <- as.numeric(strsplit(x=as.character(combo),split="")[[1]])
if (state_test==sort(state_test,decreasing = T) && length(unique(state_test))==length(state_test))	{
	sts <- 1+floor(log10(abs(combo)))
	polymorphics <- vector(length=sts)
	base <- 10^(sts-1)
	for (s in 1:sts)	{
		polymorphics[s] <- floor(abs(combo)/base)
		combo <- combo%%base
		base <- base/10
		}
	} else	{
	breakpt <- match(max(state_test),state_test);
	if (breakpt > 2)	{
		i <- 1;
		while (i < breakpt)	{
			j <- i+1;
			state_test[i] <- (10*state_test[i])+state_test[j];
			state_test[j] <- -1;
#			print(state_test);
			i <- j+1;
			}
		polymorphics <- state_test[state_test>=minst];
		} else if (sum(state_test<minst)>0)	{
		# this should happen
		i <- 1;
		while (i < length(state_test))	{
			j <- i+1;
			state_test[i] <- (10*state_test[i])+state_test[j];
			state_test[j] <- -1;
#			print(state_test);
			i <- j+1;
			}
		polymorphics <- state_test[state_test>=minst];
		} else if ((length(state_test) %% 2)==0)	{
		i <- 1;
		while (i < length(state_test))	{
			j <- i+1;
			state_test[i] <- (10*state_test[i])+state_test[j];
			state_test[j] <- -1;
#			print(state_test);
			i <- j+1;
			}
		polymorphics <- state_test[state_test>=minst];
		} else	{
		i <- 1;
		while (i < length(state_test))	{
			j <- i+1;
			state_test[i] <- (10*state_test[i])+state_test[j];
			state_test[j] <- -1;
#			print(state_test);
			i <- j+1;
			}
		polymorphics <- state_test[state_test>=minst];
		}
	}
return (polymorphics);
}

unravel_polymorph <- function(poly)	{
combo <- -1*poly
sts <- 1+floor(log10(abs(combo)))
polymorphics <- vector(length=sts)

base <- 10^(sts-1)
for (s in 1:sts)	{
	polymorphics[s] <- floor(abs(combo)/base)
	combo <- combo%%base
	base <- base/10
	}
return (polymorphics)
}

ravel_polymorph <- function(polystates)	{
polystates <- sort(polystates,decreasing = TRUE);
polym <- polystates[1];
for (st in 2:length(polystates))	polym <- (10*polym)+polystates[st]
return(-1*polym)
}

ravel_polymorph_for_file <- function(polystates)	{
polystates <- sort(polystates,decreasing = FALSE);
return(paste("(",paste(polystates,collapse=""),")",sep=""));
}

#### SUMMARIZE CHARACTER DATA ####
# count taxa scored with something other than missing or inapplicable
count_scored_characters_per_otu <- function(chmatrix,UNKNOWN=-11,INAP=-22)	{
nch <- ncol(chmatrix);
notu <- nrow(chmatrix);
#scored <- vector(length=nch)
scored <- c();
for (s in 1:notu)
	scored <- c(scored,notu - (sum(chmatrix[s,]==UNKNOWN)+sum(chmatrix[s,]==INAP)));
return(scored);
}

# count missing and/or inapplicable per otu
count_scored_otus_per_character <- function(chmatrix,UNKNOWN=-11,INAP=-22)	{
if (is.matrix(chmatrix))	{
	nchars <- ncol(chmatrix);
	notu <- nrow(chmatrix);
	} else	{
	nchars <- 1;
	notu <- length(chmatrix);
	dummy <- array(0,dim=c(length(chmatrix),1));
	dummy[,1] <- chmatrix;
	chmatrix <- dummy;
	}
#scored <- vector(length=nch)
scored <- c();
for (c in 1:nchars)
	scored <- c(scored,notu - (sum(chmatrix[,c]==UNKNOWN)+sum(chmatrix[,c]==INAP)));
return(scored);
}

# count missing and/or inapplicable per otu
count_scored_otus_per_character_state <- function(chmatrix,chstates,UNKNOWN=-11,INAP=-22)	{
nch <- ncol(chmatrix);
#notu <- nrow(chmatrix);
#scored <- vector(length=nch)
scored <- array(0,dim=c(nch,max(chstates)))
for (c in 1:nch)	{
	for (st in 1:chstates[c])	{
		stt <- st-1
		scored[c,st] <- sum(chmatrix[,c]==stt)
		}
	}
return(scored);
}

# count characters with autapomorphic taxa
count_autapomorphic_characters <- function(chmatrix,chstates,UNKNOWN=-11,INAP=-22)	{
# count number of characters with an autapomorphic state
nchars <- ncol(chmatrix);
notus_per_chstate <- count_scored_otus_per_character_state(chmatrix,chstates);
autaps <- 0
for (c in 1:nchars)	{
	if(sum(notus_per_chstate[c,1:chstates[c]]==1)>0)	{
		autaps <- autaps+1;
		}
	}
return(autaps);
}

# count states coding only one taxon
count_autapomorphic_states <- function(chmatrix,chstates,UNKNOWN=-11,INAP=-22)	{
# count number of states that are autapomorphic
nchars <- ncol(chmatrix);
notus_per_chstate <- count_scored_otus_per_character_state(chmatrix,chstates);
autaps <- 0
for (ch in 1:nchars)
	autaps <- autaps+sum(notus_per_chstate[ch,1:chstates[ch]]==1);
return(autaps);
}

# routine to list all characters with at least one autapomorphic state
list_autapomorphic_characters <- function(chmatrix,chstates,UNKNOWN=-11,INAP=-22)	{
nchars <- ncol(chmatrix);
#scored <- vector(length=nch)
notus_per_chstate <- count_scored_otus_per_character_state(chmatrix,chstates);
autaps <- c()
for (ch in 1:nchars)	{
	if(sum(notus_per_chstate[ch,1:chstates[ch]]==1)>0)	{
		autaps <- c(autaps,ch);
		}
	}
return(autaps);
}

# routine to list all character states that are autapomorphic
list_autapomorphic_states <- function(chmatrix,chstates,UNKNOWN=-11,INAP=-22)	{
nchars <- ncol(chmatrix);
#scored <- vector(length=nch)
notus_per_chstate <- count_scored_otus_per_character_state(chmatrix,chstates);
autaps <- c()
for (c in 1:nchars)	{
	if(sum(notus_per_chstate[c,1:chstates[c]]==1)>0)	{
		autaps <- rbind(autaps,c(c,chstates[c]));
		}
	}
return(autaps);
}

# count the number of scorings that are polymorphic. There can be 1 per character per taxon
count_polymorphic_scorings <- function(chmatrix,UNKNOWN=-11,INAP=-22)	{
# count number of characters with an autapomorphic state
if (!is.matrix(chmatrix))	{
	chmatrix <- data.frame(ch=chmatrix);
	}
nchars <- ncol(chmatrix);
polymorphs <- vector(length=nchars);
for (ch in 1:nchars)	{
	char_states <- chmatrix[,ch];
	char_states <- char_states[char_states!=UNKNOWN];
	char_states <- char_states[char_states!=INAP];
	polys <- char_states[char_states<0];
	polymorphs[ch] <- length(polys);
	}
return(polymorphs);
}

# count the number of characters with a polymorphic scoring
count_polymorphic_characters <- function(chmatrix,UNKNOWN=-11,INAP=-22)	{
# count number of characters with an autapomorphic state
if (!is.matrix(chmatrix))	{
	chmatrix <- data.frame(ch=chmatrix);
	}
nchars <- ncol(chmatrix);
polymorphs <- count_polymorphic_scorings(chmatrix,UNKNOWN,INAP);
return(sum(polymorphs>0));
}

# count the number of states per character that show polymorphism
count_polymorphic_states <- function(chmatrix,UNKNOWN=-11,INAP=-22)	{
# count number of characters with an autapomorphic state
if (!is.matrix(chmatrix))
	chmatrix <- data.frame(ch=chmatrix);
nchars <- ncol(chmatrix);
polymorphic_states <- vector(length=nchars);
for (ch in 1:nchars)	{
	char_states <- sort(unique(chmatrix[,ch]));
	char_states <- char_states[char_states!=UNKNOWN];
	char_states <- char_states[char_states!=INAP];
	polys <- char_states[char_states<0];
	poly_states <- c();
	pp <- 0;
	while (pp < length(polys))	{
		pp <- pp+1;
		poly_states <- unique(c(poly_states,unravel_polymorph_badass(polys[pp])));
		}
	polymorphic_states[ch] <- length(polys);
	}
return(polymorphic_states);
}

# get number of states for each character
count_states <- function(chmatrix,UNKNOWN=-11,INAP=-22,include_polymorphs=T)	{
if (!is.matrix(chmatrix))	{
	chmatrix <- data.frame(ch=chmatrix);
	}
nchars <- ncol(chmatrix);
nstates <- c();
for (ch in 1:nchars)	{
	char_states <- sort(unique(chmatrix[,ch]));
	char_states <- char_states[char_states!=UNKNOWN];
	char_states <- char_states[char_states!=INAP];
	unique_states <- char_states[char_states>=0]
	if (sum(char_states<0)>0 && include_polymorphs)	{
		polys <- char_states[char_states<0];
		polystates <- c();
		for (pp in polys)	{
			polystates <- c(polystates,unravel_polymorph(poly=pp));
			}
		unique_states <- sort(unique(c(unique_states,polystates)));
		}
	nstates <- c(nstates,length(unique_states));
	}	# pick up here!!!
return(nstates);
}

count_states_old <- function(chmatrix,UNKNOWN=-11,INAP=-22)	{
# 2020-09-01: fixed breakdown of polymorphics
nchars <- ncol(chmatrix);
nstates <- c();
for (ch in 1:nchars)	{
	char_states <- sort(unique(chmatrix[,ch]));
	char_states <- char_states[char_states!=UNKNOWN];
	char_states <- char_states[char_states!=INAP];
	if (sum(char_states<0)>0)	{
		while (char_states[1]<0)	{
			char_states <- sort(unique(c(char_states[2:length(char_states)],unravel_polymorph(char_states[1]))));
			}
		}
	nstates <- c(nstates,length(char_states));
	}
return(nstates);
}

maximum_state <- function(chmatrix,UNKNOWN=-11,INAP=-22)	{
# 2020-09-01: fixed breakdown of polymorphics
nchars <- ncol(chmatrix);
maxstates <- c();
for (ch in 1:nchars)	{
	char_states <- sort(unique(chmatrix[,ch]));
	char_states <- char_states[char_states!=UNKNOWN];
	char_states <- char_states[char_states!=INAP];
	if (length(char_states)==0)	{
		char_states <- 0;
		}	else if (sum(char_states<0)>0)	{
		while (char_states[1]<0)	{
			char_states <- sort(unique(c(char_states[2:length(char_states)],unravel_polymorph(char_states[1]))));
			}
		}
	maxstates <- c(maxstates,max(char_states));
	}
return(maxstates);
}

minimum_state <- function(chmatrix,UNKNOWN=-11,INAP=-22)	{
# 2020-09-01: fixed breakdown of polymorphics
nchars <- ncol(chmatrix);
minstates <- c();
for (ch in 1:nchars)	{
	char_states <- sort(unique(chmatrix[,ch]));
	char_states <- char_states[char_states!=UNKNOWN];
	char_states <- char_states[char_states!=INAP];
	if (length(char_states)==0)	{
		char_states <- 0;
		} else if (sum(char_states<0)>0)	{
		while (char_states[1]<0)	{
			char_states <- sort(unique(c(char_states[2:length(char_states)],unravel_polymorph(char_states[1]))));
			}
		}
	minstates <- c(minstates,min(char_states));
	}
return(minstates);
}


#### MODIFY & WRITE NEXUS FILES ####
# routine to "clean" character matrix (e.g., remove gaps in coding, standarize minimum states, etc.)
mundify_character_matrix <- function(chmatrix,minst=0,UNKNOWN=-11,INAP=-22)	{
notu <- nrow(chmatrix);	# replaces spc to standardize coding.
nchars <- ncol(chmatrix);
min_states <- minimum_state(chmatrix,UNKNOWN=UNKNOWN,INAP=INAP);
max_states <- maximum_state(chmatrix,UNKNOWN=UNKNOWN,INAP=INAP);
for (ch in 1:nchars)	{
	rem <- c((1:notu)[chmatrix[,ch]==UNKNOWN],(1:notu)[chmatrix[,ch]==INAP]);
	if (length(rem)>0)	{
		test <- chmatrix[-rem,ch]
		}	else test <- chmatrix[,ch];
	# check polymorphics for anything that needs to be changed
	if (length(rem) < notu)	{
		polys <- sum(test<0);	# taxa with polymorphic scores
		coded <- sort(unique(test[test>=0]));
		if (polys>0)	{
#			examps <- test[test<0]
			polycoded <- sort(unique(test[test<0]))
			for (i in 1:length(polycoded))	{
				polystates <- unravel_polymorph(polycoded[i])
				coded <- sort(unique(c(coded,polystates)))
#				if (min(polystates)<minstch)	minstch <- min(polystates)
				}
			} else	{
			polycoded <- c();
			}
		minstch <- min(coded);
	# eliminate gaps in states
		if (sum(!min(coded):max(coded) %in% coded)>0)	{
#			new_codes <- match(coded,coded)-(1-minst);
			new_codes <- match(coded,coded)-(1-minstch);
			for (st in 1:length(coded))	{
				if (coded[st]!=new_codes[st])	{
					rec <- (1:notu)[chmatrix[,ch]==coded[st]];
					chmatrix[rec,ch] <- new_codes[st];
					redo_poly <- 1;
					while (redo_poly <= length(polycoded))	{
						polystates <- unravel_polymorph(polycoded[redo_poly]);
						polystates[polystates==coded[st]] <- new_codes[st]
						newpolystates <- ravel_polymorph(polystates);
						testp <- (1:notu)[chmatrix[,ch]==polycoded[redo_poly]];
#						if (newpolystates != polycoded[redo_poly])	{
						chmatrix[testp,ch] <- newpolystates;
#							}
#						polycoded[redo_poly] <- chmatrix[,ch][chmatrix[,ch] %in% polycoded[redo_poly]] <- newpolystates
						polycoded[redo_poly] <- newpolystates;
						redo_poly <- redo_poly+1;
						}
					coded[st] <- new_codes[st];
					}
				}
			}
		# standardize minimum state
		# simple cheat: subtract 1111 to polymorphics
		if (minstch!=minst)	{
			adj <- minst-minstch;
			test2 <- (1:notu)[chmatrix[,ch]>=0];
			chmatrix[test2,ch] <- chmatrix[test2,ch]+adj;
			if (polys>0)	{
				examps2 <- polycoded;
				for (i in 1:length(polycoded))	{
					testp <- (1:notu)[chmatrix[,ch]==polycoded[i]];
					examps2[i] <- examps2[i]-(adj*floor(10^floor(log10(abs(polycoded[i])))*10/9));
					chmatrix[testp,ch] <- examps2[i];
					}
				polycoded <- examps2;
				} # end rescoring of polytomies
			} # end rescaling stats
		} # end case where rem < notu
	} # end search of characters;	
return(chmatrix);
}

# routine to remove invariant and/or unscored characters from matrix
remove_invariants_from_character_matrix <- function(chmatrix,minst=0,UNKNOWN=-11,INAP=-22)	{
notu <- nrow(chmatrix)	# replaces spc to standardize coding.
ncharss <- ncol(chmatrix)
rem_char <- c()
for (ch in 1:ncharss)	{
	rem <- c((1:notu)[chmatrix[,ch]==UNKNOWN],(1:notu)[chmatrix[,ch]==INAP])
	if (length(rem)>0)	{
		test <- chmatrix[-rem,ch]
		} else	test <- chmatrix[,ch]
	if (length(unique(test))<2)	rem_char <- c(rem_char,ch)
	}
return(chmatrix[,-rem_char])
}

# generate composite score from 2+ scored taxa
accersi_composite_scores <- function(mini_matrix,return_polymorph=TRUE,UNKNOWN=-11,INAP=-22)	{
notu <- nrow(mini_matrix);
nchars <- ncol(mini_matrix);
composite_score <- c();
ch <- 1;
for (ch in 1:nchars)	{
	composite_states <- unique(mini_matrix[,ch]);
	if (length(composite_states)==1)	{
		composite_score <- c(composite_score,composite_states);
		} else	{
		composite_states <- composite_states[composite_states!=INAP];
		if (length(composite_states)>1)	{
			composite_states <- composite_states[composite_states!=UNKNOWN];
			if (length(composite_states)>1)	{
				polyscored <- composite_states[composite_states<0];
				if (length(polyscored)>0)	{
					accersi_states <- sapply(polyscored,unravel_polymorph);
					composite_states <- sort(unique(c(accersi_states,composite_states[composite_states>=0])));
					}
				composite_score <- c(composite_score,ravel_polymorph(composite_states));
				} else	{
				composite_score <- c(composite_score,composite_states);
				}
			} else	{
			composite_score <- c(composite_score,composite_states);
			}
		}
#	ch <- ch+1;
#	composite_score;
	}
#rbind(mini_matrix,composite_score)
return(composite_score);
}

scribio_nexus_file_from_chmatrix <- function(ch_matrix,new_file_name,UNKNOWN=-11,INAP=-22)	{
notu <- nrow(ch_matrix);
taxon_names <- rownames(ch_matrix);
nchars <- ncol(ch_matrix);
nstates <- count_states(chmatrix = ch_matrix);

nexus_file_content <- c();
nexus_file_content <- rbind("#NEXUS","","BEGIN DATA;")
nexus_file_content <- rbind(nexus_file_content,paste("	DIMENSIONS  NTAX=",notu," NCHAR=",nchars,";",sep=""));
if (max(nstates)<10) {
	state_symbols <- " ";
	for (st in 1:max(nstates))
		state_symbols <- paste(state_symbols,st-1,sep=" ");	
	} else	{
	mxl <- max(nstates)-10;
	letter_states <- LETTERS[!LETTERS %in% c("I","O")][1:mxl]
	all_states <- c(0:9,letter_states);
	state_symbols <- paste(all_states,collapse=" ");
	}
nexus_file_content <- rbind(nexus_file_content,paste("	FORMAT DATATYPE = STANDARD RESPECTCASE GAP = - MISSING = ? SYMBOLS = \"",state_symbols,"\";"));
nexus_file_content <- rbind(nexus_file_content,"	MATRIX");

string_to_count <- taxon_names;
name_lengths <- sapply(string_to_count,count_characters_in_string);
max_name_length <- max(name_lengths);
need_quotes <- c(".","(",")","[","]");
for (nn in 1:notu)	{
	test_name <- strsplit(taxon_names[nn],split="",fixed=TRUE)[[1]]
	if (sum(test_name %in% need_quotes)==0)	{
		taxon <- gsub(" ","_",taxon_names[nn]);
		} else	{
		taxon <- paste("\"",taxon_names[nn],"\"",sep="");
		name_lengths[nn] <- name_lengths[nn]+2;
		}
	this_line <- paste("\t",taxon,paste(rep(" ",(5+(max_name_length-name_lengths[nn]))),collapse=""),sep="");
	otu_code <- c();
	for (ch in 1:nchars)	{
		if (ch_matrix[nn,ch]>=0 && ch_matrix[nn,ch]<=9)	{
			otu_code <- paste(otu_code,ch_matrix[nn,ch],sep="");
			} else if (ch_matrix[nn,ch]>9)	{
			otu_code <- paste(otu_code,all_states[1+ch_matrix[nn,ch]],sep="");	# note: we need +1 because of state 0
			} else if (ch_matrix[nn,ch]==UNKNOWN)	{
			otu_code <- paste(otu_code,"?",sep="");
			} else if (ch_matrix[nn,ch]==INAP)	{
			otu_code <- paste(otu_code,"-",sep="");
			} else if (ch_matrix[nn,ch]<0)	{
			polystates <- strsplit(as.character(ch_matrix[nn,ch]),split="",fixed=TRUE)[[1]];
			polystates <- as.numeric(polystates[polystates!="-"]);
			otu_code <- paste(otu_code,ravel_polymorph_for_file(polystates),sep="");
			}
		}
	nexus_file_content <- rbind(nexus_file_content,paste(this_line,otu_code,sep=""));
	}

nexus_file_content <- rbind(nexus_file_content,";");
nexus_file_content <- rbind(nexus_file_content,"END;");
nexus_file_content <- rbind(nexus_file_content,"begin mrbayes;");
nexus_file_content <- rbind(nexus_file_content,"	set autoclose=yes nowarn=yes;");
nexus_file_content <- rbind(nexus_file_content,"	lset nst=6 rates=invgamma;");
nexus_file_content <- rbind(nexus_file_content,"	unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all); ");
nexus_file_content <- rbind(nexus_file_content,"	prset applyto=(all) ratepr=variable;");
nexus_file_content <- rbind(nexus_file_content,"	mcmcp ngen= 100000000 relburnin=yes burninfrac=0.25 printfreq=10000  samplefreq=10000 nchains=4 savebrlens=yes;");
nexus_file_content <- rbind(nexus_file_content,"	mcmc;");
nexus_file_content <- rbind(nexus_file_content,"	sumt;");
nexus_file_content <- rbind(nexus_file_content,"end;");
write(nexus_file_content,file=new_file_name);
}

ravel_polymorph_for_file <- function(polystates)	{
polystates <- sort(polystates,decreasing = FALSE);
return(paste("(",paste(polystates,collapse=""),")",sep=""));
}

#### READ NEWICK FILES ####
#### convert (1,(2,3)) to vector_tree = 4 5 5 -1 4
read_newick_tree_from_chosen_file <- function() {
newicktree_file <- file.choose();
newick_tree <- scan(file=newicktree_file,what=character(),sep="\n");
nexus_string <- strsplit(newick_tree,split="",fixed=TRUE)[[1]]
nodes <- 0
for (i in 1:length(nexus_string))		if (nexus_string[i]=="(")	nodes <- nodes+1
# get clades
clades <- vector(length=nodes)
for (c in 1:nodes)	clades[c] <- c
# get taxa
notu <- p <- 0
for (i in 1:length(nexus_string))	{
	if (nexus_string[i]>="0" && nexus_string[i]<="9")	{
		otu <- as.numeric(nexus_string[i])+(otu * (10^p))
		p <- p+1
		if (otu>notu)	notu <- otu
		} else {
		p <- otu <- 0
		}
	}
vector_tree <- vector(length=notu+max(clades))
for (c in 1:nodes)	clades[c] <- -1
cl <- c <- 0
i <- 1
for (i in 1:length(nexus_string))	{
	if (nexus_string[i]=="(")	{
		sp <- p <- 0
		cl <- cl+1
		if (cl>1)	{
			vector_tree[notu+cl] <- clades[c]+notu
			} else vector_tree[notu+1] <- -1
		c <- c+1
		clades[c] <- cl
		} else if (nexus_string[i]==")")	{
		c <- c-1
		sp <- p <- 0
		} else if (nexus_string[i]==",")	{
		sp <- p <- 0
		} else if (nexus_string[i]>="0" && nexus_string[i]<="9")	{
		sp <- as.numeric(nexus_string[i])+(sp*10)
		p <- p+1
		if (nexus_string[i+1]<"0" || nexus_string[i]>"9")	vector_tree[sp] <- notu+clades[c]
		}
	}

return(vector_tree)
}

read_newick_tree_from_file <- function(newicktree_file) {
newick_tree <- scan(file=newicktree_file,what=character(),sep="\n")
nexus_string <- strsplit(newick_tree,split="",fixed=TRUE)[[1]]
nodes <- 0
for (i in 1:length(nexus_string))		if (nexus_string[i]=="(")	nodes <- nodes+1
# get clades
clades <- 1:nodes;
#clades <- vector(length=nodes)
#for (c in 1:nodes)	clades[c] <- c;
# get taxa
notu <- p <- 0;
for (i in 1:length(nexus_string))	{
	if (nexus_string[i]>="0" && nexus_string[i]<="9")	{
		otu <- as.numeric(nexus_string[i])+(otu * 10)
		p <- p+1
		if (otu>notu)	notu <- otu
		} else {
		p <- otu <- 0
		}
	}
vector_tree <- vector(length=notu+max(clades))
for (c in 1:nodes)	clades[c] <- -1
cl <- c <- 0
i <- 1
for (i in 1:length(nexus_string))	{
	if (nexus_string[i]=="(")	{
		sp <- p <- 0
		cl <- cl+1
		if (cl>1)	{
			vector_tree[notu+cl] <- clades[c]+notu
			} else vector_tree[notu+1] <- -1
		c <- c+1
		clades[c] <- cl
		} else if (nexus_string[i]==")")	{
		c <- c-1
		sp <- p <- 0
		} else if (nexus_string[i]==",")	{
		sp <- p <- 0
		} else if (nexus_string[i]>="0" && nexus_string[i]<="9")	{
		sp <- as.numeric(nexus_string[i])+(sp*10);
		p <- p+1
		if (nexus_string[i+1]<"0" || nexus_string[i]>"9")	vector_tree[sp] <- notu+clades[c];
		}
	}

return(vector_tree)
}

#### convert vector_tree = 4 5 5 -1 4 to (1,(2,3))
#### 	where number is the htu number of the clade to which a species or htu belong
#### does not work yet
write_newick_string_from_vector_tree <- function(vector_tree) {
mat_tree <- transform_vector_tree_to_matrix_tree(vector_tree);

nodes <- 0;
if (length(newick_string)==1)	newick_string <- strsplit(newick_string,split="",fixed=TRUE)[[1]];
for (i in 1:length(newick_string))		if (newick_string[i]=="(")	nodes <- nodes+1;
# get clades
clades <- vector(length=nodes);
for (c in 1:nodes)	clades[c] <- c;
# get taxa
notu <- p <- 0
for (i in 1:length(newick_string))	{
	if (newick_string[i]>="0" && newick_string[i]<="9")	{
		otu <- as.numeric(newick_string[i])+(otu * (10^p))
		p <- p+1
		if (otu>notu)	notu <- otu
		} else {
		p <- otu <- 0
		}
	}
vector_tree <- vector(length=notu+max(clades))
for (c in 1:nodes)	clades[c] <- -1
cl <- c <- 0
i <- 1
for (i in 1:length(newick_string))	{
	if (newick_string[i]=="(")	{
		sp <- p <- 0
		cl <- cl+1
		if (cl>1)	{
			vector_tree[notu+cl] <- clades[c]+notu
			} else vector_tree[notu+1] <- -1
		c <- c+1
		clades[c] <- cl
		} else if (newick_string[i]==")")	{
		c <- c-1
		sp <- p <- 0
		} else if (newick_string[i]==",")	{
		sp <- p <- 0
		} else if (newick_string[i]>="0" && newick_string[i]<="9")	{
		sp <- as.numeric(newick_string[i])+(sp*10)
		p <- p+1
		if (newick_string[i+1]<"0" || newick_string[i]>"9")	vector_tree[sp] <- notu+clades[c]
		}
	}

return(vector_tree)
}

read_newick_string <- function(newick_string) {
nodes <- 0;
if (length(newick_string)==1)	newick_string <- strsplit(newick_string,split="",fixed=TRUE)[[1]];
for (i in 1:length(newick_string))		if (newick_string[i]=="(")	nodes <- nodes+1;
# get clades
clades <- vector(length=nodes);
for (c in 1:nodes)	clades[c] <- c;
# get taxa
notu <- p <- 0
for (i in 1:length(newick_string))	{
	if (newick_string[i]>="0" && newick_string[i]<="9")	{
		otu <- as.numeric(newick_string[i])+(otu * (10^p))
		p <- p+1
		if (otu>notu)	notu <- otu
		} else {
		p <- otu <- 0
		}
	}
vector_tree <- vector(length=notu+max(clades))
for (c in 1:nodes)	clades[c] <- -1
cl <- c <- 0
i <- 1
for (i in 1:length(newick_string))	{
	if (newick_string[i]=="(")	{
		sp <- p <- 0
		cl <- cl+1
		if (cl>1)	{
			vector_tree[notu+cl] <- clades[c]+notu
			} else vector_tree[notu+1] <- -1
		c <- c+1
		clades[c] <- cl
		} else if (newick_string[i]==")")	{
		c <- c-1
		sp <- p <- 0
		} else if (newick_string[i]==",")	{
		sp <- p <- 0
		} else if (newick_string[i]>="0" && newick_string[i]<="9")	{
		sp <- as.numeric(newick_string[i])+(sp*10)
		p <- p+1
		if (newick_string[i+1]<"0" || newick_string[i]>"9")	vector_tree[sp] <- notu+clades[c]
		}
	}

return(vector_tree)
}

# written for Cinctan project
transform_newick_string_to_venn_tree <- function(newick_string)	{
atomized_newick <- strsplit(newick_string,"")[[1]];
l_m <- length(atomized_newick);
clade_bounds_l <- (1:l_m)[atomized_newick %in% "("];
clade_bounds_e <- clade_bounds_r <- (1:l_m)[atomized_newick %in% ")"];
nNodes <- length(clade_bounds_l);	# number of clades;
names(clade_bounds_l) <- names(clade_bounds_r) <- names(clade_bounds_e) <- 1:nNodes;
# get the first possible right paren ending this clade;
for (nn in 1:nNodes)	clade_bounds_e[nn] <- sum(clade_bounds_r>clade_bounds_l[nn])
clade_bounds_e_unq <- unique(clade_bounds_e);
clade_boundaries <- clade_boundaries_orig <- cbind(clade_bounds_l,clade_bounds_r,clade_bounds_e);
#1 <- nn <- 1
while (length(clade_bounds_e_unq)>0)	{
	this_group <- sum(clade_boundaries_orig[,3]==clade_bounds_e_unq[1]);
	if (length(clade_bounds_e_unq)>1)	{
		this_group_starts <- this_group-sum(clade_boundaries_orig[,2]<clade_boundaries_orig[this_group+1,1]);
		} else	{
		this_group_starts <- 0;
		}
	clade_boundaries_orig[1:this_group,2] <- sort(clade_boundaries_orig[1:this_group,2],decreasing=T);
	clade_boundaries[clade_boundaries[,1] %in% clade_boundaries_orig[,1],] <- clade_boundaries_orig;
	clade_boundaries_orig <- clade_boundaries_orig[!(1:nrow(clade_boundaries_orig)) %in% ((this_group_starts+1):this_group),]
	if (nrow(clade_boundaries_orig)>0)	clade_boundaries_orig[,2] <- sort(clade_boundaries_orig[,2]);
	n <- 0;
	while (n < nrow(clade_boundaries_orig))	{
		n <- n+1;
		clade_boundaries_orig[n,3] <- sum(clade_boundaries_orig[,2]>clade_boundaries_orig[n,1]);
		}
	clade_bounds_e_unq <- unique(clade_boundaries_orig[,3]);
	}

venn_tree <- array(0,dim=c(nNodes,nNodes+1));
for (nn in 1:nNodes)	{
	lp <- clade_boundaries[nn,1];
	rp <- clade_boundaries[nn,2];
	this_clade <- paste(atomized_newick[lp:rp],collapse="");
#	print(this_clade)
	this_clade <- gsub("\\(","",this_clade);
	this_clade <- gsub(")","",this_clade);
	these_prog <- sort(as.numeric(str_split(this_clade,",")[[1]]));
	venn_tree[nn,1:length(these_prog)] <- these_prog;
	}
return(venn_tree);
}

# written for Cinctan project
accersi_clade_boundaries_from_newick_string <- function(newick_string)	{
atomized_newick <- strsplit(newick_string,"")[[1]];
l_m <- length(atomized_newick);
clade_bounds_l <- (1:l_m)[atomized_newick %in% "("];
clade_bounds_e <- clade_bounds_r <- (1:l_m)[atomized_newick %in% ")"];
nNodes <- length(clade_bounds_l);	# number of clades;
names(clade_bounds_l) <- names(clade_bounds_r) <- names(clade_bounds_e) <- 1:nNodes;
# get the first possible right paren ending this clade;
for (nn in 1:nNodes)	clade_bounds_e[nn] <- sum(clade_bounds_r>clade_bounds_l[nn])
clade_bounds_e_unq <- unique(clade_bounds_e);
clade_boundaries <- clade_boundaries_orig <- cbind(clade_bounds_l,clade_bounds_r,clade_bounds_e);
#1 <- nn <- 1
while (length(clade_bounds_e_unq)>0)	{
	this_group <- sum(clade_boundaries_orig[,3]==clade_bounds_e_unq[1]);
	if (length(clade_bounds_e_unq)>1)	{
		this_group_starts <- this_group-sum(clade_boundaries_orig[,2]<clade_boundaries_orig[this_group+1,1]);
		} else	{
		this_group_starts <- 0;
		}
	clade_boundaries_orig[1:this_group,2] <- sort(clade_boundaries_orig[1:this_group,2],decreasing=T);
	clade_boundaries[clade_boundaries[,1] %in% clade_boundaries_orig[,1],] <- clade_boundaries_orig;
	clade_boundaries_orig <- clade_boundaries_orig[!(1:nrow(clade_boundaries_orig)) %in% ((this_group_starts+1):this_group),]
	if (nrow(clade_boundaries_orig)>0)	clade_boundaries_orig[,2] <- sort(clade_boundaries_orig[,2]);
	n <- 0;
	while (n < nrow(clade_boundaries_orig))	{
		n <- n+1;
		clade_boundaries_orig[n,3] <- sum(clade_boundaries_orig[,2]>clade_boundaries_orig[n,1]);
		}
	clade_bounds_e_unq <- unique(clade_boundaries_orig[,3]);
	}
clade_boundaries <- data.frame(lp=as.numeric(clade_boundaries[,1]),
							   rp=as.numeric(clade_boundaries[,2]))
return(clade_boundaries);
}

#### 	where number is the htu number of the clade to which a species or htu belong
# newick_string_ancestored <- newick_string_taxa_only_raw
# written for Cinctan project
# updated 2020-12-30 to allow outgroup to be ancestral
find_newick_ancestors <- function(newick_string_ancestored)	{
atomic_ancestral <- strsplit(newick_string_ancestored,"")[[1]];
a_a <- length(atomic_ancestral);
l_paren <- (1:a_a)[atomic_ancestral=="("];
r_paren <- (1:a_a)[atomic_ancestral==")"];
sisters <- (1:a_a)[atomic_ancestral==","];
otu_nos <- (1:a_a)[!(1:a_a) %in% c(l_paren,r_paren,sisters)];
otu_nos <- otu_nos[otu_nos!=length(atomic_ancestral)];
notu <- 1;
breaks <- c();
for (i in 2:length(otu_nos))	{
	if ((otu_nos[i]-1)>otu_nos[i-1])	{
		notu <- notu+1;
		breaks <- c(breaks,otu_nos[i]-1);
		}
	}
sampled_ancestors <- array(0,dim=notu);
ancestral_starts <- 1+breaks[atomic_ancestral[breaks]==")"];
ab <- 0;
while (ab < length(ancestral_starts))	{
	ab <- ab+1;
	dd <- ancestral_starts[ab];
	st_hr <- match(dd,otu_nos);
	this_anc <- as.numeric(atomic_ancestral[otu_nos[st_hr]]);
	while(st_hr < length(otu_nos) && otu_nos[st_hr+1]==(otu_nos[st_hr]+1))	{
		st_hr <- st_hr+1;
		this_anc <- (10*this_anc)+as.numeric(atomic_ancestral[otu_nos[st_hr]]);
		}
	sampled_ancestors[this_anc] <- 1;
	}
return(sampled_ancestors);
}

# written for Cinctan project
fix_newick_ancestors <- function(newick_string_ancestored)	{
atomic_ancestral <- strsplit(newick_string_ancestored,"")[[1]];
a_a <- length(atomic_ancestral);
l_paren <- (1:a_a)[atomic_ancestral=="("];
r_paren <- (1:a_a)[atomic_ancestral==")"];
sisters <- (1:a_a)[atomic_ancestral==","];
otu_nos <- (1:a_a)[!(1:a_a) %in% c(l_paren,r_paren,sisters)];
otu_nos <- otu_nos[otu_nos!=length(atomic_ancestral)];
for (rp in length(r_paren):1)	{
	if (!is.na(match(1,otu_nos-r_paren[rp])))	{
		an <- r_paren[rp];
		an_no <- c();
		while (atomic_ancestral[an+1] %in% as.character(0:9))	{
			an_no  <- c(an_no,an+1);
			an <- an+1;
			}
		atomic_ancestral <- c(atomic_ancestral[1:(r_paren[rp]-1)],
							  ",",
							  atomic_ancestral[an_no],
							  ")",
							  atomic_ancestral[(an+1):a_a]);
		a_a <- length(atomic_ancestral);
		l_paren <- (1:a_a)[atomic_ancestral=="("];
		r_paren <- (1:a_a)[atomic_ancestral==")"];
		sisters <- (1:a_a)[atomic_ancestral==","];
		otu_nos <- (1:a_a)[!(1:a_a) %in% c(l_paren,r_paren,sisters)];
		otu_nos <- otu_nos[otu_nos!=length(atomic_ancestral)];
		}
	}
revised_newick_string <- paste(atomic_ancestral,collapse="");
return(revised_newick_string);
}

# written for Cinctan project
# heavily modified 2020-12
#read_newick_string_mcmc <- function(newick_string_full,otu_names) {
read_newick_string_mcmc <- function(newick_string_full,otu_names) {
otu_names_nex <- gsub(" ","_",otu_names);
simple_newick_string <- molecularize <- strsplit(newick_string_full,split="")[[1]];
left_brackets <- (1:length(molecularize))[molecularize %in% "["];
right_brackets <- (1:length(molecularize))[molecularize %in% "]"];
for (br in 1:length(left_brackets))
	molecularize[left_brackets[br]:right_brackets[br]] <- "";
newick_string_taxa_only_rawwest <- newick_string_taxa_raw <- newick_string_taxa_only <- paste(molecularize[molecularize!=""],collapse="");
notu <- length(otu_names);
branch_durations <- array(0,dim=(2*notu)-1);
for (nn in 1:notu)	{
	dummy_newick <- gsub(paste(as.character(otu_names_nex[nn]),":",sep=""),"•",newick_string_taxa_only);
	dummy_newick <- strsplit(dummy_newick,split="")[[1]];
	dd <- 1+match("•",dummy_newick);
	b_d <- dummy_newick[dd];
	dd <- dd+1;
	while (dummy_newick[dd] %in% c(".",0:9))	{
		b_d <- paste(b_d,dummy_newick[dd],sep="");
		dd <- dd+1;
		}
	branch_durations[nn] <- as.numeric(b_d);
	}
for (i in 0:9)	newick_string_taxa_only <- gsub(i,"",newick_string_taxa_only);
newick_string_taxa_only <- gsub(":","",newick_string_taxa_only);
newick_string_taxa_only <- gsub("-","",newick_string_taxa_only);
newick_string_taxa_only <- gsub("\\.","",newick_string_taxa_only);
for (nn in 1:notu)	{
	newick_string_taxa_only <- gsub(otu_names_nex[nn],as.character(nn),newick_string_taxa_only);
	newick_string_taxa_raw <- gsub(otu_names_nex[nn],as.character(nn),newick_string_taxa_raw);
	}
newick_string_taxa_only_atomized <- strsplit(newick_string_taxa_only,"")[[1]];
nstoa <- length(newick_string_taxa_only_atomized);
if (newick_string_taxa_only_atomized[nstoa]!=";")	{
	newick_string_taxa_only_atomized <- c(newick_string_taxa_only_atomized,";");
	nstoa <- length(newick_string_taxa_only_atomized);
	}
if (newick_string_taxa_only_atomized[nstoa-1]!=")")	{
	newick_string_taxa_only_atomized[nstoa] <- ");";
	newick_string_taxa_only_atomized <- c("(",newick_string_taxa_only_atomized);
	}
newick_string_taxa_only <- paste(newick_string_taxa_only_atomized,collapse="");
ancestral <- find_newick_ancestors(newick_string_ancestored=newick_string_taxa_only);
names(ancestral) <- otu_names_nex;
newick_string_taxa_only_raw <- newick_string_taxa_only;
newick_string_taxa_only <- fix_newick_ancestors(newick_string_taxa_only);
vector_tree <- read_newick_string(newick_string_taxa_only);
mat_tree <- transform_vector_tree_to_matrix_tree(vector_tree);
newick_string_taxa_raw <- strsplit(newick_string_taxa_raw,split="")[[1]];
clade_ends <- (1:length(newick_string_taxa_raw))[newick_string_taxa_raw %in% ")"];
colons <- (1:length(newick_string_taxa_raw))[newick_string_taxa_raw %in% ":"];
names(clade_ends) <- length(clade_ends):1;
clade_colons <- 1+clade_ends[(clade_ends+1) %in% colons];
for (cc in 1:length(clade_colons))	{
	n_node <- notu+as.numeric(names(clade_colons)[cc]);
	dd <- clade_colons[cc]+1;
	b_d <- newick_string_taxa_raw[dd];
	dd <- dd+1;
	while (newick_string_taxa_raw[dd] %in% c(".",0:9))	{
		b_d <- paste(b_d,newick_string_taxa_raw[dd],sep="");
		dd <- dd+1;
		}
	branch_durations[n_node] <- as.numeric(b_d);
	}
nNodes <- length(clade_ends);
newick_string_taxa_raw <- paste(newick_string_taxa_raw,collapse="");
if (nNodes<10)	{
	node_names <- paste("node_",1:nNodes,sep="")
	} else if (nNodes<100)	{
	node_names <- c(paste("node_0",1:9,sep=""),
					paste("node_",10:nNodes,sep=""));
	} else	{
	node_names <- c(paste("node_00",1:9,sep=""),
					paste("node_0",10:99,sep=""),
					paste("node_",100:nNodes,sep=""));
	}
names(branch_durations) <- c(otu_names_nex,node_names);

vector_tree_raw <- read_newick_string(newick_string_taxa_only_raw);
#vector_tree <- read_newick_string(newick_string_taxa_only);
molecularize <- strsplit(newick_string_taxa_only,split="")[[1]];

venn_tree_newick <- transform_newick_string_to_venn_tree(newick_string = newick_string_taxa_only);

# have: newick_string_taxa_only,newick_string_taxa_only_raw,vector_tree,ancestral,branch_durations);
# need: clade_posteriors, prob_ancestor,hpd;
# uset mat_tree and newick_string_taxa_only to figure out which node is what number
newick_rem_info <- gsub("sampled_ancestor=","•",newick_string_full);
newick_rem_info <- gsub("age_95%_HPD=","§",newick_rem_info);
newick_rem_info <- gsub("posterior=","¶",newick_rem_info);
hpd <- data.frame(lb=as.numeric(rep(0,notu+nNodes)),ub=as.numeric(rep(0,notu+nNodes)));
molecularized <- str_split(newick_rem_info,"")[[1]];
molecules <- length(molecularized);
panc_boundaries <- (1:molecules)[molecularized=="•"];
post_boundaries <- (1:molecules)[molecularized=="¶"];
hpd_boundaries <- (1:molecules)[molecularized=="§"];
taxon_boundaries <- array(0,dim=c(notu,2));
for (nn in 1:notu)	{
	taxon_dummy <- gsub(otu_names_nex[nn],"£",newick_rem_info);
	taxon_boundaries[nn,1] <- (1:length(str_split(taxon_dummy,"")[[1]]))[str_split(taxon_dummy,"")[[1]]=="£"];
	taxon_boundaries[nn,2] <- taxon_boundaries[nn,1]+length(str_split(otu_names_nex[nn],"")[[1]])-1;
	}
clade_boundaries <- accersi_clade_boundaries_from_newick_string(newick_rem_info);
colnames(taxon_boundaries) <- colnames(clade_boundaries);
tu_boundaries <- rbind(taxon_boundaries,clade_boundaries);
rownames(tu_boundaries) <- all_names <- c(otu_names_nex,node_names);
tu_boundaries <- tu_boundaries[order(tu_boundaries$rp),];
brackets_l <- (1:molecules)[molecularized=="{"];
brackets_r <- (1:molecules)[molecularized=="}"];
for (hp in 1:length(hpd_boundaries))	{
	tx <- sum(tu_boundaries$rp<hpd_boundaries[hp]);
	txn <- match(rownames(tu_boundaries)[tx],all_names);
	i <- brackets_l[1+sum(hpd_boundaries[hp]>brackets_l)]+1;
	j <- brackets_r[1+sum(hpd_boundaries[hp]>brackets_r)]-1;
	hpd[txn,] <- as.numeric(str_split(paste(molecularized[i:j],collapse=""),",")[[1]]);
	}
rownames(hpd) <- all_names;
hpd[vector_tree[(1:notu)[ancestral==1]],] <- hpd[(1:notu)[ancestral==1],];

rownames(taxon_boundaries) <- otu_names_nex;
taxon_boundaries <- taxon_boundaries[order(taxon_boundaries[,1]),];
brackets_l <- (1:molecules)[molecularized=="["];
brackets_r <- (1:molecules)[molecularized=="]"];
prob_ancestor <- array(0,dim=notu);
for (pb in 1:length(panc_boundaries))	{
	tx <- sum(taxon_boundaries[,1]<panc_boundaries[pb]);
	txn <- match(rownames(taxon_boundaries)[tx],otu_names_nex);
	i <- panc_boundaries[pb]+1;
	pranc <- c();
	while (!molecularized[i] %in% c(",","]"))	{
		pranc <- paste(pranc,molecularized[i],sep="");
		i <- i+1;
		}
	prob_ancestor[txn] <- as.numeric(pranc);
	}
names(prob_ancestor) <- otu_names_nex;

clade_posteriors <- array(0,dim=nNodes);
rownames(clade_boundaries) <- names(clade_posteriors) <- node_names;
clade_boundaries <- clade_boundaries[order(clade_boundaries$rp),];
for (pp in 1:length(post_boundaries))	{
	cl <- sum(clade_boundaries$rp<post_boundaries[pp]);
	cld <- match(rownames(clade_boundaries)[cl],node_names);
	i <- post_boundaries[pp]+1;
	postp <- c();
	while (!molecularized[i] %in% c(",","]"))	{
		postp <- paste(postp,molecularized[i],sep="");
		i <- i+1;
		}
	clade_posteriors[cld] <- as.numeric(postp);
#	molecularized[i:i+-1:10]
	}
#sum(str_split(newick_rem_info,"")[[1]]=="•")
#sum(str_split(newick_rem_info,"")[[1]]=="§")
#sum(str_split(newick_rem_info,"")[[1]]=="¶")
output <- list(newick_string_taxa_only,newick_string_taxa_only_raw,vector_tree,clade_posteriors,ancestral,prob_ancestor,hpd,branch_durations);
names(output) <- c("newick_modified","newick","vector_tree","clade_posteriors","ancestral","prob_ancestor","hpd","branch_durations");
return(output);
}

#### MODIFY TREES STORED IN MEMORY ####
transform_matrix_tree_to_vector_tree <- function (matrix_tree)	{
Nnode <- dim(matrix_tree)[1]
notus <- max(matrix_tree)-Nnode
ttus <- max(matrix_tree)
vector_tree <- vector(length=ttus)
vector_tree[notus+1] <- -1
for (n in Nnode:1)
	vector_tree[matrix_tree[n,]] <- n+notus
return(vector_tree)
}

# routine to extract vector tree from matrix giving total progeny of a node 
transform_venn_tree_to_vector_tree <- function (venn_tree)	{
Nnode <- dim(venn_tree)[1]
notus <- dim(venn_tree)[2]
max_otus <- max(venn_tree)
base <- max_otus+1
vtree <- vector(length=(max_otus+Nnode))
otus <- sort(venn_tree[1,])
for (s in 1:notus)	{
	spc <- otus[s]
	vtree[spc] <- max_otus+sort(which(venn_tree==spc,arr.ind=TRUE)[,1],decreasing=TRUE)[1]
	}
vtree[base] <- -1
vtree[base+1] <- base
for (n in 3:Nnode)	{
	htu <- max_otus+n
	lead <- venn_tree[n,1]
	vtree[htu] <- max_otus+sort(which(venn_tree[1:(n-1),]==lead,arr.ind=TRUE)[,1],decreasing=TRUE)[1]
	}
return(vtree)
}

transform_vector_tree_to_matrix_tree <- function(vector_tree)	{
node_rosetta <- sort(unique(vector_tree[vector_tree>0]))
Nnodes <- length(node_rosetta)
maxtomy <- max((hist(vector_tree[vector_tree>1],breaks=((min(vector_tree[vector_tree>1])-1):max(vector_tree[vector_tree>1])),plot=FALSE)$counts))
#order(vector_tree)[2:length(vector_tree)]
node_rich <- vector(length=Nnodes)
matrix_tree <- matrix(0,Nnodes,maxtomy)
for (i in 1:length(vector_tree))	{
	node <- match(vector_tree[i],node_rosetta)
	if(!is.na(node))	{
		node_rich[node] <- node_rich[node]+1
		matrix_tree[node,node_rich[node]] <- i
		}
#	if (vector_tree[i]>=node_rosetta[1])	{
#		node <- match(vector_tree[i],node_rosetta)
#		node_rich[node] <- node_rich[node]+1
#		matrix_tree[node,node_rich[node]] <- i
#		}
	}
return(matrix_tree)
}

transform_vector_tree_to_venn_tree <- function(vector_tree)	{
ohtu <- length(vector_tree);
base <- match(-1,vector_tree);
otu <- base-1
htu <- ohtu-otu
venn_tree <- matrix(0,ohtu,otu)
for (i in 1:otu)	venn_tree[base,i] <- i

node_rich <- vector(length=ohtu)
for (sp in otu:1)	if (vector_tree[sp]!=0)			node_rich[vector_tree[sp]] <- node_rich[vector_tree[sp]]+1
for (nd in ohtu:(base+1))	if (vector_tree[nd]>0)	node_rich[vector_tree[nd]] <- node_rich[vector_tree[nd]]+node_rich[nd]
node_div <- vector(length=ohtu)
for (sp in 1:otu)	{
	node_div[vector_tree[sp]] <- node_div[vector_tree[sp]]+1
	venn_tree[vector_tree[sp],node_div[vector_tree[sp]]] <- sp
	}

for (nd in ohtu:(base+1))	{
	anc <- vector_tree[nd]
	for (i in 1:node_div[nd])	{
		node_div[anc] <- node_div[anc]+1
		venn_tree[anc,node_div[anc]] <- venn_tree[nd,i]
		}
	}
#venn_tree[base:ohtu,1:15]

return(venn_tree[base:ohtu,])
}

create_phylo_class_from_nexus_tree <- function(vector_tree,tip.label)	{
htu1 <- min(vector_tree[vector_tree>0])
Nnode <- 1+(max(vector_tree)-htu1)
otus <- htu1-1
edges <- matrix(0,otus+Nnode-1,2)
j <- 0
for (i in 1:length(vector_tree))	{
	if (vector_tree[i]!=-1)	{
		j <- j+1
		edges[j,1] <- vector_tree[i]
		edges[j,2] <- i
		}
	}
output_tree <- list(edges,tip.label,Nnode)
names(output_tree) <- c("edge","tip.label","Nnode")
class(output_tree) <- "phylo"
#str(output_tree)
return(output_tree)
}

create_phylo_class_from_nexus_tree_file <- function(nexustreefile,tip.label)	{
vector_tree <- read_newick_tree_from_file(nexustreefile)
htu1 <- min(vector_tree[vector_tree>0])
Nnode <- 1+(max(vector_tree)-htu1)
otus <- htu1-1
edges <- matrix(0,otus+Nnode-1,2)
j <- 0
for (i in 1:length(vector_tree))	{
	if (vector_tree[i]!=-1)	{
		j <- j+1
		edges[j,1] <- vector_tree[i]
		edges[j,2] <- i
		}
	}
output_tree <- list(edges,tip.label,Nnode)
names(output_tree) <- c("edge","tip.label","Nnode")
class(output_tree) <- "phylo"
#str(output_tree)
return(output_tree)
}