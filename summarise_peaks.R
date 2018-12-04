library(tidyverse)
library(readxl)
source("take.R")

split_path <- function(path) {
  if (dirname(path) %in% c(".", path))
	return(basename(path))
  return(c(basename(path), split_path(dirname(path))))
}

files <- list.files(pattern = "mERG_peaks.csv", recursive=TRUE)
data <- list()
for (i in seq_along(files)){
	# i=1
	print(files[i])

	data[[i]] <- read_csv(files[i]) %>%
		add_column(mouse = split_path(files[i])[3], 
		           graft = split_path(files[i])[4]) %>%
		separate(graft, into=c("graft", "after_transplant"), sep="_")
}
data <- bind_rows(data)

data <- data %>% 
	separate(file, into=c("mouse", "stimuli", "condition"), sep="_") %>%
	separate(mouse, into=c("mouse", "eye"), sep=-1)
	


data <- data %>% 
	mutate(condition=toupper(condition)) %>% #force upper case to avoid errors
	mutate(condition=factor(condition, levels = c("PRE", "AMES", "L-AP4", "WASHOUT"))) %>%
	mutate(condition = forcats::fct_recode(condition, 
                                         "Pre" = "PRE", 
                                         "AMES" = "AMES", 
                                         "L-AP4" = "L-AP4",
                                         "Washout" = "WASHOUT")) %>%
	mutate(stimuli=factor(stimuli, 
	       levels =c("allND2mAmERG", "allND10mAmERG"), 
	       # labels=c("10.56 log photons/cm2/s", "12.84 log photons/cm2/s"))) %>%
	       labels=c("weak", "strong"))) %>%
	mutate(graft = factor(graft, levels=c("WT", "B4", "ISL1"), labels=c("wt", "Bhlhb4", "Islet1")),
	       after_transplant = factor(after_transplant, levels=c("5W", "8W", "12W")))

data %>% distinct(condition) %>% pull()
data %>% distinct(stimuli) %>% pull()
data %>% distinct(graft) %>% pull()
data %>% distinct(after_transplant) %>% pull()

data <- data %>% add_column(graft_covered = NA)
ch_file <- list.files(pattern="[Ee]lectrodes")
if (length(ch_file)){
	ch_data <- read_excel(ch_file) %>%
	separate(`Mouse ID`, into = c("mouse", "eye"), sep=-1)

	for (m in seq_along(ch_data$mouse)){
		mouse_id <- ch_data$mouse[m] 
		eye_id <-  ch_data$eye[m]
		ch_id <- ch_data[m,] %>% select(X__2:X__65) %>%  as.integer()
		
		data <- data %>%	
			mutate(graft_covered = ifelse(mouse==mouse_id & eye==eye_id, channel %in% ch_id, graft_covered)) 
	}
}
data

saveRDS(data, "merg_data.rds")

