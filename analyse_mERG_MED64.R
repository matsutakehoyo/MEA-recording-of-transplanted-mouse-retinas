#analyse mERG data
library(tidyverse)
source("take.R")
source("analyse_merg_peaks.R")
library(signal)			#for signal *dplyr::filter


peak_w = 10
peak_span = 0.02
a_delay = 55
b_delay = 120
c_delay = 350
d_delay = 1000

split_path <- function(path) {
  if (dirname(path) %in% c(".", path))
	return(basename(path))
  return(c(basename(path), split_path(dirname(path))))
}

#mea filter 
pass_filter <- butter(n = 2, W = c(1/10000, 50/10000), type = "pass") #1-50Hz band pass filter for mERG


logfile <- gen_start_log("Import MEA h5 data to R and analyse mERG")           #log results in file
#files to process
#MED recordings were converted to csv files that end in -3digits.csv, 
files <- list.files(pattern = "\\.csv", recursive=TRUE)
#to avoid csv created by previous execution
files <- files[!grepl(pattern="mERG_peaks.csv", files)]

gen_log(paste("Files to process:"), logfile)
gen_log(capture.output(print(files)), logfile)


for (i in seq_along(files)){
	# i=1
	print(files[i])
	# things are saved here
	# filename <- gsub(pattern=".csv", replacement="", files[i])
	filename <- gsub(pattern=".csv", replacement="", files[i])

	# read data
	mea_data <- read_csv(files[i], skip=3, col_names=FALSE, cols(.default = col_double())) 
	colnames(mea_data) <- c("time",  seq(1:64))
	
	mea_data <- mea_data %>% 
		gather(key=channel, value=value, -time) %>% #long format
		mutate(channel = as.integer(gsub("\\D", "", channel))) %>%
		mutate(col = ((channel-1) %% 8 +1) %>% as.integer(),
		       row = ((channel-1) %/% 8 + 1) %>% as.integer()) %>%
		mutate(light_on = 2500) %>%
		add_column(filename = filename, .before=1) %>%
		group_by(channel) %>%
		mutate(filter = signal::filter(pass_filter, value)) %>%
		ungroup()

	gen_log("          mERG data loaded", logfile, datetime=TRUE)
	# add filter for mERG analsysis 1-50Hz

	# glimpse(mea_data)
	flash_time <- mea_data$light_on %>% unique()
	
	gen_log("Calculating peaks for each channel", logfile, datetime=TRUE)
	#save plots of individual channels here
	
	dir.create(file.path(filename, "channels"), showWarnings = FALSE, recursive=TRUE)
	dir.create(path=(file.path(filename, "channels_fit")), showWarnings = FALSE, recursive=TRUE)
	erg_peaks <- tibble(file = character(),
						channel = integer(), 
						row = integer(), 
						col = integer(),
						distance = double(),
						a_abs = double(), 
						b_abs = double(),
						a_amplitude = double(), 
						b_amplitude = double(),
						a_time = double(),
						b_time = double(),
						base = double(),
						noise = double())
	
	pb_ch <- progress_estimated(64)
	for (ch in seq(1:64)){
		pb_ch$tick()$print()
		erg_ch <- mea_data %>% 
			dplyr::filter(channel==ch) %>% 
			dplyr::filter(time > (flash_time - 500) & time < (flash_time + 1500))
			
		if(nrow(erg_ch)){
			base_ch <- erg_ch %>% 
				dplyr::filter(time < flash_time) %>%
				select(value) %>% pull() %>% mean()

			noise_ch <- erg_ch %>% 
				dplyr::filter(time < flash_time) %>%
				select(value) %>% pull() %>% sd()

			# test_findpeaks(x=erg_ch$time, y=erg_ch$filter, w=10, span=0.02)
			erg_ch_peak <- analyse_peaks(x=erg_ch$time, y=erg_ch$filter, flash_time=flash_time,
			                             plot=FALSE, w=peak_w, span=peak_span)
			p <- plot_peaks(x=erg_ch$time, peak=erg_ch_peak, event_pos=flash_time) + 
			scale_x_continuous(limits = c(flash_time - 500, flash_time + 1000), expand = c(0, 0))
			ggsave(paste0(filename, "/channels_fit/", sprintf("CH%02d", ch), "_ab.png"), p, width=3, height=2)

			#are detected peaks are in reasonable location?
			peak_times_neg <- erg_ch$time[erg_ch_peak$i_min[erg_ch_peak$signal_neg]]-flash_time
			peak_times_pos <- erg_ch$time[erg_ch_peak$i_max[erg_ch_peak$signal_pos]]-flash_time
		
			condition_a <- which( 0 < peak_times_neg & peak_times_neg < a_delay)
			condition_b <- which(a_delay - 25 < peak_times_pos & peak_times_pos < b_delay)

		
			if (length(condition_a)){			
				aabs <- min(erg_ch_peak$y_hat[which(erg_ch$time %in% (peak_times_neg[condition_a] + flash_time))])
				aamplitude <- aabs - base_ch
				atime <- peak_times_neg[peak_times_neg %in% (erg_ch$time[which(erg_ch_peak$y_hat==aabs)] - flash_time)]
			} else {
				aabs <- NA
				aamplitude <- NA
				atime <- NA
			}
	
			if (length(condition_b)){
				babs <- max(erg_ch_peak$y_hat[which(erg_ch$time %in% (peak_times_pos[condition_b] + flash_time))])
				if (is.na(aabs)){
					if (babs > 0){
						bamplitude <- babs	
					} else {
						bamplitude <- NA
					}		
				} else if (aabs < 0){
					bamplitude <- sum(babs, -aamplitude, na.rm=TRUE) #b wave is difference from a wave not from 0
				} else {
					bamplitude <- NA
				}
				
				btime <- peak_times_pos[peak_times_pos %in% (erg_ch$time[which(erg_ch_peak$y_hat==babs)] - flash_time)]
				
			} else {
				babs <- NA
				bamplitude <- NA
				btime <- NA
			}
		
	
			# dist <- sqrt((col_coord(cent_channel) - col_coord(ch))^2 + 
			# 			 (row_coord(cent_channel) - row_coord(ch))^2) * distance_factor
			dist <- NA	

		
			erg_peaks <- bind_rows(erg_peaks,
							tibble(
									file = filename,
								   channel = ch, 
								   row = erg_ch$row %>% unique(),
								   col = erg_ch$col %>% unique(),
								   distance = dist,
								   a_abs = aabs, 
								   b_abs = babs, 
								   a_amplitude = aamplitude, 
								   b_amplitude = bamplitude,
								   a_time = atime,
								   b_time = btime,
								   base = base_ch, 
								   noise = noise_ch))
		}#if(nrow(erg_ch)){
	}#for (channel in seq(1:64))
	#save file
	write.csv(erg_peaks, file = paste0(filename, "/mERG_peaks.csv"), row.names = FALSE)

	erg_peaks <- erg_peaks %>% 
		gather(key=key, value=value, a_abs:b_time) %>%
		separate(key, into=c("peak", "key")) %>%
		arrange(channel, peak, key) %>%
		spread(key, value)
	
	erg_peaks_summary <- erg_peaks %>%
		group_by(peak) %>% 
		summarise(amp_mean=mean(amplitude, na.rm = TRUE), 
		          abs_mean=mean(abs, na.rm = TRUE),
		          amp_sd=sd(amplitude, na.rm = TRUE), 
		          time_mean=mean(time, na.rm = TRUE), 
		          time_sd=sd(time, na.rm = TRUE))
	
	p_erg <- ggplot(mea_data, aes(x = time - flash_time, y = filter)) +
		geom_line(size = 0.2) +
		geom_point(data = erg_peaks, aes(x = time, y = abs, color = peak), size=2, alpha=1/2) +
		facet_grid(row~col) +
		xlab("time after flash(ms)") + ylab("response (uV)") + ggtitle("mERG") +
		scale_x_continuous(limits = c(-500, 1000), expand = c(0, 0)) +
		basic_style
	ggsave(paste0(filename, "/mERG overview.png"), width = 24, height = 18)

	p_average <- ggplot(data=mea_data, aes(x = time - flash_time, y = filter)) +
		stat_summary(fun.y = mean, geom = "line", color = "red") +
		geom_line(size = 0.1, color = "gray", alpha = 1/10) +
		geom_point(data = erg_peaks, aes(x = time, y = abs, color = peak)) +
		geom_point(data = erg_peaks_summary, aes(x = time_mean, y = abs_mean), color = "red", size = 3) +
		scale_x_continuous(limits = c(-500, 1000), expand = c(0, 0)) +
		xlab("time after flash(ms)") +
		ylab("response (uV)") +
		basic_style
	ggsave(paste0(filename, "/mERG average.png"), width = 12, height = 9)

	pb_ch <- progress_estimated(64)
	for (ch in seq(1:64)){	
		pb_ch$tick()$print()
		p_ch <-  ggplot(mea_data %>% dplyr::filter(channel==ch), aes(x = time - flash_time, y = filter)) +
			geom_line(size = 0.5, alpha = 1/2) +
			geom_point(data = erg_peaks %>% dplyr::filter(channel==ch), aes(x = time, y = abs, color = peak), size=3, alpha = 1/2) +
			geom_vline(aes(xintercept=0), color="gold") +
			scale_x_continuous(limits = c(-500, 1000), expand = c(0, 0)) +
			xlab("time after flash(ms)") +
			ylab("response (uV)") +
			ggtitle(paste0("CH:", ch)) +
			basic_style
		ggsave(paste0(filename, "/channels/", sprintf("CH%02d", ch), "_peaks.pdf"), p_ch, width=3, height=2)
	}

	p_summary_amplitude <- ggplot(erg_peaks, aes(x=peak, y=amplitude)) +
		# geom_boxplot(outlier.shape = NA) +
		geom_violin() +
		geom_jitter(width = 0.2, color="gray") +
		xlab("peaks") +
		ylab("Amplitude(mV)") +
		ggtitle("Micro ERG peak amplitudes") +
		basic_style
	ggsave(paste0(filename, "/mERG amplitude summary.pdf"))

	p_summary_timetopeak <- ggplot(erg_peaks, aes(x=peak, y=time)) +
		# geom_boxplot(outlier.shape = NA) +
		geom_violin() +
		geom_jitter(width = 0.2, color="gray") +
		xlab("peaks") +
		ylab("time to peak(ms)") +
		ggtitle("Micro ERG peak delay") +
		basic_style
	ggsave(paste0(filename, "/mERG time to peak sumamry.pdf"))

	#heat map 
	h <- ggplot(data = erg_peaks, aes(x=as.factor(col), y=as.factor(row), fill=amplitude)) + 
			geom_tile(color = "white")+
			geom_tile()+
			facet_wrap(~peak)+
			#scale_y_reverse()+
			scale_y_discrete(limits = rev(levels(factor(erg_peaks$col))))+
			scale_x_discrete(position = "top")+
			scale_fill_gradient2(low = "steelblue", high = "firebrick", mid = "white",            #high="coral"
													midpoint = 0, space = "Lab")+    
			coord_fixed()+
			basic_style+
			#theme_minimal()
			theme(
				axis.title.x = element_blank(),
				axis.title.y = element_blank(),
				panel.grid.major = element_blank(),
				panel.border = element_blank(),
				panel.background = element_blank(),
				axis.ticks = element_blank())
	ggsave(paste0(filename, "/mERG heatmap.pdf"), width=8, height=8)
	
	h_value <- h + geom_text(aes(as.factor(col), as.factor(row), label = sprintf("%.3f", amplitude)), color = "black", size = 2)
	ggsave(paste0(filename, "/mERG heatmap values.pdf"), width=8, height=8)

} #for (i in seq_along(files)){
