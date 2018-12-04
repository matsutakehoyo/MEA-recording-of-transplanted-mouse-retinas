#analyse mERG peaks


#time to peak bins to detect a, b, c, d peaks
#awave is negative peak between 0 and 50
#bwave is positive peak between 50 and 150
#cwave is negative peak between 150 and 450
#dwave is positive peak between 450 and 150


analyse_peaks <- function(data, peak_w = 100, peak_span = 0.005, 
                          a_delay = 55, b_delay = 120, c_delay = 350, d_delay = 1000, logfile=logfile){

}


#this finds max and min and return a list with their values and index
findpeaks <- function(x, y, w, ...) {
	require(zoo)
	n <- length(y)
	y_smooth <- loess(y ~ x, ...)$fitted
	sd <- sd(y - y_smooth)
	y_max <- rollapply(zoo(y_smooth), 2*w+1, max, align="center")
	y_min <- rollapply(zoo(y_smooth), 2*w+1, min, align="center")
	delta_max <- y_max - y_smooth[-c(1:w, n+1-1:w)]
	i_max <- which(delta_max <= 0) + w
	delta_min <- y_min - y_smooth[-c(1:w, n+1-1:w)]
	i_min <- which(delta_min >= 0) + w
	list(sd = sd, 
		 x_max = x[i_max], 
		 x_min = x[i_min], 
		 i_max = i_max, 
		 i_min = i_min, 
		 y_hat = y_smooth)
}

test_findpeaks <- function(x, y, w, span) {
	peaks <- findpeaks(x, y, w=w, span=span)  
	plot(x, y, cex=0.75, col="gray", type="l", main=paste("w = ", w, ", span = ", span, sep=""))
	lines(x, peaks$y_hat,  lwd=2) #$
	y_min <- min(y)
	y_max <- max(y)
	#peaks
	sapply(peaks$i_max, function(i) lines(c(x[i],x[i]), c(y_min, peaks$y_hat[i]),col="Red", lty=2))
	points(x[peaks$i_max], peaks$y_hat[peaks$i_max], col="Red", pch=19, cex=1.25)
	#valleys
	sapply(peaks$i_min, function(i) lines(c(x[i],x[i]), c(y_max, peaks$y_hat[i]),col="Blue", lty=2))
	points(x[peaks$i_min], peaks$y_hat[peaks$i_min], col="Blue", pch=19, cex=1.25)
}

#this function performs baseline analysis and peak detection
#analyse_peaks(x=filter(erg, channel=="CH01")$Time, y=filter(erg, channel=="CH01")$mV, w=300, span=0.005)
analyse_peaks <- function(x, y, flash_time, plot=TRUE, w, span, ...){
	y_min <- min(y)
	y_max <- max(y)
	flash_index <- which(x==flash_time)
	y_peaks <- findpeaks(x = x, y = y, w = w, span = span)
	#b/d peaks
	signal_pos <- which( (x[y_peaks$i_max] > flash_time + a_delay - 20) & (x[y_peaks$i_max] < flash_time + b_delay) |
						 (x[y_peaks$i_max] > flash_time + c_delay/2) & (x[y_peaks$i_max] < flash_time + d_delay) )
	#a/c peaks
	signal_neg <- which( (x[y_peaks$i_min] > flash_time + 10) & (x[y_peaks$i_min] < flash_time + a_delay) |
						 (x[y_peaks$i_min] > flash_time + a_delay/2) & (x[y_peaks$i_min] < flash_time + c_delay) )
	#signal_pos <- which(y_peaks$y_hat[y_peaks$i_max] > sd_tolerance * y_peaks$sd)
	#signal_neg <- which(abs(-y_peaks$y_hat[y_peaks$i_min]) > sd_tolerance * y_peaks$sd)
	
	if (plot==TRUE){
		plot(x=x, y=y, xlab="time(msec)", ylab="mV", cex=0.75, col="gray", type="l")  #original data
		lines(x=x, y=y_peaks$y_hat, lwd=2)      #smoothed trace
		#peaks
		points(x[y_peaks$i_max], y_peaks$y_hat[y_peaks$i_max], col="hotpink", pch=19, cex=0.6)
		points(x[y_peaks$i_max[signal_pos]], y_peaks$y_hat[y_peaks$i_max[signal_pos]], col="Red", pch=19, cex=1.25)
		sapply(y_peaks$i_max[signal_pos], function(i) lines(c(x[i],x[i]), c(y_min, y_peaks$y_hat[i]),col="Red", lty=2))
		#valleys
		points(x[y_peaks$i_min], y_peaks$y_hat[y_peaks$i_min], col="royalblue", pch=19, cex=0.6)
		points(x[y_peaks$i_min[signal_neg]], y_peaks$y_hat[y_peaks$i_min[signal_neg]], col="Blue", pch=19, cex=1.25)
		sapply(y_peaks$i_min[signal_neg], function(i) lines(c(x[i],x[i]), c(y_peaks$y_hat[i], y_max),col="Blue", lty=2))
		#light stimulus
		lines(c(x[flash_index],x[flash_index]), c(y_min, y_max),col="Gold", lwd=2)
	}

	list(original       = y, 
		 sd             = y_peaks$sd, 
		 x_max          = y_peaks$x_max,   #x value for local maxima
		 x_min          = y_peaks$x_min,   #x value for local minima
		 i_max          = y_peaks$i_max,   #index for local maxima
		 i_min          = y_peaks$i_min,   #index for local minima
		 y_hat          = y_peaks$y_hat,   #smoothed trace
		 signal_pos     = signal_pos, #i_max index for positive signals 
		 signal_neg     = signal_neg #i_min index for negative signals
		 )
}

plot_peaks <- function(x, peaks, event_pos){
	flash_index = which(x==event_pos)
	df <- tibble(frametime=x, original=peaks$original, smooth=peaks$y_hat)
	p <- ggplot(df) + 
		geom_line(aes(x=frametime, y=original), color="darkgray") +
		geom_line(aes(x=frametime, y=smooth), color="black") +
		geom_vline(data=tibble(x=x[flash_index]), aes(xintercept=x), color="gold") +
		xlab("time(s)") + ylab("micro Volts") +
		basic_style

	#max
	if (length(peaks$i_max)){
		p <- p +
			geom_point(data=tibble(x=peaks$x_max, y=peaks$y_hat[peaks$i_max]), aes(x=x, y=y), color="pink", alpha=1/2)
	
		if (length(peaks$signal_pos)){
			p <- p +
				geom_point(data=tibble(x=peaks$x_max[peaks$signal_pos], 
									   y=peaks$y_hat[peaks$i_max[peaks$signal_pos]]), 
							aes(x=x, y=y), color="red", alpha=1/2) +
				geom_segment(data=tibble(x=peaks$x_max[peaks$signal_pos], 
										 y=min(peaks$original), 
										 xend=peaks$x_max[peaks$signal_pos], 
										 yend=peaks$y_hat[peaks$i_max[peaks$signal_pos]]), 
							aes(x=x, y=y, xend=xend, yend=yend), col="red", linetype="dotted", size=0.5) 
		}
	}
	#min
	if (length(peaks$i_min)){
		p <- p +
			geom_point(data=tibble(x=peaks$x_min, y=peaks$y_hat[peaks$i_min]), aes(x=x, y=y), color="lightblue", alpha=1/2) 
	  
		if (length(peaks$signal_neg)){
			p <- p +
				geom_point(data = tibble(x = peaks$x_min[peaks$signal_neg], 
										 y = peaks$y_hat[peaks$i_min[peaks$signal_neg]]), 
							aes(x=x, y=y), color="blue", alpha=1/2) +
				geom_segment(data=tibble(x=peaks$x_min[peaks$signal_neg], 
										 y=max(peaks$original), 
										 xend=peaks$x_min[peaks$signal_neg], 
										 yend=peaks$y_hat[peaks$i_min[peaks$signal_neg]]), 
							aes(x=x, y=y, xend=xend, yend=yend), col="blue", linetype="dotted", size=0.5)
		}
	}
	p
}

