library(grid)
library(gridExtra)
library(wesanderson)
library(scales)

#graph styles
basic_style = theme(plot.title = element_text(lineheight=.8, face="bold"),
					plot.background = element_rect(fill = "transparent", colour = NA),
					legend.background = element_rect(fill = "transparent", colour = NA),
					legend.title=element_blank(), 
					strip.text = element_text(face = "bold"),
					axis.title = element_text(face="bold"))

# for black background
# update_geom_defaults("text",  list(colour = "black"))


#Colorblind-friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

keynote = theme(plot.title = element_text(size=20),
				axis.title = element_text(size=18, face="bold"),
				axis.text = element_text(size=18, face="bold"),
				strip.text = element_text(size=18, face="bold"),
				legend.text=element_text(size=18))

dark_background =  theme(plot.title = element_text(color="white"),
						 axis.title = element_text(color="white"), 
						 axis.text = element_text(color="white"), 
						 strip.text = element_text(color="white"),
						 strip.background = element_rect(colour = "white"), 
						 axis.ticks = element_line(color="white"))

base_size = 12
base_family = ""
dark_background <- theme(
	# Specify axis options
	axis.line = element_blank(),  
	axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
	axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
	axis.ticks = element_line(color = "white", size  =  0.2),  
	axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),
	axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),
	axis.ticks.length = unit(0.3, "lines"),   
	# Specify legend options
	legend.background = element_rect(color = NA, fill = "black"),  
	legend.key = element_rect(color = "white",  fill = "black"),  
	legend.key.size = unit(1.2, "lines"),  
	legend.key.height = NULL,  
	legend.key.width = NULL,      
	legend.text = element_text(size = base_size*0.8, color = "white"),  
	legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
	legend.position = "right",  
	legend.text.align = NULL,  
	legend.title.align = NULL,  
	legend.direction = "vertical",  
	legend.box = NULL, 
	# Specify panel options
	panel.background = element_rect(fill = "black", color  =  NA),  
	panel.border = element_rect(fill = NA, color = "white"),  
	panel.grid.major = element_line(color = "grey35"),  
	panel.grid.minor = element_line(color = "grey20"),  
	panel.spacing = unit(0.5, "lines"),   
	# Specify facetting options
	strip.background = element_rect(fill = "grey30", color = "grey10"),  
	strip.text.x = element_text(size = base_size*0.8, color = "white"),  
	strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
	# Specify plot options
	plot.background = element_rect(color = "black", fill = "black"),  
	plot.title = element_text(size = base_size*1.2, color = "white"),  
	plot.margin = unit(rep(1, 4), "lines")
	
)

transparent = theme(
	panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
	panel.grid.minor = element_line(colour = "gray90"), 
	panel.grid.major = element_line(colour = "gray80"),
	strip.background = element_rect(colour = "black", fill="transparent"),
	plot.background = element_rect(fill = "transparent",colour = NA),
	legend.background= element_rect(fill = "transparent",colour = NA))

# function to squish axis when plotting, 
# use it inside scale_y_continuous(trans = squish_trans(-2, 2, 4))
squish_trans <- function(from, to, factor) {

  trans <- function(x) {

    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to

    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)

    return(x)
  }

  inv <- function(x) {

    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor

    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))

    return(x)
  }

  # return the transformation
  return(trans_new("squished", trans, inv))
}

#log script
gen_start_log<-function(name="Debug"){
  dir.create("log", recursive=TRUE, showWarnings=FALSE) # make the folders before saving data
  local_log_file<-paste(format(Sys.time(),"%Y%m%d%H%M%S"), name, "log", sep=".")
  return(local_log_file)
}

gen_log<-function(text, log_file, datetime=FALSE){
  dir.create(file.path(".","log"), recursive=TRUE, showWarnings=FALSE)
  if (datetime) {
    write(paste(Sys.time(), text, sep=": "), file=file.path(".","log",log_file), append=TRUE)
    print(paste(Sys.time(), text, sep=": "))
  }
  else {
    write(text, file=file.path(".","log",log_file), append=TRUE)
    print(text)
  }
}

#calculate bw for histogram 
hist_bw <- function(x){
	bw <- hist(x, breaks = 'Scott', plot = FALSE)
	bw <- bw$mids %>% diff() %>% unique()
	return(bw)
}

#calculate mode of density distribution
dens_mode <- function(x){
	d <- density(x, adjust=0.1) 
	d$x[which.max(d$y)]
}

# return shape & rate parameters from mean and sd (from DBDA
# from Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition:)
gamma_shape_rate = function( mean , sd ) {
	if ( mean <=0 ) stop("mean must be > 0")
	if ( sd <=0 ) stop("sd must be > 0")
	shape = mean^2/sd^2
	rate = mean/sd^2
	return( list( shape=shape , rate=rate ) )
}

# return shape & rate parameters from mean and sd
# from Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition:)
gamma_shape_rate2 = function( mode , sd ) {
	if ( mode <=0 ) stop("mode must be > 0")
	if ( sd <=0 ) stop("sd must be > 0")
	rate = ( mode + sqrt( mode^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )
	shape = 1 + mode * rate
	return( list( shape=shape , rate=rate ) )
}

#pairs/scatter plot
#psecify color_fill as string
take_pairs <- function(df, color_fill=NA){
	require(GGally)
	ggp <- df %>%
		{if(is.na(color_fill)) 
			ggpairs(., diag = 'blank', lower = 'blank') 
		else ggpairs(., diag = 'blank', lower = 'blank', mapping = ggplot2::aes_string(color = color_fill))}
	print("customizing pairwise comparison plots: diagonal")
	#ggpair diagonal 
	for(i in seq_along(df)) {
		x <- df[,i] %>% unlist() #unlist to untibble it
		if(!is.na(color_fill)){
			p <- ggplot(eval(parse(text = paste0("data.frame(x, group=df$", color_fill, ")"))), 
						aes(x, fill=group, color=group))
		} else {
			p <- ggplot(data.frame(x), aes(x))
		}
		p <- p + theme(text=element_text(size=14), axis.text.x=element_text(angle=40, vjust=1, hjust=1))
		if (class(x) == 'factor') {
			p <- p + geom_bar(color='grey20')
		} else {
			p <- p + geom_histogram(binwidth=hist_bw(x), color='grey20')
			p <- p + geom_line(eval(bquote(aes(y=..count..*.(hist_bw(x))))), stat='density')
		}
		# p <- p + 
		# 	geom_label(data=data.frame(x=-Inf, y=Inf, label=colnames(df)[i]), 
		# 				aes(x=x, y=y, label=label), hjust=0, vjust=1)
		p <- p + scale_fill_manual(values = cbPalette)
		ggp <- putPlot(ggp, p, i, i) 
	}
	
	#upper panels change the cor display
	zcolat <- seq(-1, 1, length=81)
	zcolre <- c(zcolat[1:40]+1, rev(zcolat[41:81]))
	print("customizing pairwise comparison plots: upper")
	for(i in 1:(ncol(df)-1)) {
		for(j in (i+1):ncol(df)) {
			x <- df[,i] %>% pull() %>% as.numeric()  #because cor requires numeric
			y <- df[,j] %>% pull() #%>% as.numeric()
			if (class(y) != 'factor' ){
				y <- as.numeric(y)
				r <- cor(x, y, method='spearman', use='pairwise.complete.obs')
				zcol <- lattice::level.colors(r, at=zcolat,
											  col.regions=colorRampPalette(c(scales::muted('red'), 'white', scales::muted('blue')), 
											  							 space='rgb')(81))
				textcol <- ifelse(abs(r) < 0.4, 'grey20', 'white')
				ell <- ellipse::ellipse(r, level=0.95, type='l', npoints=50, scale=c(.2, .2), centre=c(.5, .5))
				p <- ggplot(data.frame(ell), aes(x=x, y=y))
				p <- p + theme_bw() + 
					theme(
						plot.background=element_blank(),
						panel.grid.major=element_blank(), 
						panel.grid.minor=element_blank(),
						panel.border=element_blank(), axis.ticks=element_blank()
					)
				p <- p + geom_polygon(fill=zcol, color=zcol)
				p <- p + geom_text(data=NULL, x=.5, y=.5, label=100*round(r, 2), size=6, col=textcol)
				ggp <- putPlot(ggp, p, i, j)	
			} else {
				ggp[i, j] <- ggp[i, j] + scale_fill_manual(values = cbPalette)				
			}
		}
	}
	
	#lower panels
	print("customizing pairwise comparison plots: lower")
	for(j in 1:(ncol(df)-1)) {
		for(i in (j+1):ncol(df)) {
			x <- df[,j] %>% pull() 
			y <- df[,i] %>% pull() 
			if (class(y) == 'factor' ) {
				if(!is.na(color_fill)){
					p <- ggplot(eval(parse(text = paste0("data.frame(x, y, group=df$", color_fill, ")"))), 
								aes(x, fill=group, color=group)) +
						scale_fill_manual(values = cbPalette) +
						scale_color_manual(values = cbPalette)
				} else {
					p <- ggplot(data.frame(x, y), aes(x))
				}
				p <- p + theme(text=element_text(size=14), axis.text.x=element_text(angle=40, vjust=1, hjust=1))
				if (class(x) == 'factor') {
					p <- p + geom_col(aes(y=y)) + facet_grid(y~.)
				} else {
					p <- p + geom_histogram(binwidth=hist_bw(x)) + facet_grid(y~.)
				}
			} else {
				if(!is.na(color_fill)){
					p <- ggplot(eval(parse(text = paste0("data.frame(x, y, group=df$", color_fill, ")"))), 
								aes(x, y, fill=group, color=group)) +
						scale_fill_manual(values = cbPalette) +
						scale_color_manual(values = cbPalette)
				} else {
					p <- ggplot(data.frame(x, y), aes(x=x, y=y)) 
				}
				
				p <- p + theme(text=element_text(size=14), axis.text.x=element_text(angle=40, vjust=1, hjust=1))
				if (class(x) == 'factor') {
					p <- p + geom_boxplot(alpha=3/6, outlier.size=0, fill='white')
					p <- p + geom_point(position=position_jitter(w=0.4, h=0), size=1)
				} else {
					p <- p + geom_point(size=1)
				}
			}
			ggp <- putPlot(ggp, p, i, j)
		}
	}
	ggp
}


