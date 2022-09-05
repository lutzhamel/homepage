### map-utils.r 
##################################################################################################
# NOTE: this script is no longer supported, it has been superseded by the package popsom in CRAN # 
##################################################################################################
# version 2.2
# (c) 2009-2012 Lutz Hamel, University of Rhode Island
#
# This file constitues a set of routines which are useful in constructing
# and evaluating self-organizing maps (SOMs).  The utilities are built around
# the 'som' library available from the CRAN archive (see the package installer
# and package manager).  The main utilities available in this file are:
#	map.build - constructs a map
#	map.convergence - reports the accuracy of the map in terms of modeling the 
#                     underlying data distribution (100% if all feature distributions
#                     are modeled correctly, 0% if none are)
#	map.significance - graphically reports the significance of each feature with
#                      respect to the self-organizing map model
#	map.umat - displays the "unified distance matrix" (umat) of the SOM model (lighter
#              represent strong clusters, dark colors represent weak clusters)
#	map.starburst - displays the starburst representation of the SOM model, the centers of
#                   starbursts are the centers of clusters
#	map.projection - print a table with the associations of labels with map elements
#   map.feature - compute and display the enhanced unified distance matrix for a 
#                 feature of the training data
#
### Papers that document the theoretical aspects of this software package:
#
# "Bayesian Probability Approach to Feature Significance for Infrared 
# Spectra of Bacteria", Lutz Hamel, Chris W. Brown, Applied Spectroscopy, Volume 66, Number 1, 2012.
#
# "A Population Based Convergence Criterion for Self-Organizing Maps", Lutz Hamel and Benjamin Ott. 
# Proceeding of the 8th International Conference on Data Mining (DMIN'12), to appear.
#
# "Improved Interpretability of the Unified Distance Matrix with Connected Components", Lutz Hamel and 
# Chris W. Brown. Proceeding of the 7th International Conference on Data Mining (DMIN'11), July 18-21, 2011, 
# Las Vegas Nevada, USA, ISBN: 1-60132-168-6, pp338-343, CSREA Press, 2011.
#
# (preprints of these papers are available at www.cs.uri.edu/~hamel)
#
### License
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
### Example
# source("map-utils.r")
# data(iris)
# df <- subset(iris,select=-Species)
# labels <- subset(iris,select=Species)
# m <- map.build(df,labels,xdim=15,ydim=10,train=1000)
# map.convergence(m)
# map.significance(m)
# map.starburst(m)

# load libraries
require(som)
require(fields)
require(graphics)

### map.build -- construct a SOM, returns an object of class 'map'
# parameters:
# - data - a dataframe where each row contains an unlabeled training instance
# - labels - a vector or dataframe with one label for each observation in data
# - xdim,ydim - the dimensions of the map
# - alpha - the learning rate, should be a positive non-zero real number
# - train - number of training iterations
# retuns:
# - an object of type 'map'

# Hint: if your training data does not have any labels you can construct
#       simple label vector as follows: labels <- 1:ncol(training.data)

map.build <- function(data,labels,xdim=10,ydim=5,alpha=.6,train=1000) {

	# compute the initial neighborhood radius
	r <- sqrt(xdim^2 + ydim^2)

	# train the map
	m <- som(data,
			 xdim=xdim,
			 ydim=ydim,
			 init="random",
			 alpha=c(alpha,alpha),
			 alphaType="linear",
			 neigh="bubble",
			 topol="rect",
			 radius=c(r,r),
			 rlen=c(1,train))

	# we add one more field to the list we get back from som
	# namely the labels as a data frame.
	m$labels <- data.frame(labels)
	
	# the new list is now an object of type 'map'
	class(m) <- "map"
	
	# return the map
	m
}


### map.convergence - evaluate the convergence of a map using the F-test and 
#                     a Bayesian estimate of the variance in the training data.
# parameters:
# - map is an object if type 'map'
# - conf.int is the confidence interval of the convergence test (default 95%)
#
# - return value is the convergence index of the map (variance captured by the map so far)

# Hint: the convergence index is the variance of the trainig data captured by the map;
#       maps with convergence of less than 90% are typically not trustworthy.  Of course,
#       the precise cut-off depends on the noise level in your training data. 

map.convergence <- function(map,conf.int=.95) {

	 if (class(map) != "map")
		stop("map.convergence: first argument is not a map object.")
	 
	 # map.df is a dataframe that contains the code vectors
	 # note: map$code is what the 'som' package produces
	 map.df <- data.frame(map$code)

	 # data.df is a dataframe that contain the training data
	 # note: map$data is what the 'som' package returns
	 data.df <- data.frame(map$data)

	 # do the F-test on a pair of datasets: code vectors/training data
	 l <- df.var.test(map.df,data.df,conf=conf.int)

	 # compute the variance captured by the map
	 nfeatures <- ncol(map.df)
	 prob.v <- map.significance(map,graphics=FALSE)
	 var.sum <- 0
	 for (i in 1:nfeatures) {
	 	if (l$conf.int.lo[i] <= 1.0 && l$conf.int.hi[i] >= 1.0) {
				var.sum <- var.sum + prob.v[i]
		#cat("Feature",i,":\t",as.real(l$ratio[i]),"\t(",l$conf.int.lo[i],"-",l$conf.int.hi[i],")\n")
		}
	}

	# return the variance captured by converged features
	var.sum
}


### map.starburst - compute and display the starburst representation of clusters
# parameters:
# - map is an object if type 'map'
# - explicit controls the shape of the connected components
# - smoothing controls the smoothing level of the umat (NULL,0,>0)

map.starburst <- function(map,explicit=FALSE,smoothing=2) {

	if (class(map) != "map")
		stop("map.starburst: first argument is not a map object.")

	umat <- compute.umat(map,smoothing=smoothing)
	plot.heat(map,umat,explicit=explicit,comp=TRUE)
}


### map.umat - compute and display the unified distance matrix
# parameters:
# - map is an object if type 'map'

map.umat <- function(map) {

	if (class(map) != "map")
		stop("map.umat: first argument is not a map object.")

	umat <- compute.umat(map,smoothing=NULL)
	plot.heat(map,umat,comp=FALSE)
}


### map.feature - compute and display the enhanced unified distance matrix for a 
#                 feature of the training data
# parameters:
# - map is an object if type 'map'
# - feature is an integer as the index of the feature. 
# - explicit controls the shape of the connected components
# - smoothing controls the smoothing level of the umat (NULL,0,>0)

map.feature <- function(map,feature,explicit=FALSE,smoothing=2) {

	if (class(map) != "map")
		stop("map.feature: first argument is not a map object.")

	if (feature > ncol(data.frame(map$data)) || feature < 1)
		stop("map.feature: illegal feature index.")

	heat <- compute.plane(map,feature,smoothing=smoothing)
	plot.heat(map,heat,explicit=explicit,comp=TRUE)
}


### map.projection - print the association of labels with map elements
# parameters:
# - map is an object if type 'map'
# return values:
# - a dataframe containing the projection onto the map for each observation

map.projection <- function(map) {

	if (class(map) != "map")
		stop("map.projection: first argument is not a map object.")

	labels <- map$labels
	if (is.null(labels))
		stop("map.projection: no labels available")
		
	x <- map$visual$x + 1
	y <- map$visual$y +1
	names(labels[[1]]) <- "labels"
	names(x[[1]]) <- "x"
	names(y[[1]]) <- "y"
	
	data.frame(labels,x,y)
}

### map.significance - compute the relative significance of each feature and plot it
# parameters:
# - map is an object if type 'map'
# - graphics is a switch that controls whether a plot is generated or not
# - feature.labels is a switch to allow the plotting of feature names vs feature indices
# return value:
# - a vector containing the significance for each feature

map.significance <- function(map,graphics=TRUE,feature.labels=TRUE) {

	if (class(map) != "map")
		stop("map.significance: first argument is not a map object.")

	data.df <- data.frame(map$data)
	nfeatures <- ncol(data.df)

	# Compute the variance of each feature on the map
	var.v <- array(data=1,dim=nfeatures)
	for (i in 1:nfeatures) { 
		var.v[i] <- var(data.df[[i]]);
	}

	# we use the variance of a feature as likelihood of
	# being an important feature, compute the Bayesian
	# probability of significance using uniform priors

	var.sum <- sum(var.v)
	prob.v <- var.v/var.sum

	# plot the significance
	if (graphics) {
		par.v <- map.graphics.set()	

		y <- max(prob.v)
		plot.new()
		plot.window(xlim=c(1,nfeatures),ylim=c(0,y))
		box()
		
		title(xlab="Features",ylab="Significance")
		
		xticks <- seq(1,nfeatures,1)
		yticks <- seq(0,y,y/4)
		if (feature.labels)
			xlabels <- names(data.df)
		else 
			xlabels <- seq(1,nfeatures,1)
		ylabels <- formatC(seq(0,y,y/4),digits=2)
		axis(1,at=xticks,labels=xlabels)
		axis(2,at=yticks,labels=ylabels)
		
		points(1:nfeatures,prob.v,type="h")

		map.graphics.reset(par.v)

	} else {
		prob.v
	}
}



############################### local functions #################################

# map.graphics.set -- set the graphics environment for our map utilities
#                     the return value is the original graphics param vector 

map.graphics.set <- function() {
	par.v <- par()
	par(ps=6)
	par.v
}

# map.graphics.reset -- reset the graphics environment to the original state
# parameter - a vector containing the settings for the original state

map.graphics.reset <- function(par.vector) {
	par(ps=par.vector$ps)
}

### plot.heat - plot a heat map based on a 'map', this plot also contains the connected
#               components of the map based on the landscape of the heat map
# parameters:
# - map is an object if type 'map'
# - heat is a 2D heat map of the map returned by 'map'
# - labels is a vector with labels of the original training data set
# - explicit controls the shape of the connected components
# - comp controls whether we plot the connected components on the heat map

plot.heat <- function(map,heat,explicit=FALSE,comp=TRUE) {

	labels <- map$labels
	if (is.null(labels))
		stop("plot.heat: no labels available")

	x <- map$xdim
	y <- map$ydim
	nobs <- nrow(map$data)
	count <- array(data=0,dim=c(x,y))

	# bin the heat values into 100 bins used for the 100 heat colors below
	heat.v <- as.vector(heat)
	heat.v <- cut(heat.v,breaks=100,labels=FALSE)
	heat <- array(data=heat.v,dim=c(x,y))
	
	# set up the graphics window
	par.v <- map.graphics.set()
	plot.new()
	plot.window(xlim=c(0,x),ylim=c(0,y))
	box()
	
	title(xlab="x",ylab="y")
	
	xticks <- seq(0.5,x-0.5,1)
	yticks <- seq(0.5,y-0.5,1)
	xlabels <- seq(1,x,1)
	ylabels <- seq(1,y,1)
	axis(1,at=xticks,labels=xlabels)
	axis(3,at=xticks,labels=xlabels)
	axis(2,at=yticks,labels=ylabels)
	axis(4,at=yticks,labels=ylabels)
		
	# plot heat
	colors<- heat.colors(100)
	
	for (ix in 1:x) {
		for (iy in 1:y) {
			rect(ix-1,iy-1,ix,iy,col=colors[100 - heat[ix,iy] + 1],border=NA)
		}
	}
	
	# put the connected component lines on the map
	if (comp) {

		# compute the connected components
		coords <- compute.internal.nodes(map,heat,explicit)

		for(ix in 1:x){
			for (iy in 1:y) {
				cx <- coords$xcoords[ix,iy]
				cy <- coords$ycoords[ix,iy]
				points(c(ix,cx)-.5,c(iy,cy)-.5,type="l",col="grey")
			}
		}
	}

	# put the labels on the map
	# count the labels in each map cell
	for(i in 1:nobs){
		ix <- map$visual$x[i]
		iy <- map$visual$y[i]
		count[ix+1,iy+1] <- count[ix+1,iy+1]+1
	}
	
	for(i in 1:nobs){
		ix <- map$visual$x[i]
		iy <- map$visual$y[i]
		# we only print one label per cell
		if (count[ix+1,iy+1] > 0) {
			count[ix+1,iy+1] <- 0
			ix <- ix + .5
			iy <- iy + .5
			l <- labels[i,1]
			text(ix,iy,labels=l)
		}
	}

	map.graphics.reset(par.v)
}


### compute.internal.nodes -- compute the centroid for each point on the map
# parameters:
# - map is an object if type 'map'
# - heat is a matrix representing the heat map representation
# - explicit controls the shape of the connected component
# return value:
# - a list containing the matrices with the same x-y dims as the original map containing the centroid x-y coordinates

compute.internal.nodes <- function(map,heat,explicit=FALSE) {
	xdim <- map$xdim
	ydim <- map$ydim
	x.coords <- array(data=-1,dim=c(xdim,ydim))
	y.coords <- array(data=-1,dim=c(xdim,ydim))
	max.val <- max(heat)
	
	find.internal.node <- function(ix,iy) {

		# first we check if the current position is already associated
		# with a centroid.  if so, simply return the coordinates
		# of that centroid
		if (x.coords[ix,iy] > -1 && y.coords[ix,iy] > -1) {
			list(bestx=x.coords[ix,iy],besty=y.coords[ix,iy])
		}

		# try to find a smaller value in the immediate neighborhood
		# make our current position the square with the minimum value.
		# if a minimum value other that our own current value cannot be
		# found then we are at a minimum.
		#
		# search the neighborhood; three different cases: inner element, corner element, side element
		min.val <- heat[ix,iy]
		min.x <- ix
		min.y <- iy
		
		# (ix,iy) is an inner map element
		if (ix > 1 && ix < xdim && iy > 1 && iy < ydim) {
			if (heat[ix-1,iy-1] < min.val) {
				min.val <- heat[ix-1,iy-1]
				min.x <- ix-1
				min.y <- iy-1
			}
			if (heat[ix,iy-1] < min.val) {
				min.val <- heat[ix,iy-1]
				min.x <- ix
				min.y <- iy-1
			}
			if (heat[ix+1,iy-1] < min.val) {
				min.val <- heat[ix+1,iy-1]
				min.x <- ix+1
				min.y <- iy-1
			}
			if (heat[ix+1,iy] < min.val) {
				min.val <- heat[ix+1,iy]
				min.x <- ix+1
				min.y <- iy
			}
			if (heat[ix+1,iy+1] < min.val) {
				min.val <- heat[ix+1,iy+1]
				min.x <- ix+1
				min.y <- iy+1
			}
			if (heat[ix,iy+1] < min.val) {
				min.val <- heat[ix,iy+1]
				min.x <- ix
				min.y <- iy+1
			}
			if (heat[ix-1,iy+1] < min.val) {
				min.val <- heat[ix-1,iy+1]
				min.x <- ix-1
				min.y <- iy+1
			}
			if (heat[ix-1,iy] < min.val) {
				min.val <- heat[ix-1,iy]
				min.x <- ix-1
				min.y <- iy
			}			
		} 
		
		# (ix,iy) is bottom left corner
		else if (ix == 1 && iy == 1) {
			if (heat[ix+1,iy] < min.val) {
				min.val <- heat[ix+1,iy]
				min.x <- ix+1
				min.y <- iy
			}
			if (heat[ix+1,iy+1] < min.val) {
				min.val <- heat[ix+1,iy+1]
				min.x <- ix+1
				min.y <- iy+1
			}
			if (heat[ix,iy+1] < min.val) {
				min.val <- heat[ix,iy+1]
				min.x <- ix
				min.y <- iy+1
			}
		}

		# (ix,iy) is bottom right corner
		else if (ix == xdim && iy == 1) {
			if (heat[ix,iy+1] < min.val) {
				min.val <- heat[ix,iy+1]
				min.x <- ix
				min.y <- iy+1
			}
			if (heat[ix-1,iy+1] < min.val) {
				min.val <- heat[ix-1,iy+1]
				min.x <- ix-1
				min.y <- iy+1
			}
			if (heat[ix-1,iy] < min.val) {
				min.val <- heat[ix-1,iy]
				min.x <- ix-1
				min.y <- iy
			}			
		} 

		# (ix,iy) is top right corner
		else if (ix == xdim && iy == ydim) {
			if (heat[ix-1,iy-1] < min.val) {
				min.val <- heat[ix-1,iy-1]
				min.x <- ix-1
				min.y <- iy-1
			}
			if (heat[ix,iy-1] < min.val) {
				min.val <- heat[ix,iy-1]
				min.x <- ix
				min.y <- iy-1
			}
			if (heat[ix-1,iy] < min.val) {
				min.val <- heat[ix-1,iy]
				min.x <- ix-1
				min.y <- iy
			}			
		}

		# (ix,iy) is top left corner
		else if (ix == 1 && iy == ydim) {
			if (heat[ix,iy-1] < min.val) {
				min.val <- heat[ix,iy-1]
				min.x <- ix
				min.y <- iy-1
			}
			if (heat[ix+1,iy-1] < min.val) {
				min.val <- heat[ix+1,iy-1]
				min.x <- ix+1
				min.y <- iy-1
			}
			if (heat[ix+1,iy] < min.val) {
				min.val <- heat[ix+1,iy]
				min.x <- ix+1
				min.y <- iy
			}			
		}
		
		# (ix,iy) is a left side element
		else if (ix == 1  && iy > 1 && iy < ydim) {
			if (heat[ix,iy-1] < min.val) {
				min.val <- heat[ix,iy-1]
				min.x <- ix
				min.y <- iy-1
			}
			if (heat[ix+1,iy-1] < min.val) {
				min.val <- heat[ix+1,iy-1]
				min.x <- ix+1
				min.y <- iy-1
			}
			if (heat[ix+1,iy] < min.val) {
				min.val <- heat[ix+1,iy]
				min.x <- ix+1
				min.y <- iy
			}
			if (heat[ix+1,iy+1] < min.val) {
				min.val <- heat[ix+1,iy+1]
				min.x <- ix+1
				min.y <- iy+1
			}
			if (heat[ix,iy+1] < min.val) {
				min.val <- heat[ix,iy+1]
				min.x <- ix
				min.y <- iy+1
			}
		}
		
		# (ix,iy) is a bottom side element
		else if (ix > 1 && ix < xdim && iy == 1 ) {
			if (heat[ix+1,iy] < min.val) {
				min.val <- heat[ix+1,iy]
				min.x <- ix+1
				min.y <- iy
			}
			if (heat[ix+1,iy+1] < min.val) {
				min.val <- heat[ix+1,iy+1]
				min.x <- ix+1
				min.y <- iy+1
			}
			if (heat[ix,iy+1] < min.val) {
				min.val <- heat[ix,iy+1]
				min.x <- ix
				min.y <- iy+1
			}
			if (heat[ix-1,iy+1] < min.val) {
				min.val <- heat[ix-1,iy+1]
				min.x <- ix-1
				min.y <- iy+1
			}
			if (heat[ix-1,iy] < min.val) {
				min.val <- heat[ix-1,iy]
				min.x <- ix-1
				min.y <- iy
			}			
		} 

		# (ix,iy) is a right side element
		else if (ix == xdim && iy > 1 && iy < ydim) {
			if (heat[ix-1,iy-1] < min.val) {
				min.val <- heat[ix-1,iy-1]
				min.x <- ix-1
				min.y <- iy-1
			}
			if (heat[ix,iy-1] < min.val) {
				min.val <- heat[ix,iy-1]
				min.x <- ix
				min.y <- iy-1
			}
			if (heat[ix,iy+1] < min.val) {
				min.val <- heat[ix,iy+1]
				min.x <- ix
				min.y <- iy+1
			}
			if (heat[ix-1,iy+1] < min.val) {
				min.val <- heat[ix-1,iy+1]
				min.x <- ix-1
				min.y <- iy+1
			}
			if (heat[ix-1,iy] < min.val) {
				min.val <- heat[ix-1,iy]
				min.x <- ix-1
				min.y <- iy
			}			
		} 
		
		# (ix,iy) is a top side element
		else if (ix > 1 && ix < xdim && iy == ydim) {
			if (heat[ix-1,iy-1] < min.val) {
				min.val <- heat[ix-1,iy-1]
				min.x <- ix-1
				min.y <- iy-1
			}
			if (heat[ix,iy-1] < min.val) {
				min.val <- heat[ix,iy-1]
				min.x <- ix
				min.y <- iy-1
			}
			if (heat[ix+1,iy-1] < min.val) {
				min.val <- heat[ix+1,iy-1]
				min.x <- ix+1
				min.y <- iy-1
			}
			if (heat[ix+1,iy] < min.val) {
				min.val <- heat[ix+1,iy]
				min.x <- ix+1
				min.y <- iy
			}
			if (heat[ix-1,iy] < min.val) {
				min.val <- heat[ix-1,iy]
				min.x <- ix-1
				min.y <- iy
			}			
		} 

		#if successful
		# move to the square with the smaller value, i.e., call find.internal.node on this new square
		# note the RETURNED x-y coords in the x.coords and y.coords matrix at the current location
		# return the RETURNED x-y coordinates
		if (min.x != ix || min.y != iy) {
			r.val <- find.internal.node(min.x,min.y)

			# if explicit is set show the exact connected component
			# otherwise construct a connected componenent where all
			# nodes are connected to a centrol node
			if (explicit) {
				x.coords[ix,iy] <<- min.x
				y.coords[ix,iy] <<- min.y
				list(bestx=min.x,besty=min.y)
			}
			else {
				x.coords[ix,iy] <<- r.val$bestx
				y.coords[ix,iy] <<- r.val$besty
				r.val
			}
		}
		#else
		# we have found a minimum
		# note the current x-y in the x.coords and y.coords matrix
		# return the current x-y coordinates
		else {
			x.coords[ix,iy] <<- ix
			y.coords[ix,iy] <<- iy
			list(bestx=ix,besty=iy)
		}
	} # end function find.internal.node

	# iterate over the map and find the centroid for each element
	for (i in 1:xdim) {
		for (j in 1:ydim) {
			find.internal.node(i,j)
		}
	}
	
	list(xcoords=x.coords,ycoords=y.coords)
}


### compute.umat -- compute the unified distance matrix
# parameters:
# - map is an object if type 'map'
# - smoothing is either NULL, 0, or a positive floating point value controlling the 
#         smoothing of the umat representation
# return value:
# - a matrix with the same x-y dims as the original map containing the umat values

compute.umat <- function(map,smoothing=NULL) {

	d <- dist(data.frame(map$code))
	umat <- compute.heat(map,d,smoothing)
	
	umat
}

### compute.plane -- compute the heat matrix of a plane
# parameters:
# - map is an object if type 'map'
# - plane is a plane index
# - smoothing is either NULL, 0, or a positive floating point value controlling the 
#         smoothing of the umat representation
# return value:
# - a matrix with the same x-y dims as the original map containing the umat values

compute.plane <- function(map,plane,smoothing=NULL) {

	map.df <- data.frame(map$code)

	if (plane < 1 && plane > ncol(map.df))
	   stop("compute.plane: bad plane index")

	d <- dist(map.df[[plane]])
	heat <- compute.heat(map,d,smoothing)
	
	heat
}

### compute.heat -- compute a heat value map representation of the given distance matrix
# parameters:
# - map is an object if type 'map'
# - d is a distance matrix computed via the 'dist' function
# - smoothing is either NULL, 0, or a positive floating point value controlling the 
#         smoothing of the umat representation
# return value:
# - a matrix with the same x-y dims as the original map containing the heat

compute.heat <- function(map,d,smoothing=NULL) {

	d <- as.matrix(d)
	x <- map$xdim
	y <- map$ydim
	heat <- array(data=0,dim=c(x,y))

	# this function translates our 2-dim map coordinates
	# into the 1-dim coordinates of the code vector
	xl <- function(ix,iy,xdim) {
	#cat("converting (",ix,",",iy,") to row", ix + (iy-1) *xdim,"\n")
		ix + (iy-1) * xdim
	}

	# iterate over the inner nodes and compute their umat values
	for (ix in 2:(x-1)) {
		for (iy in 2:(y-1)) {
			sum <- 
			       d[xl(ix,iy,x),xl(ix-1,iy-1,x)] +
				   d[xl(ix,iy,x),xl(ix,iy-1,x)] +
				   d[xl(ix,iy,x),xl(ix+1,iy-1,x)] +
				   d[xl(ix,iy,x),xl(ix+1,iy,x)] +
				   d[xl(ix,iy,x),xl(ix+1,iy+1,x)] +
				   d[xl(ix,iy,x),xl(ix,iy+1,x)] +
				   d[xl(ix,iy,x),xl(ix-1,iy+1,x)] +
				   d[xl(ix,iy,x),xl(ix-1,iy,x)]
			heat[ix,iy] <- sum/8
		}
		
	}
	
	# iterate over bottom x axis
	for (ix in 2:(x-1)) {
		iy <- 1
		sum <- 
			   d[xl(ix,iy,x),xl(ix+1,iy,x)] +
			   d[xl(ix,iy,x),xl(ix+1,iy+1,x)] +
			   d[xl(ix,iy,x),xl(ix,iy+1,x)] +
			   d[xl(ix,iy,x),xl(ix-1,iy+1,x)] +
			   d[xl(ix,iy,x),xl(ix-1,iy,x)]
		heat[ix,iy] <- sum/5
	}
	
	# iterate over top x axis
	for (ix in 2:(x-1)) {
		iy <- y
		sum <- 
		       d[xl(ix,iy,x),xl(ix-1,iy-1,x)] +
			   d[xl(ix,iy,x),xl(ix,iy-1,x)] +
			   d[xl(ix,iy,x),xl(ix+1,iy-1,x)] +
			   d[xl(ix,iy,x),xl(ix+1,iy,x)] +
			   d[xl(ix,iy,x),xl(ix-1,iy,x)]
		heat[ix,iy] <- sum/5
	}
	
	# iterate over the left y-axis
	for (iy in 2:(y-1)) {
		ix <- 1
		sum <- 
			   d[xl(ix,iy,x),xl(ix,iy-1,x)] +
			   d[xl(ix,iy,x),xl(ix+1,iy-1,x)] +
			   d[xl(ix,iy,x),xl(ix+1,iy,x)] +
			   d[xl(ix,iy,x),xl(ix+1,iy+1,x)] +
			   d[xl(ix,iy,x),xl(ix,iy+1,x)] 
		heat[ix,iy] <- sum/5
	}
	
	# iterate over the right y-axis
	for (iy in 2:(y-1)) {
		ix <- x
		sum <- 
		       d[xl(ix,iy,x),xl(ix-1,iy-1,x)] +
			   d[xl(ix,iy,x),xl(ix,iy-1,x)] +
			   d[xl(ix,iy,x),xl(ix,iy+1,x)] +
			   d[xl(ix,iy,x),xl(ix-1,iy+1,x)] +
			   d[xl(ix,iy,x),xl(ix-1,iy,x)]
		heat[ix,iy] <- sum/5
	}
	
	# bottem left corner
	ix <- 1
	iy <- 1
	sum <- 
		   d[xl(ix,iy,x),xl(ix+1,iy,x)] +
		   d[xl(ix,iy,x),xl(ix+1,iy+1,x)] +
		   d[xl(ix,iy,x),xl(ix,iy+1,x)] 
	heat[ix,iy] <- sum/3

	# bottem right corner
	ix <- x
	iy <- 1
	sum <- 
		   d[xl(ix,iy,x),xl(ix,iy+1,x)] +
		   d[xl(ix,iy,x),xl(ix-1,iy+1,x)] +
		   d[xl(ix,iy,x),xl(ix-1,iy,x)]
	heat[ix,iy] <- sum/3

	# top left corner
	ix <- 1
	iy <- y
	sum <- 
		   d[xl(ix,iy,x),xl(ix,iy-1,x)] +
		   d[xl(ix,iy,x),xl(ix+1,iy-1,x)] +
		   d[xl(ix,iy,x),xl(ix+1,iy,x)] 
	heat[ix,iy] <- sum/3

	# top right corner
	ix <- x
	iy <- y
	sum <- 
	       d[xl(ix,iy,x),xl(ix-1,iy-1,x)] +
		   d[xl(ix,iy,x),xl(ix,iy-1,x)] +
		   d[xl(ix,iy,x),xl(ix-1,iy,x)]
	heat[ix,iy] <- sum/3

	# smooth the heat map
	xcoords <- c()
	ycoords <- c()
	for (i in 1:y) { 
		for (j in 1:x) {
			ycoords <- c(ycoords, i)
			xcoords <- c(xcoords,j)
		}
	}
	xycoords <- data.frame(xcoords,ycoords)

	if (!is.null(smoothing)) {
		if (smoothing == 0) 
			heat <- smooth.2d(as.vector(heat),x=as.matrix(xycoords),nrow=x,ncol=y,surface=FALSE)
		else if (smoothing > 0)
			heat <- smooth.2d(as.vector(heat),x=as.matrix(xycoords),nrow=x,ncol=y,surface=FALSE,theta=smoothing)
		else
			stop("compute.heat: bad value for smoothing parameter")
	}

	heat
}

### df.var.test -- a function that applies the F-test testing the ratio
#                  of the variances of the two data frames
# parameters:
# - df1,df2 - data frames with the same number of columns
# - conf - confidence level for the F-test (default .95)

df.var.test <- function(df1,df2,conf = .95) {

	if (length(df1) != length(df2))
		stop("df.var.test: cannot compare variances of data frames")

	# init our working arrays
	var.ratio.v <- array(data=1,dim=length(df1))
	var.confintlo.v <- array(data=1,dim=length(df1))
	var.confinthi.v <- array(data=1,dim=length(df1))
	
	# compute the F-test on each feature in our populations
	for (i in 1:length(df1)) { 
		t <- var.test(df1[[i]],df2[[i]],conf.level=conf)
		var.ratio.v[i] <- t$estimate
		#cat("Feature",i,"confidence interval =",t$conf.int,"\n")
		var.confintlo.v[i] <- t$conf.int[1]
		var.confinthi.v[i] <- t$conf.int[2]
	}

	# return a list with the ratios and conf intervals for each feature
	list(ratio=var.ratio.v,conf.int.lo=var.confintlo.v,conf.int.hi=var.confinthi.v)
}



