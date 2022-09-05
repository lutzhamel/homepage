# variance based feature selection

# (c) 2009-2012 Lutz Hamel, University of Rhode Island

# This is the R code used for the Bayesian feature selection described in the paper:
#
# "Bayesian Probability Approach to Feature Significance for Infrared 
# Spectra of Bacteria", Lutz Hamel, Chris W. Brown, Applied Spectroscopy, Volume 66, Number 1, 2012.
# (preprints of these papers are available at www.cs.uri.edu/~hamel)

### License
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.

###
# select features by probability, 
# returns a data frame with the narrowed data set

select.p <- function(data,p=0.9,top=TRUE) {
	
	data.df <- as.data.frame(data)
	
	# compute the probabilities, sort them, and then
	# add them until we hit our threshold p

	prob.v <- prob.vector(data.df,decreasing=top)
	prob <- 0
	p.ix <- 1

	while (prob < p) {
		prob <- prob + prob.v$sorted[p.ix]
		p.ix <- p.ix + 1
	}
	
	# generate all the features that are part of this
	# confidence interval

	wn <- prob.v$index[1]
	data.df.narrow <- data.frame(data.df[[wn]])
	names(data.df.narrow) <- c(names(data.df)[wn])

	for (i in 2:(p.ix-1)) { 		
		wn <- prob.v$index[i]
		x <- data.frame(data.df[[wn]])
		names(x) <- c(names(data.df)[wn])
		data.df.narrow <- data.frame(data.df.narrow,x)
	}

	data.df.narrow
}


###
# select n features, 
# returns a data frame with the narrowed data set

select.n <- function(data,n=5,top=TRUE) {
	
	data.df <- as.data.frame(data)

	# compute the probabilities, sort them, and then
	# add them until we hit our threshold p

	prob.v <- prob.vector(data.df,decreasing=top)

	# generate all the features that are part of this
	# interval

	wn <- prob.v$index[1]
	data.df.narrow <- data.frame(data.df[[wn]])
	names(data.df.narrow) <- c(names(data.df)[wn])

	for (i in 2:n) { 		
		wn <- prob.v$index[i]
		x <- data.frame(data.df[[wn]])
		names(x) <- c(names(data.df)[wn])
		data.df.narrow <- data.frame(data.df.narrow,x)
	}

	data.df.narrow
}


###
# select n features at random, 
# returns a data frame with the narrowed data set

select.nrandom <- function(data,n=5) {
	
	data.df <- as.data.frame(data)

	# select n features from the set of features 
	# randomly without replacemente
	ix.v <- sample(1:ncol(data.df),n,replace=FALSE)
	
	# select the features in the data set
	data.df.narrow <- data.df[,ix.v]
	
	data.df.narrow
}


###
# compute the probability with which the n
# features represent significance

prob.n <- function(data,n,top=TRUE) {
	
	data.df <- as.data.frame(data)

	# compute the probabilities, sort them, and then
	# add them until we hit our threshold n

	prob.v <- prob.vector(data.df,decreasing=top)

	prob <- 0
	for (i in 1:n) {
		prob <- prob + prob.v$sorted[i]
	}
	
	prob
}

###
# compute a probability vector for the features in data
# returns a list of various probability values
# features - list of the probabilities of the features in the order as they 
#            appear in the original data set
# sorted -   list of the probabilities given in the order specified
# index -    gives the index to the original value

prob.vector <- function(data,decreasing=TRUE) {

	data.df <- as.data.frame(data)
	nfeatures <- ncol(data.df)
	
	# Compute the variance of each feature in the data
	
	var.v <- array(data=1,dim=nfeatures)
	
	for (i in 1:nfeatures) { 
		var.v[i] <- var(data.df[[i]]);
	}

	# we use the variance of a feature as likelihood of
	# being an important feature, turn the likelihood
	# into a probability by scaling.  this is
	# bayes rules with constant priors

	var.sum <- sum(var.v)
	prob.v <- var.v/var.sum

	# now figure out how many probabilities of top
	# features we have to add together in order to reach 
	# a confidence of p.
	# note: returns a list with two components: $x - sorted values:
	#		$ix - index to the original values

	prob.v.sorted <- sort(prob.v,decreasing=decreasing,index.return=TRUE)

	# we return a list of various probability values
	# features - list of the probabilities of the features in the order as they 
	#            appear in the original data set
	# sorted -   list of the probabilities given in the order specified
	# index -    gives the index to the original value
    list(features=prob.v,sorted=prob.v.sorted$x,index=prob.v.sorted$ix)
}

###
# plot the significance of the sorted features.

features <- function(data,graphics=TRUE,decreasing=TRUE,xlabels=TRUE,p=0.9,type="l") {

	data.df <- as.data.frame(data)
	nfeatures <- ncol(data.df)

	prob.v <- prob.vector(data.df,decreasing=decreasing)

	
	# find the features that finished the probability curve
	prob <- 0
	p.ix <- 0

	while (prob < p) {
		p.ix <- p.ix + 1
		prob <- prob + prob.v$sorted[p.ix]
	}

	# plot the significance
	if (graphics) {
		y <- max(prob.v$sorted)
		plot.new()
		plot.window(xlim=c(1,nfeatures),ylim=c(0,y))
		box()
		
		title(xlab="Features",ylab="Significance")
		
		xticks <- seq(1,nfeatures,1)
		yticks <- seq(0,y,y/4)
		if (xlabels)
			xlabels <- names(data.df)[prob.v$index]
		else 
			xlabels <- seq(1,nfeatures,1)
		ylabels <- formatC(seq(0,y,y/4),digits=2)
		axis(1,at=xticks,labels=xlabels)
		axis(2,at=yticks,labels=ylabels)
		
		points(1:nfeatures,prob.v$sorted,type=type)

		points(c(p.ix,p.ix),c(0,y),type="l",col="red")
		
		if (p.ix <= nfeatures/2)
			align <- 4
		else
			align <- 2
			
		text(p.ix,y/4*3,paste("p =",p),pos=align)
	} else {
		prob.v
	}
}


###
# plot the significance of the features.

significance <- function(data,graphics=TRUE,type="l") {

	data.df <- as.data.frame(data)
	nfeatures <- ncol(data.df)

	prob.v <- prob.vector(data.df)$features

	# plot the significance
	if (graphics) {
		y <- max(prob.v)
		plot.new()
		plot.window(xlim=c(1,nfeatures),ylim=c(0,y))
		box()
		
		title(xlab="Features",ylab="Significance")
		
		xticks <- seq(1,nfeatures,1)
		yticks <- seq(0,y,y/4)
		xlabels <- names(data.df)
		ylabels <- formatC(seq(0,y,y/4),digits=2)
		axis(1,at=xticks,labels=xlabels)
#		axis(3,at=xticks,labels=xlabels)
		axis(2,at=yticks,labels=ylabels)
#		axis(4,at=yticks,labels=ylabels)

		points(1:nfeatures,prob.v,type=type)
	} else {
		prob.v
	}
}

