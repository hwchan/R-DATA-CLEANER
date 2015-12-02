# GUI for creating Heatmaps, Violin plots, Correlation plots, Residuals to y=x, and Peak position calculations
# Written by Howie Ho Wai Chan (University of British Columbia) 
# April 17, 2014
# Please read the documentation for in-depth instructions on how to use and how to customize 

library(gWidgets)
library(gWidgetstcltk)
library(pheatmap);
library(RColorBrewer);
library(vioplot);

# Sets the main window
win <- gwindow("Plotting tool")

# Label and field for reference table
label1 <- glabel("Gene expression reference table", container = win)
gRef <- ggroup(container = win);
refFile <- gedit("Select a file...", width=75, container = gRef);
btn1 <- gbutton(text = "Browse", handler = function(h, ...){
	svalue(refFile) <- tclvalue(tkgetOpenFile(filetypes = "{{All files} *} {{Tab delimited} {.txt .dlm .tab}} "));
}, container = gRef);

# Label and field for specified gene list
label2 <- glabel("Specified gene list", container = win)
gSpec <- ggroup(container = win);
specFile <- gedit("Select a file...", width=75, container = gSpec);
btn2 <- gbutton(text = "Browse", handler = function(h, ...){
	svalue(specFile) <- tclvalue(tkgetOpenFile(filetypes = "{{All files} *} {{Tab delimited} {.txt .dlm .tab}} "));
}, container = gSpec);

# Label and field for bulk gene list
label3 <- glabel("Bulk gene list", container = win)
gBulk <- ggroup(container = win);
bulkFile <- gedit("Select a file...", width=75, container = gBulk);
btn3 <- gbutton(text = "Browse", handler = function(h, ...){
	svalue(bulkFile) <- tclvalue(tkgetOpenFile(filetypes = "{{All files} *} {{Tab delimited} {.txt .dlm .tab}} "));
}, container = gBulk);

enabled(bulkFile) <- F;
enabled(btn3) <- F;

# Label and field for output folder
label4 <- glabel("Output file", container = win)
gOut <- ggroup(container = win);
outField <- gedit("Save as...", width=75, container = gOut);
btn4 <- gbutton(text = "Browse", handler = function(h, ...){
	svalue(outField) <- tclvalue(tkgetSaveFile());
}, container = gOut);

# Radio box for log functions
gLog <- ggroup(container = win);
logRadio <- gradio(c("No log", "Log2(n)", "Log2(n+1)", "Log10(n)", "Log10(n+1)"), selected=3, horizontal=T, container=gLog);

# Checkbox for heatmap row labels
gCheck <- ggroup(container = win);
chkLabels <- gcheckbox(
  text      = "Heatmap row labels",
  checked   = F,
  container = gCheck
)

# Checkbox for using specified gene list
chkSpec <- gcheckbox(
	text	= "Plot only specified genes",
	checked	= F,
	handler = function(h, ...){
		if(svalue(chkSpec)){
			enabled(specFile) <- T;
			enabled(btn2) <- T;
		} else{
			enabled(specFile) <- F;
			enabled(btn2) <- F;
		}
	},container = gCheck
)
enabled(specFile) <- F;
enabled(btn2) <- F;
		
# Checkbox for printing residuals list
chkRes <- gcheckbox(
	text	= "Print residuals",
	checked = F,
	container = gCheck
)
enabled(chkRes) <- F;

# Radio box for plot functions
gPlot <- ggroup(container = win);
plotRadio <- gradio(c("Heatmap", "Violin Plot", "Correlation Plot", "Peak Calculations"), selected=1, horizontal=T, handler = function(h, ...){
	if(svalue(plotRadio) == "Heatmap"){
		enabled(refFile) <- T;
		enabled(btn1) <- T;
		
		enabled(bulkFile) <- F;
		enabled(btn3) <- F;
		enabled(outField) <- T;
		enabled(btn4) <- T;
		enabled(chkLabels) <- T;
		enabled(chkSpec) <- T;
		enabled(chkRes) <- F;
	}
	if(svalue(plotRadio) == "Violin Plot"){
		enabled(refFile) <- T;
		enabled(btn1) <- T;

		enabled(bulkFile) <- F;
		enabled(btn3) <- F;
		enabled(outField) <- T;
		enabled(btn4) <- T;
		enabled(chkLabels) <- F;
		enabled(chkSpec) <- T;
		enabled(chkRes) <- F;
	}
	if(svalue(plotRadio) == "Correlation Plot"){
		enabled(refFile) <- T;
		enabled(btn1) <- T;
		enabled(specFile) <- F;
		enabled(btn2) <- F;
		enabled(bulkFile) <- T;
		enabled(btn3) <- T;
		enabled(outField) <- T;
		enabled(btn4) <- T;
		enabled(chkLabels) <- F;
		enabled(chkSpec) <- F;
		svalue(chkSpec) <- F;
		enabled(chkRes) <- T;
	}
	if(svalue(plotRadio) == "Peak Calculations"){
		enabled(refFile) <- T;
		enabled(btn1) <- T;

		enabled(bulkFile) <- F;
		enabled(btn3) <- F;
		enabled(outField) <- T;
		enabled(btn4) <- T;
		enabled(chkLabels) <- F;
		enabled(chkSpec) <- T;
		enabled(chkRes) <- F;
	}
}, container=gPlot);



# Plot heatmap
plotHeatmap <- function(){
	colors <- colorRampPalette(c("blue", "white", "red"))(n = 30);
	myData <- logData(refFile);
	if(svalue(chkSpec)){
		lookUp <- logData(specFile);
		index <- match(row.names(lookUp), row.names(myData));
		myData <- myData[na.omit(index),];
	}	
	#Print - with row labels
	if(svalue(chkLabels)){
		out <- paste(svalue(outField), "(Heatmap with-labels).png");
		png(out, width = 5*300, height = 20*nrow(myData)+400, res = 300, pointsize = 8);
		pheatmap(myData, fontsize_row=4, fontsize_col=6, border_color=NA, breaks=c(-15:15), color=colors, treeheight_row=0);
		dev.off();
		svalue(statusBar) <- out;
		print(out);
		flush.console();
	}
	# Print - without row labels
	else{
		out <- paste(svalue(outField), "(Heatmap no-labels).png");
		png(out);
		pheatmap(myData, fontsize_col=9, border_color=NA, breaks=c(-15:15), color=colors, treeheight_row=0, show_rownames=F);
		dev.off();
		svalue(statusBar) <- out;
		print(out);
		flush.console();
	}
	
}

# Plot violin plot
plotViolin <- function(){
	myData <- logData(refFile);
	if(svalue(chkSpec)){
		lookUp <- logData(specFile);
		index <- match(row.names(lookUp), row.names(myData));
		myData <- myData[na.omit(index),];
	}	
	floorLim <- floor(min(myData));
	ceilLim <- ceiling(max(myData));
	i <- 1;
	while(i <= nrow(myData)){
		x <- myData[i,];
		out <- paste(svalue(outField), " (", row.names(x), ").png", sep="");
		png(out);
		vioplot(unlist(x), wex = .5, names = row.names(x), col = "white", lwd = 3, ylim = c(floorLim, ceilLim), drawRect = F);
		dev.off();
		print(out);
		flush.console();
		i <- i+1;
	}
	svalue(statusBar) <- paste("Finished printing violin plots at", svalue(outField));
}

# Plot correlation graph
plotCorr <- function(){
	x <- logData(bulkFile);	# iterates with j
	y <- logData(refFile);	# iterates with i

	# Calculate the graph size
	ceilLim <- ceiling(max(x,y));
	floorLim <- floor(min(x,y));
	
	# Print out the graphs & calculate residuals
	DF <- y;
	j <- 1;
	while(j <= ncol(x)){
		# Prints the graphs: (ncol(x) * ncol(y) /2) number of graphs so we don't have duplicates (e.g. C01-C04 & C04-C01) 
		i <- j;
		while(i <= ncol(y)){
			out <- paste(svalue(outField), " (", colnames(y)[i], "-", colnames(x)[j], ").png", sep="");
			# Print the graph with lines and labels
			png(out);
				plot(x[,j], y[,i], xlim=c(floorLim,ceilLim), ylim=c(floorLim,ceilLim), xlab=colnames(x)[j], ylab=colnames(y)[i], cex=.5, pch=20, main=paste(i,". ",colnames(y)[i]," to ",colnames(x)[j], " Correlation Plot\n(Log2[m+1])", sep="")); 
				abline(0, 1, col="red");					# y=x line
				#abline(0, 1, col="dimgray", lty="dotted");	# dotted gray y=x line
				#abline(lm(y[,i]~x[,1]), col="red");		# Regression line
				#abline(lm(y[,i]~x[,1]-1), col="blue");		# Fake regression line (y=0 intercept)
				text(floorLim+3,ceilLim-.5,paste("Pearson:", round(cor(x[,j],y[,i]),2)),cex=1.5);
			dev.off();
			if(svalue(chkRes)){
				# Residuals calculation (distance between line and point: (x-y)/sqrt(2))
				DF[,i] <- (x[,j]-y[,i]);
			}
			print(out);
			flush.console();
			i <- i+1;
		}
		
		# Calculate the second half of the residuals. This cannot be in the above loop because that loop skips over half
		i <- 1;
		if(svalue(chkRes)){
			while(i < j){
				# Residuals calculation (distance between line and point: (x-y)/sqrt(2))
				DF[,i] <- (x[,j]-y[,i]);
				i <- i+1;
			}
		}
		# Print residuals
		DF <- DF/sqrt(2);
		out2 <- paste(svalue(outField), " Residuals", " (", colnames(y)[i], ").txt", sep="");
		write.table(DF, out2, sep = "\t");
		print(out2);
		flush.console();
		j <- j+1;
	}
	if(svalue(chkRes)){
		svalue(statusBar) <- paste("Finished printing correlation plots and residuals list at", svalue(outField));
	}else{
		svalue(statusBar) <- paste("Finished printing correlation plots at", svalue(outField));
	}
}

# Calculate the x-coordinate of peak locations
calcPeak <- function(){
	myData <- logData(refFile);
	if(svalue(chkSpec)){
		lookUp <- logData(specFile);
		index <- match(row.names(lookUp), row.names(myData));
		myData <- myData[na.omit(index),];
	}	

	# Create a soon-to-be data frame with 5 columns
	DF <- c(0,0,0,0,0);

	# Test if the peak is legitimate. A legitimate peak is determined by a local maximum that has a y-value that's at least 25% of the highest y-value
	# a = x location of peak
	# d = density of gene
	test <- function(a,d){
		threshold <- max(d$y)/4
		if(d$y[a] > threshold){
			return(d$x[a]);
		}
		return(0);
	}
	# Find possible peaks (local maximum)
	i <- 1;
	while(i <= nrow(myData)){
		gene <- as.numeric(myData[i,]);
		d <- density(gene, adjust=1.5);
		peakLoc <- which(diff(sign(diff(d$y)))==-2);
		# Test if the peak is legitimate
		peaks <- sapply(peakLoc, test, d);
		DF <- rbind(DF, peaks=peaks[seq(DF)]);
		rownames(DF)[i+1] <- row.names(myData[i,]);
		i <- i+1;
	}
	DF[is.na(DF)] <- 0;
	DF <- subset(DF, rowSums(DF > 0) > 1);
	# Remove empty columns
	DF <- DF[,colSums(DF>0) > 0];
	out <- paste(svalue(outField), " (Peak positions).txt", sep="");
	write.table(DF, out, sep = "\t");
	svalue(statusBar) <- paste("Finished printing", out);
	print(out);
	flush.console();
}

# Take the logarithm of the input
logData <- function(input){
	myData <- read.table(svalue(input));
	
	if(svalue(logRadio) == "No log"){
		rankDataLog <- myData;
	}
	if(svalue(logRadio) == "Log2(n)"){
		rankDataLog <- log2(myData);
		rankDataLog[rankDataLog == -Inf] = 0;
	}
	if(svalue(logRadio) == "Log2(n+1)"){
		rankDataLog <- log2(myData+1);
	}
	if(svalue(logRadio) == "Log10(n)"){
		rankDataLog <- log10(myData);
		rankDataLog[rankDataLog == -Inf] = 0;
	}
	if(svalue(logRadio) == "Log10(n+1)"){
		rankDataLog <- log10(myData+1);
	}
	
	return(rankDataLog);
}

# Check if input in field is valid 
checkField <- function(input){
	if(svalue(input) != "" && svalue(input) != "Select a file..." && svalue(input) != "Select a folder..." && svalue(input) != "N/A" && svalue(input) != "Save as..."){
		return(TRUE);
	}else{
		return(FALSE);
	}
}

# Start button function
plotStart <- function(h, ...){
	tryCatch({
		if(svalue(plotRadio) == "Heatmap"){
			if(checkField(refFile) && checkField(outField)){
				svalue(statusBar) <- "Printing heatmap...";
				plotHeatmap();
			}else{
				svalue(statusBar) <- "Error: invalid files or fields.";
			}
		}
		if(svalue(plotRadio) == "Violin Plot"){
			if(checkField(refFile) && checkField(outField)){
				svalue(statusBar) <- "Printing violin plots...";
				plotViolin();
			}else{
				svalue(statusBar) <- "Error: invalid files or fields.";
			}
		}
		if(svalue(plotRadio) == "Correlation Plot"){
			if(checkField(refFile) && checkField(outField) && checkField(bulkFile)){
				svalue(statusBar) <- "Printing correlation plots...";
				plotCorr();
			}else{
				svalue(statusBar) <- "Error: invalid files or fields.";
			}
		}
		if(svalue(plotRadio) == "Peak Calculations"){
			if(checkField(refFile) && checkField(outField)){
				svalue(statusBar) <- "Calculating peak locations...";
				calcPeak();
			}else{
				svalue(statusBar) <- "Error: invalid files or fields.";
			}
		}
	},
	error = function(e){
		svalue(statusBar) <- "Error reading files. Check console for details.";
		print(e);
	})
}

# Status bar
statusBar <- gstatusbar("", container = win)

# Start button
gBtns <- ggroup(container = win);
startbtn <- gbutton(text = "Start", handler = plotStart, container =gBtns);
