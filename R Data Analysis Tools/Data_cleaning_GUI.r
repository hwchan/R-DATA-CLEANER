# GUI for cleaning (somewhat) raw data for data analysis:
# removes duplicates, renames genes with their alias, removes genes with low FPKM, returns only selected genes 
# Written by Howie Ho Wai Chan (University of British Columbia) 
# April 17, 2014
# Please read the documentation for in-depth instructions on how to use and how to customize 

library(gWidgets)
library(gWidgetstcltk)
library(tools)

# Sets the main window
win <- gwindow("Data cleaning tool")

# Label and field for reference table
label1 <- glabel("Gene expression reference table", container = win)
gRef <- ggroup(container = win);
refFile <- gedit("Select a file...", width=75, container = gRef);
btn1 <- gbutton(text = "Browse", handler = function(h, ...){
	svalue(refFile) <- tclvalue(tkgetOpenFile(filetypes = "{{All files} *} {{Tab delimited} {.txt .dlm .tab}} "));
}, container = gRef);

# Label and field for alias list
label2 <- glabel("Gene alias list", container = win)
gAlias <- ggroup(container = win);
aliasFile <- gedit("Select a file...", width=75, container = gAlias);
btn2 <- gbutton(text = "Browse", handler = function(h, ...){
	svalue(aliasFile) <- tclvalue(tkgetOpenFile(filetypes = "{{All files} *} {{Tab delimited} {.txt .dlm .tab}} "));
}, container = gAlias);

# Label and field for specified gene list
label3 <- glabel("Specified gene list", container = win)
gSpec <- ggroup(container = win);
specFile <- gedit("Select a file...", width=75, container = gSpec);
btn3 <- gbutton(text = "Browse", handler = function(h, ...){
	svalue(specFile) <- tclvalue(tkgetOpenFile(filetypes = "{{All files} *} {{Tab delimited} {.txt .dlm .tab}} "));
}, container = gSpec);

# Label and field for output file
label4 <- glabel("Output file", container = win)
gOut <- ggroup(container = win);
outField <- gedit("Save as...", width=75, container = gOut);
btn4 <- gbutton(text = "Browse", handler = function(h, ...){
	svalue(outField) <- tclvalue(tkgetSaveFile(filetypes="{{Tab delimited text file} {.txt}}"));
	if(file_ext(svalue(outField)) != "txt"){
		svalue(outField) <- paste(svalue(outField), ".txt", sep="");
	}
}, container = gOut);

# Disable fields
enabled(aliasFile) <- F;
enabled(btn2) <- F;
enabled(specFile) <- F;
enabled(btn3) <- F;
			
# Checkboxes for options
gChk <- ggroup(container = win);
chkRmv <- gcheckbox(
	text = "Remove duplicates",
	checked = TRUE,
	container = gChk
)
chkRnm <- gcheckbox(
	text = "Rename with alias",
	checked = FALSE,
	handler	= function(h, ...){
		if(svalue(chkRnm)){
			enabled(aliasFile) <- T;
			enabled(btn2) <- T;
		}else{
			enabled(aliasFile) <- F;
			enabled(btn2) <- F;
		}
	},
	container = gChk
)

chkSpec <- gcheckbox(
	text = "Get specified genes",
	checked = FALSE,
	handler	= function(h, ...){
		if(svalue(chkSpec)){
			enabled(specFile) <- T;
			enabled(btn3) <- T;
		}else{
			enabled(specFile) <- F;
			enabled(btn3) <- F;
		}
	},
	container = gChk
)

gChkOut <- ggroup(container = win);
chkOut <- gcheckbox(
	text = "Remove genes that have < 3 cells with > 1 FPKM",
	checked = FALSE,
	container = gChkOut
)

# Check if input in field is valid
checkField <- function(input){
	if(svalue(input) != "" && svalue(input) != "Select a file..." && svalue(input) != "Select a folder..." && svalue(input) != "N/A" && svalue(input) != "Save as..." && svalue(input) != ".txt"){
		return(TRUE);
	}else{
		return(FALSE);
	}
}

doThings <- function(h, ...){
	# Check for missing fields
	if(svalue(chkRmv)){
		if(!checkField(refFile) || !checkField(outField)){
			svalue(statusBar) <- "Error: invalid files or fields.";
			return();
		}
	}
	if(svalue(chkRnm)){
		if(!checkField(refFile) || !checkField(outField) || !checkField(aliasFile)){
			svalue(statusBar) <- "Error: invalid files or fields.";
			return();
		}
	}
	if(svalue(chkOut)){
		if(!checkField(refFile) || !checkField(outField)){
			svalue(statusBar) <- "Error: invalid files or fields.";
			return();
		}
	}
	if(svalue(chkSpec)){
		if(!checkField(refFile) || !checkField(outField) || !checkField(specFile)){
			svalue(statusBar) <- "Error: invalid files or fields.";
			return();
		}
	}

	# Do everything in one pass
	tryCatch({
		myData <- read.table(svalue(refFile), row.names = NULL);
		
		# Remove duplicates
		if(svalue(chkRmv)){
			svalue(statusBar) <- "Removing duplicates with mean...";
			myData <- aggregate(.~row.names, data = myData, mean);
		}
		
		# Reset data format so the other methods can read (set the rownames as first column)
		rownames(myData) <- myData[,1];
		myData[,1] <- NULL;
		
		# Rename with alias
		if(svalue(chkRnm)){
			svalue(statusBar) <- "Renaming with aliases...";
			lookUp <- read.table(svalue(aliasFile), header=T);
			
			index <- match(lookUp[,"From"], row.names(myData));
			loo<-which(!is.na(index));
			row.names(myData)[na.omit(index)] <- as.character(lookUp[loo,][,"To"]);
		}
		
		# Remove poorly expressed genes (keep genes with > 2 cells, > 1 expression)
		if(svalue(chkOut)){
			svalue(statusBar) <- "Removing poorly expressed genes...";
			FPKM <- 1;			# FPKM value to be greater than to be satisfactory
			numberCells <- 2;	# number of cells (with enough FPKM) to be greater than to be satisfactory
			great3 <- function(x){
				length(x[x>FPKM])>numberCells	
			}
			index <- apply(myData,1,great3);	# 1 = apply by rows
			myData<-myData[index,];
		}
		
		# Return only specified genes
		if(svalue(chkSpec)){
			svalue(statusBar) <- "Returning specified genes...";
			specified <- read.table(svalue(specFile));

			index <- match(row.names(specified), row.names(myData));
			myData <- myData[na.omit(index),];
		}
		
		write.table(myData, svalue(outField), sep="\t");
		svalue(statusBar) <- paste("File saved at", svalue(outField));
		print(svalue(outField));
		flush.console();
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
startbtn <- gbutton(text = "Start", handler = doThings, container = gBtns);

