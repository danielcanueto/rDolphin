
## Internal function for parsing Bruker acquisition files
## inFile - string; directory containing the necessary acquisition files
## params - string; desired parameters to return, if missing will return
##           relevant paramaters for 1D or 2D file
## note: all values are returned as string arguments
## returns values for the designated acquisition parameters
topspin_read_spectrum2 <- function(partname, filename, minppm, maxppm){ 
  acquPar=parseAcqus(partname)
  pars=parseProcs(filename)
  storedpars <-append(pars, acquPar)
  if (is.null(storedpars$SI)) {
    storedpars$SI     = 0
    storedpars$OFFSET    = 0
    storedpars$SW     = 0
    storedpars$real   = NaN
    storedpars$NC_proc = 0
	storedpars$XDIM = 0
    return()
  }
  
  storedpars$RG = as.numeric(storedpars$RG)
  storedpars$OFFSET=as.numeric(storedpars$OFFSET)
  storedpars$SI=as.numeric(storedpars$SI)
  storedpars$NC_proc=as.numeric(storedpars$NC_proc)
  storedpars$SW=as.numeric(storedpars$SW)
  storedpars$XDIM=as.numeric(storedpars$XDIM)

 
  minppmindex = floor((storedpars$OFFSET-minppm)/storedpars$SW*(storedpars$SI-1))
  maxppmindex = ceiling((storedpars$OFFSET-maxppm)/storedpars$SW*(storedpars$SI-1))
  
  if (maxppmindex<0) {
    maxppmindex = 0
    minppmindex = storedpars$SI-1
  }
  realfile = paste(filename, '1r', sep='/')
  
  if (file.exists(realfile)==0) {
    sprintf('Error: file %s not found', realfile)
    storedpars$SI     = 0
    storedpars$OFFSET    = 0
    storedpars$SW     = 0
    storedpars$real   = NaN
    storedpars$NC_proc = 0
	    storedpars$XDIM = 0

  }
  
  
  readCon <- file(realfile, 'rb')
  data <- try(readBin(readCon, size=4, what='integer', n=storedpars$SI, endian=storedpars$BYTORDP),
              silent=TRUE)
  storedpars$real=data[(maxppmindex+1):(minppmindex+1)]
  close(readCon)
  
  storedpars$OFFSET = (storedpars$OFFSET - storedpars$SW/(storedpars$SI-1)*maxppmindex);
  storedpars$SW = (storedpars$OFFSET - storedpars$SW/(storedpars$SI-1)*maxppmindex) - (storedpars$OFFSET - storedpars$SW/(storedpars$SI-1)*minppmindex)
  storedpars$SI=length(storedpars$real)
  return(storedpars)
}





parseAcqus <- function(inDir, params){ 

  ## Designate parameters if not provided
  if (missing(params))
    params <- c('NC', 'RG', 'OVERFLW', 'NS', 'DATE')
  paramVar <- paste('##$', params, sep='')
  
  ## Search inDir for necessary acquisition parameter files
  acqus <- list.files(inDir, full.names=TRUE, pattern='^acqus$')[1]
  if (is.na(acqus)) acqus <- list.files(inDir, full.names=TRUE, pattern='^acqu$')[1]	
  if (is.na(acqus)) {
    paste('Could not find acquisition parameter files ("acqu" or "acqus")', 
               ' in:\n"', inDir, '".', sep='')
    return()
  }
  acqu2s <- list.files(inDir, full.names=TRUE, pattern='^acqu2s')[1]
  if (is.na(acqu2s))
    acqu2s <- list.files(inDir, full.names=TRUE, pattern='^acqu2$')[1]
  if (is.na(acqu2s)) files <- acqus else files <- c(acqus, acqu2s)
  
  ## Search acquisition files for designated parameters
  acquPar <- NULL
  for (i in seq_along(files)){
    
    ## Determine paramater/value separator
    for (paramSep in c('= ', '=', ' =', ' = ')){
      splitText <- strsplit(readLines(files[i]), paramSep)
      parNames <- sapply(splitText, function(x) x[1])
      parVals <- sapply(splitText, function(x) x[2])
      matches <- match(paramVar, parNames)
      if (any(is.na(matches)))
        next
      else
        break
    }
    
    ## Return an error if any parameters can not be found
    if (any(is.na(matches)))
      stop(paste('One or more of the following parameters could not be found: ', 
                 paste("'", params[which(is.na(matches))], "'", sep='', 
                       collapse=', '), ' in:\n"', files[i], sep=''))
    acquPar <- rbind(acquPar, parVals[matches])
  }
  
  ## Format the data
  colnames(acquPar) <- params
  acquPar <- data.frame(acquPar, stringsAsFactors=FALSE)
  if (!is.null(acquPar$NUC1)){
    for (i in seq_along(acquPar$NUC1)){
      acquPar$NUC1[i] <- unlist(strsplit(unlist(strsplit(acquPar$NUC1[i], 
                                                         '<'))[2], '>'))
    }
  }
  if (!is.na(acqu2s))
    rownames(acquPar) <- c('w2', 'w1')
  return(acquPar)
}


## Internal function for parsing Bruker processing files
## inFile - string; directory containing the necessary processing files
## params - string; desired parameters to return, if missing will return
##           relevant paramaters for 1D or 2D file
## note: all values are returned as string arguments
## returns values for the designated processing parameters
parseProcs <- function(inDir, params){ 
  
  ## Designate parameters if not provided
  if (missing(params))
    params <- c('SW_p','SF','SI','OFFSET','NC_proc','BYTORDP','XDIM')
  paramVar <- paste('##$', params, sep='')
  
  ## Search inDir for necessary processing parameter files
  procs <- list.files(inDir, full.names=TRUE, pattern='^procs$')[1]
  if (is.na(procs))
    procs <- list.files(inDir, full.names=TRUE, pattern='^proc$')[1]
  if (is.na(procs)) {
    paste('Could not find processing parameter files ("proc" or "procs")', 
               ' in:\n"', inDir, '".', sep='')
  return()
}
  proc2s <- list.files(inDir, full.names=TRUE, pattern='^proc2s$')[1]
  if (is.na(proc2s))
    proc2s <- list.files(inDir, full.names=TRUE, pattern='^proc2$')[1]
  if (is.na(proc2s))
    files <- procs
  else
    files <- c(procs, proc2s)
  
  ## Search processing files for designated parameters
  pars <- NULL
  for (i in seq_along(files)){
    
    ## Determine paramater/value separator
    for (paramSep in c('= ', '=', ' =', ' = ')){
      splitText <- strsplit(readLines(files[i]), paramSep)
      parNames <- sapply(splitText, function(x) x[1])
      parVals <- sapply(splitText, function(x) x[2])
      matches <- match(paramVar, parNames)
      if (any(is.na(matches)))
        next
      else
        break
    }
    
    ## Return an error if any parameters can not be found
    if (any(is.na(matches)))
      stop(paste('One or more of the following parameters could not be found: ', 
                 paste("'", params[which(is.na(matches))], "'", sep='', 
                       collapse=', '), ' in:\n"', files[i], sep=''))
    pars <- rbind(pars, parVals[matches])
  }
  
  ## Format the data
  colnames(pars) <- params
  pars <- data.frame(pars, stringsAsFactors=FALSE)
  if (!is.null(pars$BYTORDP))
    pars$BYTORDP <- ifelse(as.numeric(pars$BYTORDP), 'big', 'little')
  if (!is.na(proc2s))
    rownames(pars) <- c('w2', 'w1')
  pars$SW=as.numeric(pars$SW_p)/as.numeric(pars$SF)
  pars=pars[3:8]
  return(pars)
}


