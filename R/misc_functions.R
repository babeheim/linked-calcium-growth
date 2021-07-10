
insert_row <- function(d, after_row, label) {
  add_row <- matrix("", ncol = ncol(d), nrow = 1)
  add_row[,1] <- label
  rbind(
    d[1:after_row,],
    add_row,
    d[(after_row + 1):nrow(d),]
  )
}

extract_sampling_time <- function(s) {
  s <- gsub("Sampling .*, ", "", s)
  s <- gsub(" total", "", s)
  units <- gsub("\\d.*\\s", "", s)
  s <- as.numeric(gsub("\\s.*$", "", s))
  tar <- which(units == "seconds")
  if (length(tar) > 0) s[tar] <- s[tar]/60
  tar <- which(units == "hours")
  if (length(tar) > 0) s[tar] <- s[tar] * 60
  s <- round(s, 1)
  return(s)
}


difftime_to_char <- function(x) {
  paste(as.character(x), attr(x, "units"))
}


say_festival <- function(message) {
  system(paste0("echo '", message, "' | festival --tts"))
}


dir_init <- function(path, verbose=FALSE, overwrite=TRUE){
  if(substr(path, 1, 2)!='./') stop('path argument must be formatted
  with "./" at beginning')
  contents <- dir(path, recursive=TRUE)
  if(dir.exists(path)){
    if(overwrite){
      if(verbose){
          if(length(contents)==0) print(paste('folder ', path, ' created.', sep=""))
          if(length(contents)>0) print(paste('folder ', path, ' wiped of ', length(contents), ' files/folders.', sep=""))
        }
      if(dir.exists(path)) unlink(path, recursive=TRUE)
      dir.create(path)
    }
  } else {
  if(verbose){
    print(paste('folder ', path, ' created.', sep=""))
  }
    dir.create(path)
  }
}

texttab <- function(input.matrix, hlines=NA){
  output <- character(nrow(input.matrix))
  for(i in 1:nrow(input.matrix)){
    add.amps <- paste(input.matrix[i,], collapse=" & ")
    output[i] <- paste(add.amps, "\\\\", sep=" ")
  }
  if(all(!is.na(hlines))){
    for(i in 1:length(hlines)) output <- append(output, "\\midrule", hlines[i]+(i-1))
  }
  return(output)
}

col_alpha <- function (acol, alpha = 0.2){
  acol <- col2rgb(acol)
  acol.red <- acol["red",]/255
  acol.green <- acol["green",]/255
  acol.blue <- acol["blue",]/255
  acol <- mapply(
    function(red, green, blue, alphas) rgb(red, green, blue, alphas), acol.red, acol.green, acol.blue, alpha
  )
  return(as.character(acol))
}

gradient_maker <- function(start=NA, stop=NA, 
  cols=c("darkorange", "white", "darkcyan"), vis=FALSE, n=1000){
  if(is.na(start) | is.na(stop)) stop("need to specify 
    start and stop points on a numerical scale")
  colfunc <- colorRampPalette(cols)
  color.list <- colfunc(n)
  color.locations <- seq(start, stop, length=n)
  names(color.locations) <- color.list
  if(vis==TRUE) plot(color.locations, rep(1, n), 
    col=color.list, pch="|", ylim=c(0.9, 1.1), cex=5)
  return(color.locations)
}

data_gradient <- function(data, colors=c("darkorange", 
  "white", "darkcyan"), my.start=NA, my.stop=NA){
  if(is.na(my.start)) my.start <- min(data, na.rm=TRUE)
  if(is.na(my.stop)) my.stop <- max(data, na.rm=TRUE)
  my.gradient <- gradient_maker(start=my.start, 
    stop=my.stop, cols=colors)
  if(any(data > max(my.gradient), na.rm=T) | 
    any(data < min(my.gradient), na.rm=T)){
    warning("data is not within gradient range")
  }
  data.colors <- rep(NA, length(data))
  for(i in 1:length(data)){
    if(!is.na(data[i])){
      data.colors[i] <- names(my.gradient)[which.min(abs(data[i]-my.gradient))]
    }
  }
  data.colors
}
