
se <- function(z) sd(z, na.rm = TRUE) / sqrt(sum(!is.na(z)))

calc_twotail_p <- function(z) 2 * pmin(pnorm(z), 1 - pnorm(z))

psign <- function(samples, n_digits = 3) {
  if (mean(samples) > 0) output <- round(mean(samples < 0), n_digits)
  if (mean(samples) < 0) output <- round(mean(samples > 0), n_digits)
  float_digits <- paste0("%.", n_digits, "f")
  output <- sprintf(float_digits, output)
  null_entry <- paste0("0.", paste0(rep("0", n_digits), collapse = ""))
  output[output  ==  null_entry] <- paste0("$<0.", paste0(rep("0", n_digits - 1), collapse = ""),"1$")
  return(output)
}

prep_latex_variables <- function(named_list) {
  out <- character()
  for (i in 1:length(named_list)) {
    out[i] <- paste0("\\newcommand{\\", names(named_list)[i], "}{", named_list[[i]], "}")
  }
  return(out)
}

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
  tar <- which(units  ==  "seconds")
  if (length(tar) > 0) s[tar] <- s[tar] / 60
  tar <- which(units  ==  "hours")
  if (length(tar) > 0) s[tar] <- s[tar] * 60
  s <- round(s, 1)
  return(s)
}

difftime_to_char <- function(x) {
  paste(as.character(x), attr(x, "units"))
}

dir_init <- function(path, verbose=FALSE, overwrite=TRUE) {
  if (substr(path, 1, 2) != './') stop('path argument must be formatted
  with "./" at beginning')
  contents <- dir(path, recursive=TRUE)
  if (dir.exists(path)) {
    if (overwrite) {
      if (verbose) {
          if (length(contents) == 0) print(paste('folder ', path, ' created.', sep = ""))
          if (length(contents)>0) print(paste('folder ', path, ' wiped of ', length(contents), ' files / folders.', sep = ""))
        }
      if (dir.exists(path)) unlink(path, recursive=TRUE)
      dir.create(path)
    }
  } else {
  if (verbose) {
    print(paste('folder ', path, ' created.', sep = ""))
  }
    dir.create(path)
  }
}

texttab <- function(input.matrix, hlines=NA) {
  output <- character(nrow(input.matrix))
  for (i in 1:nrow(input.matrix)) {
    add.amps <- paste(input.matrix[i,], collapse = " & ")
    output[i] <- paste(add.amps, "\\\\", sep = " ")
  }
  if (all(!is.na(hlines))) {
    for (i in seq_along(hlines)) output <- append(output, "\\midrule", hlines[i] + (i - 1))
  }
  return(output)
}

col_alpha <- function (acol, alpha = 0.2) {
  acol <- col2rgb(acol)
  acol.red <- acol["red",] / 255
  acol.green <- acol["green",] / 255
  acol.blue <- acol["blue",] / 255
  acol <- mapply(
    function(red, green, blue, alphas) rgb(red, green, blue, alphas),
      acol.red, acol.green, acol.blue, alpha
  )
  return(as.character(acol))
}
