
# function definition
blastKK <- function (x, database = "nr", hitListSize = "10",
                        filter = "L", expect = "10", program = "blastn",
                        attempts = 10) {
  baseUrl <- "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
  query <- paste("QUERY=", as.character(x), "&DATABASE=", database,
                 "&HITLIST_SIZE=", hitListSize, "&FILTER=", filter, "&EXPECT=",
                 expect, "&PROGRAM=", program, sep = "")
  url0 <- sprintf("%s?%s&CMD=Put", baseUrl, query)
  results <- tempfile()
  Sys.sleep(5)
  require(XML)
  post <- htmlTreeParse(url0, useInternalNodes = TRUE)
  x <- post[["string(//comment()[contains(., \"QBlastInfoBegin\")])"]]
  rid <- sub(".*RID = ([[:alnum:]]+).*", "\\1", x)
  rtoe <- as.integer(sub(".*RTOE = ([[:digit:]]+).*", "\\1",
                         x))
  url1 <- sprintf("%s?RID=%s&FORMAT_TYPE=XML&CMD=Get", baseUrl,
                  rid)
  Sys.sleep(rtoe)
  .tryParseResult <- function(url, attempts){
    for (i in 1:(attempts+1)) {
      result <- tryCatch({
        xmlTreeParse(url, useInternalNodes=TRUE,
                     error = xmlErrorCumulator(immediate=FALSE))
      }, error=function(err) NULL)
      if (!is.null(result)) return(result)
      Sys.sleep(10)
    }
    stop(paste("no results after ", attempts,
               " attempts; please try again later", sep = ""))
  }
  result <- .tryParseResult(url1, attempts)
  qseq <- xpathApply(result, "//Hsp_qseq", xmlValue)
  hseq <- xpathApply(result, "//Hsp_hseq", xmlValue)
  require(Biostrings)
  res <- list()
  for (i in seq_len(length(qseq))) {
    res[i] <- DNAMultipleAlignment(c(hseq[[i]], qseq[[i]]),
                                   rowmask = as(IRanges(), "NormalIRanges"), colmask = as(IRanges(),
                                                                                          "NormalIRanges"))
  }
  res
}

# function definition
blast <- function (x, database = "nr", hitListSize = "10",
                        filter = "L", expect = "10", program = "blastn",
                        attempts = 10) {
  baseUrl <- "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
  query <- paste("QUERY=", as.character(x), "&DATABASE=", database,
                 "&HITLIST_SIZE=", hitListSize, "&FILTER=", filter, "&EXPECT=",
                 expect, "&PROGRAM=", program, sep = "")
  url0 <- sprintf("%s?%s&CMD=Put", baseUrl, query)
  results <- tempfile()
  Sys.sleep(5)
  require(XML)
  post <- htmlTreeParse(url0, useInternalNodes = TRUE)
  x <- post[["string(//comment()[contains(., \"QBlastInfoBegin\")])"]]
  rid <- sub(".*RID = ([[:alnum:]]+).*", "\\1", x)
  rtoe <- as.integer(sub(".*RTOE = ([[:digit:]]+).*", "\\1",
                         x))
  url1 <- sprintf("%s?RID=%s&FORMAT_TYPE=XML&CMD=Get", baseUrl,
                  rid)
  Sys.sleep(rtoe)
  .tryParseResult <- function(url, attempts){
    for (i in 1:(attempts+1)) {
      result <- tryCatch({
        xmlTreeParse(url, useInternalNodes=TRUE,
                     error = xmlErrorCumulator(immediate=FALSE))
      }, error=function(err) NULL)
      if (!is.null(result)) return(result)
      Sys.sleep(10)
    }
    stop(paste("no results after ", attempts,
               " attempts; please try again later", sep = ""))
  }
  result <- .tryParseResult(url1, attempts)
  #qseq <- xpathApply(result, "//Hsp_qseq", xmlValue)
  #hseq <- xpathApply(result, "//Hsp_hseq", xmlValue)
  hspev <- xpathApply(result, "//Hsp_evalue", xmlValue)
  hitid <- xpathApply(result, "//Hit_id", xmlValue)
  hitdef <- xpathApply(result, "//Hit_def", xmlValue)
  hspnum <- xpathApply(result, "//Hsp_num", xmlValue)
  res<-NULL
  kk<-0
  for(k in 1:length(hspev))
  {
    if(hspnum[k]==1){
      kk<-kk+1
      res<-rbind(res,c(hitid[kk],hitdef[kk],hspnum[k],hspev[k]))
    } else {res<-rbind(res,c(hitid[kk],hitdef[kk],hspnum[k],hspev[k]))}
  }
  colnames(res)<-c("id","definition","hitnum","evalue")
  return(data.frame(res))
}

# function definition
blast.organism <- function (x, organism="Animalia", database = "nr", hitListSize = "10",
                        filter = "L", expect = "10", program = "blastn",
                        attempts = 10) {
  baseUrl <- "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
  query <- paste("QUERY=", as.character(x), "&DATABASE=", database,
                 "&HITLIST_SIZE=", hitListSize, "&FILTER=", filter, "&EXPECT=",
                 expect, "&PROGRAM=", program,
                 "&ENTREZ_QUERY=", organism, sep = "")
  url0 <- sprintf("%s?%s&CMD=Put", baseUrl, query)
  results <- tempfile()
  #Sys.sleep(5)
  require(XML)
  post <- htmlTreeParse(url0, useInternalNodes = TRUE)
  x <- post[["string(//comment()[contains(., \"QBlastInfoBegin\")])"]]
  rid <- sub(".*RID = ([[:alnum:]]+).*", "\\1", x)
  rtoe <- as.integer(sub(".*RTOE = ([[:digit:]]+).*", "\\1",x))
  url1 <- sprintf("%s?RID=%s&FORMAT_TYPE=XML&CMD=Get", baseUrl, rid)
  #Sys.sleep(rtoe)
  .tryParseResult <- function(url, attempts){
    for (i in 1:(attempts+1)) {
      result <- tryCatch({
        xmlTreeParse(url, useInternalNodes=TRUE,
                     error = xmlErrorCumulator(immediate=FALSE))
      }, error=function(err) NULL)
      if (!is.null(result)) return(result)
      Sys.sleep(10)
    }
    stop(paste("no results after ", attempts,
               " attempts; please try again later", sep = ""))
  }
  result <- .tryParseResult(url1, attempts)
  #qseq <- xpathApply(result, "//Hsp_qseq", xmlValue)
  #hseq <- xpathApply(result, "//Hsp_hseq", xmlValue)
  hspev <- xpathApply(result, "//Hsp_evalue", xmlValue)
  if(length(hspev)==0){
    res<-c(NA,NA,NA,NA)
    names(res)<-c("id","definition","hitnum","evalue")
    return(data.frame(res))
  } else {
    hitid <- xpathApply(result, "//Hit_id", xmlValue)
    hitdef <- xpathApply(result, "//Hit_def", xmlValue)
    hspnum <- xpathApply(result, "//Hsp_num", xmlValue)
    res<-NULL
    kk<-0
    for(k in 1:length(hspev))
    {
      if(hspnum[k]==1){
        kk<-kk+1
        res<-rbind(res,c(hitid[kk],hitdef[kk],hspnum[k],hspev[k]))
      } else {res<-rbind(res,c(hitid[kk],hitdef[kk],hspnum[k],hspev[k]))}
    }
    colnames(res)<-c("id","definition","hitnum","evalue")
    return(data.frame(res))
  }
}
