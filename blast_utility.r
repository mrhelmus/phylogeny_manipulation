#' BLAST to NCBI
#' \code{blastn.organism} Send a BLAST query to NCBI for a gene ID or string
#'
#' @param x sequence as a character vector or an integer corresponding to an entrez gene ID
#' @param organism name of NCBI clades in which to search for BLAST hits
#' @param database NCBI database to use, default “nr”
#' @param hitListSize number of hits to keep
#' @param filter sequence filter; “L” for Low Complexity, “R” for Human Repeats, “m” for Mask lookup
#' @param expect the BLAST ‘expect’ value above which matches will be returned
#' @param attempts the number of times to query the NCBI server until an stop error is trigered
#' @param snooze length of time, in seconds, to wait before another query attempt is attempted
#'
#' @details BLASTs a supplied sequence to the genbank nucleotide database within the archieved sequences of a supplied NCBI taxonomic name.
#' The NCBI API is documented at \link{http://www.ncbi.nlm.nih.gov/blast/Doc/urlapi.html}.
#' @return a data.frame with ids, definitions and e-values for all significant hits
#' @author M.R. Helmus edited Kevin Keenan's version of R. Gentleman \code{blastSequences} code in \code{annotate} \link{http://rstudio-pubs-static.s3.amazonaws.com/12097_1352791b169f423f910d93222a4c2d85.html}
#' @examples \dontrun{
#'  # Search for a gi number
#'  #gi|300791433|Daphnia mendotae isolate G2/2 cytochrome oxidase subunit I (COI) gene, partial cds; mitochondrial.
#'   blastn.organism(300791433, organism="Daphnia", database = "nr", hitListSize = "10", filter = "L", expect = "1e-10", attempts = 10, snooze=5)
#'
#'  # Search for a sequence
#'   myseq<-"TATTTTTGGAATTTGGTCTGGGATAGTCGGAACCGCTCTTAGTTTACTGATCCGGGCTGAACTTGGACAATCAGGAAGATTAATTGGGGATGACCAAATTTACAATGTAATTGTAACTGCCCACGCTTTTGTAATACTTTTT"
#'   blastn.organism(myseq, organism="Daphnia", database = "nr", hitListSize = "10", filter = "L", expect = "1e-10", attempts = 10, snooze=5)
#' }
blastn.organism <- function (x, organism="Animalia", database = "nr", hitListSize = "10",
                        filter = "L", expect = "1e-10", attempts = 10,snooze=5)
{
  program <- "blastn"
  baseUrl <- "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
  query <- paste("QUERY=", as.character(x), "&DATABASE=", database,
                 "&HITLIST_SIZE=", hitListSize, "&FILTER=", filter, "&EXPECT=",
                 expect, "&PROGRAM=", program,
                 "&ENTREZ_QUERY=", organism, sep = "")
  url0 <- sprintf("%s?%s&CMD=Put", baseUrl, query)
  results <- tempfile()
  Sys.sleep(snooze)
  require(XML)
  post <- htmlTreeParse(url0, useInternalNodes = TRUE)
  X <- post[["string(//comment()[contains(., \"QBlastInfoBegin\")])"]]
  rid <- sub(".*RID = ([[:alnum:]]+).*", "\\1", X)
  rtoe <- as.integer(sub(".*RTOE = ([[:digit:]]+).*", "\\1",X))
  url1 <- sprintf("%s?RID=%s&FORMAT_TYPE=XML&CMD=Get", baseUrl, rid)
  Sys.sleep(rtoe)
  .tryParseResult <- function(url, attempts){
    for (i in 1:(attempts+1)) {
      result <- tryCatch({
        xmlTreeParse(url, useInternalNodes=TRUE,
                     error = xmlErrorCumulator(immediate=FALSE))
      }, error=function(err) NULL)
      if (!is.null(result)) return(result)
      Sys.sleep(snooze)
    }
    alarm()
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
    return(data.frame(t(res)))
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

