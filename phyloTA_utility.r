#' \code{getSpp} parses a fasta file to return species names
#'
#' @param x a fasta read cluster from PhyloTA with define as taxon label
#' @details gives taxon labels, species names from a fasta object with annotations as given by phyloTA with 'Use define as taxon label'
#' @return a character vector
#' @examples \dontrun{
#' data(daphnia)
#' x <- daphnia@fasta
#' getSpp(x)
#' }

getSpp<-function(x){
          ans<-unlist(getAnnot(x))
          ans<-strsplit(ans,"[|]")
          hj<-function(x){x[3]}
          ans<-unlist(lapply(ans,hj))
          ans<-strsplit(ans," ")
          hj<-function(x){paste(x[1:2],collapse=" ")}
          ans<-unlist(lapply(ans,hj))
          ans
}

#' \code{getGen} parses a fasta file to return genus names
#'
#' @param x a fasta read cluster from PhyloTA with define as taxon label
#' @details gives taxon labels, species names from a fasta object with annotations as given by phyloTA with 'Use define as taxon label'
#' @return a character vector
#' @examples \dontrun{
#' data(daphnia)
#' x <- daphnia@fasta
#' getGen(x)
#' }

getGen<-function(x){
          ans<-unlist(getAnnot(x))
          ans<-strsplit(ans,"[|]")
          hj<-function(x){x[3]}
          ans<-unlist(lapply(ans,hj))
          ans<-strsplit(ans," ")
          hj<-function(x){x[1]}
          ans<-unlist(lapply(ans,hj))
          ans
}


#' \code{mti.match} returns multiple matches
#'
#' @param x a fasta read cluster from PhyloTA with define as taxon label
#' @details gives taxon labels, species names from a fasta object with annotations as given by phyloTA with 'Use define as taxon label'
#' @return a character vector
#' @examples \dontrun{
#' data(daphnia)
#' x <- daphnia@fasta
#' getSpp(x)
#' }

mti.match<-function(x,table){
  which(!is.na(match(table,x)))
}
