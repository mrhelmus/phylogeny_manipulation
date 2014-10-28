#' \code{comparative.comm} creates a community comparative ecology object
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
          require(seqinr)
          ans<-unlist(getAnnot(x))
          ans<-strsplit(ans,"[|]")
          hj<-function(x){x[3]}
          ans<-unlist(lapply(ans,hj))
          ans<-strsplit(ans," ")
          hj<-function(x){paste(x[1:2],collapse=" ")}
          ans<-unlist(lapply(ans,hj))
}

