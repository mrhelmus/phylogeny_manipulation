#' Calculate phylogenetic independent correlations
#'
#'
#' @param trait1 spp x trait
#' @param trait2 spp x trait
#' @param traits a spp x traits object
#' @param tree a phylo object that matches the spp in traits
#' @return a symetrical correlation matrix based on phylogenetic independent contrasts
#' @author M.R. Helmus

picor<-function(trait1,trait2,tree)
{
  pic(trait1,tree)%*%pic(trait2,tree)/sqrt(sum(pic(trait1,tree)^2)*sum(pic(trait2,tree)^2))# must be on a dicotomous tree
}

pic.traits<-function(traits,tree)
{
  require(miscTools)
  cors<-NULL
  ntraits<-dim(traits)[2]
  nspp<-dim(traits)[1]
  for(i in 1:(ntraits))
  {
    trait1<-traits[,i]
    for(j in (i):ntraits)
    {
      trait2<-traits[,j]
      cors<-c(cors,picor(trait1,trait2,tree))
    }
  }
  #make matrix
  mat<-symMatrix(cors)
  rownames(mat)<-colnames(traits)
  colnames(mat)<-colnames(traits)
  #return(cors)
  return(mat)
}

#p<-pic.traits(traits,tree.resolved)
#p
#round(p-cor(traits),3)
