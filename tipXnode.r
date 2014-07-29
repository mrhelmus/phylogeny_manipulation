#Grafts a tip on a tree halfway the distance between two supplied nodes, added branch length is the mean bl of the decendent tips


require(ape)
require(geiger)
tipXnode<-function(tree,addtip,where.nodes=NULL){
  emat<-tree$edge
  elen<-tree$edge.length[apply(emat,1,paste,collapse="")==paste(where.nodes,collapse="")]
  nedg<-elen/2
  tps<-tips(tree,where.nodes[2])
  ind<-sapply(tps,match,tree$tip)
  kl<-dist.nodes(tree)[where.nodes[2],ind]
  nedgadd<-mean(kl)+nedg
  txt<-paste("(",addtip,":",nedgadd,");",sep="")
  add<-read.tree(text=txt)
  nieuw<-bind.tree(tree,add,where=where.nodes[2],position=nedg)
  return(nieuw)
}

