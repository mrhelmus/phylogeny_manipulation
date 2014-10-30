#Grafts a tip on a tree halfway the distance between a supplied tip and its node

tipXtip<-function(tree,addtip,where.tip=NULL){
  emat<-tree$edge
  efoc<-which.edge(tree, where.tip)
  foctip<-match(where.tip,tree$tip)
  elen<-tree$edge.length
  nedg<-elen[efoc]/2
  txt<-paste("(",addtip,":",nedg,");",sep="")
  add<-read.tree(text=txt)
  nieuw<-bind.tree(tree,add,where=foctip,position=nedg)
  return(nieuw)
}



