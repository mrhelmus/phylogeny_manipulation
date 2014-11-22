#Grafts a tip on a tree halfway the distance between two supplied nodes, added branch length is the mean bl of the decendent tips

# simple wrap function for the bind.tip function of phytools
# to graft a tip as a polytomy to a node

require(ape)
require(phytools)
tipXpoly<-function(tree,addtip,where=NULL,edge.length=NULL,position=0)
{
  nieuw<-bind.tip(tree=tree,tip.label=addtip,edge.length=edge.length,where=where,position=position)
  return(nieuw)
}

