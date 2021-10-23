
library(ape)
library(ggtree)
library(treeio)


###  Attempt1 full tree without node support - TEST TREE
# First need to edit the tip labels
tree<- read.tree("tree_totara_peronospora_mrbayes_FINALTREEFILE.newick")
#tiplabel<-tree$tip.label
#write.csv(tiplabel, "tiplabel.csv") # add column with changed tips and column with strain name
# Now plot
info <- read.csv("tiplabel.csv")
p <- ggtree(tree) %<+% info+ xlim(NA, 0.3)
 p+  geom_tiplab(aes(label=paste0('italic(', species, ')~', strain)), parse=T, size=1.3)+
    geom_hilight(node=158, fill="grey")+
    geom_hilight(node=161, fill="grey")+
    #geom_cladelabel(node=158, label="15", fontsize=1.5, offset=0.06)+ # add cladelabels in powerpoint
    #geom_cladelabel(node=161, label="5", fontsize=1.5, offset=0.07)+
    geom_treescale(x=0.2, y=-3, fontsize=1.5)


#  ----FINAL TREES WITH SUPPORT VALUES -----
# Read in tree
tree<-read.beast("tree_totara_peronospora_mrbayes_figtree.nexus") # for node support export nexus tree with metacomments from geneious and midpointrooted in figtree
tree@data
#Read in tiplabels
info <- read.csv("tiplabel.csv")

# final full tree
#p1 is full tree midpoint-rooted with posterior probability support values 
p1 <- ggtree(tree) %<+% info + xlim(NA, 0.3)
p1+  geom_tiplab(aes(label=paste0('italic(', species, ')~', strain)), parse=T, size=1.8)+geom_rootedge(rootedge = 0.005)+
    geom_text2(aes(label=round(as.numeric(`Posterior Probability`),2), x=branch), size = 1.8, vjust=0)+
    geom_hilight(node=147, fill="grey")+
    geom_hilight(node=150, fill="grey")+
    #geom_cladelabel(node=147, label="15", fontsize=2, offset = 0.01)+ # add cladelabels in powerpoint
    #geom_cladelabel(node=150, label="5", fontsize=2, offset = 0.01)+
    geom_treescale(x=0.2, y=-3, fontsize=1.8)
# Save 
ggsave("Finaltreeforsubmission.tiff", units="in", width = 10, height = 9, device='tiff', dpi=700) # add clade label in powerpoint
ggsave("Finaltreeforsubmission.pdf", units="in", width = 10, height = 9, device='pdf', dpi=700)
dev.off()

# to look up node label etc
pdataframe<-as.data.frame(p1 %>% as.treedata %>% as_tibble) 

# final subtree with clades 5,15,3
pdataframe<-as.data.frame(p2 %>% as.treedata %>% as_tibble) 
newtree<-drop.tip(tree, c(1:34, 48:99), subtree=FALSE, trim.internal = TRUE) #remove nodes as numbered in pdataframe

p2 <- ggtree(newtree) %<+% info+ xlim(NA, 0.45)
p2+  geom_tiplab(aes(label=paste0('italic(', species, ')~', strain)), parse=T, size=2.8)+geom_rootedge(rootedge = 0.005)+
    #geom_text2(aes(label=round(as.numeric(`Posterior Probability`),2), x=branch, hjust=1.5), size = 2.3)+ # somehow not printing at the correct nodes so add in powerpoint
    geom_cladelabel(node=23, label="15", fontsize=2.8, offset= 0.16, offset.text = 0.01)+ 
    geom_cladelabel(node=26, label="5", fontsize=2.8, offset = 0.16, offset.text = 0.01)+
    geom_cladelabel(node=18, label="3", fontsize=2.8, offset = 0.16, offset.text = 0.01)+
    geom_treescale(x=0.2, y=-1, fontsize=2.5)
# Save 
ggsave("Finaltreeforsubmissioncollapsed.tiff", units="in", width = 6, height = 5, device='tiff', dpi=700) # add clade label in powerpoint
ggsave("Finaltreeforsubmissioncollapsed.pdf", units="in", width = 6, height = 5, device='pdf', dpi=700)

write.table(newtree@data$`Posterior Probability`, "posterior.txt")

dev.off()
ggsave("test.tiff", units="in", width = 6, height = 5, device='tiff', dpi=700) # add clade label in powerpoint
ggsave("Finaltreeforsubmissioncollapsed.pdf", units="in", width = 6, height = 5, device='pdf', dpi=700)

# Not used but could be useful next time
# midpointrooting tree 
#library(phangorn)
#phylotree<-as.phylo(tree) # phylotree loses branch support values
#tree_midroot <- midpoint(phylotree, node.labels = "support")
    
