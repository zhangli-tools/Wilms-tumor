library(readxl)
dev.signatures=read_xlsx("marker.gene.sets.xlsx",sheet = 4)
dev.signatures=as.data.frame(dev.signatures)
dev.signatures.list=lapply(unique(dev.signatures$"Differentiation phase"),function(x){
  dev.signatures$Gene[dev.signatures$"Differentiation phase"==x]
})
names(dev.signatures.list)=unique(dev.signatures$"Differentiation phase")
library(GSVA)
WT.RNA.expr=readRDS("../RNA/WT.RNA.expr.RDS")
dev.RNA.deg.score=gsva(log2(WT.RNA.expr+.1),gset.idx.list = dev.signatures.list,method="ssgsea")
WT.protein.expr=readRDS("../Protein/WT.protein.expr.RDS")
clinical.all=readRDS("../Clinical/clinical.all.RDS")
dev.pt.deg.score=gsva(WT.protein.expr[,clinical.all$`Pre-treated_(Chemotherapy)`[match(colnames(WT.protein.expr),clinical.all$SampleID)]==0],gset.idx.list = dev.signatures.list,method="ssgsea")
subtype=readRDS("../Figure-4/subtype.RDS")
RNA.sub=paste0("S",subtype[colnames(dev.RNA.deg.score)])
RNA.sub[clinical.all$Tissue[match(colnames(dev.RNA.deg.score),clinical.all$SampleID)]=="N"]="NAT"
RNA.sub[RNA.sub=="SNA"]=NA
pt.sub=paste0("S",subtype[colnames(dev.pt.deg.score)])
pt.sub[clinical.all$Tissue[match(colnames(dev.pt.deg.score),clinical.all$SampleID)]=="N"]="NAT"
pt.sub[pt.sub=="SNA"]=NA
cal.stat=function(score,group){
  score1=scales::rescale(score,to=c(-3,3))
  avg=tapply(score1,group,function(x){mean(x,na.rm=T)})
  cmps=combn(unique(group[!is.na(group) & group!="SNA"]),2)
  pval=apply(cmps,2,function(t){
    t.test(score~factor(group,levels = t))$p.value
  })
  names(pval)=apply(cmps,2,function(x){paste(x,collapse = ".vs.")})
  y=c(avg,-log10(pval))
  y
}
score.mat=lapply(rownames(dev.RNA.deg.score),function(x){
  RNA=cal.stat(dev.RNA.deg.score[x,],RNA.sub)
  pt=cal.stat(dev.pt.deg.score[x,],pt.sub)
  y=rbind(RNA=RNA[1:4],Protein=pt[1:4])
  rownames(y)=paste(x,rownames(y),sep="_")
  y
})
pval.mat=lapply(rownames(dev.RNA.deg.score),function(x){
  RNA=cal.stat(dev.RNA.deg.score[x,],RNA.sub)
  pt=cal.stat(dev.pt.deg.score[x,],pt.sub)
  y=rbind(RNA=RNA[5:10],Protein=pt[5:10])
  rownames(y)=paste(x,rownames(y),sep="_")
  y
})
sig.mat=lapply(pval.mat,function(mat){
  res=apply(mat,2,function(x){
    p=10^(-x)
    y=rep("",length(p))
    y[p<0.1]="."
    y[p<0.05]="*"
    y[p<0.01]="**"
    y[p<0.001]="***"
    y
  })
})
x.mat=lapply(pval.mat,function(mat){
  res=t(apply(mat,1,function(x){
    y=seq(0,1,length.out=length(x)+1)
    dif=y[2]-y[1]
    y=y[1:(length(y)-1)]+dif/2
  }))
  res=res[2:1,]
})
y.mat=lapply(pval.mat,function(mat){
  res=apply(mat,2,function(x){
    y=seq(-0.1,1,length.out=length(x)+1)
    dif=y[2]-y[1]
    y=y[1:(length(y)-1)]+dif/2
  })
  res=res[2:1,]
})

library(ComplexHeatmap)
p1=Heatmap(eval(as.call(c(rbind,score.mat))),cluster_columns  = F,cluster_rows = F,name = "Score",row_names_side = "left",border = T,row_split = factor(rep(rownames(dev.RNA.deg.score),each=2),levels = rownames(dev.RNA.deg.score)),col = RColorBrewer::brewer.pal(n = 11,name = "RdBu")[11:1],row_title=F,rect_gp = gpar(col = "white", lwd = 1.5))
p2=Heatmap(eval(as.call(c(rbind,pval.mat))),cluster_columns  = F,cluster_rows = F,name = "Pvalue",row_names_side = "left",show_row_names   = F,col  =pal_material(palette = "orange")(10),border = T,row_split = factor(rep(rownames(dev.pt.deg.score),each=2),levels = rownames(dev.pt.deg.score)),row_title=F,rect_gp = gpar(col = "white", lwd = 1.5))
ht.list=p1 + p2
draw(ht.list)
for(i in 1:length(sig.mat)){
  decorate_heatmap_body("Pvalue", {
    grid.text(c(sig.mat[[i]]), c(x.mat[[i]]), c(y.mat[[i]]), default.units = "npc")
  }, slice = i)
}
## Fig 5A end

## Fig 5B
MET.markers=c("CDH4", "CDH6", "PAX3", "FGFR2", "FGF1", "FGF7")
MET.markers.RNA=sapply(MET.markers,function(x){
  sub=paste0("S",subtype[colnames(WT.RNA.expr)])
  RNA.tissue.type=clinical.all$Tissue[match(colnames(WT.RNA.expr),clinical.all$SampleID)]
  sub[RNA.tissue.type=="N"]="NAT"
  sub[sub=="SNA"]=NA
  y=scales::rescale(log2(WT.RNA.expr[x,]+1),to=c(-3,3))
  tapply(y,sub,function(t){mean(t)})
})
MET.markers.RNA.sig=sapply(MET.markers,function(x){
  sub=paste0("S",subtype[colnames(WT.RNA.expr)])
  RNA.tissue.type=clinical.all$Tissue[match(colnames(WT.RNA.expr),clinical.all$SampleID)]
  sub[RNA.tissue.type=="N"]="NAT"
  sub[sub=="SNA"]=NA
  y=sapply(c("S1","S2","S3"),function(t){
    m=wilcox.test(log2(WT.RNA.expr[x,which(sub==t)]+1),log2(WT.RNA.expr[x,which(sub=="NAT")]+1))$p.value
    if(m<0.001){m="***"}else if(m<0.01){m="**"}else if(m<0.05){m="*"}else{m=""}
  })
  y=c("",y)
})
MET.markers.pt=sapply(MET.markers,function(x){
  sub=paste0("S",subtype[colnames(WT.protein.expr)])
  sub[clinical.all[match(colnames(WT.protein.expr),clinical.all$SampleID),]$Tissue=="N"]="NAT"
  sub[sub=="SNA"]=NA
  if(x %in% rownames(WT.protein.expr))
  {y=scales::rescale(WT.protein.expr[x,],to=c(-3,3))
  y=tapply(y,sub,function(t){mean(t)})}else{y=c(NA,NA,NA,NA)}
  y
})
MET.markers.pt.sig=sapply(MET.markers,function(x){
  sub=paste0("S",subtype[colnames(WT.protein.expr)])
  sub[clinical.all[match(colnames(WT.protein.expr),clinical.all$SampleID),]$Tissue=="N"]="NAT"
  sub[sub=="SNA"]=NA
  if(x %in% rownames(WT.protein.expr))
  {y=sapply(c("S1","S2","S3"),function(t){
    m=wilcox.test(log2(WT.protein.expr[x,which(sub==t)]+1),log2(WT.protein.expr[x,which(sub=="NAT")]+1))$p.value
    if(m<0.001){m="***"}else if(m<0.01){m="**"}else if(m<0.05){m="*"}else{m=""}
  })
  }else{
    y=c("","","")
  }
  y=c("",y)
})
ComplexHeatmap::pheatmap(cbind(t(MET.markers.RNA),t(MET.markers.pt)),cluster_rows = F,cluster_cols = F,color=RColorBrewer::brewer.pal(11,name = "RdBu")[11:1],display_numbers = cbind(t(MET.markers.RNA.sig),t(MET.markers.pt.sig)),column_split=factor(rep(c("RNA","Protein"),each=4),levels=c("RNA","Protein")),cellwidth = 20,cellheight = 25,border_gp = gpar(col = "black"),name="Z-score",border_color = NA,column_title="Mesenchymal-epithelial transition (MET) signatures")
### Figure 5B end








####### Fig 5C
limma.deg=function(expr,label)
{
  library(limma)
  design.subtype=cbind(yes=label,no=1-label)
  subtype.fit = lmFit(expr, design.subtype)
  subtype.fit = eBayes(subtype.fit)
  subtype.contrast.matrix = makeContrasts(CancervNormal =yes - 
                                            no, levels = design.subtype)
  subtype.fit2 = contrasts.fit(subtype.fit, subtype.contrast.matrix)
  subtype.fit2 = eBayes(subtype.fit2)
  subtype.results = topTable(subtype.fit2, number = nrow(expr), sort.by = "p", 
                             adjust.method = "BH")
  
}

library(data.table)
TFs=fread("Science.Table S5.txt")
TFs.by.stage=TFs[TFs$qval<0.05,]
TFs.df=data.frame(TF=TFs.by.stage$Symbol,group=colnames(TFs.by.stage[,6:8])[apply(as.matrix(TFs.by.stage[,6:8]),1,which.max)])
library(clusterProfiler)
subtype.vs.nat.deg.RNA=lapply(1:3,function(i){
  RNA.tissue.type=clinical.all$Tissue[match(colnames(WT.RNA.expr),clinical.all$SampleID)]
  nat=WT.RNA.expr[,RNA.tissue.type=="N"]
  wt=WT.RNA.expr[,names(which(subtype==i))]
  cls=rep(c(0,1),c(ncol(nat),ncol(wt)))
  y=limma.deg(log2(cbind(nat,wt)+.1),cls)
})
subtype.vs.nat.deg.pt=lapply(1:3,function(i){
  protein.tissue.type=clinical.all$Tissue[match(colnames(WT.protein.expr),clinical.all$SampleID)]
  nat=WT.protein.expr[,protein.tissue.type=="N"]
  wt=WT.protein.expr[,names(which(subtype==i))]
  cls=rep(c(0,1),c(ncol(nat),ncol(wt)))
  y=limma.deg(cbind(nat,wt),cls)
})
obj.tf.RNA.subtype.vs.nat=lapply(subtype.vs.nat.deg.RNA,function(x){
  stat=x$t
  names(stat)=rownames(x)
  obj.tf.RNA=GSEA(sort(stat,decreasing = T),TERM2GENE = as.matrix(TFs.df[,2:1]),pAdjustMethod = "none",pvalueCutoff = 1)
})
obj.tf.pt.subtype.vs.nat=lapply(subtype.vs.nat.deg.pt,function(x){
  stat=x$t
  names(stat)=x$ID
  obj.tf.RNA=GSEA(sort(stat,decreasing = T),TERM2GENE = as.matrix(TFs.df[,2:1]),pAdjustMethod = "none",pvalueCutoff = 1)
})

subtype.pairwise.deg.RNA=lapply(1:3,function(i){
  #nat=WT.RNA.expr[,RNA.tissue.type=="N"]
  ij=combn(1:3,2)[,i]
  wt=WT.RNA.expr[,names(subtype)]
  ij.subtype=subtype[subtype %in% ij]
  y=limma.deg(log2(wt[,subtype %in% ij]+.1),label = ifelse(ij.subtype==c(ij[1]),1,0))
})
subtype.pairwise.deg.pt=lapply(1:3,function(i){
  #nat=WT.RNA.expr[,RNA.tissue.type=="N"]
  ij=combn(1:3,2)[,i]
  wt=WT.protein.expr[,names(subtype)]
  ij.subtype=subtype[subtype %in% ij]
  y=limma.deg(wt[,subtype %in% ij],label = ifelse(ij.subtype==c(ij[1]),1,0))
})
obj.tf.RNA.subtype.pairwise=lapply(subtype.pairwise.deg.RNA,function(x){
  stat=x$t
  names(stat)=rownames(x)
  obj.tf.RNA=GSEA(sort(stat,decreasing = T),TERM2GENE = as.matrix(TFs.df[,2:1]),pAdjustMethod = "none",pvalueCutoff = 1)
})
obj.tf.pt.subtype.pairwise=lapply(subtype.pairwise.deg.pt,function(x){
  stat=x$t
  names(stat)=x$ID
  obj.tf.RNA=GSEA(sort(stat,decreasing = T),TERM2GENE = as.matrix(TFs.df[,2:1]),pAdjustMethod = "none",pvalueCutoff = 1)
})
RNA.NES.pval=lapply(c(obj.tf.RNA.subtype.pairwise,obj.tf.RNA.subtype.vs.nat),function(x){
  x=x@result[,c("Description","NES","pvalue")]
})
pt.NES.pval=lapply(c(obj.tf.pt.subtype.pairwise,obj.tf.pt.subtype.vs.nat),function(x){
  x=x@result[,c("Description","NES","pvalue")]
})
names(RNA.NES.pval)=c("S1 vs S2","S1 vs S3","S2 vs S3","S1 vs NAT","S2 vs NAT","S3 vs NAT")
names(pt.NES.pval)=c("S1 vs S2","S1 vs S3","S2 vs S3","S1 vs NAT","S2 vs NAT","S3 vs NAT")
df.NES.pval=rbind(cbind(Comp=rep(names(RNA.NES.pval),each=3),
                        eval(as.call(c(rbind,RNA.NES.pval))),data="RNA"),
                  cbind(Comp=rep(names(pt.NES.pval),each=3),
                        eval(as.call(c(rbind,pt.NES.pval))),data="Protein"))
df.NES.pval$pvalue=ifelse(df.NES.pval$pvalue<0.05,1,0)
df.NES.pval$Comp=factor(df.NES.pval$Comp,levels=names(pt.NES.pval)[c(4:6,1:3)])
df.NES.pval$data=factor(df.NES.pval$data,levels = c("RNA","Protein"))
ggscatter(data = df.NES.pval,x = "Comp",y="Description",color = "NES",size = "pvalue",facet.by = "data")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_color_gsea()+
  theme(panel.grid.major = element_line(colour = "lightgrey",linetype = 3,size = 0.5),
        strip.text = element_text(size = 12))+xlab("Comparisons")+ylab("Developmental phases")
#### Fig 5C end

### Fig 5D

#######
RNA.genes=lapply(c(subtype.vs.nat.deg.RNA,subtype.pairwise.deg.RNA),function(x){
  down=rownames(x)[x$logFC< log2(2/4) & x$adj.P.Val<0.05]
  up=rownames(x)[x$logFC> log2(2) & x$adj.P.Val<0.05]
  y1=intersect(down,TFs.df$TF[TFs.df$group=="UB"])
  y2=intersect(up,TFs.df$TF[TFs.df$group=="PV"])
  y3=intersect(up,TFs.df$TF[TFs.df$group=="CM"])
  y=list(y1,y2,y3)
})
pt.genes=lapply(c(subtype.vs.nat.deg.pt,subtype.pairwise.deg.pt),function(x){
  down=x$ID[x$logFC<0 & x$adj.P.Val<0.05]
  up=x$ID[x$logFC>0 & x$adj.P.Val<0.05]
  y1=intersect(down,TFs.df$TF[TFs.df$group=="UB"])
  y2=intersect(up,TFs.df$TF[TFs.df$group=="PV"])
  y3=intersect(up,TFs.df$TF[TFs.df$group=="CM"])
  y=list(y1,y2,y3)
})

RNA.genes.res=unique(unlist(c(RNA.genes[[1]],
                       RNA.genes[[2]][c(1,3)],
                       RNA.genes[[3]],
                       RNA.genes[[4]],
                       RNA.genes[[5]][c(1,3)]
)))
pt.genes.res=unique(unlist(c(unlist(pt.genes[[1]]),
                      unlist(pt.genes[[2]][c(1,3)]),
                      unlist(pt.genes[[3]]),
                      unlist(pt.genes[[4]]),
                      unlist(pt.genes[[5]][c(1,3)])
)))

overlap.genes.tf=intersect(RNA.genes.res,pt.genes.res)
cor.TF.vs.target=cor(t(WT.RNA.expr[,intersect(colnames(WT.RNA.expr),colnames(WT.protein.expr))]),t(WT.protein.expr[overlap.genes.tf,intersect(colnames(WT.RNA.expr),colnames(WT.protein.expr))]),method="spearman")
library(clusterProfiler)
load("TF.database.Rdata")
db.tf.all=readRDS("db.tf.all.RDS")
obj.list=lapply(c(intersect(colnames(cor.TF.vs.target),c(names(CHEA),names(JASPAR),names(which(sapply(overlap.genes.tf,function(x){length(grep(x,unique(db.tf.all[,1])))})>0))))),function(x){
  stat=cor.TF.vs.target[,x]
  obj=GSEA(sort(stat,decreasing = T),TERM2GENE = db.tf.all[grep(x,db.tf.all[,1]),],maxGSSize = 2000,pvalueCutoff = 1,pAdjustMethod = "none")
  res=obj@result
  if(length(grep(x,res$Description))>0)
  {
    res=cbind(TF=x,res[grep(x,res$Description),])
    res=res[res$NES>0 & res$pvalue<0.05,]
  }else{res}
  res
})
obj.list=obj.list[sapply(obj.list,nrow)>0]
#obj.list[[2]]=obj.list[[2]][1,]
#obj.list[[8]]=obj.list[[8]][1,]
#obj.list[[16]]=obj.list[[16]][3,]
#obj.list[[18]]=obj.list[[18]][1,]
obj.list.mat=eval(as.call(c(rbind,obj.list)))
obj.list.mat$TF=factor(obj.list.mat$TF,levels = unique(obj.list.mat$TF[order(obj.list.mat$NES)]))
## remove some duplicated TFs manually
ggscatter(obj.list.mat,y = "NES",x = "TF",color = "pvalue",size = 5)+scale_color_material("orange",reverse = T)+
  theme(axis.line=element_blank(),axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "right",panel.border = element_rect(fill = NA,colour = "black"),panel.grid.major.x  = element_line(colour = "lightgrey",linetype = 3,linewidth = 0.25))+
  ylab("Normalized enrichment score\nfor TF target genes")+xlab("Transcription factors (TFs)")
###  Fig.5D end






