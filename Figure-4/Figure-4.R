WT.protein.expr=readRDS("../Protein/WT.protein.expr.RDS")
WT.RNA.expr=readRDS("../RNA/WT.RNA.expr.RDS")
clinical.all=readRDS("../Clinical/clinical.all.RDS")
RNA.pt.shared=intersect(colnames(WT.RNA.expr),colnames(WT.protein.expr))
RNA.pt.shared.tumors=RNA.pt.shared[clinical.all$Tissue[match(RNA.pt.shared,clinical.all$SampleID)]=="T"]
library(CancerSubtypes)
data1=FSbyVar(log2(.1+WT.RNA.expr)[,RNA.pt.shared.tumors], cut.type="topk",value=5000)
data2=FSbyVar(WT.protein.expr[,RNA.pt.shared.tumors], cut.type="topk",value=4000)

WT=list(GeneExp=data1,miRNAExp=data2)
result=ExecuteSNF.CC(WT, clusterNum=3, K=55, alpha=0.5, t=20,
                      maxK = 5, pItem = 0.9,reps=500, 
                      title = "Subtypes", plot = "pdf",
                      finalLinkage ="average")
subtype=result$group
names(subtype)=colnames(data1)
phosphoL6=readRDS("../Phosphoprotein/phosphoL6.RDS")

#### Figure 4A
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
subtype.deg.RNA=lapply(1:3,function(i){
  limma.deg(log2(0.1+WT.RNA.expr[,names(subtype)]),ifelse(subtype==i,1,0))
})
subtype.deg.protein=lapply(1:3,function(i){
  limma.deg(WT.protein.expr[,names(subtype)],ifelse(subtype==i,1,0))
})
subtype.deg.ph=lapply(1:3,function(i){
  overlap.sample.subtype=intersect(names(subtype),colnames(phosphoL6))
  limma.deg(phosphoL6[,overlap.sample.subtype],ifelse(subtype[overlap.sample.subtype]==i,1,0))
})
subtype.deg.RNA.list=lapply(subtype.deg.RNA,function(x){rownames(x)[which(x$logFC>log2(1.5) & x$adj.P.Val<0.05)]})
subtype.deg.protein.list=lapply(subtype.deg.protein,function(x){x$ID[which(x$logFC>log2(1.2) & x$adj.P.Val<0.05)]})
subtype.deg.ph.list=lapply(subtype.deg.ph,function(x){
  y=x$ID[which(x$logFC>log2(1.5)  & x$P.Value<0.05)]
  y=unique(sapply(strsplit(y,";"),function(t){t[1]}))
})
df.subtype.all=rbind(data.frame(subtype=rep(c("S1","S2","S3"),sapply(subtype.deg.RNA.list,length)),gene=unlist(subtype.deg.RNA.list),data="RNA"),
                     data.frame(subtype=rep(c("S1","S2","S3"),sapply(subtype.deg.protein.list,length)),gene=unlist(subtype.deg.protein.list),data="Protein"),
                     data.frame(subtype=rep(c("S1","S2","S3"),sapply(subtype.deg.ph.list,length)),gene=unlist(subtype.deg.ph.list),data="Phosphoprotein"))
db=msigdbr::msigdbr()
db.pathway=rbind(cbind(db$gs_description,db$gene_symbol)[grep("CP:",db$gs_subcat),],
                 cbind(db$gs_name,db$gene_symbol)[grep("HALLMARK",db$gs_name),])

library(clusterProfiler)
obj.subtype.all=compareCluster(gene~subtype+data,fun = "enricher",data = df.subtype.all,TERM2GENE=db.pathway,pAdjustMethod="none",qvalueCutoff=1)
obj.subtype.all@compareClusterResult=obj.subtype.all@compareClusterResult[obj.subtype.all@compareClusterResult$qvalue<0.25,]
shared.pathways=lapply(c("S1","S2","S3"),function(x){
  names(which(table(obj.subtype.all@compareClusterResult$Description[obj.subtype.all@compareClusterResult$subtype==x & obj.subtype.all@compareClusterResult$data %in% c("RNA","Protein")])==2))
})
RNA.scaledata=t(apply(WT.RNA.expr[unlist(subtype.deg.RNA.list),names(subtype)],1,function(x){x=log2(0.1+x);y=scales::rescale(x,to=c(-3,3))}))
protein.scaledata=t(apply(WT.protein.expr[unlist(subtype.deg.protein.list),names(subtype)],1,function(x){x=log2(x);y=scales::rescale(x,to=c(-3,3))}))
library(ComplexHeatmap)
library(ggsci)
library(circlize)
d3.cols=pal_d3()(10)
anno.for.heatmap=readRDS("anno.for.heatmap.RDS")
annot.colors=list(Gender=c(Male=d3.cols[1],Female=d3.cols[7]),
                  Age=colorRamp2(c(0, 13), c("white", "brown")),
                  Risk=c("Intermediate risk"=d3.cols[2],"High risk"=d3.cols[4]),
                  Stage=c("1"=gplots::colorpanel(5,"white",d3.cols[4])[1],
                          "2"=gplots::colorpanel(5,"white",d3.cols[4])[2],
                          "3"=gplots::colorpanel(5,"white",d3.cols[4])[3],
                          "4"=gplots::colorpanel(5,"white",d3.cols[4])[4],
                          "5"=gplots::colorpanel(5,"white",d3.cols[4])[5]),
                  Matastasis=c(Yes="black",No="grey"),
                  Site=c(Left=d3.cols[1],Right=d3.cols[3],Bilateral=d3.cols[2]),
                  Event=c("Event-free"="grey","Recurrence"="grey50",Dead="black")
)
mut.cols=lapply(colnames(anno.for.heatmap)[8:13],function(x){y=c(Yes="black",No="grey")})
names(mut.cols)=colnames(anno.for.heatmap)[8:13]
annot.colors=c(annot.colors,mut.cols)
ha1 = HeatmapAnnotation(df =anno.for.heatmap[order(subtype),], col = annot.colors,na_col="white",border = T,gp=gpar(col="lightgrey"),
                        annotation_name_side="right")
RNA.heatmap.col=gplots::bluered(11)
pt.heatmap.col=RColorBrewer::brewer.pal(n = 11,name = "RdBu")[11:1]

p1=Heatmap(RNA.scaledata[,order(subtype)],cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,use_raster = T,name = "RNA",
           row_split = rep(c("S1","S2","S3"),sapply(subtype.deg.RNA.list,length)),column_split = factor(sort(paste0("S",subtype)),levels=c("S1","S2","S3")),height =unit(8,"cm"),
           col = pt.heatmap.col,raster_by_magick = F,border = "black")
p2=Heatmap(protein.scaledata[,order(subtype)],cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,use_raster = T,name = "Protein",
           row_split = rep(c("S1","S2","S3"),sapply(subtype.deg.protein.list,length)),column_split = factor(sort(paste0("S",subtype)),levels=c("S1","S2","S3")),
           heatmap_height = unit(6,"cm"),col=pt.heatmap.col,raster_by_magick = F,
           border = "black")
ht_list = ha1 %v% p1 %v% p2
heatmap2ggplot=function(heatmap){
  grob = grid::grid.grabExpr(ComplexHeatmap::draw(heatmap)) 
  p=ggpubr::as_ggplot(grob)
  return(p)
}
p=heatmap2ggplot(draw(ht_list))
#### Figure 4A end
HE.data=readRDS("../Clinical/HE.data.RDS")
HE.tab=table(sapply(strsplit(HE.data$HE.0524,""),function(x){x=x[1];ifelse(x %in% LETTERS,x,NA)}),subtype)
library(ggsci)
pval=sapply(c(1,3,2),function(i){
  chisq.test(rbind(HE.tab[i,],apply(HE.tab[-i,],2,sum)))$p.value
})
par(mai=c(0.5,0.5,0.5,0.5))
m=barplot(t(apply(HE.tab[c(1,3,2),],2,function(x){x/sum(x)})),beside = T,col =pal_lancet(alpha = 0.5)(3) ,ylim=c(0,1.1),names.arg = c("Blastemal","Stromal","Epithelial"),ylab="Proportion (%)",axes = F)
rect(xleft = 0.75,ybottom = 0,xright = max(m)+0.85,ytop = 1.1)
text(m[c(1,5,9)],t(apply(HE.tab[c(1,3,2),],2,function(x){x/sum(x)}))[c(1,5,9)]+0.025,sapply(pval,function(p){if(p>0.05){y="NS"}else if(p>0.01){y="*"}else if(p>0.001){y="**"}else{y="***"};y}))
axis(side=2,at = seq(0,0.8,by=0.2),labels = seq(0,80,by=20),pos = 0.75)
legend("top",legend = c("S1","S2","S3"),fill=pal_lancet(alpha = 0.5)(3),box.col = NA,horiz = T,bg = NA)
## Fig.4B end

###
estimate.score=fread("../Deconvolution/estimate.out.gct",skip=2)[,-1]
estimate.score=as.matrix(estimate.score[,-1],rownames = estimate.score$Description)
library(ggpubr)

p1=ggboxplot(data.frame(score=estimate.score[1,names(subtype)],subtype=paste0("S",subtype)),x = "subtype",y = "score",color = "subtype",add = "jitter")+stat_compare_means(comparisons = list(c("S1","S2"),c("S1","S3"),c("S2","S3")),label = "p.signif",method = "wilcox.test")+scale_color_lancet(alpha = 0.5)+theme(axis.line = element_blank(),panel.border = element_rect(fill = NA,colour = "black"),legend.position = "none")+xlab("")+ylab("Stromal score")
p2=ggboxplot(data.frame(score=estimate.score[2,names(subtype)],subtype=paste0("S",subtype)),x = "subtype",y = "score",color = "subtype",add = "jitter")+stat_compare_means(comparisons = list(c("S1","S2"),c("S1","S3"),c("S2","S3")),label = "p.signif",method = "wilcox.test")+scale_color_lancet(alpha = 0.5)+theme(axis.line = element_blank(),panel.border = element_rect(fill = NA,colour = "black"),legend.position = "none")+xlab("")+ylab("Immune score")
p1+p2
### Fig 4C end


subtype.deg.overlap=lapply(1:3,function(i){
  RNA=subtype.deg.RNA[[i]]
  pt=subtype.deg.protein[[i]]
  y=intersect(rownames(RNA)[which(RNA$logFC>log2(1) & RNA$P.Value<0.05)],
              pt$ID[which(pt$logFC>0  & pt$P.Value<0.05)])
  
})

calc.value.zscore=function(gene){
  if(gene %in% rownames(WT.RNA.expr))
  {RNA=log2(0.1+WT.RNA.expr[gene,names(subtype)])
  RNA=(RNA-mean(RNA))/sd(RNA)}else{
    RNA=rep(NA,length(subtype))
  }
  if(gene %in% rownames(WT.protein.expr))
  {pt=WT.protein.expr[gene,names(subtype)]
  pt=(pt-mean(pt))/sd(pt)}else{
    pt=rep(NA,length(subtype))
  }
  y=rbind(tapply(RNA,subtype,mean),tapply(pt,subtype,mean))
  colnames(y)=paste0("S",1:3)
  # ComplexHeatmap::Heatmap(y,cluster_columns = F,cluster_rows = F,rect_gp = gpar(col = "white", lwd = 2),heatmap_height=unit(2.5,"cm"),heatmap_width =unit(3,"cm"),show_heatmap_legend = F,column_title=gene,col = RColorBrewer::brewer.pal(    n = 11,name = "RdBu")[11:1])
  col.fun=circlize::colorRamp2(seq(-1,1,length.out=100), colors = gplots::colorpanel(100,low =  RColorBrewer::brewer.pal(n = 11,name = "RdBu")[10],mid =  RColorBrewer::brewer.pal(n = 11,name = "RdBu")[6],
                                                                                     high =  RColorBrewer::brewer.pal(n = 11,name = "RdBu")[2]))
ComplexHeatmap::Heatmap(y,cluster_columns = F,cluster_rows = F,rect_gp = gpar(col = "white", lwd = 2),heatmap_height
                        =unit(2.5,"cm"),heatmap_width =unit(3,"cm"),show_heatmap_legend = F,column_title=gene,col=col.fun,border = T,
                        border_gp= gpar(col="black"),show_column_names = F)
}

S1=subtype.deg.overlap[[1]]
S2=subtype.deg.overlap[[2]]
S3=subtype.deg.overlap[[3]]
dir.create("gene.heatmap/S1",recursive = T)
dir.create("gene.heatmap/S2",recursive = T)
dir.create("gene.heatmap/S3",recursive = T)
for(i in S1){
  pdf(file=paste0("gene.heatmap/S1/",i,".pdf"),width = 5,height = 5)
  p=calc.value.zscore(i)
  print(p)
  dev.off()
}
for(i in S2){
  pdf(file=paste0("gene.heatmap/S2/",i,".pdf"),width = 5,height = 5)
  p=calc.value.zscore(i)
  print(p)
  dev.off()
}
for(i in S3){
  pdf(file=paste0("gene.heatmap/S3/",i,".pdf"),width = 5,height = 5)
  p=calc.value.zscore(i)
  print(p)
  dev.off()
}
### Fig. 4E end



subtype.deg.RNA.list=lapply(subtype.deg.RNA,function(x){rownames(x)[which(x$logFC>log2(1.5) & x$adj.P.Val<0.25 & x$P.Value<0.05)]})
subtype.deg.protein.list=lapply(subtype.deg.protein,function(x){x$ID[which(x$logFC>log2(1) & x$adj.P.Val<0.25 & x$P.Value<0.05)]})


subtype.sigatures=lapply(1:3,function(i){
  intersect(subtype.deg.RNA.list[[i]],
            subtype.deg.protein.list[[i]])
})
source("../TARGET/NTP.R")
template=matrix(0,nrow=length(unique(unlist(subtype.sigatures))),ncol=3,dimnames = list(unique(unlist(subtype.sigatures)),c("S1","S2","S3")))
for(i in 1:3){template[subtype.sigatures[[i]],i]=1}
target.expr=readRDS("../TARGET/gene.expr.RDS")
expr.genes=rownames(target.expr)[apply(target.expr>1,1,sum)>0]
template.hc=template[rownames(template) %in% expr.genes,]
template.hc1=template.hc[unlist(tapply(rownames(template.hc),apply(template.hc,1,which.max),head,n=34)),]
pred.label=NTP_Fuc(log2(target.expr[rownames(template.hc1),]+1),template = template.hc1)
pred.label.target=pred.label
library(survival)
library(survminer)
#pred.label=pred.label[as.numeric(pred.label[,"BH.FDR"])<0.25,]
df.clinical=readRDS("../TARGET/clinical.RDS")
df.clinical$subtype=paste0("S",pred.label[match(df.clinical$Sample_ID,pred.label[,1]),2])
df.clinical$FDR=as.numeric(pred.label[match(df.clinical$Sample_ID,pred.label[,1]),"BH.FDR"])
df.clinical1=df.clinical[df.clinical$FDR<0.05,]
fit.efs=survfit(Surv(EFS_time,EFS)~subtype,data = df.clinical1)
p1=ggsurvplot(fit.efs,data = df.clinical1,pval = T,palette = "aaas",legend.title="Subtype",legend.labs=c("S1","S2","S3"),risk.table = T,xlab="Event-free survival (days)")
p1
fit.os=survfit(Surv(OS_time,OS)~subtype,data = df.clinical1)
p2=ggsurvplot(fit.os,data = df.clinical,pval = T,palette = "aaas",legend.title="Subtype",legend.labs=c("S1","S2","S3"),risk.table = T,xlab="Overall survival (days)")
p2
ggarrange(plotlist = list(p1$plot+theme(axis.line = element_blank(),panel.border = element_rect(fill = NA,colour = "black")),
                          p2$plot+theme(axis.line = element_blank(),panel.border = element_rect(fill = NA,colour = "black"))),nrow=1,ncol=2,common.legend = T,legend = "bottom")




############
kinase.substrate.mat=readRDS("../Phosphoprotein/kinase.substrate.mat.RDS")
subtype.deg.ph=lapply(1:3,function(i){
  overlap.sample.subtype=intersect(names(subtype),colnames(phosphoL6))
  limma.deg(phosphoL6[,overlap.sample.subtype],ifelse(subtype[overlap.sample.subtype]==i,1,0))
})
stat.ph.list=lapply(1:3,function(i){
  x=subtype.deg.ph[[i]]
  stat=x$t
  y=subtype.deg.protein.list[[i]]
  names(stat)=x$ID
  stat=stat[!sapply(strsplit(names(stat),";"),function(t){t[1]}) %in% y]
  y=GSEA(sort(stat,decreasing = T),TERM2GENE=kinase.substrate.mat,minGSSize = 4,pAdjustMethod = "none",maxGSSize = 1000)
})
subtype.ph.gsea=lapply(stat.ph.list,function(x){x@result[x@result$NES>0 & x@result$pvalue<0.05 & x@result$qvalue<0.15,]})
subtype.ph.gsea=lapply(1:3,function(i){
  x=subtype.ph.gsea[[i]]
  total.sites=lapply(subtype.deg.ph,function(x){
    y=x$ID[which(x$logFC>log2(2)  & x$P.Value<0.1)]
  })[[i]]
  total=sapply(strsplit(x$core_enrichment,"\\/"),length)
  sites=sapply(strsplit(x$core_enrichment,"\\/"),function(t){paste(intersect(t,total.sites),collapse = "/")})
  num=sapply(strsplit(x$core_enrichment,"\\/"),function(t){length(intersect(t,total.sites))})
  x$sig.sites=sites
  x$count=num
  x$total=total
  y=x
})
protein.results=limma.deg(WT.protein.expr,ifelse(clinical.all$Tissue[match(colnames(WT.protein.expr),clinical.all$SampleID)]=="T",1,0))

subtype.ph.gsea.sig=lapply(subtype.ph.gsea,function(x){
  x[x$count>=3,]
})

S1.kinase.substrate=strsplit(subtype.ph.gsea.sig[[1]]$sig.sites,"\\/")
names(S1.kinase.substrate)=subtype.ph.gsea.sig[[1]]$ID

S2.kinase.substrate=strsplit(subtype.ph.gsea.sig[[2]]$sig.sites,"\\/")
names(S2.kinase.substrate)=subtype.ph.gsea.sig[[2]]$ID

S3.kinase.substrate=strsplit(subtype.ph.gsea.sig[[3]]$sig.sites,"\\/")
names(S3.kinase.substrate)=subtype.ph.gsea.sig[[3]]$ID

tmp.expr1=t(apply(phosphoL6[unlist(c(S1.kinase.substrate,S2.kinase.substrate,S3.kinase.substrate)),clinical.all$Tissue
                            [match(colnames(phosphoL6),clinical.all$SampleID)]=="T"],1,function(x){
                              scales::rescale(x,to=c(-3,3))
                            }))
rownames(tmp.expr1)=paste0("s",1:nrow(tmp.expr1))
anno.sites1=data.frame(kinase=rep(names(c(S1.kinase.substrate,S2.kinase.substrate,S3.kinase.substrate)),sapply(c(S1.kinase.substrate,S2.kinase.substrate,S3.kinase.substrate),length)),row.names = rownames(tmp.expr1))
cols.kinase1=pal_d3(palette = "category20")(20)[3:(length(unique(anno.sites1$kinase))+2)]
names(cols.kinase1)=unique(anno.sites1$kinase)
ha2=rowAnnotation(df = anno.sites1,col=list(kinase=cols.kinase1),show_legend=F)
lbs1=unlist(c(S1.kinase.substrate,S2.kinase.substrate,S3.kinase.substrate))
lbs1=gsub(";","_",gsub(";$","",lbs1))
ha3=HeatmapAnnotation(Subtype=paste0("S",subtype[colnames(tmp.expr1)][order(subtype[colnames(tmp.expr1)])]),col = list(Subtype=c(S1=pal_d3()(10)[3],S2=pal_d3()(10)[4],S3=pal_d3()(10)[5])),show_legend = F)
tmp.expr1=tmp.expr1[,!(colnames(tmp.expr1) %in% c("WT1","WT2"))]
pdf(file="ph.heatmap.subtype.pdf",width = 7,height = 7)
Heatmap(tmp.expr1[,order(subtype[colnames(tmp.expr1)])],cluster_columns = F,cluster_rows = F,show_column_names = F,left_annotation = ha2,row_labels = lbs1,row_names_gp = gpar(fontsize = 8),top_annotation = ha3,row_split = factor(anno.sites1$kinase,levels=unique(anno.sites1$kinase)),row_title_rot = c(0),name = "Phosphorylated levels",column_split = paste0("S",subtype[colnames(tmp.expr1)])[order(subtype[colnames(tmp.expr1)])],col = RColorBrewer::brewer.pal(n = 11,name = "RdBu")[11:1],border = "black")
dev.off()
### Fig.4G end
GSE31403.expr.mat=readRDS("../GSE31403/gene.expr.RDS")
pred.label=NTP_Fuc(GSE31403.expr.mat[rownames(GSE31403.expr.mat) %in% rownames(template),],template = template[rownames(GSE31403.expr.mat[rownames(GSE31403.expr.mat) %in% rownames(template),]),])
pred.label.GSE31403=pred.label[as.numeric(pred.label[,"BH.FDR"])<0.25,]
GSE31403.clinical=readRDS("../GSE31403/clinical.RDS")

chi.p=chisq.test(table(GSE31403.clinical[,1],pred.label[,2]))$p.value
barplot(apply(table(GSE31403.clinical[,1],pred.label[,2]),2,function(x){x/sum(x)}),col=pal_d3()(6),
        main=paste("Chi-square test, P =",round(chi.p,3),"\nGSE31403"),xlim=c(0,4.5),name=c("S1","S2","S3"))
legend("right",legend = c("I","II","III","IV","V"),fill = pal_d3()(5),box.col = NA,bg = NA)
### Fig. 4I end

