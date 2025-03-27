library(data.table)
TR4.df=fread("../Deconvolution/CIBERSORTx_TR4.txt")
TR4.df=as.data.frame(TR4.df)
cohort=readRDS("../Deconvolution/cohort.RDS")
TR4.df$cohort=cohort$group[match(TR4.df$Mixture,cohort$sample)]
TR4.df$cohort[TR4.df$cohort=="TARGET"]="TARGET-WT"
TR4.df$cohort[TR4.df$cohort=="Xinhua"]="Xinhua-WT"
TR4.df$cohort[TR4.df$cohort=="KIRP"]="TCGA-KIRP"
TR4.df$cohort=factor(TR4.df$cohort,levels = c("Xinhua-WT","TARGET-WT","TCGA-KICH","TCGA-KIRC","TCGA-KIRP"))
library(ggplot2);library(ggpubr);library(ggsci)
clinical.all=readRDS("../Clinical/clinical.all.RDS")
TR4.df=TR4.df[!(TR4.df$Mixture %in% c("WT1","WT2")),]
TR4.df1=TR4.df[TR4.df$cohort %in% c("TARGET-WT", "TCGA-KICH","TCGA-KIRC", "TCGA-KIRP") | TR4.df$Mixture %in% clinical.all$SampleID[clinical.all$Tissue=="T"],]
p=ggdensity(TR4.df1,x = "CD45",fill = "cohort",palette = "npg",xlab="Immune cell proportion",ylab = "Density")+
  facet_wrap(~cohort,scales = "free_y",ncol=1,strip.position = "right")+
  theme(legend.position = "none",strip.background = element_rect(fill = NA,colour = NA),
        strip.text.y = element_text(angle = 0),axis.text.y=element_blank(),
        axis.ticks.y = element_blank(),axis.line = element_blank()
  )
### Fig.6A
TR4.df.Xinhua=TR4.df[TR4.df$cohort=="Xinhua-WT",]
library(maftools)
laml=readRDS("../Mutations/laml.RDS")
laml.tmb=tmb(laml,captureSize = 60)
tmb.res=laml.tmb
tmb.res$Immune=TR4.df.Xinhua$CD45[match(tmb.res$Tumor_Sample_Barcode,TR4.df.Xinhua$Mixture)]
tmb.res.noNA=tmb.res[!is.na(tmb.res$Immune),]
ggscatter(tmb.res.noNA,x = "Immune",y = "total_perMB",xlab="Immune cell proportion",ylab="Tumor mutation burden (TMB)",color = "#008B457F",add = "reg.line",size = 0.5,cor.coef = T)+ylim(0,0.5)+
  theme(axis.line = element_blank(),panel.border = element_rect(fill = NA,colour = "black"))
### Fig 6B end


LM22.df=as.data.frame(fread("../Deconvolution/CIBERSORTx_LM22.txt"))
LM22.mat=as.matrix(LM22.df[,c(2:23)])
rownames(LM22.mat)=LM22.df$Mixture
LM22.mat=t(LM22.mat)
subtype=readRDS("../Figure-4/subtype.RDS")
LM22.mat=LM22.mat[,names(subtype)]
library(scales)
LM22.mat.adjust=sapply(colnames(LM22.mat),function(x){
  y=LM22.mat[,x]
  y=y*TR4.df.Xinhua[TR4.df.Xinhua$Mixture==x,"CD45"]
})
LM22.mat.avg=apply(LM22.mat.adjust,1,function(x){x=scales::rescale(x,to=c(0,1));tapply(x,subtype,mean)})
immune.pval=apply(LM22.mat.adjust,1,function(x){kruskal.test(x~subtype)$p.value})
immune.pval1=apply(LM22.mat.adjust,1,function(x){anova(lm(x~subtype))$"Pr(>F)"[1]})
pheatmap::pheatmap(t(LM22.mat.avg),cluster_cols = F,labels_row = paste(colnames(LM22.mat.avg)," (",scientific(immune.pval1,digits = 3),")",sep=""),
                   color = RColorBrewer::brewer.pal(n = 11,name = "RdBu")[11:1],treeheight_row = 0,labels_col = c("S1","S2","S3"),
                   border_color = NA,rect_gp=gpar(col="black"),border=T,border_gp=gpar(col="black"))
### Fig 6C
library(readxl)
cytokines=read_xlsx("cytokines.xlsx",sheet = 1,col_names = F)
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
WT.RNA.expr=readRDS("../RNA/WT.RNA.expr.RDS")
WT.protein.expr=readRDS("../Protein/WT.protein.expr.RDS")
subtype.deg.RNA=lapply(1:3,function(i){
  limma.deg(log2(WT.RNA.expr[,names(subtype)]+0.1),ifelse(subtype==i,1,0))
})
subtype.deg.pt=lapply(1:3,function(i){
  limma.deg(WT.protein.expr[,names(subtype)],ifelse(subtype==i,1,0))
})
subtype.deg.RNA.list=lapply(subtype.deg.RNA,function(x){
  rownames(x)[which(x$adj.P.Val<0.25 & x$logFC>log2(1) & x$P.Value<0.05)]
})
subtype.deg.pt.list=lapply(subtype.deg.pt,function(x){
  x$ID[which(x$adj.P.Val<0.25 & x$logFC>log2(1) & x$P.Value<0.05)]
})
cytokines.RNA=t(apply(as.matrix(cytokines),1,function(x){
  y=sapply(subtype.deg.RNA.list,function(t){sum(t==x[2])})
  if(sum(y)>0){
    y=c(x,paste(paste0("S",which(y>0)),collapse = ";"))
  }else{
    y=c(x,NA)
  }
  y
}))

cytokines.pt=t(apply(as.matrix(cytokines),1,function(x){
  y=sapply(subtype.deg.pt.list,function(t){sum(t==x[2])})
  if(sum(y)>0){
    y=c(x,paste(paste0("S",which(y>0)),collapse = ";"))
  }else{
    y=c(x,NA)
  }
  y
}))
df.subtype.deg=cytokines.RNA[!is.na(cytokines.RNA[,3]),]
rownames(df.subtype.deg)=df.subtype.deg[,2]
df.subtype.deg=as.data.frame(df.subtype.deg)
colnames(df.subtype.deg)=c("type","gene","subtype")
df.subtype.deg=df.subtype.deg[order(df.subtype.deg$type,df.subtype.deg$subtype),]
df.subtype.dep=cytokines.pt[!is.na(cytokines.pt[,3]),]
rownames(df.subtype.dep)=df.subtype.dep[,2]
df.subtype.dep=as.data.frame(df.subtype.dep)
colnames(df.subtype.dep)=c("type","gene","subtype")
df.subtype.dep=df.subtype.dep[order(df.subtype.dep$type,df.subtype.dep$subtype),]
anno.subtype=RColorBrewer::brewer.pal(n = 3,name = "Set2")
names(anno.subtype)=c("S1","S2","S3")
p1=ComplexHeatmap::pheatmap(t(apply(WT.RNA.expr[,names(subtype)],1,function(x){x=scales::rescale(log2(x+1),to=c(-3,3))}))[rownames(df.subtype.deg),order(subtype)],
                            cluster_rows = F,cluster_cols = F,color = RColorBrewer::brewer.pal(n = 11,name = "RdBu")[11:1],
                            annotation_col = data.frame(subtype=paste0("S",subtype[order(subtype)]),row.names = names(subtype)[order(subtype)]),
                            annotation_colors = list(subtype=anno.subtype),
                            row_split=df.subtype.deg$type,
                            row_title_rot=0,show_colnames = F,annotation_row = df.subtype.deg[,c(3,1)],column_split=paste0("S",subtype)[order(subtype)],
                            cellheight = 6,cellwidth = 3,
                            border_color = NA,border="black",fontsize_row = 7)  

p2=ComplexHeatmap::pheatmap(t(apply(WT.protein.expr[rownames(df.subtype.dep),names(subtype)],1,function(x){x=scales::rescale(x,to=c(-3,3))}))[,order(subtype)],
                            cluster_rows = F,cluster_cols = F,color = RColorBrewer::brewer.pal(n = 11,name = "RdBu")[11:1],
                            #annotation_col = data.frame(subtype=paste0("S",subtype[order(subtype)]),row.names = names(subtype)[order(subtype)]),
                            annotation_colors = list(subtype=anno.subtype),
                            row_split=df.subtype.dep$type,
                            row_title_rot=0,show_colnames = F,
                            annotation_row = df.subtype.dep[,c(3,1)],
                            column_split=paste0("S",subtype)[order(subtype)],
                            cellheight = 6,cellwidth = 3,border_color = NA,
                            border="black",fontsize_row = 7)
library(ComplexHeatmap)
p1 %v% p2
### Fig.6D


Co.inhibitors=c("CTLA4","TIGIT","BTLA","CD48","PDCD1","LAG3","CD274","HAVCR2","BTN2A2","LAIR2","BTN3A1","PDCD1LG2","BTN1A1","VTCN1","BTNL2")
WT.RNA.expr=readRDS("../RNA/WT.RNA.expr.RDS")
ICB.by.subtype=t(apply(log2(WT.RNA.expr+1)[Co.inhibitors,names(subtype)],1,function(x){
  x1=(x-mean(x))/sd(x)
  y=c(tapply(x1,subtype,mean),kruskal.test(x~subtype)$p.value)
  
}))
WT.protein.expr=readRDS("../Protein/WT.protein.expr.RDS")
ICB.by.subtype.pt=t(apply(WT.protein.expr[rownames(WT.protein.expr) %in% Co.inhibitors,names(subtype)],1,function(x){
  x1=(x-mean(x))/sd(x)
  y=c(tapply(x1,subtype,mean),kruskal.test(x~subtype)$p.value)
  
}))
signif2star=function(x){
  if(x>0.05)
  {y="NS"}else if(x<0.05 & x>0.01)
  {y="*"}else if(x<0.01 & x>0.001)
  {y="**"}else if(x<0.001){y="***"}
  y
}

library(ComplexHeatmap)
p2=ComplexHeatmap::pheatmap(t(ICB.by.subtype[,1:3]),cluster_rows = F,cluster_cols = F,
                            color = RColorBrewer::brewer.pal(n = 11,name = "RdBu")[11:1],
                            breaks = seq(-1,1,length.out=11),
                            labels_row = c("S1","S2","S3"),name="Cell abundance\nZ-score",cellwidth = 10,cellheight = 15,
                            labels_col = paste(rownames(ICB.by.subtype)," (",sapply(ICB.by.subtype[,4],signif2star),")",sep=""),
                            border_color = NA,border="black"
)  


p3=ComplexHeatmap::pheatmap(t(ICB.by.subtype.pt[,1:3]),cluster_rows = F,cluster_cols = F,
                            color = RColorBrewer::brewer.pal(n = 11,name = "RdBu")[11:1],
                            breaks = seq(-1,1,length.out=11),
                            labels_row = c("S1","S2","S3"),name="Cell abundance\nZ-score",cellwidth = 10,cellheight = 15,labels_col = paste(rownames(ICB.by.subtype.pt)," (",sapply(ICB.by.subtype.pt[,4],signif2star),")",sep=""),
                            border_color = NA,border="black"
)  
p2+p3
### Fig 6E








