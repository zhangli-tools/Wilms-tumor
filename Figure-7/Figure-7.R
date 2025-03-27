WT.RNA.expr=readRDS("../RNA/WT.RNA.expr.RDS")
WT.protein.expr=readRDS("../Protein/WT.protein.expr.RDS")
clinical.all=readRDS("../Clinical/clinical.all.RDS")
protein.tissue.type=clinical.all$Tissue[match(colnames(WT.protein.expr),clinical.all$SampleID)]
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
protein.results=limma.deg(WT.protein.expr,label = ifelse(protein.tissue.type=="T",1,0))
RNA.results=limma.deg(log2(0.1+WT.RNA.expr),label = ifelse(clinical.all$Tissue[match(colnames(WT.RNA.expr),clinical.all$SampleID)]=="T",1,0))
phosphoL6=readRDS("../Phosphoprotein/phosphoL6.RDS")
ph.result=limma.deg(phosphoL6,label = ifelse(clinical.all$Tissue[match(colnames(phosphoL6),clinical.all$SampleID)]=="T",1,0))
library(data.table)
drug.target.db=fread("../dgidb/interactions.tsv")
drug.db=fread("../dgidb/drugs.tsv")
WT.vs.NAT.deglist=list(rownames(RNA.results)[RNA.results$adj.P.Val<0.25 & RNA.results$P.Value<0.05 & RNA.results$logFC>log2(2)],
                       rownames(RNA.results)[RNA.results$adj.P.Val<0.25 & RNA.results$P.Value<0.05 & RNA.results$logFC< -log2(2)])
WT.vs.NAT.deplist=list(
  pt.up=protein.results$ID[which(protein.results$adj.P.Val<0.25 & protein.results$P.Value<0.05 & protein.results$logFC>log2(1.2))],
  pt.down=protein.results$ID[which(protein.results$adj.P.Val<0.25 & protein.results$P.Value<0.05& protein.results$logFC< -log2(1.2))]
)

WT.vs.NAT.dephlist.site=list(unique(ph.result$ID[which(ph.result$adj.P.Val<0.05 & ph.result$logFC>1)]),
                             unique(ph.result$ID[which(ph.result$adj.P.Val<0.05 & ph.result$logFC< -1)]))
WT.vs.NAT.dephlist.symbol=sapply(WT.vs.NAT.dephlist.site,function(x){unique(sapply(strsplit(x,";"),function(t){t[1]}))})
subtype=readRDS("../Figure-4/subtype.RDS")
subtype.deg.RNA=lapply(1:3,function(i){
  limma.deg(log2(0.1+WT.RNA.expr[,names(subtype)]),ifelse(subtype==i,1,0))
})
subtype.deg.protein=lapply(1:3,function(i){
  limma.deg(WT.protein.expr[,names(subtype)],ifelse(subtype==i,1,0))
})

subtype.shared.degs=lapply(1:3,function(i){
  RNA=subtype.deg.RNA[[i]]
  pt=subtype.deg.protein[[i]]
  x=rownames(RNA)[RNA$P.Value<0.05 & RNA$logFC>0.5]
  y=pt$ID[pt$P.Value<0.05 & pt$logFC>0]
  intersect(x,y)
})
WT.vs.NAT.drug.target.RNA.pt=intersect(intersect(WT.vs.NAT.deglist[[1]],WT.vs.NAT.deplist[[1]]),drug.target.db$gene_name[drug.target.db$interaction_types %in% c("inhibitor","antibody")])
S1.drug.target.RNA.pt=intersect(subtype.shared.degs[[1]],drug.target.db$gene_name[drug.target.db$interaction_types %in% c("inhibitor","antibody")])
#subtype.drug.target=sapply(subtype.shared.degs,intersect,drug.target.db$gene_name[drug.target.db$interaction_types %in% c("inhibitor","antibody")])
S1.drug.target.kinase=intersect(c("CDK1","CDK2","CDK5","MAPK1","MAPK8","GSK3A","GSK3B","CSNK1A1"),
                                drug.target.db$gene_name[drug.target.db$interaction_types %in% c("inhibitor","antibody")])
WT.vs.NAT.kinases=strsplit("ATM
AURKB
CDK1
CDK2
MAPK10
MAPK8
MAPK9
MAPKAPK2
PRKDC","\n")[[1]]

WT.kinase.drug.target=intersect(WT.vs.NAT.kinases,drug.target.db$gene_name[drug.target.db$interaction_types %in% c("inhibitor","antibody")])

shared.cis.amp1=readRDS("../Figure-3/shared.cis.amp1.RDS")
cis.CNV.drug.target=intersect(drug.target.db$gene_name[drug.target.db$interaction_types %in% c("inhibitor","antibody")],shared.cis.amp1)
laml=readRDS("../Mutations/laml.RDS")
mutation.drug.target=intersect(drug.target.db$gene_name[drug.target.db$interaction_types %in% c("inhibitor","antibody")], c("CTNNB1","WT1","AMER1","BRAF","TP53","CCNE1","PIK3CB"))
total.drugs=rbind(data.frame(gene=WT.vs.NAT.drug.target.RNA.pt,group="RNA/protein-WT"),
                  data.frame(gene=S1.drug.target.RNA.pt,group="RNA/protein-S1"),
                  data.frame(gene=WT.kinase.drug.target,group="Kinase-WT"),
                  data.frame(gene=S1.drug.target.kinase,group="Kinase-S1"),
                  cis.CNV=data.frame(gene=cis.CNV.drug.target,group="cis-CNA"),
                  mutation=data.frame(gene=mutation.drug.target,group="Mutation"))
total.drugs$group=factor(total.drugs$group,levels = unique(total.drugs$group)[c(9,8,1,5,2:4,6:7)])
total.drugs.sig=total.drugs[total.drugs$gene %in% names(which(table(total.drugs$gene)>=2)),]
library(ggpubr)
#x=load("/Volumes/Elements/0710/0825/TARGET-survival/drug.target.cox.out.Rdata")
#load("/Volumes/Elements/0710/0825/TARGET-survival/drug.target.cox.out.OS.Rdata")
target.mRNA=readRDS("../TARGET/gene.expr.RDS")
pdata=readRDS("../TARGET/clinical.RDS")
target.mRNA=target.mRNA[,intersect(colnames(target.mRNA),pdata$Sample_ID)]
mRNA.pdata=pdata[match(colnames(target.mRNA),pdata$Sample_ID),]

cox.out.EFS=t(sapply(unique(total.drugs.sig$gene),function(x){
  print(x)
  y=t(sapply(30:70,function(t){x1=ifelse(target.mRNA[x,]>quantile(target.mRNA[x,],t/100),1,0)
  c(t/100,summary(coxph(Surv(EFS_time,EFS)~x1,data = cbind(mRNA.pdata,x1)))$coefficients)
  }))
  y1=matrix(y[y[,2]>0,],ncol=6)
  if(nrow(y1)>1 & !is.null(nrow(y1)))
  {y=y1[which.min(y1[,6]),]}else{
    y=y[which.min(y[,6]),]
  }
  y
}))
cox.out.OS=t(sapply(unique(total.drugs.sig$gene),function(x){
  print(x)
  y=t(sapply(30:70,function(t){x1=ifelse(target.mRNA[x,]>quantile(target.mRNA[x,],t/100),1,0)
  c(t/100,summary(coxph(Surv(OS_time,OS)~x1,data = cbind(mRNA.pdata,x1)))$coefficients)
  }))
  y1=matrix(y[y[,2]>0,],ncol=6)
  if(nrow(y1)>1 & !is.null(nrow(y1)))
  {y=y1[which.min(y1[,6]),]}else{
    y=y[which.min(y[,6]),]
  }
  y
}))

#total.drugs.sig=total.drugs[total.drugs$gene %in% unique(total.drugs$gene[total.drugs$group %in% c("Mutation","cis-CNA","Kinase-WT")]),]
#total.drugs.sig=total.drugs[total.drugs$group %in% c("Mutation","cis-CNA","RNA/protein-WT","Kinase-WT","RNA/protein-S1","Kinase-S1"),]
#drug.genes=unique(total.drugs.sig$gene)
target.levels=names(sort(table(total.drugs.sig$gene),decreasing = F))
#target.levels=drug.genes[order(table(total.drugs.sig$gene)[drug.genes],1-apply(cbind(cox.out.EFS[,6],cox.out.OS[,6]),1,min),decreasing = F)]
#target.levels=
idx1=which(target.levels %in% c("EHMT2","TRIM24"))
idx2=idx1[c(2,1)]
target.levels[idx1]=target.levels[idx2]
total.drugs.sig$gene=factor(total.drugs.sig$gene,levels = target.levels)
#total.drugs.sig$pval=cox.out[as.matrix(total.drugs.sig$gene),6]
p1=ggscatter(data = total.drugs.sig,x = "group",y = "gene",color = "group",size=3)+
  theme(axis.text.x=element_text(angle = 45,hjust = 1),legend.position = "none",axis.text.y = element_text(size = 10),axis.line = element_blank(),panel.grid.major = element_line(colour = "lightgrey",linetype = 3,linewidth = 0.1),panel.border = element_rect(fill = NA,colour = "black"))+
  xlab("")+ylab("Drug targets")+scale_color_d3(alpha = 0.5)
forestdf=data.frame(labels=factor(levels(total.drugs.sig$gene),levels = levels(total.drugs.sig$gene)),
                    rr=2^cox.out.EFS[levels(total.drugs.sig$gene),2],
                    rrlow=2^(cox.out.EFS[levels(total.drugs.sig$gene),2]-qnorm(0.975)*cox.out.EFS[levels(total.drugs.sig$gene),4]),
                    rrhigh=2^(cox.out.EFS[levels(total.drugs.sig$gene),2]+qnorm(0.975)*cox.out.EFS[levels(total.drugs.sig$gene),4]))
p2 <- ggplot(forestdf, aes(x = rr, y = labels, xmin = rrlow, xmax = rrhigh,color="lightgrey"),size=1) +
  geom_hline(aes(yintercept = labels, colour = NA), size = 7) + 
  geom_pointrange(shape = 22, fill = pal_d3()(10)[4]) +
  geom_vline(xintercept = 1, linetype = 3) +
  xlab("Hazard ratio (EFS)") +
  ylab("Adjusted Relative Risk with 95% Confidence Interval") +
  theme_classic() +
  scale_colour_identity() +
  scale_y_discrete(limits = levels(rev(forestdf$labels))) +
  scale_x_log10(limits = c(0.35, 4), 
                breaks = c(0.35, 1, 2, 4), 
                labels = c( "0.35", "1", "2", "4"), expand = c(0,0)) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.line = element_blank(),panel.border = element_rect(fill = NA,colour = "black"))
forestdf1=data.frame(labels=factor(levels(total.drugs.sig$gene),levels = levels(total.drugs.sig$gene)),
                     rr=2^cox.out.OS[levels(total.drugs.sig$gene),2],
                     rrlow=2^(cox.out.OS[levels(total.drugs.sig$gene),2]-qnorm(0.975)*cox.out.OS[levels(total.drugs.sig$gene),4]),
                     rrhigh=2^(cox.out.OS[levels(total.drugs.sig$gene),2]+qnorm(0.975)*cox.out.OS[levels(total.drugs.sig$gene),4]))

p3 <- ggplot(forestdf1, aes(x = rr, y = labels, xmin = rrlow, xmax = rrhigh,color="lightgrey"),size=1) +
  geom_hline(aes(yintercept = labels, colour = NA), size = 7) + 
  geom_pointrange(shape = 22, fill = pal_d3()(10)[4]) +
  geom_vline(xintercept = 1, linetype = 3) +
  xlab("Hazard ratio (OS)") +
  ylab("Adjusted Relative Risk with 95% Confidence Interval") +
  theme_classic() +
  scale_colour_identity() +
  scale_y_discrete(limits = levels(rev(forestdf$labels))) +
  scale_x_log10(limits = c(0.25, 4), 
                breaks = c(0.25, 1, 2, 4), 
                labels = c( "0.25", "1", "2", "4"), expand = c(0,0)) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.line = element_blank(),panel.border = element_rect(fill = NA,colour = "black"))

p1+p2+p3
### Fig.7A


library(clusterProfiler)
db=msigdbr::msigdbr()
db.pathway=rbind(cbind(db$gs_description,db$gene_symbol)[grep("CP:",db$gs_subcat),],
                 cbind(db$gs_name,db$gene_symbol)[grep("HALLMARK",db$gs_name),])

obj=enricher(unique(total.drugs$gene),TERM2GENE=db.pathway,pAdjustMethod = "none",pvalueCutoff = 1)
path=c("Cell cycle","Regulation of TP53 Activity through Phosphorylation",
       "Wnt signaling pathway",
       "Regulation of TP53 Activity through Methylation",
       "MAPK signaling pathway","Transcriptional regulation by RUNX3",
       "Chromatin modifying enzymes","ErbB signaling pathway")
obj.res=obj@result
obj.res$geneID=sapply(strsplit(obj.res$geneID,"\\/"),function(x){
  paste(intersect(x,total.drugs.sig$gene),collapse = "/")
})
p=barplot(obj,showCategory = path)+scale_fill_material(palette = "orange",reverse = T)+
  geom_text(data =obj.res[obj.res$Description %in% path,] ,aes(x=0,label=geneID,hjust=0))
p=p+theme(panel.grid.major = element_blank(),panel.grid = element_blank())
p+scale_y_discrete(position = "right")
####Fig.7B
phosphoL6=readRDS("../Phosphoprotein/phosphoL6.RDS")
kinase.substrate.mat=readRDS("../Phosphoprotein/kinase.substrate.mat.RDS")
library(GSVA)
hallmark=cbind(db$gs_name,db$gene_symbol)[grep("HALLMARK",db$gs_name),]
cdk1and2.acitivity=GSVA::gsva(phosphoL6,sapply(c("CDK1","CDK2"),function(x){kinase.substrate.mat[kinase.substrate.mat[,1]==x,2]}),method="ssgsea")
RNA.e2f.act=gsva(log2(1+WT.RNA.expr),gset.idx.list = list(hallmark[hallmark[,1]=="HALLMARK_E2F_TARGETS",2]))
pt.e2f.act=gsva(WT.protein.expr,gset.idx.list = list(hallmark[hallmark[,1]=="HALLMARK_E2F_TARGETS",2]))
overlap.samples.across.multics=intersect(intersect(colnames(WT.RNA.expr),colnames(WT.protein.expr)),colnames(phosphoL6))
CDK.mat=rbind(cdk1and2.acitivity[,overlap.samples.across.multics],
              WT.protein.expr[c("CDK1","CDK2"),overlap.samples.across.multics],
              log2(WT.RNA.expr[c("CDK1","CDK2"),overlap.samples.across.multics]+1),
              phosphoL6[c("RB1;S249;","RB1;S37;"),overlap.samples.across.multics],
              RNA.e2f.act[,overlap.samples.across.multics],
              pt.e2f.act[,overlap.samples.across.multics])
rownames(CDK.mat)=c("CDK1 activity","CDK2 activity","CDK1 (protein)","CDK2 (protein)","CDK1 (RNA)","CDK1 (RNA)","RB1 pS249","RB1 pS37","E2F target (RNA)","E2F target (protein)")

CDK.mat.sample.anno=data.frame(Tissue=ifelse(clinical.all$Tissue[match(overlap.samples.across.multics,clinical.all$SampleID)]=="T","WT","NAT"),
                               Subtype=sapply(subtype[overlap.samples.across.multics],function(x){if(is.na(x)){y=x}else{y=paste0("S",x)};y}),row.names = overlap.samples.across.multics)[,2:1]
cls=CDK.mat.sample.anno$Subtype[order(CDK.mat.sample.anno$Tissue,CDK.mat.sample.anno$Subtype)]
cls[is.na(cls)]="NAT"
cls=factor(cls,levels = c("NAT","S1","S2","S3"))
ComplexHeatmap::pheatmap(CDK.mat[,order(clinical.all$Tissue[match(overlap.samples.across.multics,clinical.all$SampleID)],subtype[overlap.samples.across.multics])],
                         cluster_rows = F,scale="row",cluster_cols = F,
                         annotation_col = CDK.mat.sample.anno,
                         color = RColorBrewer::brewer.pal(n = 11,name = "RdBu")[11:1],
                         show_rownames = T,show_colnames = F,
                         annotation_colors = list(Tissue=c("WT"=pal_d3()(2)[2],"NAT"=pal_d3()(2)[1]),
                                                  Subtype=c(S1=pal_lancet()(3)[1],
                                                            S2=pal_lancet()(3)[2],
                                                            S3=pal_lancet()(3)[3])),
                         cellwidth = 7,cellheight = 12,border_gp = gpar(col = "black",size=0.5),
                         border_color = NA,column_split=cls,name="Z-score")



load("../CellLine/cell.line.fpkm.Rdata")

deg=limma.deg(log2(0.1+cell.line.fpkm),label = c(1,1,1,0,0,0))
up=rownames(deg)[which(deg$logFC>log2(1.1) & deg$P.Value<0.05)]
down=rownames(deg)[which(deg$P.Value<0.05 & deg$logFC< -log2(1.1))]
df=data.frame(group=rep(c("siEHMT2-up","siEHMT2-down"),c(length(up),length(down))),gene=c(up,down))

library(clusterProfiler)
obj=compareCluster(gene~group,fun = "enricher",data = df,TERM2GENE=db.pathway,qvalueCutoff=1,pvalueCutoff=1,pAdjustMethod="none")
obj@compareClusterResult=obj@compareClusterResult[obj@compareClusterResult$pvalue<0.05,]
#pathways=readxl::read_xlsx("../../Figures/cell_line/pathway.xlsx",sheet = "Sheet1")
library(enrichplot)
pathways=obj@compareClusterResult
show.pathways=c("Genes controlling nephrogenesis",unlist(tapply(pathways$Description,pathways$Cluster,head,n=10))[-10])
names(show.pathways)=NULL
edox2 <- pairwise_termsim(obj,showCategory = 1e10)
treeplot(edox2,nCluster=3,geneClusterPanel = "pie",hclust_method="complete",nWords=0,showCategory=show.pathways)+scale_fill_d3()
glist.show=strsplit(obj@compareClusterResult[obj@compareClusterResult$Description %in% show.pathways[c(1,2,4)],]$geneID,"\\/")
names(glist.show)=obj@compareClusterResult[obj@compareClusterResult$Description %in% show.pathways[c(1,2,4)],]$Description
glist.show
glist.show.df=data.frame(gene=unlist(glist.show),pathway=rep(names(glist.show),sapply(glist.show,length)))
glist.show.df=glist.show.df[order(glist.show.df$gene),]
glist.show.df.sort=glist.show.df[order(glist.show.df$pathway),]
ComplexHeatmap::pheatmap(t(apply(cell.line.fpkm[glist.show.df.sort$gene,6:1],1,function(x){x=log2(x+1);scales::rescale(x,to=c(-3,3))})),cluster_rows = T,cluster_cols = F,row_split=glist.show.df.sort$pathway,name="RNA abudance",row_title_rot=90,
                         annotation_col = data.frame(Group=rep(c("siEHMT2","siNC"),each=3),row.names = colnames(cell.line.fpkm)),annotation_colors = list(Group=c(siEHMT2="lightblue",siNC="brown1")),
                         show_colnames = F,color = RColorBrewer::brewer.pal(n = 11,name = "RdBu")[10:2],border_color = NA,border="black",row_title_side="right")




