library(ggpubr)
library(ggVennDiagram)
library(ggsci)
library(ggplot2)
library(clusterProfiler)
library(xlsx)
library(GSVA)
library(ComplexHeatmap)
library(scales)
library(survminer)
library(survival)
library(enrichplot)

## Load datasets
clinical.all=readRDS("../Clinical/clinical.all.RDS")
WT.RNA.expr=readRDS("../RNA/WT.RNA.expr.RDS")
WT.protein.expr=readRDS("../Protein/WT.protein.expr.RDS")
laml=readRDS("../Mutations/laml.RDS")
phosphoL6=readRDS("../Phosphoprotein/phosphoL6.RDS")

### Figure 2AB
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
RNA.results=limma.deg(log2(0.1+WT.RNA.expr),ifelse(clinical.all$Tissue=="T",1,0)[match(colnames(WT.RNA.expr),clinical.all$SampleID)])
protein.results=limma.deg(WT.protein.expr,label = ifelse(clinical.all$Tissue=="T",1,0)[match(colnames(WT.protein.expr),clinical.all$SampleID)])
WT.vs.NAT.deplist=list(
  pt.up=protein.results$ID[which(protein.results$adj.P.Val<0.05 & protein.results$logFC>log2(1.2))],
  pt.down=protein.results$ID[which(protein.results$adj.P.Val<0.05 & protein.results$logFC< -log2(1.2))]
)
ph.result=limma.deg(phosphoL6,label = ifelse(clinical.all$Tissue=="T",1,0)[match(colnames(phosphoL6),clinical.all$SampleID)])

WT.vs.NAT.dephlist.site=list(unique(ph.result$ID[which(ph.result$adj.P.Val<0.05 & ph.result$logFC>1)]),
                             unique(ph.result$ID[which(ph.result$adj.P.Val<0.05 & ph.result$logFC< -1)]))
WT.vs.NAT.dephlist.symbol=sapply(WT.vs.NAT.dephlist.site,function(x){unique(sapply(strsplit(x,";"),function(t){t[1]}))})

WT.vs.NAT.deglist=list(rownames(RNA.results)[RNA.results$adj.P.Val<0.05 & RNA.results$P.Value<0.05 & RNA.results$logFC>log2(2)],
                       rownames(RNA.results)[RNA.results$adj.P.Val<0.05 & RNA.results$P.Value<0.05 & RNA.results$logFC< -log2(2)])
library(ggVennDiagram)
p.up=ggVennDiagram(list(RNA=WT.vs.NAT.deglist[[1]],
                        Protein=WT.vs.NAT.deplist[[1]],
                        Phosphoprotein=WT.vs.NAT.dephlist.symbol[[1]])[c(2,1,3)],label="both",label_percent_digit=2,label_alpha = 0)

p.up=p.up+scale_fill_material(palette = "orange")+scale_colour_manual(values = rep("white",3))
p.up

p.down=ggVennDiagram(list(RNA=WT.vs.NAT.deglist[[2]],
                          Protein=WT.vs.NAT.deplist[[2]],
                          Phosphoprotein=WT.vs.NAT.dephlist.symbol[[2]])[c(2,1,3)],label="both",label_percent_digit=2,edge_color="black",edge_size = 0.5,label_alpha = 0)
p.down=p.down+scale_fill_material(palette = "blue")+scale_colour_manual(values = rep("white",3))
p.down
ggarrange(plotlist = list(p.up,p.down),nrow=1,ncol=2)
### Figure 2AB end
df.up=data.frame(gene=unlist(list(RNA=WT.vs.NAT.deglist[[1]],
                                  Protein=WT.vs.NAT.deplist[[1]],
                                  Phosphoprotein=WT.vs.NAT.dephlist.symbol[[1]])),
                 group=rep(c("RNA","Protein","Phosphoprotein"),sapply(list(RNA=WT.vs.NAT.deglist[[1]],
                                                                           Protein=WT.vs.NAT.deplist[[1]],
                                                                           Phosphoprotein=WT.vs.NAT.dephlist.symbol[[1]]),length)))
df.down=data.frame(gene=unlist(list(RNA=WT.vs.NAT.deglist[[2]],
                                    Protein=WT.vs.NAT.deplist[[2]],
                                    Phosphoprotein=WT.vs.NAT.dephlist.symbol[[2]])),
                   group=rep(c("RNA","Protein","Phosphoprotein"),sapply(list(RNA=WT.vs.NAT.deglist[[2]],
                                                                             Protein=WT.vs.NAT.deplist[[2]],
                                                                             Phosphoprotein=WT.vs.NAT.dephlist.symbol[[2]]),length)))
### Figure 2D
library(clusterProfiler)
db=msigdbr::msigdbr()
db.pathway=rbind(cbind(db$gs_description,db$gene_symbol)[grep("CP:",db$gs_subcat),],
                 cbind(db$gs_name,db$gene_symbol)[grep("HALLMARK",db$gs_name),])

obj.up=compareCluster(gene~group,fun = "enricher",data = df.up,TERM2GENE=db.pathway)
obj.down=compareCluster(gene~group,fun = "enricher",data = df.down,TERM2GENE=db.pathway)

path.list=list(RNA.processing=c("mRNA Processing","mRNA Splicing","RNA Polymerase II Transcription Termination","Processing of Intronless Pre-mRNAs"),
               WNT=c("Signaling by WNT","TCF dependent signaling in response to WNT","Repression of WNT target genes"),
               DNA.damage.repair=c("DNA Damage Response","DNA Damage Bypass","G2/M DNA damage checkpoint","SUMOylation of DNA damage response and repair proteins"),
               Cell.proliferation.and.cell.cycle=c("DNA Replication","Mitotic Metaphase and Anaphase","Cell cycle","Retinoblastoma Gene in Cancer"),
               embryonic.development=c("Development of ureteric collection system","Genes controlling nephrogenesis","Nephrogenesis"),
               epigenetic.regulation=c("SUMOylation of DNA methylation proteins","DNA methylation","Epigenetic regulation of gene expression","ERCC6 (CSB) and EHMT2 (G9a) positively regulate rRNA expression"),
               Renal.function=c("Proximal tubule transport","Proximal tubule bicarbonate reclamation"),
               metabolism=c("ABC transporters","SLC transporter disorders","Disorders of transmembrane transporters","Amino Acid metabolism","Gluconeogenesis","Fatty acid metabolism","Oxidative phosphorylation","Glycolysis / Gluconeogenesis"),
               NOTCH=c("Notch-HLH transcription pathway","Signaling by NOTCH1","Canonical  and Non-canonical Notch signaling"),
               TME=c("Complement cascade","FCGR activation","L1CAM interactions","MHC class II antigen presentation","Interferon gamma signaling","CD22 mediated BCR regulation","Focal Adhesion","Antigen processing and presentation","Cell-Cell communication","Leukocyte transendothelial migration"))
xx.up<- pairwise_termsim(obj.up,showCategory = unlist(path.list[c(1:4,5,6,9)]))  
xx.down<- pairwise_termsim(obj.down,showCategory = unlist(path.list[c(7:8,10)]))  

p.up=emapplot(xx.up,showCategory=unlist(path.list[c(1:4,5,6,9)]),label_style="ggforce",shadowtext=F,min_edge=0.1)+scale_fill_d3()
p.down=emapplot(xx.down,showCategory=unlist(path.list[c(7:8,10)]),label_style="ggforce",shadowtext=F,min_edge=0.05)+scale_fill_d3()
ggarrange(plotlist = list(p.up,p.down),nrow=1,ncol=2,common.legend = T,legend="right")
### Figure 2D end

### Figure 2C
boxplot.RNA.tn=function(gene){
  RNA.tissue.type=clinical.all$Tissue[match(colnames(WT.RNA.expr),clinical.all$SampleID)]
  if(length(gene)==1)
  {
    df=data.frame(score=log2(WT.RNA.expr[gene,]+1),group=factor(ifelse(RNA.tissue.type=="T","WT","NAT"),levels = c("NAT","WT")),gene=rep(gene,each=length(RNA.tissue.type)))
    df$score=(df$score-mean(df$score))/sd(df$score)
  }else{
    df=data.frame(score=c(apply(log2(1+WT.RNA.expr[gene,]),1,function(x){
      (x-mean(x))/sd(x)
    })),
    group=factor(ifelse(RNA.tissue.type=="T","WT","NAT"),levels = c("NAT","WT")),
    gene=factor(rep(gene,each=length(RNA.tissue.type)),levels = gene) )
  }
  
  cmps=matrix(c("WT","NAT"),nrow=2)
  ggboxplot(data = df,x = "group",y = "score",add.params=list(size=0.5),color = "group",palette = "d3",add = "jitter",error.plot="mean_linerange",facet.by = "gene",nrow=1)+
    stat_compare_means(comparisons = lapply(1:ncol(cmps),function(j){cmps[,j]}),method = "wilcox.test" )+ylim(min(df$score),max(df$score)+(max(df$score)-min(df$score))*0.7)+xlab("")+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5),
          strip.background = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.text = element_text(size = 14))+ylab("RNA abundance")
}
boxplot.protein.tn=function(gene){
  protein.tissue.type=clinical.all$Tissue[match(colnames(WT.protein.expr),clinical.all$SampleID)]
  if(length(gene)==1)
  {
    df=data.frame(score=WT.protein.expr[gene,],group=factor(ifelse(protein.tissue.type=="T","WT","NAT"),levels = c("NAT","WT")),gene=rep(gene,each=length(protein.tissue.type)))
    df$score=(df$score-mean(df$score))/sd(df$score)
  }else{
    df=data.frame(score=c(apply(WT.protein.expr[gene,],1,function(x){
      (x-mean(x))/sd(x)
    })),
    group=factor(ifelse(protein.tissue.type=="T","WT","NAT"),levels = c("NAT","WT")),
    gene=factor(rep(gene,each=length(protein.tissue.type)),levels = gene) )
  }
  
  cmps=matrix(c("WT","NAT"),nrow=2)
  ggboxplot(data = df,x = "group",y = "score",color = "group",add.params=list(size=0.5),palette = "d3",add = "jitter",error.plot="mean_linerange",facet.by = "gene",nrow=1)+
    stat_compare_means(comparisons = lapply(1:ncol(cmps),function(j){cmps[,j]}),method = "wilcox.test" )+ylim(min(df$score),max(df$score)+(max(df$score)-min(df$score))*0.7)+xlab("")+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5),
          strip.background = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.text = element_text(size = 14))+ylab("Protein abundance")
}
p1=boxplot.RNA.tn(c("WT1","VIM","NCAM1","GPC3","DES"))+ylim(-2.2,4.5)
p2=boxplot.protein.tn(c("WT1","VIM","NCAM1","GPC3","DES"))+ylim(-2.2,3)
ggarrange(plotlist = list(p1,p2),nrow=2,ncol=1)
### Figure 2C end

### Figure 2E
pathway.genes=readxl::read_xlsx("pathways-selected.xlsx",sheet = "Sheet2")
overlap.up=intersect(WT.vs.NAT.deglist[[1]],WT.vs.NAT.deplist[[1]])
overlap.down=intersect(WT.vs.NAT.deglist[[2]],WT.vs.NAT.deplist[[2]])
pathway.genes=pathway.genes[pathway.genes$Gene %in% c(overlap.up,overlap.down),]
pathway.genes=pathway.genes[order(match(pathway.genes$Path,
                                        c("Cell proliferation","DNA damage repair","Epigenetic regulation",
                                          "NOTCH signaling","RNA processing","Wnt signaling","Stroma",
                                          "Metabolism","Renal function","Immune"))),]
library(ComplexHeatmap)
library(ggsci)
ha.RNA=HeatmapAnnotation(df = data.frame(Tissue=ifelse(clinical.all$Tissue=="T","WT","NAT")[match(colnames(WT.RNA.expr),clinical.all$SampleID)][order(clinical.all$Tissue[match(colnames(WT.RNA.expr),clinical.all$SampleID)])]),col = list(Tissue=c(WT=pal_d3()(2)[2],NAT=pal_d3()(2)[1])),border = T,show_annotation_name = F)
ods.RNA=order(clinical.all$Tissue[match(colnames(WT.RNA.expr),clinical.all$SampleID)],
              apply(t(apply(log2(1+WT.RNA.expr)[pathway.genes$Gene[1:61],],1,function(x){scales::rescale(x,to=c(-3,3))})),2,mean),
              apply(t(apply(log2(1+WT.RNA.expr)[pathway.genes$Gene[-c(1:61)],],1,function(x){scales::rescale(x,to=c(-3,3))})),2,mean))
p.RNA.path.genes=Heatmap(t(apply(log2(1+WT.RNA.expr[pathway.genes$Gene,ods.RNA]),1,function(x){scales::rescale(x,to=c(-3,3))})),show_column_names = F,row_names_gp = gpar(cex=0.5),cluster_rows = F,top_annotation = ha.RNA,cluster_columns = F,name='RNA abundance',row_split = factor(pathway.genes$Path,levels = unique(pathway.genes$Path)),row_title_rot=0,border = T,col = RColorBrewer::brewer.pal(n = 11,name = "RdBu")[11:1],use_raster = F)
ha.pt=HeatmapAnnotation(df = data.frame(Tissue=ifelse(clinical.all$Tissue=="T","WT","NAT")[match(colnames(WT.protein.expr),clinical.all$SampleID)][order(clinical.all$Tissue[match(colnames(WT.protein.expr),clinical.all$SampleID)])]),col = list(Tissue=c(WT=pal_d3()(2)[2],NAT=pal_d3()(2)[1])),border = T)
ods.pt=order(clinical.all$Tissue[match(colnames(WT.protein.expr),clinical.all$SampleID)],
             apply(t(apply(WT.protein.expr[pathway.genes$Gene[1:61],],1,function(x){scales::rescale(x,to=c(-3,3))})),2,mean),
             apply(t(apply(WT.protein.expr[pathway.genes$Gene[-c(1:61)],],1,function(x){scales::rescale(x,to=c(-3,3))})),2,mean))
glist=lapply(unique(pathway.genes$Path),function(x){
  data.frame(genes=pathway.genes$Gene[pathway.genes$Path==x],fontsize=8)
})
names(glist)=unique(pathway.genes$Path)
row.ha =  rowAnnotation(textbox = anno_textbox(align_to=factor(pathway.genes$Path,levels=unique(pathway.genes$Path)), glist,background_gp = gpar(fill=NA)))
p.pt.path.genes=Heatmap(t(apply(log2(1+WT.protein.expr[pathway.genes$Gene,ods.pt]),1,function(x){scales::rescale(x,to=c(-3,3))})),show_column_names = F,row_names_gp = gpar(cex=0.35),cluster_rows = F,top_annotation = ha.pt,cluster_columns = F,name="Protein abundance",border = T,col=RColorBrewer::brewer.pal(n = 11,name = "RdBu")[11:1],use_raster = F,right_annotation = row.ha,show_row_names = F,row_split = factor(pathway.genes$Path,levels=unique(pathway.genes$Path)))
p.RNA.path.genes+p.pt.path.genes
### Figure 2E end

### Figure 2F
stat.ph=ph.result$t
names(stat.ph)=ph.result$ID
stat.ph.exc.up=stat.ph[!sapply(strsplit(names(stat.ph),"\\;"),function(x){x[1]}) %in% protein.results$ID[which(protein.results$logFC>log2(1.5) & protein.results$adj.P.Val<0.25 & protein.results$P.Value<0.05) ]]
kinase.substrate.mat=readRDS("../Phosphoprotein/kinase.substrate.mat.RDS")
obj.up=GSEA(sort(stat.ph.exc.up,decreasing = T),TERM2GENE = kinase.substrate.mat,minGSSize = 3,pAdjustMethod = "none",maxGSSize = 2000,pvalueCutoff = 1,exponent=1)
obj.up@result=obj.up@result[obj.up@result$pvalue<0.05 & obj.up@result$qvalue<0.25 & obj.up@result$NES>0,]
stat.ph.exc.down=stat.ph[!sapply(strsplit(names(stat.ph),"\\;"),function(x){x[1]}) %in% protein.results$ID[which(protein.results$logFC< -log2(1.5) & protein.results$adj.P.Val<0.25 & protein.results$P.Value<0.05) ]]
obj.down=GSEA(sort(stat.ph.exc.down,decreasing = T),TERM2GENE = kinase.substrate.mat,minGSSize = 3,pAdjustMethod = "none",maxGSSize = 1000)
obj.down@result=obj.down@result[obj.down@result$pvalue<0.05 & obj.down@result$qvalue<0.25 & obj.down@result$NES<0,]
up.ph.sites=ph.result$ID[ph.result$adj.P.Val<0.05 & ph.result$logFC>1]
down.ph.sites=ph.result$ID[ph.result$adj.P.Val<0.05 & ph.result$logFC< -1]
up.kinase.substrate=sapply(strsplit(obj.up@result$core_enrichment,"\\/")[obj.up@result$pvalue<0.05 & obj.up@result$qvalue<0.25 & obj.up@result$NES>0],intersect,up.ph.sites)
names(up.kinase.substrate)=obj.up@result$ID[obj.up@result$pvalue<0.05 & obj.up@result$qvalue<0.25 & obj.up@result$NES>0]

down.kinase.substrate=sapply(strsplit(obj.down@result$core_enrichment,"\\/")[obj.down@result$pvalue<0.05 & obj.down@result$qvalue<0.25 & obj.down@result$NES<0],intersect,down.ph.sites)
names(down.kinase.substrate)=obj.down@result$ID[obj.down@result$pvalue<0.05 & obj.down@result$qvalue<0.25 & obj.down@result$NES<0]
library(scales)
tmp.expr=t(apply(phosphoL6[unlist(c(up.kinase.substrate,down.kinase.substrate)),],1,function(x){
  rescale(x,to=c(-3,3))
}))
rownames(tmp.expr)=paste0("s",1:nrow(tmp.expr))
library(ComplexHeatmap)
library(ggsci)

de.kinase.substrate.list=sapply(unique(c(sort(names(up.kinase.substrate)),sort(names(down.kinase.substrate)))),function(x){
  y1=up.kinase.substrate[[x]]
  y2=down.kinase.substrate[[x]]
  y=c(y1,y2)
  y=y[!is.na(y)]
})
kinase.mat=matrix(0,nrow=nrow(tmp.expr),ncol=length(de.kinase.substrate.list))
rownames(kinase.mat)=rownames(tmp.expr)
colnames(kinase.mat)=names(de.kinase.substrate.list)
for(i in 1:ncol(kinase.mat)){
  kinase.mat[(cumsum(c(0,sapply(de.kinase.substrate.list,length)))[i]+1):cumsum(sapply(de.kinase.substrate.list,length))[i],i]=1
}
heat1=Heatmap(kinase.mat,cluster_columns = F,cluster_rows = F,col = c("white",RColorBrewer::brewer.pal(11,name = "RdBu")[4]),show_row_names = F,border = "lightgrey",column_names_side="top",name="Kinase-substrate pair",rect_gp=gpar(col="lightgrey",lwd=0.5))
ha1=HeatmapAnnotation(Tissue=rep(c("NAT","WT"),each=23),col = list(Tissue=c(NAT=pal_d3()(2)[1],WT=pal_d3()(2)[2])),show_legend = F)
sites=grep("MCM|HDAC|RB1|MLH1|MYC|EZH2|VIM|RPS6",unlist(de.kinase.substrate.list),value=T)
#ha2=HeatmapAnnotation()
ha2 = rowAnnotation(foo = anno_mark(at = grep("MCM|HDAC|RB1|MLH1|MYC|EZH2|VIM|RPS6",unlist(de.kinase.substrate.list)), 
                                    labels = grep("MCM|HDAC|RB1|MLH1|MYC|EZH2|VIM|RPS6",unlist(de.kinase.substrate.list),value=T)))
heat1+Heatmap(tmp.expr[,order(clinical.all$Tissue[match(colnames(tmp.expr),clinical.all$SampleID)])],cluster_columns = F,cluster_rows = F,show_column_names = F,row_names_gp = gpar(fontsize = 6),row_title_rot = c(0),row_title_gp=gpar(fontsize=10),top_annotation = ha1,column_split = rep(c("NAT","WT"),each=23),name = "Phosphorylated levels",col = RColorBrewer::brewer.pal(n = 11,name = "RdBu")[11:1],right_annotation = ha2,show_row_names = F)
## Figure 2F end

