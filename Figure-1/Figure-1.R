## Load datasets
clinical.all=readRDS("../Clinical/clinical.all.RDS")
WT.RNA.expr=readRDS("../RNA/WT.RNA.expr.RDS")
WT.protein.expr=readRDS("../Protein/WT.protein.expr.RDS")
laml=readRDS("../Mutations/laml.RDS")
phosphoL6=readRDS("../Phosphoprotein/phosphoL6.RDS")

## Figure 1B
mat=matrix(0,nrow=9,ncol=length(unique(clinical.all$Patient_ID)),dimnames = list(NULL,unique(clinical.all$Patient_ID)))
mat[1,]=1
mat[2,clinical.all$Patient_ID[match(laml@clinical.data$Tumor_Sample_Barcode,clinical.all$SampleID)]]=2
mat[3,clinical.all$Patient_ID[match(laml@clinical.data$Tumor_Sample_Barcode,clinical.all$SampleID)]]=3
mat[4,clinical.all$Patient_ID[match(colnames(WT.RNA.expr)[clinical.all$Tissue[match(colnames(WT.RNA.expr),clinical.all$SampleID)]=="T"],clinical.all$SampleID)]]=4
mat[5,clinical.all$Patient_ID[match(colnames(WT.RNA.expr)[clinical.all$Tissue[match(colnames(WT.RNA.expr),clinical.all$SampleID)]=="N"],clinical.all$SampleID)]]=5
mat[6,clinical.all$Patient_ID[match(colnames(WT.protein.expr)[clinical.all$Tissue[match(colnames(WT.protein.expr),clinical.all$SampleID)]=="T"],clinical.all$SampleID)]]=6
mat[7,clinical.all$Patient_ID[match(colnames(WT.protein.expr)[clinical.all$Tissue[match(colnames(WT.protein.expr),clinical.all$SampleID)]=="N"],clinical.all$SampleID)]]=7
mat[8,clinical.all$Patient_ID[match(colnames(phosphoL6)[clinical.all$Tissue[match(colnames(phosphoL6),clinical.all$SampleID)]=="T"],clinical.all$SampleID)]]=8
mat[9,clinical.all$Patient_ID[match(colnames(phosphoL6)[clinical.all$Tissue[match(colnames(phosphoL6),clinical.all$SampleID)]=="N"],clinical.all$SampleID)]]=9

mat1=mat[,order(mat[1,],mat[2,],mat[3,],mat[4,],mat[5,],mat[6,],mat[7,],mat[8,],mat[9,],decreasing = T)]
df.mat=data.frame(omics=rep(c("Total","Genome-T","Genome-N","Transcriptome-T","Transcriptome-N",
                              "Proteome-T","Proteome-N","Phosphoproteome-T","Phosphoproteome-N"),ncol(mat)),
                  sample=rep(colnames(mat),each=9),value=c(mat))
library(ggplot2);library(ggsci)
df.mat$value1=as.character(df.mat$value)
angles=seq(0,360,length.out=ncol(mat))
df.mat$sample.name=df.mat$sample
#df.mat$sample.name=paste("WT",sapply(match(df.mat$sample,colnames(mat1)),function(x){paste(c(rep(0,3-nchar(x)),x),collapse = "")}),sep="")
#clinical.all$publish.id=df.mat$sample.name[match(clinical.all$Patint_ID,df.mat$sample)]
ggplot(df.mat, aes(x=sample.name, y=factor(omics,levels = unique(omics)), fill=value1)) +
  geom_tile(colour="white") +
  scale_fill_manual(labels=c("NA","Total","Exome (WT)","Exome (NAT)",
                             "Transcriptome (WT)","Transcriptome (NAT)",
                             "Proteome (WT)","Proteome (NAT)",
                             "Phosphoproteome (WT)","Phosphoproteome (NAT)"),
                    values = c("grey95",pal_d3(palette = "category20",alpha = 0.5)(10)[-8] )) +
  #scale_y_discrete(breaks=y_breaks, labels=y_labels) +
  coord_polar(theta="x") +
  theme(panel.background=element_blank(),
        axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y  = element_blank(),
        axis.text.x = element_text(vjust = 1,angle = 270-(angles[1:95]+angles[2:96])/2)
        )+
  scale_x_discrete()+guides(fill=guide_legend(title = "Omics"))
#### Figure 1B end

#### Figure 1C
vars=c("Patient_ID","Gender","Pre-treated_tumors_histology","Risk_group","Site","Pre-treated_(Chemotherapy)","Matastasis","Stage","fulstate")
patient.clinical=unique(as.matrix(clinical.all[clinical.all$Tissue=="T",..vars]))
patient.clinical[is.na(patient.clinical)]="NA"
patient.clinical[patient.clinical[,5]=="bilateral",5]="Bilateral"
patient.clinical[patient.clinical[,5]=="",5]="NA"
patient.clinical[,7]=ifelse(patient.clinical[,7]=="","No","Yes")
df.pie=lapply(colnames(patient.clinical)[-1],function(x){
  y1=table(patient.clinical[,x])/nrow(patient.clinical)
  y=data.frame(var=x,group=names(y1),prop=c(y1))
})
df.pie=eval(as.call(c(rbind,df.pie)))
library(scales)
df.pie$group=paste(df.pie$group," (",percent(df.pie$prop),")",sep="")
library(ggrepel)
#  geom_text_repel(aes(y =c(cumsum(c(0,prop[match(levels(factor(group)),group)]))[-length(prop)]+prop[match(levels(factor(group)),group)])/2 , 
#label = percent(prop[match(levels(factor(group)),group)])), size=5)

plist=lapply(colnames(patient.clinical)[-1],function(t){
  p=ggplot(df.pie[df.pie$var==t,], aes(x=var, y=prop,fill=group,group=group))+
    geom_bar(width = 1, stat = "identity")+coord_polar("y",direction = -1)+theme_void()
  #p=p+scale_fill_discrete(labels=paste(p$data$group," (",percent(p$data$prop),")",sep=""))  
  if(t=="Pre-treated_tumors_histology")
  {p=p+scale_fill_lancet(alpha = 0.6)}else if(t=="Gender")
  {
    p=p+scale_fill_nejm(alpha = 0.6)
  }else if(t=="Risk_group")
  {
    p=p+scale_fill_aaas(alpha=0.6)
  }else if(t=="Site"){
    p=p+scale_fill_jco(alpha=0.6)
  }else if(t=="fulstate"){
    p=p+scale_fill_jco(alpha=0.6)
  }else{
    p=p+scale_fill_jama(alpha=0.6)
  }
  
  p=p+guides(fill=guide_legend(title=t))
})
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
library(ggplot2)
library(gridExtra)
mylegend=lapply(plist,g_legend)
library(ggpubr)
plist1=lapply(colnames(patient.clinical)[-1],function(t){
  p=ggplot(df.pie[df.pie$var==t,], aes(x=var, y=prop,fill=group,group=group))+
    geom_bar(width = 1, stat = "identity")+coord_polar("y",direction = -1)+theme_void()
  if(t=="Pre-treated_tumors_histology")
  {p=p+scale_fill_lancet(alpha = 0.6)}else if(t=="Gender")
  {
    p=p+scale_fill_nejm(alpha = 0.6)
  }else if(t=="Risk_group")
  {
    p=p+scale_fill_aaas(alpha=0.6)
  }else if(t=="Site"){
    p=p+scale_fill_jco(alpha=0.6)
  }else if(t=="fulstate"){
    p=p+scale_fill_jco(alpha=0.6)
  }else {
    p=p+scale_fill_jama(alpha=0.6)
  }
  
  p+theme(text = element_blank(),legend.position = "none")
})
library(ggplot2);library(ggsci);library(ggpubr)
ggarrange(plotlist = list(ggarrange(plotlist =plist1,nrow=1),ggarrange(plotlist =mylegend,nrow=1,align="hv")),nrow=2,heights = c(0.7,0.3)[2:1],align = "h")
#### Figure 1C end

#### Figure 1D
raw.phosph.data=readRDS("../Phosphoprotein/raw.phosph.data.RDS")
raw.phosph.data.sites=readRDS("../Phosphoprotein/raw.phosph.data.sites.RDS")
WT.protein.expr.raw=readRDS("../Protein/WT.protein.expr.raw.RDS")
RNA.genes=apply(WT.RNA.expr>1,2,sum)
pt.genes=apply(WT.protein.expr.raw[,colnames(WT.protein.expr)]>0,2,sum)
ph.sites=apply(raw.phosph.data,2,function(x){sum(!is.na(x))})
ph.genes=apply(raw.phosph.data,2,function(x){g=raw.phosph.data.sites[!is.na(x)];
g=length(unique(sapply(strsplit(g,";"),function(t){t[1]})))
})

df.count=data.frame(
  count=c(RNA.genes,pt.genes,ph.sites,ph.genes),
  omics=factor(rep(c("RNAs","Proteins","Phosphorylation sites","Phosphoproteins"),c(length(RNA.genes),length(pt.genes),length(ph.sites),length(ph.genes))),levels = c("RNAs","Proteins","Phosphorylation sites","Phosphoproteins")),
  Tissue=ifelse(clinical.all$Tissue[match(names(c(RNA.genes,pt.genes,ph.sites,ph.genes)),clinical.all$SampleID)]=="T","WT","NAT"),
  sample=c(names(RNA.genes),names(pt.genes),names(ph.sites),names(ph.genes))
)
library(ggplot2);library(ggsci);library(ggpubr)
plist=lapply(levels(df.count$omics),function(x){
  p=ggplot(df.count[df.count$omics==x,], aes(y=count,fill=Tissue,x=Tissue)) +
    geom_boxplot(outlier.color = NA)+scale_fill_d3(alpha = 0.5)+
    theme(panel.background = element_rect(fill = NA,colour = "black"),
          panel.grid.major = element_line(colour = "grey",linetype = 3),
          axis.text.x =  element_blank(),axis.ticks.x = element_blank(),
          strip.background = element_blank(),strip.text =  element_blank(),
          axis.text.y  = element_text(size=12))+guides(color=NULL)+
    xlab("")+stat_compare_means(method = "wilcox.test",comparisons = list(c("WT","NAT")))+
    geom_jitter(data = df.count[df.count$omics==x,],aes(color=Tissue),shape=16, position=position_jitter(0.2))+scale_color_d3(alpha = 0.5)
  p=p+ylim(min(p$data$count)*0.95,max(p$data$count)*1.05)
  ylabs=c("# of RNAs","# of proteins","# of Phosphorylated sites","# of phosphorylated proteins")[match(unique(p$data$omics),levels(df.count$omics))]
  p=p+ylab(ylabs)
})
ggarrange(plotlist = plist,nrow=1,ncol=4,common.legend = T)
#### Figure 1D end

#### Figure 1E
overlap.samples=intersect(colnames(WT.RNA.expr),colnames(WT.protein.expr))
overlap.tumor.samples=overlap.samples[clinical.all$Tissue[match(overlap.samples,clinical.all$SampleID)]=="T"]
overlap.normal.samples=overlap.samples[clinical.all$Tissue[match(overlap.samples,clinical.all$SampleID)]=="N"]
k=0.1
RNA.pt.overlap.genes.n=intersect(names(which(apply(WT.protein.expr.raw[rownames(WT.protein.expr),overlap.normal.samples]==0,1,sum)<=(34*k))),
                                 names(which(apply(WT.RNA.expr[,overlap.normal.samples]<=1,1,sum)<=(34*k))))
RNA.pt.overlap.genes.t=intersect(names(which(apply(WT.protein.expr.raw[rownames(WT.protein.expr),overlap.tumor.samples]==0,1,sum)<=(59*k))),
                                 names(which(apply(WT.RNA.expr[,overlap.tumor.samples]<=1,1,sum)<=(59*k))))

cor.RNA.vs.pt.n=sapply(RNA.pt.overlap.genes.n,function(x){
  cor(WT.protein.expr[x,overlap.normal.samples],WT.RNA.expr[x,overlap.normal.samples],method="spearman")
})
cor.RNA.vs.pt.t=sapply(RNA.pt.overlap.genes.t,function(x){
  cor(WT.protein.expr[x,overlap.tumor.samples],WT.RNA.expr[x,overlap.tumor.samples],method="spearman")
})

library(clusterProfiler)
db=msigdbr::msigdbr()
db.pathway=rbind(cbind(db$gs_description,db$gene_symbol)[grep("CP:",db$gs_subcat),],
                 cbind(db$gs_name,db$gene_symbol)[grep("HALLMARK",db$gs_name),])
cor.RNA.vs.pt.t.pathway=GSEA(sort(cor.RNA.vs.pt.t,decreasing = T),TERM2GENE=db.pathway,pAdjustMethod = "none",pvalueCutoff = 1)
cor.RNA.vs.pt.n.pathway=GSEA(sort(cor.RNA.vs.pt.n,decreasing = T),TERM2GENE=db.pathway,pAdjustMethod = "none",pvalueCutoff = 1)

pathways=c("CDK Regulation of DNA Replication","Genes controlling nephrogenesis","Tumor suppressor activity of SMARCB1","Development of ureteric collection system",
           "Respiratory electron transport","Fatty acid metabolism","IL-2 signaling pathway","B cell receptor signaling pathway")

cor.pathway.plot1=function(corr,gsea.obj,col="red",pathways){
  df.cor=data.frame(cor=corr,
                    rank=rank(-corr),
                    gene=names(corr))
  p1=ggplot(df.cor, aes(x=cor)) +
    geom_histogram(aes(y=..density..),fill=col,color=col, position="identity")+
    geom_density(color=col)+
    geom_vline(data=data.frame(x=median(corr)), aes(xintercept=x),
               linetype="dashed")+
    labs(x="Spearman correlation", y = "Density")+
    theme_classic()

  p1=p1+geom_text(data = data.frame(x=median(corr),y=layer_scales(p1)$y$range$range[2]*1.1,label=paste("median =",round(median(corr),2))),aes(x=x,y=y,label=label))+
    coord_cartesian(xlim = c(-2,1))
 
  p2=ggplot(data = df.cor)+xlim(min(df.cor$cor),max(df.cor$cor))
  
  for(i in 1:length(pathways)){
    print(i)
    p2=p2+geom_segment(data =cbind(df.cor[ df.cor$gene %in% gsea.obj@geneSets[[pathways[i]]],],y=length(pathways)-i+1) ,aes(x = cor,y = y-1,xend = cor,yend = y-0.2),size=0.2)+
      theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.background = element_rect(fill = NA))
    
  }
  
  p2=p2+geom_text(data = data.frame(x=min(df.cor$cor),y=seq(0.4,length(pathways),by=1)[length(pathways):1],label=pathways),aes(x=x,y=y,label=label),hjust=1)+xlim(-2,1)+
    xlab("")+ylab("")+coord_cartesian(xlim = c(-2,1))
  plot_grid(plotlist = list(p1,p2),nrow=2,ncol=1,align = "v")
}
library(cowplot)
p.normal=cor.pathway.plot1(cor.RNA.vs.pt.n,gsea.obj = cor.RNA.vs.pt.n.pathway,col = pal_d3(alpha = 0.5)(10)[1],pathways = pathways)
p.tumor=cor.pathway.plot1(cor.RNA.vs.pt.t,gsea.obj = cor.RNA.vs.pt.t.pathway,col = pal_d3(alpha = 0.5)(10)[2],pathways = pathways)
ggarrange(plotlist = list(p.normal,p.tumor),nrow=1,ncol=2)
#### Figure 1E end


