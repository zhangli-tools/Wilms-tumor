laml=readRDS("../Mutations/laml.RDS")
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
### Figure 3A
oncoplot(laml,genes =  c("CTNNB1","WT1","AMER1","BRAF","TP53","CCNE1","PIK3CB"),removeNonMutated = F,colors = vc_cols)
### Figure 3A end

plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)
library("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
laml.tnm = trinucleotideMatrix(maf = laml, prefix = '', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
library('NMF')
laml.sign = estimateSignatures(mat = laml.tnm, nTry = 6)
plotCophenetic(res = laml.sign)
laml.sig = extractSignatures(mat = laml.tnm, n = 2)
laml.og30.cosm = compareSignatures(nmfRes = laml.sig, sig_db = "legacy")
### Figure 3B
maftools::plotSignatures(nmfRes = laml.sig, title_size = 1.2, sig_db = "legacy")
### Figure 3B end

OncogenicPathways(laml)
PlotOncogenicPathways(laml,pathways = "WNT")

mut.mat=t(sapply(unique(laml@data$Hugo_Symbol),function(x){
  y=rep(0,length(unique(laml@data$Tumor_Sample_Barcode)))
  names(y)=c(as.matrix(unique(laml@data$Tumor_Sample_Barcode)))
  y[c(as.matrix(laml@data$Tumor_Sample_Barcode[laml@data$Hugo_Symbol==x]))]=1
  y
}))
wnt.genes=c( "CTNNB1","AMER1")
db=msigdbr::msigdbr()
db.pathway=rbind(cbind(db$gs_description,db$gene_symbol)[grep("CP:",db$gs_subcat),],
                 cbind(db$gs_name,db$gene_symbol)[grep("HALLMARK",db$gs_name),])

WT.RNA.expr=readRDS("../RNA/WT.RNA.expr.RDS")
WT.protein.expr=readRDS("../Protein/WT.protein.expr.RDS")
clinical.all=readRDS("../Clinical/clinical.all.RDS")
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
mut.pathway.deg=function(gs){
  tumor.samples=intersect(colnames(mut.mat),colnames(WT.RNA.expr))
  normal.samples=colnames(WT.RNA.expr)[clinical.all$Tissue[match(colnames(WT.RNA.expr),clinical.all$SampleID)]=="N"]
  expr1=log2(WT.RNA.expr[,tumor.samples]+0.1)
  mut.label=apply(mut.mat[rownames(mut.mat) %in% gs,tumor.samples],2,sum)
  print(table(mut.label>0))
  deg1=limma.deg(expr1,ifelse(mut.label>0,1,0))
  expr2=log2(WT.RNA.expr[,c(normal.samples,tumor.samples[mut.label>0])]+0.1)
  deg2=limma.deg(expr2,rep(c(0,1),c(length(normal.samples),sum(mut.label>0))))
  y=list(deg1,deg2)
}
mut.pathway.dep=function(gs){
  tumor.samples=intersect(colnames(mut.mat),colnames(WT.protein.expr))
  normal.samples=colnames(WT.protein.expr)[clinical.all$Tissue[match(colnames(WT.protein.expr),clinical.all$SampleID)]=="N"]
  expr1=WT.protein.expr[,tumor.samples]
  mut.label=apply(mut.mat[rownames(mut.mat) %in% gs,tumor.samples],2,sum)
  print(table(mut.label>0))
  deg1=limma.deg(expr1,ifelse(mut.label>0,1,0))
  expr2=WT.protein.expr[,c(normal.samples,tumor.samples[mut.label>0])]
  deg2=limma.deg(expr2,rep(c(0,1),c(length(normal.samples),sum(mut.label>0))))
  y=list(deg1,deg2)
}

wnt.deg=mut.pathway.deg(wnt.genes)
wnt.dep=mut.pathway.dep(wnt.genes)
wnt.RNA.genes=intersect(rownames(wnt.deg[[1]])[wnt.deg[[1]]$P.Value<0.05 & wnt.deg[[1]]$logFC>log2(1)],
                        rownames(wnt.deg[[1]])[wnt.deg[[2]]$P.Value<0.05 & wnt.deg[[2]]$logFC>log2(1)])
obj.wnt.RNA=enricher(wnt.RNA.genes,TERM2GENE=db.pathway)
wnt.pt.genes=intersect(wnt.dep[[1]]$ID[wnt.dep[[1]]$P.Value<0.1 & wnt.dep[[1]]$logFC>0],
                        wnt.dep[[1]]$ID[wnt.dep[[2]]$P.Value<0.05 & wnt.dep[[2]]$logFC>0])
obj.wnt.pt=enricher(wnt.pt.genes,TERM2GENE=db.pathway)

wnt.cls.RNA=ifelse(clinical.all$Tissue[match(colnames(WT.RNA.expr),clinical.all$SampleID)]=="T","WT","NAT")
names(wnt.cls.RNA)=colnames(WT.RNA.expr)
wnt.cls.RNA[names(which(apply(mut.mat[wnt.genes,],2,sum)>0))]="Mutation"
wnt.cls.RNA[!(names(wnt.cls.RNA) %in% colnames(mut.mat)) & wnt.cls.RNA=="WT"]=NA
wnt.cls.RNA=wnt.cls.RNA[!is.na(wnt.cls.RNA)]
Fig.3C.data.RNA=t(apply(log2(0.1+WT.RNA.expr[c("MYC","WNT5A","DKK2","NKD1","FZD1","SERPINF1","DKK1","WNT10A","PORCN","CAMK2A","NKD2","NFATC2","WNT16","PLCB3","PPP3CC","PLCB4","NOTUM","WIF1"),names(wnt.cls.RNA)]),1,function(x){
  x=scales::rescale(x,to=c(-3,3))
  tapply(x,wnt.cls.RNA,mean)
}))
library(ComplexHeatmap)
ComplexHeatmap::pheatmap(Fig.3C.data.RNA[,c(2,3,1)],cluster_cols = F,cluster_rows = F,color = gplots::colorpanel(100,low = "blue",mid = "white",high = "red"),name="RNA expression",border_color = NA,border_gp=gpar(col="black"))
wnt.cls.pt=ifelse(clinical.all$Tissue[match(colnames(WT.protein.expr),clinical.all$SampleID)]=="T","WT","NAT")
names(wnt.cls.pt)=colnames(WT.protein.expr)
wnt.cls.pt[names(which(apply(mut.mat[wnt.genes,],2,sum)>0))]="Mutation"
wnt.cls.pt[!(names(wnt.cls.pt) %in% colnames(mut.mat)) & wnt.cls.pt=="WT"]=NA
wnt.cls.pt=wnt.cls.pt[!is.na(wnt.cls.pt)]

Fig.3C.data.pt=t(apply(WT.protein.expr[c("ACTN2","MYBPC1","MYH3","MYH8","MYL1","MYL2","MYL4","MYOM1","TNNC2","TNNT3","TTN","VIM"),names(wnt.cls.pt)],1,function(x){
  x=scales::rescale(x,to=c(-3,3))
  tapply(x,wnt.cls.pt,mean)
}))

Fig.3C.r2=ComplexHeatmap::pheatmap(Fig.3C.data.pt[,c(2,3,1)],cluster_cols = F,cluster_rows = F,color = RColorBrewer::brewer.pal(11,name = "RdBu")[11:1],name="Protein expression",border_color = NA,border_gp=gpar(col="black"),)
Fig.3C.data.RNA2=t(apply(log2(0.1+WT.RNA.expr[c("ACTN2","MYBPC1","MYH3","MYH8","MYL1","MYL2","MYL4","MYOM1","TNNC2","TNNT3","TTN","VIM"),names(wnt.cls.RNA)]),1,function(x){
  x=scales::rescale(x,to=c(-3,3))
  tapply(x,wnt.cls.RNA,mean)
}))
Fig.3C.r1=ComplexHeatmap::pheatmap(Fig.3C.data.RNA2[,c(2,3,1)],cluster_cols = F,cluster_rows = F,color =gplots::colorpanel(n = 11,low = "blue",mid = "white",high = "red"),name="RNA expression",border_color = NA,border_gp=gpar(col="black"),)
Fig.3C.r1+Fig.3C.r2
### Fig. 3C end



### Figure 3D
library(estimate)
outputGCT(as.data.frame(log2(1+WT.RNA.expr[,clinical.all$Tissue[match(colnames(WT.RNA.expr),clinical.all$SampleID)]=="T"])),output.f = "WT.tumor.gct")
estimateScore("WT.tumor.gct",output.ds = "estimate.out.gct",platform = "illumina")
library(data.table)
estimate.score=fread("estimate.out.gct",skip = 2)
estimate.score=as.matrix(estimate.score[,-c(1,2)],rownames=estimate.score$Description)

df.fig3d=data.frame(group=wnt.cls.RNA[intersect(names(wnt.cls.RNA),colnames(estimate.score))],
                    score=estimate.score[1,intersect(names(wnt.cls.RNA),colnames(estimate.score))])
library(ggpubr)
df.fig3d$group[df.fig3d$group=="WT"]="Wild"
ggboxplot(df.fig3d,x = "group",y = "score",color = "group",palette = "d3",add = "jitter",add.params = list(color="lightgrey",size=0.5),
          ylab="Stromal score",xlab="Genotype of WNT genes\n(CTNNB1/AMER1)")+
  stat_compare_means(comparisons = list(c("Wild","Mutation")))+
  theme(legend.position = "none",axis.line = element_blank(),panel.border = element_rect(fill = NA,colour = "black"))
### Figure 3D end

### Figure 3E
laml.gistic = readGistic(gisticAllLesionsFile = "../CNV/632363/WT.36.all_lesions.conf_90.txt",
                         gisticAmpGenesFile = "../CNV/632363/WT.36.amp_genes.conf_90.txt",
                         gisticDelGenesFile = "../CNV/632363/WT.36.del_genes.conf_90.txt",
                         gisticScoresFile = "../CNV/632363/WT.36.scores.gistic",cnLevel = "all")
gisticChromPlot(gistic = laml.gistic, markBands = "all",ref.build = "hg19")
### Figure 3E end


### Figure 3G
library(data.table)
cnv.mat=fread("../CNV/632363/WT.36.all_data_by_genes.txt")
cnv.mat=as.matrix(cnv.mat[,-c(1:3)],rownames = cnv.mat$`Gene Symbol`)
#colnames(cnv.mat)=clinical.all$SampleID[clinical.all$Tissue=="T"][match(colnames(cnv.mat),clinical.all$Patint_ID[clinical.all$Tissue=="T"])]
source("../CNV/multiOmicsViz.R")
sig.cnv.genes=names(which(apply(cnv.mat>0.4,1,sum)>=4 | apply(cnv.mat< -0.4,1,sum)>=4))
overlap.genes=intersect(intersect(sig.cnv.genes,rownames(WT.protein.expr)),rownames(WT.RNA.expr))
datacache=readRDS("../CNV/datacache.Rdata")
multiOmicsViz1(sourceOmics = as.data.frame(cnv.mat[overlap.genes,]),
               sourceOmicsName = "CNA",
               chrome_sourceOmics = "All",
               targetOmicsList = list(as.data.frame(log2(WT.RNA.expr[overlap.genes,]+0.1)),as.data.frame(WT.protein.expr[overlap.genes,])),
               targetOmicsName = c("mRNA","Protein"),
               chrome_targetOmics = "All",
               fdrThr = 0.25,
               outputfile = "Figure 3F",
               nThreads = 1)
calculateCorForTwoMatrices=function (matrix1, matrix2, fdr) 
{
  if (class(matrix1)[1] == "SummarizedExperiment") {
    matrix1 <- assays(matrix1, n = 1)
  }
  else {
    if (!is.matrix(matrix1) && !is.data.frame(matrix1)) {
      stop("Matrix1 should be a R matrix, data.frame \n      or SummarizedExperiment object.")
    }
  }
  if (class(matrix2)[1] == "SummarizedExperiment") {
    matrix2 <- assays(matrix2, n = 1)
  }
  else {
    if (!is.matrix(matrix2) && !is.data.frame(matrix2)) {
      stop("Matrix2 should be a R matrix, data.frame or \n        SummarizedExperiment object.")
    }
  }
  ov_sample <- intersect(colnames(matrix1), colnames(matrix2))
  if (length(ov_sample) < 6) {
    stop("The number of overlapping samples of two matrices is less than 6 \n    and thus can not perform correlation analysis. Please check the data.")
  }
  matrix1 <- matrix1[, ov_sample]
  matrix2 <- matrix2[, ov_sample]
  suppressWarnings(corrArray <- cor(t(matrix1), t(matrix2), 
                                    method = "spearman"))
  corrArray[is.na(corrArray)] <- 0
  n <- t(!is.na(t(matrix1))) %*% (!is.na(t(matrix2)))
  suppressWarnings(t <- (corrArray * sqrt(n - 2))/sqrt(1 - 
                                                         corrArray^2))
  suppressWarnings(corrP <- 2 * (1 - pt(abs(t), (n - 2))))
  corrP[is.na(corrP)] <- 0
  corrP[corrP > 1] <- 1
  corrSig <- apply(corrP, 1, function(x){
    ifelse(x<fdr,1,0)})
  corrSig <- t(corrSig)
  corrArray <- sign(corrArray) * corrSig
  return(corrArray)
}
library(data.table)

sig.RNA=calculateCorForTwoMatrices(matrix1 = as.data.frame(cnv.mat[overlap.genes,]),as.data.frame(log2(WT.RNA.expr[overlap.genes,]+1)),fdr = 0.25)
sig.protein=calculateCorForTwoMatrices(matrix1 = as.data.frame(cnv.mat[overlap.genes,]),as.data.frame(WT.protein.expr[overlap.genes,]),fdr = 0.25)

RNA.results=limma.deg(log2(0.1+WT.RNA.expr),ifelse(clinical.all$Tissue[match(colnames(WT.RNA.expr),clinical.all$SampleID)]=="T",1,0))
WT.vs.NAT.deglist=list(rownames(RNA.results)[RNA.results$adj.P.Val<0.25  & RNA.results$logFC>
                                               log2(1.5)],
                       rownames(RNA.results)[RNA.results$adj.P.Val<0.05  & RNA.results$logFC<
                                               -log2(1.5)])

protein.results=limma.deg(WT.protein.expr,ifelse(clinical.all$Tissue[match(colnames(WT.protein.expr),clinical.all$SampleID)]=="T",1,0))
WT.vs.NAT.deplist=list(protein.results$ID[protein.results$adj.P.Val<0.25 & protein.results$P.Value<0.05 & protein.results$logFC>log2(1.25)],
                       protein.results$ID[protein.results$adj.P.Val<0.25 & protein.results$P.Value<0.05 & protein.results$logFC< -log2(1.25)])
shared.cis=names(which(diag(sig.RNA[overlap.genes,])>0 & diag(sig.protein[overlap.genes,])>0))
shared.cis.amp=names(which(apply(cnv.mat[shared.cis,]>0.4,1,sum)>=4 & apply(cnv.mat[shared.cis,]>0.4,1,sum)>apply(cnv.mat[shared.cis,]< -0.4,1,sum)))
shared.cis.del=names(which(apply(cnv.mat[shared.cis,]< -0.4,1,sum)>=4 & apply(cnv.mat[shared.cis,]>0.4,1,sum)<apply(cnv.mat[shared.cis,]< -0.4,1,sum)))
shared.cis.amp1=intersect(shared.cis.amp,intersect(WT.vs.NAT.deglist[[1]],WT.vs.NAT.deplist[[1]]))
shared.cis.del1=intersect(shared.cis.del,intersect(WT.vs.NAT.deglist[[2]],WT.vs.NAT.deplist[[2]]))
db.C1=cbind(db$gs_name,db$gene_symbol)[db$gs_cat=="C1",]
db.C1.1=db.C1
db.C1.1[,1]=gsub("[0-9]+$","",db.C1.1[,1])
obj.chr=compareCluster(gene~group,fun = "enricher",data = data.frame(gene=c(shared.cis.amp1,shared.cis.del1),group=rep(c("Gain","Deletion"),c(length(shared.cis.amp1),length(shared.cis.del1)))),TERM2GENE=db.C1.1,minGSSize=4,maxGSSize=2500,pAdjustMethod="fdr")
obj.chr@compareClusterResult=obj.chr@compareClusterResult[obj.chr@compareClusterResult$Count>=5 & obj.chr@compareClusterResult$qvalue<0.05,]
dotplot(obj.chr,showCategory=10,size="count")+
  theme(panel.grid.major = element_line(linetype = 3,colour = "lightgrey",size = 0.5),axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  scale_color_material(palette = "red",reverse = T,alpha = 0.65)+
  ggtitle("cis-regulatory CNAs")
### Figure 3G end


### Figure 3H
go.term.del=c("electron transport chain","carboxylic acid catabolic process",
              "fatty acid oxidation","hexose metabolic process",
              "organic acid metabolic process")
go.term.amp=c("RNA splicing","regulation of G0 to G1 transition",
              "DNA replication","mitotic sister chromatid segregation",
              "chromatin remodeling","G1/S transition of mitotic cell cycle")
df.cis.genes=data.frame(gene=c(shared.cis.amp1,shared.cis.del1),group=rep(c("Gain","Deletion"),c(length(shared.cis.amp1),length(shared.cis.del1))))
library(org.Hs.eg.db)
obj.path=compareCluster(gene~group,fun = "enrichGO",data = df.cis.genes[df.cis.genes$gene %in% unlist(strsplit(obj.chr@compareClusterResult$geneID[1:13],"\\/")),],keyType="SYMBOL",OrgDb=org.Hs.eg.db,ont="BP",minGSSize=4,maxGSSize=2500,pAdjustMethod="none")
obj.path@compareClusterResult=obj.path@compareClusterResult[obj.path@compareClusterResult$Count>=2,]
dotplot(obj.path,showCategory=c(go.term.amp,go.term.del),label_format=120)
df.cnv.go=data.frame(Cluster=rep(obj.path@compareClusterResult$Cluster,obj.path@compareClusterResult$Count),
                     Description=rep(obj.path@compareClusterResult$Description,obj.path@compareClusterResult$Count),
                     gene=unlist(strsplit(obj.path@compareClusterResult$geneID,"\\/")))
df.cnv.go=df.cnv.go[df.cnv.go$Description %in% c(go.term.amp,go.term.del),]
df.cnv.go$chr=db.C1.1[match(df.cnv.go$gene,db.C1.1[,2]),1]
cnv.status=apply(rbind(cnv.mat[unique(df.cnv.go$gene[df.cnv.go$Cluster=="Deletion"]),]< -0.4,cnv.mat[unique(df.cnv.go$gene[df.cnv.go$Cluster=="Gain"]),]>0.4),1,function(x){ifelse(x,1,0)})
CNV.cor.mat=1-as.matrix(dist(t(cnv.status),method = "binary"))
CNV.gene.pairs=which(CNV.cor.mat>0.35,arr.ind = T)
CNV.gene.pairs=CNV.gene.pairs[CNV.gene.pairs[,1]<CNV.gene.pairs[,2],]
CNV.gene.pairs=t(apply(CNV.gene.pairs,1,function(x){y=c(rownames(CNV.cor.mat)[x],CNV.cor.mat[x[1],x[2]])}))
CNV.gene.pairs=as.data.frame(CNV.gene.pairs)
CNV.gene.pairs$chr1=df.cnv.go$chr[match(CNV.gene.pairs$V1,df.cnv.go$gene)]
CNV.gene.pairs$chr2=df.cnv.go$chr[match(CNV.gene.pairs$V2,df.cnv.go$gene)]
CNV.gene.pairs.test=apply(as.matrix(CNV.gene.pairs),1,function(x){
  glist=strsplit(obj.path@compareClusterResult[obj.path@compareClusterResult$Description %in% c(go.term.del,go.term.amp),]$geneID,"\\/")
  test=sapply(glist,function(t){sum(t %in% x)})
  sum(test==2)>0
})
library(igraph)
g=graph_from_data_frame(data.frame(CNV.gene.pairs[,1],CNV.gene.pairs[,2],dist=as.numeric(CNV.gene.pairs[,3]))[CNV.gene.pairs.test,],directed = F)
deg <- degree(g, mode="all")
cols.pt=c(RColorBrewer::brewer.pal(n = 14,name = "Paired"),RColorBrewer::brewer.pal(n = 8,name = "Set2"))
gene.cls=paste(df.cnv.go$Cluster[match(V(g)$name,df.cnv.go$gene)],df.cnv.go$chr[match(V(g)$name,df.cnv.go$gene)],sep=":")
V(g)$color=cols.pt[match(gene.cls,unique(gene.cls))]
V(g)$size <- log2(deg)*2
V(g)$label.cex=0.75
E(g)$width <- as.numeric(apply(as_edgelist(g),1,function(x){CNV.gene.pairs[CNV.gene.pairs$V1==x[1] & CNV.gene.pairs$V2==x[2],"V3"]}))
E(g)$color="grey"
V(g)$size=4
plot(g, layout=layout_with_kk,edge.arrow.size=.25,edge.curved=.05,edge.label.cex=0.5,vertex.label.cex=0.5,vertex.frame.color="grey",vertex.label.color="black",
     vertex.label.cex=0.5)                                                                                                           
### Figure 3H end


                                                                                                              
                                                                                                              