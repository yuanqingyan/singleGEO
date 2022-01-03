## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = FALSE, 
  message = FALSE,
  comment = "#>",
  fig.width=6, 
  fig.height=6,
  fig.align = "center",
  echo=TRUE
)

## ---- eval=FALSE--------------------------------------------------------------
#  library(devtools)
#  devtools::install_github("yuanqingyan/singleGEO")

## ----eval=TRUE----------------------------------------------------------------
library(singleGEO)
##Download the GEO meta database. Remove "Demo=TRUE" or set "Demo=FALSE" when you start your own project.
MyDataBase<-GetGeoMetaDatabase(Sqlfile=NULL,Demo=TRUE)

## ----eval=TRUE----------------------------------------------------------------
# Will study mouse and human data
MyOrganism<-c("Mus musculus","Homo sapiens")
# Will study lung. Some study may use "pulmonary" instead of "lung", so two options are specified
search_kw1<-c("lung","pulmonary")
KeyList<-list(kw1=search_kw1)
MySearch<-Get_Keyword_Meta(GeoDataBase=MyDataBase,
                           Organism=MyOrganism, 
                           InputSearchList=KeyList)
colnames(MySearch)
unique(MySearch$gse)

## ----eval=FALSE---------------------------------------------------------------
#  head(MySearch[,c("title","gse")])

## ----eval=TRUE,echo=FALSE-----------------------------------------------------
head(MySearch[,c("title","gse")])

## ----eval=FALSE---------------------------------------------------------------
#  MySearch[MySearch$gse=="GSE158127","summary"][1]

## ----eval=TRUE,echo=F---------------------------------------------------------
MySearch[MySearch$gse=="GSE158127","summary"][1]

## ----eval=TRUE----------------------------------------------------------------
Organism_hs<-c("Homo sapiens")
Data_platform<-c("scRNAseq","scATACseq")
search_kw_ade<-c("adenocarcinoma")
KeyList_ade<-list(kw_ad=search_kw_ade)
MySearch_ade<-Get_Keyword_Meta(GeoDataBase=MyDataBase,
                               Organism=Organism_hs,
                               InputSearchList=KeyList_ade,
                               Platform=Data_platform)
unique(MySearch_ade$gse)

## ----eval=FALSE---------------------------------------------------------------
#  MySearch_ade[MySearch_ade$gse=="GSE142285","title"][1]

## ----eval=TRUE,echo=F---------------------------------------------------------
MySearch_ade[MySearch_ade$gse=="GSE142285","title"][1]

## ----eval=TRUE----------------------------------------------------------------
MyOrganism_ms<-c("Mus musculus")
search_brain<-c("brain")
KeyList_brain<-list(kw_brain=search_brain)

MySearch_brain_scRNA<-Get_Keyword_Meta(GeoDataBase=MyDataBase,
                                       Organism=MyOrganism_ms,
                                       InputSearchList=KeyList_brain,
                                       Platform="scRNAseq")
MySearch_brain_scATAC<-Get_Keyword_Meta(GeoDataBase=MyDataBase,
                                        Organism=MyOrganism_ms,
                                        InputSearchList=KeyList_brain,
                                        Platform="scATACseq")
(Both_RNA_ATAC<-intersect(unique(MySearch_brain_scRNA$gse),
                          unique(MySearch_brain_scATAC$gse)))

## ----eval=FALSE---------------------------------------------------------------
#  MySearch_brain_scRNA[MySearch_brain_scRNA$gse==Both_RNA_ATAC[1],"overall_design"][1]

## ----eval=TRUE,echo=F---------------------------------------------------------
MySearch_brain_scRNA[MySearch_brain_scRNA$gse==Both_RNA_ATAC[1],"overall_design"][1]

## ----eval=TRUE----------------------------------------------------------------
search_kw2<-c("fibrosis")
search_kw3<-c("new","novel")
search_kw4<-c("subtype","subpopulation")
KeyList_new<-list(kw1=search_kw1,
                  kw2=search_kw2,
                  kw3=search_kw3,
                  kw4=search_kw4)
MySearch_new<-Get_Keyword_Meta(GeoDataBase=MyDataBase,
                               Organism=MyOrganism_ms,
                               InputSearchList=KeyList_new)
unique(MySearch_new$gse)

## ----eval=FALSE---------------------------------------------------------------
#  MySearch_new[1,"title"]

## ----eval=TRUE,echo=F---------------------------------------------------------
MySearch_new[1,"title"]

## ----eval=FALSE---------------------------------------------------------------
#  MySearch_new$summary[1]

## ----eval=TRUE,echo=F---------------------------------------------------------
MySearch_new$summary[1]

## ----eval=TRUE,echo=F---------------------------------------------------------
MySearch_new$overall_design[1]

## ----eval=TRUE,echo=T---------------------------------------------------------
gsm_info<-Get_GSMFromGSE(GeoDataBase=MyDataBase,
                         GSE.ID=MySearch_new$gse[1])
colnames(gsm_info)

## ----eval=FALSE---------------------------------------------------------------
#  gsm_info[,c("title","gsm")]

## ----eval=TRUE,echo=F---------------------------------------------------------
gsm_info[,c("title","gsm")]

## ----eval=TRUE----------------------------------------------------------------
search_kw5<-c("fas-signaling","fas signaling","fas pathway","fas-pathway")
KeyList_fas<-list(kw1=search_kw1,
                  kw2=search_kw2,
                  kw5=search_kw5)
MySearch_fas<-Get_Keyword_Meta(GeoDataBase=MyDataBase,
                               Organism=MyOrganism_ms,
                               InputSearchList=KeyList_fas)
unique(MySearch_fas$gse)

## ----eval=FALSE---------------------------------------------------------------
#  MySearch_fas[1,"title"]

## ----eval=TRUE,echo=F---------------------------------------------------------
MySearch_fas[1,"title"]

## ----eval=TRUE,echo=F---------------------------------------------------------
MySearch_fas$overall_design[1]

## ----eval=FALSE---------------------------------------------------------------
#  MySearch_fas$summary[1]

## ----eval=TRUE,echo=F---------------------------------------------------------
MySearch_fas$summary[1]

## ----eval=FALSE---------------------------------------------------------------
#  MySearch_fas$supplementary_file[1]

## ----eval=TRUE,echo=F---------------------------------------------------------
MySearch_fas$supplementary_file[1]

## ----eval=TRUE----------------------------------------------------------------
myGeoID<-c("GSE140032")
##Set DecompressFile=TRUE if you want to download and decompress the data in the same time
DownloadFileInfo<-Get_Geo_Data(DownGSEID=myGeoID,DecompressFile=FALSE)
DownloadFileInfo
ExtractFile(ExtractID=myGeoID,DecompressType=c("tar","gz"))
list.files("GSE140032")

## ----eval=TRUE----------------------------------------------------------------
###The following code generate an example dataset from GSE134174 (The same data has been built and can be loaded by data(testData_GSE134174)). Note that only 4 samples with 800 cells were selected for illustration purpose.
# Get_Geo_Data(DownGSEID="GSE134174",DecompressFile="TRUE",DecompressType="All")
# selPatient<-c("T101","T85","T153","T164")
# raw_GSE134174<-read.delim("./GSE134174/GSE134174_Processed_invivo_raw.txt",header=T,sep="\t")
# meta_GSE134174<-read.delim("./GSE134174/GSE134174_Processed_invivo_metadata.txt",header=T,sep="\t")
# meta_sel0<-meta_GSE134174[meta_GSE134174$Donor %in% selPatient,]
# raw_sel0<-raw_GSE134174[,colnames(raw_GSE134174) %in% meta_sel0$Cell]
# set.seed(12345);sampleIndex<-sample(1:nrow(meta_sel0),size=800)
# meta_sel<-meta_sel0[sampleIndex,];raw_sel<-raw_sel0[,paste(meta_sel$Cell)]
# raw_sel<-raw_sel[rowSums(raw_sel)>0,]
# testData_GSE134174<-list(TwoRawData=raw_sel,TwoMetaData=meta_sel)

data(testData_GSE134174)
test_meta<-testData_GSE134174$TwoMetaData
test_data<-testData_GSE134174$TwoRawData
list_GSE134174<-Splitdata_MakeDataList(InputData=test_data,Group=test_meta$Donor)
##The following code setting a big number of maxFeature,maxCount,maxMT and small number ##of minFeature,minCount will keep all the cells
seu_GSE134174<-MakeSeuObj_FromRawRNAData(RawList=list_GSE134174,
                                      GSE.ID="GSE134174",
                                      MinFeature=1,
                                      MaxFeature=750000,
                                      MinCount=1,
                                      MaxCount=4000000,
                                      MaxMT=100)
names(seu_GSE134174)
Seurat::DimPlot(seu_GSE134174[[1]],repel = TRUE)
####Use sctransform normalization and make Seurat object
seu_GSE134174_sct<-MakeSeuObj_FromRawRNAData(RawList=list_GSE134174,
                                          GSE.ID="GSE134174",
                                          MinFeature=1,
                                          MaxFeature=750000,
                                          MinCount=1,
                                          MaxCount=4000000,
                                          MaxMT=100,
                                          Norm.method = "sct")
names(seu_GSE134174_sct)

## ----eval=TRUE----------------------------------------------------------------
#####to make Seurat object from normalized by LogNormalized method
#Norm_GSE134174<-lapply(seu_GSE134174,function(x) {Seurat::GetAssayData(x,slot = "data")})
data(Norm_GSE134174)
seu_Norm_GSE134174<-MakeSeuObj_FromNormalizedRNAData(NormList=Norm_GSE134174,GSE.ID="GSE134174")
Seurat::DimPlot(seu_Norm_GSE134174[[1]],repel = TRUE)

#####to make Seurat object from normalized by sctransform method
list_Norm.sct_GSE134174<-lapply(seu_GSE134174_sct,function(x) Seurat::GetAssayData(x,slot="data"))
seu_Norm.sct_GSE134174<-MakeSeuObj_FromNormalizedRNAData(NormList=list_Norm.sct_GSE134174,Frow.which.Norm="sct",GSE.ID="GSE134174")
Seurat::DimPlot(seu_Norm.sct_GSE134174[[1]],repel = TRUE)


## ----eval=TRUE----------------------------------------------------------------
Int_GSE134174<-SeuObj_integration(Object.list=seu_GSE134174)
###Add additional meta data
row.names(test_meta)<-test_meta$Cell
test_meta<-test_meta[row.names(Int_GSE134174@meta.data),]
Int_GSE134174<-Seurat::AddMetaData(object=Int_GSE134174,metadata =as.data.frame(test_meta))
Seurat::DimPlot(Int_GSE134174,group.by="cluster_ident",label=TRUE,repel = TRUE)+ NoLegend()

#####Data integration with one single dataset normalized with sctransform method
Int_GSE134174_sct<-SeuObj_integration(Object.list=seu_GSE134174_sct,Frow.which.Norm="sct")
test_meta_sct<-test_meta[row.names(Int_GSE134174_sct@meta.data),]
Int_GSE134174_sct<-Seurat::AddMetaData(object=Int_GSE134174_sct,metadata =as.data.frame(test_meta_sct))
Seurat::DimPlot(Int_GSE134174_sct,group.by="cluster_ident",label=TRUE,repel = TRUE)+ NoLegend()

## ----eval=TRUE,fig.height = 4, fig.width = 8, fig.align = "center"------------
# ##The codes were used to generate the example data list of GSE104154 (only a small subset of the data were used)
#   myGeoList<-c("GSE104154", "GSE161648")
#   myAllF<-Get_Geo_Data(DownGSEID=myGeoList,DecompressFile=TRUE,DecompressType="All")
#   dat_104_raw<-read.csv("GSE104154/GSE104154_d0_d21_sma_tm_Expr_raw.csv")
#   meta_104_2<-as.data.frame(readxl::read_xlsx("GSE104154/GSE104154_cell_type_annotation_d0_d21.xlsx",col_names=T,sheet=2))
#   meta_104_3<-as.data.frame(readxl::read_xlsx("GSE104154/GSE104154_cell_type_annotation_d0_d21.xlsx",col_names=T,sheet=3))
#   meta_104_2_3<-rbind(meta_104_2,meta_104_3)
#   meta_104_2_3$ID<-sapply(strsplit(meta_104_2_3$Barcode,split="\\."),function(x) x[2])
#   sleSize=300
#   sampleInd<-unlist(lapply(unique(meta_104_2_3$ID),function(x) {temp<-meta_104_2_3[meta_104_2_3$ID %in% x,]
#     set.seed(1234);sleSize<-min(nrow(temp),sleSize);SelectCellInd<-sample(1:nrow(temp),size=sleSize,replace=F)
#     returnBarcode<-temp$Barcode[SelectCellInd]}))
#   meta_23<-meta_104_2_3[meta_104_2_3$Barcode %in% sampleInd,]
#   dat_104Raw_withAno<-dat_104_raw[,paste(meta_23$Barcode)]
#   rowSum_dat104Raw<-rowSums(dat_104Raw_withAno)
#   dat_104Raw_f<-dat_104Raw_withAno[rowSum_dat104Raw>0,]
#   removeDupGene<-(!duplicated(dat_104_raw[rowSum_dat104Raw>0,2]))
#   dat_104Raw_f2<-dat_104Raw_f[removeDupGene,]
#   row.names(dat_104Raw_f2)<-dat_104_raw[rowSum_dat104Raw>0,2][removeDupGene]
#   t_datRaw<-data.frame(barcode=sapply(strsplit(colnames(dat_104Raw_f2),split="\\."),function(x) x[2]),t(dat_104Raw_f2))
#   t_datRaw$barcode<-dplyr::recode(t_datRaw$barcode,"1"="d0_sma.pos_tm.pos","2"="d0_sma.neg_tm.pos",
#   "3"="d0_sma.neg_tm.neg", "4"="d21_sma.pos_tm.pos","5"="d21_sma.neg_tm.pos","6"="d21_sma.neg_tm.neg")
#   RawDataList_GSE104154<-Splitdata_MakeDataList(InputData=t(t_datRaw[,2:ncol(t_datRaw)]),Group=t_datRaw$barcode)
#   meta_GSE104154<-meta_23
#
# ##The codes were used to generate the example data list of GSE161648
#   GEO_ID<-"GSE161648"
#   list_barcodefile<-list.files(GEO_ID,pattern="*_barcodes.tsv.gz")
#   (sampleName<-gsub("_barcodes.tsv.gz","",list_barcodefile))
#   sapply(sampleName,function(x){
#     NewSampleFolder<-sprintf("%s/%s",GEO_ID,x)
#     if (!file.exists(NewSampleFolder)){dir.create(NewSampleFolder)}
#     file.copy(sprintf("%s/%s_barcodes.tsv.gz",GEO_ID,x), sprintf("%s/barcodes.tsv.gz",NewSampleFolder))
#     file.copy(sprintf("%s/%s_matrix.mtx.gz",GEO_ID,x), sprintf("%s/matrix.mtx.gz",NewSampleFolder))
#     file.copy(sprintf("%s/%s_features.tsv.gz",GEO_ID,x), sprintf("%s/features.tsv.gz",NewSampleFolder))
#   })
#   GSE161648_10xRawData<-sapply(sampleName,function(x){rawdata<- Seurat::Read10X(data.dir =sprintf("%s/%s",GEO_ID,x))})
#   sleSize1=500
#   RawDataList_GSE161648<-lapply(GSE161648_10xRawData,function(x){sleSize1<-min(ncol(as.matrix(x)),sleSize1)
#     set.seed(1234);SelectCellInd<-sample(1:ncol(as.matrix(x)),size=sleSize1,replace=F)
#     returnCell<-(as.matrix(x))[,SelectCellInd]})

  ##load the example data
  data(Dat_GSE54_GSE48)
  RawDataList_GSE104154<-Dat_GSE54_GSE48$dat_GSE104154$RawDataList
  meta_GSE104154<-Dat_GSE54_GSE48$dat_GSE104154$MetaData
  row.names(meta_GSE104154)<-meta_GSE104154$Barcode
  RawDataList_GSE161648<-Dat_GSE54_GSE48$dat_GSE161648$RawDataList
  ###make Seurat object fro the raw data
  seu_GSE104154_raw<-MakeSeuObj_FromRawRNAData(RawList=RawDataList_GSE104154,
                                               GSE.ID="GSE104154",
                                               MtPattern='^mt-',
                                               MinFeature=1,
                                               MaxFeature=750000,
                                               MinCount=1, 
                                               MaxCount=4000000,
                                               MaxMT=100)
  seu_GSE161648_raw<-MakeSeuObj_FromRawRNAData(RawList=RawDataList_GSE161648,GSE.ID="GSE161648",MtPattern='^mt-')
  ###Data integration
  Int_GSE54_GSE48<-SeuObj_integration(Object.list=seu_GSE104154_raw,Object.list2=seu_GSE161648_raw)
  
  ###Evaluate the potential batch effect
  Seurat::DimPlot(Int_GSE54_GSE48,group.by="orig.ident",split.by="GSE",repel = TRUE)
  Seurat::DimPlot(Int_GSE54_GSE48,group.by="seurat_clusters",split.by="GSE",label=TRUE,repel = TRUE)+ NoLegend()
  
  ###Add annotation to GSE104154 as provided by the data submitter
  Int_GSE54_GSE48$rowname<-row.names(Int_GSE54_GSE48@meta.data)
  meta_subsetOfInt_GSE54<-subset(Int_GSE54_GSE48,GSE=="GSE104154")@meta.data
  meta_subsetOfInt_GSE54$Barcode<-sapply(strsplit(meta_subsetOfInt_GSE54$rowname,split="\\_"),function(x) x[1])
  ano_meta_subsetOfInt_GSE54<-dplyr::left_join(meta_subsetOfInt_GSE54,meta_GSE104154,by="Barcode")
  row.names(ano_meta_subsetOfInt_GSE54)<-ano_meta_subsetOfInt_GSE54$rowname
  ano_meta_subsetOfInt_GSE54<-ano_meta_subsetOfInt_GSE54[,c("rowname","defined")]
  Int_GSE54_GSE48@meta.data<-dplyr::left_join(Int_GSE54_GSE48@meta.data,ano_meta_subsetOfInt_GSE54,by="rowname")
  row.names(Int_GSE54_GSE48@meta.data)<-Int_GSE54_GSE48@meta.data$rowname
  Int_GSE54_GSE48$defined<-ifelse(is.na(Int_GSE54_GSE48$defined),
                                  paste(Int_GSE54_GSE48$seurat_clusters),Int_GSE54_GSE48$defined)

  Seurat::DimPlot(Int_GSE54_GSE48,group.by="defined",split.by="GSE",label=TRUE,repel = TRUE)+ NoLegend()
  ###Test whether loss Fas signaling will affect PDGFrb cell fraction in GSE161648 dataset
  Subset_GSE161648FromInt<-subset(Int_GSE54_GSE48,GSE=="GSE161648")
  Subset_GSE161648FromInt$FasStatus<-ifelse(Subset_GSE161648FromInt$ID %in% c("GSM4911963_sixwkfasdel","GSM4911965_threewkfasdel"),"Loss","NoLoss")
  tab_dat<-as.matrix(as.data.frame(unclass(table(Subset_GSE161648FromInt@meta.data[,c("FasStatus","defined")]))))
  ###PDGFrb is located in the cluster of 10. No obvious different is observed and loss of Fas signaling probably won't affect PDGFrb cell fraction.
  (tab_pop<-prop.table(tab_dat, margin = 1))
  dat_barplot<-data.frame(FasStatus=rep(c("Loss","NoLoss"),ncol(tab_pop)),
                        defined=rep(colnames(tab_pop),each=2),
                        Percentage=as.numeric(tab_pop))

## ----eval=TRUE,fig.align = "center"-------------------------------------------
  ###Barplot. Note that PDGFrb is located in the cluster of 7.
  pl<-ggplot2::ggplot(dat_barplot, aes(fill=defined, y=Percentage, x=FasStatus))+
    geom_bar(position="stack", stat="identity") + 
    theme_bw()
  print(pl)

## ----eval=TRUE----------------------------------------------------------------
   data(Data_Ref_transfer)
   data(Data_Query_transfer)

  Transfer_Result<-SeuObj_transfer(RefObj=Data_Ref_transfer,
                                   QueryObj=Data_Query_transfer, 
                                   RefLabel="defined")
  Seurat::DimPlot(Transfer_Result, 
                  group.by="predicted.RefLabel", 
                  reduction = "umap",
                  label=TRUE,
                  repel = TRUE)+ NoLegend()

## -----------------------------------------------------------------------------
sessionInfo()

