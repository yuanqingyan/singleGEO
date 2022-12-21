#' @title MakeSeuObj_FromRawRNAData
#'
#' @description Make Seurat object from the raw count data which were downloaded from GEO database. To use this function, a list containing the raw data should be provided. This function is re-written from Seurat package.
#'
#'
#' @param RawList A list containing the raw data. Each element of the list is for one raw dataset. The name of the element will be used as the ID of the raw dataset. Required.
#' @param MtPattern The pattern to recognize mitochondria gene. Different species can have different labels. For example, for human data, MtPattern='^MT-' is appropriate, while for mouse it should be MtPattern='^mt-'. Optional. Default: MtPattern='^MT-'
#' @param GSE.ID The GSE ID for this project. Default: Test
#' @param MinFeature Minimal feature count to keep the cell. Default: 200
#' @param MaxFeature Maximal feature count to keep the cell. Default: 7500
#' @param MinCount Minimal RNA count to keep the cell. Default: 400
#' @param MaxCount Maximal RNA count to keep the cell. Default: 40000
#' @param MaxMT Percentage of mitochondria gene to filter out the cell. Default: 10
#' @param Norm.method The method to normalize data. Two normalization methods are allowed and supported: lognorm (For LogNormalization) and sct(For sctransform). Default: lognorm
#' @param Do.scale Whether to scale the data or not. Binary data. This will feed to ScaleData function. Default: TRUE
#' @param Do.center Whether to center the data or not. Binary data. This will feed to ScaleData function. Default: TRUE
#' @param Scale.factor The factor to scale up the data. Default: 10000
#' @param Feature.selection.method The method for the top variable feature selection. This will feed to FindVariableFeatures function. Default: vst
#' @param Nfeatures The number of top variable features for the FindVariableFeatures function. Default: 2000
#' @param Dims The number of top dimensions of reduction to use for the functions of FindNeighbors and RunUMAP. Default: 1:30
#' @param Npcs The number of pc to use for the functions of RunPCA. Default: 50
#' @param Resolution The resolution value for FindClusters function. Default: 0.8
#' @param Algorithm The algorithm to be used in FindClusters. Default: 1
#'
#' @return A list with each element corresponding to a Seurat object of a raw data.
#' @export
#' @examples
#'
#' library(singleGEO)
#'
#' #
#'
#' data(testData_GSE134174)
#'
#' test_dat<-testData_GSE134174$TwoRawData
#' test_meta<-testData_GSE134174$TwoMetaData
#' list_GSE134174<-Splitdata_MakeDataList(InputData=test_dat,Group=test_meta$Donor)
#'
#' ##The following code setting big number of MaxFeature,MaxCount,MaxMT and small number
#' ##of MinFeature,MinCount will keep all the cells
#'
#' seu_GSE134174<-MakeSeuObj_FromRawRNAData(RawList=list_GSE134174,GSE.ID="GSE134174",MinFeature=1,MaxFeature=750000,
#'                                        MinCount=1, MaxCount=4000000,MaxMT=100)



MakeSeuObj_FromRawRNAData<-function(RawList=RawList,MtPattern='^MT-',GSE.ID="Test",MinFeature=200,MaxFeature=7500,
                                 MinCount=400, MaxCount=40000,MaxMT=10,
                                 Norm.method = "lognorm", Scale.factor = 10000,
                                 Feature.selection.method = "vst", Nfeatures = 2000,Npcs =50,Dims = 1:30,
                                 Resolution=0.8,Algorithm=1,
                                 Do.scale = TRUE,Do.center = TRUE){
  name_RawList<-names(RawList)
  listOut<-sapply(name_RawList,function(x){
    temp_data<-CheckRowName(DataInput=RawList[[x]])
    MtPattern<-gsub("-",".",MtPattern);MtPattern<-gsub("_",".",MtPattern);MtPattern<-gsub("\\.","\\\\.",MtPattern)
    seu_raw<- Seurat::CreateSeuratObject(counts =temp_data,
                                         project =x,
                                         min.cells = 1,
                                         min.features = 1)
    seu_raw[['percent.mt']] <- Seurat::PercentageFeatureSet(seu_raw, pattern=MtPattern)
    seu_raw[['ID']] <- x
    seu_raw[['GSE']] <- GSE.ID
    seu_sub <- subset(seu_raw,
                      subset = nFeature_RNA > MinFeature &
                        nFeature_RNA < MaxFeature &
                        nCount_RNA > MinCount &
                        nCount_RNA < MaxCount &
                        percent.mt < MaxMT)

    seu_sub <-Seurat_normalization(Seu_obj_input=seu_sub,
                                   Norm.method=Norm.method,
                                   Scale.factor=Scale.factor,
                                   Feature.selection.method=Feature.selection.method,
                                   Nfeatures = Nfeatures,
                                   Do.scale = Do.scale,
                                   Do.center = Do.center)

    seu_sub <-Seurat_dimensionReduc(Seu_obj_input=seu_sub,
                                    npcs =Npcs,
                                    Dims = Dims,
                                    Resolution=Resolution,
                                    Algorithm=Algorithm)

  })
  return(listOut)
}

