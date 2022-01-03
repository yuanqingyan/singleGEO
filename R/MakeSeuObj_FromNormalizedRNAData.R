#' @title MakeSeuObj_FromNormalizedRNAData
#'
#' @description Make Seurat object from the normalized data which were downloaded from GEO database. To use this function, a list containing the normalized data should be provided. This function is re-written from Seurat package.
#'
#' @param NormList A list containing the normalized data. Each element of the list is for one normalized dataset. The name of the element will be used as the ID of the normalized dataset. Required.
#' @param GSE.ID The GSE ID for this project. Default: Test
#' @param Frow.which.Norm Which normalization method the data are obtained. Two normalization methods are allowed and supported: lognorm (For LogNormalization) and sct(For sctransform). Default: lognorm
#' @param Feature.selection.method The method for the top variable feature selection. This will feed to FindVariableFeatures function. Default: vst
#' @param Nfeatures The number of top variable features for the FindVariableFeatures function. Default: 2000
#' @param Do.scale Whether to scale the data or not. Binary data. This will feed to ScaleData function. Default: TRUE
#' @param Do.center Whether to center the data or not. Binary data. This will feed to ScaleData function. Default: TRUE
#' @param Dims The number of top dimensions of reduction to use for the functions of FindNeighbors and RunUMAP. Default: 1:30
#' @param Resolution The resolution value for FindClusters function. Default: 0.8
#' @param Algorithm The algorithm to be used in FindClusters. Default: 1
#'
#' @return A list with each element corresponding to a Seurat object of a normalized data. The object is ready for the downstream analysis by Seurat package.
#' @export
#' @examples
#'
#' library(singleGEO)
#'
#' #
#' data(Norm_GSE134174)
#' seu_Norm_GSE134174<-MakeSeuObj_FromNormalizedRNAData(NormList=Norm_GSE134174,GSE.ID="GSE134174")

MakeSeuObj_FromNormalizedRNAData<-function(NormList=NormList,
                                        GSE.ID="Test",
                                        Frow.which.Norm="lognorm",
                                        Feature.selection.method = "vst",
                                        Nfeatures = 2000,
                                        Do.scale = TRUE,
                                        Do.center = TRUE,
                                        Resolution=0.8,
                                        Algorithm=1,
                                        Dims = 1:30){
  name_NormList<-names(NormList)
  listOut<-sapply(name_NormList,function(x){
    temp_data<-CheckRowName(DataInput=NormList[[x]])
    seu_Nor<- Seurat::CreateSeuratObject(counts =temp_data,
                                         project =x,
                                         min.cells = 1,
                                         min.features = 1)
    seu_Nor[['ID']] <- x
    seu_Nor[['GSE']] <- GSE.ID

    seu_Nor <-CheckNormMethodForNormData(Seu_obj_input=seu_Nor,
                                         Feature.selection.method=Feature.selection.method,
                                         Nfeatures = Nfeatures,Dims = Dims,Do.scale = Do.scale,Do.center = Do.center)

    seu_Nor <-Seurat_dimensionReduc(Seu_obj_input=seu_Nor,
                                    Dims = Dims,
                                    Resolution=Resolution,
                                    Algorithm=Algorithm)
  })
  return(listOut)
}
