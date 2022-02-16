#' @title SeuObj_integration
#'
#' @description Make Seurat object from the normalized data which were downloaded from GEO database. To use this function, a list containing the normalized data should be provided. This function is re-written from Seurat package.
#'
#' @param Object.list A list containing the seurat object of each samples. This should be returned by MakeSeuObj_FromNormalizedData or MakeSeuObj_FromRawData. Required.
#' @param Object.list2 An additional list containing the seurat object of each samples for data integration. This should be returned by MakeSeuObj_FromNormalizedData or MakeSeuObj_FromRawData. Optional: NULL.
#' @param Frow.which.Norm Which normalization method the data are obtained. Two normalization methods are allowed and supported: lognorm (For LogNormalization) and sct(For sctransform). Default: lognorm
#' @param SampleNameAsReference A string of names to be used as the reference for data integration. The names should match the object names of the input list. Try the command of "names(object.list)" and/or "names(object.list2)" to check. Optional: NULL.
#' @param NumberOfSampleForReference The number of objects will be used as the reference during integration. This parameter will be masked if SampleNameAsReference is not NULL. Default: 4
#' @param NumberOfSampleForReference2 The number of objects in Object.list2 will be used as the reference during integration. If this parameter is null and object.list2 is not null, it is set to be equal of SampleNameAsReference. Default: NULL
#' @param Nfeatures The number of top variable features for the SelectIntegrationFeatures function. Default: 2000
#' @param Anchor.reduction Dimensional reduction method to find anchors. Default: rpca
#' @param Do.scale Whether to scale the data or not. Binary data. This will feed to ScaleData function. Default: TRUE
#' @param Do.center Whether to center the data or not. Binary data. This will feed to ScaleData function. Default: TRUE
#' @param Dims.anchor The number of top dimensions for the IntegrateData function. Default: 1:50
#' @param Dims.umap The number of top dimensions of reduction to use for the functions of FindNeighbors and RunUMAP. Default: 1:30
#' @param K.weight Number of neighbors to consider when weighting anchors in "IntegrateData" function. Default: 100
#' @param Resolution The resolution value for FindClusters function. Default: 0.8
#' @param Algorithm The algorithm to be used in FindClusters. Default: 1
#' @param Future.globals.maxSize Maximum allowed total size of global variables. Default: 12000*1024^2 (12Gb)
#'
#' @return A Seurat object with data integrated.
#' @export
#' @examples
#'
#' library(singleGEO)
#'
#' #
#' data(listDat_GSE134174)
#' listDat_obj<-listDat_GSE134174$data
#' listDat_meta<-listDat_GSE134174$meta
#' Int_GSE134174<-SeuObj_integration(Object.list=listDat_obj)
#' row.names(listDat_meta)<-listDat_meta$Cell
#' listDat_meta<-listDat_meta[row.names(Int_GSE134174@meta.data),]
#' Int_GSE134174<-Seurat::AddMetaData(object=Int_GSE134174,metadata =as.data.frame(listDat_meta))
#' Seurat::DimPlot(Int_GSE134174,group.by="cluster_ident",label=TRUE)

SeuObj_integration<-function(Object.list = object.list,
                             Object.list2 = NULL,
                             Frow.which.Norm="lognorm",
                             SampleNameAsReference=NULL,
                             NumberOfSampleForReference=4,
                             NumberOfSampleForReference2=NULL,
                             Nfeatures = 2000,
                             Do.scale = TRUE,
                             Do.center = TRUE,
                             Anchor.reduction = 'rpca',
                             Dims.anchor=1:50,
                             Dims.umap=1:30,
                             K.weight = 100,
                             Resolution=0.8,
                             Algorithm=1,
                             Future.globals.maxSize = 12000*1024^2){
  options(future.globals.maxSize = Future.globals.maxSize)

  ###if neither NumberOfSampleForReference or SampleNameAsReference
  if(is.null(NumberOfSampleForReference) & is.null(SampleNameAsReference)){
    NumberOfSampleForReference=5000
    print("NumberOfSampleForReference is NULL; Set to 5000!")}

  ###if provide NumberOfSampleForReference
  if(!is.null(NumberOfSampleForReference)){
    referenceSampleN<-SelectObjWithTopNumCluster(ObjList=Object.list,TopN=NumberOfSampleForReference)

    object.list.all<-Object.list
    res_checkObj2<-checkObj2_forInt(Obj1=object.list.all,
                                    Obj2=Object.list2,
                                    RefSamN=referenceSampleN,
                                    NumOfSamForRef=NumberOfSampleForReference,
                                    NumOfSamForRef2=NumberOfSampleForReference2)
    object.list.all<-res_checkObj2$Obj.all
    referenceSampleN<-res_checkObj2$RefSamN.all

    ReferenceIndex<-checkRefIndex(obj.all=object.list.all,
                                  refSamN=referenceSampleN)
  }

  ###if provide SampleNameAsReference
  if(!is.null(SampleNameAsReference)){
    object.list.all<-Object.list
    if(!is.null(Object.list2)){
      object.list.all<-c(Object.list,Object.list2)
    }
    ReferenceIndex<-checkRefIndex(obj.all=object.list.all,
                                  refSamN=SampleNameAsReference)
  }



  features <- Seurat::SelectIntegrationFeatures(object.list = object.list.all, nfeatures = Nfeatures)

  if(Frow.which.Norm %in% c("lognorm","Lognorm","LogNorm","LogNormalize")){
    object.list.all <- lapply(object.list.all, FUN = function(x) {
      x <- Seurat::ScaleData(x, features = features,do.scale = Do.scale,do.center = Do.center, verbose = FALSE)
      x <- Seurat::RunPCA(x, features = features, verbose = FALSE)
    })

    anchors <- Seurat::FindIntegrationAnchors(object.list = object.list.all,
                                              reference =ReferenceIndex,
                                              normalization.method = "LogNormalize",
                                              anchor.features = features,
                                              reduction = Anchor.reduction,
                                              dims = Dims.anchor,
                                              verbose = FALSE)
    seu_int <- Seurat::IntegrateData(anchorset = anchors,
                                     dims = Dims.anchor,
                                     k.weight = K.weight,
                                     normalization.method ="LogNormalize",
                                     verbose = FALSE)

    Seurat::DefaultAssay(seu_int) <- "integrated"
    seu_int <- Seurat::ScaleData(object = seu_int,do.scale = Do.scale,do.center = Do.center,verbose = FALSE)
    seu_int <- Seurat::RunPCA(object = seu_int,verbose = FALSE)
    seu_int <- Seurat::FindNeighbors(seu_int, dims = Dims.umap,verbose = FALSE)
    seu_int <- Seurat::FindClusters(seu_int,resolution=Resolution,algorithm=Algorithm,verbose = FALSE)
    seu_int <- Seurat::RunUMAP(object = seu_int, dims = Dims.umap,return.model =TRUE,verbose = FALSE)
  }else if(Frow.which.Norm %in% c("sct","sctransform","SCT","Sctransform","Sct")){
    sct.list <- Seurat::PrepSCTIntegration(object.list =object.list.all, anchor.features = features)
    sct.anchors <- Seurat::FindIntegrationAnchors(object.list = sct.list,reference =ReferenceIndex, normalization.method = "SCT",anchor.features = features, reduction = Anchor.reduction, dims = Dims.anchor,verbose = FALSE)
    seu_int <- Seurat::IntegrateData(anchorset = sct.anchors,dims = Dims.anchor, k.weight = K.weight, normalization.method = "SCT",verbose = FALSE)
    seu_int <- Seurat::RunPCA(seu_int, verbose = FALSE)
    seu_int <- Seurat::RunUMAP(seu_int, dims = Dims.umap,return.model =TRUE,verbose = FALSE)
  }else{
    stop("Norm.method should be 'lognorm' or 'sct'")
  }

  return(seu_int)
}
