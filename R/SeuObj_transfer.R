#' @title SeuObj_transfer
#'
#' @description Transfer reference data information onto query data.
#'
#' @param RefObj A reference Seurat object. Required.
#' @param QueryObj The query Seurat object. Required.
#' @param RefLabel The column name of meta data from RefObj. This label will transfer to query data. Required.
#' @param Norm.method Normalization method. The standard LogNormalize or sctransform. Options of: "lognorm" or "sct". Default: "lognorm"
#' @param Future.globals.maxSize Maximum allowed total size of global variables. Default: 12000*1024^2 (12Gb)
#'
#'
#' @return The Seurat object of query data with predicted labels from reference data.
#' @export
#' @examples
#'
#' library(singleGEO)
#' ###
#'
#' data(Data_Ref_transfer)
#' data(Data_Query_transfer)
#'
#' Transfer_Result<-SeuObj_transfer(RefObj=Data_Ref_transfer,
#'                                 QueryObj=Data_Query_transfer,
#'                                  RefLabel="defined")
#' Seurat::DimPlot(Transfer_Result,
#'                 group.by="predicted.RefLabel",
#'                 reduction = "umap",
#'                 label=TRUE,
#'                 repel = TRUE)+ NoLegend()


SeuObj_transfer<-function(RefObj=RefObj,
                          QueryObj=QueryObj,
                          RefLabel=RefLabel,
                          Norm.method = "lognorm",
                          Future.globals.maxSize = 12000*1024^2){
  options(future.globals.maxSize = Future.globals.maxSize)
  if(Norm.method %in% c("lognorm","Lognorm","LogNorm","LogNormalize")){
    normalization.method<-"LogNormalize"
  }else if(Norm.method %in% c("sct","sctransform","SCT","Sctransform","Sct")){
    normalization.method<-"SCT"
  }else{
    stop("Norm.method should be 'lognorm' or 'sct'")
  }

   FindAnchor <- Seurat::FindTransferAnchors(reference =RefObj,
                                            query = QueryObj,
                                            normalization.method = normalization.method,
                                            dims = 1:30,
                                            reference.reduction ="pca",
                                            verbose = FALSE)

  TransferObj <- Seurat::MapQuery(anchorset = FindAnchor,
                                  reference = RefObj,
                                  query =QueryObj,
                                  refdata = list(RefLabel= RefLabel),
                                  reference.reduction ="pca",
                                  reduction.model ="umap",
                                  verbose = FALSE)
  return(TransferObj)
}
