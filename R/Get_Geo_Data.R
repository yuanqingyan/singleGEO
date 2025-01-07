#' @title Get_Geo_Data
#'
#' @description Download all supplementary data of GSE IDs provided by the user
#'
#' @param DownGSEID The GSE IDs for data downloading. Required
#' @param DecompressFile Whether decompress the files or not. Binary variable with TRUE or FALSE. Optional. Default: FALSE
#' @param DecompressType Which types of files to be decompressed if DecompressFile is set to TRUE. Options: "gz"(for .gz),"tar"(for .tar or .tar.gz),"zip"(for .zip) or "All"(for .tar, .tar.gz, .gz, .zip). Optional. Default: "All"
#'
#' @return The list of files downloaded
#' @export
#' @examples
#'
#' library(singleGEO)
#' ###not run
#' ##myGeoID<-c("GSE131800")
#' ##DownloadFileInfo<-Get_Geo_Data(DownGSEID=myGeoID,DecompressFile=FALSE)


Get_Geo_Data<-function(DownGSEID=DownGSEID,DecompressFile=FALSE,DecompressType="All"){
  sapply(DownGSEID,function(iGSEID){
    # AlreadyFile<-CheckFolderExist(checkList=iGSEID)
    # if(length(AlreadyFile)>0){
    #   removeFolder(ToBeRemove=AlreadyFile)
    # }

    downloadGSESupplFile(GSEID=iGSEID,destfolder=getwd())

    if(DecompressFile){
      ExtractFile(ExtractID=iGSEID,DecompressType=DecompressType)
    }
  })


  ListOfAllFile<-NULL

  if(length(DownGSEID)>0){
    ListOfAllFile<-sapply(DownGSEID,function(x){
      list.files(path=paste(x))
    })
  }
  return(ListOfAllFile)
}


