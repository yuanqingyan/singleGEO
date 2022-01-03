#' @title ExtractFile
#'
#' @description Decompress the files
#'
#' @param ExtractID The list of GSE IDs. Required
#' @param DecompressType Which types of files to be decompressed if DecompressFile is set to TRUE. Options: "gz"(for .gz),"tar"(for .tar or .tar.gz),"zip"(for .zip) or "All"(for .tar, .tar.gz, .gz, .zip). Optional. Default: "All"
#'
#' @return NULL
#' @export
#' @examples
#'
#' library(singleGEO)
#' ###not run
#' ##myGeoID<-c("GSE131800")
#' ##DownloadFileInfo<-Get_Geo_Data(DownGSEID=myGeoID,DecompressFile=FALSE)
#' ##ExtractFile(ExtractID=myGeoID,DecompressType="tar")


ExtractFile<-function(ExtractID=ExtractID,DecompressType=decompressType){
  if(sum(DecompressType %in% c("All","all","ALL"))>=1){
    extractAll(ExtractID=ExtractID)
  }
  if(sum(DecompressType %in% c("Tar","tar","TAR","Tar.gz","tar.gz","Tar.Gz","TAR.GZ"))>=1){
    extractTar(ExtractID=ExtractID)
  }
  if(sum(DecompressType %in% c("Gz","gz","GZ"))>=1){
    extractGZ(ExtractID=ExtractID)
  }
  if(sum(DecompressType %in% c("Zip","zip","ZIP"))>=1){
    extractZip(ExtractID=ExtractID)
  }
}


