#' @title Get_GSMFromGSE
#'
#' @description Get GSM information matching to the GSE IDs
#'
#'
#' @param GeoDataBase The object returned by getGeoMetaDatabase function. Required
#' @param GSE.ID The GSE IDs
#'
#' @return The GSM information
#' @export
#' @examples
#'
#' library(singleGEO)
#' MyDataBase<-GetGeoMetaDatabase(Sqlfile=NULL,Demo=TRUE)
#' GSEID<-c("GSE122960","GSE158127")
#' GSM_info<-Get_GSMFromGSE(GeoDataBase=MyDataBase,GSE.ID=GSEID)
#'


Get_GSMFromGSE<-function(GeoDataBase=GeoDataBase,
                         GSE.ID=c("GSE122960","GSE158127")){
  output<-Obtain_GSM_fromGSE(Gseid=GSE.ID,db=GeoDataBase)
  return(output)
}


