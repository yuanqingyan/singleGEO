#' @title GetGeoMetaDatabase
#'
#' @description Download and connect to a SQLite database file
#'
#' @param Sqlfile Using existing GEO meta database or downloading the new one. If "GEOmetadb.sqlite" file has been downloaded but it is not in current folder, the path of the file should be provided. If file of "GEOmetadb.sqlite" can be detected under current folder, the program will use the file. If no "GEOmetadb.sqlite"" file is detected, the program will download the new one and save it in current folder. Optional. Default: NULL(check file under current folder or automatically download the new one)
#' @param Demo Whether to use Demo database or not. Default: FALSE
#'
#' @return The meta database
#' @export
#' @examples
#'
#' library(singleGEO)
#' GeoDataBase<-GetGeoMetaDatabase(Sqlfile=NULL,Demo=TRUE)

GetGeoMetaDatabase<-function(Sqlfile=NULL,Demo=FALSE){
  GeoDataBase=NULL
  if(Demo==TRUE & is.null(Sqlfile)){
    print("Using demo GEO meta database")
    DemoFileLocation<-system.file("extdata","demo.GEOmetadb.sqlite", package = "singleGEO")
    GeoDataBase<-RSQLite::dbConnect(RSQLite::SQLite(), DemoFileLocation)
  }else{
    if(!is.null(Sqlfile)) {
      print("Using GEO meta database provided by user")
      GEOSqlfile<-Sqlfile
    } else if (file.exists("GEOmetadb.sqlite")){
      print("Using GEO meta database within current folder")
      GEOSqlfile<- "GEOmetadb.sqlite"
    }else{
      file.remove("GEOmetadb.sqlite.gz")
      print("Downloading metadb from GEO....")
      #GEOSqlfile<- GEOmetadb::getSQLiteFile()
      GEOSqlfile<- singleGEOgetSQLiteFile()
    }
    GeoDataBase<-RSQLite::dbConnect(RSQLite::SQLite(), GEOSqlfile)
  }
  return(GeoDataBase)
}

