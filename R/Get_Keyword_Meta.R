#' @title Get_Keyword_Meta
#'
#' @description Filter the meta database to obtain the study using single cell RNAseq data and containing the key information provided by user
#'
#' @param GeoDataBase The object returned by getGeoMetaDatabase function. Required
#' @param Organism The organisms to search for. One or multiple organisms can be provided. Optional. Default: c("Homo sapiens")
#' @param InputSearchList The list of keywords which are described in GEO database. Multiple keywords can be provided and the keywords should be in a list. Required. Default: c("single")
#' @param Platform The platform of single cell sequencing. Currently only support single cell RNAseq(scRNAseq) and single cell ATACseq(scATACseq). Options: scRNAseq, scATACseq. Default: "scRNAseq"
#'
#' @return The filtered meta data
#' @export
#' @examples
#'
#' library(singleGEO)
#' MyDataBase<-GetGeoMetaDatabase(Sqlfile=NULL,Demo=TRUE)
#' MyOrganism<-c("Mus musculus","Homo sapiens")
#' search_kw1<-c("lung","pulmonary")
#' search_kw2<-c("fibrosis")
#' KeyList<-list(kw1=search_kw1,kw2=search_kw2)
#' MySearch<-Get_Keyword_Meta(GeoDataBase=MyDataBase,Organism=MyOrganism, InputSearchList=KeyList)
#'


Get_Keyword_Meta<-function(GeoDataBase=GeoDataBase,
                           Organism="Homo sapiens",
                           InputSearchList=list(kw1="single"),
                           Platform="scRNAseq"){
  get_organism<-Try_Filter_organism(keyword=Organism)
  get_scPlatform<-Try_GetPlatform(Platform=Platform)
  get_userKW<-Try_Filter_userKW(kwlist=InputSearchList)

  print("Database searching")
  ##step 1## filter out organism; obtain gsm#####
  CM_Obtain_GSM<-sprintf("SELECT gsm FROM gsm WHERE %s",get_organism)
  Obtain_GSM<-RSQLite::dbGetQuery(GeoDataBase,CM_Obtain_GSM)
  Input_gsm<-paste(sapply(Obtain_GSM$gsm,function(x) sprintf("'%s'",x)),collapse=",")
  ##step 2## obtain unique gse#####
  CM_Obtain_GSE_FromGSM<-sprintf("SELECT gse FROM gse_gsm WHERE gsm IN (%s)",Input_gsm)
  Obtain_GSE_FromGSM<-RSQLite::dbGetQuery(GeoDataBase,CM_Obtain_GSE_FromGSM)
  Input_gse<-paste(sapply(unique(Obtain_GSE_FromGSM$gse),function(x) sprintf("'%s'",x)),collapse=",")

  Unique_GSE_FromGSm_OrganismFilter<-sprintf("(gse IN (%s))",Input_gse)
  ##step 3## Final step to obtain gse information#####
  CM_Obtain_GSE_FromGSE<-sprintf("SELECT %s FROM gse WHERE %s AND %s AND %s",select_OutGSES_item(), Unique_GSE_FromGSm_OrganismFilter,get_scPlatform,get_userKW)
  res_query<-RSQLite::dbGetQuery(GeoDataBase,CM_Obtain_GSE_FromGSE)
  return(res_query)
}


