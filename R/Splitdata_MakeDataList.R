#' @title Splitdata_MakeDataList
#'
#' @description Split the data by a factor and make a list
#'
#' @param InputData The input matrix or data frame. Each row is a gene and each column is a cell ID. Required
#' @param Group The factor to split the matrix/data frame
#'
#' @return A list of matrix split by Group variable
#' @export
#' @examples
#'
#' test.matrix<-matrix(sample(0:3,size=900,replace=TRUE),nrow=100,ncol=90)
#' colnames(test.matrix)<-paste0("Cell",1:90,sep="")
#' row.names(test.matrix)<-paste0("Gene",1:100,sep="")
#' sparseMatrix<-as(test.matrix,"sparseMatrix")
#' pdata<-data.frame(CellID=paste0("Cell",1:90,sep=""),CellGroup=rep(c("A","B","C"),c(20,30,40)))
#' list_test<-Splitdata_MakeDataList(InputData=sparseMatrix,Group=pdata$CellGroup)


Splitdata_MakeDataList<-function(InputData,Group=group){
  if(!length(Group)==ncol(InputData)){
    stop("The length of 'Group' factor is not equal to column length of input data!")
  }else{
    uniqG<-unique(Group)
    OutList<-list()
    for(i in uniqG){
      whichCol<-which(Group %in% i)
      OutList[[i]]<-InputData[,whichCol]}
  }
  return(OutList)
}

