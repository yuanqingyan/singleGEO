# MakeSinSqlite<-function(FilePath=FilePath){
#   geometadbfile <- GEOmetadb::getSQLiteFile()
#   FullData<-GetGeoMetaDatabase(Sqlfile="GEOmetadb.sqlite")
#
#   get_organism<-Filter_organism(keyword=c("Mus musculus","Homo sapiens"))
#   get_platform<-GetPlatform(Platform=c("scRNAseq","scATACseq"))
#   cmd_query<-sprintf("SELECT %s FROM %s WHERE %s AND %s",
#                      select_gsmgse_item(),join_gse_gsm(),get_organism, get_platform)
#   res_query<-RSQLite::dbGetQuery(FullData,cmd_query)
#
#   unique_gse<-unique(res_query$gse)
#   length_gse<-length(unique_gse)
#   try_length<-900
#   len_ceiling<-ceiling(length_gse/try_length)
#
#   mydb<-RSQLite::dbConnect(RSQLite::SQLite(),dbname=sprintf("%s/singleGEO.sqlite",FilePath))
#   lapply(1:len_ceiling,function(x){
#     print("Build singleGEO.sqlite ......This step can take >1hr! Be patient!")
#     start<-(x-1)*try_length+1
#     end<-min(x*try_length,length_gse)
#
#     temp_gse<-unique_gse[start:end]
#     sleGSE<-as.character(sapply(temp_gse,function(y){sprintf("(gse.gse like '%s')",y)}))
#     sleGSE_collap<-paste(sleGSE,collapse=" OR ")
#     queryCMD_selDB_gse<-sprintf("SELECT gse.* FROM %s WHERE %s",join_gse_gsm(),sleGSE_collap)
#     queryCMD_selDB_gsm<-sprintf("SELECT gsm.* FROM %s WHERE %s",join_gse_gsm(),sleGSE_collap)
#     queryCMD_selDB_gse_gsm<-sprintf("SELECT gse.gse,gsm.gsm FROM %s WHERE %s",join_gse_gsm(),sleGSE_collap)
#
#
#     subset_gse<-RSQLite::dbGetQuery(FullData,queryCMD_selDB_gse)
#     subset_gsm<-RSQLite::dbGetQuery(FullData,queryCMD_selDB_gsm)
#     subset_gse_gsm<-RSQLite::dbGetQuery(FullData,queryCMD_selDB_gse_gsm)
#
#     RSQLite::dbWriteTable(conn = mydb, name = "gsm", subset_gsm, append=T,row.names=FALSE)
#     RSQLite::dbWriteTable(conn = mydb, name = "gse_gsm", subset_gse_gsm, append=T,row.names=FALSE)
#     RSQLite::dbWriteTable(conn = mydb, name = "gse", subset_gse, append=T,row.names=FALSE)
#   })
#
#   #print(sprintf("%s GSE records were found", nrow(RSQLite::dbGetQuery(mydb,cmd_query))))
#   RSQLite::dbDisconnect(mydb)
#   file.remove("GEOmetadb.sqlite")
# }

select_OutGSES_item<-function(){
  gse_item<-c("title", "gse", "summary","type" ,"overall_design","supplementary_file")
  return(paste(gse_item,collapse=","))
}



Try_gse_search_item<-function(){
  gseItem<-c("title","summary", "overall_design")
  return(gseItem)
}

Try_Filter_organism<-function(keyword){
  Fil_organism<-sprintf("(organism_ch1 IN (%s))",collapse_keywork(keyword))
  return(Fil_organism)
}

Try_Filter_sequencing<-function(){
  Fil_seq<-"(type LIKE '%high%throughput%sequencing%')"
  return(Fil_seq)
}

Try_Filter_SingleCell<-function(){
  singleCellItem<-keyword_singlecell_search()
  gseItem<-Try_gse_search_item()

  single_gse<-sapply(gseItem,function(y){
    single_like<-sapply(singleCellItem,function(x){sprintf("(%s LIKE '%%%s%%')",y,x)})
    single_collap<-paste(single_like,collapse=" OR ")
  })

  single_gse_add<-sapply(single_gse,function(z){sprintf("(%s)",z)})
  out_single_gse0<-paste(single_gse_add,collapse=" OR ")
  out_single_gse<-sprintf("(%s)",out_single_gse0) #####################(((title LIKE '%single-cell%') OR (title LIKE '%single cell%')
  return(out_single_gse)
}

Try_Filter_RNAseq<-function(){
  transcripItem<-keyword_RNAseq_search()
  gseItem<-Try_gse_search_item()

  RNA_gse<-sapply(gseItem,function(y){
    RNA_like<-sapply(transcripItem,function(x){sprintf("(%s LIKE '%%%s%%')",y,x)})
    RNA_collap<-paste(RNA_like,collapse=" OR ")
  })

  RNA_gse_add<-sapply(RNA_gse,function(z){sprintf("(%s)",z)})
  out_RNA_gse0<-paste(RNA_gse_add,collapse=" OR ")
  out_RNA_gse<-sprintf("(%s)",out_RNA_gse0)  #############(((title LIKE '%transcriptomes%') OR (title LIKE '%transcriptome%')
  return(out_RNA_gse)
}

Try_Filter_scRNAseq<-function(){
  fSingle<-Try_Filter_SingleCell()
  fRNA<-Try_Filter_RNAseq()
  fExp<-Try_Filter_sequencing()
  f_SingleRNA<-sprintf("%s AND %s AND %s",fExp,fSingle,fRNA)
  return(f_SingleRNA)
}

Try_Filter_ATACseq<-function(){
  ATAC_kw<-keyword_ATACseq_search()
  gseItem<-Try_gse_search_item()

  ATAC_gse<-sapply(gseItem,function(y){
    ATAC_like<-sapply(ATAC_kw,function(x){sprintf("(%s LIKE '%%%s%%')",y,x)})
    ATAC_collap<-paste(ATAC_like,collapse=" OR ")
  })

  ATAC_gse_add<-sapply(ATAC_gse,function(z){sprintf("(%s)",z)})
  out_ATAC_gse0<-paste(ATAC_gse_add,collapse=" OR ")
  out_ATAC_gse<-sprintf("(%s)",out_ATAC_gse0)
  return(out_ATAC_gse)
}

Try_Filter_scATACseq<-function(){
  fSingle<-Try_Filter_SingleCell()
  fATAC<-Try_Filter_ATACseq()
  fExp<-Try_Filter_sequencing()
  f_SingleATAC<-sprintf("%s AND %s AND %s",fExp,fSingle,fATAC)
  return(f_SingleATAC)
}

Try_GetPlatform<-function(Platform=Platform){
  RNA_platform<-NULL
  ATAC_platform<-NULL

  PlatformNotIn<-Platform[!Platform %in% c("scRNAseq","scATACseq")]
  if(length(PlatformNotIn)>0){
    stop(sprintf("Only support scRNAseq and/or scATACseq. Your %s is wrong!",PlatformNotIn))
  }

  if("scRNAseq" %in% Platform){
    RNA_platform<-Try_Filter_scRNAseq()
  }
  if("scATACseq" %in% Platform){
    ATAC_platform<-Try_Filter_scATACseq()
  }
  com_platform<-c(RNA_platform,ATAC_platform)
  com_platform<-com_platform[!is.null(com_platform)]
  if(length(com_platform)>1){
    com_platform_add<-sapply(com_platform,function(z){sprintf("(%s)",z)})
    out_platform0<-paste(com_platform_add,collapse=" OR ")
    out_platform<-sprintf("(%s)",out_platform0)
  }else if(length(com_platform)==1){
    out_platform<-com_platform
  }else{
    stop("No valid Platform is specified")
  }
  return(out_platform)
}


Try_Filter_userKW<-function(kwlist=kwlist){
  checkKW<-check_KW_list(user_kw=kwlist)
  nkw<-checkKW$nkw
  user_keyword<-checkKW$user_kw

  kw_name<-sprintf("user_kw%s",1:nkw)
  for(i in 1:nkw){assign(sprintf("user_kw%s",i),user_keyword[[i]])} #eval(parse(text=sprintf("kw%s",i)))
  gseItem<-Try_gse_search_item()

  user_wd_string<-sapply(kw_name,function(x1){
    temp_x1<-eval(parse(text=x1))
    temp_x2<-sapply(temp_x1,function(x2){
      temp_x2_1<-sapply(gseItem,function(y){
        temp_x2_2<-sapply(temp_x1,function(x3){sprintf("(%s LIKE '%%%s%%')",y,x3)})
        temp_x2_3<-paste(temp_x2_2,collapse=" OR ")})
    })
    temp_x2_add<-sapply(temp_x2,function(z){sprintf("(%s)",z)})
    out_temp_x20<-paste(temp_x2_add,collapse=" OR ")
    out_temp_x2<-sprintf("(%s)",out_temp_x20)
  })

  user_wd_string_add<-sapply(user_wd_string,function(z){sprintf("(%s)",z)})
  out_user_wd_string0<-paste(user_wd_string_add,collapse=" AND ")
  out_user_wd_string<-sprintf("(%s)",out_user_wd_string0)
  return(out_user_wd_string)
}


#
# select_gsmgse_item<-function(){
#   gsmgse_item<-c("gsm.title", "gsm.gsm", "gsm.source_name_ch1", "gsm.characteristics_ch1",
#                  "gsm.treatment_protocol_ch1", "gsm.treatment_protocol_ch2",
#                  "gse.title", "gse.gse", "gse.summary", "gse.overall_design","gse.pubmed_id","gse.supplementary_file")
#   return(paste(gsmgse_item,collapse=","))
# }
collapse_keywork<-function(keywork){
  temp<-sapply(keywork,function(x) sprintf("'%s'",x))
  temp_collap<-paste(temp,collapse=",")
  return(temp_collap)
}

keyword_singlecell_search<-function(){
  singleCellSearchItem<-c('single-cell','single cell','single cells',
                    'single-nuclei','single nuclei','single-nucleus','single nucleus')
  return(singleCellSearchItem)
}

keyword_RNAseq_search<-function(){
  RNAseqSearchItem<-c('transcriptomes','transcriptome','transcriptomic','transcriptomics',
                      'RNA-Seq','RNA-sequencing','RNAseq','RNAsequencing','RNA seq','RNA sequencing',
                      'scRNA-seq','scRNA-sequencing','scRNA sequencing','snRNA-seq',
                      'snRNA-sequencing','snRNA sequencing','sNuc-Seq','sNuc-Sequencing','sNuc Sequencing',
                      '(scRNA-seq)','[scRNA-seq]','(snRNA-seq)','[snRNA-seq]')
  return(RNAseqSearchItem)
}
#
# gse_search_item<-function(){
#   gseItem<-c("gse.title","gse.summary", "gse.overall_design")
#   return(gseItem)
# }

# Filter_organism<-function(keyword){
#   Fil_organism<-sprintf("(gsm.organism_ch1 IN (%s))",collapse_keywork(keyword))
#   return(Fil_organism)
# }
#
# Filter_sequencing<-function(){
#   Fil_seq<-"(gse.type LIKE '%high%throughput%sequencing%')"
#   return(Fil_seq)
# }
#
# Filter_SingleCell<-function(){
#   singleCellItem<-keyword_singlecell_search()
#   gseItem<-gse_search_item()
#
#   single_gse<-sapply(gseItem,function(y){
#     single_like<-sapply(singleCellItem,function(x){sprintf("(%s LIKE '%%%s%%')",y,x)})
#     single_collap<-paste(single_like,collapse=" OR ")
#   })
#
#   single_gse_add<-sapply(single_gse,function(z){sprintf("(%s)",z)})
#   out_single_gse0<-paste(single_gse_add,collapse=" OR ")
#   out_single_gse<-sprintf("(%s)",out_single_gse0)
#   return(out_single_gse)
# }
#
# Filter_RNAseq<-function(){
#   transcripItem<-keyword_RNAseq_search()
#   gseItem<-gse_search_item()
#
#   RNA_gse<-sapply(gseItem,function(y){
#     RNA_like<-sapply(transcripItem,function(x){sprintf("(%s LIKE '%%%s%%')",y,x)})
#     RNA_collap<-paste(RNA_like,collapse=" OR ")
#   })
#
#   RNA_gse_add<-sapply(RNA_gse,function(z){sprintf("(%s)",z)})
#   out_RNA_gse0<-paste(RNA_gse_add,collapse=" OR ")
#   out_RNA_gse<-sprintf("(%s)",out_RNA_gse0)
#   return(out_RNA_gse)
# }
# Filter_scRNAseq<-function(){
#   fSingle<-Filter_SingleCell()
#   fRNA<-Filter_RNAseq()
#   fExp<-Filter_sequencing()
#   f_SingleRNA<-sprintf("%s AND %s AND %s",fExp,fSingle,fRNA)
#   return(f_SingleRNA)
# }

keyword_ATACseq_search<-function(){
  ATACseqSearchItem<-c('assay of transposase accessible chromatin sequencing',
                       'ATAC','ATAC seq','ATAC-Seq','ATAC-sequencing','ATACseq','ATACsequencing',
                       'ATAC seq','ATAC sequencing','(scATAC-seq)','[scATAC-seq]',
                      'scATAC','scATAC-seq','scATAC-sequencing','scATAC sequencing',
                      'sciATAC-seq','sciATAC-Sequencing','(sciATAC-seq)','[sciATAC-seq]')
  return(ATACseqSearchItem)
}

#
# Filter_ATACseq<-function(){
#   ATAC_kw<-keyword_ATACseq_search()
#   gseItem<-gse_search_item()
#
#   ATAC_gse<-sapply(gseItem,function(y){
#     ATAC_like<-sapply(ATAC_kw,function(x){sprintf("(%s LIKE '%%%s%%')",y,x)})
#     ATAC_collap<-paste(ATAC_like,collapse=" OR ")
#   })
#
#   ATAC_gse_add<-sapply(ATAC_gse,function(z){sprintf("(%s)",z)})
#   out_ATAC_gse0<-paste(ATAC_gse_add,collapse=" OR ")
#   out_ATAC_gse<-sprintf("(%s)",out_ATAC_gse0)
#   return(out_ATAC_gse)
# }
#
# Filter_scATACseq<-function(){
#   fSingle<-Filter_SingleCell()
#   fATAC<-Filter_ATACseq()
#   fExp<-Filter_sequencing()
#   f_SingleATAC<-sprintf("%s AND %s AND %s",fExp,fSingle,fATAC)
#   return(f_SingleATAC)
# }
#
# GetPlatform<-function(Platform=Platform){
#   RNA_platform<-NULL
#   ATAC_platform<-NULL
#
#   PlatformNotIn<-Platform[!Platform %in% c("scRNAseq","scATACseq")]
#   if(length(PlatformNotIn)>0){
#     stop(sprintf("Only support scRNAseq and/or scATACseq. Your %s is wrong!",PlatformNotIn))
#   }
#
#   if("scRNAseq" %in% Platform){
#     RNA_platform<-Filter_scRNAseq()
#   }
#   if("scATACseq" %in% Platform){
#     ATAC_platform<-Filter_scATACseq()
#   }
#   com_platform<-c(RNA_platform,ATAC_platform)
#   com_platform<-com_platform[!is.null(com_platform)]
#   if(length(com_platform)>1){
#     com_platform_add<-sapply(com_platform,function(z){sprintf("(%s)",z)})
#     out_platform0<-paste(com_platform_add,collapse=" OR ")
#     out_platform<-sprintf("(%s)",out_platform0)
#   }else if(length(com_platform)==1){
#     out_platform<-com_platform
#   }else{
#     stop("No valid Platform is specified")
#   }
#   return(out_platform)
# }


check_KW_list<- function(user_kw=user_kw){
  if(!is.list(user_kw)){
    stop("Error: User keywords should be in a list. Try list(kw1=kw1,kw2=kw2).")
  }else{
    nkw<-length(user_kw)
    return(list(nkw=nkw,user_kw=user_kw))
  }
}
#
# Filter_userKW<-function(kwlist=kwlist){
#   checkKW<-check_KW_list(user_kw=kwlist)
#   nkw<-checkKW$nkw
#   user_keyword<-checkKW$user_kw
#
#   kw_name<-sprintf("user_kw%s",1:nkw)
#   for(i in 1:nkw){assign(sprintf("user_kw%s",i),user_keyword[[i]])} #eval(parse(text=sprintf("kw%s",i)))
#   gseItem<-gse_search_item()
#
#   user_wd_string<-sapply(kw_name,function(x1){
#     temp_x1<-eval(parse(text=x1))
#     temp_x2<-sapply(temp_x1,function(x2){
#       temp_x2_1<-sapply(gseItem,function(y){
#         temp_x2_2<-sapply(temp_x1,function(x3){sprintf("(%s LIKE '%%%s%%')",y,x3)})
#         temp_x2_3<-paste(temp_x2_2,collapse=" OR ")})
#     })
#     temp_x2_add<-sapply(temp_x2,function(z){sprintf("(%s)",z)})
#     out_temp_x20<-paste(temp_x2_add,collapse=" OR ")
#     out_temp_x2<-sprintf("(%s)",out_temp_x20)
#   })
#
#   user_wd_string_add<-sapply(user_wd_string,function(z){sprintf("(%s)",z)})
#   out_user_wd_string0<-paste(user_wd_string_add,collapse=" AND ")
#   out_user_wd_string<-sprintf("(%s)",out_user_wd_string0)
#   return(out_user_wd_string)
# }

select_OutGSM_item<-function(){
  gsm_item<-c("gsm.title", "gsm.gsm","gse_gsm.gse", "gsm.source_name_ch1", "gsm.characteristics_ch1",
              "gsm.treatment_protocol_ch1", "gsm.treatment_protocol_ch2")
  return(paste(gsm_item,collapse=","))
}

join_gse_gsm<-function(){
  join_gsegsm<-"gsm JOIN gse_gsm ON gsm.gsm=gse_gsm.gsm"
  return(join_gsegsm)
}

Obtain_GSM_fromGSE<-function(Gseid=Gseid,db=db){
  gseid_colla<-paste(sapply(unique(Gseid),function(x) sprintf("'%s'",x)),collapse=",")
  cmd_obtainGSM<-sprintf("SELECT %s FROM %s WHERE gse_gsm.gse IN (%s)",select_OutGSM_item(),join_gse_gsm(),gseid_colla)
  GsmOut<-RSQLite::dbGetQuery(db,cmd_obtainGSM)
  return(GsmOut)
}

downloadSupplFile<-function(downGSEList=downGSEList){
  sapply(downGSEList,function(x){
    GEOquery::getGEOSuppFiles(x)
    })
}

isEmptyFile<-function(inputString){
  isEmptyF<-ifelse(length(inputString)==0,TRUE,FALSE)
  return(isEmptyF)
}

extractZip<-function(ExtractID=ExtractID){
  sapply(ExtractID,function(x){
    allzipfile<-list.files(path=paste(x),pattern=".zip")
    if(!isEmptyFile(allzipfile)){
      sapply(allzipfile,function(y){unzip(zipfile=sprintf("%s/%s",x,y),exdir=sprintf("%s",x))})}
  })
}

extractGZ<-function(ExtractID=ExtractID){
  sapply(ExtractID,function(x){
    allGzfile<-list.files(path=paste(x),pattern=".gz")
    if(!isEmptyFile(allGzfile)){
      allTarGzfile<-list.files(path=paste(x),pattern="tar.gz")
      if(isEmptyFile(allTarGzfile)){
        sapply(allGzfile,function(y){GEOquery::gunzip(filename=sprintf("%s/%s",x,y),remove=FALSE,overwrite =TRUE)})
        }else{
        allTarGzfile<-list.files(path=paste(x),pattern="tar.gz")
        allGzfile_unique<-allGzfile[!allGzfile %in% allTarGzfile]
        if(!isEmptyFile(allGzfile_unique)){
          sapply(allGzfile_unique,function(y){GEOquery::gunzip(filename=sprintf("%s/%s",x,y),remove=FALSE,overwrite =TRUE)})
        }}
      }
    })
}

extractTar<-function(ExtractID=ExtractID){ ##can also handle .tar.gz
  sapply(ExtractID,function(x){
    allTarfile<-list.files(path=paste(x),pattern=".tar")
    if(!isEmptyFile(allTarfile)){
      sapply(allTarfile,function(y){untar(tarfile=sprintf("%s/%s",x,y),exdir=sprintf("%s",x))})}
  })
}

extractAll<-function(ExtractID=ExtractID){
  extractTar(ExtractID=ExtractID)
  extractGZ(ExtractID=ExtractID)
  extractZip(ExtractID=ExtractID)
}


CheckRowName<-function(DataInput=DataInput){
  rowN0<-row.names(DataInput)
  rowN1<-gsub("-",".",rowN0)
  rowN<-gsub("_",".",rowN1)
  row.names(DataInput)<-rowN
  return(DataInput)
}

CheckFolderExist<-function(checkList=checkList){
  currentFile<-list.files()
  NewCheckList<-checkList
  if(sum(checkList %in% currentFile)>0){
    NewCheckList<-checkList[!checkList %in% currentFile]
    InCheckList<-checkList[checkList %in% currentFile]
  #  print(sprintf("%s already exist!",InCheckList))
  }
  return(NewCheckList)
}

removeFolder<-function(ToBeRemove=ToBeRemove){
  sapply(ToBeRemove,function(x){unlink(x, recursive = TRUE)})
}

Seurat_normalization<-function(Seu_obj_input=seu_obj_input,
                               Norm.method=Norm.method,
                               Scale.factor=Scale.factor,
                               Feature.selection.method=Feature.selection.method,
                               Nfeatures = Nfeatures,Do.scale = Do.scale,Do.center = Do.center){

  if(Norm.method %in% c("lognorm","Lognorm","LogNorm","LogNormalize")){
    seu_obj_input <- Seurat::NormalizeData(Seu_obj_input, normalization.method ="LogNormalize", scale.factor = Scale.factor,verbose = FALSE)
    seu_obj_input <- Seurat::FindVariableFeatures(seu_obj_input, selection.method = Feature.selection.method, nfeatures = Nfeatures,verbose = FALSE)
    seu_obj_input <- Seurat::ScaleData(seu_obj_input,do.scale = Do.scale,do.center = Do.center,verbose = FALSE)
  }else if(Norm.method %in% c("sct","sctransform","SCT","Sctransform","Sct")){
    seu_obj_input <- Seurat::SCTransform(Seu_obj_input, method = "glmGamPoi",verbose = FALSE,variable.features.n =Nfeatures,do.scale = Do.scale,do.center = Do.center)
  }else{
    stop("Norm.method should be 'lognorm' or 'sct'")
  }
  return(seu_obj_input)
}


Seurat_dimensionReduc<-function(Seu_obj_input=seu_obj_input,npcs =npcs,Dims = dims,Resolution=Resolution,Algorithm=Algorithm){
  seu_obj_input <- Seurat::RunPCA(Seu_obj_input,npcs =npcs,verbose = FALSE)
  seu_obj_input <- Seurat::FindNeighbors(seu_obj_input, dims = Dims,verbose = FALSE)
  seu_obj_input <- Seurat::FindClusters(seu_obj_input,resolution=Resolution,algorithm = Algorithm,verbose = FALSE)
  seu_obj_input <- Seurat::RunUMAP(seu_obj_input, dims = Dims,verbose = FALSE)
  return(seu_obj_input)
}

CheckNormMethodForNormData<-function(Seu_obj_input=seu_obj_input,
                                     Feature.selection.method=Feature.selection.method,
                                     Nfeatures = Nfeatures,Dims = Dims,Do.scale = Do.scale,Do.center = Do.center){
    seu_obj_input <- Seurat::FindVariableFeatures(Seu_obj_input, selection.method = Feature.selection.method, nfeatures = Nfeatures,verbose = FALSE)
    seu_obj_input <- Seurat::ScaleData(seu_obj_input,do.scale = Do.scale,do.center = Do.center,verbose = FALSE)
  return(seu_obj_input)

}


CheckUniqueClusterNumber<-function(SeuObj=SeuObj){
  ClusterCount<-length(unique(SeuObj@meta.data$seurat_clusters))
  return(ClusterCount)
}



SelectObjWithTopNumCluster<-function(ObjList=objList,TopN=topN){
  if(length(ObjList)<TopN){TopN<-length(ObjList)}
  countUniqNClust<-lapply(ObjList,function(x) CheckUniqueClusterNumber(x))
  topN_Name<-names(ObjList[order(-unlist(countUniqNClust))][1:TopN])
  return(topN_Name)
}

checkObj2_forInt<-function(Obj1=Obj1,
                           Obj2=Obj2,
                           RefSamN=referenceSampleN,
                           NumOfSamForRef=NumberOfSampleForReference,
                           NumOfSamForRef2=NumberOfSampleForReference2){
  Obj.all<-Obj1
  RefSamN.all<-RefSamN
  if(!is.null(Obj2)){
    if(sum(names(Obj2) %in% names(Obj1))>0){
      stop("Object.list2 has the item name as Object.list1! Change the item name.")}
    RefSamN2<-SelectObjWithTopNumCluster(ObjList=Obj2,TopN=NumOfSamForRef)
    if(!is.null(NumOfSamForRef2)){
      RefSamN2<-SelectObjWithTopNumCluster(ObjList=Obj2,TopN=NumOfSamForRef2)
    }
    RefSamN.all<-c(RefSamN,RefSamN2)
    Obj.all<-c(Obj1,Obj2)
  }
  Out<-list(Obj.all=Obj.all,RefSamN.all=RefSamN.all)
  return(Out)
}

checkRefIndex<-function(obj.all=object.list.all,
                        refSamN=referenceSampleN){
  refIndex<-which(names(obj.all) %in% refSamN)
  return(refIndex)
}


