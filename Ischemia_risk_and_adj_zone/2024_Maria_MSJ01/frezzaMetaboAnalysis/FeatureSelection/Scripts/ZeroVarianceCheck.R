

ZeroVarCheck <- function(Result_Folder,
                         scriptpath = scriptpath,
                         OutputName = OutputName,
                         Input_data, 
                         RemoveFeatures = TRUE){
  
# Check metabolite variance
metabolite_var <-  as.data.frame(apply(Input_data, 2, function(x) var(x, na.rm = TRUE)) |> t()) # calculate each metabolites variance
metabolite_zero_var_list <- colnames(metabolite_var)[which(metabolite_var[1,]==0)]# takes the names of metabolites with zero variance and puts them in list
  
# Remove the zero variance metabolites
if(RemoveFeatures==TRUE){
  # Print a warning if Zero var metabolites were identified
  if(length(colnames(metabolite_var)[which(metabolite_var[1,]==0)]) > 0 ){
    message(paste0(length(metabolite_zero_var_list)), " Metabolites with zero variance have been identified:", paste0(metabolite_zero_var_list, by=", "))
  }
  Input_data_filtered <- Input_data |> select(-all_of(metabolite_zero_var_list))
}
# Save resulting table
write.table(metabolite_zero_var_list, row.names = FALSE,col.names = FALSE, file =  paste(Result_Folder,"/Zero_variance_metabolites", OutputName,".csv",sep =  "")) # save zero var metabolite list
  
return(list("Input_data_filtered"=Input_data_filtered,"ZeroVarMetabolites" = metabolite_zero_var_list) )
}
