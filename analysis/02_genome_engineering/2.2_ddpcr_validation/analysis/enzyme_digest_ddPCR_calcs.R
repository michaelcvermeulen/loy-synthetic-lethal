# Function to calculate DNA and water volumes
calculate_volumes <- function(C_sample,sample_name) {
  # Constants
  total_volume <- 20          # Total reaction volume (in µL)
  RE_volume <- 2              # Restriction Enzyme volume (in µL)
  buffer_volume <- 2          # Buffer volume (in µL)
  DNA_required <- 1000        # DNA required (in ng)
  
  # Calculate V_DNA (µL)
  V_DNA <- DNA_required / C_sample
  
  # Calculate V_water (µL)
  V_water <- total_volume - (V_DNA + RE_volume + buffer_volume)
  
  # If more DNA is needed, double the reaction
  if(V_water < 0){
    
    message(C_sample,"ng/uL is too low for a 20uL reaction - DOUBLING to 40uL")
    total_volume <- total_volume * 2
    RE_volume <-  RE_volume * 2              
    buffer_volume <- buffer_volume * 2 
    
  }
  
  V_water <- total_volume - (V_DNA + RE_volume + buffer_volume)
  
  # 10 µL of 2X ddPCR Supermix for Probes.
  # 1 µL of 20X TaqMan Probe Assay.
  # Template DNA/RNA (amount varies based on experiment).  50ng
  # Nuclease-free water to bring the total volume to 20 µL.
  
  (DNA_required / total_volume) -> DNA_conc
  
  
  
  
  # Return the results
  data.frame(
       Sample = sample_name,
       DNA_concentration = paste0(C_sample,"ng/uL"),
       Total_Reaction_volume = paste0(total_volume ,"uL"),
       RE_volume = paste0(RE_volume,"uL"),
       buffer_volume = paste0(buffer_volume,"uL"),
       DNA_sample_volume = paste0(round(V_DNA,3),"uL"), 
       Water_volume = paste0(round(V_water,3),"uL"),
       DNA_conc = paste0(DNA_conc, "ng/uL")) -> a 
  
  DNA_conc <- as.numeric(stringr::str_split(string = a$DNA_conc, pattern = "ng/uL", simplify = T)[,1])
  Reaction_volume <- 20
  Taqman_Primers <- 1
  ddPCR_master_mix <- 10
  ddPCR_DNA_Volume <- 50 / DNA_conc
  ddPCR_Water_Volume <- 20 - (ddPCR_DNA_Volume + ddPCR_master_mix + Taqman_Primers)
  
  data.frame(
    Sample = sample_name,
    DNA_concentration = paste0(C_sample,"ng/uL"),
    Total_Reaction_volume = paste0(Reaction_volume ,"uL"),
    ddPCR_master_mix = paste0(ddPCR_master_mix,"uL"),
    Taqman_Primers = paste0(Taqman_Primers,"uL"),
    ddPCR_DNA_Volume = paste0(ddPCR_DNA_Volume,"uL"), 
    Water_volume = paste0(ddPCR_Water_Volume,"uL"),
    Total_DNA = paste0(ddPCR_DNA_Volume * DNA_conc, "ng")) -> b 
  
  return(list(a,b))
  
}

# Example: Assume sample concentration is 200 ng/µL
C_sample <- 502.9
results <- calculate_volumes(C_sample, "")


## input the filename of the input sheet in the same format as the example sheet
## output_dir is where the output files go
## experiment name is appended to the outputfile
make_ddPCR_experiment_table <- function(input_sheet = "A:/Dropbox/LOY_paralogs/code/Verification/example_input_sheet.csv",
                                        output_dir = "A:/Dropbox/LOY_paralogs/code/Verification/", experiment_name){
 
  data.table::fread(input_sheet) -> input_sheet
  lapply(X = 1:nrow(input_sheet), FUN = function(z){
    
    input_sheet[z,] -> q
    
    results <- calculate_volumes(q$Conc, q$Sample)
    return(results)
    
  }) -> l
  lapply(l, function(x) x[[1]]) %>% data.table::rbindlist() -> digest
  lapply(l, function(x) x[[2]]) %>% data.table::rbindlist() -> ddPCR
  
  
  data.table::fwrite(file = paste0(output_dir,"/",experiment_name,"_restriction_enzyme.csv"), x = digest )
  data.table::fwrite(file = paste0(output_dir,"/",experiment_name,"_ddCPR.csv"), x = ddPCR )
  
  return(list(digest,ddPCR))
  
}


make_ddPCR_experiment_table(input_sheet = "A:/Dropbox/LOY_paralogs/code/Verification/NOV28.csv",
                            experiment_name = "NOV28")


