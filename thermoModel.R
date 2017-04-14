#thermoModel.R

# Sig2GRN 
# The simulation generates the dynamical trends of the signaling proteins using a generalized thermodynamic model
# Thermodynamic model is derived under the assumption that the system is in the thermodynamic equilibrium. 
# Therefore the gene expression level is defined as a function of the activity levels of the bound transcription factors

thermoModel <- function(thermoParametersMap){
  #thermoParametersMap: R hash environment which the keys are the names of the transcription factor(node) and the values are time-series data of the transcription factor activities (value ∈[0,1] ] at each simulation iteration) 
  #global variable FinalDataGraph: R hash environment which consists of time-series data of the transcription factor activities 
  #edges: list of edges that directed graph G as an input of Sig2GRN 
  #resultsPos: list of positive combinations from G
  # resultAll: list of all combinations from G
  # ==== SETUP All input nodes from thermoParametersMap ======================================================  
  
  FinalDataGraph <<- new.env(hash=T)
  edges <- vector(mode="character")
  resultsPos <-list()
  resultsAll <-list()
  numOfedges = length(edges)
  inputNodes = ls(thermoParametersMap)
  numberOfTime = length(FinalDataGraph[[inputNode]])
  thermoFinalResult<- vector(mode="double", length = numberOfTime)

  
  for(edge in edges){
    edgeSplit = strsplit(edge, "[-]")
    source = edgeSplit[1]
    target = edgeSplit[2]
  }
  activators <- vector(mode="character")
  inhibitors <- vector(mode="character")
    
  for(i in 1:length(inputNodes)){
    checkEdge <- paste(inputNodes[i], "target", sep="-")
    if(checkEdge %in% edges){
      activators[i] <- inputNodes[i]
    }
    else{
      inhibitors[i] <- inputNodes[i]
    }
  }
  
  if(length(activators) == 0){
    for(i in 1:numberOfTime){
      thermoFinalResult[i] = 0
    }
  }
  else{
    for (time in numberOfTime){
      totalPositive = 0
      totalCombination = 0
      
      #all Positive
      for(i in 1:length(resultsPos)){
      array = list()
      j = 1
      for(inputNode in resultsPos){
        activityLevel = FinalDataGraph[[inputNode]][time]
        concentration = activityLevel * thermoParametersMap[[inputNode]]
        array[j] <- concentration
        j = j + 1
      }
      totalConcentration = 1.0
      for(num in array){
        if(num != 0){
          totalConcentration = totalConcentration * num
        }
      }
        
      totalPos = totalPos + totalConcentration
    }  
      
    #all combination
    for(i in 1:length(resultsAll)){  
      array = list()
      j = 1
      for(inputNode in resultsAll){
        activityLevel = FinalDataGraph[[inputNode]][time]
        concentration = activityLevel * thermoParametersMap[[inputNode]]
        array[j] <- concentration
        j = j + 1
      }
      totalConcentration = 1.0
      for(num in array){
        if(num != 0){
          totalConcentration = totalConcentration * num
        }
      }
      
      totalCombination = totalCombination + totalConcentration
    }
      
    thermoFinalResult[i] = totalPos / (totalConcentration + 1)
    }
  }
    return (thermoFinalResult)

}