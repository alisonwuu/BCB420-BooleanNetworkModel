#booleanModel.R

# Sig2GRN 
# The simulation generates the dynamical trends of the signaling proteins using a generalized logical model. 
# The goal of the upstream simulation is to generate the dynamical trends of the transcription factors under the perturbations. After generating the time-series data of the transcription factor activities (value ∈[0,1] at each simulation iteration), the simulated transcription factor activities are employed to predict the gene expression patterns over time using the Boolean Model.
# The AND logical relation is assigned to the transcription factors that have the same transcriptional regulation type (e.g., activation or inhibition) for a gene, so that the gene will be switched ON (or OFF) when the maximum activity level of activating (or inhibiting) transcription factors surpasses a user-defined threshold (value ∈(0,1)). 
# When both activation and inhibition regulations are present on the same gene, the inhibition is assumed to precede the activation. The simulation result of the Boolean model is a list of 0s and 1s, over the course of time.

booleanModel <- function(inhibitionActivationNodes){
  #inhibitionActivationNodes: R hash environment which the keys are the names of the transcription factor(node) and the values are time-series data of the transcription factor activities (value ∈[0,1] ] at each simulation iteration) 
  #global variable FinalDataGraph: R hash environment which consists of time-series data of the transcription factor activities 

  #==== SETUP All input nodes from inhibitionActionNodes ======================================================  
  # Check if all input nodes are activation

  FinalDataGraph <<- new.env(hash=T)
  inputNodes = ls(inhibitionActivationNodes)
  
  allActivation = TRUE
  for(inputNode in inputNodes){
    allActivation = inhibitionActivationNodes[[inputNode]]
  }
  
  # Create a hash that holds the name of all node and their time series data 
  allRequiredInputNodes <- new.env(hash=T)
  if(allActivation){
    for(inputNode in inputNodes){
      assign(inputNode, FinalDataGraph[[inputNode]], allRequiredInputNodes)
    }
  }
  else{
    for(inputName in inputNodes){
      if(!inhibitionActivationNodes[[inputNode]]){
        assign(inputNode, FinalDataGraph[[inputNode]], allRequiredInputNodes)
      }
    }
  }
  
  # Check the threshold of each time series data of each gene
  # For all activation, if greater than user-defined threshold, the expression level will be 1, otherwise 0
  # For all inhibition, sum up the total threshold of all nodes of each time series. if the sum is greater than the multiple of the user-defined threshold and the number of nodes, then the expression level will be 0, otherwise 1.
  numberOfTime = length(FinalDataGraph[[inputNode]])
  threshold = 0.3
  expressionResult<- vector(mode="double", length = numberOfTime)
  if(allActivation){
    for(time in 1:numberOfTime){
      allGreaterThanThreshold = FALSE
      allRequiredInputNodeNames = ls(allRequiredInputNodes)
      for(inputNode in allRequiredInputNodeNames){
        allGreaterThanThreshold = allRequiredInputNodes[[allRequiredInputNodeNames]][time] > threshold
        if(allGreaterThanThreshold){
          break
        }
      }
    }
   if(allGreaterThanThreshold){
     expressionResult[time] = 1
   }
   else{
     expressionResult[time] = 0
   }
  }else{
    for(time in 1:numberOfTime){
      allRequiredInputNodeNames = ls(allRequiredInputNodes)
      sumOfValues = 0.0
      for(inputNode in allRequiredInputNodeNames){
        sumOfValues = sumOfValues + allRequiredInputNodes[[allRequiredInputNodeNames]][time]
      }
      totalNodes = length(allRequiredInputNodeNames)
      if(sumOfValues >= threshold * totalNodes){
        expressionResult[time] = 0
      }
      else{
        expressionResult[time] = 1
      }
    }
  }
  
  return (expressionResult)
  

  
}