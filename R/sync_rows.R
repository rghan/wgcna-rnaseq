sync.rows <- function(dataset.1,dataset.2){
  # function that takes in two datasets and synchronizes
  # their rows for the same name, so that at the end both
  # datasets contain the same rownames
  # prerequisite: a rowname can be present only once for each dataset
  
  # create output variables
  output.dataset.1 <- matrix(ncol=ncol(dataset.1),nrow=0)
  output.dataset.2 <- matrix(ncol=ncol(dataset.2),nrow=0)
  
  # set the colnames to input parameters
  colnames(output.dataset.1) <- colnames(dataset.1) 
  colnames(output.dataset.2) <- colnames(dataset.2) 
  
  # now go over all rows of dataset 1 and compare to dataset 2
  for(i in 1:nrow(dataset.1)){
    # read entries
    x <- which(tolower(row.names(dataset.2))==tolower(row.names(dataset.1))[i])
    
    # see if entry exists
    if(length(x)>0) { # it exists in second dataset
      output.dataset.1 <- rbind(output.dataset.1,dataset.1[i,])
      output.dataset.2 <- rbind(output.dataset.2,dataset.2[x,])
    } # end if
  } # end for loop
  
  # return variable
  y <- cbind(output.dataset.1,output.dataset.2)
  y
} # end function

sync.cols <- function(dataset.1,dataset.2){
  # function that takes in two datasets and synchronizes
  # their cols for the same name, so that at the end both
  # datasets contain the same rownames
  # prerequisite: a colname can be present only once for each dataset
  
  # create output variables
  output.dataset.1 <- matrix(nrow=nrow(dataset.1),ncol=0)
  output.dataset.2 <- matrix(nrow=nrow(dataset.2),ncol=0)
  
  # set the rownames to input parameters
  rownames(output.dataset.1) <- rownames(dataset.1) 
  rownames(output.dataset.2) <- rownames(dataset.2) 
  
  # now go over all cols of dataset 1 and compare to dataset 2
  for(i in 1:ncol(dataset.1)){
    # read entries
    x <- which(tolwer(colnames(dataset.2))==tolower(colnames(dataset.1))[i])
    
    # see if entry exists
    if(length(x)>0) { # it exists in second dataset
      output.dataset.1 <- cbind(output.dataset.1,dataset.1[,i])
      colnames(output.dataset.1)[ncol(output.dataset.1)] <- colnames(dataset.1)[i]
      output.dataset.2 <- cbind(output.dataset.2,dataset.2[,x])
    } # end if
  } # end for loop
  
  # return variable
  y <- rbind(output.dataset.1,output.dataset.2)
  y
} # end function