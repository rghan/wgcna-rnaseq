divide.genotypes.color <- function(data, accs){
  # function dividing the datasets accoring to the groups in accs
  # we have four potential groups
  ret.list <- vector("list",2)
  
  group.5.white <- which(accs[,2]==5 & tolower(accs[,3])=="w"); 
  group.5.red <- which(accs[,2]==5 & tolower(accs[,3])=="r"); 
  
  rows.5.w <- match(tolower(accs[group.5.white,1]),tolower(rownames(data))); #rows.5 <- rows.5[-(which(is.na(rows.5)))]
  if(length(rows.5.w>0)){ 
    temp <- vector("integer",0)
    for(i in 1:length(rows.5.w)) if(is.na(rows.5.w[i])) temp <- c(temp,i)
    if(length(temp)>0) rows.5.w <- rows.5.w[-temp]
  }
  rows.5.r <- match(tolower(accs[group.5.red,1]),tolower(rownames(data))); #rows.4 <- rows.4[-(which(is.na(rows.4)))]
  if(length(rows.5.r>0)){
    temp <- vector("integer",0)
    for(i in 1:length(rows.5.r)) if(is.na(rows.5.r[i])) temp <- c(temp,i)
    if(length(temp)>0) rows.5.r <- rows.5.r[-temp]
  }
  ret.list[[1]] <- data[rows.5.w,]
  ret.list[[2]] <- data[rows.5.r,]
  return(ret.list)
} # end function

divide.genotypes <- function(data, accs){
  # function dividing the datasets accoring to the groups in accs
  # we have four potential groups
  ret.list <- vector("list",4)
  
  group.5 <- which(accs[,2]==5); 
  group.4 <- which(accs[,2]==4); 
  group.3 <- which(accs[,2]==3);
  group.1 <- which(accs[,2]==1); 
  
  rows.5 <- match(tolower(accs[group.5,1]),tolower(rownames(data))); #rows.5 <- rows.5[-(which(is.na(rows.5)))]
  if(length(rows.5>0)){ 
    temp <- vector("integer",0)
    for(i in 1:length(rows.5)) if(is.na(rows.5[i])) temp <- c(temp,i)
    if(length(temp)>0) rows.5 <- rows.5[-temp]
  }
  rows.4 <- match(tolower(accs[group.4,1]),tolower(rownames(data))); #rows.4 <- rows.4[-(which(is.na(rows.4)))]
  if(length(rows.4>0)){
    temp <- vector("integer",0)
    for(i in 1:length(rows.4)) if(is.na(rows.4[i])) temp <- c(temp,i)
    if(length(temp)>0) rows.4 <- rows.4[-temp]
  }
  rows.3 <- match(tolower(accs[group.3,1]),tolower(rownames(data))); #rows.3 <- rows.3[-(which(is.na(rows.3)))]
  if(length(rows.3>0)){
    temp <- vector("integer",0)
    for(i in 1:length(rows.3)) if(is.na(rows.3[i])) temp <- c(temp,i)
    if(length(temp)>0) rows.3 <- rows.3[-temp]
  }
  rows.1 <- match(tolower(accs[group.1,1]),tolower(rownames(data))); #rows.1 <- rows.1[-(which(is.na(rows.1)))]
  if(length(rows.1>0)){
    temp <- vector("integer",0)
    for(i in 1:length(rows.1)) if(is.na(rows.1[i])) temp <- c(temp,i)
    if(length(temp)>0) rows.1 <- rows.1[-temp]
  }
  ret.list[[1]] <- data[rows.1,]
  ret.list[[2]] <- data[rows.3,]
  ret.list[[3]] <- data[rows.4,]
  ret.list[[4]] <- data[rows.5,]
  return(ret.list)
} # end function

network.difference <- function(sig.cor.ls){
  # function determine the network difference for each network submitted to other networks submitted
  # with parameter of type list
  # prepare a return variable of type list
  ret.list <- sig.cor.ls
  
  for(i in 1:length(sig.cor.ls)){
    # determine sequence for network comparison within list vector
    seq <- NULL
    if(i==1)  seq  <- c(2:length(sig.cor.ls))
    else if(i==length(sig.cor.ls)){
      if(i==2) seq <- 1
      else seq <- c(1:(length(sig.cor.ls)-1))
    }
    else seq <- c(1:(i-1),(i+1):length(sig.cor.ls))  
    
    links <- which(sig.cor.ls[[i]]!=0,arr.ind=T)
    
    # now make the comparison
    for(j in 1:nrow(links)){
      row.n <- tolower(rownames(sig.cor.ls[[i]]))[links[j,1]]
      col.n <- tolower(colnames(sig.cor.ls[[i]]))[links[j,2]]
      for(k in seq){
        row.num <- which(tolower(rownames(sig.cor.ls[[k]]))==row.n)
        col.num <- which(tolower(colnames(sig.cor.ls[[k]]))==col.n)
        # if exists then set corresponding cell value of ret.list to zero
        if(length(row.num)>=1 & length(col.num)>=1){
          # check for the value
          if(sig.cor.ls[[k]][row.num,col.num]!=0){
            # it does exist
            ret.list[[i]][links[j,1],links[j,2]] <- 0
            break; # no need to look into the other networks
          } # end if
        } # end if
      } # end nested nested loop
    } # end nested loop
  } # end for loop
  return(ret.list)
} # end function
################################################################################################################################################################
network.intersect <- function(sig.cor.ls){
  # function determining the network intersect between all supplied networks in paratmeter sig.cor.ls of type list
  
  # first convert all matrices omittiing the actual r values
  for(i in 1:length(sig.cor.ls)){
    sig.cor.ls[[i]][sig.cor.ls[[i]]<0] <- 2 # for negative correlation
    sig.cor.ls[[i]][sig.cor.ls[[i]]>0] <- 1 # for positive correlation
  } # end for loop
  
  # find all connections in the first matrix and use it as originator
  ret.mat <- sig.cor.ls[[1]]
  links <- which(ret.mat!=0,arr.ind=T)
  
  # now look if these links are also present in the other matrices
  # !!!! very important the nodes must have same names
  
  for(i in 1:nrow(links)){
    row.n <- tolower(rownames(ret.mat)[links[i,1]])
    col.n <- tolower(colnames(ret.mat)[links[i,2]])
    
    for(j in 2:length(sig.cor.ls)){
      row.num <- which(tolower(rownames(sig.cor.ls[[j]]))==row.n)
      col.num <- which(tolower(colnames(sig.cor.ls[[j]]))==col.n)
      # check if the connection is existent in the other matrix
      if(length(row.num)>=1 & length(col.num)>=1 ){
        # nodes exist
        # now check if the connection between nodes exist
        con <- sig.cor.ls[[j]][row.num,col.num]
        if(con==0){
          # it doesn't exist so set connection to 0 in originator matrix and break the loop
          ret.mat[links[i,1],links[i,2]] <- 0
          break;
        }
        else{
          if(con==1)
            ret.mat[links[i,1],links[i,2]] <- (ret.mat[links[i,1],links[i,2]])+1
          else
            ret.mat[links[i,1],links[i,2]] <- (ret.mat[links[i,1],links[i,2]])+2
        }
      } # end if
      else{
        # connection does not exist
        # set connection to 0
        ret.mat[links[i,1],links[i,2]] <- 0
        break;
      } # end else
    } # end nested loop
  } # end for loop
  
  # some explanation for the coming code:
  # the dataset returned contains only those connection that are present in all networks, i.e. positive and negative
  # the positive connections are denoted by value 1 and the negative by value 2
  # if we add up the values and the same connection was positive in all networks, we will get a value of 
  # datasets provided. Therefore, if value/number of datasets = 1 - only under this scenario we have the same positive connections
  # if value/number of datasets = 2 then we have all negative correlation in all datasets
  # any value inbetween 1 and 2 indicates that the connetion exist in all datasets but they are not the same
  # positive and negtive
  # thus
  ret.mat[ret.mat != 0] <- (ret.mat[ret.mat != 0])/(length(sig.cor.ls))
  # set diagonal to zero - just to make sure
  for(i in nrow(ret.mat)) ret.mat[i,i] <- 0
  
  return(ret.mat)
  
} # end function
################################################################################################################################################################
compare.igraph.communites <- function(igr.list){
  # function comparing the communities of all graphs provided in the parameter igr.list of type list containing igraph objects
  # hardcoded is the walktrap.community to create communities and compare them
  # the return variable is a symmetric matrix
  library(igraph)
  ret.mat <- matrix(nrow=length(igr.list),ncol=length(igr.list),data=0)
  rownames(ret.mat) <- c(1:length(igr.list))
  colnames(ret.mat) <- c(1:length(igr.list))
  for(i in 1:length(igr.list)){
    colnames(ret.mat)[i] <- igr.list[[i]]$name
    rownames(ret.mat)[i] <- colnames(ret.mat)[i]
  } # end for loop
  # fill in value
  for(i in 1:(length(igr.list)-1)){
    for(j in i:length(igr.list)){
      comp <- compare(membership(walktrap.community(igr.list[[i]])),
                      membership(walktrap.community(igr.list[[j]])),method ="rand")
      ret.mat[i,j] <- comp; ret.mat[j,i] <- comp
    } # end nested loop
  } # end for loop
  # return variable
  return(ret.mat)
} # end function
################################################################################################################################################################

get.igraph <- function(data){
  # function taking in the correlation adjacency matrix and returning the network as an igraph object for subsequent analysis
  # load library
  library(igraph)
  sig.netw.zero <- data
  sig.netw.zero[which(sig.netw.zero!=0)]  <- 1
  # convert diagonal to 0's
  for(i in 1:nrow(sig.netw.zero)){sig.netw.zero[i,i] <- 0 }
  # make network graph
  gr <- graph.adjacency(sig.netw.zero,mode="lower",diag=F,add.colnames=T)
  gr$names <- row.names(sig.netw.zero)
  # find nonconnected nodes
  zero <- which(degree(gr)==0)
  # delete nonconnected nodes from graph
  gr <- delete.vertices(gr,zero)
  if(length(zero)!=0) gr$names <- gr$names[-zero]
  return(gr) # undirected
} # end function
################################################################################################################################################################

normalzation_pietro <- function(data){
  # making sure data is in matrix format and numeric
  data <- as.matrix(data)
  row.names(data) <- data[,1]
  data <- data[,-1]
  data <- convert_to_num(data)
  
  data <- remove_na_col(data)
  data <- remove_na_row(data)
  data <- remove_zero_col(data)
  data <- remove_zero_row(data)
  library(pcaMethods)
  data <- pca(data,method="ppca")@completeObs
  # finding all rows of blank and standard mix and deleting them
  data <- data[-(which(tolower(row.names(data))=="blank")),]
  data <- data[-(which(tolower(row.names(data))=="stdmix")),]
  
  # finding all QCs and storing them
  QC <- data[which(tolower(row.names(data))=="qc"),]
  data <- data[-(which(tolower(row.names(data))=="qc")),]
  avg.QC <- apply(QC,2,mean,trim=0.05)
  
  # do the normalization by the mean of the QC
  for(i in 1:ncol(data)) data[,i] <- data[,i]/avg.QC[i]
  
  data[data=="Inf"] <- NA
  
  data <- remove_na_col(data)
  data <- remove_na_row(data)
  data <- remove_zero_col(data)
  data <- remove_zero_row(data)
  # return
  return(data)
} # end function


################################################################################################################################################################

determine.best.clustering <- function(data){
  # Function using the pam - kmedeoid clustering iteratively using the sliouhette widhts to determine the 
  # best number of groups in a given dataset= data
  # load library
  library(cluster)
  
  sil <- numeric(30) # vector holding the sillouhette widths
  for(k in 2:30){ # lets suppose we have a maximum of 30 clusters
    sil[k] <- pam(data,k)$silinfo$avg.width
  } # end for loop
  k.best <- which.max(sil)
  cat("silhouette-optimal number of clusters:", k.best, "\n")
  
  # print graphical output
  plot(1:30, sil, type= "l", main = "pam() clustering assessment",
       xlab= "k  (# clusters)", ylab = "average silhouette width")
  axis(1, k.best, paste("best",k.best,sep="\n"), col = "red", col.axis = "red")
  
  return(k.best)
  
}# end function
################################################################################################################################################################
sort.cluster <- function(data,k){
  # function returns a dataset = data into number of clusters specified = k
  # to partition the data the pam function of library cluster is used
  
  # load library 
  library(cluster)
  pam <- pam(data,k)
  groups <- pam$clustering
  # do grouping !!! by rows
  ret.var <- vector("list",k)
  for(i in 1:k){
    ret.var[[i]] <- data[which(groups==i),]
  } # end for loop
  return(ret.var)
  
} # end function
################################################################################################################################################################
# another time
# prepare.clust.list <- function(groups){
#   # groups of type list
#   rows <- nrow(groups[[1]])
#   if(length(groups)>1){
#     for(i in 2:length(groups)){
#       if(nrow(groups[[i]])>rows) rows <- nrow(groups[[i]])
#     }
#   }
#     
#   ret.val <- matrix(nrow=rows,ncol=length(groups))
#   for(i in 1:length(groups)){
#     ret.val[1:nrow(groups[[i]]),i] <- rownames(groups[[i]])
#     colnames(ret.val)[i] <- paste("group",i)
#   }
#   return(ret.val)
# }
################################################################################################################################################################
triangle <- function(data){
  # get only the upper triangle of a symmetric matrix and store in a vector
  ret.vec <- NULL
  first.col <- 2
  
  for(k in 1:nrow(data)){
    for(i in first.col:ncol(data)){
      ret.vec <- rbind(ret.vec,data[k,i])
    }# end nested loop
    first.col <- first.col+1
    if(first.col==ncol(data)){break}
  } # end for loop
  
  return(ret.vec)  
} # end function

################################################################################################################################################################
# FDR function showing p values, the frequency, and the actual number of false positives
# within the significantly identified instances

fdr.fun <- function(p.value){
  # load library
  library(qvalue)
  q.vec <- qvalue(p.value) # determining qvalue
  # preparing return matrix
  fdr.matrix <- matrix(ncol=5,nrow=0)
  colnames(fdr.matrix) <- c("pvalue","qvalue","frequency","False positives","%")
  
  seq <- seq(0.001,0.05,0.001)  # for fdr
  
  for(i in seq){
    fdr.matrix <- rbind(fdr.matrix,c(i,max(q.vec$qvalues[q.vec$pvalues <= i]),length(which(q.vec$pvalues<=i)),0))
    fdr.matrix[nrow(fdr.matrix),4] <- fdr.matrix[nrow(fdr.matrix),2]*fdr.matrix[nrow(fdr.matrix),3]
    fdr.matrix[nrow(fdr.matrix),5] <- fdr.matrix[nrow(fdr.matrix),4]/fdr.matrix[nrow(fdr.matrix),3]
  }
  
  return(fdr.matrix)
  
} # end function
################################################################################################################################################################
# removing all columns with NAs only

remove_zero_col <- function(matrix){
  x  <- NULL
  for(i in 1:ncol(matrix)){
    # print(colnames(matrix)[i])
    for(j in 1:nrow(matrix)){
      if(matrix[j,i]!=0){
        break;    	# break the loop if the entry is not 0
      }
      # this if will only be entered if the last entry was 0
      if(j==nrow(matrix)){
        x <- rbind(x,i)
      }
    } # end nested for loop
  } # end for loop
  # remove all columns with NAs only
  
  if(length(x)>0){
    matrix <- matrix[,-x]
  }
  
  matrix   # return variable
  
} # end function

# removing all rows with NAs only

remove_zero_row <- function(matrix){
  x  <- NULL
  for(i in 1:nrow(matrix)){
    for(j in 1:ncol(matrix)){
      if(matrix[i,j]!=0){
        break;			# break the loop if the entry is not 0
      }
      # this if will only be entered if the last entry was 0
      if(j==ncol(matrix)){
        x <- rbind(x,i)
      }
    } # end nested for loop
  } # end for loop
  # remove all rows with NAs only
  if(length(x)>0){
    matrix <- matrix[-x,]
  }
  matrix   # return variable
} # end function
################################################################################################################################################################
# removing all columns with NAs only

remove_na_col <- function(matrix){
  x  <- NULL
  for(i in 1:ncol(matrix)){
    for(j in 1:nrow(matrix)){
      if(is.na(matrix[j,i])==F){
        break;			# break the loop if the entry is not NA
      }
      # this if will only be entered if the last entry was NA
      if(j==nrow(matrix)){
        x <- rbind(x,i)
      }
    } # end nested for loop
  } # end for loop
  # remove all columns with NAs only
  
  if(length(x)>0){
    matrix <- matrix[,-x]
  }
  
  matrix   # return variable
  
} # end function

# removing all rows with NAs only

remove_na_row <- function(matrix){
  x  <- NULL
  for(i in 1:nrow(matrix)){
    for(j in 1:ncol(matrix)){
      if(is.na(matrix[i,j])==F){
        break;			# break the loop if the entry is not NA
      }
      # this if will only be entered if the last entry was NA
      if(j==ncol(matrix)){
        x <- rbind(x,i)
      }
    } # end nested for loop
  } # end for loop
  # remove all rows with NAs only
  if(length(x)>0){
    matrix <- matrix[-x,]
  }
  matrix   # return variable
} # end function
################################################################################################################################################################
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
    x <- which(tolower(colnames(dataset.2))==tolower(colnames(dataset.1))[i])
    
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
################################################################################################################################################################
# function converting all entries in a matrix to numeric instead of character

# matrix parameter is the matrix to convert
convert_to_num <- function(matrix){
  
  # empty matrix
  conv_mat <- matrix(nrow=nrow(matrix),ncol=ncol(matrix),dimnames=list(rownames(matrix),colnames(matrix)))
  
  # fill matrix con_mat with values from matrix
  for(i in 1:nrow(matrix)){
    for(j in 1:ncol(matrix)){
      conv_mat[i,j] <- as.numeric(as.character(matrix[i,j]))
    } # end nested for loop
    
  } # end for loop
  
  # return variable
  conv_mat
} # end function