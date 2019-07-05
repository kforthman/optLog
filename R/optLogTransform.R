#' Optimize log transformation
#' 
#' This function finds the optimal transformation for normalization of each of the variables and outputs a matrix of the transformed data.
#' @param mydata The dataset you would like to transform. Must be in vector or matrix form. If given a matrix, the function will transform each column seprately. Works best if columns are named, particularly if you are exporting plots.
#' @param skew_thresh The threshold skew value required for transformation. If the skew of the variable is less than skew_thresh, it will be considered normal and will not be transformed.
#' @param scale If set to TRUE, the resulting transformation will have zero mean and unit variance.
#' @param hist_raw_folder The name of the folder where you would like to save a histogram showing the distribution of the raw data. If you do not wish to save these plots, set to NA.
#' @param hist_trans_folder The name of the folder where you would like to save a histogram showing the distribution of the transformed data. If you do not wish to save these plots, set to NA.
#' @param kurt_folder The name of the folder where you would like to save a plot showing the optimal kurtosis with respect to the transformation variable. If you do not wish to save these plots, set to NA.
#' @param n_trans_val The number of gridpoints representing different strengths of transformation we want to test for getting the most normal curve. The higher this number is, the better nomalcy you will get, but the function will also take longer to run.
#' @export
#' @examples 
#' data("USArrests")
#' optLogTransform(USArrests$Murder)

optLogTransform <- function(mydata, skew_thresh = 1, n_trans_val = 50, scale = T, hist_raw_folder = NA,  hist_trans_folder = NA, kurt_folder = NA){
  
  mydata <- as.matrix(mydata)
  
  if(sum(!is.numeric(mydata))>0){
    stop("Please input a numeric dataset.")
  }
  if(sum(mydata<0)>0){
    stop("The transformation only currently works with positive data.")
  }
  
  # Define the number of variables.
  n <- dim(mydata)[1]
  nvar <- dim(mydata)[2]

  # Create an empty matrix to store the transformed variables.
  trans_data <- matrix(nrow = n, ncol = 0)
  
  # Data is bent so that is is curved in a particular direction.
  # Right skewed data will be bent so that they curve concave down.
  # Left skewed data will be bent so that they curve concave up.
  for(i in 1:nvar){
    var_name <- names(mydata)[i]
    var_obs <- mydata[,i]
    
    var_obs_skewness <- round(skewness(var_obs, na.rm = T), 1)
    
    message( paste0(i, '_', var_name, ' ... Skewness: ', var_obs_skewness) )
    
    m <- max(var_obs, na.rm = T)
    # The transformation values alter the strength of the transformation. The larger the 
    # trans_val is, the less strong the transformation will be.
    trans_val <- seq(0, m, length.out = n_trans_val+1)
    trans_val <- trans_val[2:(n_trans_val+1)]
    # The trans values are given a logarithmic step scale. This is because
    # smaller values have a larger incremental effect on the transformation
    # strength.
    trans_val <- (m*(exp(5*trans_val/m) - 1)) / (exp(5) - 1)
    
    # Now, the function does one transformation with variables skewed to the right
    # and the opposite transformation for variables skewed to the left.
    if(var_obs_skewness > skew_thresh){
      
      # The function that transforms the variables has different effects for
      # different values of the trans_val. For smaller values of trans_val,
      # the transformation function will have a stronger effect on the variable's
      # distribution.
      trans_func <- function(trans_val){
        h <- m/(log(m+trans_val)-log(trans_val))
        var_obs_trans <- h*(log(var_obs + trans_val) - log(trans_val))
        return(var_obs_trans)
      }
      
      # Each column of the following matrix is a transformation of a different strength.
      var_obs_trans_mat <- sapply(trans_val, trans_func)
      colnames(var_obs_trans_mat) <- trans_val #round(trans_val, 1)
      # The kurtosis is found for each strength of transformation.
      kurt <- sapply(seq(1,n_trans_val), function(x){kurtosis(var_obs_trans_mat[,x], na.rm = T)})
      
      # The trans_val that results in the best kurtosis is found
      kurt_min <- min(abs(kurt), na.rm = T)
      kurt_min_loc <- which.min(abs(kurt))
      
      # Record the optimal trans val and the resulting kurtosis of its distribution
      mytrans <- trans_val[kurt_min_loc]
      mykurt <- kurt[kurt_min_loc]
      
      # Record the transformed variable.
      var_obs_trans <- as.matrix(var_obs_trans_mat[,kurt_min_loc])
      h <- m/(log(m + mytrans)-log(mytrans))
      mytext <- paste0('right skewed; transformation = ', round(h,3), ' * (log(x + ',
                       round(mytrans,3),
                       ') - log(', round(mytrans,3), ')) ; kurtosis = ',
                       round(mykurt,3))
    }else if(var_obs_skewness < -skew_thresh){
      
      # This function is the same as the one above, but flipped on the y=x axis.
      trans_func <- function(trans_val){
        h <- m/(log(m+trans_val)-log(trans_val))
        var_obs_trans <-  exp(var_obs/h + log(trans_val)) - trans_val
        return(var_obs_trans)
      }
      
      var_obs_trans_mat <- sapply(trans_val, trans_func)
      colnames(var_obs_trans_mat) <- trans_val #round(trans_val, 1)
      kurt <- sapply(seq(1,n_trans_val), function(x){kurtosis(var_obs_trans_mat[,x], na.rm = T)})
      
      kurt_min <- min(abs(kurt), na.rm = T)
      kurt_min_loc <- which.min(abs(kurt))
      
      mytrans <- trans_val[kurt_min_loc]
      mykurt <- kurt[kurt_min_loc]
      
      var_obs_trans <- as.matrix(var_obs_trans_mat[,kurt_min_loc])
      h <- m/(log(m + mytrans)-log(mytrans))
      mytext <- paste0('left skewed; transformation = e^(x/',
                       round(h,3), ' + log(',
                       round(mytrans,3), ')) - ', round(mytrans,3), ' ; kurtosis = ',
                       round(mykurt,3))
    }else{
      message("\tVariable is already normal; no need for transformation.")
      # If the variable is not that skewed, it is not transformed.
      var_obs_trans <- var_obs
      kurt <- seq(0, 0, length.out = n_trans_val)
      mykurt <- kurtosis(var_obs_trans, na.rm = T)
      mytext <- paste0('normal, no transformation; kurtosis = ', round(mykurt,3))
    }
    
    # If the observations are supposed to be of zero mean and unit variance, they are scaled.
    if(scale){
      var_obs_trans <- scale(var_obs_trans)
    }
    # # The following is to assure that the breaks are uniform accross variables.
    # max_var_obs_trans <- max(var_obs_trans, na.rm = T)
    # min_var_obs_trans <- min(var_obs_trans, na.rm = T)
    # breakpoints_var_obs_trans <- seq(min_var_obs_trans, max_var_obs_trans, length.out = 50)
    
    # Create histograms for the transformed variables.
    if(!is.na(hist_trans_folder)){
      png(paste0(hist_trans_folder, "/", i, "_", var_name, "_hist.png"))
      hist(var_obs_trans, main = paste(var_name, "Transformed"), 
           breaks = 50)
      mtext(mytext, cex = 0.9)
      dev.off()
    }
    
    # Create histograms for the raw variables.
    if(!is.na(hist_raw_folder)){
      png(paste0(hist_raw_folder, "/", i, "_", var_name, "_hist.png"))
      hist(var_obs, main = paste(var_name), 
           breaks = 50)
      #mtext(mytext, cex = 0.9)
      dev.off()
    }
    
    if(!is.na(kurt_folder)){
      # The kurtosis for each transformation value.
      png(paste0(kurt_folder, "/", i, "_", var_name, "_kurt.png"))
      plot(trans_val, kurt, main = var_name)
      dev.off()
    }
    
    # Binds the transformed variable to a new matrix.
    trans_data <- cbind(trans_data, matrix(var_obs_trans))
  }
  
  # trans_data <- knnImputation(trans_data)
  
  colnames(trans_data) <- names(mydata)
  rownames(trans_data) <- rownames(mydata)
  return(trans_data)
}