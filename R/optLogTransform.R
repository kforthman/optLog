#' Optimize log transformation
#' 
#' This function finds the optimal transformation for normalization of each of the variables and outputs a matrix of the transformed data. If "scaled" is set to false, unless "retain_domain" is set to FALSE, the domain will be [0,1].
#' @param mydata The dataset you would like to transform. Must be in vector or matrix form. If given a matrix, the function will transform each column seprately. Works best if columns are named, particularly if you are exporting plots.
#' @param type The type of transformation can be either logarithmic or power; "log" and "power" respectively.
#' @param skew_thresh The threshold skew value required for transformation. If the skew of the variable is less than skew_thresh, it will be considered normal and will not be transformed.
#' @param scaled If set to TRUE, the resulting transformation will have zero mean and unit variance.
#' @param hist_raw_folder The name of the folder where you would like to save a histogram showing the distribution of the raw data. If you do not wish to save these plots, set to NA.
#' @param hist_trans_folder The name of the folder where you would like to save a histogram showing the distribution of the transformed data. If you do not wish to save these plots, set to NA.
#' @param skew_folder The name of the folder where you would like to save a plot showing the optimal skew with respect to the transformation variable. If you do not wish to save these plots, set to NA.
#' @param n_trans_val The number of gridpoints representing different strengths of transformation we want to test for getting the most normal curve. The higher this number is, the better nomalcy you will get, but the function will also take longer to run.
#' @export
#' @examples 
#' library(optLog)
#' 
#' # First generate a random normal dataset.
#' mydata <- rnorm(100, mean = 0, sd = 1)
#' hist(mydata)
#' # Add skew to the dataset.
#' mydata_skew <- 1-(1-mydata)^2
#' hist(mydata_skew)
#' 
#' # Use optLogTransform to remove the skew.
#' mydata_transformed <- optLogTransform(as.matrix(mydata_skew), type = "power", scaled = F)
#' hist(mydata_transformed)

optLogTransform <- function(mydata, type = 'log', skew_thresh = 1, n_trans_val = 50, scaled = T, retain_domain = F, hist_raw_folder = NA,  hist_trans_folder = NA, skew_folder = NA){
  
  # ---- Conditions for function use ----
  
  if(scaled && retain_domain){
    stop("You cannot set both scaled and retain_domain to TRUE.")
  }
  if(!is.matrix(mydata)){
    stop("Please input dataset in matrix form.")
  }
  
  if(sum(!is.numeric(mydata))>0){
    stop("Please input a numeric dataset.")
  }
  
  
  # ---- Setup  ----
  
  # Define the number of variables.
  n <- dim(mydata)[1]
  nvar <- dim(mydata)[2]
  
  # Create an empty matrix to store the transformed variables.
  trans_data <- matrix(nrow = n, ncol = 0)
  
  
  # ---- Transform the data  ----
  
  # Data is bent so that is is curved in a particular direction.
  # Right skewed data will be bent so that they curve concave down.
  # Left skewed data will be bent so that they curve concave up.
  for(i in 1:nvar){
    var_name <- names(mydata)[i]
    if(is.null(var_name)){ var_name <- paste0("var_", i)}
    
    var_obs <- mydata[,i]
    
    var_obs_skewness <- round(e1071::skewness(var_obs, na.rm = T), 1)
    
    message( paste0(i, '_', var_name, ' ... Skewness: ', var_obs_skewness) )
    
    m <- max(var_obs, na.rm = T)
    u <- min(var_obs, na.rm = T)
    
    # Change domain to be between 0 and 1.
    var_obs <- (var_obs - u)/(m - u)
    
    # The transformation values alter the strength of the transformation. The larger the 
    # trans_val is, the less strong the transformation will be.
    trans_val <- seq(0, 1, length.out = n_trans_val+1)
    trans_val <- trans_val[2:(n_trans_val+1)]
    
    if(abs(var_obs_skewness) > skew_thresh){
      
      # Now, the function does one transformation with variables skewed to the right
      # and the opposite transformation for variables skewed to the left.
      if(var_obs_skewness < 0){ var_obs <- -var_obs + 1 }
      
      # The function that transforms the variables has different effects for
      # different values of the trans_val. For smaller values of trans_val,
      # the transformation function will have a stronger effect on the variable's
      # distribution.
      
      # Each column of the following matrix is a transformation of a different strength.
      if(type == 'power'){
        skew_val <- sapply(trans_val, function(x){e1071::skewness(var_obs^x, na.rm = T)^2})
      }else if(type == 'log'){
        skew_val <-  sapply(trans_val, function(x){e1071::skewness(log(var_obs+trans_val), na.rm = T)^2})
      }
      
      # The trans_val that results in the best skew is found
      trans_opt <- trans_val[which.min(skew_val)]
      
      # Record the optimal trans val and the resulting kurtosis of its distribution
      if(type == 'power'){
        var_trans <- var_obs^trans_opt
      }else if(type == 'log'){
        var_trans <- (log(var_obs+trans_opt)-log(trans_opt))/(log(1+trans_opt)-log(trans_opt))
      }
      
      if(var_obs_skewness < 0){ var_trans <- -var_trans+1 }
      
      
      if(retain_domain){ var_trans <-  var_trans * (m - u) + u}
      
      if(var_obs_skewness > 0){ 
        if(type == 'power'){
          mytext <- paste0("x^", trans_opt)
        }else if(type == 'log'){
          mytext <- paste0("(log(x + ", trans_opt, ") - log(", trans_opt, "))/(log(1+", trans_opt, ")-log(", trans_opt, "))")
        }
      }else if(var_obs_skewness < 0){ 
        if(type == 'power'){
          mytext <- paste0("1-(1-x)^", trans_opt)
        }else if(type == 'log'){
          mytext <- paste0("1 - (log(1- x + ", trans_opt, ") - log(", trans_opt, "))/(log(1+", trans_opt, ")-log(", trans_opt, "))")
        } }
      
    }else{
      # If the variable is not that skewed, it is not transformed.
      message("\tVariable is already normal; no need for transformation.")
      var_trans <- var_obs
      skew_val <- seq(0, 0, length.out = n_trans_val)
      mytext <- paste0('normal, no transformation')
    }
    
    # If the observations are supposed to be of zero mean and unit variance, they are scaled.
    if(scaled){
      var_trans <- scale(var_trans)
    }
    
    # ---- Create Histograms ----
    
    # Create histograms for the transformed variables.
    if(!is.na(hist_trans_folder)){
      png(paste0(hist_trans_folder, "/", i, "_", var_name, "_hist.png"))
      hist(var_trans, main = paste(var_name, "Transformed"), 
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
    
    if(!is.na(skew_folder)){
      # The kurtosis for each transformation value.
      png(paste0(skew_folder, "/", i, "_", var_name, "_kurt.png"))
      plot(trans_val, skew_val, main = var_name)
      dev.off()
    }
    
    # Binds the transformed variable to a new matrix. ----
    trans_data <- cbind(trans_data, matrix(var_trans))
  }
  
  # trans_data <- knnImputation(trans_data)
  
  # ---- Return transformed dataset ----
  colnames(trans_data) <- names(mydata)
  rownames(trans_data) <- rownames(mydata)
  return(trans_data)
}