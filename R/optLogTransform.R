#' Optimize log transformation
#'
#' This function finds the optimal transformation for normalization of each of the variables and outputs a matrix of the transformed data. If "scaled" is set to false, unless "retain_domain" is set to FALSE, the domain will be [0,1].
#' The output is a list containing four objects:
#' \enumerate{
#' \item a string vector listing the function used to transform each variable,
#' \item a numeric vector giving the skew of each transformed variable,
#' \item a numeric vector giving the optimal transform value for each variable,
#' \item and a matrix of the transformed data.
#' }
#' @param mydata The dataset you would like to transform. Must be in vector or matrix form. If given a matrix, the function will transform each column seprately. Works best if columns are named, particularly if you are exporting plots.
#' @param type The type of transformation can be either logarithmic or power; "log" and "power" respectively.
#' @param skew_thresh The threshold skew value required for transformation. If the skew of the variable is less than skew_thresh, it will be considered normal and will not be transformed.
#' #' @param n_trans_val The number of gridpoints representing different strengths of transformation we want to test for getting the most normal curve. The higher this number is, the better noimportFrom("grDevices", "dev.off", "png")
#' @param scaled If set to TRUE, the resulting transformation will have zero mean and unit variance.
#' @param retain_domain Set to TRUE if you would like the transformed data to have the same domain as the original dataset (not recommended).
#' @param hist_raw_folder The name of the folder where you would like to save a histogram showing the distribution of the raw data. If you do not wish to save these plots, set to NA.
#' @param hist_trans_folder The name of the folder where you would like to save a histogram showing the distribution of the transformed data. If you do not wish to save these plots, set to NA.
#' @param skew_folder The name of the folder where you would like to save a plot showing the optimal skew with respect to the transformation variable. If you do not wish to save these plots, set to NA.
#' @importFrom grDevices dev.off png
#' @importFrom graphics hist mtext plot
#' @export
#' @examples
#' library(optLog)
#'
#' # First generate a random normal dataset.
#' mydata <- rnorm(100, mean = 0, sd = 1)
#' hist(mydata)
#' # Add skew to the dataset.
#' mydata_skew <- cbind(1-(1-mydata)^2, (1-mydata)^2, mydata)
#' colnames(mydata_skew) <- c("Variable 1", "Variable 2", "Variable 3")
#' for(i in 1:3){hist(mydata_skew[,i])}
#'
#' # Use optLogTransform to remove the skew.
#' mydata_transformed <- optLogTransform(mydata_skew, type = "power", scaled = FALSE)
#' for(i in 1:3){hist(mydata_transformed$data[,i])}

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
  trans_funcs <- matrix(nrow = 1, ncol = 0)
  trans_skews <- matrix(nrow = 1, ncol = 0)
  trans_opts <- matrix(nrow = 1, ncol = 0)


  # ---- Transform the data  ----

  # Data is bent so that is is curved in a particular direction.
  # Right skewed data will be bent so that they curve concave down.
  # Left skewed data will be bent so that they curve concave up.
  for(i in 1:nvar){
    var_name <- colnames(mydata)[i]
    if(is.null(var_name)){ var_name <- paste0("var_", i)}

    var_obs <- mydata[,i]

    var_obs_skewness <- round(e1071::skewness(var_obs, na.rm = T), 1)

    message( paste0(var_name, ' ... Skewness: ', var_obs_skewness) )



    trans_func <- "x"

    # The transformation values alter the strength of the transformation. The larger the
    # trans_val is, the less strong the transformation will be.
    trans_val <- seq(0, 1, length.out = n_trans_val+1)
    trans_val <- trans_val[2:(n_trans_val+1)]

    if(abs(var_obs_skewness) > skew_thresh){

      # Now, the function does one transformation with variables skewed to the right
      # and the opposite transformation for variables skewed to the left.
      if(var_obs_skewness < 0){
        var_obs <- -var_obs
        trans_func <- paste0("-", trans_func)
      }


      m <- max(var_obs, na.rm = T)
      u <- min(var_obs, na.rm = T)
      if(type == 'power'){
        var_obs <- (var_obs-u)/(m-u)
        trans_func <- paste0(round(1/(m-u),3), "(", trans_func, " - ", round(u,3), ")")
        }
      if(type == 'log'){trans_val <- trans_val*m - trans_val*u - u}

      # The function that transforms the variables has different effects for
      # different values of the trans_val. For smaller values of trans_val,
      # the transformation function will have a stronger effect on the variable's
      # distribution.

      # Each column of the following matrix is a transformation of a different strength.
      if(type == 'power'){
        skew_val <- sapply(trans_val, function(x){e1071::skewness(var_obs^x, na.rm = T)^2})
      }else if(type == 'log'){
        skew_val <-  sapply(trans_val, function(x){e1071::skewness(log(var_obs+x), na.rm = T)^2})
      }

      # The trans_val that results in the best skew is found
      trans_opt <- trans_val[which.min(skew_val)]

      # Record the optimal trans val and the resulting kurtosis of its distribution
      if(type == 'power'){
        var_trans <- var_obs^trans_opt
        trans_func <- paste0("(", trans_func, ")^", trans_opt)
      }else if(type == 'log'){
        #var_trans <- (log(var_obs+trans_opt)-log(trans_opt))/(log(1+trans_opt)-log(trans_opt))
        var_trans <- log(var_obs + trans_opt)
        #trans_func <- paste0("( log(", trans_func, " + ", trans_opt, ") - log(", trans_opt, ") ) / ( log(1 + ", trans_opt, ") - log(", trans_opt, ") )")
        trans_func <- paste0("log(", trans_func, " + ", trans_opt, ")")
      }

      if(var_obs_skewness < 0){
        var_trans <- -var_trans
        trans_func <- paste0("-", trans_func)
      }


      if(retain_domain){
        var_trans <-  var_trans * (m - u) + u
        trans_func <- paste0("( ", trans_func, " ) * (", round(m, 3), " - ", round(u, 3), ") + ", round(u, 3))
      }


      if(!is.na(skew_folder)){
        # The kurtosis for each transformation value.
        png(paste0(skew_folder, "/", i, "_", var_name, "_skew.png"))
        plot(trans_val, skew_val, main = var_name, type = 'l')
        abline(a = 0, b = 0, col = "red")
        dev.off()
      }

    }else{
      # If the variable is not that skewed, it is not transformed.
      message("\tVariable is already normal; no need for transformation.")
      var_trans <- var_obs
      skew_val <- seq(0, 0, length.out = n_trans_val)
    }

    # If the observations are supposed to be of zero mean and unit variance, they are scaled.
    if(scaled){
      var_trans <- scale(var_trans)
      trans_func <- paste0("scale(", trans_func, ")")
    }

    if(is.na(trans_func)){
      mytext <- paste0('normal, no transformation')
    }else{
      mytext <- trans_func
    }
    # ---- Create Histograms ----

    # Create histograms for the transformed variables.
    if(!is.na(hist_trans_folder)){
      png(paste0(hist_trans_folder, "/", i, "_", var_name, "_hist.png"))
      hist(var_trans, main = paste(var_name, "Transformed"),
           breaks = 25)
      mtext(mytext, cex = 0.9)
      dev.off()
    }

    # Create histograms for the raw variables.
    if(!is.na(hist_raw_folder)){
      png(paste0(hist_raw_folder, "/", i, "_", var_name, "_hist.png"))
      hist(mydata[,i], main = paste(var_name),
           breaks = 25)
      #mtext(mytext, cex = 0.9)
      dev.off()
    }



    # Binds the transformed variable to a new matrix. ----
    trans_data <- cbind(trans_data, matrix(var_trans))
    trans_funcs <- cbind(trans_funcs, trans_func)
    trans_opts <- cbind(trans_opts, trans_opt)
    trans_skews <- cbind(trans_skews, e1071::skewness(var_trans))
  }

  # trans_data <- knnImputation(trans_data)

  # ---- Return transformed dataset ----
  colnames(trans_data) <- names(mydata)
  rownames(trans_data) <- rownames(mydata)

  output <- list(names = colnames(mydata), func = trans_funcs, skew = trans_skews, trans_val = trans_opts, data = trans_data)

  return(output)
}
