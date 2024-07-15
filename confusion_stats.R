
# if you need to check these later, look at obsidian ROC curve note
# obsidian://open?vault=simulations&file=ROC%20curves
get_TPR <- function(data, mode = "ignore_sign", return_components = FALSE) {
    if (mode == "cluster") {
        TP <- data$TP_cluster 
        FN <- data$FN_cluster
    } else if (mode == "ignore_sign") {
        TP <- data$TP_ignoreSign
        FN <- data$FN_ignoreSign
    } else if (mode == "ignore_direction") {
        TP <- data$TP_ignoreDirection
        FN <- data$FN_ignoreDirection
    } else if (mode == "sign") {
        # tp/(tp+fn)
        TP <- data$TP_sign
        FN <- data$FN_sign
    } else if (mode == "cluster_cov") {
        # tp/(tp+fn)
        TP <- data$TP_clusterCov
        FN <- data$FN_clusterCov
    } else {
        print("mode get_TPR not implemented")
    }
    TPR <- TP / (TP + FN)
    if (return_components) {
      return(list(TPR = TPR, TP = TP, FN = FN))
    } else {
      return(TPR)
    }
}

get_FPR <- function(data, mode = "ignore_sign", return_components = FALSE) {
    if (mode == "cluster") {
      FP <- data$FP_cluster
      TN <- data$TN_cluster
    } else if (mode == "ignore_sign") {
      FP <- data$FP_ignoreSign
      TN <- data$TN_ignoreSign
    } else if (mode == "ignore_direction") {
      FP <- data$FP_ignoreDirection
      TN <- data$TN_ignoreDirection
    } else if (mode == "sign") {
      FP <- data$FP_sign
      TN <- data$TN_sign
    } else if (mode == "cluster_cov") {
      FP <- data$FP_clusterCov
      TN <- data$TN_clusterCov
    } else {
        print("mode get_FPR not implemented")
    }
    FPR <- FP / (FP + TN)
    if (return_components) {
      return(list(FPR = FPR, FP = FP, TN = TN))
    } else {
      return(FPR)
    }
}

get_precision <- function(data, mode = "ignore_sign", return_components = FALSE) {
    if (mode == "cluster") {
        TP <- data$TP_cluster 
        FP <- data$FP_cluster
    } else if (mode == "ignore_sign") {
        TP <- data$TP_ignoreSign
        FP <- data$FP_ignoreSign
    } else if (mode == "ignore_direction") {
        TP <- data$TP_ignoreDirection
        FP <- data$FP_ignoreDirection
    } else if(mode == "sign") {
      TP <- data$TP_sign
      FP <- data$FP_sign 
    } else if (mode == "cluster_cov") {
      TP <- data$TP_clusterCov
      FP <- data$FP_clusterCov
    } else {
        print("mode get_precision not implemented")
    }
    prec <- TP / (TP + FP)
    if (return_components) {
      return(list(prec = prec, TP = TP, FP = FP))
    } else {
      return(prec)
    }
}

get_recall <- function(data, mode = "ignore_sign", return_components = FALSE) {
    if (mode == "cluster") {
      TP <- data$TP_cluster
      FN <- data$FN_cluster
    } else if (mode == "ignore_sign") {
      TP <- data$TP_ignoreSign
      FN <- data$FN_ignoreSign 
    } else if (mode == "ignore_direction") {
      TP <- data$TP_ignoreDirection
      FN <- data$FN_ignoreDirection
    } else if (mode == "sign") {      
      TP <- data$TP_sign
      FN <- data$FN_sign 
    } else if (mode == "cluster_cov") {
      TP <- data$TP_clusterCov
      FN <- data$FN_clusterCov
    } else {
        print("mode get_recall not implemented")
    }
    recall <- TP / (TP + FN)
    if (return_components) {
      return(list(recall = recall, TP = TP, FN = FN))
    } else {
      return(recall)
    }
}

get_falseDiscovery <- function(data, mode = "ignore_sign", return_components = FALSE) {
    if (mode == "cluster") {
      TP <- data$TP_cluster
      FP <- data$FP_cluster
    } else if (mode == "ignore_sign") {
      TP <- data$TP_ignoreSign
      FP <- data$FP_ignoreSign 
    } else if (mode == "ignore_direction") {
      TP <- data$TP_ignoreDirection
      FP <- data$FP_ignoreDirection
    } else if (mode == "sign") {      
      TP <- data$TP_sign
      FP <- data$FP_sign 
    } else if (mode == "cluster_cov") {
      TP <- data$TP_clusterCov
      FP <- data$FP_clusterCov
    } else {
        print("mode get_falseDiscovery not implemented")
    }
    FDR <- FP / (FP + TP)
    if (return_components) {
      return(list(FDR = FDR, FP = FP, TP = TP))
    } else {
      return(FDR)
    }
}
