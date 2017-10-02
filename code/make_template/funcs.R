## functions used to create / modify template structure

## creates a new band with name newb
## newb is created as a copy of oldb
AddBand <- function(tem,newb,oldb){
    tem$betas <- cbind(tem$betas,tem$betas[,oldb])
    colnames(tem$betas)[ncol(tem$betas)] <- newb
    tem$betas <- tem$betas[,order(colnames(tem$betas))]
    ## dust
    tem$dust[newb] <- tem$dust[oldb]
    tem$dust <- tem$dust[order(names(tem$dust))]
    ## model error
    tem$model_error[newb] <- tem$model_error[oldb]
    tem$model_error <- tem$model_error[order(names(tem$model_error))]
    ## templates
    te <- matrix(0,nrow=nrow(tem$templates)+1,ncol=ncol(tem$templates))
    rownames(te) <- c(rownames(tem$templates),newb)
    te[1:nrow(tem$templates),] <- tem$templates
    te[newb,] <- te[oldb,]
    te <- te[order(rownames(te)),]
    tem$templates <- te
    ## templatesd
    te <- matrix(0,nrow=nrow(tem$templatesd)+1,ncol=ncol(tem$templatesd))
    rownames(te) <- c(rownames(tem$templatesd),newb)
    te[1:nrow(tem$templatesd),] <- tem$templatesd
    te[newb,] <- te[oldb,]
    te <- te[order(rownames(te)),]
    tem$templatesd <- te
    ## tempalate_funcs
    tem$template_funcs[newb] <- tem$template_funcs[oldb]
    tem$template_funcs <- tem$template_funcs[order(names(tem$template_funcs))]
    ## tempalated_funcs
    tem$templated_funcs[newb] <- tem$templated_funcs[oldb]
    tem$templated_funcs <- tem$templated_funcs[order(names(tem$templated_funcs))]
    return(tem)
}


RemoveBand <- function(tem,rb){
    bands <- colnames(tem$betas)
    bands <- bands[!(bands %in% rb)]
    tem$betas <- tem$betas[,bands]
    tem$dust <- tem$dust[bands]
    tem$templates <- tem$templates[bands,]
    tem$templatesd <- tem$templatesd[bands,]
    tem$template_funcs <- tem$template_funcs[bands]
    tem$templated_funcs <- tem$templated_funcs[bands]
    tem$model_error <- tem$model_error[bands]
    return(tem)
}
