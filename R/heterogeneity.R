#' Calculate Age-related Expression and Heterogeneity Changes for Gene
#'
#' @param exprs a numeric vector with the expression values, ordered the same as age vector
#' @param age a numeric vector, where the names correspond to samples (the same as colnames of the given matrix).
#' @param modex expression change calculation method. 'linear' or 'loess', defaults to 'linear'
#' @param het_change_met heterogeneity change calculation method. 'LR', for 'linear regression', or and correlation method accepted by cor.test() function.
#'
#' @return a list with i) summary results with the effect size and p value for expression level and heterogeneity changes, ii) residuals from expr~age model, which corresponds to unexplained variance, iii) expression values passed into the function.
#' @export
#'
calc.het.feat <- function(exprs, age, modex = 'linear',
                          het_change_met) {
    if ( modex == 'linear') {
        lmx <- lm(exprs ~ age)
        expr_change <- unname(lmx$coef[2])
        expr_change.p <- unname(summary(lmx)$coef[2, 4])
    } else if ( modex == 'loess') {
        lmx <- loess(exprs ~ age)
        expr_change <- NA
        expr_change.p <- NA
    }
    resids <- lmx$resid
    if ( het_change_met == "LR") {
        lmxx <- lm(abs(resids) ~ age)
        het_change <- unname(lmxx$coefficients[2])
        het_change.p <- unname((summary(lmxx)$coef)[2,4])
    } else {
        co_het <- cor.test(abs(resids), age, method = het_change_met)
        het_change <- unname(co_het$est)
        het_change.p <- unname(co_het$p.val)
    }
    res <- list(summary = c(level_change = expr_change,
                            level_change.p = expr_change.p,
                            heterogeneity_change = het_change,
                            heterogeneity_change.p = het_change.p))
    res$residuals <- resids
    res$expression <- exprs
    return(res)
}

#' Calculate Age-related Expression and Heterogeneity Changes
#'
#' @param exprmat a numeric matrix with the expression values, where columns are the samples and rows are probesets, transcripts, or genes.
#' @param age a numeric vector, where the names correspond to samples (the same as colnames of the given matrix).
#' @param age_in type of the age vector, allowed values are days or years. defaults to 'days'
#' @param age_to final format of age vector. allowed values are 'years', 'days', 'pw-N', and 'lg-N', where N is any number. 'pw' means power. e.g. pw-0.5 means sqrt(age), and 'lg' means log, e.g. lg-2 means log2(age).defaults to 'pw-0.25'
#' @param batch_corr batch correction strategy. available values are i) NC: No Correction, ii) QN: Quantile Normalization, iii) LR: Linear regression (requires covariates), iv) SVA: surrogate variable analysis, v) LR+QN: Linear regression followed by quantile normalization, and vi) SVA+QN: SVA followed by quantile normalization. Defaults to 'NC'.
#' @param modex expression change calculation method. 'linear' or 'loess', defaults to 'linear'
#' @param tr_log2 logical to set log2 transforming expression matrix. defaults to TRUE.
#' @param sc_features logical to set whether to scale features. defaults to TRUE.
#' @param covariates a list of covarietes where each element is a vector with sampleIDs as names
#' @param het_change_met heterogeneity change calculation method. 'LR', for 'linear regression', or and correlation method accepted by cor.test() function.
#' @param padj_met method for multiple test correction. value is passed to 'p.adjust' function. defaults to 'fdr'.
#'
#' @return a list object with summary results including expression level and heterogeneity changes, input values, intermediate values, and heterogeneity matrix.
#'
#' @importFrom tibble as.tibble
#' @importFrom dplyr select
#' @importFrom dplyr everything
#'
#' @export
#'
calc.het <- function(exprmat, age, age_in = 'days',
                     age_to = 'pw-0.25',
                     batch_corr = "NC",
                     modex = "linear",
                     tr_log2 = T, sc_features = T,
                     covariates = NA,
                     het_change_met = "spearman",
                     padj_met = "fdr") {
    retval <- list()
    snms <- intersect(colnames(exprmat), names(age))
    retval$sampleID <- snms
    exprmat <- exprmat[, snms]
    retval$input_expr <- exprmat
    age <- age[snms]
    if ( !all(is.na(covariates)) ) { covariates = lapply(covariates,
                                                    function(cov) cov[snms] )}
    retval$input_age <- age
    age <- transform_age(age, type = age_in, to = age_to)
    retval$usedAge <- age
    exprmat <- exprmat[complete.cases(exprmat),]
    if (tr_log2 == T) {
        exprmat <- trans_log2(exprmat)
    }
    if (batch_corr == "NC") {
        exprmat <- exprmat
    } else if (batch_corr == "QN") {
        exprmat <- cfr_quantileNorm(exprmat)
    } else if (batch_corr == "LR") {
        LR_res <- cfr_linReg(exprmat, covariates)
        exprmat <- LR_res$correctedExp
        retval$LR_res <- LR_res
    } else if (batch_corr == "LR+QN") {
        LR_res <- cfr_linReg(exprmat, covariates)
        exprmat <- LR_res$correctedExp
        exprmat <- cfr_quantileNorm(exprmat)
        retval$LR_res <- LR_res
    } else if (batch_corr == "SVA") {
        sva_res <- cfr_sva(exprmat, age, covariates)
        exprmat <- sva_res$correctedExp
        retval$sva_res <- sva_res
    } else if (batch_corr == "SVA+QN") {
        sva_res <- cfr_sva(exprmat, age, covariates)
        exprmat <- sva_res$correctedExp
        exprmat <- cfr_quantileNorm(exprmat)
        retval$sva_res <- sva_res
    }
    exprmat <- feature_scale(exprmat)
    exprmat <- exprmat[complete.cases(exprmat),]
    resx <- apply(exprmat, 1, function(expx) {
        calc.het.feat(expx, age, modex = modex,
                      het_change_met = het_change_met)
    })
    feature_het <- tibble::as.tibble(t(sapply(resx, function(x) x$summary)))
    feature_het$level_change.p.adj <- p.adjust(feature_het$level_change.p,
                                               method = padj_met)
    feature_het$heterogeneity_change.p.adj <- p.adjust(feature_het$heterogeneity_change.p,
                                              method = padj_met)
    feature_het$feature <- rownames(exprmat)
    feature_het <- dplyr::select(feature_het, 'feature', dplyr::everything())
    retval$residMat <- t(sapply(resx,function(x)x$resid))
    retval$usedMat <- t(sapply(resx,function(x)x$expr))
    retval$feature_result <- feature_het
    return(retval)
}
