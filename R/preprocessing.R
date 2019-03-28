#' Log2 Transformation
#'
#' @param mat a numeric matrix with the expression values, where columns are the samples and rows are probesets, transcripts, or genes.
#'
#' @return log2 transformed matrix
#'
#' @export
#'
#' @examples
#' myexp <- sapply(1:10, function(i) rnbinom(10000, mu = 4, size = 1))
#' myexp_l2 <- trans_log2(myexp)
#' par(mfrow = c(1,2))
#' hist(myexp)
#' hist(myexp_l2)
#'
trans_log2 <- function(mat){
    return(log2(mat + 1))
}

#' Quantile Normalize Expression Values
#'
#' A function to apply quantile normalization. It is build on preprocessCore::normalize.quantiles function, but keeps the dimnames as is.
#'
#' @param mat a numeric matrix with the expression values, where columns are the samples and rows are probesets, transcripts, or genes.
#'
#' @return quantile normalized expression values as a matrix
#'
#' @importFrom preprocessCore normalize.quantiles
#'
#' @export
#' @examples
#' myexp <- sapply( 1:10, function(i){
#' rnorm(n = 10000, mean = sample(1:3, 1), sd = sample(c(1, 3), 1))
#' })
#' myexp_qn <- cfr_quantileNorm(myexp)
#' par(mfrow = c(1,2))
#' boxplot(myexp, main = 'Before QN')
#' boxplot(myexp_qn, main = 'After QN')
#'

cfr_quantileNorm <- function(mat) {
    res <- preprocessCore::normalize.quantiles(mat)
    colnames(res) <- colnames(mat)
    rownames(res) <- rownames(mat)
    return(res)
}

#' Construct model matrix
#'
#' Using a list of covarietes, constructs a model matrix to be used in a model
#'
#' @param covs a list of covarietes
#'
#' @return a dataframe that can be used as a design matrix in a model
#'
#' @export
#' @examples
#' mycovs <- list(batch = sample(c(0, 1), 6, replace = TRUE),
#' array = as.factor(sample(c('a', 'b', 'c'), 6, replace = TRUE)),
#' day = sample(c(0, 1, 3), 6, replace = TRUE))
#' des_mat <- const.mod(mycovs)
#' print(mycovs)
#' print(des_mat)
#'
const.mod <- function(covs) {
    res <- lapply(covs, function(nm) {
        model.matrix(~0 + nm)
    })
    res2 <- matrix(ncol = 0, nrow = length(covs[[1]]))
    for (i in names(res)) {
        if (ncol(res[[i]]) > 1) {
            colnames(res[[i]]) <- paste(i, levels(covs[[i]]), sep = ".")
        } else {
            colnames(res[[i]]) <- i
        }
        res2 <- cbind(res2, res[[i]])
    }
    return(res2)
}

#' Removing Confounding Factor Effects using Linear Regression
#'
#' @param mat a numeric matrix with the expression values, where columns are the samples and rows are probesets, transcripts, or genes.
#' @param cov a list with covariates to regress out from the expression matrix
#'
#' @return a list with three elements, i) corrected expression matrix, ii) coefficients for the covariates and iii) p values for each covariate
#'
#' @export
#'

cfr_linReg <- function(mat, cov){
    cov <- const.mod(cov)
    temp_mat = apply(mat, 1, function(x){
        lmx = lm(x ~ cov)
    })
    corr_mat = t(sapply(temp_mat, function(x) x$resid ))
    colnames(corr_mat) <- colnames(mat)
    rownames(corr_mat) <- rownames(mat)
    cov_coef <- t(sapply(temp_mat, function(x) x$coeff ))
    colnames(cov_coef) <- gsub('covsx', '', colnames(cov_coef))
    cov_p <- t(sapply(temp_mat, function(x) (summary(x)$coef)[,4] ))
    colnames(cov_p) <- gsub('covsx', '', colnames(cov_p))
    res <- list(correctedExp = corr_mat,
                cov_coef = cov_coef,
                cov_p = cov_p)
    return(res)
}

#' Removing Confounding Factor Effects using SVA
#'
#' @param mat a numeric matrix with the expression values, where columns are the samples and rows are probesets, transcripts, or genes.
#' @param age a numeric vector, where the names correspond to samples (the same as colnames of the given matrix).
#' @param cov covariates to analyse SVs.
#'
#' @return a list with i) a numeric matrix with the expression values, ii) SVs for each sample, and iii) correlation between SVs and covariates
#'
#' @importFrom sva sva
#' @importFrom reshape2 melt
#'
#' @export


cfr_sva <- function(mat, age, cov = NA) {
    samples <- unique(intersect(colnames(mat), names(age)))
    mat <- mat[, samples]
    age <- age[samples]
    mod <- model.matrix(~age, data = as.data.frame(colnames(mat)))
    mod0 <- model.matrix(~1, data = as.data.frame(colnames(mat)))
    svobj <- sva::sva(mat, mod, mod0)
    if (!all(is.na(cov))) {
        sv_cov_cors <- apply(svobj$sv, 2, function(x) {
            cx <- const.mod(cov)
            lmx <- lm(x ~ cx)
            estx <- lmx$coef[-1]
            names(estx) <- substr(names(estx),3,nchar(names(estx)))
            px <- summary(lmx)$coef[-1, 4]
            names(px) <- substr(names(px),3,nchar(names(px)))
            resultx <- data.frame(estimate = estx[colnames(cx)],
                       p = px[colnames(cx)])
            resultx$p.adj <- p.adjust(resultx$p,'fdr')
            resultx$cov <- colnames(cx)
            resultx
        })
        sv_cov_cors <- reshape2::melt(sv_cov_cors,id.vars = c('cov','estimate',
                                                              'p', 'p.adj'))
        colnames(sv_cov_cors)[5] = 'SV'
    } else {
        sv_cov_cors <- NULL
    }
    corr.mat <- t(apply(mat, 1, function(genexp) lm(genexp ~ svobj$sv)$resid))
    SVs <- svobj$sv
    colnames(SVs) = paste('SV', 1:ncol(SVs), sep = '_')
    rownames(SVs) = samples
    res <- list(correctedExp = corr.mat, SVs = SVs, SV_cov_corr = sv_cov_cors)
    return(res)
}

#' Scale Expression Values for Features (probsets, transcripts, genes)
#'
#' @param mat a numeric matrix with the expression values, where columns are the samples and rows are probesets, transcripts, or genes.
#'
#' @return scaled matrix
#' @export
#'
feature_scale <- function(mat) {
    res <- t(apply(mat, 1, scale))
    colnames(res) <- colnames(mat)
    rownames(res) <- rownames(mat)
    return(res)
}

#' Age Transforming Function
#'
#' @param age a numeric vector with ages
#' @param type type of the age vector, allowed values are days or years. defaults to 'days'
#' @param to final format of age vector. allowed values are 'years', 'days', 'pw-N', and 'lg-N', where N is any number. 'pw' means power. e.g. pw-0.5 means sqrt(age), and 'lg' means log, e.g. lg-2 means log2(age). Defaults to 'years'.
#'
#' @return a vector with ages in the desired scale
#' @export
#'
#' @examples
#' agevec <- seq(1:20)
#' transform_age(agevec, type = 'years', to = 'days')
#' transform_age(agevec, type = 'years', to = 'pw-0.25')
#' transform_age(agevec, type = 'years', to = 'lg-10')
#'
transform_age <- function(age, type = 'days',
                          to = 'years') {
    if (type == 'years') {
        age <- age * 365
    }
    if (to == 'days') {
        return(age)
    } else if (to == 'years') {
        return(age/365)
    } else if (substr(to, 1, 2) == 'pw') {
        trval <- as.numeric(strsplit(to, '-')[[1]][2])
        return(age^trval)
    } else if (substr(to, 1, 2) == 'lg') {
        trval <- as.numeric(strsplit(to, '-')[[1]][2])
        return(log(age,base = trval))
    }
}
