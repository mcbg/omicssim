#' @import data.table distr
#' @export
sim_data <- function(n, group_no, rp, m, m_bm, ref, FUN = rnorm, debug=FALSE) {
  # simulate data
  k <- 2 * n  / group_no # number of samples in each group
  i <- rep(seq_len(group_no), k)

  ds <- .sim_data_helper(n = n, m = m, m_biomarkers = m_bm, seed = rp, mus = ref, FUN = FUN)
  ds$a <- i

  # check balance between trt / plates
  if(!ds[, sum(trt == 1) == sum(trt == 0) , a][, all(V1)]) {
    TMP <<- ds
    stop('balance between trt and groups')
  }
  return(ds)
}


.sim_data_helper <- function(n, # samples per a arm
                     m, # total analytes
                     m_biomarkers = 1, # biomarkers
                     seed = 123,
                     FUN,
                     mus = NULL) {
  set.seed(seed)

  # generate mu
  if (is.null(mus)) {
    mus <- rpois(m, lambda = 1) + 1
  }

  # check if mu = 0
  if(any(mus == 0)) stop('mu must be > 0')

  # generate analytes
  # first m_biomarkers are biomarkers, rest have no treatment effect
  ll <- lapply(seq_along(mus), \(i) {
    mu <- mus[i]
    d <- fifelse(i <= m_biomarkers, 1/2, 0)
    x <- c(FUN(n) - d, d + FUN(n))
    exp(x) * mu
  })
  ds <- as.data.table(ll)
  names(ds)[seq_len(m_biomarkers)] <- sub('^V', 'BM', names(ds)[seq_len(m_biomarkers)])

  # trt
  ds[, trt := rep(c(1, 0) %>% as.character(), each=n)]
  ds[, id := seq_len(nrow(ds))]

  return(copy(ds))
}
