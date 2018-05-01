plot_factors <- function(scores, type = "Umean",
            clust = if (exists('clust', where=scores)) scores$clust else NULL,
            limx = NULL, limy = NULL) {

p_gg <- requireNamespace("GGally", quietly = TRUE) &&
  requireNamespace("ggplot2", quietly = TRUE)

if (type == "Umean") {
  mat <- scores$Umean
}
if (type == "Uclust") {
  mat <- scores$Uclust
}
if (type == "Uscores") {
  mat <- scores$Uscores
}
if (is.null(clust)) {
  if (type == "Uscores") {
    q <- dim(mat)[2]
    g <- dim(mat)[3]
      for (i in 1 : g) {
        mat <- as.matrix(scores$Uscores[,, i])
        colnames(mat) <- paste0("u_", 1 : q)

        if (q == 1) {
          plot(mat, 1 : length(mat), xlim = range(c(mat, limx)), axes = FALSE,
                 xlab = expression(widehat(u)[1]), ylab = "", type = "p",
                 pch = 20, bg = 1)
          axis(side=1)
        }

        if (p_gg) {
          df <- data.frame(mat)
          print(GGally::ggpairs(df, columns = 1 : q))

        } else {
          if (q == 2)
            plot(mat[, c(1, 2)], ylim = range(c(limy, mat[, 2])),
               pch = 20, bg =  1,  xlim = range(c(limx, mat[, 1])),
               xlab = expression(widehat(u)[1]),
               ylab = expression(widehat(u)[2]))

          if (q > 2)
            pairs(mat, bg =  1, pch = 20)
        }

        readline(prompt = "Press [enter] to continue")
      }

  } else if((type == "Umean") || (type == "Uclust")) {
        stop('For type= "Umean" or "Uclust", clust needs to be specified')
  }

} else {
  q <- dim(mat)[2]
  g <- length(unique(clust))
  if ((type == "Umean") || (type == "Uclust"))
    it  <- 1

  if (type == "Uscores") {
    it <- g <- dim(mat)[3]
  }

  for (i in 1:it) {

    if (type == "Uscores")
      mat <- as.matrix(scores$Uscores[,, i])

    if (q == 1) {
      plot(mat, 1 : length(mat),  xlim = range(c(limx, mat)), axes = FALSE,
           xlab = expression(widehat(u)[1]), ylab = "", type = "p",
           pch = if (g <= 5) {20 + as.numeric(clust)} else{as.numeric(clust)},
           col = as.numeric(clust), bg =  as.numeric(clust))
      axis(side = 1)
    }

    if (p_gg) {
      colnames(mat) <- paste0("u_", 1 : q)
      df <- data.frame(mat, clust=factor(clust))
      print(GGally::ggpairs(df, columns = 1 : q,
                      mapping = ggplot2::aes(colour = clust)))

    } else {
      if (q == 2)
        plot(mat[, c(1, 2)], col = as.numeric(clust),
           ylim = range(c(limy, mat[, 2])),
           pch = if (g <= 5) {20 + as.numeric(clust)} else{as.numeric(clust)},
           bg = as.numeric(clust), xlim = range(c(limx, mat[, 1])),
           xlab = expression(widehat(u)[1]), ylab=expression(widehat(u)[2]))

      if (q > 2)
        pairs(mat, col = as.numeric(clust), bg = as.numeric(clust),
          pch = if (g <= 5) {20 + as.numeric(clust)} else {as.numeric(clust)})
    }

    if (i < it)
      readline(prompt="Press [enter] to continue")
  }
}
}
