try_udag2pdagRelaxed = function (gInput, verbose = FALSE, unfVect = NULL, solve.confl = FALSE, 
          orientCollider = TRUE, rules = rep(TRUE, 3)) 
{
  orientConflictCollider <- function(pdag, x, y, z) {
    if (pdag[x, y] == 1) {
      pdag[y, x] <- 0
    }
    else {
      pdag[x, y] <- pdag[y, x] <- 2
    }
    if (pdag[z, y] == 1) {
      pdag[y, z] <- 0
    }
    else {
      pdag[z, y] <- pdag[y, z] <- 2
    }
    pdag
  }
  rule1 <- function(pdag, solve.confl = FALSE, unfVect = NULL) {
    search.pdag <- pdag
    ind <- which(pdag == 1 & t(pdag) == 0, arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      isC <- which(search.pdag[b, ] == 1 & search.pdag[, 
                                                       b] == 1 & search.pdag[a, ] == 0 & search.pdag[, 
                                                                                                     a] == 0)
      if (length(isC) > 0) {
        for (ii in seq_along(isC)) {
          c <- isC[ii]
          if (!solve.confl | (pdag[b, c] == 1 & pdag[c, 
                                                     b] == 1)) {
            if (!is.null(unfVect)) {
              if (!any(unfVect == triple2numb(p, a, b, 
                                              c), na.rm = TRUE) && !any(unfVect == 
                                                                        triple2numb(p, c, b, a), na.rm = TRUE)) {
                pdag[b, c] <- 1
                pdag[c, b] <- 0
              }
            }
            else {
              pdag[b, c] <- 1
              pdag[c, b] <- 0
            }
            if (verbose) 
              cat("\nRule 1':", a, "->", b, " and ", 
                  b, "-", c, " where ", a, " and ", c, 
                  " not connected and ", a, b, c, " faithful triple: ", 
                  b, "->", c, "\n")
          }
          else if (pdag[b, c] == 0 & pdag[c, b] == 1) {
            if (!is.null(unfVect)) {
              if (!any(unfVect == triple2numb(p, a, b, 
                                              c), na.rm = TRUE) && !any(unfVect == 
                                                                        triple2numb(p, c, b, a), na.rm = TRUE)) {
                pdag[b, c] <- 2
                pdag[c, b] <- 2
                if (verbose) 
                  cat("\nRule 1':", a, "->", b, "<-", 
                      c, " but ", b, "->", c, "also possible and", 
                      a, b, c, " faithful triple: ", a, 
                      "->", b, "<->", c, "\n")
              }
            }
            else {
              pdag[b, c] <- 2
              pdag[c, b] <- 2
              if (verbose) 
                cat("\nRule 1':", a, "->", b, "<-", c, 
                    " but ", b, "->", c, "also possible and", 
                    a, b, c, " faithful triple: ", a, "->", 
                    b, "<->", c, "\n")
            }
          }
        }
      }
      if (!solve.confl) 
        search.pdag <- pdag
    }
    pdag
  }
  rule2 <- function(pdag, solve.confl = FALSE) {
    search.pdag <- pdag
    ind <- which(search.pdag == 1 & t(search.pdag) == 1, 
                 arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      isC <- which(search.pdag[a, ] == 1 & search.pdag[, 
                                                       a] == 0 & search.pdag[, b] == 1 & search.pdag[b, 
                                                                                                     ] == 0)
      for (ii in seq_along(isC)) {
        c <- isC[ii]
        if (!solve.confl | (pdag[a, b] == 1 & pdag[b, 
                                                   a] == 1)) {
          pdag[a, b] <- 1
          pdag[b, a] <- 0
          if (verbose) 
            cat("\nRule 2: Chain ", a, "->", c, "->", 
                b, ":", a, "->", b, "\n")
        }
        else if (pdag[a, b] == 0 & pdag[b, a] == 1) {
          pdag[a, b] <- 2
          pdag[b, a] <- 2
          if (verbose) 
            cat("\nRule 2: Chain ", a, "->", c, "->", 
                b, ":", a, "<->", b, "\n")
        }
      }
      if (!solve.confl) 
        search.pdag <- pdag
    }
    pdag
  }
  rule3 <- function(pdag, solve.confl = FALSE, unfVect = NULL) {
    search.pdag <- pdag
    ind <- which(search.pdag == 1 & t(search.pdag) == 1, 
                 arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      c <- which(search.pdag[a, ] == 1 & search.pdag[, 
                                                     a] == 1 & search.pdag[, b] == 1 & search.pdag[b, 
                                                                                                   ] == 0)
      if (length(c) >= 2) {
        cmb.C <- combn(c, 2)
        cC1 <- cmb.C[1, ]
        cC2 <- cmb.C[2, ]
        for (j in seq_along(cC1)) {
          c1 <- cC1[j]
          c2 <- cC2[j]
          if (search.pdag[c1, c2] == 0 && search.pdag[c2, 
                                                      c1] == 0) {
            if (!is.null(unfVect)) {
              if (!any(unfVect == triple2numb(p, c1, 
                                              a, c2), na.rm = TRUE) && !any(unfVect == 
                                                                            triple2numb(p, c2, a, c1), na.rm = TRUE)) {
                if (!solve.confl | (pdag[a, b] == 1 & 
                                    pdag[b, a] == 1)) {
                  pdag[a, b] <- 1
                  pdag[b, a] <- 0
                  if (!solve.confl) 
                    search.pdag <- pdag
                  if (verbose) 
                    cat("\nRule 3':", a, c1, c2, "faithful triple: ", 
                        a, "->", b, "\n")
                  break
                }
                else if (pdag[a, b] == 0 & pdag[b, a] == 
                         1) {
                  pdag[a, b] <- pdag[b, a] <- 2
                  if (verbose) 
                    cat("\nRule 3':", a, c1, c2, "faithful triple: ", 
                        a, "<->", b, "\n")
                  break
                }
              }
            }
            else {
              if (!solve.confl | (pdag[a, b] == 1 & pdag[b, 
                                                         a] == 1)) {
                pdag[a, b] <- 1
                pdag[b, a] <- 0
                if (!solve.confl) 
                  search.pdag <- pdag
                if (verbose) 
                  cat("\nRule 3':", a, c1, c2, "faithful triple: ", 
                      a, "->", b, "\n")
                break
              }
              else if (pdag[a, b] == 0 & pdag[b, a] == 
                       1) {
                pdag[a, b] <- pdag[b, a] <- 2
                if (verbose) 
                  cat("\nRule 3':", a, c1, c2, "faithful triple: ", 
                      a, "<->", b, "\n")
                break
              }
            }
          }
        }
      }
    }
    pdag
  }
  if (numEdges(gInput@graph) == 0) 
    return(gInput)
  g <- as(gInput@graph, "matrix")
  p <- nrow(g)
  pdag <- g
  #Assumptions (Before finding the v structures)
   # pdag[5,4]<- 0
     #pdag[4,3]<- 0 
    #  pdag[2,1]<- 0
    #   pdag[3,1]<- 0
    # pdag[5,3] = 0
  #End of assumptions
    if (orientCollider) {
    ind <- which(g == 1, arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      x <- ind[i, 1]
      y <- ind[i, 2]
      allZ <- setdiff(which(g[y, ] == 1), x)
      for (z in allZ) {
        if (g[x, z] == 0 && !((y %in% gInput@sepset[[x]][[z]]) || 
                              (y %in% gInput@sepset[[z]][[x]]))) {
          if (length(unfVect) == 0) {
            if (!solve.confl) {
              pdag[x, y] <- pdag[z, y] <- 1
              pdag[y, x] <- pdag[y, z] <- 0
            }
            else {
              pdag <- orientConflictCollider(pdag, x, 
                                             y, z)
            }
          }
          else {
            if (!any(unfVect == triple2numb(p, x, y, 
                                            z), na.rm = TRUE) && !any(unfVect == triple2numb(p, 
                                                                                             z, y, x), na.rm = TRUE)) {
              if (!solve.confl) {
                pdag[x, y] <- pdag[z, y] <- 1
                pdag[y, x] <- pdag[y, z] <- 0
              }
              else {
                pdag <- orientConflictCollider(pdag, 
                                               x, y, z)
              }
            }
          }
        }
      }
    }
  }
  
  # #Assumptions (after finding the v-structures )
   # pdag[4,3] <- 0
   #  pdag[2,5] <- 0
   
   repeat {
    old_pdag <- pdag
    if (rules[1]) {
      pdag <- rule1(pdag, solve.confl = solve.confl, unfVect = unfVect)
    }
    if (rules[2]) {
      pdag <- rule2(pdag, solve.confl = solve.confl)
    }
    if (rules[3]) {
      pdag <- rule3(pdag, solve.confl = solve.confl, unfVect = unfVect)
    }
    if (all(pdag == old_pdag)) 
      break
  }
  gInput@graph <- as(pdag, "graphNEL")
  gInput
}