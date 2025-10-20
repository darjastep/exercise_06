rm(list=ls())
if (!requireNamespace("gtools", quietly = TRUE)) {
  install.packages("gtools")
}
library(gtools)

#task 1
DoubleDigestProblem <- function(A_fragments, B_fragments, AB_fragments) {
  
  # sorting
  AB_fragments <- sort(AB_fragments)
  total_length <- sum(AB_fragments)
  
  
  A_perms <- permutations(n = length(A_fragments), r = length(A_fragments), v = A_fragments)
  B_perms <- permutations(n = length(B_fragments), r = length(B_fragments), v = B_fragments)
  #prazdny list kam se budou ukladat vysledky
  valid_results <- list()
  
  for (i in 1:nrow(A_perms)) {
    A_perm <- A_perms[i,] # vyber konkretni usporadani 
    A_cuts <- cumsum(A_perm) # cum sum co nam najde "Cuts"
    A_sites <- c(0, A_cuts) # prida nulu na zacatek, mame sites, 
    
    # udelame to same pro B
    for (j in 1:nrow(B_perms)) {
      B_perm <- B_perms[j, ]
      B_cuts <- cumsum(B_perm)
      B_sites <- c(0, B_cuts)
      
      # kontrola delek A a B..
      #tail(A_sites, 1) vrati posledni cislo , to posledni cislo znamena celkovy pocet BP v danem DNA
      if (tail(A_sites, 1) != total_length || tail(B_sites, 1) != total_length) {
        next
      }
      
      
      
      # spojeni "strihu" (cuts) do jednoho vektoru , sort je seradi podle velikosti :) +unique odstrani duplicitni strihy
      all_sites <- sort(unique(c(A_sites, B_sites)))
      
      # vypocet rozdilu mezi sousednimi strihy - diff(c(0, 2, 3, 6))
      #                                          => c(2, 1, 3)
      combined_fragments <- diff(all_sites)
      
      #je potreba dat sort na comb_frag, jelikoz to se jeste nesortovalo
      #porovnani jestli sedi nase combined fragments s AB???:
      if (identical(sort(combined_fragments), AB_fragments)) {
        result <- list(
          A_sites = A_sites[-c(1, length(A_sites))],
          B_sites = B_sites[-c(1, length(B_sites))],
          combined_fragments = combined_fragments
        )
        #pridani results do listu
        valid_results[[length(valid_results) + 1]] <- result
      }
  
    }
   }
    return(valid_results)
  }
    
A <- c(3,5,2,10)
B <- c(7,3,10)
AB <- c(1, 2, 2, 5, 5, 5)

DoubleDigestProblem(A,B,AB)


#task 2
#pomocna funkce delta:
Delta <- function(y, X) {
  abs(y - X)  # vrací vektor vzdáleností mezi y a každým prvkem X
}

#pomocna funkce remove
Remove <- function(deltaX, diffs) {
  for (d in diffs) {
    i <- match(d, deltaX)
    if (!is.na(i)) {
      deltaX <- deltaX[-i]
    } else {
      return(NULL)
    }
  }
  return(deltaX)
}

#funkce place
Place <- function(deltaX, X, width) {
  if (length(deltaX) ==0){
    return(X)
  }
  y <- max(deltaX)
  # první pokus: přidat y jako řez
  deltas1 <- Delta(y, X)
  if (all(sapply(unique(deltas1), function(x) sum(deltas1 == x) <= sum(deltaX == x)))) {
    new_deltaX <- Remove(deltaX, deltas1)
    if (!is.null(new_deltaX)) {
      X_new <- c(X, y)
      result <- Place(new_deltaX, X_new, width)
      if (!is.null(result)) return(result)
    }
  }
  # druhý pokus: přidat (width - y)
  y2 <- width - y
  deltas2 <- Delta(y2, X)
  if (all(sapply(unique(deltas2), function(x) sum(deltas2 == x) <= sum(deltaX == x)))) {
    new_deltaX2 <- Remove(deltaX, deltas2)
    if (!is.null(new_deltaX2)) {
      X_new <- c(X, y2)
      result <- Place(new_deltaX2, X_new, width)
      if (!is.null(result)) return(result)
    }
  }
  return(NULL) #zadne platne reseni
}

deltaX <- c(2, 2, 3, 3, 4, 5, 6, 7, 8, 10)
width <- 10
X <- c(0, width)

#kontrola
solution <- Place(deltaX, X, width)
