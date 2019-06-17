#' Divisive hierarchcal clustering
#'
#' This function recursively splits an n x p matrix into smaller and smaller subsets,
#' returning a "dendrogram" object.
#'
#' @param x a matrix
#' @param algo character string giving the partitioning algorithm to be used
#'   to split the data. Currently only "kmeans" is supported.
#' @param stand logical indicating whether the matrix should be standardised
#'   prior to the recursive partitioning procedure. Defaults to FALSE.
#' @param ... further arguments to be passed to splitting methods (not including
#'   \code{centers} if \code{algo = kmeans}).
#' @return Returns an object of class \code{"dendrogram"}.
#'
#' @details This function creates a dendrogram by successively splitting
#'   the dataset into smaller and smaller subsets (recursive
#'   partitioning). This is a divisive, or "top-down" approach to tree-building,
#'   as opposed to agglomerative "bottom-up" methods such as neighbor joining
#'   and UPGMA. It is particularly useful for large large datasets with many records
#'   (\emph{n} > 10,000) since the need to compute a large \emph{n} * \emph{n}
#'   distance matrix is circumvented.
#'
#'   If a more accurate tree is required, users can increase the value
#'   of \code{nstart} passed to \code{kmeans} \emph{via} the \code{...} argument.
#'   While this can increase computation time, it can improve accuracy
#'   considerably.
#'
#'
#' @author Shaun Wilkinson
#'
#' @references TBA
#'
#'
#' @examples
#' \dontrun{
#' ## Cluster a subsample of the iris dataset
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(999)
#' iris50 <- iris[sample(x = 1:150, size = 50, replace = FALSE),]
#' x <- as.matrix(iris50[, 1:4])
#' rownames(x) <- iris50[, 5]
#' dnd <- dclust(x, nstart = 20)
#' plot(dnd, horiz = TRUE, yaxt = "n")
#'
#' ## Color labels according to species
#' rectify_labels <- function(node, x){
#'   newlab <- factor(rownames(x))[unlist(node, use.names = FALSE)]
#'   attr(node, "label") <- newlab
#'   return(node)
#' }
#' dnd <- dendrapply(dnd, rectify_labels, x = x)
#'
#' ## Create a color palette as a data.frame with one row for each species
#' uniqspp <- as.character(unique(iris50$Species))
#' colormap <- data.frame(Species = uniqspp, color = rainbow(n = length(uniqspp)))
#' colormap[, 2] <- c("red", "blue", "green")
#'
#' ## Color the inner dendrogram edges
#' color_dendro <- function(node, colormap){
#'   if(is.leaf(node)){
#'     nodecol <- colormap$color[match(attr(node, "label"), colormap$Species)]
#'     attr(node, "nodePar") <- list(pch = NA, lab.col = nodecol)
#'     attr(node, "edgePar") <- list(col = nodecol)
#'   }else{
#'     spp <- attr(node, "label")
#'     dominantspp <- levels(spp)[which.max(tabulate(spp))]
#'     edgecol <- colormap$color[match(dominantspp, colormap$Species)]
#'     attr(node, "edgePar") <- list(col = edgecol)
#'   }
#'   return(node)
#' }
#' dnd <- dendrapply(dnd, color_dendro, colormap = colormap)
#'
#' ## Plot the dendrogram
#' plot(dnd, horiz = TRUE, yaxt = "n")
#' }
################################################################################
dclust <- function(x, algo = "kmeans", stand = FALSE, ...){
  x <- as.matrix(x)
  stopifnot(is.matrix(x))
  if(stand) x <- scale(x)
  nrec <- nrow(x)
  if(nrec == 1){# singleton tree (leaf)
    tree <- 1
    attr(tree, "leaf") <- TRUE
    attr(tree, "height") <- 0
    attr(tree, "midpoint") <- 0
    attr(tree, "label") <- names(x)[1]
    attr(tree, "members") <- 1
    class(tree) <- "dendrogram"
    return(tree)
  }
  if(is.null(rownames(x))) rownames(x) <- paste0("RECORD", 1:nrec)
  catchnames <- rownames(x)
  charstrings <- apply(x, 1, paste0, collapse = "")
  hashes <- openssl::md5(charstrings)
  duplicates <- duplicated(hashes)
  nurec <- sum(!duplicates)
  point <- function(h){
    uh <- unique(h)
    pointers <- seq_along(uh)
    names(pointers) <- uh
    unname(pointers[h])
  }
  pointers <- point(hashes)
  x <- x[!duplicates, ]
  if(sum(!duplicates) == 1){
    tree <- vector(mode = "list", length = nrow(x))
    attr(tree, "height") <- 0
    attr(tree, "midpoint") <- 0.5
    attr(tree, "members") <- nrow(x)
    for(i in seq_len(nrow(x))){
      tree[[i]] <- i
      attr(tree[[i]], "height") <- 0
      attr(tree[[i]], "midpoint") <- 0
      attr(tree[[i]], "members") <- 1
      attr(tree[[i]], "leaf") <- TRUE
      attr(tree[[i]], "label") <- catchnames[i]
    }
    class(tree) <- "dendrogram"
    return(tree)
  }
  #kcounts <- kcount(x, k = k, residues = residues, gap = gap, named = FALSE)
  kcounts <- x
  tree <- 1
  attr(tree, "leaf") <- TRUE
  attr(tree, "records") <- 1:nurec
  attr(tree, "height") <- 0
  #attr(tree, "kvector") <- apply(kcounts, 2, mean)
  ## define recursive splitting functions
  clustern <- function(node, kcs, ...){
    if(!is.list(node) & length(attr(node, "records")) > 1){
      ## fork leaves only
      recs <- kcs[attr(node, "records"), , drop = FALSE]
      errfun <- function(er){## used when >3 uniq hashes but kmeans throws error
        cat("yikes!\n")
        out <- list()
        nrs <- nrow(recs)
        cls <- rep(1, nrs)
        cls[sample(1:nrs, 1)] <- 2 ## peel randomly selected one off
        out$cluster <- cls
        out$centers <- rbind(apply(recs[cls == 1, , drop = FALSE], 2, mean),
                             apply(recs[cls == 2, , drop = FALSE], 2, mean))
        out$totss <- 0.0001########################3 arbitrary small value
        return(out)
      }
      km <- if(nrow(recs) > 2){
        tryCatch(kmeans(recs, centers = 2, ... = ...),
                 error = errfun, warning = errfun)
      }else{
        ## totss = sum(apply(Z, 2, var) * (nrow(Z) - 1))
        list(cluster = 1:2, centers = recs, totss = sum(apply(recs, 2, var)))
      }
      tmpattr <- attributes(node)
      node <- vector(mode = "list", length = 2)
      attributes(node) <- tmpattr
      attr(node, "leaf") <- NULL
      attr(node, "height") <- km$totss
      for(i in 1:2){
        node[[i]] <- 1

        # attr(node[[i]], "kvector") <- km$centers[i, ]
        # kmatrix <- rbind(attr(node, "kvector"), attr(node[[i]], "kvector"))
        # rownames(kmatrix) <- paste(1:2)
        # diffheight <- sqrt(sum((apply(kmatrix, 2, function(v) v[1] - v[2]))^2))
        # attr(node[[i]], "height") <- attr(node, "height") - diffheight ## cleaned up later

        attr(node[[i]], "height") <- 0 ## replaced later if not a leaf
        attr(node[[i]], "leaf") <- TRUE
        attr(node[[i]], "records") <- attr(node, "records")[km$cluster == i]

      }
    }
    return(node)
  }
  clusterr <- function(tree, kcs, ...){ # kcs is the kfreq matrix
    tree <- clustern(tree, kcs, ...)
    if(is.list(tree)) tree[] <- lapply(tree, clusterr, kcs = kcs, ...)
    return(tree)
  }
  ##  build tree recursively
  tree <- clusterr(tree, kcs = kcounts, ... = ...)
  tree <- phylogram::remidpoint(tree)
  class(tree) <- "dendrogram"
  tree <- phylogram::reposition(tree)
  if(any(duplicates)){
    reduplicate <- function(node, pointers){
      attr(node, "records") <- which(pointers %in% attr(node, "records"))
      if(is.leaf(node)){
        lams <- length(attr(node, "records"))
        if(lams > 1){
          labs <- attr(node, "label")
          hght <- attr(node, "height")
          recs <- attr(node, "records")
          node <- vector(mode = "list", length = lams)
          attr(node, "height") <- hght
          attr(node, "records") <- recs
          for(i in 1:lams){
            node[[i]] <- 1
            attr(node[[i]], "height") <- hght
            attr(node[[i]], "label") <- labs[i]
            attr(node[[i]], "records") <- recs[i]
            attr(node[[i]], "leaf") <- TRUE
          }
        }
      }
      return(node)
    }
    tree <- dendrapply(tree, reduplicate, pointers)
    tree <- phylogram::remidpoint(tree)
  }
  label <- function(node, labs){
    if(is.leaf(node)) attr(node, "label") <- labs[attr(node, "records")]
    return(node)
  }
  tree <- dendrapply(tree, label, labs = catchnames)
  rmrecs <- function(node){
    if(is.leaf(node)){
      tmpattr <- attributes(node)
      node[] <- tmpattr$records
      tmpattr$records <- NULL
      #tmpattr$kvector <- NULL
      attributes(node) <- tmpattr
    }else{
      attr(node, "records") <- NULL
      #attr(node, "kvector") <- NULL
    }
    return(node)
  }
  tree <- dendrapply(tree, rmrecs)
  return(tree)
}
