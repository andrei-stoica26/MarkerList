#' @importClassesFrom Seurat Seurat
#' @importFrom MatrixGenerics rowMaxs
#' @importFrom qs2 qs_save
#' @importFrom Seurat FindMarkers
#'
NULL

#' Order markers decreasingly by ratio
#'
#' This function calculates pct.ratio (pct.1 / pct.2) for a data frame generated
#' using FindMarkers, then used it to order the data frame decreasingly.
#'
#' @param markers A data frame generated using FindMarkers
#'
#' @return A data frame of markers ordered decreasingly by pct.ratio (pct.1 / pct.2)
#'
#' @export
#'
ratioOrder <- function(markers){
  markers$pct.ratio <- markers$pct.1 / markers$pct.2
  markers <- markers[order(markers$pct.ratio, decreasing=TRUE), ]
  return(markers)
}

#' Generate markers of Seurat identity classes and order them decreasingly
#'
#' This function is a wrapper around FindMarkers from Seurat that relaxes the
#' criteria for marker selection as default option, calculates pct.ratio (pct. 1
#'  / pct.2), and used it to order decreasingly the data frame.
#'
#' @param seuratObj A Seurat object
#' @inheritParams Seurat::FindMarkers
#'
#' @return A data frame of markers ordered decreasingly by pct.ratio (pct.1 / pct.2)
#'
#' @export
#'
orderedFindMarkers <- function(seuratObj, min.pct = 0, logfc.threshold = 0, min.diff.pct = 0, only.pos = T, densify = T,
                               min.cells.group = 1, ...){
  markers <- FindMarkers(seuratObj, min.pct=min.pct, logfc.threshold=logfc.threshold, min.diff.pct=min.diff.pct,
                         only.pos=only.pos, densify=densify, min.cells.group = min.cells.group, ...)
  if(nrow(markers))
  {
    markers <- subset(markers, p_val_adj < 0.05)
    markers <- ratioOrder(markers)
  }
  return(markers)
}

#' Generate markers of Seurat identity classes and order them decreasingly
#'
#' This function is a wrapper around FindMarkers from Seurat that relaxes the
#' criteria for marker selection as default option, calculates pct.ratio (pct. 1
#'  / pct.2), and used it to order decreasingly the data frame.
#'
#' @param seuratObj A Seurat object
#' @inheritParams Seurat::FindMarkers
#' @param objectName Character to be used in the file name when saving the markers.
#' Default is NULL (the markers will not be saved)
#'
#' @export
#'
#' @return A data frame of markers ordered decreasingly by pct.ratio (pct.1 / pct.2)
#'
allMarkers <- function(seuratObj, objectName = NULL, group.by = 'seurat_clusters', ...){
  allMarkers <- lapply(levels(seuratObj), function(x) {
    message(paste0('Finding markers for identity class: ', x, '...'))
    return(orderedFindMarkers(seuratObj, ident.1 = x, ...))
  })
  if (!is.null(objectName))
    qs_save(allMarkers, paste0(objectName, 'AllMarkers.qs'))
  return(allMarkers)
}


#' Rank markers of Seurat identity classes
#'
#' This function ranks markers generated using FindMarkers or orderedFindMarkers
#' using input criteria, and filters top markers using a percentage argument.
#'
#' @param markers A data frame generated using FindMarkers or orderedFindMarkers
#' The latter is required if pct.ratio is listed among the criteria
#' @param criteria Columns of the input data frame that will be involved in the
#' ranking. Must be character objects
#' @param signs A vector of the same lengths as criteria, in which -1 and 1 are
#' the only allowed values; -1 indicates that the corresponding criterion will
#' be used to rank the markers decreasingly, while 1 indicates that it will be
#' used to rank the markers increasingly
#' Default is NULL (the markers will not be saved)
#' @param percMarkers Percentage of top markers retained. If NULL, all markers
#' will be retained
#' @param rankFun The function used to rank markers
#'
#' @return A data frame of markers ordered decreasingly by pct.ratio (pct.1 / pct.2)
#'
#' @export
#'
rankMarkers <- function(markers, criteria = c('avg_log2FC', 'pct.1', 'p_val_adj', 'pct.ratio'), signs = c(-1, -1, 1, -1),
                        percMarkers = 2, rankFun = MatrixGenerics::rowMaxs){
  for (i in 1:length(criteria))
    markers[, paste0(criteria[i], '_rank')] <- rank(signs[i] * markers[, criteria[i]], ties.method = 'min')
  markers$rank <- rankFun(as.matrix(markers[, (ncol(markers) - length(criteria) + 1):ncol(markers)]))
  markers <- markers[order(markers$rank), ]
  if (!is.null(percMarkers))
    markers <- markers[1:round(percMarkers / 100 * nrow(markers)), ]
  return(markers)
}

#' Generate list of marker names from list of marker data frames
#'
#' This function takes a list of marker data frames generated using FindMarkers
#' or orderedFindMarkers, ranks and filters the markers in each data frame,
#' and returns a list of marker names.
#'
#' @param markerList A data frame generated using FindMarkers or orderedFindMarkers
#' The latter is required if pct.ratio is listed among the criteria
#' @param removeRep If true, remove the markers appearing in more than one data
#' frame. Default is FALSE
#' @inheritParams rankMarkers
#' @param nMarkersThr The lowest acceptable number of markers in a data frame;
#' if a data frame has fewer marker than this number, its markers will be
#' excluded from the output
#'
#' @return A data frame of markers ordered decreasingly by pct.ratio (pct.1 / pct.2)
#'
#' @export
#'
markerListNames <- function(markerList, removeRep = FALSE, percMarkers = 2, nMarkersThr = 5){
  if(removeRep)
    markerList <- removeRepeatedMarkers(markerList)
  names(markerList) <- 0:(length(markerList) - 1)
  res <- lapply(markerList, function(x) rankMarkers(x, percMarkers=percMarkers))
  res <- res[which(sapply(res, function(x) nrow(x) >= nMarkersThr))]
  res <- lapply(res, rownames)
  return(res)
}

#' Remove non-exclusive markers
#'
#' This function removes all markers found in more than one data frame from an
#' input list
#'
#' @inheritParams markerListNames
#'
#' @return A list of marker data frames containing only exclusive markers
#'
#' @export
#'
removeRepeatedMarkers <- function(markerList){
  df <- data.frame(table(unlist(sapply(markerList, rownames))))
  df <- subset(df, Freq == 1)
  exclusiveMarkers <- df$Var1
  res <- lapply(markerList, function(x) x[intersect(rownames(x), exclusiveMarkers), ])
  return (res)
}
