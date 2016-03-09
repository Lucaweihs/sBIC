toEdgeListMat = function(edgeListVec) {
  if (is.matrix(edgeListVec)) {
    return(edgeListVec)
  }
  return(t(matrix(edgeListVec, 2, length(edgeListVec) / 2)))
}

toEdgeListVec = function(edgeListMat) {
  if (is.vector(edgeListMat)) {
    return(edgeListMat)
  }
  return(as.numeric(t(edgeListMat)))
}

getSubModelSupports = function(m, edges) {
  # This function computes possible getSubModelSupports (q-subforests)
  # the output is in the {0,1}-vectors indicating which edges are in.
  # m = number of leaves
  # edges = |E|x2 matrix listing the edges
  edges = toEdgeListMat(edges)

  ne = nrow(edges)
  if (ne == 0) {
    return(matrix(numeric(0), ncol = 0, nrow = 1))
  }

  subs = combinat::hcube(rep(2, ne), 1,-1) # indicators of all possible subsets of E
  subms = rep(1, 2 ^ ne)

  # Compute degree of each node
  nodeDegree = numeric(max(edges))
  for (i in 1:nrow(edges)) {
    for (j in 1:ncol(edges)) {
      nodeDegree[edges[i, j]] = nodeDegree[edges[i, j]] + 1
    }
  }

  maxVertex = max(edges)
  if (maxVertex > m) {
    for (i in 1:2 ^ ne) {
      indicatorMatrix = matrix(c(subs[i,], subs[i,]), ne, 2)
      for (j in (m + 1):max(edges)) {
        if (sum(indicatorMatrix * edges == j) == 1) {
          subms[i] = 0
        }
      }
    }
  }

  if (sum(subms) == 1) {
    return(t(subs[subms == 1, ]))
  } else if (ne == 1) {
    return(matrix(subs[subms == 1, ]))
  } else {
    return(subs[subms == 1, ])
  }
}

getSubmodelEdges = function(subs, edges) {
  #returns the set of edges for a given indicator function
  return(toEdgeListMat(edges)[subs == 1, , drop = F])
}

subModelsToDAG = function(subModels) {
  numSubModels = nrow(subModels)
  if (ncol(subModels) == 0) {
    return(graph.empty(1))
  }

  lessThanMatrix = matrix(numeric(numSubModels ^ 2), ncol = numSubModels)

  for (i in 1:(numSubModels - 1)) {
    curSupport = subModels[i, ]
    curNumEdges = sum(curSupport)
    for (j in (i + 1):numSubModels) {
      compSupport = subModels[j, ]
      compNumEdges = sum(compSupport)
      intersectSum = sum(curSupport * compSupport)
      if (intersectSum == curNumEdges) {
        lessThanMatrix[i, j] = 1
      } else if (intersectSum == compNumEdges) {
        lessThanMatrix[j, i] = 1
      }
    }
  }

  greaterThanMatrix = 1 - lessThanMatrix
  trueLessThanMatrix = lessThanMatrix
  for (i in 1:numSubModels) {
    for (j in which(lessThanMatrix[i, ] == 1)) {
      trueLessThanMatrix[i, ] = trueLessThanMatrix[i, ] * greaterThanMatrix[j, ]
    }
  }

  numEdgesInLattice = sum(trueLessThanMatrix)
  edgeMatrix = matrix(numeric(numEdgesInLattice * 2), ncol = 2)

  k = 1
  for (i in 1:numSubModels) {
    for (j in which(trueLessThanMatrix[i, ] == 1)) {
      edgeMatrix[k, ] = c(i, j)
      k = k + 1
    }
  }

  return(graph.edgelist(toEdgeListMat(edgeMatrix)))
}

forestDim = function(edgeList, support, numLeaves) {
  if (numLeaves == 0) {
    return(0)
  }
  edgeList = toEdgeListMat(edgeList[support == 1,])
  g = graph.edgelist(edgeList, directed = F)
  return(numLeaves + nrow(edgeList) - sum(degree(g) == 2))
}

forestFromTreeEdgesAndSupport = function(tree, edges, support) {
  if (is.vector(edges)) {
    edges = t(matrix(edges, 2, length(edges) / 2))
  }
  edgesToDelete = as.vector(t(edges[support == 0, ]))
  return(delete.edges(tree, get.edge.ids(tree, edgesToDelete)))
}

edgesInRootedDAGOrder = function(edgeList, root = "max") {
  if (is.vector(edgeList)) {
    edgeList = t(matrix(edgeList, 2, length(edgeList) / 2))
  }
  if (nrow(edgeList) == 0) {
    return(edgeList)
  }
  if (root == "max") {
    root = max(edgeList)
  }

  g = graph.empty(max(edgeList), directed = F)
  g = add.edges(g, as.vector(t(edgeList)))
  dfsResults = graph.dfs(g, max(edgeList), father = T, unreachable = T)
  dfsFather = as.numeric(dfsResults$father)

  for (i in 1:nrow(edgeList)) {
    v1 = edgeList[i, 1]
    v2 = edgeList[i, 2]
    if (!is.na(dfsFather[v1]) && dfsFather[v1] == v2) {
      edgeList[i,] = c(v2, v1)
    }
  }

  return(edgeList)
}

removeDeg2Nodes = function(forest) {
  deg2Nodes = which(degree(forest) == 2)

  while (length(deg2Nodes) > 0) {
    deg2Node = deg2Nodes[1]
    adjNodes = as.numeric(get.adjlist(forest)[[deg2Node]])
    edgeIdsToDelete = as.numeric(get.edge.ids(forest, c(deg2Node, adjNodes[1], deg2Node, adjNodes[2])))
    forest = delete.edges(forest, edgeIdsToDelete)
    forest = add.edges(forest, adjNodes)
    deg2Nodes = which(degree(forest) == 2)
  }
  return(forest)
}

###
# Takes as input two vectors (pBig, pSmall) of the same length n which
# represent partitions of 1,...,n. Here i,j of pBig (pSmall) are
# considered to be in the same partition if pBig[i] == pBig[j]
# (pSmall[i] == pSmall[j]). The function returns true if pSmall is a
# refinement (sub-partition) of bBig.
###
isSubPartition = function(pBig, pSmall) {
  bigInds = sort(unique(pBig))
  for (bigInd in bigInds) {
    smallClustInds = sort(unique(pSmall[pBig == bigInd]))
    if (any(smallClustInds %in% pSmall[pBig != bigInd])) {
      return(F)
    }
  }
  return(T)
}

###
# Requires that either f1 or f2 is a sub-forest of the other and
# that the leaves are the first numLeaves vertices of the two forests.
###
forestDistance = function(f1, f2, numLeaves) {
  if (length(V(f1)) == 0) {
    return(0)
  }
  leaves = 1:numLeaves
  f1 = removeDeg2Nodes(f1)
  f2 = removeDeg2Nodes(f2)
  f1Clusters = clusters(f1)$membership[leaves]
  f2Clusters = clusters(f2)$membership[leaves]

  # There are fewer, bigger, clusters in one forest
  forestFlipped = F
  if (max(f1Clusters) <= max(f2Clusters)) {
    fBig = f1
    bigClusters = f1Clusters
    fSmall = f2
    smallClusters = f2Clusters
  } else {
    forestFlipped = T
    fBig = f2
    bigClusters = f2Clusters
    fSmall = f1
    smallClusters = f1Clusters
  }

  if (!isSubPartition(bigClusters, smallClusters)) {
    return(NA)
  }

  numBigClusters = length(unique(bigClusters))
  numSmallClusters = length(unique(smallClusters))
  dist = numSmallClusters - numBigClusters
  while (numBigClusters != numSmallClusters) {
    edgeList = get.edgelist(fBig)
    didChange = F
    for (i in 1:nrow(edgeList)) {
      edge = edgeList[i, ]
      edgeId = get.edge.ids(fBig, edge)
      newF = removeDeg2Nodes(delete.edges(fBig, edgeId))
      newClusters = clusters(newF)$membership[leaves]
      if (isSubPartition(newClusters, smallClusters)) {
        didChange = T
        fBig = newF
        bigClusters = newClusters
        numBigClusters = length(unique(bigClusters))
        break
      }
    }
    if (!didChange) {
      return(NA)
    }
  }

  bigForestString = forestToString(fBig, numLeaves)
  smallForestString = forestToString(fSmall, numLeaves)
  if (bigForestString == smallForestString) {
    return(dist * (-1) ^ forestFlipped)
  } else {
    return(NA)
  }
}

###
# Forest with no edges is considered to be at depth 0 while a full
# tree has depth numLeaves-1
###
forestDepth = function(forest, numLeaves) {
  numClusters = length(unique(clusters(forest)$membership[1:numLeaves]))
  return((numLeaves - 1) - (numClusters - 1))
}


# Helper function
isRootedDag = function(g) {
  if (!is.simple(g)) {
    return(F)
  }
  clusts = clusters(g, mode = "weak")
  membership = clusts$membership
  numClusts = clusts$no

  for (i in 1:numClusts) {
    sourceFound = F
    vertices = which(membership == i)
    for (v in vertices) {
      if (length(incident(g, v, "in")) == 0) {
        if (sourceFound == T) {
          return(F)
        }
        sourceFound = T
      }
    }
  }
  return(T)
}

is.support = function(support) {
  return(is.vector(support) && is.atomic(support)
         &&
           (length(support) == 0 || all(support == 1 || support == 0)))
}

is.atomic.vec = function(vec) {
  return(is.atomic(vec) && is.vector(vec))
}

is.edgelist.mat = function(edgeList) {
  return(is.matrix(edgeList) &&
           (ncol(edgeList) == 2) &&
           all(edgeList %% 1 == 0) &&
           all(edgeList > 0))
}

is.edgelist.vec = function(edgeList) {
  return(
    is.vector(edgeList) &&
      is.atomic(edgeList) &&
      (length(edgeList) %% 2 == 0) &&
      all(edgeList %% 1 == 0) &&
      all(edgeList > 0)
  )
}

is.edgelist.vec.or.mat = function(edgeList) {
  return(((is.matrix(edgeList) &&
             ncol(edgeList) == 2) ||
            (
              is.vector(edgeList) &&
                is.atomic(edgeList) &&
                length(edgeList) %% 2 == 0
            )
  ) &&
    all(edgeList %% 1 == 0) &&
    all(edgeList > 0))
}

isBinaryEdgelistMat = function(edgeList, numLeaves) {
  if (!is.edgelist.mat(edgeList)) {
    return(F)
  }
  g = graph.edgelist(edgeList, directed = F)
  v = length(V(g))
  if (v < numLeaves) {
    g = add.vertices(g, numLeaves - v)
  }
  return(is.binary.forest(g, numLeaves))
}

isBinaryEdgelistVecOrMat = function(edgeList, numLeaves) {
  return(isBinaryEdgelistMat(toEdgeListMat(edgeList), numLeaves))
}

is.pos.int = function(num) {
  return(length(num) == 1 &&
           (is.numeric(num) || is.integer(num)) && (num %% 1 == 0 && num > 0))
}

is.nonneg.int = function(num) {
  return(length(num) == 1 &&
           (is.numeric(num) ||
              is.integer(num)) && (num %% 1 == 0 && num >= 0))
}

is.data.matrix = function(mat) {
  return(is.matrix(mat) && nrow(mat) > 0)
}

is.forest = function(f) {
  minSpanTree = minimum.spanning.tree(f)
  return(!is.directed(f) && length(E(f)) == length(E(minSpanTree)))
}

is.binary.forest = function(f, numLeaves) {
  if (numLeaves == 0) {
    return(length(V(f)) == 0)
  }
  degs = degree(f)
  return(
    is.forest(f) &&
      numLeaves %% 1 == 0 &&
      numLeaves > 0 &&
      numLeaves <= length(V(f)) &&
      all(degs[1:numLeaves] %in% c(0, 1)) &&
      all(degs[-(1:numLeaves)] %in% c(0, 2, 3))
  )
}

###
# Creates a unique string representation of the bifucating forest.
##
forestToString = function(forest, numLeaves) {
  subTreeMembership = clusters(forest)$membership

  numSubTrees = max(subTreeMembership)
  subTreeVertexList = vector("list", numSubTrees)

  subTreeOrder = numeric(numSubTrees)
  subTreeContainsLeaf = logical(numSubTrees)
  for (i in 1:numSubTrees) {
    subTreeVertexList[[i]] = sort(which(subTreeMembership == i))
    subTreeOrder[i] = subTreeVertexList[[i]][1]
    if (subTreeOrder[i] <= numLeaves) {
      subTreeContainsLeaf[i] = T
    }
  }
  subTreePermutation = order(subTreeOrder)

  newSubTreeVertexList = vector("list", sum(subTreeContainsLeaf))
  k = 1
  for (i in subTreePermutation) {
    if (subTreeContainsLeaf[i] == T) {
      newSubTreeVertexList[[k]] = subTreeVertexList[[i]]
      k = k + 1
    }
  }
  subTreeVertexList = newSubTreeVertexList
  numSubTrees = length(subTreeVertexList)

  adjList = get.adjlist(forest, "all")
  ancestors = vector("list", length(V(forest)))

  subTreeStrings = character(numSubTrees)
  for (i in 1:numSubTrees) {
    subTreeVertices = subTreeVertexList[[i]]

    if (length(subTreeVertices) == 1) {
      subTreeStrings[i] = paste("(", subTreeVertices[1], ")", sep = "")
    } else if (sum(subTreeVertices <= numLeaves) == 2) {
      subTreeLeaves = sort(subTreeVertices[which(subTreeVertices <= numLeaves)])
      subTreeStrings[i] = paste("((", subTreeLeaves[1], ")(", subTreeLeaves[2], "))", sep =
                                  "")
    } else {
      # The smallest node in the subtree, note that this
      # is the first in subTreeVertices by the above ordering
      # We know that this smallest node is a leaf by the above pruning.
      smallestObserved = subTreeVertices[1]

      smallestObservedDegree3Parent = adjList[[smallestObserved]]

      last = smallestObserved
      while (length(adjList[[smallestObservedDegree3Parent]]) != 3) {
        curNeighbors = adjList[[smallestObservedDegree3Parent]]
        current = smallestObservedDegree3Parent
        smallestObservedDegree3Parent = curNeighbors[which(curNeighbors != last)]
        last = current
      }

      dfsResult = graph.dfs(
        forest,
        root = smallestObservedDegree3Parent,
        unreachable = F,
        father = T,
        order.out = T
      )

      orderOut = as.numeric(dfsResult$order.out)
      order = as.numeric(dfsResult$order)
      parents = as.numeric(dfsResult$father)

      orderOut = orderOut[!is.na(orderOut)]
      order = order[!is.na(order)]

      for (node in orderOut) {
        if (!is.na(parents[node]) && parents[node] != 0) {
          parent = parents[node]
          ancestors[[parent]] = sort(c(ancestors[[parent]], ancestors[[node]], node))
        }
      }

      subTreeStrings[i] = treeToString(smallestObservedDegree3Parent,
                                       parents,
                                       ancestors,
                                       adjList,
                                       numLeaves)
    }
  }

  return(paste(subTreeStrings, collapse = "|"))
}

###
# Helper function, do not call.
###
getSubForestString = function(forest, edgeList, support, numLeaves) {
  edgesToDelete = edgeList[support == 0, ]
  if (!is.vector(edgesToDelete)) {
    edgesToDelete = as.vector(t(edgesToDelete))
  }
  edgesToDelete = get.edge.ids(forest, edgesToDelete)
  return(forestToString(delete.edges(forest, edgesToDelete), numLeaves))
}

###
# Helper function, do not call.
###
treeToString = function(root,
                        parents,
                        ancestors,
                        adjList,
                        numLeaves) {
  rootNeighbors = adjList[[root]]
  if (length(rootNeighbors) == 1 || length(rootNeighbors) == 0) {
    if (root <= numLeaves) {
      return(paste("(", as.character(root), ")", sep = ""))
    } else {
      return("")
    }
  } else if (length(rootNeighbors) == 2) {
    print("ERROR: root should never have degree 2.")
    return(F)
  } else if (length(rootNeighbors) == 3) {
    if (!is.na(parents[root]) && parents[root] != 0) {
      rootNeighbors = rootNeighbors[rootNeighbors != parents[root]]
    }

    neighborPermutation = numeric(length(rootNeighbors))
    for (i in 1:length(rootNeighbors)) {
      neighborPermutation[i] = min(rootNeighbors[i], ancestors[[rootNeighbors[i]]][1])
    }
    rootNeighbors = rootNeighbors[order(neighborPermutation)]

    treeString = character(length(rootNeighbors))
    for (i in 1:length(rootNeighbors)) {
      newRoot = rootNeighbors[i]
      while (length(adjList[[newRoot]]) == 2) {
        newRoot = adjList[[newRoot]][which(adjList[[newRoot]] != parents[newRoot])]
      }
      treeString[i] = treeToString(newRoot, parents, ancestors, adjList, numLeaves)
    }
    treeString = treeString[treeString != ""]

    if (length(treeString) == 0) {
      return("")
    } else if (length(treeString) == 1) {
      return(treeString[1])
    } else {
      return(paste("(", paste(treeString, collapse = ""), ")", sep = ""))
    }
  } else {
    print("Error: input tree is not bivariate.")
    return(F)
  }
}

generateRandomBinaryTree = function(numLeaves, directedOrdering = F) {
  if (numLeaves %in% c(0, 1)) {
    return(matrix(numeric(0), ncol = 2))
  }

  edgeMatrix = t(c(1, 2))
  n = 2
  while (n < numLeaves) {
    maxVertex = max(edgeMatrix)
    edgeMatrix = edgeMatrix + (edgeMatrix > n)

    edgeToSplit = sample(1:(dim(edgeMatrix)[1]), 1)

    newEdges = matrix(
      c(edgeMatrix[edgeToSplit, 1], maxVertex + 2, edgeMatrix[edgeToSplit, 2], maxVertex + 2),
      ncol = 2,
      byrow = T
    )
    newEdgeMatrix = edgeMatrix
    newEdgeMatrix[edgeToSplit, ] = c(maxVertex + 2, n + 1)
    newEdgeMatrix = rbind(newEdgeMatrix, newEdges)

    edgeMatrix = newEdgeMatrix
    n = n + 1
  }

  if (directedOrdering) {
    g = graph.edgelist(edgeMatrix, directed = F)
    dfsResults = graph.dfs(g, max(edgeMatrix), father = T)

    edgeMatrix = edgeMatrix * 0

    for (j in 1:(length(dfsResults$order) - 1)) {
      node = dfsResults$order[j + 1]
      father = dfsResults$father[node]
      edgeMatrix[j, 1] = father
      edgeMatrix[j, 2] = node
    }
  }

  return(edgeMatrix)
}

generateRandomBinaryForest = function(numLeaves, depth) {
  if (numLeaves < 0 || depth < 0 || depth > numLeaves - 1) {
    stop("Invalid number of leaves or depth")
  }
  if (numLeaves == 1) {
    return(graph.empty(1, directed = F))
  }
  edgeList = generateRandomBinaryTree(numLeaves, T)

  numVertices = 2 * numLeaves - 2
  support = rep(1, nrow(edgeList))
  subForestIndicator = support

  if ((numLeaves - 1) - depth > 0) {
    for (i in 1:((numLeaves - 1) - depth)) {
      supports = generateOneStepSubForestSupports(subForestIndicator, edgeList, numLeaves, numVertices)
      subForestIndicator = supports[sample(nrow(supports), 1),]
    }
  }

  tree = graph.edgelist(edgeList, directed = F)
  edgesToDelete = edgeList[subForestIndicator == 0, ]
  if (!is.vector(edgesToDelete)) {
    edgesToDelete = as.vector(t(edgesToDelete))
  }
  edgesToDelete = get.edge.ids(tree, edgesToDelete)
  return(removeDeg2Nodes(delete.edges(tree, edgesToDelete)))
}

###
# Generates the subforests F' that are strictly less than the given
# forest F and for which there exists no other subforest Q such that
# F' subset Q subset F (deg 2 nodes have been removed).
#
# Assumes the support comes from a valid minimal subforest F.
###
generateOneStepSubForestSupports = function(support, edgeList, numLeaves, numVertices) {
  edgeList = toEdgeListMat(edgeList)

  numPossibleSubForests = sum(support)
  if (numPossibleSubForests == 0) {
    return(numeric(0))
  }

  subForestSupports = matrix(numeric(numPossibleSubForests * length(support)), ncol =
                               length(support))
  edgesToRemove = which(support != 0)

  seenSupports = hash()
  k = 1
  for (i in 1:numPossibleSubForests) {
    newSupport = support
    newSupport[edgesToRemove[i]] = 0
    newSupport = pruneEdges(newSupport, edgeList, numLeaves, numVertices)
    newSupportString = paste(newSupport, collapse = "")
    if (!has.key(newSupportString, seenSupports)) {
      seenSupports[[newSupportString]] = 1
      subForestSupports[k,] = newSupport
      k = k + 1
    }
  }

  if (k - 1 == 1) {
    return(t(subForestSupports[1:(k - 1),]))
  } else {
    return(subForestSupports[1:(k - 1),])
  }
}

pruneEdges = function(support, edgeList, numLeaves, numVertices) {
  if (nrow(edgeList) == 0) {
    return(numeric(0))
  }
  return(pruneEdgesHelper(support, edgeList, numLeaves, numVertices))
}