library(foreign)
root = "G:/gradient"
outShedEdgesFile = "C:/Users/ruesca/Dropbox/gradient/shedEdges.RData"
# outShedEdgesFile = "C:/Users/ruesca/Dropbox/gradient/shedEdges_milwaukee.RData" # Milwaukee Data
source(paste(root, "/code/gradientFunctions.r", sep=""))
options(warn=2)
 
# args = commandArgs(trailingOnly=TRUE)
# 
# print(args)
# 
# edgeFile = args[1]
# nodeFile = args[2]
# shedFile = args[3]
# nodeRelFile = args[4]
# relFile = args[5]
# damFile = args[6]
# outCsv = args[7]
# outDbf = args[8]
# useSavedEdges = as.logical(args[9])
 
# edgeFile = paste(root, "/data/flowlines.dbf", sep="") # Milwaukee Data
edgeFile = paste(root, "/stateData/flowlines.dbf", sep="")
nodeFile = paste(root, "/stateData/nodes_rawElev.dbf", sep="")
shedFile = paste(root, "/stateData/sheds_rawElevMin.csv", sep="")
nodeRelFile = paste(root, "/relationshipFiles/noderelationships.csv", sep="")
relFile = paste(root, "/relationshipFiles/relationships.csv", sep="")
damFile = paste(root, "/stateData/dams.csv", sep="")
# outCsv = paste(root, "/temp/gradient_milwaukee.csv", sep="") # Milwaukee Data
# outDbf = paste(root, "/temp/gradient_milwaukee.dbf", sep="") # Milwaukee Data
outCsv = paste(root, "/gradient.csv", sep="")
outDbf = paste(root, "/gradient.dbf", sep="")
useSavedEdges = TRUE

# edgeCols = c(7,8,10,11) # Milwaukee Data
edgeCols = c(6,7,9,10)

if (!useSavedEdges) {
    shedEdges = formatShedEdges(edgeFile
                                , edgeCols
                                , nodeFile
                                , shedFile
                                , nodeRelFile
                                , relFile
                                , damFile
                                , outShedEdgesFile)
} else {
    load(outShedEdgesFile)
}

# Iterate over outlets
outlets = which(is.na(shedEdges$TOTRACEID))

print(" ")
print("Executing gradient-correction algorithm...")
print(" ")
for (outlet in outlets) {
    print(paste("Fixing upstream from REACHID", shedEdges$REACHID[outlet]))
    outletComplete=FALSE
    flagHeadwater=FALSE
    flagConfluence=TRUE
    row = outlet
    singleFeature = length(which(shedEdges$TOTRACEID == shedEdges$TRACEID[row])) == 0
    if (singleFeature) {
        print(paste(row, shedEdges[row,"seedtype"]))
        negGrad = shedEdges$maxElevFix1[row] < shedEdges$minElevFix1[row]
        if (negGrad) {
            shedEdges[row, c("minElevFix2","maxElevFix2")] = shedEdges[row, "minElevFix2"]
        } else {
            shedEdges[row, c("minElevFix2","maxElevFix2")] = shedEdges[row, c("minElevFix1","maxElevFix1")]
        }
        shedEdges$elevCheck[row] = 1
        next
    }
    shedEdges$minElevFix2[outlet] = shedEdges$minElevFix1[outlet]
    while (!outletComplete) {
        if (flagHeadwater) {
            flagHeadwater=FALSE
            flagConfluence=TRUE
            # Find the next highest drainage area where elevation has not yet been checked
            notFixed = shedEdges[is.na(shedEdges$elevCheck) & (shedEdges$outlet == outlet),]
            if (nrow(notFixed) == 0) {
                outletComplete = TRUE
                next
            }
            maxDrainArea = notFixed[which(notFixed$cellCount == max(notFixed$cellCount)),]
            # Select the TRACEID that flows into something other than itself (i.e., the downstream-most)
            maxDrainArea = maxDrainArea[!(maxDrainArea$TOTRACEID %in% maxDrainArea$TRACEID),]
            maxDrainArea = maxDrainArea[1,]
            row = which(shedEdges$TRACEID == maxDrainArea$TRACEID)
        }
        print(paste(row, shedEdges[row,"seedtype"]))
        if (is.na(shedEdges$minElevFix2[row])) {
            shedEdges$minElevFix2[row] = shedEdges$minElevFix1[row]
        }
        negGrad = shedEdges$maxElevFix1[row] < shedEdges$minElevFix2[row]
        if (!negGrad) {
            to = shedEdges$TOTRACEID[row]
            if (is.na(to)) {
                minElev = shedEdges[row, "minElevFix2"]
            } else {
                minElev = shedEdges[shedEdges$TRACEID==to,"maxElevFix2"]
            }
            is.stream = !(shedEdges[row, "seedtype"] %in% c("lake", "dangle"))
            if (is.stream) { 
                shedEdges$maxElevFix2[row] = shedEdges$maxElevFix1[row]
                if (is.na(shedEdges$minElevFix2[row])) {
                    shedEdges$minElevFix2[row] = minElev
                }
                shedEdges$elevCheck[row] = 1
                froms = shedEdges[which(shedEdges$TOTRACEID == shedEdges$TRACEID[row]),]
                # Anchor all minimum elevations of upstream connected nodes
                if (nrow(froms) > 0) {
                    shedEdges$minElevFix2[which(shedEdges$TOTRACEID == shedEdges$TRACEID[row])] = shedEdges$maxElevFix2[row]
                }
            } else {
                # Define existence of dam as a feature with a dam on it, or having a downstream feature with a dam on it
                if (is.na(to)) {
                    damBool = shedEdges$DAM[row]
                } else {
                    damBool = shedEdges$DAM[row] | shedEdges$DAM[shedEdges$TRACEID == to]
                }                
                # If a reservoir exists, we must use the raw elevation
                # Otherwise, minumum elevations will propogate due to 
                # lake gradients being coerced to zero.
                if (damBool) {
                    # First, check if raw gradient is negative
                    rawNegBool = shedEdges$minElevFix1[row] > shedEdges$maxElevFix1[row]
                    if (!rawNegBool) {
                        minElev = shedEdges$minElevFix1[row]
                    }
                }
                lakeRows = shedEdges[shedEdges$REACHID == shedEdges$REACHID[row],]
                shedEdges[which(shedEdges$TRACEID %in% lakeRows$TRACEID),c("minElevFix2", "maxElevFix2")] = minElev
                shedEdges$elevCheck[shedEdges$REACHID == shedEdges$REACHID[row]] = 1
                froms = shedEdges[which(shedEdges$TOTRACEID %in% lakeRows$TRACEID),]
                froms = froms[which(froms$REACHID != shedEdges$REACHID[row]),]
                # Anchor all minimum elevations of upstream connected nodes
                if (nrow(froms) > 0) {
                    shedEdges[which(shedEdges$TRACEID %in% froms$TRACEID), "minElevFix2"] = minElev
                }
            }
            # Find the row of the next upstream mainstem reach           
            if (nrow(froms) == 0) { 
                flagHeadwater = TRUE
                print("headwater")
            } else {
                from = froms[which(froms$cellCount == max(froms$cellCount))[1],]
                row = which(shedEdges$TRACEID == from$TRACEID)
            }
        } else {
            # If a negative gradient occurs at a confluence, do not select the downstream
            # segment as the low anchor.
            if (is.na(shedEdges$TOTRACEID[row])) {
                lakeImmediatelyDownstream = FALSE
            } else {
                lakeImmediatelyDownstream = shedEdges$seedtype[shedEdges$TRACEID == shedEdges$TOTRACEID[row]] %in% c("lake", "dangle")
            }
            if (flagConfluence | lakeImmediatelyDownstream) {
                anchorLow = shedEdges[row,]
                if (is.na(anchorLow$minElevFix2)) {
                    anchorLow$minElevFix2 = anchorLow$minElevFix1
                    shedEdges$minElevFix2[shedEdges$TRACEID == anchorLow$TRACEID] = 
                        anchorLow$minElevFix1
                }
                flagConfluence = FALSE
            } else {
                anchorLow = shedEdges[shedEdges$TRACEID == shedEdges$TOTRACEID[row],]
            }
            # Find the next segment where the minimum elevation is greater than the low anchor point
            betweenData = findHighAnchor(shedEdges, anchorLow, row)
            anchorHigh = betweenData$anchorHigh; interpolationIds = betweenData$interpolationIds
            shedEdges = interpolateElevations(shedEdges, anchorLow, anchorHigh, interpolationIds)
            shedEdges$elevCheck[shedEdges$TRACEID %in% interpolationIds] = 1
            # If the high anchor was irreconcilable, flag it as a headwater so the algorithm skips it
            # in the next iteration.
            highAnchorCheck = shedEdges$elevCheck[which(shedEdges$TRACEID == anchorHigh$TRACEID)]
            highAnchorChecked = !is.na(highAnchorCheck)
            if (highAnchorChecked) {
                flagHeadwater = TRUE
                print("headwater")
            } else {
                row = which(shedEdges$TRACEID == anchorHigh$TRACEID)
            }
        }
    }
}

recovery=shedEdges
shedEdges[c("minElevFix1", "maxElevFix1", "minElevFix2", "maxElevFix2")] =
              shedEdges[c("minElevFix1", "maxElevFix1", "minElevFix2", "maxElevFix2")] / 1000
gradient = ((shedEdges$maxElevFix2 - shedEdges$minElevFix2) / shedEdges$length) * 100
shedEdges$gradient = gradient

options(scipen=500)
shedEdges = shedEdges[order(shedEdges$gradient, decreasing=TRUE),]
shedEdges = shedEdges[c("REACHID", "minElevFix1", "maxElevFix1"
                        , "minElevFix2", "maxElevFix2","gradient")]
shedEdges = aggregate(shedEdges[,2:6], list(shedEdges[,1]), mean)
names(shedEdges) = c("REACHID","minElevRaw","maxElevRaw","minElevFix","maxElevFix","gradient")
write.dbf(shedEdges, file=outDbf)
write.csv(shedEdges, file=outCsv, row.names=F, na="")
