findLakeTribs = function (shedEdges, froms) {
    upstreamLakeRows = froms
    lakeIDs = froms$REACHID
    end = F
    while (!end) {
        froms = shedEdges[shedEdges$TOTRACEID %in% froms$TRACEID,]
        froms = froms[(froms$REACHID %in% lakeIDs) & !(froms$TRACEID %in% upstreamLakeRows$TRACEID),]
        if (nrow(froms) > 0 ) {
            upstreamLakeRows = rbind(upstreamLakeRows, froms)
        } else {
            end = T
        }
    }
    tribs = shedEdges[shedEdges$TOTRACEID %in% upstreamLakeRows$TRACEID
        & !(shedEdges$REACHID %in% lakeIDs),] # All tribs of lake
}



findHighAnchor = function (shedEdges, anchorLow, row) {
    confluence = shedEdges$TRACEID[row] == anchorLow$TRACEID
    if (!confluence) {
        to = shedEdges[row,]
        upstreamids = c(anchorLow$TRACEID, to$TRACEID)
    } else {
        to = anchorLow
        upstreamids = anchorLow$TRACEID
    }
    foundHighAnchor = FALSE
    while (!foundHighAnchor) {
        froms = shedEdges[which(shedEdges$TOTRACEID == to$TRACEID),]
        if (nrow(froms) == 0) {
            from = shedEdges[shedEdges$TRACEID == upstreamids[length(upstreamids)],]
            foundHighAnchor = TRUE
        } else {
            # Check if froms are segments in a lake
            # If yes, find the largest contributor to all TRACEIDs comprising lake
            # by selecting tribs and then re-selecting downstream
            if (all(froms$seedtype %in% c("lake", "dangle"))) {
                tribs = findLakeTribs(shedEdges,froms) # All tribs of lake
                if (nrow(tribs) == 0) { # Headwater lake
                    from = froms
                    foundHighAnchor = TRUE
                } else {
                    from = tribs[which(tribs$cellCount == max(tribs$cellCount))[1],]
                    # Trace back down to low anchor
                    upstreamids = from$TRACEID
                    down = from
                    while (!(anchorLow$TRACEID %in% upstreamids)) {
                        down = shedEdges[shedEdges$TRACEID == down$TOTRACEID,]
                        upstreamids = c(upstreamids, down$TRACEID)
                    }
                    # reverse upstreamids excluding the largest contributor (index 1)
                    upstreamids = upstreamids[length(upstreamids):2]
                    foundHighAnchor = from$minElevFix1 > anchorLow$minElevFix2
                }
            } else {
                from = froms[which(froms$cellCount == max(froms$cellCount))[1],]
                foundHighAnchor = from$minElevFix1 > anchorLow$minElevFix2
            }
        }
        if (!foundHighAnchor) {
            upstreamids = c(upstreamids, from$TRACEID)
        }
        to=from
    }
    return(list(anchorHigh=from, interpolationIds=upstreamids))
}

interpolateElevations = function (shedEdges, anchorLow, anchorHigh, interpolationIds) {
    flagIrreconcilable = FALSE
    # Because gradients of lakes are coerced to zero, we only calculate slope
    # over the length of streams in interpolation set
    lStreams = sum(shedEdges[shedEdges$TRACEID %in% interpolationIds
                             & shedEdges$seedtype != "lake", "length"])
    if (anchorHigh$minElevFix1 > anchorLow$minElevFix2) {
        slope = (anchorHigh$minElevFix1 - anchorLow$minElevFix2) / lStreams
    } else if (anchorHigh$maxElevFix1 > anchorLow$minElevFix2) {
        slope = (anchorHigh$maxElevFix1 - anchorLow$minElevFix2) / lStreams
    } else {
        flagIrreconcilable = TRUE
    }
    for (interpolationId in interpolationIds) {
        if (flagIrreconcilable) {
            shedEdges$maxElevFix2[which(shedEdges$TRACEID==interpolationId)] = NA
            shedEdges$minElevFix2[which(shedEdges$TOTRACEID==interpolationId)] = NA
        } else {
            is.stream = !(shedEdges[shedEdges$TRACEID==interpolationId, "seedtype"] %in% c("lake", "dangle"))
            if (is.stream) {
                l = shedEdges[shedEdges$TRACEID==interpolationId,"length"]
                to = shedEdges$TOTRACEID[shedEdges$TRACEID == interpolationId]
                if (is.na(to)) {
                    minElev = shedEdges[row, "minElevFix2"]
                } else {
                    minElev = shedEdges[shedEdges$TRACEID==to,"maxElevFix2"]
                }
                maxElev = minElev + (slope * l)
                # Anchor max elevation
                shedEdges$maxElevFix2[which(shedEdges$TRACEID==interpolationId)] = maxElev
                shedEdges$minElevFix2[which(shedEdges$TRACEID==interpolationId)] = minElev
                # Anchor min elevation of upstream tribs
                shedEdges$minElevFix2[which(shedEdges$TOTRACEID==interpolationId)] = maxElev
            } else {
                reachid = shedEdges$REACHID[shedEdges$TRACEID==interpolationId]
                lakeRows = shedEdges[shedEdges$REACHID == reachid,]
                # Anchor min elevations of all lake segments
                shedEdges[which(shedEdges$REACHID == reachid), c("minElevFix2", "maxElevFix2")] =
                    shedEdges$minElevFix2[which(shedEdges$TRACEID==interpolationId)]
                shedEdges$elevCheck[which(shedEdges$REACHID == reachid)] = 1
                # Anchor min elevation of upstream tribs
                shedEdges$minElevFix2[shedEdges$TOTRACEID %in% lakeRows$TRACEID] =
                    shedEdges$minElevFix2[which(shedEdges$TRACEID==interpolationId)]
            }
        }
    }
    if (flagIrreconcilable) {
        shedEdges[which(shedEdges$TRACEID==anchorHigh$TRACEID), c("minElevFix2", "maxElevFix2")] = NA
    }
    return(shedEdges)
}

elevProfile=function(up, down, shedEdges, minCol, maxCol){
    traceids=array(up)
    rownum=which(shedEdges$TRACEID==up)
    minelv=array(shedEdges[rownum, minCol])
    maxelv=array(shedEdges[rownum, maxCol])
    lengths=array(shedEdges$length[rownum])
    nextrow=which(shedEdges$TRACEID==shedEdges$TOTRACEID[rownum])
    x=1
    end = FALSE
    while (!end) {
        x=x+1
        traceids[x]=shedEdges$TRACEID[nextrow]
        maxelv[x]=shedEdges[nextrow, maxCol]
        minelv[x]=shedEdges[nextrow, minCol]
        lengths[x]=shedEdges$length[nextrow]
        nextrow=which(shedEdges$TRACEID==shedEdges$TOTRACEID[nextrow])
        atDown = shedEdges$TRACEID[nextrow] == down
        atOutlet = length(nextrow) == 0
        if (atDown | atOutlet) {end = TRUE}
        print(shedEdges$TRACEID[nextrow])
    }
    out=data.frame(traceids, maxelv, minelv, lengths)
    out$l=NA
    for (l in 1:nrow(out)){
        out$l[l]=sum(out$lengths[1:l])
    }
    return(out)
}

formatShedEdges = function() {
    #Steps in getting data ready for gradient analysis
    #In GIS
    #1) Extract Values to Points with node points and raw elevation
    #2) Zonal Statistics As Table for minimum raw elevation value for watersheds
    #3) Export table from #2 into a csv file so it can be easily imported into R
    
    print("Reading input data...")
    #nodes with elevations
    nodes=read.dbf(nodeFile)
    #sheds with elevations
    sheds=read.csv(shedFile)
    edges=read.dbf(edgeFile)
    #read in node relationships and edge relationships (from main database)
    noderels=read.csv(nodeRelFile)
    rels=read.csv(relFile)
    print("Combining input data...")
    #combine edges with sheds to get minimum elevations for each REACHi
    edgeMinElev=merge(edges[,edgeCols], sheds, by.x="REACHID", by.y="key", all.x=TRUE)
    colnames(edgeMinElev)[6]="minElevRaw"
    #get maxElev for features that are sources
    sources=nodes[which(nodes$node_cat=="Source"),c(4:6)]
    sources2=merge(sources[,c(2:3)], noderels[,c(2:3)], by.x="pointid", by.y="fromnode")
    edgeElev=merge(edgeMinElev, sources2, by.x="TRACEID", by.y="TRACEID", all.x=TRUE)
    colnames(edgeElev)[8]="maxElevRaw"
    #get to/from info for each segment
    finalEdges=merge(edgeElev, rels[, 2:3], by.x="TRACEID", by.y="FROM_TRACEID", all.x=TRUE)
    finalEdges=finalEdges[,c(1:6, 8:9)]
    colnames(finalEdges)[4]="length"
    
    #table needs these columns: TRACEID, REACHID, minElevRaw, maxElevRaw, TO_TRACEID, minElevC1
    
    #change up topology to get from-to relationships that exclude those segments without sheds
    #NA in new_TO for features that do not go anywhere (either becuase isolated or because go to features without sheds (ie across Lake Michigan/Superior HUC
    print("Excluding features without an associated watershed...")
    new_TO_TRACEID=array(NA, nrow(finalEdges))
    for (f in 1:nrow(finalEdges)){
        flag="up"
        rownum=f
        while(flag=="up"){
            rownum=which(finalEdges$TRACEID==finalEdges[["TO_TRACEID"]][rownum])
            if (length(rownum)==0){ #basically, if at end of trace, no "to" feature
                flag="down"
            } else if (is.na(finalEdges$minElevRaw[rownum])==FALSE){ #if the min elevation is not NA, i.e. we got to a segment with a shed
                flag="down"
            }
        }
        new_TO_TRACEID[f]=ifelse(length(rownum)==0, NA, finalEdges$TRACEID[rownum])
    }
    finalEdges$new_TO=new_TO_TRACEID
    
    # any headwater feature too small to have its own shed needs (thus has NA for minimum shed elevation)to have its 
    # source elevation assigned to first feature downstream with shed
    # if this occurs at a confluence, we might be in trouble with script so far...
    print("Assigning downstream ID to features lacking a watershed...")
    finalEdges$maxElevFix1=finalEdges$maxElevRaw
    fixes=finalEdges[which(is.na(finalEdges$minElevRaw)==TRUE & is.na(finalEdges$maxElevRaw)==FALSE),]
    for (f in 1:nrow(fixes)){
        rownum=which(finalEdges$TRACEID==fixes$new_TO[f])
        finalEdges$maxElevFix1[rownum]=fixes$maxElevFix1[f]
    }
    #subset out features we are interested in, namely those with watersheds (which means they should have minimum shed elevation data
    shedEdges=finalEdges[which(is.na(finalEdges$minElevRaw)==FALSE),]
    
    
    #rectify values at confluences
    print("Rectifying elevations at confluences...")
    shedEdges$minElevFix1=shedEdges$minElevRaw #create column to store new elevations
    to_freq=as.data.frame(table(shedEdges$new_TO)) 
    confluence=subset(to_freq, to_freq$Freq>1) #select out features that more than one feature goes to
    for (v in 1:nrow(confluence)){
        rows=which(shedEdges$new_TO==confluence$Var1[v])
        shedEdges$minElevFix1[rows]=min(shedEdges$minElevFix1[rows])
    }
    
    #assign maximum elevations to features based on minimum elevations of upstream features
    #I take the minimum but values should be the same because of confluence rectifying above
    print("Assigning maximum elevations based on minimum elevations of upstream features...")
    for (f in 1:nrow(shedEdges)){
        if (is.na(shedEdges$maxElevFix1)[f]==TRUE){
            upsegs=which(shedEdges$new_TO==shedEdges$TRACEID[f])
            shedEdges$maxElevFix1[f]=min(shedEdges$minElevFix1[upsegs])
        }
    }
    
    print("Formatting main table...")
    shedEdges = shedEdges[,c("TRACEID", "new_TO", "cellCount", "REACHID", "seedtype", "length", "minElevFix1", "maxElevFix1")]
    shedEdges[c("minElevFix2", "maxElevFix2", "elevCheck")] = NA
    names(shedEdges)[2] = "TOTRACEID"
    shedEdges$minElevFix1 = as.integer(round(shedEdges$minElevFix1 * 1000))
    shedEdges$maxElevFix1 = as.integer(round(shedEdges$maxElevFix1 * 1000))
    save(shedEdges,file=outShedEdgesFile)
    return(shedEdges)
}
