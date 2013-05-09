source("G:/gradient/code/gradientFunctions.r")

up = 207658
down = 245393

new = elevProfile(up, down, shedEdges, "minElevFix2", "maxElevFix2")
old = elevProfile(up, down, shedEdges, "minElevFix1", "maxElevFix1")


plot(minelv ~ l
     , ylim = c(min(c(old$minelv,new$minelv),na.rm=T)
                , max(c(old$minelv,new$minelv),na.rm=T))
     , data = old
     , xlab = "distance downstream (m)"
     , ylab = "minimum segment elevation (m)"
     , type = "l"
     , col = "blue"
     , lwd = 2)
lines(minelv ~ l, data = new, col = "red", lwd = 2)
