setwd("C:/Users/SIMILIER/Desktop/All document/SIMILIEN/PHD FILES/PHD_FILE/Data/DHS 2020")
library(rgdal)
library(geoR)
library(haven)
stuntingdata <- read_dta("childrendata2020.dta")
Rwanda <- readOGR(dsn = "RWGE81FL", layer = "RWGE81FL")
Rwanda
head(Rwanda@data, n = 5)
library(ggmap)
library(tmap)
library(rgdal)
library(dplyr)
library(lme4)
stuntingdataset<-merge(stuntingdata,Rwanda,by="DHSCLUST", all.x=TRUE,all.y=TRUE )
save(stuntingdataset, file = "stuntingdataset.RData")
stuntingdataset$stunting
dim(unique(stuntingdataset[, c("LONGNUM","LATNUM")]))
library(dplyr) 
d <- group_by(stuntingdataset, LONGNUM,LATNUM,DHSREGNA) %>% 
  summarize( total = n(), Stunted = sum(stunting), 
             prev = Stunted / total ) 
library(sp) 
library(rgdal) 
sps <- SpatialPoints(d[, c("LONGNUM","LATNUM")],
                     proj4string = CRS("+proj=utm +zone=36") )
spst <- spTransform(sps, CRS("+proj=longlat +datum=WGS84"))
d[, c("LONG","LAT")] <- coordinates(spst) 
d
library(leaflet) 
library(viridis)
library(mapview)
library(maps)
pal <- colorBin(c("green","#0018F9","#FFD300","red"), bins = c(0, 0.25, 0.5, 0.75, 1)) 
leaflet(d) %>% 
  addProviderTiles(providers$OpenStreetMap.HOT) %>% 
  addCircles(lng = ~LONGNUM, lat = ~LATNUM, color = ~ pal(prev)) %>% 
  addLegend("bottomright", pal = pal, values = ~prev, title = "Stunting prevalence" ) %>% 
  addScaleBar(position = c("bottomleft")) 
library(raster) 
r<- getData('alt', country = 'RWA', mask = TRUE)
head(r)
pal <- colorNumeric(c("green", "#EFFD5F","yellow","red"), values(r), na.color = "transparent" )
leaflet(d) %>% addProviderTiles(providers$OpenStreetMap.HOT) %>% 
  addRasterImage(r, colors = pal, opacity = 0.5) %>% 
  addLegend("bottomright", pal = pal, values = values(r),
            title = "Altitude" ) %>% addScaleBar(position = c("bottomleft"))
d$alt <- raster::extract(r, d[, c("LONGNUM","LATNUM")])
head(d)
library(INLA) 
coo <- cbind(d$LONGNUM, d$LATNUM) 
mesh <- inla.mesh.2d( 
  loc = coo, max.edge = c(0.1, 5), 
  cutoff = 0.1 )
mesh$n
plot(mesh)
points(coo, col = "red")
spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)
indexs <- inla.spde.make.index("s", spde$n.spde)
lengths(indexs)
dp <- rasterToPoints(r)
dim(dp)
ra <- aggregate(r, fact = 5, fun = mean)
A <- inla.spde.make.A(mesh = mesh, loc = coo)
dp <- rasterToPoints(ra)
dim(dp)
coop <- dp[, c("x", "y")]
Ap <- inla.spde.make.A(mesh = mesh, loc = coop)
# stack for estimation stk.e
stk.e <- inla.stack(
  tag = "est",
  data = list(y = d$Stunted, numtrials = d$total),
  A = list(1, A),
  effects = list(data.frame(b0 = 1, altitude = d$alt), s = indexs)
)
# stack for prediction stk.p
stk.p <- inla.stack(
  tag = "pred",
  data = list(y = NA, numtrials = NA),
  A = list(1, Ap),
  effects = list(data.frame(b0 = 1, altitude = dp[, 3]),
                 s = indexs
  )
)

# stk.full has stk.e and stk.p
stk.full <- inla.stack(stk.e, stk.p)
formula <-y ~ 0 + b0 + altitude+ f(s, model = spde)
res <- inla(formula,
            family = "binomial", Ntrials = numtrials,
            control.family = list(link = "logit"),
            data = inla.stack.data(stk.full),
            control.predictor = list(
              compute = TRUE, link = 1,
              A = inla.stack.A(stk.full)
            ))
res$summary.fitted.values
res$marginals.fitted.values
summary(res)
res$summary.hyperpar
#stunting prevalence predictions
index <- inla.stack.index(stack = stk.full, tag = "pred")$data
prev_mean <- res$summary.fitted.values[index, "mean"]
prev_ll <- res$summary.fitted.values[index, "0.025quant"]
prev_ul <- res$summary.fitted.values[index, "0.975quant"]
pal <- colorNumeric(c("green", "#EFFD5F","yellow","red"), c(0, 1), na.color = "transparent")
leaflet() %>%
  addProviderTiles(providers$OpenStreetMap.HOT) %>%
  addCircles(
    lng = coop[, 1], lat = coop[, 2],
    color = pal(prev_mean)
  ) %>%
  addLegend("bottomright",
            pal = pal, values = prev_mean,
            title = "Prev."
  ) %>%
  addScaleBar(position = c("bottomleft"))
r_prev_mean <- rasterize(
  x = coop, y = ra, field = prev_mean,
  fun = mean
)
#Map of stunting prevalence predictions in The Rwanda created with leaflet() and addRasterImage().
pal <- colorNumeric( c("green","yellow","#BF0A30"), c(0,1), na.color = "transparent")
leaflet() %>%
  addProviderTiles(providers$OpenStreetMap.HOT) %>%
  addRasterImage(r_prev_mean, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright",
            pal = pal,
            values = values(r_prev_mean), title = "Prediction prevalence"
  ) %>%
  addScaleBar(position = c("bottomleft"))
#Approach to create maps with the lower limits of the prevalence predictions.
r_prev_ll <- rasterize(
  x = coop, y = ra, field = prev_ll,
  fun = mean
)
leaflet(d) %>%
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addRasterImage(r_prev_ll, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright",
            pal = pal,
            values = values(r_prev_ll), title = "Lower Limits"
  ) %>%
  addScaleBar(position = c("bottomleft"))
#Approach to create maps with the upper limits of the prevalence predictions.
r_prev_ul <- rasterize(
  x = coop, y = ra, field = prev_ul,
  fun = mean
)
leaflet() %>%
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addRasterImage(r_prev_ul, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright",
            pal = pal,
            values = values(r_prev_ul), title = "Upper Limits"
  ) %>%
  addScaleBar(position = c("bottomleft"))
index <- inla.stack.index(stack = stk.full, tag = "pred")$data
marg <- res$marginals.fitted.values[index][[1]]
#k<-inla.pmarginal(q =0.3, marginal = marg)
#k
1 - inla.pmarginal(q =0.3, marginal = marg)
excprob <- sapply(res$marginals.fitted.values[index],
                  FUN = function(marg){1-inla.pmarginal(q = 0.3, marginal = marg)})
head(excprob)
#PLOT EXCEDENCE PROBABILITY 
r_excprob <- rasterize(
  x = coop, y = ra, field = excprob,
  fun = mean
)
pal <- colorNumeric(c("green","yellow","red"), c(0, 1), na.color = "transparent")
leaflet() %>%
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addRasterImage(r_excprob, colors = pal, opacity = 0.5) %>%
  addLegend("bottomright",
            pal = pal,
            values = values(r_excprob), title = "P(p>0.3)"
  ) %>%
  addScaleBar(position = c("bottomleft"))
save(d, file = "d.RData")
