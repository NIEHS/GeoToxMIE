######################################################
# Manuscript Figure 3 
# By: Kristin Eccles
# Date: March 7th, 2022
# Written in R Version 4.0.2
######################################################
# Libraries
library(ggplot2)
library(viridis)
library(dplyr)
library(tigris)
library(sf)
library(reshape2)
library(sjPlot)
library(stringr)
library(maps)
library(ggpubr)
library(forcats)

# Import data
hill2.fit <- get(load("/Volumes/SHAG/GeoTox/data/Hill_2param_model_fit.RData"))
hill2.fit$chnm <- rownames(hill2.fit)

hill2.fit$chnm<- str_replace_all(hill2.fit$chnm, "C.I. Disperse Black 6", "3,3'-Dimethoxybenzidine")
hill2.fit$chnm <- str_replace_all(hill2.fit$chnm, "C.I. Azoic Diazo Component", "Benzidine")

#NATA
nata_df<- read.csv("/Volumes/SHAG/GeoTox/data/2014_NATA_CONCS.csv")
nata_chemicals <- read.csv("/Volumes/SHAG/GeoTox/data/NATA_pollutant_names_casrn.csv")


ice_data <- get(load("/Volumes/SHAG/GeoTox/data//210105_ICE_cHTS_invitrodbv33.Rdata"))
ice_data <- subset(ice_data, new_hitc == 1)
states <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))

county_2014 <-st_read("/Volumes/SHAG/GeoTox/data/cb_2014_us_county_5m/cb_2014_us_county_5m.shp")
county_2014$GEOID <- as.numeric(county_2014$GEOID)

#limit to continental USA
county_2014<- subset(county_2014,  STATEFP != "02" )
county_2014<- subset(county_2014,  STATEFP != "15" )
county_2014<- subset(county_2014,  STATEFP != "60" )
county_2014<- subset(county_2014,  STATEFP != "66" )
county_2014<- subset(county_2014,  STATEFP != "69" )
county_2014<- subset(county_2014,  STATEFP != "72" )
county_2014<- subset(county_2014,  STATEFP != "78" )
county_2014$countyid <-as.numeric(paste0(county_2014$STATEFP, county_2014$COUNTYFP))

###############################################################################
# Determine what the overlap between NATA chemicals and tox21 chemicals
#simplify names
nata_chemicals$web_name <- str_to_title(nata_chemicals$web_name, locale = "en")
nata_chemicals$web_name <- str_replace(nata_chemicals$web_name , " \\s*\\([^\\)]+\\)", "")

nata_tox21 <- merge(ice_data, nata_chemicals,  by.x = "casn", by.y = "casrn") # n = 3533
nata_tox21$presence <- "1" 
length(unique(nata_tox21$aenm)) # 643 assays

nata_tox21_unique <- nata_tox21[!duplicated(nata_tox21$chnm), ]
length(nata_tox21_unique$chnm) # n= 83 unique chemicals

# summarized number of chemicals by assay
nata_tox21_count <- as.data.frame(nata_tox21) %>%
  group_by(nata_tox21$aenm) %>%
  summarise(count = n())
nata_tox21_count <- as.data.frame(nata_tox21_count[,1:2])
colnames(nata_tox21_count) <- cbind("assay name", "chemical_count")
nata_tox21_count
#write.csv(nata_tox21_count, "nata_tox21_count.csv")

# aggregate to county
nata_county <- aggregate(nata_df[,7:ncol(nata_df)] , by=list(nata_df$STCOFIPS) , FUN=mean) 

#melt
nata_county_stack <- melt(nata_county, id=1)
colnames(nata_county_stack) <- cbind("FIPS", "chemical", "concentration")

# add casrn
nata_county_stack <- left_join(nata_county_stack, nata_chemicals, by=c("chemical" = "smoke_name"), keep= FALSE)

#remove any missing data
nata_county_stack <- na.omit(nata_county_stack)

# range normalize air concentrations
fac = interaction(factor(nata_county_stack$chemical), drop = T)
fnorm.var = function(y){ return((y - min(y)) / (max(y) - min(y)))}
fnorm.tab = function(tab) {apply(tab, 2, fnorm.var)}
concentration_norm=unsplit(lapply(split(as.data.frame(nata_county_stack$concentration), fac), 
                                  FUN = function(x) as.data.frame(fnorm.tab(x))), fac)
colnames(concentration_norm) <- "concentration_norm"
nata_county_stack <- cbind(nata_county_stack, concentration_norm)


# limit to chemicals avaliable in httk

in.chems <-  c("98-86-2","92-87-5","92-52-4","117-81-7","133-06-2","532-27-4","133-90-4","57-74-9","510-15-6","94-75-7" ,
               "64-67-5","132-64-9","106-46-7","111-44-4","79-44-7","131-11-3","77-78-1","119-90-4","121-14-2","534-52-1", 
               "51-28-5","121-69-7","107-21-1","51-79-6","76-44-8","822-06-0","77-47-4","123-31-9","72-43-5" , 
               "101-77-9","56-38-2","82-68-8","87-86-5","1120-71-4", "114-26-1","91-22-5","96-09-3","95-80-7","584-84-9" ,
               "95-95-4","1582-09-8")

nata_county_stack <-nata_county_stack[(nata_county_stack$casrn %in% in.chems),]

# subset tox21 data
ice_cyp1a1_up <- subset(ice_data, aenm == "LTEA_HepaRG_CYP1A1_up" )

nata_county_cyp1a1_up <- left_join(nata_county_stack, ice_cyp1a1_up, by=c("casrn" = "casn"), keep= FALSE)
nata_county_cyp1a1_up<- nata_county_cyp1a1_up[!is.na(nata_county_cyp1a1_up$ACC), ]
nata_county_cyp1a1_up<- nata_county_cyp1a1_up[!is.na(nata_county_cyp1a1_up$concentration), ]
nata_county_cyp1a1_up <-  mutate(nata_county_cyp1a1_up, 
                                 count = ifelse(concentration >0, 1, 0))

nata_county_cyp1a1_up$chnm<- str_replace_all(nata_county_cyp1a1_up$chnm, "C.I. Disperse Black 6", "3,3'-Dimethoxybenzidine")
nata_county_cyp1a1_up$chnm <- str_replace_all(nata_county_cyp1a1_up$chnm, "C.I. Azoic Diazo Component", "Benzidine")

#### Heatmap ####
heatmap_cyp1a1_df<- left_join(nata_county_cyp1a1_up, hill2.fit, by= c("casrn" = "casn"), keep=FALSE )
heatmap_df<- as.data.frame(heatmap_cyp1a1_df) %>%
  group_by(aenm, chnm.y) %>%
  summarise(mean=mean(10^(logAC50)))%>%
  as.data.frame()
summary(heatmap_df$mean)

heatmap_cyp1a1_up <- ggplot(heatmap_df, aes(x = fct_reorder(chnm.y, mean),y = aenm,  fill= mean)) + 
  geom_tile()+
  theme_bw()+
  scale_fill_viridis_c(direction = 1,option = "A") +
  labs(y="Assay", x="Chemical", fill = "EC50 (μM)")+
  theme(axis.text.x = element_text(angle = 45, hjust=1), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        legend.title=element_text(size=14),
        text = element_text(size = 18), 
        aspect.ratio=1.5/10,
        plot.margin = margin(10, 10, 10, 100))
heatmap_cyp1a1_up
#save_plot("heatmap_cyp1a1_up.tif", heatmap_cyp1a1_up, width = 40, height = 30, dpi = 200)

##### Facet Wrap Map ####
nata_county_cyp1a1_up<- left_join(nata_county_cyp1a1_up, hill2.fit, by= c("casrn" = "casn"), keep=FALSE )
cyp1a1_up_nata_county_sp <- left_join(county_2014, nata_county_cyp1a1_up, by=c("countyid" = "FIPS"), keep=FALSE)
cyp1a1_up_nata_county_sf <-st_as_sf(cyp1a1_up_nata_county_sp)

cyp1a1_up_all_map <- ggplot(data = cyp1a1_up_nata_county_sf, aes(fill=concentration_norm)) +
  geom_sf(lwd = 0)+
  facet_wrap(vars(fct_reorder(chnm.x, logAC50)), ncol=7)+
  theme_bw()+
  geom_sf(data = subset(cyp1a1_up_nata_county_sf, concentration == 0), aes(fill=concentration), fill = "light grey", size=0)+
  geom_sf(data = states, fill = NA, size=0.15)+
  scale_fill_viridis_c(name = "Normalized Concentration", direction = -1,option = "A") +
  theme(axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(),
        legend.text=element_text(angle = 45),
        legend.position="bottom", 
        strip.text = element_text(size = 10),
        text = element_text(size = 18))
#save_plot("facet_cyp1a1_up_nata.tif", cyp1a1_up_all_map, width = 40, height = 30, dpi = 200)


#### Count Map ####
# cyp1a1_up composite count
count_county_cyp1a1_up <- as.data.frame(nata_county_cyp1a1_up) %>%
  group_by(nata_county_cyp1a1_up$FIPS) %>%
  summarise(sum = sum(count))
count_county_cyp1a1_up <- as.data.frame(count_county_cyp1a1_up)
colnames(count_county_cyp1a1_up) <- cbind("FIPS", "count")
summary(count_county_cyp1a1_up)

count_county_cyp1a1_up_sp<- left_join(county_2014, count_county_cyp1a1_up, by=c("countyid" = "FIPS"), keep=FALSE)
count_county_cyp1a1_up_sf <-st_as_sf(count_county_cyp1a1_up_sp)

count_cyp1a1_up_map <- ggplot(data = count_county_cyp1a1_up_sf, aes(fill=count)) +
  geom_sf(lwd = 0)+
  theme_bw()+
  geom_sf(data = states, fill = NA, size=0.15)+
  scale_fill_viridis_c(name="# Chemicals", direction = -1,option = "A") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), 
        legend.position="right", 
        text = element_text(size = 18), 
        aspect.ratio=5/10)
count_cyp1a1_up_map
#save_plot("count_cyp1a1_up_uM_sumAC50.tif", count_cyp1a1_up_map, width = 40, height = 30, dpi = 200)

#### Combine Plots V ####
cyp1a1_v=ggarrange(cyp1a1_up_all_map, count_cyp1a1_up_map, heatmap_cyp1a1_up, 
                     labels = c("A", "B", "C"),
                     #vjust = -0.005,
                     #hjust = -0.5,
                     widths = c(1, 0.65, 0.5),
                     ncol = 1, nrow = 3,
                     heights = c(1, 0.65, 0.5), 
                     align = "v",
                     font.label = list(size = 20, color = "black", face = "bold"),
                     common.legend = FALSE)



#Plot figures with dpi=300
save_plot("/Volumes/SHAG/GeoTox/data/plots/Figure3_v20220420.tif", cyp1a1_v, width = 45, height = 55, dpi = 300)





