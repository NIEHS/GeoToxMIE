######################################################
# NATA and Tox21
# By: Kristin Eccles
# Date: May7th, 2021
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

# Import data
hill2.fit <- get(load("/Volumes/SHAG/GeoTox/data/Hill_2param_model_fit.RData"))

#NATA
nata_df<- read.csv("/Volumes/SHAG/GeoTox/data/2014_NATA_CONCS.csv")
nata_chemicals <- read.csv("/Volumes/SHAG/GeoTox/data/NATA_pollutant_names_casrn.csv")
county_2014 <-st_read("/Volumes/SHAG/GeoTox/data/cb_2014_us_county_5m/cb_2014_us_county_5m.shp")

# There were some repeats of reading in data?!
#simplify names
nata_chemicals$web_name <- str_to_title(nata_chemicals$web_name, locale = "en")
nata_chemicals$web_name <- str_replace(nata_chemicals$web_name , " \\s*\\([^\\)]+\\)", "")

# aggregate to county
nata_county_mean <- aggregate(nata_df[,7:ncol(nata_df)] , by=list(nata_df$STCOFIPS) , FUN=mean) 
nata_county_sd <- aggregate(nata_df[,7:ncol(nata_df)] , by=list(nata_df$STCOFIPS) , FUN=sd) 

#melt
county_stack_mean <- melt(nata_county_mean, id=1)
colnames(county_stack_mean) <- cbind("FIPS", "chemical", "concentration_mean")

county_stack_sd <- melt(nata_county_sd, id=1)
colnames(county_stack_sd) <- cbind("FIPS", "chemical", "concentration_sd")

county_stack <- left_join(county_stack_mean, county_stack_sd, by=c("FIPS","chemical"), keep= FALSE)

# add casrn
county_stack <- left_join(county_stack, nata_chemicals, by=c("chemical" = "smoke_name"), keep= FALSE)

# Create a STATE ID 
county_stack$STATE <- county_stack$FIPS %>% as.character()
county_stack$STATE[nchar(county_stack$STATE)==4] <- str_pad(substr(county_stack$STATE[nchar(county_stack$STATE)==4],1,1),
                                                            width = 2,side = "left",pad = "0")
county_stack$STATE[nchar(county_stack$STATE)==5] <- substr(county_stack$STATE[nchar(county_stack$STATE)==5],1,2)


#limit to continental USA
county_stack<- subset(county_stack,  STATE != "02" )
county_stack<- subset(county_stack,  STATE != "15" )
county_stack<- subset(county_stack,  STATE != "60" )
county_stack<- subset(county_stack,  STATE != "66" )
county_stack<- subset(county_stack,  STATE != "69" )
county_stack<- subset(county_stack,  STATE != "72" )
county_stack<- subset(county_stack,  STATE != "78" )

### filter county_stack by the 41 chemicals in the analysis 

county_stack_filter <- filter(county_stack,casrn %in% hill2.fit$casn)

################################################################################
#### Join with the CYP1A1 2-parameter hill model fit data ####
county_cyp1a1_up <- left_join(county_stack_filter, hill2.fit, by=c("casrn" = "casn"), keep= FALSE)

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

################################################################################
#### CYP1A1 ####
hill2.fit$chemical_name <- rownames(hill2.fit)
hill2.fit$aenm <- "LTEA_HepaRG_CYP1A1_up"

#### Heatmap ####
heatmap_df<- as.data.frame(hill2.fit) %>%
  group_by(aenm, chemical_name)%>%
  summarise(mean=mean(10^logAC50))%>%
  as.data.frame()
heatmap_df
summary(heatmap_df$mean)

heatmap_cyp1a1_up <- ggplot(heatmap_df, aes(x = fct_reorder(chemical_name, mean),y = aenm,  fill= mean)) + 
  geom_tile()+
  theme_bw()+
  scale_fill_viridis_c(direction = 1,option = "A") +
  labs(y="Assay", x="Chemical", fill = "AC50 (μM)")+
  theme(axis.text.x = element_text(angle = 90, hjust=1), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.ticks.x = element_blank())+
  theme(text = element_text(size = 8), aspect.ratio=3/10)
heatmap_cyp1a1_up
#save_plot("heatmap_cyp1a1_up.tif", heatmap_cyp1a1_up, width = 40, height = 30, dpi = 200)

##### Facet Wrap Map ####
cyp1a1_up_nata_county_sp <- left_join(county_2014, nata_county_stack, by=c("countyid" = "FIPS"), keep=FALSE)
cyp1a1_up_nata_county_sf <-st_as_sf(cyp1a1_up_nata_county_sp)

cyp1a1_up_all_map <- ggplot(data = cyp1a1_up_nata_county_sf, aes(fill=concentration_norm)) +
  geom_sf(lwd = 0)+
  facet_wrap(vars(web_name), ncol=6, nrow=7)+
  theme_bw()+
  geom_sf(data = subset(cyp1a1_up_nata_county_sf, concentration == 0), aes(fill=concentration), fill = "light grey", size=0)+
  geom_sf(data = states, fill = NA, size=0.15)+
  scale_fill_viridis_c(name = "Normalized Concentration", direction = -1,option = "A") +
  theme(axis.ticks.x=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.y=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(),
        legend.text=element_text(angle = 45), legend.title=element_text(size=12), 
        legend.position="bottom", 
        strip.text = element_text(size = 6))
cyp1a1_up_all_map
#cyp1a1_up_all_map
#save_plot("facet_cyp1a1_up_nata.tif", cyp1a1_up_all_map, width = 40, height = 30, dpi = 200)


#### Count Map ####
nata_county_stack <-  mutate(nata_county_stack, 
                                 count = ifelse(concentration >0, 1, 0))
# cyp1a1_up composite count
count_county_cyp1a1_up <- as.data.frame(nata_county_stack) %>%
  group_by(FIPS) %>%
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
  theme(text = element_text(size = 12), axis.ticks.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), 
        legend.position="bottom")
count_cyp1a1_up_map
#save_plot("count_cyp1a1_up_uM_sumACC.tif", count_cyp1a1_up_map, width = 40, height = 30, dpi = 200)

#### Combine Plots ####
cyp1a1_up1=ggarrange(count_cyp1a1_up_map, heatmap_cyp1a1_up, 
                     labels = c("B", "C"),
                     #vjust = -0.005,
                     #hjust = -0.5,
                     ncol = 1, nrow = 2,
                     #heights = c(1, 0.25), 
                     font.label = list(size = 16, color = "black", face = "bold"),
                     common.legend = FALSE)

cyp1a1_up=ggarrange(cyp1a1_up_all_map, cyp1a1_up1,
                    labels = c("A"),
                    vjust = 1,
                    #hjust = -0.5,
                    ncol = 2, nrow = 1,
                    widths = c(1, 0.5),
                    font.label = list(size = 16, color = "black", face = "bold"),
                    common.legend = FALSE)

#Plot figures with dpi=300
save_plot("/Volumes/SHAG/GeoTox/data/plots/cyp1a1_up_figure2.tif", cyp1a1_up, width = 32, height = 20, dpi = 300)
