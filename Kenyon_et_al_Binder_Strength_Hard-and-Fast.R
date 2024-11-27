# Proc B Submission Code

# Submission 28 Nov 2024

# # Libraries ----

library(tidyverse)
library(glmmTMB)
library(emmeans)
library(DHARMa)
library(MuMIn)
library(lme4)
library(effects)
library(gridExtra)
library(lubridate)
library(car)
library(flextable) # for saving tables
library(officer) # for saving tables
library(broom.mixed)
library(vegan)

# Writing a function for saving the tables:

save_tables_to_docx <- function(emmeans_obj, 
                                emmeans_filename, 
                                contrasts_filename) {
  # Extract emmeans and contrasts
  emmeans_df <- emmeans_obj$emmeans %>% as.data.frame()
  contrasts_df <- emmeans_obj$contrasts %>% as.data.frame()
  
  # Set flextable defaults
  set_flextable_defaults(font.family = "Times New Roman", font.size = 9,
                         padding.bottom = 0, padding.top = 0)
  
  # Create flextables
  emmeans_table <- flextable(emmeans_df) %>% 
    colformat_double(big.mark = ",", digits = 3, na_str = "N/A")
  
  contrasts_table <- flextable(contrasts_df) %>% 
    colformat_double(big.mark = ",", digits = 3, na_str = "N/A")
  
  # Save flextables to Word documents
  save_as_docx(emmeans_table, path = emmeans_filename)
  save_as_docx(contrasts_table, path = contrasts_filename)
}

#save_tables_to_docx(glmmTemp_2comps, "glmmTemp_2compsT.docx", "glmmTemp_2compsT2.docx")


# Cleaning & formatting -----

gridid <- read.csv("gridid3.csv", header =TRUE)
pairid <- read.csv("pairid2.csv", header =TRUE)
binderid <- read.csv("binderid3.csv", header = TRUE)
pieceid <- read.csv("pieceid2.csv", header =TRUE)
gapid <- read.csv("gapid2.csv", header =TRUE) 
contactid <- read.csv("contactid2.csv", header = TRUE) 

cols1 <- c("grid_id", "org_grid_id", "pair_arrangement", "experiment", "location", "site", "depth", "exposure", 
           "layer")
gridid <- gridid %>% mutate_at(cols1, factor)
sapply(gridid, class)

gridid$site_unique <- paste(gridid$site, gridid$depth, sep = "_")
gridid$site_unique <- as.factor(gridid$site_unique)
levels(gridid$site_unique) # 16 levels of site_unique (16 sites according to exposure and depth)

gridid$site_unique <- factor(gridid$site_unique, levels = c("Humpy Sheltered_Shallow slope", 
                                                            "Halfway Sheltered_Shallow slope", 
                                                            "Clam Bay_Shallow slope", 
                                                            "Humpy Exposed_Shallow slope", 
                                                            "Halfway Exposed_Shallow slope", 
                                                            "Red Beach_Shallow slope", 
                                                            "Coral Gardens_Reef flat",
                                                            "Halfway_Reef flat", 
                                                            "Coral Gardens_Shallow slope", 
                                                            "Coral Gardens_Deep slope", 
                                                            "Halfway_Shallow slope",
                                                            "Halfway_Deep slope", 
                                                            "Eco1_Shallow slope" ,"Eco1_Deep slope",
                                                            "Eco2_Shallow slope" , "Eco2_Deep slope"))


cols2 <- c("grid_id", "bound", "unique_pair_id", "force_type", "experiment")
cols3 <- c("break_force", "turf_height_cm")

pairid <- pairid %>%  mutate_at(cols2, factor) %>% mutate_at(cols3, as.numeric)
sapply(pairid, class)

str(gridid)

cols4 <- c("grid_id", "unique_pair_id", "experiment")

cols6 <- c("bind_length_x_thickness","bind_length", "bind_thickness")

binderid <- binderid %>% mutate_at(cols4, factor)  %>%
  mutate_at(cols6, as.numeric)

cols7 <- c("grid_id", "unique_pair_id", "experiment", "side2")

gapid <- gapid %>% mutate_at(cols7, factor)
contactid <- contactid %>% mutate_at(cols7, factor)

cols5 <- c("grid_id", "piece_num", "unique_pair_id", "experiment")
cols8 <- c("axial_length", "width_1", "width_2", "width_3")
pieceid <- pieceid %>% mutate_at(cols5, factor) %>%
  mutate_at(cols8, as.numeric)

# Add a days_deployed column to gridid

gridid <- gridid %>%  
  mutate(date_deployed = dmy(date_deployed),
         date_collected = dmy(date_collected),
         days_deployed = as.numeric(difftime(date_collected, date_deployed, units = "days"))) 

gridid <- gridid %>% 
  mutate(date_deployed30 = months_deployed*30)

# MAIN EXPERIMENT DATAFRAME -----

levels(gridid$experiment)

# gridid & pairid

gridmain <- gridid %>% filter(experiment == "main") %>% droplevels()
pairmain <- pairid %>% filter(experiment == "main") %>% droplevels()

names(gridmain)
names(pairmain)

gridpair <- inner_join(gridmain, pairmain, by="grid_id")

write.csv(gridpair, file = "gridpair3.csv") # checked this

# gridid & binderid

bindermain <- binderid %>% filter(experiment == "main") %>% droplevels()
gridbinder <- inner_join(bindermain, gridmain, by="grid_id")


# Supp. 4 N to m/s calculation ----

# Calculate rubble piece measurements ----

pieceid$experiment <- as.factor(pieceid$experiment)
piecemain <- pieceid %>% filter(experiment == "main") %>% droplevels()

str(piecemain)

# Join gridid and pieceid so that the rubble length information is availabile

gridpiecemain <- inner_join(gridmain, piecemain, by="grid_id")
write.csv(gridpiecemain, file = "gridpiecemain2.csv") # checked this

piecemain <- piecemain %>% mutate(avgdiameterperpiece = (width_1 + width_2 + width_3)/3)

pieceavg <- piecemain %>% group_by(unique_pair_id) %>%
  summarise(avgrubblelength= mean(axial_length),
            maxrubblelength= max(axial_length),
            avgrubblediameter = mean(avgdiameterperpiece),
            maxrubblediameter = max(avgdiameterperpiece)) %>% as.data.frame() 

write.csv(pieceavg, "pieceavg2.csv")

range(pieceavg$avgrubblediameter, na.rm = TRUE) # diameter from 0.8 to 2.5, seems right
range(pieceavg$avgrubblelength, na.rm = TRUE) # length of 4 to 19, looks right

# Calculating bind measurements ----

names(bindermain)

binderidsum <- bindermain %>% group_by(unique_pair_id) %>%
  summarise(pairbindareacm_LT_org = sum(bind_length_x_thickness)) %>%  as.data.frame() # this is the chosen metric for calculating bind area

# now join this with the pairmain sheet  (which is just pairid for the main experiment only)

pairbindmain <- inner_join(pairmain, binderidsum, by="unique_pair_id")

head(pieceavg)

# Add in the rubble piece measurements

pairbindpiecemain <- inner_join(pairbindmain, pieceavg, by="unique_pair_id")

# Join this with the site data

gridpairbindpiecemain <- inner_join(gridmain, pairbindpiecemain, by="grid_id")

# Formatting

gridpairbindpiecemain <- gridpairbindpiecemain %>% mutate(boundB = ifelse(gridpairbindpiecemain$bound == "Yes", 1, 0))

gridpairbindpiecemain$layer <- factor(gridpairbindpiecemain$layer ,
                                      levels = c('Top', 'Bottom'),
                                      labels = c('Surface grids', 'Buried grids'))
gridpairbindpiecemain$location <- factor(gridpairbindpiecemain$location,
                                         levels = c('Keppels', 'Heron'),
                                         labels = c('Keppels (inshore)', 'Heron (offshore)'))

gridpairbindpiecemain$exposuredepth <- paste(gridpairbindpiecemain$exposure, gridpairbindpiecemain$depth, sep = "_")
gridpairbindpiecemain$exposuredepth <- as.factor(gridpairbindpiecemain$exposuredepth)

gridpairbindpiecemain$exposuredepth <- factor(gridpairbindpiecemain$exposuredepth, levels = c("Low_Reef flat",
                                                                                              "Low_Shallow slope",
                                                                                              "Low_Deep slope",
                                                                                              "High_Shallow slope", 
                                                                                              "High_Deep slope"),
                                              labels = c("Sheltered Reef Flat", "Sheltered Shallow",
                                                         "Sheltered Deep", "Exposed Shallow", "Exposed Deep"))

gridpairbindpiecemain$locexpdepth <- paste(gridpairbindpiecemain$location, gridpairbindpiecemain$exposuredepth, sep = "_")
levels(gridpairbindpiecemain$locexpdepth)
gridpairbindpiecemain$locexpdepth <- as.factor(gridpairbindpiecemain$locexpdepth)

gridpairbindpiecemain$locexpdepth <- factor(gridpairbindpiecemain$locexpdepth, levels = c("Keppels (inshore)_Sheltered Shallow",
                                                                                          "Keppels (inshore)_Exposed Shallow",
                                                                                          "Heron (offshore)_Sheltered Reef Flat",
                                                                                          "Heron (offshore)_Sheltered Shallow" , 
                                                                                          "Heron (offshore)_Sheltered Deep",
                                                                                          "Heron (offshore)_Exposed Shallow",
                                                                                          "Heron (offshore)_Exposed Deep"),
                                            labels = c("In Shelt Shal", "In Exp Shal", "Off Shelt RF", "Off Shelt Shal",
                                                       "Off Shelt Deep", "Off Exp Shal", "Off Exp Deep"))

levels(gridpairbindpiecemain$exposuredepth)
levels(gridpairbindpiecemain$location)
levels(gridpairbindpiecemain$layer)
levels(gridpairbindpiecemain$bound)

# Now, we only want break force to be 0 if the binding was 'No', if break force is NA
# for another reason (pair missing, couldn't measure properly) then we want break force to be NA still
gridpairbindpiecemain2 <- gridpairbindpiecemain %>%
  mutate(break_force = ifelse(bound == 'No' & is.na(break_force), 0, break_force))

write.csv(gridpairbindpiecemain2, "gridpairbindpiecemain2_3.csv") # checked v 2_2 for pair 1704 as example

gridpairbindpiecemain3 <- gridpairbindpiecemain2 %>% 
  mutate(Tbind_Nm_LT_org = break_force/(pairbindareacm_LT_org/10000),
         Tbind_Ncm_LT_org = break_force/(pairbindareacm_LT_org)) %>% as.data.frame() # sum of all bind areas on pair (length x thickness); Tom Baldock advised to use this

# Change the NAs to 0 if there is no binding

gridpairbindpiecemain3 <- gridpairbindpiecemain3 %>%
  mutate(Tbind_Nm_LT_org = ifelse(bound == 'No' & is.na(Tbind_Nm_LT_org), 0, Tbind_Nm_LT_org),
         Tbind_Ncm_LT_org = ifelse(bound == 'No' & is.na(Tbind_Ncm_LT_org), 0, Tbind_Ncm_LT_org))

# Calculate the velocity based on Tau, bind area, and rubble dimensions
# Variables to use in the break velocity equation are below.
# Tau here is already in N/m2, but bind area, rubble diameter 
# and length are in ‘cm’ so are being adjusted to ‘m’ here.
# alpha is 1-gamma - most conservative where the rubble piece's diameter is fully exposed, is alpha = 1 (because gamma = 0).
# theta is cos(degrees), most conservative scenario where flow is perpendicular to rubble peice, is cos(0) = 1.

gridpairbindpiecemain4 <- gridpairbindpiecemain3 %>% 
  mutate(break_velocity_BINDAREA_LT_org = sqrt((2*(Tbind_Nm_LT_org)*(pairbindareacm_LT_org/10000))/
                                                 (1020*1*((avgrubblediameter/100)*
                                                            (avgrubblelength/100)))),
         break_velocity_NOarea_Cons = sqrt((2*break_force)/
                                             (1020*1*((1*(avgrubblelength/100))*
                                                        (1*(avgrubblediameter/100))))),# testing for theta/orientation of cos(0)=1 (straight on) and alpha of 1 (no sheltering, meaning gamma = 0)
         break_velocity_NOarea_Avg = sqrt((2*break_force)/
                                            (1020*1*((0.7071067812*(avgrubblelength/100))*
                                                       ((1-habitat_meangamma)*(avgrubblediameter/100)))))) %>%  # this one uses average gamma for each habitat, and 45 degrees for theta
  as.data.frame() 


# Now, again we only want break force to be 0 if the binding was 'No', if break force is NA
# for another reason (pair missing, couldn't measure properly) then we want break force to be NA still
gridpairbindpiecemain5 <- gridpairbindpiecemain4 %>%
  mutate(break_velocity_BINDAREA_LT_org = ifelse(bound == 'No' & is.na(break_velocity_BINDAREA_LT_org), 
                                                 0, break_velocity_BINDAREA_LT_org),
         break_velocity_NOarea_Cons = ifelse(bound == 'No' & is.na(break_velocity_NOarea_Cons), 
                                             0, break_velocity_NOarea_Cons),
         break_velocity_NOarea_Avg = ifelse(bound == 'No' & is.na(break_velocity_NOarea_Avg), 
                                            0, break_velocity_NOarea_Avg))

# Remove 'rolled' tests as we don't have reliable strength data for them (can estimate but it only removes 3 reps so not a huge data loss)

gridpairbindpiecemainS <- gridpairbindpiecemain5 %>% 
  filter(force_type != "rolled") %>% droplevels()

# Display the levels of the modified factor
levels(gridpairbindpiecemainS$force_type) # Rolled has been dropped

# Remember that grids 297, 298, and 299 are not included (extra reef flat grids), and that rolled force is dropped

write.csv(gridpairbindpiecemainS, 'gridpairbindpiecemainS_29May.csv') # 1892 pairs in the surface grids dataset
# Checked this dataframe on 11 April 2024 before analyses continued
# Note: Use gridpairbindpiecemainS2 (created later) if you want the updated total_binds column but also want both layers.

# Total binds ----
# Making sure the total binds column is right (creating it from the bindermain sheet):

bindermain2 <- bindermain %>%
  group_by(unique_pair_id) %>%
  mutate(
    total_binds2 = ifelse(any(is.na(binder_num)) & any(is.na(binder_DNA_broadest_cat)), NA, 
                          ifelse(any(is.na(binder_num)) | any(binder_num > 0), n(), 0))) %>%
  distinct(unique_pair_id, .keep_all = TRUE)

# total_binds2 will have either 0 if # binds are O, or the count of rows for each unique_pair_id, if binder_num is NA or a number greater than 0
# and the ones where pairs are NA, i.e, binder_num is NA and category is NA, need to be NA in total_binds2 too.

# Now join with gridpairbindpiecemainS

gridpairbindpiecemainS2 <- inner_join(gridpairbindpiecemainS, bindermain2, by="unique_pair_id")

gridpairbindpiecemainSS2 <- gridpairbindpiecemainS2 %>% 
  filter(layer == "Surface grids") %>% droplevels()

str(gridpairbindpiecemainS)
str(gridpairbindpiecemainSS2)

range(gridpairbindpiecemainSS2$months_deployed)
gridpairbindpiecemainSS2$months_deployedFO <- ordered(gridpairbindpiecemainSS2$months_deployed)
levels(gridpairbindpiecemainSS2$months_deployedFO)

gridpairbindpiecemainSS2$months_deployedF <- factor(gridpairbindpiecemainSS2$months_deployed)
levels(gridpairbindpiecemainSS2$months_deployedF)

# Adding in Gap and Contact ID ----

# To investigate whether we need to standardise the number of binds, and strength, to
# And then for bind strength regression, also use bind area for each binder as well as number of binds for each binder as you've done

# Contact ID

# Delete rows where number is NA first:
names(contactid)
contactid <- contactid[,-c(7,8,9)] # removing 7th to 9th columns, keeping all rows
contactid <- contactid %>%
  filter(!is.na(number))

str(contactid)
contactid$side2 <- as.factor(contactid$side2)
contactid$contact_length <- as.numeric(contactid$contact_length)

# Also not doing the subtraction if both or one of the values is NA for contact_length.
contactid2 <- contactid %>%
  group_by(unique_pair_id, number) %>%
  mutate(same_as_other_side = case_when(
    !is.na(contact_length[side2 == "topside"]) & !is.na(contact_length[side2 == "underside"]) ~ 
      ifelse(abs(contact_length[side2 == "topside"] - contact_length[side2 == "underside"]) < 0.5, "yes", "no"),
    TRUE ~ NA_character_)) %>%    # If any of the values are NA, return NA
  ungroup()

write.csv(contactid2, "contactid3.csv") # Checked pair 2829

# OK, let's now get ONE value per pair, rather than topside and underside

contactid2 <- contactid2 %>%
  filter(!is.na(number))

contactid3 <- contactid2 %>%
  group_by(unique_pair_id, number) %>%
  summarise(contact_length_topunder = max(contact_length[side2 == "topside"],  # note, this is taking the max value between topside 1 and underside 1 and pasting it just once (because often it's the same gap/contact point)
                                          contact_length[side2 == "underside"], # you could take the average instead but I think max is fine
                                          na.rm = TRUE), 
            .groups = "drop") %>%  # Remove the grouping
  mutate(contact_length_topunder = ifelse(contact_length_topunder == -Inf, NA, contact_length_topunder)) %>%  # Replace -Inf with NA
  filter(!is.na(contact_length_topunder))  # Remove rows where both were NA

# Use the mutate line if you want the NAs to stay in. 
# Not every pair is there now. But if no gap/contact value, then we can't standardise it... so, OK to not be there

view(contactid3) # there are warnings but only because of the -Inf that are created, but have replaced those with NA and then removed.


# Do the same for gapid:

gapid$side2 <- as.factor(gapid$side2)
gapid$gap_length <- as.numeric(gapid$gap_length)
gapid$max_gap_span <- as.numeric(gapid$max_gap_span)
str(gapid)

# Also not doing the subtraction if both or one of the values is NA for contact_length.
# Check what this does
gapid2 <- gapid %>%
  group_by(unique_pair_id, number) %>%
  mutate(same_as_other_side = case_when(
    !is.na(gap_length[side2 == "topside"]) & !is.na(gap_length[side2 == "underside"]) ~ 
      ifelse(abs(gap_length[side2 == "topside"] - gap_length[side2 == "underside"]) < 0.5, "yes", "no"),
    TRUE ~ NA_character_)) %>%    # If any of the values are NA, return NA
  ungroup()

write.csv(gapid2, "gapid3.csv")

# OK, let's now get ONE value per pair, rather than topside and underside

gapid2 <- gapid2 %>%
  filter(!is.na(number))

gapid3 <- gapid2 %>%
  group_by(unique_pair_id, number) %>%
  summarise(gap_length_topunder = max(gap_length[side2 == "topside"],  # note, this is taking the max value between topside 1 and underside 1 and pasting it just once (because often it's the same gap/contact point)
                                      gap_length[side2 == "underside"], # you could take the average instead but I think max is fine
                                      na.rm = TRUE), 
            .groups = "drop") %>%  # Remove the grouping
  mutate(gap_length_topunder = ifelse(gap_length_topunder == -Inf, NA, gap_length_topunder)) %>%  # Replace -Inf with NA
  filter(!is.na(gap_length_topunder))  # Remove rows where both were NA
#view(gapid3) # there are warnings but because of the -Inf that are created, but have replaced those with NA and then removed.

# Add in the max_gap_span - which is per pair, not per gap.

gapid2_2 <- gapid2 %>%
  group_by(unique_pair_id) %>%
  summarise(max_gap_span = max(max_gap_span))  

#view(gapid2_2)

gapid4 <- gapid3 %>%
  left_join(gapid2_2 %>% select(unique_pair_id, max_gap_span), by = "unique_pair_id")

# View the resulting dataframe
#view(gapid4)

unique(gapid$unique_pair_id)
unique(contactid$unique_pair_id)
unique(gapid4$unique_pair_id) # 3222 pairs
unique(contactid3$unique_pair_id) # 3222 pairs

# Now add the gap and contact lengths together to get a 'binding opportunity length'

head(contactid3)
head(gapid4)

gapcontact <- full_join(contactid3, gapid4, by = c("unique_pair_id", "number"))
head(gapcontact)

gapcontactSUM <- gapcontact %>% group_by(unique_pair_id) %>%
  summarise(sumcontact = sum(contact_length_topunder, na.rm=TRUE),
            sumgap = sum(gap_length_topunder, na.rm=TRUE))

# Add max gap span

gapcontact_2 <- gapcontact %>%
  group_by(unique_pair_id) %>%
  summarise(max_gap_span = max(max_gap_span))

gapcontactSUM2 <- gapcontactSUM %>%
  left_join(gapcontact_2 %>% select(unique_pair_id, max_gap_span), by = "unique_pair_id")

#view(gapcontactSUM2)

# Also create a dataframe that tells  the number of contact points and number of gaps:

gapcontactCOUNT <- gapcontact %>%
  group_by(unique_pair_id) %>%
  summarise(
    num_contacts = sum(!is.na(contact_length_topunder)),  # Count of non-NA values in contact_length_topunder
    num_gaps = sum(!is.na(gap_length_topunder)))       # Count of non-NA values in gap_length_topunder

head(gapcontactCOUNT)


# Now join this with the sum dataframe.

gapcontactFULL <- full_join(gapcontactSUM2, gapcontactCOUNT, by = c("unique_pair_id"))

#view(gapcontactFULL) # this one has 3197 pairs, down from 3222, so 25 are NA for all of the gap/contact measurements


# Dataframe with gaps/contacts for looking at TOTAL BINDS and BIND STRENGTH per habitat:

# MAIN EXPERIMENT

head(gridpairbindpiecemainSS2) # gridpairbindpiecemainSS2 has updated total_binds column but only the surface layer.
# gridpairbindpiecemainS2 if you want the updated total_binds column but also want both layers.
head(gapcontactFULL)

gridpairbindpiecemainSS2_gc <- left_join(gridpairbindpiecemainSS2, 
                                         gapcontactFULL, by = "unique_pair_id") 

names(gridpairbindpiecemainSS2_gc)

write.csv(gridpairbindpiecemainSS2_gc, "gridpairbindpiecemainSS2_gcCHECK.csv")


# Adding in cumulative bind length ----

gridbinderB_cum2 <- gridbinder %>%
  group_by(unique_pair_id, grid_id, 
           location, site_unique, depth, months_deployed, 
           exposure, layer) %>% 
  summarise(bind_length_cum = sum(bind_length, na.rm = TRUE)) %>%
  mutate(bind_length_cum = ifelse(is.na(bind_length_cum), 0, bind_length_cum)) %>%
  ungroup()

gridpairbindpiecemainSS2_gc <- gridpairbindpiecemainSS2_gc %>%
  left_join(gridbinderB_cum2 %>% select(unique_pair_id, bind_length_cum), by = "unique_pair_id")

str(gridpairbindpiecemainSS2_gc) # this has bind area, based on length x thickness which is used to calculate velocity, and cumbindlength

# Plot cumulative bind length vs break force

p1 <- ggplot(gridpairbindpiecemainSS2_gc, aes(x = bind_length_cum, y = break_force)) +
  geom_point() + geom_smooth() + ggtitle("All pairs") + xlab ("Cumulative bind length per pair")

p2 <- ggplot(gridpairbindpiecemainSS2_gc, aes(x = pairbindareacm_LT_org, y = break_force)) +
  geom_point() + geom_smooth() + ggtitle("All pairs") + xlab ("Cumulative bind area per pair")


grid.arrange(p1,p2)


str(gridpairbindpiecemainSS2_gc)

write.csv(gridpairbindpiecemainSS2_gc, "gridpairbindpiecemainSS2_gc.csv")

















# CROSS EXPERIMENT DATAFRAME ----

# gridid & pairid

gridcross <- gridid %>% filter(experiment == "cross") %>% droplevels()
paircross <- pairid %>% filter(experiment == "cross") %>% droplevels()

str(gridcross)
str(paircross)

gridpaircross <- inner_join(gridcross, paircross, by="grid_id")

write.csv(gridpaircross, file = "gridpaircross.csv") # checked this

# gridid & binderid

bindercross <- binderid %>% filter(experiment == "cross") %>% droplevels()
gridbindercross <- inner_join(bindercross, gridcross, by="grid_id")

#view(gridbindercross)

# N to m/s calculation ----

# Calculate rubble piece measurements ----

pieceid$experiment <- as.factor(pieceid$experiment)
piececross <- pieceid %>% filter(experiment == "cross") %>% droplevels()

str(piececross)

# Join gridid and pieceid so that we have the rubble length information.

gridpiececross <- inner_join(gridcross, piececross, by="grid_id")
write.csv(gridpiececross, file = "gridpiececross.csv") # checked this

piececross <- piececross %>% mutate(avgdiameterperpiece = (width_1 + width_2 + width_3)/3)

piececrossavg <- piececross %>% group_by(unique_pair_id) %>%
  summarise(avgrubblelength= mean(axial_length),
            maxrubblelength= max(axial_length),
            avgrubblediameter = mean(avgdiameterperpiece),
            maxrubblediameter = max(avgdiameterperpiece)) %>% as.data.frame() # don't need average circumference anymore
#View(piececrossavg)
write.csv(piececrossavg, "piececrossavg.csv")

range(piececrossavg$avgrubblediameter, na.rm = TRUE) # diameter from 1 to 1.8, seems right
range(piececrossavg$avgrubblelength, na.rm = TRUE) # length of 6 to 16, looks right

# Calculate bind measurements ----

view(bindercross)

bindercrosssum <- bindercross %>% group_by(unique_pair_id) %>%
  summarise(pairbindareacm_LT_org = sum(bind_length_x_thickness)) %>%  as.data.frame() # this is the chosen metric after the comparisons done (sent to Tom/Fikri)

# now join this with the paircross sheet  (which is just pairid for the main experiment only)

pairbindcross <- inner_join(paircross, bindercrosssum, by="unique_pair_id")
#View(pairbindcross)

head(pairbindcross)

# Add in the rubble piece measurements (the average of the pair, cannot consider the lengths of each piece separately when doing the N to m/s conversion)

pairbindpiececross <- inner_join(pairbindcross, piececrossavg, by="unique_pair_id")

# Join this with the site data

gridpairbindpiececross <- inner_join(gridcross, pairbindpiececross, by="grid_id")
head(gridpairbindpiececross)

# Formatting

levels(gridpairbindpiececross$bound)

gridpairbindpiececross$bound <- factor(gridpairbindpiececross$bound, levels = c("No","No ","Yes"),
                                       labels = c("No", "No", "Yes"))

gridpairbindpiececross <- gridpairbindpiececross %>% mutate(boundB = ifelse(gridpairbindpiececross$bound == "Yes", 1, 0))

gridpairbindpiececross$exposuredepth <- paste(gridpairbindpiececross$exposure, gridpairbindpiececross$depth, sep = "_")
gridpairbindpiececross$exposuredepth <- as.factor(gridpairbindpiececross$exposuredepth)

levels(gridpairbindpiececross$exposuredepth)

gridpairbindpiececross$exposuredepth <- factor(gridpairbindpiececross$exposuredepth, levels = c("Low_Shallow slope",
                                                                                                "Low_Deep slope"),
                                               labels = c("Sheltered Shallow",
                                                          "Sheltered Deep" ))

gridpairbindpiececross$locexpdepth <- paste(gridpairbindpiececross$location, gridpairbindpiececross$exposuredepth, sep = "_")

gridpairbindpiececross$locexpdepth <- as.factor(gridpairbindpiececross$locexpdepth)
levels(gridpairbindpiececross$locexpdepth)

gridpairbindpiececross$locexpdepth <- factor(gridpairbindpiececross$locexpdepth, levels = c("Heron_Sheltered Shallow" , 
                                                                                            "Heron_Sheltered Deep"),
                                             labels = c("Off Shelt Shal",
                                                        "Off Shelt Deep"))

# Now, we only want break force to be 0 if the binding was 'No', if break force is NA
# for another reason (pair missing, couldn't measure properly) then we want break force to be NA still
gridpairbindpiececross2 <- gridpairbindpiececross %>%
  mutate(break_force = ifelse(bound == 'No' & is.na(break_force), 0, break_force))

write.csv(gridpairbindpiececross2, "gridpairbindpiececross2.csv")

gridpairbindpiececross3 <- gridpairbindpiececross2 %>% 
  mutate(Tbind_Nm_LT_org = break_force/(pairbindareacm_LT_org/10000),
         Tbind_Ncm_LT_org = break_force/(pairbindareacm_LT_org)) %>% as.data.frame() # sum of all bind areas on pair (length x thickness), this is what Tom advised to use

#View(gridpairbindpiececross3)

# Change the NAs to 0 if there is no binding

gridpairbindpiececross4 <- gridpairbindpiececross3 %>%
  mutate(Tbind_Nm_LT_org = ifelse(bound == 'No' & is.na(Tbind_Nm_LT_org), 0, Tbind_Nm_LT_org),
         Tbind_Ncm_LT_org = ifelse(bound == 'No' & is.na(Tbind_Ncm_LT_org), 0, Tbind_Ncm_LT_org))

# Calculate the velocity based on Tau, bind area, and rubble dimensions
# Variables to use in the break velocity equation are below.
# Tau here is already in N/m2, but bind area, rubble diameter 
# and length are in ‘cm’ so are being adjusted to ‘m’ here.
# alpha is 1-gamma - most conservative is alpha = 1 (because gamma = 0).
# theta is cos(degrees), most conservative scenario is cos(0) = 1.

gridpairbindpiececross5 <- gridpairbindpiececross4 %>% 
  mutate(break_velocity_BINDAREA_LT_org = sqrt((2*(Tbind_Nm_LT_org)*(pairbindareacm_LT_org/10000))/
                                                 (1020*1*((avgrubblediameter/100)*
                                                            (avgrubblelength/100)))),
         break_velocity_NOarea_Cons = sqrt((2*break_force)/
                                             (1020*1*((1*(avgrubblelength/100))*
                                                        (1*(avgrubblediameter/100))))),# testing for theta/orientation of cos(0)=1 (straight on) and alpha of 1 (no sheltering, meaning gamma = 0)
         break_velocity_NOarea_Avg = sqrt((2*break_force)/
                                            (1020*1*((0.7071067812*(avgrubblelength/100))*
                                                       ((1-habitat_meangamma)*(avgrubblediameter/100)))))) %>%  # this one uses average gamma for each habitat, and 45 degrees for theta
  as.data.frame() 


# Now, again we only want break force to be 0 if the binding was 'No', if break force is NA
# for another reason (pair missing, couldn't measure properly) then we want break force to be NA still
gridpairbindpiececross6 <- gridpairbindpiececross5 %>%
  mutate(break_velocity_BINDAREA_LT_org = ifelse(bound == 'No' & is.na(break_velocity_BINDAREA_LT_org), 
                                                 0, break_velocity_BINDAREA_LT_org),
         break_velocity_NOarea_Cons = ifelse(bound == 'No' & is.na(break_velocity_NOarea_Cons), 
                                             0, break_velocity_NOarea_Cons),
         break_velocity_NOarea_Avg = ifelse(bound == 'No' & is.na(break_velocity_NOarea_Avg), 
                                            0, break_velocity_NOarea_Avg))

levels(gridpairbindpiececross6$force_type)

# Remove 'rolled' tests 
# Cutting out 6 of 30 normal pairs in shallow and 1 of 30 normal pairs in deep. Ask Wen for velocity estimates and add in.

gridpairbindpiececrossS <- gridpairbindpiececross6 %>% 
  filter(force_type != "rolled") %>% droplevels()

# Display the levels of the modified factor
levels(gridpairbindpiececrossS$force_type) # Rolled has been dropped

write.csv(gridpairbindpiececrossS, 'gridpairbindpiececrossS.csv') 

# Total binds
# Making sure the total binds column is right (creating it from the bindercross sheet):

bindercross2 <- bindercross %>%
  group_by(unique_pair_id) %>%
  mutate(
    total_binds2 = ifelse(any(is.na(binder_num)) & any(is.na(binder_DNA_broadest_cat)), NA, 
                          ifelse(any(is.na(binder_num)) | any(binder_num > 0), n(), 0))) %>%
  distinct(unique_pair_id, .keep_all = TRUE)

# total_binds2 will have either 0 if binds is O, or the count of rows for each unique_pair_id, if binder_num is NA or a number greater than 0
# and the ones where pairs are NA, i.e, binder_num is NA and category is NA, need to be NA in total_binds2 too.

#view(bindercross2)

# Now join with gridpairbindpiececrossS

gridpairbindpiececrossS2 <- inner_join(gridpairbindpiececrossS, bindercross2, by="unique_pair_id")
# There are 173 rows meaning 173 pairs (180 minus the rolled ones) 10*3*3*2

str(gridpairbindpiececrossS2)

gridpairbindpiececrossS2 <- gridpairbindpiececrossS2 %>% 
  mutate(pair_arrangement = factor(pair_arrangement, levels = c("parallel",   "crossed", "parallel_loose"),
                                   labels = c("Parallel_Stable", 
                                              "Crossed_Stable",
                                              "Parallel_Loose")))
# Now add gapid and contactid information in, as we did with the main experiment

gridpairbindpiececrossS2_gc <- left_join(gridpairbindpiececrossS2, 
                                         gapcontactFULL, by = "unique_pair_id") 

str(gridpairbindpiececrossS2_gc)

gridpairbindpiececrossS2_gc$months_deployedFO <- as.ordered(gridpairbindpiececrossS2_gc$months_deployed)



























# Supp. 2 Average rubble size -----

# Probability density curves

legend_titleT <- "Habitat"

colsT <- c("In Shelt Shal" = "#6CA6CD",
           "In Exp Shal" = "#CD8C95", 
           "Off Shelt RF" = "#AEE2F5",
           "Off Shelt Shal" = "#6CA6CD",
           "Off Shelt Deep" = "#104E8B",
           "Off Exp Shal" = "#CD8C95",    
           "Off Exp Deep" = "#A52A2A")

linesT <- c("In Shelt Shal" = 2,
            "In Exp Shal" = 2, 
            "Off Shelt RF" = 1,
            "Off Shelt Shal" = 1,
            "Off Shelt Deep" = 1,
            "Off Exp Shal" = 1,    
            "Off Exp Deep" = 1)

# Rubble size for main experiment:
M1 <- ggplot(gridpairbindpiecemainSS2_gc, 
             aes(x=avgrubblelength, colour = locexpdepth,
                 linetype = locexpdepth)) + 
  geom_density(size = 0.8) + 
  geom_vline(aes(xintercept=mean(avgrubblelength, na.rm = TRUE)),
             color="black", linetype="dashed", size=0.5) +
  #geom_hline(yintercept = seq(0, 20, by = 4), linetype = "dotted", color = "black", size = 0.3) +
  scale_colour_manual(legend_titleT, values=colsT) +
  scale_linetype_manual(legend_titleT, values = linesT) +
  theme_classic() + 
  scale_x_continuous(limits = c(4,20),breaks = c(3,6,9,12,15,18)) +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +   # Bold the legend title
  theme(legend.position="bottom") +
  ylab("Probability density") + 
  xlab("Average rubble length per pair")

M1

# Rubble size for cross pairs:
C1 <- ggplot(gridpairbindpiececrossS2_gc, 
             aes(x=avgrubblelength, colour = locexpdepth,
                 linetype = locexpdepth)) + 
  geom_density(size = 0.8) + 
  geom_vline(aes(xintercept=mean(avgrubblelength, na.rm = TRUE)),
             color="black", linetype="dashed", size=0.5) +
  #geom_hline(yintercept = seq(0, 20, by = 4), linetype = "dotted", color = "black", size = 0.3) +
  scale_colour_manual(legend_titleT, values=colsT) +
  scale_linetype_manual(legend_titleT, values = linesT) +
  scale_x_continuous(limits = c(4,20),breaks = c(3,6,9,12,15,18)) +
  theme_classic() + 
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +   # Bold the legend title
  theme(legend.position="bottom") +
  ylab("Probability density") + 
  xlab("Average rubble length per pair")

C1

grid.arrange(M1,C1)

# Does rubble length vary between sites?

names(gridpair)
names(piecemain)

gridpiecemain <- inner_join(piecemain, gridpair, by = "unique_pair_id")

names(gridpiecemain)

# Formatting

gridpiecemain$layer <- factor(gridpiecemain$layer ,
                              levels = c('Top', 'Bottom'),
                              labels = c('Surface grids', 'Buried grids'))
gridpiecemain$location <- factor(gridpiecemain$location,
                                 levels = c('Keppels', 'Heron'),
                                 labels = c('Keppels (inshore)', 'Heron (offshore)'))

gridpiecemain$exposuredepth <- paste(gridpiecemain$exposure, gridpiecemain$depth, sep = "_")
gridpiecemain$exposuredepth <- as.factor(gridpiecemain$exposuredepth)

gridpiecemain$exposuredepth <- factor(gridpiecemain$exposuredepth, levels = c("Low_Reef flat",
                                                                              "Low_Shallow slope",
                                                                              "Low_Deep slope",
                                                                              "High_Shallow slope", 
                                                                              "High_Deep slope"),
                                      labels = c("Sheltered Reef Flat", "Sheltered Shallow",
                                                 "Sheltered Deep", "Exposed Shallow", "Exposed Deep"))

gridpiecemain$locexpdepth <- paste(gridpiecemain$location, gridpiecemain$exposuredepth, sep = "_")
levels(gridpiecemain$locexpdepth)
gridpiecemain$locexpdepth <- as.factor(gridpiecemain$locexpdepth)

gridpiecemain$locexpdepth <- factor(gridpiecemain$locexpdepth, levels = c("Keppels (inshore)_Sheltered Shallow",
                                                                          "Keppels (inshore)_Exposed Shallow",
                                                                          "Heron (offshore)_Sheltered Reef Flat",
                                                                          "Heron (offshore)_Sheltered Shallow" , 
                                                                          "Heron (offshore)_Sheltered Deep",
                                                                          "Heron (offshore)_Exposed Shallow",
                                                                          "Heron (offshore)_Exposed Deep"),
                                    labels = c("In Shelt Shal", "In Exp Shal", "Off Shelt RF", "Off Shelt Shal",
                                               "Off Shelt Deep", "Off Exp Shal", "Off Exp Deep"))

str(gridpiecemain)
gridpiecemainS <- gridpiecemain %>% filter(layer == "Surface grids") %>% droplevels()

# Model:

range(gridpiecemainS$axial_length, na.rm = TRUE)
mean(gridpiecemainS$axial_length, na.rm = TRUE)
sd(gridpiecemainS$axial_length, na.rm = TRUE)

hist(gridpiecemainS$axial_length) # pretty normal

len1_1 <- glmmTMB(axial_length ~ site_unique + (1|grid_id.y), data = gridpiecemainS)

plot(simulateResiduals(fittedModel = len1_1))

len1_1_comps <- emmeans(len1_1, pairwise ~ site_unique, type = "response")

save_tables_to_docx(len1_1_comps, "len1_1_comps1.docx", "len1_1_comps2.docx")


















# Supp. 3 Tensile vs shear testing -----


# Which grid ids have both tensile and shear testing done?

# Check levels of grid.id.x with both tensile and shear
result <- gridpairbindpiecemainSS2_gc %>%
  filter(!is.na(force_type)) %>%
  group_by(grid_id.x) %>%
  summarize(force_count = n_distinct(force_type)) %>%
  filter(force_count > 1)

View(result)

filtered_ids <- result$grid_id.x

# Select rows from the original dataset where grid.id.x is in filtered_ids
filtered_sheardat <- gridpairbindpiecemainSS2_gc %>%
  filter(grid_id.x %in% filtered_ids)

# Display the filtered rows
view(filtered_sheardat)

ggplot(filtered_sheardat, aes(x = force_type, y = break_force)) + geom_point()

# Model it: 

# This is for 45 grids and 447 pairs
# We did 118 shear and 329 tensile (average of 3 pairs per grid shear)

filtered_sheardatSUM <- filtered_sheardat %>%
  group_by(grid_id.x, force_type) %>%
  summarise(count = n(), .groups = "drop")
filtered_sheardatSUM

filtered_sheardatSUM2 <- filtered_sheardat %>%
  group_by(force_type) %>%
  summarise(count = n(), .groups = "drop")
filtered_sheardatSUM2

hist(filtered_sheardat$break_force)

shearmod <- glmmTMB(break_force ~ force_type +
                      (1|site_unique/grid_id.y), family = ziGamma(link= "log"),
                    ziformula = ~1,
                    data = filtered_sheardat)
car::Anova(shearmod) # no difference between shear and tensile in terms of the resulting break force
# Also checked for interaction with habitat and there was none.

plot(simulateResiduals(shearmod)) # good

# What about with grid as a fixed effect, is it different for any of the grids? 

filtered_sheardat$grid_id.x <- as.factor(filtered_sheardat$grid_id.x)

shearmod2 <- glmmTMB(break_force ~ force_type * grid_id.x +
                       (1|site_unique), family = ziGamma(link= "log"),
                     ziformula = ~1,
                     data = filtered_sheardat)

car::Anova(shearmod2) # no interaction between force type and grid


emmeans(shearmod2, pairwise ~ force_type | grid_id.x, type ="response")
# Plot the grid-force_type model:

shearmodgrid <- with(filtered_sheardat, list(force_type = levels(force_type),
                                             grid_id.x = levels(grid_id.x)))

shearmodpreds <- emmeans(shearmod2, ~ force_type | grid_id.x, at = shearmodgrid, 
                         type = "response") %>% as.data.frame() # error

pd <- position_dodge(0.5)

ggplot() + 
  geom_pointrange(aes(x = grid_id.x, y=response, 
                      colour = force_type,
                      fill = force_type,
                      ymin = asymp.LCL, 
                      ymax = asymp.UCL),
                  data = shearmodpreds,
                  position=pd, size =0.8) + # predicted average data points +
  geom_errorbar(aes(y=response, x=grid_id.x, 
                    ymin=asymp.LCL, ymax=asymp.UCL, 
                    colour = force_type), 
                data = shearmodpreds, width=0.3, position=pd) +
  geom_point(aes(y=break_force, x=grid_id.x, colour = force_type), 
             data = filtered_sheardat, 
             position = position_jitter(width=0.2, height=0), 
             size =0.5, alpha= 0.1) + # raw data
  geom_point() +
  geom_hline(yintercept = seq(0, 100, by = 20), linetype = "dotted", color = "black", size = 0.3) +
  theme_classic() + 
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +   # Bold the legend title
  ylab("Force to break binds (N)") + 
  xlab("Grid ID") + 
  coord_cartesian(ylim=c(0,100))






# Supp. 7 Bind density across time and space -----


freqmod2TB <- glmmTMB(total_binds2 ~ months_deployedFO * locexpdepth + 
                        (1|site_unique/grid_id.y), family=poisson, 
                      ziformula=~1, data = gridpairbindpiecemainSS2_gc) 

car::Anova(freqmod2TB) # sig interaction

#Plotting from freqmod2TB

freqmod2TBcomps <- emmeans(freqmod2TB, pairwise ~ locexpdepth | months_deployedFO, type = "response")

freqmod2TBcomps2 <- emmeans(freqmod2TB, pairwise ~  months_deployedFO | locexpdepth, type = "response")

# Add this table to Supplementary Material

save_tables_to_docx(freqmod2TBcomps, "freqmod2TBcomps1.docx", "freqmod2TBcomps2.docx")
save_tables_to_docx(freqmod2TBcomps2, "freqmod2TBcomps21.docx", "freqmod2TBcomps22.docx")

freqmod2TB_grid <- with(gridpairbindpiecemainSS2, 
                        list(locexpdepth = levels(locexpdepth), 
                             months_deployedFO = levels(months_deployedFO)))

freqmod2TB_preds <- emmeans(freqmod2TB, ~ locexpdepth | months_deployedFO, 
                            at = freqmod2TB_grid, type = "response") %>% as.data.frame()

head(freqmod2TB_preds) 

# Plot it:

shapesT <- c("In Shelt Shal" = 8,
             "In Exp Shal" = 8, 
             "Off Shelt RF" = 17,
             "Off Shelt Shal" = 17,
             "Off Shelt Deep" = 17,
             "Off Exp Shal" = 17,    
             "Off Exp Deep" = 17)

(freqmod2TB_plot <- ggplot() +
    geom_pointrange(aes(y=rate, x=months_deployedFO, 
                        colour = locexpdepth,
                        fill = locexpdepth,
                        shape = locexpdepth,
                        ymin = asymp.LCL, 
                        ymax = asymp.UCL),
                    data = freqmod2TB_preds,
                    position=pd, size =0.8) + # predicted average data points +
    geom_errorbar(aes(y=rate, x=months_deployedFO, 
                      ymin=asymp.LCL, ymax=asymp.UCL, 
                      colour = locexpdepth), 
                  data = freqmod2TB_preds, width=0.3, position=pd) +
    geom_point(aes(y=total_binds2, x=months_deployedFO, colour = locexpdepth), 
               data = gridpairbindpiecemainSS2_gc, 
               position = position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data
    geom_line(aes(y=rate, x=months_deployedFO, 
                  group = locexpdepth,
                  colour = locexpdepth), data = freqmod2TB_preds, 
              position=pd) +
    scale_colour_manual(legend_titleT, values=colsT) +
    scale_fill_manual(legend_titleT, values = colsT) +
    scale_shape_manual(legend_titleT, values = shapesT) +
    theme_classic() + 
    theme(axis.title.x = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"),
          legend.title = element_text(face = "bold")) +   # Bold the legend title
    ylab("Number binds per rubble pair") + 
    xlab("Months since deployment") +
    ylim(0,10) +
    theme(panel.grid.major.y = element_line(color = "black", linetype = "dotted", size = 0.3)))









# Supp. 8 Bind density across sumcontact -----


# tried with interactions with habitat and time and not significant

freqmodSC2 <- glmmTMB(total_binds2 ~ sumcontact + 
                        (1|site_unique/grid_id.y), family=poisson, 
                      ziformula=~1, data = gridpairbindpiecemainSS2_gc)

car::Anova(freqmodSC2) 

summary(freqmodSC2)

freqmodSC2_grid <- with(gridpairbindpiecemainSS2_gc, 
                        list(sumcontact = seq(min(sumcontact,na.rm=TRUE),
                                              max(sumcontact,na.rm = TRUE),len=100)))

freqmodSC2_preds <- emmeans(freqmodSC2, ~ sumcontact, 
                            at = freqmodSC2_grid, type = "response") %>% as.data.frame()

head(freqmodSC2_preds) 

freqmod2TB_C3_preds2$months_deployedFO <- as.ordered(freqmod2TB_C3_preds2$months_deployedFO)

(freqmodSC2_plot <- ggplot() +
    geom_point(aes(y=total_binds2, x=sumcontact), 
               data = gridpairbindpiecemainSS2_gc, 
               position = position_jitter(width=0.2, height=0), 
               size =0.5, alpha= 0.1) + # raw data
    geom_ribbon(aes(y=rate, x=sumcontact,
                    ymin=asymp.LCL,ymax = asymp.UCL), 
                data = freqmodSC2_preds, 
                size =0.5, alpha= 0.1) + # raw data
    geom_line(aes(y=rate, x=sumcontact), data = freqmodSC2_preds) +
    #scale_colour_manual(legend_titleT, values=colsT) +
    #scale_fill_manual(legend_titleT, values = colsT) +
    #scale_shape_manual(legend_titleT, values = shapesT) +
    theme_classic() + 
    theme(axis.title.x = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"),
          legend.title = element_text(face = "bold")) +   # Bold the legend title
    ylab("Number binds per rubble pair") + 
    xlab("Sum contact length per pair") +
    #ylim(0,10) +
    theme(panel.grid.major.y = element_line(color = "black", linetype = "dotted", size = 0.3)))







# Supp. 9 Bind length across time and space -----


hist(gridpairbindpiecemainSS2_gc$bind_length_cum) # ziGamma should be good

bindLmod <- glmmTMB(bind_length_cum ~ months_deployedFO * locexpdepth + 
                      (1|site_unique/grid_id.y), family=ziGamma(link = "log"), 
                    ziformula = ~1,
                    data = gridpairbindpiecemainSS2_gc) 

car::Anova(bindLmod) # sig. interaction

plot(simulateResiduals(fittedModel = bindLmod)) # decent

bindLmodcomps <- emmeans(bindLmod, pairwise ~ locexpdepth | months_deployedFO, type = "response")

# Add this table to Supplementary Material

save_tables_to_docx(bindLmodcomps, "bindLmodcomps1.docx", "bindLmodcomps2.docx")

bindLmod_grid <- with(gridpairbindpiecemainSS2_gc, 
                      list(locexpdepth = levels(locexpdepth), 
                           months_deployedFO = levels(months_deployedFO)))

bindLmod_preds <- emmeans(bindLmod, ~ locexpdepth | months_deployedFO, 
                          at = bindLmod_grid, type = "response") %>% as.data.frame()

head(bindLmod_preds) 

# Plot it:

(bindLmod_plot <- ggplot() +
    geom_pointrange(aes(y=response, x=months_deployedFO, 
                        colour = locexpdepth,
                        fill = locexpdepth,
                        shape = locexpdepth,
                        ymin = asymp.LCL, 
                        ymax = asymp.UCL),
                    data = bindLmod_preds,
                    position=pd, size =0.8) + # predicted average data points +
    geom_errorbar(aes(y=response, x=months_deployedFO, 
                      ymin=asymp.LCL, ymax=asymp.UCL, 
                      colour = locexpdepth), 
                  data = bindLmod_preds, width=0.3, position=pd) +
    geom_point(aes(y=bind_length_cum, x=months_deployedFO, colour = locexpdepth), 
               data = gridpairbindpiecemainSS2_gc, 
               position = position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data
    geom_line(aes(y=response, x=months_deployedFO, 
                  group = locexpdepth,
                  colour = locexpdepth), data = bindLmod_preds, 
              position=pd) +
    scale_colour_manual(legend_titleT, values=colsT) +
    scale_fill_manual(legend_titleT, values = colsT) +
    scale_shape_manual(legend_titleT, values = shapesT) +
    theme_classic() + 
    theme(axis.title.x = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"),
          legend.title = element_text(face = "bold")) +   # Bold the legend title
    ylab("Bind length (cm) per rubble pair") + 
    xlab("Months since deployment") +
    ylim(0,10) +
    theme(panel.grid.major.y = element_line(color = "black", linetype = "dotted", size = 0.3)))






# Supp. 10 Bind length across sumcontact -----


bindLsum <- glmmTMB(bind_length_cum ~ sumcontact + 
                      (1|site_unique/grid_id.y), family=ziGamma(link = "log"), 
                    ziformula = ~1,
                    data = gridpairbindpiecemainSS2_gc) 

car::Anova(bindLsum)
summary(bindLsum)
inverse(0.093158)

bindLsum_grid <- with(gridpairbindpiecemainSS2_gc, 
                      list(sumcontact = seq(min(sumcontact,na.rm=TRUE),
                                            max(sumcontact,na.rm = TRUE),len=100)))

bindLsum_preds <- emmeans(bindLsum, ~ sumcontact, 
                          at = bindLsum_grid, type = "response") %>% as.data.frame()

head(bindLsum_preds) 

(bindLsum_plot <- ggplot() +
    geom_point(aes(y=bind_length_cum, x=sumcontact), 
               data = gridpairbindpiecemainSS2_gc, 
               position = position_jitter(width=0.2, height=0), 
               size =0.5, alpha= 0.1) + # raw data
    geom_ribbon(aes(y=response, x=sumcontact,
                    ymin=asymp.LCL,ymax = asymp.UCL), 
                data = bindLsum_preds, 
                size =0.5, alpha= 0.1) + # raw data
    geom_line(aes(y=response, x=sumcontact), data = bindLsum_preds) +
    #scale_colour_manual(legend_titleT, values=colsT) +
    #scale_fill_manual(legend_titleT, values = colsT) +
    #scale_shape_manual(legend_titleT, values = shapesT) +
    theme_classic() + 
    theme(axis.title.x = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"),
          legend.title = element_text(face = "bold")) +   # Bold the legend title
    ylab("Cumulative bind length (cm) per rubble pair") + 
    xlab("Sum contact length per pair") +
    #ylim(0,10) +
    theme(panel.grid.major.y = element_line(color = "black", linetype = "dotted", size = 0.3)))







# Supp. 11 Bind area across time and space -----


gridpairbindpiecemainSS2_gc$pairbindareacm_LT_org <- as.numeric(gridpairbindpiecemainSS2_gc$pairbindareacm_LT_org)
hist(na.omit(gridpairbindpiecemainSS2_gc$pairbindareacm_LT_org)) # very right skewed
range(gridpairbindpiecemainSS2_gc$pairbindareacm_LT_org, na.rm = TRUE) # 0.0001 to 9.4 cm2

areamod <- glmmTMB(pairbindareacm_LT_org ~ months_deployedFO * locexpdepth + 
                     (1|site_unique/grid_id.y), family=Gamma(link = "log"), 
                   data = gridpairbindpiecemainSS2_gc) 

car::Anova(areamod)

plot(simulateResiduals(fittedModel = areamod))

areamodcomps <- emmeans(areamod, pairwise ~ locexpdepth | months_deployedFO, type = "response")

# Add this table to Supplementary Material

save_tables_to_docx(areamodcomps, "areamodcomps1.docx", "areamodcomps2.docx")

areamod_grid <- with(gridpairbindpiecemainSS2, 
                     list(locexpdepth = levels(locexpdepth), 
                          months_deployedFO = levels(months_deployedFO)))

areamod_preds <- emmeans(areamod, ~ locexpdepth | months_deployedFO, 
                         at = areamod_grid, type = "response") %>% as.data.frame()

head(areamod_preds) 

# Plot it:

(areamod_plot <- ggplot() +
    geom_pointrange(aes(y=response, x=months_deployedFO, 
                        colour = locexpdepth,
                        fill = locexpdepth,
                        shape = locexpdepth,
                        ymin = asymp.LCL, 
                        ymax = asymp.UCL),
                    data = areamod_preds,
                    position=pd, size =0.8) + # predicted average data points +
    geom_errorbar(aes(y=response, x=months_deployedFO, 
                      ymin=asymp.LCL, ymax=asymp.UCL, 
                      colour = locexpdepth), 
                  data = areamod_preds, width=0.3, position=pd) +
    geom_point(aes(y=pairbindareacm_LT_org, x=months_deployedFO, colour = locexpdepth), 
               data = gridpairbindpiecemainSS2_gc, 
               position = position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data
    geom_line(aes(y=response, x=months_deployedFO, 
                  group = locexpdepth,
                  colour = locexpdepth), data = areamod_preds, 
              position=pd) +
    scale_colour_manual(legend_titleT, values=colsT) +
    scale_fill_manual(legend_titleT, values = colsT) +
    scale_shape_manual(legend_titleT, values = shapesT) +
    theme_classic() + 
    theme(axis.title.x = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"),
          legend.title = element_text(face = "bold")) +   # Bold the legend title
    ylab("Bind area per rubble pair (cm2)") + 
    xlab("Months since deployment") +
    ylim(0,2.5) +
    theme(panel.grid.major.y = element_line(color = "black", linetype = "dotted", size = 0.3)))







# Supp. 12. Bind area across sumcontact -----


bindAsum <- glmmTMB(pairbindareacm_LT_org ~ sumcontact + 
                      (1|site_unique/grid_id.y), family=Gamma(link = "log"), 
                    data = gridpairbindpiecemainSS2_gc) 

car::Anova(bindAsum)
summary(bindAsum)

bindAsum_grid <- with(gridpairbindpiecemainSS2_gc, 
                      list(sumcontact = seq(min(sumcontact,na.rm=TRUE),
                                            max(sumcontact,na.rm = TRUE),len=100)))

bindAsum_preds <- emmeans(bindAsum, ~ sumcontact, 
                          at = bindAsum_grid, type = "response") %>% as.data.frame()

head(bindAsum_preds) 

(bindAsum_plot <- ggplot() +
    geom_point(aes(y=pairbindareacm_LT_org, x=sumcontact), 
               data = gridpairbindpiecemainSS2_gc, 
               position = position_jitter(width=0.2, height=0), 
               size =0.5, alpha= 0.1) + # raw data
    geom_ribbon(aes(y=response, x=sumcontact,
                    ymin=asymp.LCL,ymax = asymp.UCL), 
                data = bindAsum_preds, 
                size =0.5, alpha= 0.1) + # raw data
    geom_line(aes(y=response, x=sumcontact), data = bindAsum_preds) +
    #scale_colour_manual(legend_titleT, values=colsT) +
    #scale_fill_manual(legend_titleT, values = colsT) +
    #scale_shape_manual(legend_titleT, values = shapesT) +
    theme_classic() + 
    theme(axis.title.x = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"),
          legend.title = element_text(face = "bold")) +   # Bold the legend title
    ylab("Cumulative bind area (cm2) per rubble pair") + 
    xlab("Sum contact length per pair") +
    #ylim(0,10) +
    theme(panel.grid.major.y = element_line(color = "black", linetype = "dotted", size = 0.3)))









# Supp. 12 Bind area across sumcontact -----

bindAsum <- glmmTMB(pairbindareacm_LT_org ~ sumcontact + 
                      (1|site_unique/grid_id.y), family=Gamma(link = "log"), 
                    data = gridpairbindpiecemainSS2_gc) 

car::Anova(bindAsum)
summary(bindAsum)

bindAsum_grid <- with(gridpairbindpiecemainSS2_gc, 
                      list(sumcontact = seq(min(sumcontact,na.rm=TRUE),
                                            max(sumcontact,na.rm = TRUE),len=100)))

bindAsum_preds <- emmeans(bindAsum, ~ sumcontact, 
                          at = bindAsum_grid, type = "response") %>% as.data.frame()

head(bindAsum_preds) 

(bindAsum_plot <- ggplot() +
    geom_point(aes(y=pairbindareacm_LT_org, x=sumcontact), 
               data = gridpairbindpiecemainSS2_gc, 
               position = position_jitter(width=0.2, height=0), 
               size =0.5, alpha= 0.1) + # raw data
    geom_ribbon(aes(y=response, x=sumcontact,
                    ymin=asymp.LCL,ymax = asymp.UCL), 
                data = bindAsum_preds, 
                size =0.5, alpha= 0.1) + # raw data
    geom_line(aes(y=response, x=sumcontact), data = bindAsum_preds) +
    #scale_colour_manual(legend_titleT, values=colsT) +
    #scale_fill_manual(legend_titleT, values = colsT) +
    #scale_shape_manual(legend_titleT, values = shapesT) +
    theme_classic() + 
    theme(axis.title.x = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"),
          legend.title = element_text(face = "bold")) +   # Bold the legend title
    ylab("Cumulative bind area (cm2) per rubble pair") + 
    xlab("Sum contact length per pair") +
    #ylim(0,10) +
    theme(panel.grid.major.y = element_line(color = "black", linetype = "dotted", size = 0.3)))












# Supp. 13	The effect of bind density on break force -----


# Already tested for interaction between locexpdepth and total_binds on break_force, and it was not significant.
# Also tested for months and not sig interaction between number binds and months

TBBFcorr <- glmmTMB(break_force ~ 0 + total_binds2 +
                      (1|site_unique/grid_id.y), family=ziGamma(link = "log"), 
                    ziformula=~1, data = gridpairbindpiecemainSS2_gc)

summary(TBBFcorr)
car::Anova(TBBFcorr)
plot(simulateResiduals(fittedModel = TBBFcorr)) # decent

tidy_results <- broom::tidy(TBBFcorr)

write.csv(tidy_results, "TBBFcorr_summary.csv")

TBBFcorr_grid <- with(gridpairbindpiecemainSS2_gc, 
                      list(total_binds2 = seq(0, 12, len =100)))

TBBFcorr_preds <- emmeans(TBBFcorr, ~ total_binds2, at = TBBFcorr_grid, type = "response") %>% as.data.frame()

head(TBBFcorr_preds)

TBBFcorr_plot <- ggplot() +
  geom_point(aes(y=break_force, x=total_binds2), 
             data = gridpairbindpiecemainSS2_gc, 
             position = position_jitter(width=0.2, height=0), colour = "grey", size =0.5, alpha= 0.3) + # raw data
  geom_ribbon(aes(y=response, x=total_binds2, 
                  ymin=asymp.LCL, ymax=asymp.UCL), 
              data = TBBFcorr_preds, fill = "black", alpha = 0.2) +
  geom_line(aes(y=response, x=total_binds2), size = 1, colour = "black", data = TBBFcorr_preds) +
  theme_classic() + 
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +   # Bold the legend title
  ylab("Force to break binds (N)") + 
  xlab("Number binds per rubble pair") +
  theme(panel.grid.major.y = element_line(color = "black", linetype = "dotted", size = 0.3)) +
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12)) + 
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 20))  # Zoom in on both x and y axes

TBBFcorr_plot








# Supp. 14	The effect of bind length on break force -----


# Already tested for interaction between locexpdepth and total_binds on break_force, and it was not significant.
# what about including months? Also tested for and not sig interaction between number binds and months

# Dataframe:

# Want to put 0 if no binds (to match to total)

BLBFcorr <- glmmTMB(break_force ~ 0 + bind_length_cum +
                      (1|locexpdepth/site_unique/grid_id.y), family=ziGamma(link = "log"), 
                    ziformula = ~1, data = gridpairbindpiecemainSS2_gc)

summary(BLBFcorr)
car::Anova(BLBFcorr)

plot(simulateResiduals(fittedModel = BLBFcorr))

tidy_results <- broom::tidy(BLBFcorr)

write.csv(tidy_results, "BLBFcorr_summary.csv")

range(gridpairbindpiecemainSS2_gc$bind_length_cum,na.rm = TRUE)

BLBFcorr_grid <- with(gridpairbindpiecemainSS2_gc, 
                      list(bind_length_cum = seq(0, 14.9, len =100)))

BLBFcorr_preds <- emmeans(BLBFcorr, ~ bind_length_cum, 
                          at = BLBFcorr_grid, type = "response") %>% 
  as.data.frame()

head(BLBFcorr_preds)

BLBFcorr_plot <- ggplot() +
  geom_point(aes(y=break_force, x=bind_length_cum), 
             data = gridpairbindpiecemainSS2_gc, 
             position = position_jitter(width=0.2, height=0), colour = "grey", size =0.5, alpha= 0.3) + # raw data
  geom_ribbon(aes(y=response, x=bind_length_cum,
                  ymin=asymp.LCL, ymax=asymp.UCL), 
              data = BLBFcorr_preds, fill = "black", alpha = 0.2) +
  geom_line(aes(y=response, x=bind_length_cum), color = "black", size = 1, data = BLBFcorr_preds) +
  theme_classic() + 
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +   # Bold the legend title
  ylab("Force to break binds (N)") + 
  xlab("Cumulative bind length (cm) per rubble pair") +
  theme(panel.grid.major.y = element_line(color = "black", linetype = "dotted", size = 0.3))+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12)) + 
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 20))  # Zoom in on both x and y axes

BLBFcorr_plot









# Supp. 15	The effect of bind area on break force -----


BABFcorr <- glmmTMB(break_force ~ 0 + pairbindareacm_LT_org +
                      (1|locexpdepth/site_unique/grid_id.y), family=ziGamma(link = "log"), 
                    ziformula = ~1, data = gridpairbindpiecemainSS2_gc)

summary(BABFcorr)
car::Anova(BABFcorr)

plot(simulateResiduals(fittedModel = BABFcorr)) # tweedie no better...

tidy_results <- broom::tidy(BABFcorr)

write.csv(tidy_results, "BABFcorr_summary.csv")

range(gridpairbindpiecemainSS2_gc$pairbindareacm_LT_org,na.rm = TRUE)

BABFcorr_grid <- with(gridpairbindpiecemainSS2_gc, 
                      list(pairbindareacm_LT_org = seq(0, 10, len =100)))

BABFcorr_preds <- emmeans(BABFcorr, ~ pairbindareacm_LT_org, 
                          at = BABFcorr_grid, type = "response") %>% 
  as.data.frame()

head(BABFcorr_preds)

BABFcorr_plot <- ggplot() +
  geom_point(aes(y=break_force, x=pairbindareacm_LT_org), 
             data = gridpairbindpiecemainSS2_gc, 
             position = position_jitter(width=0.2, height=0), 
             colour = "grey", size =0.5, alpha= 0.3) + # raw data
  geom_ribbon(aes(y=response, x=pairbindareacm_LT_org,
                  ymin=asymp.LCL, ymax=asymp.UCL), 
              data = BABFcorr_preds, fill = "black", alpha = 0.2) +
  geom_line(aes(y=response,
                x=pairbindareacm_LT_org), 
            color = "black", size = 1, data = BABFcorr_preds) +
  theme_classic() + 
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +   # Bold the legend title
  ylab("Force to break binds (N)") + 
  xlab("Cumulative bind area (cm2) per rubble pair") +
  theme(panel.grid.major.y = element_line(color = "black", linetype = "dotted", size = 0.3))+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12)) + 
  coord_cartesian(xlim = c(0, 4.5), ylim = c(0, 20))  # Zoom in on both x and y axes

BABFcorr_plot










# Supp. 16	Break force across time and space -----


freqmod1 <- glmmTMB(break_force ~ months_deployedFO * locexpdepth + 
                      (1|site_unique/grid_id.y), family=ziGamma(link = "log"), 
                    ziformula=~1, data = gridpairbindpiecemainSS2_gc) # don't need to use gridpairbindpiecemainSS3 for frequentist because it can take out the NAs by itself


summary(freqmod1)
car::Anova(freqmod1) # There is an interaction

freqmod1comps <- emmeans(freqmod1, pairwise ~  locexpdepth| months_deployedFO, type = "response")

freqmod2comps <- emmeans(freqmod1, pairwise ~  months_deployedFO | locexpdepth, type = "response")

save_tables_to_docx(freqmod2TBcomps, "freqmod1comps1.docx", "freqmod1comps2.docx")

save_tables_to_docx(freqmod2comps, "freqmod2comps1.docx", "freqmod2comps2.docx")


freq1_grid <- with(gridpairbindpiecemainSS2, list(locexpdepth = levels(locexpdepth),
                                                  months_deployedFO = levels(months_deployedFO)))

freqmod1_preds <- emmeans(freqmod1, ~ locexpdepth | months_deployedFO, 
                          at = freq1_grid, type = "response") %>% as.data.frame()

legend_titleT <- "Habitat"
colsT <- c("In Shelt Shal" = "#6CA6CD",
           "In Exp Shal" = "#CD8C95", 
           "Off Shelt RF" = "#AEE2F5",
           "Off Shelt Shal" = "#6CA6CD",
           "Off Shelt Deep" = "#104E8B",
           "Off Exp Shal" = "#CD8C95",    
           "Off Exp Deep" = "#A52A2A")

shapesT <- c("In Shelt Shal" = 8,
             "In Exp Shal" = 8, 
             "Off Shelt RF" = 17,
             "Off Shelt Shal" = 17,
             "Off Shelt Deep" = 17,
             "Off Exp Shal" = 17,    
             "Off Exp Deep" = 17)

linesT <- c("In Shelt Shal" = 2,
            "In Exp Shal" = 2, 
            "Off Shelt RF" = 1,
            "Off Shelt Shal" = 1,
            "Off Shelt Deep" = 1,
            "Off Exp Shal" = 1,    
            "Off Exp Deep" = 1)

pd <- position_dodge(0.55)

# Months deployed

freqmod1_plot <- ggplot() +
  geom_pointrange(aes(y=response, x=months_deployedFO, 
                      colour = locexpdepth,
                      fill = locexpdepth,
                      shape = locexpdepth,
                      ymin = asymp.LCL, 
                      ymax = asymp.UCL),
                  data = freqmod1_preds,
                  position=pd, size =0.8) + # predicted average data points +
  geom_errorbar(aes(y=response, x=months_deployedFO, 
                    ymin=asymp.LCL, ymax=asymp.UCL, 
                    colour = locexpdepth), 
                data = freqmod1_preds, width=0.3, position=pd) +
  geom_point(aes(y=break_force, x=months_deployedFO, colour = locexpdepth), 
             data = gridpairbindpiecemainSS2, 
             position = position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data
  geom_line(aes(y=response, x=months_deployedFO, 
                group = locexpdepth,
                colour = locexpdepth,
                shape = locexpdepth), data = freqmod1_preds, 
            position=pd) +
  geom_hline(yintercept = seq(0, 20, by = 4), linetype = "dotted", color = "black", size = 0.3) +
  scale_colour_manual(legend_titleT, values=colsT) +
  scale_fill_manual(legend_titleT, values = colsT) +
  scale_shape_manual(legend_titleT, values = shapesT) +
  theme_classic() + 
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +   # Bold the legend title
  theme(legend.position="right") +
  ylab("Force to break binds (N)") + 
  xlab("Months since deployment") + scale_y_continuous(breaks = c(0,4,8,12,16,20), limits = c(0,20))

freqmod1_plot







# Supp. 17	Tau (force over area) across time and space ----

Taumod1 <- glmmTMB(Tbind_Ncm_LT_org ~ months_deployedFO * locexpdepth + 
                     (1|site_unique/grid_id.y), family=ziGamma(link = "log"), 
                   ziformula=~1, data = gridpairbindpiecemainSS2_gc) # don't need to use gridpairbindpiecemainSS3 for frequentist because it can take out the NAs by itself

summary(Taumod1)

# So for the ordered factor, you end up with a linear component, a quadratic component (which is the curvature), and a cubic component
# months_deployedFO.L is the linear component: it is positive (0.9), and significant (P = 2.24e-16), it means that the line is going up (as it is in our raw data)
# months_deployedFO.Q is the quadratic, which is the curvature and it is negative and significant.
# months_deployedFO.C is the cubic and it is positive but not significant.

car::Anova(Taumod1) # Interaction

plot(simulateResiduals(fittedModel = Taumod1)) # decent

emmeans(Taumod1, pairwise ~ locexpdepth | months_deployedFO, type = "response")

Taumod1comps <- emmeans(Taumod1, pairwise ~  locexpdepth| months_deployedFO, type = "response")

Taumod12comps <- emmeans(Taumod1, pairwise ~  months_deployedFO | locexpdepth, type = "response")

save_tables_to_docx(Taumod1comps, "Taumod1comps1.docx", "Taumod1comps12.docx")

save_tables_to_docx(Taumod12comps, "Taumod12comps1.docx", "Taumod12comps2.docx")


# Plot it:

Taumod11_grid <- with(gridpairbindpiecemainSS2_gc, list(locexpdepth = levels(locexpdepth),
                                                        months_deployedFO = levels(months_deployedFO)))

Taumod11_preds <- emmeans(Taumod1, ~ locexpdepth | months_deployedFO, 
                          at = Taumod11_grid, type = "response") %>% as.data.frame()

ggplot() +
  geom_pointrange(aes(y=response, x=months_deployedFO, 
                      colour = locexpdepth,
                      fill = locexpdepth,
                      shape = locexpdepth,
                      ymin = asymp.LCL, 
                      ymax = asymp.UCL),
                  data = Taumod11_preds,
                  position=pd, size =0.8) + # predicted average data points +
  geom_errorbar(aes(y=response, x=months_deployedFO, 
                    ymin=asymp.LCL, ymax=asymp.UCL, 
                    colour = locexpdepth), 
                data = Taumod11_preds, width=0.3, position=pd) +
  geom_point(aes(y=Tbind_Ncm_LT_org, x=months_deployedFO, colour = locexpdepth), 
             data = gridpairbindpiecemainSS2, 
             position = position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data
  geom_line(aes(y=response, x=months_deployedFO, 
                group = locexpdepth,
                colour = locexpdepth,
                shape = locexpdepth), data = Taumod11_preds, 
            position=pd) +
  #geom_hline(yintercept = seq(0, 20, by = 4), linetype = "dotted", color = "black", size = 0.3) +
  scale_colour_manual(legend_titleT, values=colsT) +
  scale_fill_manual(legend_titleT, values = colsT) +
  scale_shape_manual(legend_titleT, values = shapesT) +
  theme_classic() + 
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +   # Bold the legend title
  ylab("Tau - Force (N) / bind area (cm2)") + 
  xlab("Months since deployment") +
  scale_y_continuous(breaks = c(0,20,40,60,80), limits = c(0,80))






# Supp. 18	Break velocity across time and space ----


velmod <- glmmTMB(break_velocity_BINDAREA_LT_org ~ months_deployedFO *
                    locexpdepth +
                    (1|site_unique/grid_id.y), family=ziGamma(link = "log"), 
                  ziformula=~1, data = gridpairbindpiecemainSS2_gc) 

car::Anova(velmod) # sig interaction

plot(simulateResiduals(fittedModel = velmod)) # good

velmodcomps <- emmeans(velmod, pairwise ~ locexpdepth | months_deployedFO, type = "response")

velmodcomps2 <- emmeans(velmod, pairwise ~  months_deployedFO | locexpdepth, type = "response")

# Add this table to Supplementary Material

save_tables_to_docx(velmodcomps, "velmodcomps1.docx", "velmodcomps2.docx")
save_tables_to_docx(velmodcomps2, "velmodcomps2_1.docx", "velmodcomps2_2.docx")

velmod_grid <- with(gridpairbindpiecemainSS2_gc, 
                    list(locexpdepth = levels(locexpdepth), 
                         months_deployedFO = levels(months_deployedFO)))

velmod_preds <- emmeans(velmod, ~ locexpdepth | months_deployedFO, 
                        at = velmod_grid, type = "response") %>% as.data.frame()

head(velmod_preds) 

# Plot it:

(velmod_plot <- ggplot() +
    geom_pointrange(aes(y=response, x=months_deployedFO, 
                        colour = locexpdepth,
                        fill = locexpdepth,
                        shape = locexpdepth,
                        ymin = asymp.LCL, 
                        ymax = asymp.UCL),
                    data = velmod_preds,
                    position=pd, size =0.8) + # predicted average data points +
    geom_errorbar(aes(y=response, x=months_deployedFO, 
                      ymin=asymp.LCL, ymax=asymp.UCL, 
                      colour = locexpdepth), 
                  data = velmod_preds, width=0.3, position=pd) +
    geom_point(aes(y=break_velocity_BINDAREA_LT_org, x=months_deployedFO, colour = locexpdepth), 
               data = gridpairbindpiecemainSS2_gc, 
               position = position_jitter(width=0.2, height=0), size =0.5, alpha= 0.1) + # raw data
    geom_line(aes(y=response, x=months_deployedFO, 
                  group = locexpdepth,
                  colour = locexpdepth), data = velmod_preds, 
              position=pd) +
    scale_colour_manual(legend_titleT, values=colsT) +
    scale_fill_manual(legend_titleT, values = colsT) +
    scale_shape_manual(legend_titleT, values = shapesT) +
    theme_classic() + 
    theme(axis.title.x = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"),
          legend.title = element_text(face = "bold")) +   # Bold the legend title
    ylab("Velocity to break apart bound rubble pairs (m/s)") + 
    xlab("Months since deployment") +
    ylim(0,4.5) +
    theme(panel.grid.major.y = element_line(color = "black", linetype = "dotted", size = 0.3)))






# Supp. 19	Varying pair configuration - Bind density, bind length, break force and break velocity -----

# Number of binds -----

crossdat <- gridpairbindpiececrossS2_gc

hist(crossdat$total_binds2) # again right-skewed but data is discrete

freqmodcrsN1 <- glmmTMB(total_binds2 ~ pair_arrangement * locexpdepth + 
                          (1|grid_id.y),  
                        ziformula=~1, family=poisson(), data = crossdat) 

car::Anova(freqmodcrsN1) # interaction

plot(simulateResiduals(fittedModel = freqmodcrsN1)) # looks good

freqmodcrsN1comps <- emmeans(freqmodcrsN1, pairwise ~ pair_arrangement | locexpdepth, type = "response")

save_tables_to_docx(freqmodcrsN1comps, "freqmodcrsN1comps1.docx", "freqmodcrsN1comps2.docx")

emmeans(freqmodcrsN1, pairwise ~  locexpdepth | pair_arrangement, type = "response")


# Bind length -----

# Have to first add a column for cumulative bind length to cross data:

gridbindercross$binder_DNA_broadest_cat <- as.factor(gridbindercross$binder_DNA_broadest_cat)
levels(gridbindercross$binder_DNA_broadest_cat)

gridbindercross <- gridbindercross %>% mutate(binder_DNA_broadest_cat = factor(binder_DNA_broadest_cat,
                                                                               levels = c("Turf algae", "Lobophora", "Other macroalgae",
                                                                                          "Peyssonnelia",
                                                                                          "Sponge", "Colonial ascidian", "Solitary ascidian", 
                                                                                          "Colonial bryozoan","Serpulid worm", "Vermetid snail", 
                                                                                          "Bivalve","Other")))


gridbindercross_cum <- gridbindercross %>%
  group_by(unique_pair_id, grid_id, 
           location, site_unique, depth, months_deployed, 
           exposure, layer) %>% 
  summarise(bind_length_cum = sum(bind_length, na.rm = TRUE)) %>%
  mutate(bind_length_cum = ifelse(is.na(bind_length_cum), 0, bind_length_cum)) %>%
  ungroup()

crossdat <- crossdat %>%
  left_join(gridbindercross_cum %>% select(unique_pair_id, bind_length_cum), by = "unique_pair_id")

str(crossdat) # this has bind area, based on length x thickness which is used to calculate velocity,


# Model:
hist(crossdat$bind_length_cum, na.rm=TRUE)
range(crossdat$break_force, na.rm=TRUE)

lengthmodcrsN1 <- glmmTMB(bind_length_cum ~ pair_arrangement * locexpdepth + 
                            (1|grid_id.y),  
                          ziformula=~1, family=ziGamma(link = "log"), data = crossdat) 

car::Anova(lengthmodcrsN1) # no interaction

lengthmodcrsN1 <- glmmTMB(bind_length_cum ~ pair_arrangement + locexpdepth + 
                            (1|grid_id.y),  
                          ziformula=~1, family=ziGamma(link = "log"), data = crossdat) 

car::Anova(lengthmodcrsN1) # effect of pair type but not depth

lengthmodcrsN1comps <- emmeans(lengthmodcrsN1, pairwise ~ pair_arrangement, type = "response")

save_tables_to_docx(lengthmodcrsN1comps, "lengthmodcrsN1comps1.docx", "lengthmodcrsN1comps2.docx")

velmodcrs_grid = with(crossdat, list(pair_arrangement = levels(pair_arrangement),
                                     locexpdepth = levels(locexpdepth)))

lengthmodcrsN1preds <- emmeans(lengthmodcrsN1,  ~ pair_arrangement,
                               at = velmodcrs_grid,
                               type = "response") %>% as.data.frame()

head(lengthmodcrsN1preds)

pd <- position_dodge(0.5) 

plotBL <- ggplot() +
  geom_point(aes(y = bind_length_cum, x = pair_arrangement, colour = site_unique), 
             data = crossdat, 
             position = position_jitter(width = 0.2, height = 0), size = 2, alpha = 0.5) + # raw data
  geom_pointrange(aes(y = response, x = pair_arrangement,
                      ymin = asymp.LCL, 
                      ymax = asymp.UCL),
                  data = lengthmodcrsN1preds,
                  position = pd, size = 0.8) + # predicted average data points
  geom_errorbar(aes(y = response, x = pair_arrangement,
                    ymin = asymp.LCL, ymax = asymp.UCL),  # Add ymin and ymax here
                data = lengthmodcrsN1preds, width = 0.3, position = pd) +
  theme_classic() + 
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +   # Bold the legend title
  ylab("Cumulative bind length per pair (N)") +
  xlab("Pair configuration") +
  # ylim(0,10) +
  theme(panel.grid.major.y = element_line(color = "black", 
                                          linetype = "dotted", size = 0.3)) +
  theme(legend.position = "right")
plotBL

# Break force ----

forcemodcrs1 <- glmmTMB(break_force ~ pair_arrangement * locexpdepth + 
                          (1|grid_id.y), family=ziGamma(link = "log"), 
                        ziformula=~1, data = crossdat) #

plot(simulateResiduals(fittedModel = forcemodcrs1)) #looks good!

car::Anova(forcemodcrs1) # interaction not significant but leave in as have for others

forcemodcrs1comps <- emmeans(forcemodcrs1, pairwise ~ pair_arrangement | locexpdepth,type = "response") 

save_tables_to_docx(forcemodcrs1comps, "forcemodcrs1comps1.docx", 
                    "forcemodcrs1comps2.docx")



# Break velocity -----


velmodcrs1 <- glmmTMB(break_velocity_BINDAREA_LT_org ~ pair_arrangement * locexpdepth + 
                        (1|grid_id.y), family=ziGamma(link = "log"), 
                      ziformula=~1, data = crossdat)


plot(simulateResiduals(fittedModel = velmodcrs1)) #looks good!

car::Anova(velmodcrs1) # interaction significant

velmodcrs1comps <- emmeans(velmodcrs1, pairwise ~ pair_arrangement | locexpdepth,
                           type = "response")

save_tables_to_docx(velmodcrs1comps, "velmodcrs1comps1.docx", 
                    "velmodcrs1comps2.docx")

velmodcrs_grid = with(crossdat, list(pair_arrangement = levels(pair_arrangement),
                                     locexpdepth = levels(locexpdepth)))

velmodcrs1preds <- emmeans(velmodcrs1,  ~ pair_arrangement | locexpdepth,
                           at = velmodcrs_grid,
                           type = "response") %>% as.data.frame()

head(velmodcrs1preds)

pd <- position_dodge(0.5) 

colsT <- c("Off Shelt Shal" = "#6CA6CD",
           "Off Shelt Deep" = "#104E8B")
shapesT <- c("Off Shelt Shal" = 17,
             "Off Shelt Deep" = 17)

legend_titleT <- "Habitat"

plotBV <- ggplot() +
  geom_pointrange(aes(y = response, x = pair_arrangement,
                      colour = locexpdepth,
                      shape = locexpdepth,
                      ymin = asymp.LCL, 
                      ymax = asymp.UCL),
                  data = velmodcrs1preds,
                  position = pd, size = 0.8) + # predicted average data points
  geom_errorbar(aes(y = response, x = pair_arrangement,
                    colour = locexpdepth,
                    ymin = asymp.LCL, ymax = asymp.UCL),  # Add ymin and ymax here
                data = velmodcrs1preds, width = 0.3, position = pd) +
  geom_point(aes(y = break_velocity_BINDAREA_LT_org, x = pair_arrangement, 
                 colour = locexpdepth), 
             data = crossdat, 
             position = position_jitter(width = 0.2, height = 0), size = 0.5, alpha = 0.1) + # raw data
  scale_colour_manual(name = legend_titleT, values = colsT) +
  scale_shape_manual(name = legend_titleT, values = shapesT) +
  theme_classic() + 
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +   # Bold the legend title
  ylab("Velocity to break binding per pair (N)") +
  xlab("Pair configuration") +
  # ylim(0,10) +
  theme(panel.grid.major.y = element_line(color = "black", 
                                          linetype = "dotted", size = 0.3)) +
  theme(legend.position = "none")
plotBV

# Figure 3C ----
# Break force & velocity plot combined:

grid <- with(crossdat, list(pair_arrangement = levels(pair_arrangement),
                            locexpdepth = levels(locexpdepth)))

predsBF <- emmeans(forcemodcrs1, ~ pair_arrangement | locexpdepth, type = "response") %>% as.data.frame()
predsTB <- emmeans(freqmodcrsN1, ~ pair_arrangement | locexpdepth, type = "response") %>% as.data.frame()

scale_factor <- 1
velocity_offset <- 2.2  # Adjust this value as needed

combined_plotPAIRS <- ggplot() +
  # Break force plot (primary y-axis)
  geom_pointrange(aes(y = response, x = pair_arrangement,
                      colour = locexpdepth,
                      shape = locexpdepth,
                      ymin = asymp.LCL, ymax = asymp.UCL),
                  data = predsBF,
                  position = pd, size = 0.8) +
  geom_errorbar(aes(y = response, x = pair_arrangement,
                    colour = locexpdepth,
                    ymin = asymp.LCL, ymax = asymp.UCL),
                data = predsBF, width = 0.3, position = pd) +
  geom_point(aes(y = break_force, x = pair_arrangement, 
                 colour = locexpdepth), 
             data = crossdat, 
             position = position_jitter(width = 0.2, height = 0), size = 0.5, alpha = 0.1) +
  
  # Break velocity plot (secondary y-axis), adjusted for offset
  geom_pointrange(aes(y = response / scale_factor + velocity_offset, x = pair_arrangement, 
                      colour = locexpdepth,
                      shape = locexpdepth,
                      ymin = asymp.LCL / scale_factor + velocity_offset, 
                      ymax = asymp.UCL / scale_factor + velocity_offset),
                  data = velmodcrs1preds,
                  position = pd, size = 0.8) +
  geom_errorbar(aes(y = response / scale_factor + velocity_offset, x = pair_arrangement,
                    colour = locexpdepth,
                    ymin = asymp.LCL / scale_factor + velocity_offset, 
                    ymax = asymp.UCL / scale_factor + velocity_offset),
                data = velmodcrs1preds, width = 0.3, linetype = "dashed", position = pd) +
  geom_point(aes(y = break_velocity_BINDAREA_LT_org + velocity_offset, x = pair_arrangement, 
                 colour = locexpdepth), 
             data = crossdat, 
             position = position_jitter(width = 0.2, height = 0), size = 0.5, alpha = 0.1) +
  
  # Scales
  scale_colour_manual(name = legend_titleT, values = colsT) +
  scale_shape_manual(name = legend_titleT, values = shapesT) +
  
  # Primary y-axis (for break force)
  ylab("Force to break binding per pair (N)") +
  
  # Secondary y-axis (for break velocity, scaled)
  scale_y_continuous(
    limits = c(0, 4),  # Adjust limits as necessary to accommodate the offset
    sec.axis = sec_axis(~ . * scale_factor - velocity_offset, name = "Velocity to break binding per pair (m/s)")) +
  
  # x-axis
  xlab("Pair configuration") +
  
  # Theme
  theme_classic() + 
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.title.y.right = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.grid.major.y = element_line(color = "black", 
                                          linetype = "dotted", size = 0.3),
        legend.position = "none")

# Display the combined plot
combined_plotPAIRS










# Supp. 20 Binding community across time and space -----

gridbinderB <- gridbinder %>%
  filter(!is.na(binder_DNA_broadest_cat)) # removing those where there is no binding.
str(gridbinderB)

view(gridbinderB)

# Formatting

gridbinderB$layer <- factor(gridbinderB$layer ,
                            levels = c('Top', 'Bottom'),
                            labels = c('Surface grids', 'Buried grids'))
gridbinderB$location <- factor(gridbinderB$location,
                               levels = c('Keppels', 'Heron'),
                               labels = c('Keppels (inshore)', 'Heron (offshore)'))

gridbinderB$exposuredepth <- paste(gridbinderB$exposure, gridbinderB$depth, sep = "_")
gridbinderB$exposuredepth <- as.factor(gridbinderB$exposuredepth)

gridbinderB$exposuredepth <- factor(gridbinderB$exposuredepth, levels = c("Low_Reef flat",
                                                                          "Low_Shallow slope",
                                                                          "Low_Deep slope",
                                                                          "High_Shallow slope", 
                                                                          "High_Deep slope"),
                                    labels = c("Sheltered Reef Flat", "Sheltered Shallow",
                                               "Sheltered Deep", "Exposed Shallow", "Exposed Deep"))

gridbinderB$locexpdepth <- paste(gridbinderB$location, gridbinderB$exposuredepth, sep = "_")
gridbinderB$locexpdepth <- as.factor(gridbinderB$locexpdepth)
levels(gridbinderB$locexpdepth)

gridbinderB$locexpdepth <- factor(gridbinderB$locexpdepth, levels = c("Keppels (inshore)_Sheltered Shallow",
                                                                      "Keppels (inshore)_Exposed Shallow",
                                                                      "Heron (offshore)_Sheltered Reef Flat",
                                                                      "Heron (offshore)_Sheltered Shallow" , 
                                                                      "Heron (offshore)_Sheltered Deep",
                                                                      "Heron (offshore)_Exposed Shallow",
                                                                      "Heron (offshore)_Exposed Deep"),
                                  labels = c("In Shelt Shal", "In Exp Shal", "Off Shelt RF", "Off Shelt Shal",
                                             "Off Shelt Deep", "Off Exp Shal", "Off Exp Deep"))

levels(gridbinderB$exposuredepth)
levels(gridbinderB$location)
levels(gridbinderB$layer)
levels(gridbinderB$locexpdepth)

gridbinderB$binder_DNA_broadest_cat <- factor(gridbinderB$binder_DNA_broadest_cat, 
                                              levels = c("Turf algae","Lobophora" ,"Other macroalgae" ,"CCA",
                                                         "Peyssonnelia","Sponge" , "Colonial ascidian",
                                                         "Solitary ascidian","Colonial bryozoan",
                                                         "Serpulid worm","Vermetid snail","Bivalve","Hard coral","Other"))

# Make months an ordered factor again:

gridbinderB$months_deployedFO <- ordered(gridbinderB$months_deployed)
levels(gridbinderB$months_deployedFO)


# Let's look at surface grids only now:

gridbinderBS <- gridbinderB %>% filter(layer != "Buried grids") %>% droplevels()
str(gridbinderBS) # now there are 1661 rows
levels(gridbinderBS$binder_DNA_broadest_cat)

# Bind density

typecols3 <- c("Turf algae" = "#B4EEB4",
               "Lobophora" = "#B8860B", 
               "Other macroalgae" = "#79B07A",
               "CCA" =  "#EEB4B4",
               "Peyssonnelia" = "#F5BC49",
               "Sponge" = "#FFF68F",     
               "Colonial ascidian" = "#C799C6",
               "Solitary ascidian" = "#B452CD",
               "Colonial bryozoan" = "lightblue2",
               "Serpulid worm" = "#FFDEAD",
               "Vermetid snail" = "red",
               "Bivalve" = "darkgrey",
               "Hard coral" = "peachpuff3",
               "Other" = "#919DBD")


gridbinderBS_summ1 <- gridbinderBS %>%
  group_by(locexpdepth, months_deployedFO, binder_DNA_broadest_cat) %>%
  summarise(total_binds_per_type = n(), na.rm = TRUE)

(gridbinderBSplot <- ggplot(gridbinderBS_summ1, 
                            aes(x = locexpdepth, y = total_binds_per_type, fill = binder_DNA_broadest_cat)) +
    geom_bar(stat = "identity", position = "fill", colour = "black") +  # Stacks the bars to show proportions
    facet_wrap(~months_deployedFO) +  # Facet by months_deployedFO
    scale_fill_manual(values = typecols3) +  # Set colors to your predefined palette
    labs(y = "Proportion bind density", fill = "Binder Type") +  # Update y-axis label and legend title
    theme_classic() +
    theme(text = element_text(face="bold", size = 12)) + 
    theme(axis.text = element_text(colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position = "right") +
    theme(axis.text.x.bottom = element_text(size=12)) +
    theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=0.8)))

(gridbinderBSplot <- ggplot(gridbinderBS_summ1, 
                            aes(x = locexpdepth, y = total_binds_per_type, fill = binder_DNA_broadest_cat)) +
    geom_bar(stat = "identity", position = "stack", colour = "black") +  # Stacks the bars to show proportions
    facet_wrap(~months_deployedFO) +  # Facet by months_deployedFO
    scale_fill_manual(values = typecols3) +  # Set colors to your predefined palette
    labs(y = "Proportion bind density", fill = "Binder Type") +  # Update y-axis label and legend title
    theme_classic() +
    theme(text = element_text(face="bold", size = 12)) + 
    theme(axis.text = element_text(colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position = "right") +
    theme(axis.text.x.bottom = element_text(size=12)) +
    theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=0.8)))

# Bind length

gridbinderBS_summ <- gridbinderBS %>%
  group_by(locexpdepth, months_deployedFO, binder_DNA_broadest_cat) %>%
  summarise(total_bind_length = sum(bind_length, na.rm = TRUE))

(gridbinderBLplot <- ggplot(gridbinderBS_summ, 
                            aes(x = locexpdepth, y = total_bind_length, fill = binder_DNA_broadest_cat)) +
    geom_bar(stat = "identity", position = "fill", colour = "black") +  # Stacks the bars to show proportions
    facet_wrap(~months_deployedFO) +  # Facet by months_deployedFO
    scale_fill_manual(values = typecols3) +  # Set colors to your predefined palette
    labs(y = "Proportion bind length", fill = "Binder Type") +  # Update y-axis label and legend title
    theme_classic() +
    theme(text = element_text(face="bold", size = 12)) + 
    theme(axis.text = element_text(colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position = "right") +
    theme(axis.text.x.bottom = element_text(size=12)) +
    theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=0.8)))

(gridbinderBLplot <- ggplot(gridbinderBS_summ, 
                            aes(x = locexpdepth, y = total_bind_length, fill = binder_DNA_broadest_cat)) +
    geom_bar(stat = "identity", position = "fill", colour = "black") +  # Stacks the bars to show proportions
    facet_wrap(~months_deployedFO) +  # Facet by months_deployedFO
    scale_fill_manual(values = typecols3) +  # Set colors to your predefined palette
    labs(y = "Proportion bind length", fill = "Binder Type") +  # Update y-axis label and legend title
    theme_classic() +
    theme(text = element_text(face="bold", size = 12)) + 
    theme(axis.text = element_text(colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position = "right") +
    theme(axis.text.x.bottom = element_text(size=12)) +
    theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=0.8)))

# Bind area

gridbinderBS_summ_BA <- gridbinderBS %>%
  group_by(locexpdepth, months_deployedFO, binder_DNA_broadest_cat) %>%
  summarise(total_bind_area = sum(bind_length_x_thickness, na.rm = TRUE))

(gridbinderBAplot <- ggplot(gridbinderBS_summ_BA, 
                            aes(x = locexpdepth, y = total_bind_area, fill = binder_DNA_broadest_cat)) +
    geom_bar(stat = "identity", position = "fill", colour = "black") +  # Stacks the bars to show proportions
    facet_wrap(~months_deployedFO) +  # Facet by months_deployedFO
    scale_fill_manual(values = typecols3) +  # Set colors to your predefined palette
    labs(y = "Proportion bind area", fill = "Binder Type") +  # Update y-axis label and legend title
    theme_classic() +
    theme(text = element_text(face="bold", size = 12)) + 
    theme(axis.text = element_text(colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position = "right") +
    theme(axis.text.x.bottom = element_text(size=12)) +
    theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=0.8)))

(gridbinderBAplot2 <- ggplot(gridbinderBS_summ_BA, 
                             aes(x = locexpdepth, y = total_bind_area, fill = binder_DNA_broadest_cat)) +
    geom_bar(stat = "identity", position = "stack", colour = "black") +  # Stacks the bars to show proportions
    facet_wrap(~months_deployedFO) +  # Facet by months_deployedFO
    scale_fill_manual(values = typecols3) +  # Set colors to your predefined palette
    labs(y = "Total bind area", fill = "Binder Type") +  # Update y-axis label and legend title
    theme_classic() +
    theme(text = element_text(face="bold", size = 12)) + 
    theme(axis.text = element_text(colour = "black")) + 
    theme(panel.background =element_rect(colour = "black", size=1)) + 
    theme(axis.ticks.length=unit(.2,"cm")) +
    theme(legend.position = "right") +
    theme(axis.text.x.bottom = element_text(size=12)) +
    theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=0.8)))







# Supp. 21 Top contributors to binding strength – bind density -----

# Make dataframe:

gridbinderB <- gridbinder %>%
  filter(!is.na(binder_DNA_broadest_cat)) # removing those where there is no binding.
str(gridbinderB)

# Remove microhabitat, we want the binders per pair to relate to strength
gridbinderB_count <- gridbinderB %>%
  group_by(unique_pair_id, grid_id, binder_DNA_broadest_cat,
           location, site_unique, depth, months_deployed, 
           exposure, layer) %>% # removing microhabitat for my purposes (comparing to strength categories)
  summarise(Count = n()) %>%
  ungroup()

# Convert from long to wide, including '0' for missing combinations
gridbinderB_wide <- gridbinderB_count %>%
  pivot_wider(names_from = binder_DNA_broadest_cat,
              values_from = Count,
              values_fill = 0) 

# Also add a total_binds column

gridbinderB_wideTot <- left_join(gridbinderB_wide, gridpairbindpiecemainS2 %>% select(unique_pair_id, break_force), by = "unique_pair_id") # There are 1957 observations, does this sound right? pairmain has 3075 rows.
sum(pairmain$bound == "Yes", na.rm = TRUE) # 1958 rows where binding was Yes. And, one of these must have an NA break force, thus 1957 rows now.

gridbinderB_wideTot <- gridbinderB_wideTot %>%
  rowwise() %>%
  mutate(total_binds = sum(c_across("Turf algae":"Vermetid snail")))

# Remove the break_force = NA

gridbinderB_wideTotNoNA <- gridbinderB_wideTot %>%
  filter(!is.na(break_force)) 

names(gridbinderB_wideTotNoNA) <- gsub(" ", "_", names(gridbinderB_wideTotNoNA))
str(gridbinderB_wideTotNoNA)

# Now let's format, create the locexpdepth category etc.

# Formatting

gridbinderB_wideTotNoNA$layer <- factor(gridbinderB_wideTotNoNA$layer ,
                                        levels = c('Top', 'Bottom'),
                                        labels = c('Surface grids', 'Buried grids'))
gridbinderB_wideTotNoNA$location <- factor(gridbinderB_wideTotNoNA$location,
                                           levels = c('Keppels', 'Heron'),
                                           labels = c('Keppels (inshore)', 'Heron (offshore)'))

gridbinderB_wideTotNoNA$exposuredepth <- paste(gridbinderB_wideTotNoNA$exposure, gridbinderB_wideTotNoNA$depth, sep = "_")

gridbinderB_wideTotNoNA$exposuredepth <- as.factor(gridbinderB_wideTotNoNA$exposuredepth)

levels(gridbinderB_wideTotNoNA$exposuredepth)

gridbinderB_wideTotNoNA$exposuredepth <- factor(gridbinderB_wideTotNoNA$exposuredepth, levels = c("Low_Reef flat",
                                                                                                  "Low_Shallow slope",
                                                                                                  "Low_Deep slope",
                                                                                                  "High_Shallow slope", 
                                                                                                  "High_Deep slope"),
                                                labels = c("Sheltered Reef Flat", "Sheltered Shallow",
                                                           "Sheltered Deep", "Exposed Shallow", "Exposed Deep"))

gridbinderB_wideTotNoNA$locexpdepth <- paste(gridbinderB_wideTotNoNA$location, gridbinderB_wideTotNoNA$exposuredepth, sep = "_")
gridbinderB_wideTotNoNA$locexpdepth <- as.factor(gridbinderB_wideTotNoNA$locexpdepth)
levels(gridbinderB_wideTotNoNA$locexpdepth)

gridbinderB_wideTotNoNA$locexpdepth <- factor(gridbinderB_wideTotNoNA$locexpdepth, levels = c("Keppels (inshore)_Sheltered Shallow",
                                                                                              "Keppels (inshore)_Exposed Shallow",
                                                                                              "Heron (offshore)_Sheltered Reef Flat",
                                                                                              "Heron (offshore)_Sheltered Shallow" , 
                                                                                              "Heron (offshore)_Sheltered Deep",
                                                                                              "Heron (offshore)_Exposed Shallow",
                                                                                              "Heron (offshore)_Exposed Deep"),
                                              labels = c("In Shelt Shal", "In Exp Shal", "Off Shelt RF", "Off Shelt Shal",
                                                         "Off Shelt Deep", "Off Exp Shal", "Off Exp Deep"))

levels(gridbinderB_wideTotNoNA$exposuredepth)
levels(gridbinderB_wideTotNoNA$location)
levels(gridbinderB_wideTotNoNA$layer)
levels(gridbinderB_wideTotNoNA$locexpdepth)

# Make months an ordered factor again:

gridbinderB_wideTotNoNA$months_deployedFO <- ordered(gridbinderB_wideTotNoNA$months_deployed)
levels(gridbinderB_wideTotNoNA$months_deployedFO)


# Primer Dataframe ----

# Surface grids only:

levels(gridbinderB_wideTotNoNA$layer)

gridbinderB_wideTotNoNAS <- gridbinderB_wideTotNoNA %>% 
  filter(layer == "Surface grids")
dim(gridbinderB_wideTotNoNA) # 1910 with surface and buried grids
dim(gridbinderB_wideTotNoNAS) # now there are 1661 rows


# Bind density model

modB2 <- glmmTMB(break_force ~ Turf_algae +  Lobophora + Other_macroalgae +
                   CCA + Peyssonnelia + Sponge + Solitary_ascidian + Colonial_ascidian + 
                   Colonial_bryozoan + Serpulid_worm + Vermetid_snail +
                   Bivalve + Hard_coral + Other + 
                   (1|site_unique/grid_id), family = Gamma(link = "log"), 
                 data = gridbinderB_wideTotNoNAS)

AICc(modB, modB2) # Gamma is better

plot(simulateResiduals(fittedModel = modB2)) # looks ok

car::Anova(modB2)

summary(modB2)

tidy_results <- broom::tidy(modB2)

write.csv(tidy_results, "modB2summary.csv")

# Plotting the coefficients as well:

coefs2 <- tidy(modB2, effects = "fixed") %>% as.data.frame()

# Convert back to be on the response scale:
coefs2 <- coefs2 %>%
  mutate(
    response_scale_estimate = exp(estimate),
    response_scale_std_error = exp(estimate) * std.error,  # Adjust standard errors
    response_scale_lower = exp(estimate - 1.96 * std.error),
    response_scale_upper = exp(estimate + 1.96 * std.error))

# Reorder the term column in coefs2
coefs2 <- coefs2 %>%
  mutate(term = factor(term, levels = rev(c("(Intercept)","Turf_algae", "Lobophora", "Other_macroalgae",
                                            "CCA", "Peyssonnelia",
                                            "Sponge", "Colonial_ascidian", "Solitary_ascidian", 
                                            "Colonial_bryozoan","Serpulid_worm", "Vermetid_snail", "Bivalve", "Hard_coral",         
                                            "Other"))))

# Plot:

(coefsplot2 <- ggplot(coefs2, aes(x = term, y = response_scale_estimate, 
                                  ymin = response_scale_lower, ymax = response_scale_upper)) +
    geom_pointrange() +
    coord_flip() +
    xlab("Binder Category") +
    ylab("Coefficient Estimate") +
    theme_minimal())

# Get the variables in the same order as this plot and then put side by side.

predictor_vars <- all.vars(formula(modB2))[-1] 
predictor_vars <- setdiff(predictor_vars, c("grid_id", "site_unique"))

results_list <- list()
# For loop to iterate over predictor variables and perform the model
for (pred in predictor_vars) {
  # Create a grid for the current predictor variable
  grid <- with(gridbinderB_wideTotNoNAS, 
               list(seq(min(gridbinderB_wideTotNoNAS[[pred]]), max(gridbinderB_wideTotNoNAS[[pred]]), length.out = 50))) # is this enough preds.
  names(grid) <- pred
  
  # Perform emmeans for the current predictor variable
  emmeans_result <- emmeans(modB2, reformulate(pred), at = grid, type = "response") %>%
    as.data.frame()
  
  # Store the result in the list with the name of the predictor
  results_list[[pred]] <- emmeans_result
}

results_list <- as.data.frame(results_list)

# Will have to change results_list to long form.
# But first need to have the same predictor for each one (e.g., number binds from 0 to 15)
# Change column names so I can reformat:

names(results_list)

results_list2 <- results_list %>%
  rename(
    Turf_algae.x = Turf_algae.Turf_algae,
    Other_macroalgae.x = Other_macroalgae.Other_macroalgae,
    Lobophora.x = Lobophora.Lobophora,
    CCA.x = CCA.CCA,
    Peyssonnelia.x = Peyssonnelia.Peyssonnelia,
    Sponge.x = Sponge.Sponge,
    Colonial_ascidian.x = Colonial_ascidian.Colonial_ascidian,
    Solitary_ascidian.x = Solitary_ascidian.Solitary_ascidian,
    Colonial_bryozoan.x = Colonial_bryozoan.Colonial_bryozoan,
    Serpulid_worm.x = Serpulid_worm.Serpulid_worm,
    Vermetid_snail.x = Vermetid_snail.Vermetid_snail,
    Bivalve.x = Bivalve.Bivalve,
    Hard_coral.x = Hard_coral.Hard_coral,
    Other.x = Other.Other)

str(results_list2)

results_list3 <- results_list2 %>%
  rename_with(~ gsub("\\.LCL$", "_LCL", .x), ends_with(".LCL"))

results_list3 <- results_list3 %>%
  rename_with(~ gsub("\\.UCL$", "_UCL", .x), ends_with(".UCL"))

results_list3_long <- results_list3 %>%
  pivot_longer(
    cols = everything(),  # Select all columns
    names_to = c("Binder", "Metric"),  # Create 'Binder' and 'Metric' columns
    names_pattern = "(.*)\\.(x|response|SE|df|asymp_LCL|asymp_UCL)",  # Regex to capture Binder and Metric
    values_to = "Value"  # Place values in 'Value' column
  ) %>% as.data.frame()

str(results_list3_long)
results_list3_long <- results_list3_long %>% mutate_if(is.character, as.factor)
head(results_list3_long)

# Create a replicate identifier within each Binder-Matrix group
results_list3_long <- results_list3_long %>%
  group_by(Binder, Metric) %>%
  mutate(Replicate = row_number()) %>%
  ungroup()

# Pivot the dataframe to wide format, keeping replicates as separate rows
results_list3_wide <- results_list3_long %>%
  pivot_wider(
    names_from = Metric,  # Use Matrix levels as column names
    values_from = Value   # Fill columns with Value data
  )

# View the transformed dataframe
head(results_list3_wide)

results_list3_wide <-  results_list3_wide %>% as.data.frame()
head(results_list3_wide)

levels(results_list3_wide$Binder)

results_list3_wide <- results_list3_wide %>% mutate(Binder = factor(Binder,
                                                                    levels = c("Turf_algae", "Lobophora", "Other_macroalgae",
                                                                               "CCA", "Peyssonnelia",
                                                                               "Sponge", "Colonial_ascidian", "Solitary_ascidian", 
                                                                               "Colonial_bryozoan","Serpulid_worm", "Vermetid_snail", "Bivalve", "Hard_coral",         
                                                                               "Other")))


typecols2 <- c("Turf_algae" = "#B4EEB4",
               "Lobophora" = "#B8860B",
               "Other_macroalgae" = "#79B07A",
               "CCA" =  "#EEB4B4",
               "Peyssonnelia" = "#F5BC49",
               "Sponge" = "#FFF68F",     
               "Colonial_ascidian" = "#C799C6",
               "Solitary_ascidian" = "#B452CD",
               "Colonial_bryozoan" = "lightblue2",
               "Serpulid_worm" = "#FFDEAD",
               "Vermetid_snail" = "red",
               "Bivalve" = "darkgrey",
               "Hard_coral" = "peachpuff3",
               "Other" = "#919DBD")

(all.plot2 <- ggplot(results_list3_wide, aes(x = x, y = response, fill = Binder)) + 
    scale_y_continuous(limits = c(0,310), breaks = seq(0, 310, by = 10)) +  # Customize y-axis breaks
    scale_x_continuous(limits = c(0,7.2), breaks = seq(0, 7, by = 1)) +  # Customize y-axis breaks
    geom_hline(yintercept = seq(0, 310, by = 10), linetype = "dotted", color = "gray") + # Add dotted horizontal lines)
    geom_ribbon(aes(ymin = asymp_LCL, ymax = asymp_UCL), alpha = 0.15) +
    geom_line(aes(x = x, y = response, colour = Binder), size = 0.9) + 
    scale_colour_manual(name = "Binder", values = typecols2) +
    scale_fill_manual(name = "Binder", values = typecols2) +
    labs(x = "Number of binds", y = "Force to break bound rubble pair (N)") +
    theme_classic() +
    coord_cartesian(xlim = c(1, 7), ylim = c(0, 100)))  # Zoom in on both x and y axes






# Supp. 22	Top contributors to binding strength – bind length -----

gridbinderB_cum <- gridbinderB %>%
  group_by(unique_pair_id, grid_id, binder_DNA_broadest_cat,
           location, site_unique, depth, months_deployed, 
           exposure, layer) %>% # removing microhabitat for my purposes (comparing to strength categories)
  summarise(pair_bind_length_per_binder = sum(bind_length),
            pair_bind_area_per_binder_BW = sum(bind_length_x_width),
            pair_bind_area_per_binder_BT = sum(bind_length_x_thickness)) %>%
  ungroup()

# Bind length wide form
gridbinderB_cum_wideBL <- gridbinderB_cum %>%
  pivot_wider(id_cols = c("unique_pair_id", "grid_id", 
                          "location", "site_unique", "depth", 
                          "months_deployed", "exposure", "layer"),
              names_from = binder_DNA_broadest_cat,
              values_from = pair_bind_length_per_binder,
              values_fill = 0)

# Bind area length*thickness wide form
gridbinderB_cum_wideBA <- gridbinderB_cum %>%
  pivot_wider(id_cols = c("unique_pair_id", "grid_id", 
                          "location", "site_unique", "depth", 
                          "months_deployed", "exposure", "layer"),
              names_from = binder_DNA_broadest_cat,
              values_from = pair_bind_area_per_binder_BT,
              values_fill = 0)

# Join with gridpairbindpiecemainS2 which doesn't have rolled tests, but still has bottom grids, and remove later

gridbinderB_cum_wideBL_FORCE <- left_join(gridbinderB_cum_wideBL, gridpairbindpiecemainS2 %>% select(unique_pair_id, break_force), by = "unique_pair_id") # There are 1957 observations, does this sound right? pairmain has 3075 rows.

gridbinderB_cum_wideBA_FORCE <- left_join(gridbinderB_cum_wideBA, gridpairbindpiecemainS2 %>% select(unique_pair_id, break_force), by = "unique_pair_id") # There are 1957 observations, does this sound right? pairmain has 3075 rows.

# Remove the break_force = NA

gridbinderB_cum_wideBL_FORCENoNA <- gridbinderB_cum_wideBL_FORCE %>%
  filter(!is.na(break_force)) 

names(gridbinderB_cum_wideBL_FORCENoNA) <- gsub(" ", "_", names(gridbinderB_cum_wideBL_FORCENoNA))

gridbinderB_cum_wideBA_FORCENoNA <- gridbinderB_cum_wideBA_FORCE %>%
  filter(!is.na(break_force)) 

names(gridbinderB_cum_wideBA_FORCENoNA) <- gsub(" ", "_", names(gridbinderB_cum_wideBA_FORCENoNA))

# Now let's format, create the locexpdepth category etc.

# Formatting

gridbinderB_cum_wideBL_FORCENoNA$layer <- factor(gridbinderB_cum_wideBL_FORCENoNA$layer ,
                                                 levels = c('Top', 'Bottom'),
                                                 labels = c('Surface grids', 'Buried grids'))
gridbinderB_cum_wideBL_FORCENoNA$location <- factor(gridbinderB_cum_wideBL_FORCENoNA$location,
                                                    levels = c('Keppels', 'Heron'),
                                                    labels = c('Keppels (inshore)', 'Heron (offshore)'))

gridbinderB_cum_wideBL_FORCENoNA$exposuredepth <- paste(gridbinderB_cum_wideBL_FORCENoNA$exposure, gridbinderB_cum_wideBL_FORCENoNA$depth, sep = "_")

gridbinderB_cum_wideBL_FORCENoNA$exposuredepth <- as.factor(gridbinderB_cum_wideBL_FORCENoNA$exposuredepth)

levels(gridbinderB_cum_wideBL_FORCENoNA$exposuredepth)

gridbinderB_cum_wideBL_FORCENoNA$exposuredepth <- factor(gridbinderB_cum_wideBL_FORCENoNA$exposuredepth, levels = c("Low_Reef flat",
                                                                                                                    "Low_Shallow slope",
                                                                                                                    "Low_Deep slope",
                                                                                                                    "High_Shallow slope", 
                                                                                                                    "High_Deep slope"),
                                                         labels = c("Sheltered Reef Flat", "Sheltered Shallow",
                                                                    "Sheltered Deep", "Exposed Shallow", "Exposed Deep"))

gridbinderB_cum_wideBL_FORCENoNA$locexpdepth <- paste(gridbinderB_cum_wideBL_FORCENoNA$location, gridbinderB_cum_wideBL_FORCENoNA$exposuredepth, sep = "_")
gridbinderB_cum_wideBL_FORCENoNA$locexpdepth <- as.factor(gridbinderB_cum_wideBL_FORCENoNA$locexpdepth)
levels(gridbinderB_cum_wideBL_FORCENoNA$locexpdepth)

gridbinderB_cum_wideBL_FORCENoNA$locexpdepth <- factor(gridbinderB_cum_wideBL_FORCENoNA$locexpdepth, levels = c("Keppels (inshore)_Sheltered Shallow",
                                                                                                                "Keppels (inshore)_Exposed Shallow",
                                                                                                                "Heron (offshore)_Sheltered Reef Flat",
                                                                                                                "Heron (offshore)_Sheltered Shallow" , 
                                                                                                                "Heron (offshore)_Sheltered Deep",
                                                                                                                "Heron (offshore)_Exposed Shallow",
                                                                                                                "Heron (offshore)_Exposed Deep"),
                                                       labels = c("In Shelt Shal", "In Exp Shal", "Off Shelt RF", "Off Shelt Shal",
                                                                  "Off Shelt Deep", "Off Exp Shal", "Off Exp Deep"))

levels(gridbinderB_cum_wideBL_FORCENoNA$exposuredepth)
levels(gridbinderB_cum_wideBL_FORCENoNA$location)
levels(gridbinderB_cum_wideBL_FORCENoNA$layer)
levels(gridbinderB_cum_wideBL_FORCENoNA$locexpdepth)

gridbinderB_cum_wideBA_FORCENoNA$layer <- factor(gridbinderB_cum_wideBA_FORCENoNA$layer ,
                                                 levels = c('Top', 'Bottom'),
                                                 labels = c('Surface grids', 'Buried grids'))
gridbinderB_cum_wideBA_FORCENoNA$location <- factor(gridbinderB_cum_wideBA_FORCENoNA$location,
                                                    levels = c('Keppels', 'Heron'),
                                                    labels = c('Keppels (inshore)', 'Heron (offshore)'))

gridbinderB_cum_wideBA_FORCENoNA$exposuredepth <- paste(gridbinderB_cum_wideBA_FORCENoNA$exposure, gridbinderB_cum_wideBA_FORCENoNA$depth, sep = "_")

gridbinderB_cum_wideBA_FORCENoNA$exposuredepth <- as.factor(gridbinderB_cum_wideBA_FORCENoNA$exposuredepth)

levels(gridbinderB_cum_wideBA_FORCENoNA$exposuredepth)

gridbinderB_cum_wideBA_FORCENoNA$exposuredepth <- factor(gridbinderB_cum_wideBA_FORCENoNA$exposuredepth, levels = c("Low_Reef flat",
                                                                                                                    "Low_Shallow slope",
                                                                                                                    "Low_Deep slope",
                                                                                                                    "High_Shallow slope", 
                                                                                                                    "High_Deep slope"),
                                                         labels = c("Sheltered Reef Flat", "Sheltered Shallow",
                                                                    "Sheltered Deep", "Exposed Shallow", "Exposed Deep"))

gridbinderB_cum_wideBA_FORCENoNA$locexpdepth <- paste(gridbinderB_cum_wideBA_FORCENoNA$location, gridbinderB_cum_wideBA_FORCENoNA$exposuredepth, sep = "_")
gridbinderB_cum_wideBA_FORCENoNA$locexpdepth <- as.factor(gridbinderB_cum_wideBA_FORCENoNA$locexpdepth)
levels(gridbinderB_cum_wideBA_FORCENoNA$locexpdepth)

gridbinderB_cum_wideBA_FORCENoNA$locexpdepth <- factor(gridbinderB_cum_wideBA_FORCENoNA$locexpdepth, levels = c("Keppels (inshore)_Sheltered Shallow",
                                                                                                                "Keppels (inshore)_Exposed Shallow",
                                                                                                                "Heron (offshore)_Sheltered Reef Flat",
                                                                                                                "Heron (offshore)_Sheltered Shallow" , 
                                                                                                                "Heron (offshore)_Sheltered Deep",
                                                                                                                "Heron (offshore)_Exposed Shallow",
                                                                                                                "Heron (offshore)_Exposed Deep"),
                                                       labels = c("In Shelt Shal", "In Exp Shal", "Off Shelt RF", "Off Shelt Shal",
                                                                  "Off Shelt Deep", "Off Exp Shal", "Off Exp Deep"))

levels(gridbinderB_cum_wideBA_FORCENoNA$exposuredepth)
levels(gridbinderB_cum_wideBA_FORCENoNA$location)
levels(gridbinderB_cum_wideBA_FORCENoNA$layer)
levels(gridbinderB_cum_wideBA_FORCENoNA$locexpdepth)

# Make months an ordered factor again:

gridbinderB_cum_wideBL_FORCENoNA$months_deployedFO <- ordered(gridbinderB_cum_wideBL_FORCENoNA$months_deployed)
levels(gridbinderB_cum_wideBL_FORCENoNA$months_deployedFO)

gridbinderB_cum_wideBA_FORCENoNA$months_deployedFO <- ordered(gridbinderB_cum_wideBA_FORCENoNA$months_deployed)
levels(gridbinderB_cum_wideBA_FORCENoNA$months_deployedFO)

# Let's look at surface grids only now:

gridbinderB_cum_wideBL_FORCENoNAS <- gridbinderB_cum_wideBL_FORCENoNA %>% filter(layer != "Buried grids") %>% droplevels()

gridbinderB_cum_wideBA_FORCENoNAS <- gridbinderB_cum_wideBA_FORCENoNA %>% filter(layer != "Buried grids") %>% droplevels()

dim(gridbinderB_cum_wideBL_FORCENoNAS) # now there are 1661 rows
dim(gridbinderB_cum_wideBA_FORCENoNAS) # now there are 1661 rows

# Which organisms are contributing the most to the break force over time?
# By bind length & then bind area

# Add rubble length in case it is included as a covariate in the Primer analyses

str(gridbinderB_cum_wideBA_FORCENoNAS)

gridbinderB_cum_wideBL_FORCENoNAS_L <- left_join(gridbinderB_cum_wideBL_FORCENoNAS, 
                                                 gridpairbindpiecemainS2 %>% select(unique_pair_id, avgrubblelength), 
                                                 by = "unique_pair_id") 


gridbinderB_cum_wideBA_FORCENoNAS_L <- left_join(gridbinderB_cum_wideBA_FORCENoNAS, 
                                                 gridpairbindpiecemainS2 %>% select(unique_pair_id, avgrubblelength), 
                                                 by = "unique_pair_id")


# Write to csv for Primer
write.csv(gridbinderB_cum_wideBL_FORCENoNAS_L, "gridbinderB_cum_wideBL_FORCENoNAS_L.csv")
write.csv(gridbinderB_cum_wideBA_FORCENoNAS_L, "gridbinderB_cum_wideBA_FORCENoNAS_L.csv")


# Bind length model:

hist(gridbinderB_cum_wideBL_FORCENoNAS$break_force)

modBL2 <- glmmTMB(break_force ~ Turf_algae +  Lobophora + Other_macroalgae +
                    CCA + Peyssonnelia + Sponge + Solitary_ascidian + Colonial_ascidian + 
                    Colonial_bryozoan + Serpulid_worm + Vermetid_snail +
                    Bivalve + Hard_coral + Other + 
                    (1|site_unique/grid_id), family = Gamma(link = "log"), 
                  data = gridbinderB_cum_wideBL_FORCENoNAS)

plot(simulateResiduals(fittedModel = modBL2)) # looks ok

car::Anova(modBL2)
plot(allEffects(modBL2))
summary(modBL2)

tidy_resultsBL <- broom::tidy(modBL2)

write.csv(tidy_resultsBL, "modBL2summary.csv")

# Figure 5b -----

coefsBL2 <- tidy(modBL2, effects = "fixed") %>% as.data.frame()

# Convert back to be on the response scale:
coefsBL2 <- coefsBL2 %>%
  mutate(
    response_scale_estimate = exp(estimate),
    response_scale_std_error = exp(estimate) * std.error,  # Adjust standard errors
    response_scale_lower = exp(estimate - 1.96 * std.error),
    response_scale_upper = exp(estimate + 1.96 * std.error))

# Reorder the term column in coefs2
coefsBL2 <- coefsBL2 %>%
  mutate(term = factor(term, levels = rev(c("(Intercept)","Turf_algae", "Lobophora", "Other_macroalgae",
                                            "CCA", "Peyssonnelia",
                                            "Sponge", "Colonial_ascidian", "Solitary_ascidian", 
                                            "Colonial_bryozoan","Serpulid_worm", "Vermetid_snail", "Bivalve", "Hard_coral",         
                                            "Other"))))

view(coefsBL2)
(coefsBL2plot <- ggplot(coefsBL2, aes(x = term, y = response_scale_estimate, 
                                      ymin = response_scale_lower, ymax = response_scale_upper)) +
    geom_pointrange() +
    coord_flip() +
    xlab("Binder Category") +
    ylab("Coefficient Estimate") +
    theme_minimal())


predictor_vars <- all.vars(formula(modBL2))[-1] 
predictor_vars <- setdiff(predictor_vars, c("grid_id", "site_unique"))

# Create a list to store the results
results_list <- list()
# For loop to iterate over predictor variables and perform the model
for (pred in predictor_vars) {
  # Create a grid for the current predictor variable
  grid <- with(gridbinderB_cum_wideBL_FORCENoNAS, 
               list(seq(min(gridbinderB_cum_wideBL_FORCENoNAS[[pred]]), max(gridbinderB_cum_wideBL_FORCENoNAS[[pred]]), length.out = 50))) # is this enough preds.
  names(grid) <- pred
  
  # Perform emmeans for the current predictor variable
  emmeans_result <- emmeans(modBL2, reformulate(pred), at = grid, type = "response") %>%
    as.data.frame()
  
  # Store the result in the list with the name of the predictor
  results_list[[pred]] <- emmeans_result
}

view(results_list)
results_list <- as.data.frame(results_list)

names(results_list)

results_list2 <- results_list %>%
  rename(
    Turf_algae.x = Turf_algae.Turf_algae,
    Other_macroalgae.x = Other_macroalgae.Other_macroalgae,
    Lobophora.x = Lobophora.Lobophora,
    CCA.x = CCA.CCA,
    Peyssonnelia.x = Peyssonnelia.Peyssonnelia,
    Sponge.x = Sponge.Sponge,
    Colonial_ascidian.x = Colonial_ascidian.Colonial_ascidian,
    Solitary_ascidian.x = Solitary_ascidian.Solitary_ascidian,
    Colonial_bryozoan.x = Colonial_bryozoan.Colonial_bryozoan,
    Serpulid_worm.x = Serpulid_worm.Serpulid_worm,
    Vermetid_snail.x = Vermetid_snail.Vermetid_snail,
    Bivalve.x = Bivalve.Bivalve,
    Hard_coral.x = Hard_coral.Hard_coral,
    Other.x = Other.Other)

str(results_list2)

results_list3 <- results_list2 %>%
  rename_with(~ gsub("\\.LCL$", "_LCL", .x), ends_with(".LCL"))

results_list3 <- results_list3 %>%
  rename_with(~ gsub("\\.UCL$", "_UCL", .x), ends_with(".UCL"))

results_list3_long <- results_list3 %>%
  pivot_longer(
    cols = everything(),  # Select all columns
    names_to = c("Binder", "Metric"),  # Create 'Binder' and 'Metric' columns
    names_pattern = "(.*)\\.(x|response|SE|df|asymp_LCL|asymp_UCL)",  # Regex to capture Binder and Metric
    values_to = "Value"  # Place values in 'Value' column
  ) %>% as.data.frame()

str(results_list3_long)
results_list3_long <- results_list3_long %>% mutate_if(is.character, as.factor)
head(results_list3_long)

# Create a replicate identifier within each Binder-Matrix group
results_list3_long <- results_list3_long %>%
  group_by(Binder, Metric) %>%
  mutate(Replicate = row_number()) %>%
  ungroup()

# Pivot the dataframe to wide format, keeping replicates as separate rows
results_list3_wide <- results_list3_long %>%
  pivot_wider(
    names_from = Metric,  # Use Matrix levels as column names
    values_from = Value   # Fill columns with Value data
  )

# View the transformed dataframe
head(results_list3_wide)

results_list3_wide <-  results_list3_wide %>% as.data.frame()
head(results_list3_wide)

levels(results_list3_wide$Binder)

results_list3_wide <- results_list3_wide %>% mutate(Binder = factor(Binder,
                                                                    levels = c("Turf_algae", "Lobophora", "Other_macroalgae",
                                                                               "CCA", "Peyssonnelia",
                                                                               "Sponge", "Colonial_ascidian", "Solitary_ascidian", 
                                                                               "Colonial_bryozoan","Serpulid_worm", "Vermetid_snail", "Bivalve", "Hard_coral",         
                                                                               "Other")))


# Figure 5a ----

typecols2 <- c("Turf_algae" = "#B4EEB4",
               "Lobophora" = "#B8860B",
               "Other_macroalgae" = "#79B07A",
               "CCA" =  "#EEB4B4",
               "Peyssonnelia" = "#F5BC49",
               "Sponge" = "#FFF68F",     
               "Colonial_ascidian" = "#C799C6",
               "Solitary_ascidian" = "#B452CD",
               "Colonial_bryozoan" = "lightblue2",
               "Serpulid_worm" = "#FFDEAD",
               "Vermetid_snail" = "red",
               "Bivalve" = "darkgrey",
               "Hard_coral" = "peachpuff3",
               "Other" = "#919DBD")

(all.plot2 <- ggplot(results_list3_wide, aes(x = x, y = response, fill = Binder)) + 
    scale_y_continuous(limits = c(0,310), breaks = seq(0, 310, by = 10)) +  # Customize y-axis breaks
    scale_x_continuous(limits = c(0,7.2), breaks = seq(0, 7, by = 1)) +  # Customize y-axis breaks
    geom_hline(yintercept = seq(0, 310, by = 10), linetype = "dotted", color = "gray") + # Add dotted horizontal lines)
    geom_ribbon(aes(ymin = asymp_LCL, ymax = asymp_UCL), alpha = 0.15) +
    geom_line(aes(x = x, y = response, colour = Binder), size = 0.9) + 
    scale_colour_manual(name = "Binder", values = typecols2) +
    scale_fill_manual(name = "Binder", values = typecols2) +
    labs(x = "Bind length (cm)", y = "Force to break binding per pair (N)") +
    theme_classic() +
    coord_cartesian(xlim = c(1, 7), ylim = c(0, 100)))  # Zoom in on both x and y axes









# Supp. 23	Top contributors to binding strength – bind area -----


modBA2 <- glmmTMB(break_force ~ Turf_algae +  Lobophora + Other_macroalgae +
                    CCA + Peyssonnelia + Sponge + Solitary_ascidian + Colonial_ascidian + 
                    Colonial_bryozoan + Serpulid_worm + Vermetid_snail +
                    Bivalve + Hard_coral + Other + 
                    (1|site_unique/grid_id), 
                  family = Gamma(link = "log"), 
                  data = gridbinderB_cum_wideBA_FORCENoNAS)

str(gridbinderB_cum_wideBA_FORCENoNAS)

plot(simulateResiduals(fittedModel = modBA2)) # looks ok

car::Anova(modBA2)
plot(allEffects(modBA2))
summary(modBA2)

tidy_resultsBA2 <- broom::tidy(modBA2)

write.csv(tidy_resultsBA2, "modBA2summary2.csv")

# Plotting the coefficients as well:

coefsBA2 <- tidy(modBA2, effects = "fixed") %>% as.data.frame()

# Convert back to be on the response scale:

coefsBA2 <- coefsBA2 %>%
  mutate(
    response_scale_estimate = exp(estimate),
    response_scale_std_error = exp(estimate) * std.error,  # Adjust standard errors
    response_scale_lower = exp(estimate - 1.96 * std.error),
    response_scale_upper = exp(estimate + 1.96 * std.error))

# Reorder the term column in coefs2
coefsBA2 <- coefsBA2 %>%
  mutate(term = factor(term, levels = rev(c("(Intercept)","Turf_algae", "Lobophora", "Other_macroalgae",
                                            "CCA", "Peyssonnelia",
                                            "Sponge", "Colonial_ascidian", "Solitary_ascidian", 
                                            "Colonial_bryozoan","Serpulid_worm", "Vermetid_snail", "Bivalve", "Hard_coral",         
                                            "Other"))))

# Plot:

(coefsBA2plot <- ggplot(coefsBA2, aes(x = term, y = response_scale_estimate, 
                                      ymin = response_scale_lower, ymax = response_scale_upper)) +
    geom_pointrange() +
    coord_flip() +
    xlab("Binder Category") +
    ylab("Coefficient Estimate") +
    theme_minimal()) # why is this coming out weird now? CCA is the highest 

# Bind length is easy to measure.
# the bind 'area' can change depending on the binder style, bindlength x thickness will
# work for sponges that bind in between pieces, but not for a CCA encrusting that bridges between but
# also encrusts a lot of the surface of each rubble piece, providing more strength that way
# i.e., bindlength x width (or half width) might correlate better in some cases, e.g., where half of an encrusting ascidian
# is binding across two rubble pairs, and peels off one of the piece (so the area of bind length x half width)

predictor_vars <- all.vars(formula(modBA2))[-1] 
predictor_vars <- setdiff(predictor_vars, c("grid_id", "site_unique"))

results_list <- list()
# For loop to iterate over predictor variables and perform the model
for (pred in predictor_vars) {
  # Create a grid for the current predictor variable
  grid <- with(gridbinderB_cum_wideBA_FORCENoNAS, 
               list(seq(min(gridbinderB_cum_wideBA_FORCENoNAS[[pred]]), max(gridbinderB_cum_wideBA_FORCENoNAS[[pred]]), length.out = 50))) # is this enough preds.
  names(grid) <- pred
  
  # Perform emmeans for the current predictor variable
  emmeans_result <- emmeans(modBA2, reformulate(pred), at = grid, type = "response") %>%
    as.data.frame()
  
  # Store the result in the list with the name of the predictor
  results_list[[pred]] <- emmeans_result
}

results_list <- as.data.frame(results_list)

results_list2 <- results_list %>%
  rename(
    Turf_algae.x = Turf_algae.Turf_algae,
    Other_macroalgae.x = Other_macroalgae.Other_macroalgae,
    Lobophora.x = Lobophora.Lobophora,
    CCA.x = CCA.CCA,
    Peyssonnelia.x = Peyssonnelia.Peyssonnelia,
    Sponge.x = Sponge.Sponge,
    Colonial_ascidian.x = Colonial_ascidian.Colonial_ascidian,
    Solitary_ascidian.x = Solitary_ascidian.Solitary_ascidian,
    Colonial_bryozoan.x = Colonial_bryozoan.Colonial_bryozoan,
    Serpulid_worm.x = Serpulid_worm.Serpulid_worm,
    Vermetid_snail.x = Vermetid_snail.Vermetid_snail,
    Bivalve.x = Bivalve.Bivalve,
    Hard_coral.x = Hard_coral.Hard_coral,
    Other.x = Other.Other)

str(results_list2)

results_list3 <- results_list2 %>%
  rename_with(~ gsub("\\.LCL$", "_LCL", .x), ends_with(".LCL"))

results_list3 <- results_list3 %>%
  rename_with(~ gsub("\\.UCL$", "_UCL", .x), ends_with(".UCL"))

results_list3_long <- results_list3 %>%
  pivot_longer(
    cols = everything(),  # Select all columns
    names_to = c("Binder", "Metric"),  # Create 'Binder' and 'metric' columns
    names_pattern = "(.*)\\.(x|response|SE|df|asymp_LCL|asymp_UCL)", 
    values_to = "Value"  # Then place the values in 'Value' column
  ) %>% as.data.frame()

str(results_list3_long)
results_list3_long <- results_list3_long %>% mutate_if(is.character, as.factor)
head(results_list3_long)

results_list3_long <- results_list3_long %>%
  group_by(Binder, Metric) %>%
  mutate(Replicate = row_number()) %>%
  ungroup()

results_list3_wide <- results_list3_long %>%
  pivot_wider(
    names_from = Metric,  
    values_from = Value   
  )

head(results_list3_wide)

results_list3_wide <-  results_list3_wide %>% as.data.frame()
head(results_list3_wide)

levels(results_list3_wide$Binder)

results_list3_wide <- results_list3_wide %>% mutate(Binder = factor(Binder,
                                                                    levels = c("Turf_algae", "Lobophora", "Other_macroalgae",
                                                                               "CCA", "Peyssonnelia",
                                                                               "Sponge", "Colonial_ascidian", "Solitary_ascidian", 
                                                                               "Colonial_bryozoan","Serpulid_worm", "Vermetid_snail", "Bivalve", "Hard_coral",         
                                                                               "Other")))

typecols2 <- c("Turf_algae" = "#B4EEB4",
               "Lobophora" = "#B8860B",
               "Other_macroalgae" = "#79B07A",
               "CCA" =  "#EEB4B4",
               "Peyssonnelia" = "#F5BC49",
               "Sponge" = "#FFF68F",     
               "Colonial_ascidian" = "#C799C6",
               "Solitary_ascidian" = "#B452CD",
               "Colonial_bryozoan" = "lightblue2",
               "Serpulid_worm" = "#FFDEAD",
               "Vermetid_snail" = "red",
               "Bivalve" = "darkgrey",
               "Hard_coral" = "peachpuff3",
               "Other" = "#919DBD")

(all.plot2 <- ggplot(results_list3_wide, aes(x = x, y = response, fill = Binder)) + 
    scale_y_continuous(limits = c(0,310), breaks = seq(0, 310, by = 10)) +  # Customize y-axis breaks
    scale_x_continuous(limits = c(0,7.2), breaks = seq(0, 7, by = 1)) +  # Customize y-axis breaks
    geom_hline(yintercept = seq(0, 310, by = 10), linetype = "dotted", color = "gray") + # Add dotted horizontal lines)
    geom_ribbon(aes(ymin = asymp_LCL, ymax = asymp_UCL), alpha = 0.15) +
    geom_line(aes(x = x, y = response, colour = Binder), size = 0.9) + 
    scale_colour_manual(name = "Binder", values = typecols2) +
    scale_fill_manual(name = "Binder", values = typecols2) +
    labs(x = "Bind area (cm2)", y = "Force to break bound rubble pair (N)") +
    theme_classic() +
    coord_cartesian(xlim = c(1, 7), ylim = c(0, 100)))  # Zoom in on both x and y axes








# Supp. 24  umax across 65-year time series at each site ----

umax <- read.csv("umax_heron_keppel.csv", header =TRUE)
str(umax_long)

umax$date <- as.Date(umax$date, format = "%d-%b-%Y")

# Reshape the data to long format with "site" and "value" columns
umax_long <- umax %>%
  pivot_longer(cols = 2:17, names_to = "site", values_to = "umax")

umax_exp_period <- umax_long %>%
  filter(date >= as.Date("2021-09-01") & date <= as.Date("2023-04-30"))
umax_long$site <- as.factor(umax_long$site)

exp_period_umax_summ <- umax_exp_period %>% group_by(site) %>%
  summarise(maxumax = max(umax, na.rm = TRUE), 
            meanumax = mean(umax, na.rm = TRUE))

view(exp_period_umax_summ)

ggplot(exp_period_umax_summ, aes(x = site, 
                                 y = meanumax)) + 
  geom_point() + theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=0.8))

ggplot(exp_period_umax_summ, aes(x = site, 
                                 y = maxumax)) + 
  geom_point() + theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=0.8))

# do probability density curves for each habitat

total_period_umax_summ <- umax_long %>% group_by(site) %>%
  summarise(maxumax = max(umax, na.rm = TRUE), 
            meanumax = mean(umax, na.rm = TRUE))

ggplot(total_period_umax_summ, aes(x = site, 
                                   y = meanumax)) + 
  geom_point() + theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=0.8))

ggplot(total_period_umax_summ, aes(x = site, 
                                   y = maxumax)) + 
  geom_point() + theme(axis.text.x = element_text(angle = 30, vjust = 0.9, hjust=0.8))

thresholds <- c(0.3000, 1.0000, 2.0000, 3.0000)

result <- umax_long %>%
  group_by(site) %>%
  summarise(
    total = n(),
    above_0.3 = sum(umax >= 0.3, na.rm = TRUE),
    above_1 = sum(umax >= 1, na.rm = TRUE),
    above_2 = sum(umax >= 2, na.rm = TRUE),
    above_3 = sum(umax >= 3, na.rm = TRUE),
    .groups = "drop"
  )

result <- result %>%
  mutate(
    prop_above_0.3 = above_0.3 / total,
    prop_above_1 = above_1 / total,
    prop_above_2 = above_2 / total,
    prop_above_3 = above_3 / total
  )

result <- result %>%
  mutate(
    pc_above_0.3 = prop_above_0.3*100,
    pc_above_1 = prop_above_1*100,
    pc_above_2 = prop_above_2*100,
    pc_above_3 = prop_above_3*100,
    days_above_1_in_4months = prop_above_1*120,
  )

print(result)

# On 30% of day = 1 in every 3 days it goes above 1 m/s
# On 6% of days = 7.2 times in 4 months










# Figure 3a ------

# BREAK FORCE & VELOCITY SAME PLOT

range(gridpairbindpiecemainSS2_gc$break_force, na.rm=TRUE)
range(gridpairbindpiecemainSS2_gc$break_velocity_BINDAREA_LT_org, na.rm=TRUE)

# Define a scaling factor to align the y-values if needed (adjust based on your data).

# New scaling factor for smaller velocity range
scale_factor <- 0.085

legend_titleT <- "Habitat"

colsT <- c("In Shelt Shal" = "#6CA6CD",
           "In Exp Shal" = "#CD8C95", 
           "Off Shelt RF" = "#AEE2F5",
           "Off Shelt Shal" = "#6CA6CD",
           "Off Shelt Deep" = "#104E8B",
           "Off Exp Shal" = "#CD8C95",    
           "Off Exp Deep" = "#A52A2A")

linesT <- c("In Shelt Shal" = 2,
            "In Exp Shal" = 2, 
            "Off Shelt RF" = 1,
            "Off Shelt Shal" = 1,
            "Off Shelt Deep" = 1,
            "Off Exp Shal" = 1,    
            "Off Exp Deep" = 1)

shapesT <- c("In Shelt Shal" = 8,
             "In Exp Shal" = 8, 
             "Off Shelt RF" = 17,
             "Off Shelt Shal" = 17,
             "Off Shelt Deep" = 17,
             "Off Exp Shal" = 17,    
             "Off Exp Deep" = 17)

combined_plot <- ggplot() +
  
  # Plot for Break Force (Primary y-axis)
  geom_pointrange(aes(y = response, x = months_deployedFO, 
                      colour = locexpdepth, fill = locexpdepth, shape = locexpdepth,
                      ymin = asymp.LCL, ymax = asymp.UCL),
                  data = freqmod1_preds, position = pd, size = 0.8) +
  geom_point(aes(y = break_force, x = months_deployedFO, colour = locexpdepth),
             data = gridpairbindpiecemainSS2, position = position_jitter(width = 0.2, height = 0), 
             size = 0.5, alpha = 0.1) +
  geom_line(aes(y = response, x = months_deployedFO, group = locexpdepth, colour = locexpdepth),
            data = freqmod1_preds, position = pd) +
  
  # Plot for Velocity (Secondary y-axis)
  geom_pointrange(aes(y = response / scale_factor, x = months_deployedFO, 
                      colour = locexpdepth, fill = locexpdepth, shape = locexpdepth,
                      ymin = asymp.LCL / scale_factor, ymax = asymp.UCL / scale_factor),
                  data = velmod_preds, position = pd, size = 0.8, linetype = "dashed") +
  geom_point(aes(y = break_velocity_BINDAREA_LT_org / scale_factor, x = months_deployedFO, colour = locexpdepth),
             data = gridpairbindpiecemainSS2_gc, position = position_jitter(width = 0.2, height = 0), 
             size = 0.5, alpha = 0.1, shape = 17) +
  geom_line(aes(y = response / scale_factor, x = months_deployedFO, group = locexpdepth, colour = locexpdepth),
            data = velmod_preds, position = pd, linetype = "dotted") +
  
  # Customizing the primary and secondary y-axes
  scale_y_continuous(name = "Force to break binds (N)",
                     limits = c(0, 50),
                     sec.axis = sec_axis(~ . * scale_factor, name = "Velocity to break apart bound rubble pairs (m/s)")) +
  
  # Customize the x-axis and add common aesthetics
  xlab("Months since deployment") +
  scale_colour_manual(legend_titleT, values = colsT) +
  scale_fill_manual(legend_titleT, values = colsT) +
  scale_shape_manual(legend_titleT, values = shapesT) +
  
  # Styling
  theme_classic() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.title.y.right = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.grid.major.y = element_line(color = "black", linetype = "dotted", size = 0.3),
        legend.position = "none")

combined_plot







# Figure 4 ------

# Using dataframe from Top contributors analysis

#gridbinderB_wideTotNoNAS

names(gridbinderB_wideTotNoNAS)
class(gridbinderB_wideTotNoNAS$months_deployed)
class(gridbinderB_wideTotNoNAS$months_deployedFO)

#Makes the community matrix (the response variables) for surface grids only
commat2 <- gridbinderB_wideTotNoNAS %>% 
  dplyr::select("Turf_algae", "Lobophora", "Other_macroalgae",
                "CCA", "Peyssonnelia",
                "Sponge", "Colonial_ascidian", "Solitary_ascidian", 
                "Colonial_bryozoan","Serpulid_worm", "Vermetid_snail", "Bivalve", "Hard_coral",         
                "Other")

#Makes the explanatory variable matrix - remember this doesn't have topside and underside
commat2met <- gridbinderB_wideTotNoNAS %>% 
  dplyr::select(unique_pair_id, grid_id, locexpdepth, site_unique, months_deployedFO, months_deployed)

tr.commat2 <- sqrt(commat2) # transform the data with sq root

commat2.dist <-vegdist(tr.commat2, method = "bray") # now calculate the distance matrix using bray curtis

# All months

#This produce the coordinates for each point

allmonth <- metaMDS(tr.commat2, distance = "bray") # metaMDS on the transformed data, not the raw data

allmonth$stress

stressplot(allmonth) # stress is 0.18 which is below 0.2 so not a poor fit, but not a great one.

sites.scores<- as.data.frame(scores(allmonth, display = 'sites'))
sites.scores

sites.scores<- data.frame(sites.scores, gridbinderB_wideTotNoNAS)
head(sites.scores)

species.scores<- as.data.frame(scores(allmonth, display = 'species'))

species.scores$Binders<- row.names(species.scores)
head(species.scores)

# Want plot where I move from month to month, and then also show the community vectors in each month

sites.scores$locexpdepthMonth <- paste(sites.scores$locexpdepth, 
                                       sites.scores$months_deployedFO, sep="_")

sites.scores$locexpdepthMonth <- as.factor(sites.scores$locexpdepthMonth)

levels(sites.scores$locexpdepthMonth)

sites.scores$locexpdepthMonth <- 
  factor(sites.scores$locexpdepthMonth, 
         levels = c("In Shelt Shal_4", "In Shelt Shal_12", "In Shelt Shal_18", 
                    "In Exp Shal_4", "In Exp Shal_12","In Exp Shal_18",
                    "Off Shelt RF_4","Off Shelt RF_8","Off Shelt RF_12","Off Shelt RF_18",
                    "Off Shelt Shal_4","Off Shelt Shal_8","Off Shelt Shal_12","Off Shelt Shal_18",
                    "Off Shelt Deep_4","Off Shelt Deep_8","Off Shelt Deep_12","Off Shelt Deep_18",
                    "Off Exp Shal_4","Off Exp Shal_8","Off Exp Shal_12","Off Exp Shal_18",  
                    "Off Exp Deep_4", "Off Exp Deep_8", "Off Exp Deep_12","Off Exp Deep_18"),
         labels = c("ISS4","ISS12","ISS18",
                    "IES4","IES12","IES18",
                    "ORF4","ORF8","ORF12","ORF18",
                    "OSS4","OSS8","OSS12","OSS18",
                    "OSD4","OSD8","OSD12","OSD18",
                    "OES4","OES8","OES12","OES18",
                    "OED4","OED8","OED12","OED18"))


colsT <- c("ISS4" = "#6CA6CD",
           "ISS12" = "#6CA6CD",
           "ISS18" = "#6CA6CD",
           "IES4" = "#CD8C95", 
           "IES12" = "#CD8C95",
           "IES18" = "#CD8C95",
           "ORF4" = "#AEE2F5",
           "ORF8" = "#AEE2F5",
           "ORF12" = "#AEE2F5",
           "ORF18" = "#AEE2F5",
           "OSS4" = "#6CA6CD",
           "OSS8" = "#6CA6CD",
           "OSS12" = "#6CA6CD",
           "OSS18" = "#6CA6CD",
           "OSD4" = "#104E8B",
           "OSD8" = "#104E8B", 
           "OSD12" = "#104E8B",
           "OSD18" = "#104E8B",
           "OES4" = "#CD8C95",
           "OES8" = "#CD8C95",
           "OES12" = "#CD8C95",
           "OES18" = "#CD8C95",
           "OED4" = "#A52A2A",
           "OED8" = "#A52A2A",
           "OED12" = "#A52A2A",
           "OED18" = "#A52A2A")

shapesT <- c("ISS4" = 8,
             "ISS12" = 8,
             "ISS18" = 8,
             "IES4" = 8, 
             "IES12" = 8,
             "IES18" = 8,
             "ORF4" = 17,
             "ORF8" = 17,
             "ORF12" = 17,
             "ORF18" = 17,
             "OSS4" = 17,
             "OSS8" = 17,
             "OSS12" = 17,
             "OSS18" = 17,
             "OSD4" = 17,
             "OSD8" = 17, 
             "OSD12" = 17,
             "OSD18" = 17,
             "OES4" = 17,
             "OES8" = 17,
             "OES12" = 17,
             "OES18" = 17,
             "OED4" = 17,
             "OED8" = 17,
             "OED12" = 17,
             "OED18" = 17)

centroidsA <- aggregate(cbind(NMDS1, NMDS2) ~ locexpdepthMonth, data = sites.scores, FUN = mean)
centroidsA <- data.frame(centroidsA) # pulling out just the averages which is what I will plot
legend_title <- "Habitat & Month"

levels(sites.scores$locexpdepthMonth)

# plot for manuscript

library(ggrepel)

g1<-ggplot()+
  geom_segment(data=NULL, aes(y=-Inf, x=0, yend=Inf, xend=0), linetype='dotted')+
  geom_segment(data=NULL, aes(y=0, x=-Inf, yend=0, xend=Inf), linetype='dotted')+
  #geom_segment(data=topsideshallnmds.species.scores, aes(y=0, x=0, yend=NMDS2, xend=NMDS1), colour = "grey", arrow=arrow(length=unit(0.3, 'lines')))+ #adding the binding organism drivers
  #geom_text(data=topsideshallnmds.species.scores, aes(y=NMDS2, x=NMDS1, label = Binders, hjust="inward",vjust="inward"), show.legend=FALSE, colour = "grey") +
  geom_point(data=centroidsA, aes(y=NMDS2, x=NMDS1, color=locexpdepthMonth, shape = locexpdepthMonth, fill = locexpdepthMonth), size = 5)+
  scale_colour_manual(legend_title, values=colsT) +
  scale_fill_manual(legend_title, values = colsT) +
  scale_shape_manual(legend_title, values = shapesT) +
  geom_text_repel(data=centroidsA, aes(y=NMDS2, x=NMDS1, label = locexpdepthMonth,
                                       hjust="inward",vjust="inward", color=locexpdepthMonth), show.legend=FALSE) +
  theme_classic() +
  ylab("NMDS2") +
  xlab("NMDS1") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_vline(xintercept=0, linetype="dashed") + 
  geom_hline(yintercept=0, linetype="dashed") +
  theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position = "none") 
g1

# nMDS segment plot

gtopsideSEG <-ggplot()+
  geom_segment(data=NULL, aes(y=-Inf, x=0, yend=Inf, xend=0), linetype='dotted')+
  geom_segment(data=NULL, aes(y=0, x=-Inf, yend=0, xend=Inf), linetype='dotted')+
  geom_segment(data=species.scores, aes(y=0, x=0, yend=NMDS2, xend=NMDS1), 
               colour = "grey", arrow=arrow(length=unit(0.3, 'lines'))) + 
  geom_text_repel(data=species.scores, aes(y=NMDS2, x=NMDS1, 
                                           #hjust="inward",vjust="inward",
                                           label = Binders), colour = "gray27") + 
  #geom_point(data=topsideshallnmds.sites.scores, aes(y=NMDS2, x=NMDS1), color="white") +
  #geom_point(data=centroidsA, aes(y=NMDS2, x=NMDS1, colour = aspectmonthdepth)) +
  #scale_colour_manual(legend_title, values=cols) +
  theme_classic() +
  ylab("NMDS2") +
  xlab("NMDS1") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  geom_vline(xintercept=0, linetype="dashed") + 
  geom_hline(yintercept=0, linetype="dashed") +
  theme(text = element_text(size=15, face="bold")) + 
  theme(axis.text = element_text(colour = "black")) + 
  theme(panel.background =element_rect(colour = "black", size=1)) + 
  theme(axis.ticks.length=unit(.2,"cm")) +
  theme(legend.position = "none") # don't need the legend
gtopsideSEG

grid.arrange (g1, gtopsideSEG, nrow = 2)
# 10 x 8.14 inches with and without legend - adjust in Illustrator







# Results first section: Summary stats -----

str(gridpairbindpiecemainSS2_gc)

gridpairbindpiecemainSS2_gcSUM <- gridpairbindpiecemainSS2_gc %>%
  summarise(avgsumcont = mean(sumcontact, na.rm=TRUE),
            sdsumcont = sd(sumcontact, na.rm=TRUE),
            sesumcont = sdsumcont/sqrt(n()))

gridpairbindpiecemainSS2_gcSUM3 <- gridpairbindpiecemainSS2_gc %>%
  summarise(avgsumgap = mean(sumgap, na.rm=TRUE),
            sdsumgap = sd(sumgap, na.rm=TRUE),
            sesumgap = sdsumgap/sqrt(n()))

gridpairbindpiecemainSS2_gc <- gridpairbindpiecemainSS2_gc %>%
  mutate(propsumcont = sumcontact/avgrubblelength)

gridpairbindpiecemainSS2_gc <- gridpairbindpiecemainSS2_gc %>%
  mutate(propsumgap = sumgap/avgrubblelength)

gridpairbindpiecemainSS2_gcSUM2<- gridpairbindpiecemainSS2_gc %>%
  summarise(avgpropsumcont = mean(propsumcont, na.rm=TRUE),
            sdpropsumcont = sd(propsumcont, na.rm=TRUE),
            sepropsumcont = sdpropsumcont/sqrt(n()))

gridpairbindpiecemainSS2_gcSUM2<- gridpairbindpiecemainSS2_gc %>%
  summarise(count = mean(propsumcont, na.rm=TRUE),
            sdpropsumcont = sd(propsumcont, na.rm=TRUE),
            sepropsumcont = sdpropsumcont/sqrt(n()))

gridpairbindpiecemainSS2_gcSUM4<- gridpairbindpiecemainSS2_gc %>%
  summarise(avgpropsumgap = mean(propsumgap, na.rm=TRUE),
            sdpropsumgap = sd(propsumgap, na.rm=TRUE),
            sepropsumgap = sdpropsumgap/sqrt(n()))

hist(gridpairbindpiecemainSS2_gc$max_gap_span)

gridpairbindpiecemainSS2_gcSUM6<- gridpairbindpiecemainSS2_gc %>%
  summarise(avgmax_gap_span = mean(max_gap_span, na.rm=TRUE),
            maxmax_gap_span = max(max_gap_span, na.rm=TRUE),
            medianmax_gap_span = median(max_gap_span, na.rm=TRUE),
            sdmax_gap_span = sd(max_gap_span, na.rm=TRUE),
            semax_gap_span = sdmax_gap_span/sqrt(n()))

sum(gridpairbindpiecemainSS2_gc$boundB == 1, na.rm=TRUE)

sum(gridpairbindpiecemainSS2_gc$boundB == 0, na.rm=TRUE)










