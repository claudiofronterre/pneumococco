# Load required packages and functions
if (!require("pacman")) install.packages("pacman")
pkgs = c("lubridate", "dplyr") # package names
pacman::p_load(pkgs, character.only = T)
  
# Load and clead data and generate variable of interest
cocco <- readr::read_csv("data/PCVPA_31dec_adj_claudio_v2_16feb19_min.csv")
cocco <- cocco[cocco$typ_cat7 != "adult", ] # remove adults

# Generate dummy variabile for vt carriage (1 = carriage, 0 = no carriage)
cocco$vt_01 <- ifelse(cocco$vt_01 == "vt", 1, 0)

# Generate dummy variabile for vaccination status (1 = vaccinated, 0 = no vaccinated)
cocco$pcv_vaccd <- ifelse(cocco$pcv_vaccd == "PCV13-vacc'd", 1, 0)
cocco$pcv_vaccd[is.na(cocco$pcv_vaccd)] <- 0

# Calculate time since vaccination and age at vaccination
cocco <- cocco %>% 
  mutate_at(vars(contains("dat")), as.Date, format = "%d%b%Y") %>% # change to data format
  mutate(t1 = as.numeric(collection_date - pcv1dat), # time since first vacc
         age_y2 = as.numeric(collection_date - as.Date(dob, "%d%b%Y")) / 365, # child age at collection
         agev = age_y2 - t1 / 365) # age at vaccination


# Recode risk season
cocco$risk <- as.factor(ifelse(month(cocco$collection_date) > 4 & 
                               month(cocco$collection_date) < 9,
                               yes = "High", no = "Low"))


# Save the cleaned data set
saveRDS(cocco, "data/cocco_clean.rds")

