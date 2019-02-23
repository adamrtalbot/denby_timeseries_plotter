library(tidyr)
library(dplyr)
library(tibble)
library(readr)
library(RSQLite)

# Open connection to new datbase

timeseries_db <- dbConnect(drv = RSQLite::SQLite(), "timeseries-db.sqlite")

## Import botrytis data ##

let.data <- read.csv("data/lettuce_gene_count_sizefactored_wideformat.csv") %>%
  gather(key = treatment, value = counts, -X)
bot.data <- read.csv("data/bot_gene_count_sizefactored_wideformat.csv")[,1:49] %>%
  gather(key = treatment, value = counts, -X)

let.vst <- read.csv("data/lettuce_gene_count_sizefactored_wideformat_vst.csv") %>%
  gather(key = treatment, value = vst, -X) %>% mutate(vst = vst-min(vst))
bot.vst <- read.csv("data/bot_gene_count_sizefactored_wideformat_vst.csv")[,1:49] %>%
  gather(key = treatment, value = vst, -X) %>% mutate(vst = vst-min(vst))

vst.data <- rbind(let.vst, bot.vst) %>%
  separate(col = treatment, into = c("treatment", "timepoint", "replicate"), sep = "_") %>%
  select(gene = X, everything()) %>%
  mutate(timepoint = as.numeric(gsub(pattern = "tp", replacement = "", x = .$timepoint))) %>%
  mutate(replicate = gsub(pattern = "rep", replacement = "", x = .$replicate))

let.bot.data <- rbind(let.data, bot.data) %>%
  separate(col = treatment, into = c("treatment", "timepoint", "replicate"), sep = "_") %>%
  select(gene = X, everything()) %>%
  mutate(timepoint = as.numeric(gsub(pattern = "tp", replacement = "", x = .$timepoint))) %>%
  mutate(replicate = gsub(pattern = "rep", replacement = "", x = .$replicate)) %>%
  left_join(vst.data, by = c("gene", "treatment", "timepoint", "replicate")) %>%
  mutate(log2 = log(counts+1, base = 2)) %>%
  mutate(pathosystem = "let-bot")

remove(let.vst)
remove(bot.vst)
remove(let.data)
remove(bot.data)
remove(vst.data)

dbWriteTable(timeseries_db, "letbot", let.bot.data, overwrite = TRUE)

## Import sclerotinia-lettuce data ##

# Lettuce VST data

let.sclero.vst.data <- read.csv("data/lettuce_only_factored_vst.csv") %>%
  select("gene" = "X", everything()) %>% 
  gather(key = treatment, value = vst, -gene) %>%
  separate(col = treatment, into = c("treatment", "timepoint", "replicate"), sep = "_") %>%
  mutate(timepoint = as.numeric(gsub(pattern = "T", replacement = "", x = .$timepoint))) %>%
  mutate(replicate = gsub(pattern = "Rep", replacement = "", x = .$replicate)) %>%
  mutate(treatment = gsub(pattern = "Control", replacement = "Mock", x = .$treatment)) %>%
  mutate(pathosystem = "let-sclero") %>%
  mutate(vst = vst-min(vst))

# Lettuce count data

let.sclero.data <- read.csv(file = "timeseries_plotter/data/lettuce_sclero_only_factored_new.csv") %>%
  select("gene" = "X", everything()) %>% 
  gather(key = treatment, value = counts, -gene) %>%
  separate(col = treatment, into = c("treatment", "timepoint", "replicate"), sep = "_") %>%
  mutate(timepoint = as.numeric(gsub(pattern = "T", replacement = "", x = .$timepoint))) %>%
  mutate(replicate = gsub(pattern = "Rep", replacement = "", x = .$replicate)) %>%
  mutate(treatment = gsub(pattern = "Control", replacement = "Mock", x = .$treatment)) %>%
  mutate(pathosystem = "let-sclero") %>%
  left_join(y = let.sclero.vst.data, by = c("gene", "treatment", "timepoint", "replicate", "pathosystem"))


# Sclerotinia

sclero.data <- read.csv(file = "timeseries_plotter/data/sclerotinia_only_sizefactored.csv") %>%
  select("gene" = "X", everything()) %>% 
  gather(key = treatment, value = counts, -gene) %>%
  separate(col = treatment, into = c("treatment", "timepoint", "replicate"), sep = "_") %>%
  mutate(timepoint = as.numeric(gsub(pattern = "T", replacement = "", x = .$timepoint))) %>%
  mutate(replicate = gsub(pattern = "Rep", replacement = "", x = .$replicate)) %>%
  mutate(treatment = gsub(pattern = "Control", replacement = "Mock", x = .$treatment)) %>%
  mutate(pathosystem = "let-sclero") %>%
  mutate(vst = NA)

let.sclero.data <- rbind.data.frame(let.sclero.data, sclero.data) %>% 
  mutate(log2 = log(counts+1, base = 2))

remove(sclero.data)
remove(let.sclero.vst.data)

dbWriteTable(timeseries_db, "letsclero", let.sclero.data, overwrite = TRUE)

## Import Arabidopsis data ##

athal.data <- read_delim("~/Documents/HAPI Lettuce Work Folder/Experiments/Experiment 17 - Botrytis Timeseries/plotter/timeseries_plotter/data/Botrytis-Combo3.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  select(gene = X1, everything()) %>%
  gather(key = treatment, value = expression, -gene) %>%
  separate(col = treatment, into = c("timepoint", "treatment", "replicate"), sep = ":") %>%
  mutate(timepoint = as.numeric(gsub(pattern = " hours", replacement = "", x = .$timepoint))) %>%
  mutate(pathosystem = "athal-bot") %>%
  left_join(y = read_csv("~/Documents/HAPI Lettuce Work Folder/Experiments/Experiment 17 - Botrytis Timeseries/plotter/timeseries_plotter/data/CATMA_probes_v4.csv"), by = c("gene" = "CATMA")) %>% 
  select(-gene) %>% 
  select(gene = ATG, everything())

dbWriteTable(timeseries_db, "athalbot", athal.data, overwrite = TRUE)

# SQL queries:
# dbGetQuery(timeseries_db, 'SELECT * FROM athalbot WHERE gene="AT5G39610"')
# dbGetQuery(timeseries_db, 'SELECT * FROM athalbot WHERE gene in ("AT5G39610", "AT1G61110")')

# test <- dbGetQuery(timeseries_db, 'SELECT * FROM athalbot WHERE gene=:x', 
#         params = list(x = "AT5G39610"))

# test <- dbGetQuery(timeseries_db, 'SELECT * FROM athalbot WHERE gene in (:x)', 
#         params = list(x = c("AT5G39610", "AT1G61110")))

# Parameterised SQL queries;
pathosystem <- "letbot"
rs <- dbSendQuery(timeseries_db, paste0('SELECT * FROM ', pathosystem, ' WHERE gene in (:x)'))
dbBind(rs, param = list(x = "Lsat_1_v5_gn_4_1201"))
test <- dbFetch(rs)

# rs <- dbSendQuery(timeseries_db, paste('SELECT * FROM', 'athalbot', 'WHERE gene in (:x)'))
# dbBind(rs, param = list(x = "AT5G39610"))
# test <- dbFetch(rs)

## Combine data tables ##

#data <- rbind.data.frame(let.bot.data, let.sclero.data, let.data)
#rm(let.bot.data)
#rm(let.data)
#rm(let.sclero.data)
#rm(athal.data)

## Save workspace for importing ##

#dir.create("timeseries_plotter/data")
#save(list = "data", file = "timeseries_plotter/data/timeseries.Rdata", compress = FALSE)

#dbDisconnect()
