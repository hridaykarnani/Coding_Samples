# ---------------------------------------------------------------------------- #
# ----------- Panel regression by country: Local projection Method ----------- #
# ---------------------------------------------------------------------------- #

# Hriday Karnani 
# Date: March-2022

# Introduction to this code:
# In this code we estimate the interest rate response to fiscal stimulus (IRRF).
# This is estimated following a local projection methodology using panel data of 
# countries. We start by cleaning the data and generating some new variables, 
# later we estimate the IRRF, Policy RRF and Consumption RRF. We do some plots
# to show the main results.
# There are different specifications for the local projection regression, this 
# can be defined in "locals" section. We also create different folders to store
# results.
# This code is one of the main codes of a paper published on the EER, 
# "Saving Constraints, inequality, and the credit market response to fiscal stimulus"
# ---------------------------------------------------------------------------- #
options(scipen=999)
rm(list=ls())

# Install and load packages
list.of.packages <- c("haven","dplyr", "readxl", "data.table","mFilter", "plm","texreg",
                      "zoo","sandwich","lmtest","ggplot2","ggthemes","doBy","estimatr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
rm(list.of.packages,new.packages) # con esto borro los objetos

# load libraries
Packages <- c("haven","dplyr", "readxl", "data.table","mFilter", "plm","texreg",
                      "zoo","sandwich","lmtest","ggplot2","ggthemes","doBy","estimatr")
lapply(Packages, library, character.only = TRUE)

# ---------------------------------------------------------------------------- #
# Set working directory
wd <- "C:/Users/Hriday/Dropbox/Hriday/IRRF/Replication"
setwd(wd)
rdata1 <- paste0(wd,"/data/")
results1 <- paste0(wd,"/results/")
figures1 <- paste0(wd,"/tablas_figuras/")
tex1 <- paste0(wd,"/tex/")

# ---------------------------------------------------------------------------- #

# ------------------------------ Set "locals" ---------------------------------- #
# Filter
#  "HP" / "BK"
filter <- "BK"

# Lags for BP method (AG method will use lags_BP/2)
# 1:2 = from 1 to 2 / 1:4 = from 1 to 4 / 1:6 = from 1 to 6
lags <- 1:4

# Definition of IRRF
# 4 = X periods of impact (depends on lags) / 0 = only impact period
def_irrf <- 4

# Determine if we include policy rate as an indep. variable in regressions
# yes = yes in AG and BP / no = no in AG and BP, / BP = yes in BP no in AG
inc.polrate <- "BP"
# ---------------------------------------------------------------------------- #
# Create directories to save results for different possible specifications in locals
rdata <- paste0(rdata1,filter,"_",length(lags),"lags_",def_irrf,"per/")
results <- paste0(results1,filter,"_",length(lags),"lags_",def_irrf,"per/")
figures <- paste0(figures1,filter,"_",length(lags),"lags_",def_irrf,"per/")
tex <- paste0(tex1,filter,"_",length(lags),"lags_",def_irrf,"per/")

dir.create(rdata)
dir.create(results)
dir.create(figures)
dir.create(tex)

rdata <- paste0(rdata1,filter,"_",length(lags),"lags_",def_irrf,"per/",inc.polrate,"polrate/")
results <- paste0(results1,filter,"_",length(lags),"lags_",def_irrf,"per/",inc.polrate,"polrate/")
figures <- paste0(figures1,filter,"_",length(lags),"lags_",def_irrf,"per/",inc.polrate,"polrate/")
tex <- paste0(tex1,filter,"_",length(lags),"lags_",def_irrf,"per/",inc.polrate,"polrate/")
rowname <- paste0(filter,"_",length(lags),"lags_",def_irrf,"per_",inc.polrate,"polrate")

dir.create(rdata)
dir.create(results)
dir.create(figures)
dir.create(tex)

write.table(rdata, file = paste0(rdata1,"/rdata.txt"), sep = "",row.names = F)
write.table(results, file = paste0(results1,"/results.txt"), sep = "",row.names = F)
write.table(tex, file = paste0(tex1,"/tex.txt"), sep = "",row.names = F)
write.table(figures, file = paste0(figures1,"/figures.txt"), sep = "",row.names = F)
write.table(rowname, file = paste0(results1,"/rownames.txt"), sep = "",row.names = F)

lags.1 <- 1/length(lags)
aglags <- 1:(length(lags)/2)
ag.lags1 <- 1/length(aglags)
# ---------------------------------------------------------------------------- #
# Clean and work data

paneldata <- read_xlsx("data/panel_data_july21.xlsx")

# change variable name
setnames(paneldata,old="bis_cons_credit_gdp",new="credit_gdp")

# set variables as numeric
columns <- c("quarter","year","real_credit","rgdp_cyc","policy_rate","bond_yields","rgdp","inflation","c_real")
paneldata <- paneldata %>% mutate_at(columns, as.numeric)

# set NA when value is 0
paneldata$policy_rate <- ifelse(paneldata$policy_rate == 0 | is.nan(paneldata$policy_rate) & 
                                    paneldata$year < 2000, NA,paneldata$policy_rate)
paneldata$credit_gdp <- ifelse(paneldata$credit_gdp == 0, NA,paneldata$credit_gdp)
paneldata$real_credit <- ifelse(paneldata$real_credit == 0, NA,paneldata$real_credit)
paneldata$real_credit_cyc <- ifelse(paneldata$real_credit_cyc == 0, NA,paneldata$real_credit_cyc)
paneldata$rgdp[is.nan(paneldata$rgdp)] <- NA

# generate g spread
paneldata$gspread <- paneldata$bond_yields - paneldata$policy_rate

# order by country
paneldata <- paneldata[order(paneldata$country),]
# keep data until 2007
paneldata <- paneldata %>% filter(year<= 2007)
# generate copy (it will be used for AG later)
AGpaneldata <- paneldata
# read some country variables
topanel <- read_dta("data/topanel.dta")
# merge df
mpanel <- merge(paneldata,topanel,by = "country")

# drop estonia
mpanel <- mpanel %>% filter(country != "Estonia")

# generate new variables
mpanel$c_real <- mpanel$c_real/1000
mpanel$lng <- log(mpanel$G)
mpanel$lngdp <- log(mpanel$rgdp)
mpanel$lnc <- log(mpanel$c_real)

# Apply HP filter and BK filter to Real GDP. Keep trend
countries <- unique(mpanel$country)
mpanelsub3 <- list();mpanelsub4<-NA
for (i in 1:28) {
# test; i=1
    mpanelsub <- mpanel %>% filter(country == countries[i])

    # Fix NA problem
    df1 <- mpanelsub %>% filter(!is.na(mpanelsub$rgdp)) %>% 
            select(country,year,quarter,rgdp)
    df2 <-  hpfilter(df1$rgdp,freq=1600,drift=TRUE)
    df3 <- data.frame(df1, cycle_hp = df2$cycle)
    df3$trend_rgdp <- df3$rgdp- df3$cycle_hp
    
    mpanelsub2 <- merge(mpanelsub,df3,by = c("country","year","quarter","rgdp"),all = T,)
    
    df2 <- bkfilter(df1$rgdp,pl = 6,pu=32,nfix = 12,drift = T)
    df3 <- data.frame(df1, cycle_bk = df2$cycle)
    df3$trend_rgdpBK <- df3$rgdp - df3$cycle_bk
    mpanelsub3[[i]] <- merge(mpanelsub2,df3,by = c("country","year","quarter","rgdp"),all = T)
    mpanelsub4 <- rbind(mpanelsub4,mpanelsub3[[i]])
    
}

mpanelsub4$trend_rgdpBK <- ifelse(is.na(mpanelsub4$trend_rgdpBK), mpanelsub4$trend_rgdp, mpanelsub4$trend_rgdpBK)
# Create some news var based on trend real GDP
mpanelsub4 <-mpanelsub4 %>% mutate(dlng = log((G)/trend_rgdpBK), dlnc = log((c_real)/trend_rgdpBK),
                                   dlngdp = log((rgdp)/trend_rgdpBK))

if (filter == "HP") {mpanelsub4 <- mpanelsub4 %>% mutate(dlng = log((G)/trend_rgdp), dlnc = log((c_real)/trend_rgdp),
                                        dlngdp = log((rgdp)/trend_rgdp))
} else {mpanelsub4 <- mpanelsub4 %>% mutate(dlng = log((G)/trend_rgdpBK), dlnc = log((c_real)/trend_rgdpBK),
                                        dlngdp = log((rgdp)/trend_rgdpBK))}

# Create function to compute lags
L. <- function(x, k) {
    res <- as.matrix(dplyr::bind_cols(lapply(k, function(k) dplyr::lag(x, k))))
    colnames(res) <- paste0("_lag", seq_along(1:ncol(res)))
    res
}

# Create function that determines the variables of the model
var.model <- function(dep.lags, pol.rate) {
    if (pol.rate == "no") {
        model.var <- shift(bond_yields,dep.lags) ~ dlng + L.(bond_yields,lags) + L.(dlngdp,lags) + L.(dlng,lags)
    } else {
        model.var <- shift(bond_yields,dep.lags) ~ dlng + L.(bond_yields,lags) + L.(policy_rate,lags) + L.(dlngdp,lags) + L.(dlng,lags)
    }}

# Generate IRRF by country
irrf_byc_f0 <- list();irrf_byc_f1 <- list();irrf_byc_f2 <- list();irrf_byc_f3 <- list()
for (i in 1:28) {
    # BOND YIELD F0 (impact on the same period)
   data <- mpanelsub4 %>% filter(country == countries[i]) # keep data of country i
   reg <- lm_robust(var.model(0,inc.polrate),data = data,se_type = "HC1") # reg for country i
   reg2 <- c(reg$coefficients[2],reg$std.error[2]) # store results and keep uper and lower bound (next line of code)
   irrf_byc_f0[[i]] <- data.frame(estimates = reg2,low_bound = reg2[1] - reg2[2]*1.645,upp_bound = reg2[1] + reg2[2]*1.645) 
   # BOND YIELD F1   (impact on the next period, quarter)
   reg <- lm_robust(var.model(-1,inc.polrate),data = data,se_type = "HC1")
   reg2 <- c(reg$coefficients[2],reg$std.error[2])
   irrf_byc_f1[[i]] <- data.frame(estimates = reg2,low_bound = reg2[1] - reg2[2]*1.645,upp_bound = reg2[1] + reg2[2]*1.645) 
   # BOND YIELD F2
   reg <- lm_robust(var.model(-2,inc.polrate),data = data,se_type = "HC1")
   reg2 <- c(reg$coefficients[2],reg$std.error[2])
   irrf_byc_f2[[i]] <- data.frame(estimates = reg2,low_bound = reg2[1] - reg2[2]*1.645,upp_bound = reg2[1] + reg2[2]*1.645) 
   # BOND YIELD F3
   reg <- lm_robust(var.model(-3,inc.polrate),data = data,se_type = "HC1")
   reg2 <- c(reg$coefficients[2],reg$std.error[2])
   irrf_byc_f3[[i]] <- data.frame(estimates = reg2,low_bound = reg2[1] - reg2[2]*1.645,upp_bound = reg2[1] + reg2[2]*1.645) 
}


# Save in one dataframe the IRRF of every country. Weight equally each response (F0, F1, F2 and F3)
IRRF <- NA;IRRF_lb <- NA;IRRF_ub <- NA;IRRF_at0 <- NA;IRRF_at0_lb <- NA;IRRF_at0_ub <- NA
IRRF_adj <- NA;IRRFat0_adj <- NA;weights <- NA;weightsat0 <- NA;IRRF_f <- matrix(NA,28,4)
for (i in 1:28) {
    if (def_irrf == 4) {
        # Weight equally each response (F0,F1,F2,F3), upper and lower bound
        IRRF[i] <- lags.1*(irrf_byc_f0[[i]][1,1]+irrf_byc_f1[[i]][1,1]+irrf_byc_f2[[i]][1,1]+irrf_byc_f3[[i]][1,1])
        IRRF_lb[i] <- lags.1*(irrf_byc_f0[[i]][1,2]+irrf_byc_f1[[i]][1,2]+irrf_byc_f2[[i]][1,2]+irrf_byc_f3[[i]][1,2])
        IRRF_ub[i] <- lags.1*(irrf_byc_f0[[i]][1,3]+irrf_byc_f1[[i]][1,3]+irrf_byc_f2[[i]][1,3]+irrf_byc_f3[[i]][1,3])
    } else {
        IRRF[i] <- lags.1*irrf_byc_f0[[i]][1,1]
        IRRF_lb[i] <- lags.1*irrf_byc_f0[[i]][1,2]
        IRRF_ub[i] <- lags.1*irrf_byc_f0[[i]][1,3]
    }

    IRRF_at0[i] <- lags.1*irrf_byc_f0[[i]][1,1]
    IRRF_at0_lb[i] <- lags.1*irrf_byc_f0[[i]][1,2]
    IRRF_at0_ub[i] <- lags.1*irrf_byc_f0[[i]][1,3]

    IRRF_adj[i] <- IRRF[i]/(IRRF_ub[i] - IRRF_lb[i])
    weights[i] <- 1/(IRRF_ub[i] - IRRF_lb[i]) # compute weight in case WLS is required later
    IRRFat0_adj[i] <- IRRF_at0[i]/(IRRF_at0_ub[i] - IRRF_at0_lb[i])
    weightsat0[i] <- 1/(IRRF_at0_ub[i] - IRRF_at0_lb[i])
    IRRF_f[i,] <- cbind(irrf_byc_f0[[i]][1,1],irrf_byc_f1[[i]][1,1],irrf_byc_f2[[i]][1,1],irrf_byc_f3[[i]][1,1])
}

IRRF_df <- data.frame(cbind(IRRF, IRRF_lb, IRRF_ub, IRRF_at0, IRRF_at0_lb, IRRF_at0_ub, IRRF_adj,
              weights, IRRFat0_adj, weightsat0))
IRRF_df$country <- countries

rm(list=c("IRRF_at0", "IRRF_at0_lb", "IRRF_at0_ub", "IRRF_adj",
          "weights", "IRRFat0_adj", "weightsat0"))

# IRRF_df has all the important variables


# ------------------------------------- PLOT --------------------------------- #
# Create DF to generate plot
IRRF_finland <- data.frame(IRRF = rbind(irrf_byc_f0[[7]][1,1],irrf_byc_f1[[7]][1,1],irrf_byc_f2[[7]][1,1],irrf_byc_f3[[7]][1,1]),
                           IRRF_lb = rbind(irrf_byc_f0[[7]][1,2],irrf_byc_f1[[7]][1,2],irrf_byc_f2[[7]][1,2],irrf_byc_f3[[7]][1,2]),
                           IRRF_ub = rbind(irrf_byc_f0[[7]][1,3],irrf_byc_f1[[7]][1,3],irrf_byc_f2[[7]][1,3],irrf_byc_f3[[7]][1,3]),
                           h = c(0,1:3))
# plot of US' and Findlands IRRF
IRRF_us<- data.frame(IRRF = rbind(irrf_byc_f0[[28]][1,1],irrf_byc_f1[[28]][1,1],irrf_byc_f2[[28]][1,1],irrf_byc_f3[[28]][1,1]),
                           IRRF_lb = rbind(irrf_byc_f0[[28]][1,2],irrf_byc_f1[[28]][1,2],irrf_byc_f2[[28]][1,2],irrf_byc_f3[[28]][1,2]),
                           IRRF_ub = rbind(irrf_byc_f0[[28]][1,3],irrf_byc_f1[[28]][1,3],irrf_byc_f2[[28]][1,3],irrf_byc_f3[[28]][1,3]),
                           h = c(0,1:3))

# -----------------------------------------------------------------------------
# FIGURE 1 OF THE FINAL PAPER
# -----------------------------------------------------------------------------
# PLOT IRRF across h
# FINLAND
ggplot(IRRF_finland) +  
    # geom_line(aes(x = h, y = IRRF_lb)) + geom_line(aes(x = h, y = IRRF_ub)) + 
    ggtitle("IRRF Finland") +ylab("Bond yield") +
    geom_ribbon(aes(ymin=IRRF_lb, ymax=IRRF_ub,x=h),fill = "#003366", alpha=0.7) + 
    geom_line(aes(x = h,y=IRRF),color= "red",size= 1, alpha = 0.5) + 
    ylim(c(-10,70)) +theme_stata(scheme = "s1mono") 

ggsave(paste0(figures,"FIG_1_IRRF_Finland.png"))

# USA
ggplot(IRRF_us) +  
    #geom_line(aes(x = h, y = IRRF_lb)) + geom_line(aes(x = h, y = IRRF_ub)) + 
    ggtitle("IRRF USA")+ ylab("Bond yield") +
     geom_ribbon(aes(ymin=IRRF_lb, ymax=IRRF_ub,x=h), fill = "#003366", alpha=0.7)+ 
    geom_line(aes(x = h,y=IRRF),color= "red",size=1,alpha = 0.5) + 
    ylim(c(-30,10)) +theme_stata(scheme = "s1mono") 

ggsave(paste0(figures,"FIG_1_IRRF_USA.png"))

# ---------------------------------------------------------------------------- #
# -------------------------- PRRF -------------------------------------------- #
# ---------------------------------------------------------------------------- #

# we use the same DF, just make some changes
mpanel5 <- mpanel
mpanel5$rgdp <- mpanel$rgdp/1000
mpanel5$c_real <- mpanel$c_real*1000
mpanel5$lngdp <- log(mpanel5$rgdp)
mpanel5$lnc <- log(mpanel5$c_real)


mpanelsub7 <- list();mpanelsub4<-NA
for (i in 1:28) {
    # test; i=1
    mpanelsub5 <- mpanel5 %>% filter(country == countries[i])
    
    # fix NA problem
    df4 <- mpanelsub5 %>% filter(!is.na(mpanelsub5$rgdp)) %>% 
        select(country,year,quarter,rgdp)
    df5 <-  hpfilter(df4$rgdp,freq=1600,drift=TRUE)
    df6 <- data.frame(df4, cycle_hp = df5$cycle)
    df6$trend_rgdp <- df6$rgdp- df6$cycle_hp
    
    mpanelsub6 <- merge(mpanelsub5,df6,by = c("country","year","quarter","rgdp"),all = T,)
    
    df5 <- bkfilter(df4$rgdp,pl = 6,pu=32,nfix = 12,drift = T)
    df6 <- data.frame(df4, cycle_bk = df5$cycle)
    df6$trend_rgdpBK <- df6$rgdp - df6$cycle_bk
    mpanelsub7[[i]] <- merge(mpanelsub6,df6,by = c("country","year","quarter","rgdp"),all = T)
    mpanelsub4 <- rbind(mpanelsub4,mpanelsub3[[i]])
}

rm(list=c("mpanelsub5","mpanelsub6","mpanelsub7","df4","df5","df6"))

mpanelsub4$trend_rgdpBK <- ifelse(is.na(mpanelsub4$trend_rgdpBK), mpanelsub4$trend_rgdp, mpanelsub4$trend_rgdpBK)


if (filter=="HP"){mpanelsub4 <- mpanelsub4 %>% mutate(dlng = log((G)/trend_rgdp), dlnc = log((c_real)/trend_rgdp),
                                                         dlngdp = log((rgdp)/trend_rgdp))
} else {mpanelsub4 <-mpanelsub4 %>% mutate(dlng = log((G)/trend_rgdpBK), dlnc = log((c_real)/trend_rgdpBK),
                                           dlngdp = log((rgdp)/trend_rgdpBK))}

# Create function that determines the variables of the model
var.model <- function(dep.lags, pol.rate) {
    if (pol.rate == "no") {
        model.var <- shift(policy_rate,dep.lags) ~ dlng + L.(dlngdp,lags) + L.(dlng,lags) + L.(bond_yields,lags)  
    } else {
        model.var <- shift(policy_rate,dep.lags) ~ dlng + L.(dlngdp,lags) + L.(dlng,lags) + L.(policy_rate,lags) + L.(bond_yields,lags)
    }}    

prrf_byc_f0 <- list();prrf_byc_f1 <- list();prrf_byc_f2 <- list();prrf_byc_f3 <- list()
for (i in 1:28) {
    # i = 1
    # BOND YIELD F0   
    #summary(reg)
    data <- mpanelsub4 %>% filter(country == countries[i])
    reg <- lm_robust(var.model(0,inc.polrate),data = data, se_type = "HC1")
    reg2 <- c(reg$coefficients[2],reg$std.error[2])
    prrf_byc_f0[[i]] <- data.frame(estimates = reg2,low_bound = reg2[1] - reg2[2]*1.645,upp_bound = reg2[1] + reg2[2]*1.645) 
    # BOND YIELD F1   
    reg <- lm_robust(var.model(-1,inc.polrate),data = data, se_type = "HC1")
    reg2 <- c(reg$coefficients[2],reg$std.error[2])
    prrf_byc_f1[[i]] <- data.frame(estimates = reg2,low_bound = reg2[1] - reg2[2]*1.645,upp_bound = reg2[1] + reg2[2]*1.645) 
    # BOND YIELD F2
    reg <- lm_robust(var.model(-2,inc.polrate),data = data, se_type = "HC1")
    reg2 <- c(reg$coefficients[2],reg$std.error[2])
    prrf_byc_f2[[i]] <- data.frame(estimates = reg2,low_bound = reg2[1] - reg2[2]*1.645,upp_bound = reg2[1] + reg2[2]*1.645) 
    # BOND YIELD F3
    reg <- lm_robust(var.model(-3,inc.polrate),data = data, se_type = "HC1")
    reg2 <- c(reg$coefficients[2],reg$std.error[2])
    prrf_byc_f3[[i]] <- data.frame(estimates = reg2,low_bound = reg2[1] - reg2[2]*1.645,upp_bound = reg2[1] + reg2[2]*1.645) 
}

PRRF <- NA;PRRF_lb <- NA;PRRF_ub <- NA;PRRF_at0 <- NA;PRRF_at0_lb <- NA;PRRF_at0_ub <- NA
PRRF_adj <- NA;PRRFat0_adj <- NA;Pweights <- NA;Pweightsat0 <- NA;PRRF_f <- matrix(NA,28,4)
for (i in 1:28) {
    if (def_irrf == 4) {
        PRRF[i] <- lags.1*(prrf_byc_f0[[i]][1,1]+prrf_byc_f1[[i]][1,1]+prrf_byc_f2[[i]][1,1]+prrf_byc_f3[[i]][1,1])
        PRRF_lb[i] <- lags.1*(prrf_byc_f0[[i]][1,2]+prrf_byc_f1[[i]][1,2]+prrf_byc_f2[[i]][1,2]+prrf_byc_f3[[i]][1,2])
        PRRF_ub[i] <- lags.1*(prrf_byc_f0[[i]][1,3]+prrf_byc_f1[[i]][1,3]+prrf_byc_f2[[i]][1,3]+prrf_byc_f3[[i]][1,3])
    } else {
        PRRF[i] <- lags.1*prrf_byc_f0[[i]][1,1]
        PRRF_lb[i] <- lags.1*prrf_byc_f0[[i]][1,2]
        PRRF_ub[i] <- lags.1*prrf_byc_f0[[i]][1,3]
    }
    
    PRRF_at0[i] <- lags.1*prrf_byc_f0[[i]][1,1]
    PRRF_at0_lb[i] <- lags.1*prrf_byc_f0[[i]][1,2]
    PRRF_at0_ub[i] <- lags.1*prrf_byc_f0[[i]][1,3]
    
    PRRF_adj[i] <- PRRF[i]/(PRRF_ub[i] - PRRF_lb[i])
    Pweights[i] <- 1/(PRRF_ub[i] - PRRF_lb[i])
    PRRFat0_adj[i] <- PRRF_at0[i]/(PRRF_at0_ub[i] - PRRF_at0_lb[i])
    Pweightsat0[i] <- 1/(PRRF_at0_ub[i] - PRRF_at0_lb[i])
    PRRF_f[i,] <- cbind(prrf_byc_f0[[i]][1,1],prrf_byc_f1[[i]][1,1],prrf_byc_f2[[i]][1,1],prrf_byc_f3[[i]][1,1])
}

PRRF_df <- data.frame(cbind(PRRF, PRRF_lb, PRRF_ub, PRRF_at0, PRRF_at0_lb, PRRF_at0_ub, PRRF_adj,
                            Pweights, PRRFat0_adj, Pweightsat0))
PRRF_df$country <- countries

rm(list=c("PRRF_at0", "PRRF_at0_lb", "PRRF_at0_ub", "PRRF_adj",
          "Pweights", "PRRFat0_adj", "Pweightsat0"))


# ----------------------------- CRF ---------------------------------------- #

# Create function that determines the variables of the model
var.model <- function(dep.lags, pol.rate) {
    if (pol.rate == "no") {
        model.var <- shift(dlnc,dep.lags) ~ dlng + L.(dlnc,lags) + L.(dlngdp,lags) + L.(dlng,lags)
    } else {
        model.var <- shift(dlnc,dep.lags) ~ dlng + L.(dlnc,lags) + L.(policy_rate,lags) + L.(dlngdp,lags) + L.(dlng,lags)
    }}

# Generate CRF
crf_byc_f0 <- list();crf_byc_f1 <- list();crf_byc_f2 <- list();crf_byc_f3 <- list()
for (i in 1:28) {
    # BOND YIELD F0   
    data <- mpanelsub4 %>% filter(country == countries[i])
    reg <- lm_robust(var.model(0,inc.polrate),data = data,se_type = "HC1")
    reg2 <- c(reg$coefficients[2],reg$std.error[2])
    crf_byc_f0[[i]] <- data.frame(estimates = reg2,low_bound = reg2[1] - reg2[2]*1.645,upp_bound = reg2[1] + reg2[2]*1.645) 
    # BOND YIELD F1   
    reg <- lm_robust(var.model(-1,inc.polrate),data = data,se_type = "HC1")
    reg2 <- c(reg$coefficients[2],reg$std.error[2])
    crf_byc_f1[[i]] <- data.frame(estimates = reg2,low_bound = reg2[1] - reg2[2]*1.645,upp_bound = reg2[1] + reg2[2]*1.645) 
    # BOND YIELD F2
    reg <- lm_robust(var.model(-2,inc.polrate),data = data,se_type = "HC1")
    reg2 <- c(reg$coefficients[2],reg$std.error[2])
    crf_byc_f2[[i]] <- data.frame(estimates = reg2,low_bound = reg2[1] - reg2[2]*1.645,upp_bound = reg2[1] + reg2[2]*1.645) 
    # BOND YIELD F3
    reg <- lm_robust(var.model(-3,inc.polrate),data = data,se_type = "HC1")
    reg2 <- c(reg$coefficients[2],reg$std.error[2])
    crf_byc_f3[[i]] <- data.frame(estimates = reg2,low_bound = reg2[1] - reg2[2]*1.645,upp_bound = reg2[1] + reg2[2]*1.645) 
}



CRF <- NA;CRF_lb <- NA;CRF_ub <- NA;CRF_at0 <- NA;CRF_at0_lb <- NA;CRF_at0_ub <- NA
CRF_adj <- NA;CRFat0_adj <- NA;Cweights <- NA;Cweightsat0 <- NA;CRF_f <- matrix(NA,28,4)
for (i in 1:28) {
    if (def_irrf == 4) {
        CRF[i] <- lags.1*(crf_byc_f0[[i]][1,1]+crf_byc_f1[[i]][1,1]+crf_byc_f2[[i]][1,1]+crf_byc_f3[[i]][1,1])
        CRF_lb[i] <- lags.1*(crf_byc_f0[[i]][1,2]+crf_byc_f1[[i]][1,2]+crf_byc_f2[[i]][1,2]+crf_byc_f3[[i]][1,2])
        CRF_ub[i] <- lags.1*(crf_byc_f0[[i]][1,3]+crf_byc_f1[[i]][1,3]+crf_byc_f2[[i]][1,3]+crf_byc_f3[[i]][1,3])
    } else {
        CRF[i] <- lags.1*crf_byc_f0[[i]][1,1]
        CRF_lb[i] <- lags.1*crf_byc_f0[[i]][1,2]
        CRF_ub[i] <- lags.1*crf_byc_f0[[i]][1,3]
    }
    
    # question: why is it multiplied by 1/lags?
    CRF_at0[i] <- lags.1*crf_byc_f0[[i]][1,1]
    CRF_at0_lb[i] <- lags.1*crf_byc_f0[[i]][1,2]
    CRF_at0_ub[i] <- lags.1*crf_byc_f0[[i]][1,3]
    
    CRF_adj[i] <- CRF[i]/(CRF_ub[i] - CRF_lb[i])
    Cweights[i] <- 1/(CRF_ub[i] - CRF_lb[i])
    CRFat0_adj[i] <- CRF_at0[i]/(CRF_at0_ub[i] - CRF_at0_lb[i])
    Cweightsat0[i] <- 1/(CRF_at0_ub[i] - CRF_at0_lb[i])
    CRF_f[i,] <- cbind(crf_byc_f0[[i]][1,1],crf_byc_f1[[i]][1,1],crf_byc_f2[[i]][1,1],crf_byc_f3[[i]][1,1])
    
}

CRF_df <- data.frame(cbind(CRF, CRF_lb, CRF_ub, CRF_at0, CRF_at0_lb, CRF_at0_ub, CRF_adj,
                            Cweights, CRFat0_adj, Cweightsat0))
CRF_df$country <- countries

rm(list=c("CRF_at0", "CRF_at0_lb", "CRF_at0_ub", "CRF_adj",
          "Cweights", "CRFat0_adj", "Cweightsat0"))



# ---------------------------- EXPORT ---------------------------------------- #
# PRRF_df has all the important variables
PRRF_df[,-11]
IRRF_df[,-11]
CRF_df[,-11]
# export df 

writexl::write_xlsx(IRRF_df,paste0(results,"irrf_lp_country.xlsx"))
writexl::write_xlsx(PRRF_df,paste0(results,"prrf_lp_country.xlsx"))
writexl::write_xlsx(CRF_df,paste0(results,"crf_lp_country.xlsx"))


# ---------------------------------------------------------------------------- #
# --------------------------- IIRF ------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Now we follow Auerbach and Gorodnichenko's methodology to estimate the IRRF
# To do this we use semester variables instead of quarterly data
# generate semester variable
AGpaneldata$year_half <- ifelse(AGpaneldata$quarter==1|AGpaneldata$quarter==2,1,2)
# drop Estonia
AGpaneldata <- AGpaneldata[!(AGpaneldata$country == "Estonia"),]
countries <- unique(AGpaneldata$country)
# Create contry code
cty_code <- c("AUS","AUT","BEL","CAN","CZE","DNK","FIN","FRA","DEU","GRC","HUN",
              "ISL","IRL","ITA","JPN","KOR","NLD","NZL","NOR","POL","PRT","SVK",
              "SVN","ESP","SWE","CHE","GBR","USA")
AGpaneldata$location <- 0
# AGpaneldata$location_id <- 0
for (i in 1:28) {
    AGpaneldata$location <- ifelse(AGpaneldata$country==countries[i],cty_code[i],AGpaneldata$location)
#    AGpaneldata$location_id <- ifelse(AGpaneldata$country==countries[i],i,AGpaneldata$location_id)
}

# We use the mean between two quarters
mean_na <- function(x)base::mean(x,na.rm=T)
gshocks <- read_dta(paste0(rdata1,"g_shocks_oecd_forecast_AGAERPP13.dta"))
AGpanel <- summaryBy(bond_yields + policy_rate + gspread + rgdp ~ location + year + year_half,
                  keep.names = T,data = AGpaneldata,FUN = mean_na)

AGpanel <- merge(x= AGpanel,y = gshocks, by=c("location","year","year_half"))
AGpanel <- AGpanel[order(AGpanel$location,AGpanel$year),]

AGpanel <- AGpanel %>% group_by(location) %>% mutate(rgdp_g0 = 100*(rgdp-shift(rgdp))/shift(rgdp))

AGpanel <- AGpanel %>% group_by(location) %>% mutate(dropshock = mean(pure2007_FE_0,na.rm=T))
AGpanel <- AGpanel[!(is.nan(AGpanel$dropshock)),]
AGpanel <- AGpanel[!(AGpanel$location == "POL"),]
AGpanel <- AGpanel[order(AGpanel$location,AGpanel$year,AGpanel$year_half),]

ineq_AG <- readxl::read_excel(paste0(rdata1,"new_ineq_toAG.xlsx"))
ineq_AG <- ineq_AG %>% rename(inequality = p90_p10)

AGpanel <- merge(AGpanel, ineq_AG,by = "location")

# ---------------------------------------------------------------------------- #
# --------------------- Bond yields by country ------------------------------- #
# ---------------------------------------------------------------------------- #

# Create var that 
var.model <- function(dep.vari,dep.lags, pol.rate) {
    if (pol.rate == "no" | pol.rate== "BP") {
        model.var <- shift(dep.vari,dep.lags) ~ pure2007_FE_1  
    } else if (pol.rate == "yes") {
        if (dep.vari == "data1$bond_yields") {
            model.var <- shift(dep.vari,dep.lags) ~ pure2007_FE_1 + L.(policy_rate,aglags)
        } else {
            model.var <- shift(dep.vari,dep.lags) ~ pure2007_FE_1 + L.(bond_yields,aglags)
        }
    }}    


countries <- unique(AGpanel$country)
IIRF_bond1_<-NA;IIRF_bond1_lb_<-NA;IIRF_bond1_ub_<-NA;PRRF_bond1_<-NA;PRRF_bond1_lb_<-NA;PRRF_bond1_ub_<-NA
for (i in 1:length(unique(AGpanel$location))) {
# i = 1
#summary(breg)
    data1 <- AGpanel %>% filter(country == countries[i])
    # IRRF
    breg <- lm_robust(var.model(data1$bond_yields,0,inc.polrate),data=data1,se_type = "HC1")       
    breg2 <- c(breg$coefficients[2],breg$std.error[2])
    # IIRF F1
    breg <- lm_robust(var.model(data1$bond_yields,-1,inc.polrate),data=data1,se_type = "HC1")       
    breg3 <- c(breg$coefficients[2],breg$std.error[2])
    if (def_irrf == 4) {
        IIRF_bond1_[i] <- ag.lags1*(breg2[1] + breg3[1])
        IIRF_bond1_lb_[i] <- ag.lags1*((breg2[1] - breg2[2]*1.645) + (breg3[1] - breg3[2]*1.645))
        IIRF_bond1_ub_[i] <- ag.lags1*((breg2[1] + breg2[2]*1.645) + (breg3[1] + breg3[2]*1.645))
        
    } else {
        IIRF_bond1_[i] <- ag.lags1*(breg2[1])
        IIRF_bond1_lb_[i] <- ag.lags1*((breg2[1] - breg2[2]*1.645))
        IIRF_bond1_ub_[i] <- ag.lags1*((breg2[1] + breg2[2]*1.645))
    }
    # --------------------- PRRF ------------------------------------------ #
    # PRRF
    breg <- lm_robust(var.model(data1$policy_rate,0,inc.polrate),data=data1,se_type = "HC1")       
    breg2 <- c(breg$coefficients[2],breg$std.error[2])
    # PRRF F1
    breg <- lm_robust(var.model(data1$policy_rate,-1,inc.polrate),data=data1,se_type = "HC1")       
    breg3 <- c(breg$coefficients[2],breg$std.error[2])
    if (def_irrf == 4) {
        PRRF_bond1_[i] <- ag.lags1*(breg2[1] + breg3[1])
        PRRF_bond1_lb_[i] <- ag.lags1*((breg2[1] - breg2[2]*1.645) + (breg3[1] - breg3[2]*1.645))
        PRRF_bond1_ub_[i] <- ag.lags1*((breg2[1] + breg2[2]*1.645) + (breg3[1] + breg3[2]*1.645))
    } else { 
        PRRF_bond1_[i] <- ag.lags1*(breg2[1])
        PRRF_bond1_lb_[i] <- ag.lags1*((breg2[1] - breg2[2]*1.645))
        PRRF_bond1_ub_[i] <- ag.lags1*((breg2[1] + breg2[2]*1.645))
        }
}

IIRF_df1 <- data.frame(IIRF_bond1_,IIRF_bond1_lb_,IIRF_bond1_ub_)
PRRF_df1 <- data.frame(PRRF_bond1_,PRRF_bond1_lb_,PRRF_bond1_ub_)
IIRF_df1$country_code <- unique(AGpanel$location)
PRRF_df1$country_code <- unique(AGpanel$location)

write.csv(IIRF_df1,paste0(rdata,"iirf_oneatatime.csv"),row.names = F)
write.csv(PRRF_df1,paste0(rdata,"prrf_oneatatime.csv"),row.names = F)

# source("IRRF_Rstudio/codes/IRRF_inequality.R")
