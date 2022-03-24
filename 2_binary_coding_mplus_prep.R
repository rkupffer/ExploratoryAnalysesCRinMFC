####________________________ binary coding and preparation of data for TIRT-analyses in Mplus __________________________####
# data from OSF
#data <- readOSF("HERE", header = TRUE, stringsAsFactors = FALSE) 
df <- read.csv("../OSF/DataExploratoryAnalysesCRinMFC.csv")[-1]

# which item belongs to which trait and how is it keyed
all_items <- read.csv2("../OSF/2_All_Items_Coding.csv", sep = ",")

# packages
library(CRinMFC)
library(TirtAutomation)
library(multiplex)

# ------------------------------------------------- bft ---> bt ------------------------------------------------- 
# recode triplets
df_bt <- reTri(subset(df, select = MT01_01:MT20_03), fbname = "bt")           

# merge with df
df <- cbind(df, df_bt)       

# select all binary coded variables and ID
df_bt <- subset(df, select = c(ID, b.bt.1:b.bt.60))

# recode missings for mplus
df_bt[,2:61] <- lapply(df_bt[,2:61], function(x) car::recode(x, recodes = "NA=-99; 1=1; 0=0"))

write.dat(df_bt, "3_Mplus")

# matrix (items x traits) with design loadings
m_dload <- designLoadings(all_items, quest = "BT", no.traits = 5)
colnames(m_dload) <- c("TNeu", "TExt", "TOpe", "TAgr", "TCon")
rownames(m_dload) <- ifelse((1:nrow(m_dload))<10,paste0("T","0",1:nrow(m_dload)), paste0("T",1:nrow(m_dload)))

tirt_bt <- tirtMplusSyntax(design.load = m_dload, 
                           names.pairs = NULL,
                           item.short = "T",
                           id.var = "ID",
                           file.data = "df_bt.dat", 
                           title = "03_TIRT_BFT", 
                           out.command = "sampstat standardized;",
                           fscores.file = "fscores_bt.dat",
                           missings.as = "-99")

# save as input-file
cat(paste(tirt_bt, collapse="\n\n"), file="3_Mplus/03_TIRT_BT.inp")

# ------------------------------------------------- bfi ---> bi ------------------------------------------------- 
# recode triplets
df_bi <- reTri(subset(df, select = BF01_01:BF20_03), fbname = "bi")     

# merge with df
df <- cbind(df, df_bi)  

# select all binary coded variables and ID
df_bi <- subset(df, select = c(ID, b.bi.1:b.bi.60))

# recode missings for mplus
df_bi[,2:61] <- lapply(df_bi[,2:61], function(x) car::recode(x, recodes = "NA=-99; 1=1; 0=0"))

write.dat(df_bi, "3_Mplus")

# matrix (items x traits) with design loadings
m_dload <- designLoadings(all_items, quest = "BI", no.traits = 5)
colnames(m_dload) <- c("INeu", "IExt", "IOpe", "IAgr", "ICon")
rownames(m_dload) <- ifelse((1:nrow(m_dload))<10,paste0("I","0",1:nrow(m_dload)), paste0("I",1:nrow(m_dload)))
head(m_dload)

tirt_bi <- tirtMplusSyntax(design.load = m_dload, 
                             names.pairs = NULL,
                             item.short = "I",
                             id.var = "ID",
                             file.data = "df_bi.dat", 
                             title = "03_TIRT_BFI", 
                             out.command = "sampstat standardized;",
                             fscores.file = "fscores_bi.dat",
                             missings.as = "-99")

# save as input-file
cat(paste(tirt_bi, collapse="\n\n"), file="3_Mplus/03_TIRT_BI.inp")



# ------------------------------------------------- SD3 ---> sd ------------------------------------------------- 
# recode triplets
df_sd <- reTri(subset(df, select = DT01_01:DT09_03), fbname = "sd")     

# merge with df
df <- cbind(df, df_sd)     

# select all binary coded variables and ID
df_sd <- subset(df, select = c(ID, b.sd.1:b.sd.27))

# recode missings for mplus
df_sd[,2:28] <- lapply(df_sd[,2:28], function(x) car::recode(x, recodes = "NA=-99; 1=1; 0=0"))

write.dat(df_sd, "3_Mplus")

# matrix (items x traits) with design loadings
m_dload <- designLoadings(all_items, quest = "SD", no.traits = 3)
colnames(m_dload) <- c("SMac", "SNar", "SPsy")
rownames(m_dload) <- ifelse((1:nrow(m_dload))<10,paste0("S","0",1:nrow(m_dload)), paste0("S",1:nrow(m_dload)))
head(m_dload)

tirt_sd <- tirtMplusSyntax(design.load = m_dload, 
                             names.pairs = NULL,
                             item.short = "S",
                             id.var = "ID",
                             file.data = "df_sd.dat", 
                             title = "03_TIRT_SD3", 
                             out.command = "sampstat standardized;",
                             fscores.file = "fscores_sd.dat",
                             missings.as = "-99")

# save as input-file
cat(paste(tirt_sd, collapse="\n\n"), file="3_Mplus/03_TIRT_SD.inp")


# ------------------------------------------------- HEXACO ---> ho ------------------------------------------------- 
# recode triplets
df_ho <- reTri(subset(df, select = HM01_01:HM20_03), fbname = "ho")     

# merge with df
df <- cbind(df, df_ho)       

# select all binary coded variables and ID
df_ho <- subset(df, select = c(ID, b.ho.1:b.ho.60))

# recode missings for mplus
df_ho[,2:61] <- lapply(df_ho[,2:61], function(x) car::recode(x, recodes = "NA=-99; 1=1; 0=0"))

write.dat(df_ho, "3_Mplus")


# matrix (items x traits) with design loadings
m_dload <- designLoadings(all_items, quest = "HO", no.traits = 6)
colnames(m_dload) <- c("HHon", "HEmo", "HeXt", "HAng", "HCon", "HOpe")
rownames(m_dload) <- ifelse((1:nrow(m_dload))<10,paste0("H","0",1:nrow(m_dload)), paste0("H",1:nrow(m_dload)))
head(m_dload)

tirt_ho <- tirtMplusSyntax(design.load = m_dload, 
                             names.pairs = NULL,
                             item.short = "H",
                             id.var = "ID",
                             file.data = "df_ho.dat", 
                             title = "03_TIRT_HEXACO", 
                             out.command = "sampstat standardized;",
                             fscores.file = "fscores_ho.dat",
                             missings.as = "-99")

# save as input-file
cat(paste(tirt_ho, collapse="\n\n"), file="3_Mplus/03_TIRT_HO.inp")

# ------------------------------------------------- IPIP ---> ip ------------------------------------------------- 
# recode triplets
df_ip <- reTri(subset(df, select = IP01_01:IP20_03), fbname = "ip")     

# merge with df
df <- cbind(df, df_ip)       

# select all binary coded variables and ID
df_ip <- subset(df, select = c(ID, b.ip.1:b.ip.60))

# recode missings for mplus
df_ip[,2:61] <- lapply(df_ip[,2:61], function(x) car::recode(x, recodes = "NA=-99; 1=1; 0=0"))

write.dat(df_ip, "3_Mplus")

# matrix (items x traits) with design loadings
m_dload <- designLoadings(all_items, quest = "IP", no.traits = 6)
colnames(m_dload) <- c("ICon", "IDis", "IImp", "IJoy", "IRec", "ISel")
rownames(m_dload) <- ifelse((1:nrow(m_dload))<10,paste0("I","0",1:nrow(m_dload)), paste0("I",1:nrow(m_dload)))
head(m_dload)

tirt_ip <- tirtMplusSyntax(design.load = m_dload, 
                             names.pairs = NULL,
                             item.short = "I",
                             id.var = "ID",
                             file.data = "df_ip.dat", 
                             title = "03_TIRT_IPIP", 
                             out.command = "sampstat standardized;",
                             fscores.file = "fscores_ip.dat",
                             missings.as = "-99")

# save as input-file
cat(paste(tirt_ip, collapse="\n\n"), file="3_Mplus/03_TIRT_IP.inp")

# ------------------------------------------------- ORVIS ---> os ------------------------------------------------- 
# due to estimation problems, triplet 12 will be excluded
os2 <- subset(df, select = c(OR01_01:OR11_03,OR13_01:OR30_03)) 

df_os <- reTri(os2, "os") # recode triplets of ORVIS again but without triplet 12
# this way the items will be named throughout

df <- cbind(df, df_os)

df_os <- subset(df, select = c(ID, b.os.1:b.os.87))

df_os[,2:88] <- lapply(df_os[,2:88], function(x) car::recode(x, recodes = "NA=-99; 1=1; 0=0"))

write.dat(df_os, "3_Mplus")

#new all_items table without triplet 12
all_itemsOS_without12 <- all_items[which(all_items$questionnaire == "OS" & all_items$triplet != 12), ]

# matrix (items x traits) with design loadings
m_dload <- designLoadings(all_itemsOS_without12, quest = "OS", no.traits = 8)
colnames(m_dload) <- c("OLea", "OOrg", "OAlt", "OCre", "OAna", "OPro", "OAdv", "OEru")
rownames(m_dload) <- ifelse((1:nrow(m_dload))<10,paste0("O","0",1:nrow(m_dload)), paste0("O",1:nrow(m_dload)))
head(m_dload)

tirt_os <- tirtMplusSyntax(design.load = m_dload, 
                             names.pairs = NULL,
                             item.short = "O",
                             id.var = "ID",
                             file.data = "df_os.dat", 
                             title = "03_TIRT_ORVIS", 
                             out.command = "sampstat standardized;",
                             fscores.file = "fscores_os.dat",
                             missings.as = "-99")

# modifications in ORVIS TIRT model due to convergence problems
## fix correlations to specific values
## add some comments to the file
## fix eO63 (and not eO61)


# add a comment to the Mplus input file:
tirt_os$data <- "DATA: FILE IS 'df_os.dat';
    ! without triplet 12 and fixed correlations
    ! to values from Pozzebon et al. (2010)
    ! Table 2 (upper diagonale)"

# add fixed correlations 
tirt_os$model <- "MODEL: \nOLea BY \nO01O02*1 (L_O01) \nO01O03*1 (L_O01) \nO04O06*-1 (L_O06_n) 
\nO05O06*-1 (L_O06_n) \nO07O08*-1 (L_O08_n) \nO08O09*1 (L_O08) \nO10O11*1 (L_O10) 
\nO10O12*1 (L_O10) \nO13O14*-1 (L_O14_n) \nO14O15*1 (L_O14) \nO16O18*-1 (L_O18_n) 
\nO17O18*-1 (L_O18_n) \nO19O20*1 (L_O19) \nO19O21*1 (L_O19) \nO22O23*-1 (L_O23_n) 
\nO23O24*1 (L_O23) \nO25O27*-1 (L_O27_n) \nO26O27*-1 (L_O27_n) \nO28O29*1 (L_O28) 
\nO28O30*1 (L_O28) \nO31O32*-1 (L_O32_n) \nO32O33*1 (L_O32); \n\nOOrg BY 
\nO01O02*-1 (L_O02_n) \nO02O03*1 (L_O02) \nO04O05*1 (L_O04) \nO04O06*1 (L_O04) 
\nO07O09*-1 (L_O09_n) \nO08O09*-1 (L_O09_n) \nO34O35*-1 (L_O35_n) \nO35O36*1 (L_O35) 
\nO37O39*-1 (L_O39_n) \nO38O39*-1 (L_O39_n) \nO40O41*1 (L_O40) \nO40O42*1 (L_O40) 
\nO43O44*-1 (L_O44_n) \nO44O45*1 (L_O44) \nO46O48*-1 (L_O48_n) \nO47O48*-1 (L_O48_n) 
\nO49O51*-1 (L_O51_n) \nO50O51*-1 (L_O51_n) \nO52O53*-1 (L_O53_n) \nO53O54*1 (L_O53) 
\nO55O57*-1 (L_O57_n) \nO56O57*-1 (L_O57_n) \nO58O59*1 (L_O58) \nO58O60*1 (L_O58); 
\n\nOAlt BY \nO01O03*-1 (L_O03_n) \nO02O03*-1 (L_O03_n) \nO10O11*-1 (L_O11_n) 
\nO11O12*1 (L_O11) \nO13O15*-1 (L_O15_n) \nO14O15*-1 (L_O15_n) \nO34O35*1 (L_O34) 
\nO34O36*1 (L_O34) \nO37O38*1 (L_O37) \nO37O39*1 (L_O37) \nO61O62*-1 (L_O62_n) 
\nO62O63*1 (L_O62) \nO64O66*-1 (L_O66_n) \nO65O66*-1 (L_O66_n) \nO67O68*1 (L_O67) 
\nO67O69*1 (L_O67) \nO70O71*-1 (L_O71_n) \nO71O72*1 (L_O71) \nO73O75*-1 (L_O75_n) 
\nO74O75*-1 (L_O75_n) \nO76O77*1 (L_O76) \nO76O78*1 (L_O76) \nO79O80*-1 (L_O80_n) 
\nO80O81*1 (L_O80); \n\nOCre BY \nO04O05*-1 (L_O05_n) \nO05O06*1 (L_O05) 
\nO10O12*-1 (L_O12_n) \nO11O12*-1 (L_O12_n) \nO19O20*-1 (L_O20_n) \nO20O21*1 (L_O20) 
\nO40O41*-1 (L_O41_n) \nO41O42*1 (L_O41) \nO43O44*1 (L_O43) \nO43O45*1 (L_O43) 
\nO52O53*1 (L_O52) \nO52O54*1 (L_O52) \nO61O63*-1 (L_O63_n) \nO62O63*-1 (L_O63_n) 
\nO64O65*-1 (L_O65_n) \nO65O66*1 (L_O65) \nO67O69*-1 (L_O69_n) \nO68O69*-1 (L_O69_n) 
\nO70O71*1 (L_O70) \nO70O72*1 (L_O70) \nO79O80*1 (L_O79) \nO79O81*1 (L_O79) 
\nO82O83*-1 (L_O83_n) \nO83O84*1 (L_O83) \nO85O87*-1 (L_O87_n) \nO86O87*-1 (L_O87_n);
\n\nOAna BY \nO19O21*-1 (L_O21_n) \nO20O21*-1 (L_O21_n) \nO28O29*-1 (L_O29_n) 
\nO29O30*1 (L_O29) \nO31O32*1 (L_O31) \nO31O33*1 (L_O31) \nO34O36*-1 (L_O36_n) 
\nO35O36*-1 (L_O36_n) \nO40O42*-1 (L_O42_n) \nO41O42*-1 (L_O42_n) \nO46O47*1 (L_O46) 
\nO46O48*1 (L_O46) \nO49O50*-1 (L_O50_n) \nO50O51*1 (L_O50) \nO61O62*1 (L_O61) 
\nO61O63*1 (L_O61) \nO70O72*-1 (L_O72_n) \nO71O72*-1 (L_O72_n) \nO85O86*1 (L_O85)
\nO85O87*1 (L_O85); \n\nOPro BY \nO07O08*1 (L_O07) \nO07O09*1 (L_O07) 
\nO16O17*1 (L_O16) \nO16O18*1 (L_O16) \nO22O24*-1 (L_O24_n) \nO23O24*-1 (L_O24_n) 
\nO28O30*-1 (L_O30_n) \nO29O30*-1 (L_O30_n) \nO43O45*-1 (L_O45_n) 
\nO44O45*-1 (L_O45_n) \nO46O47*-1 (L_O47_n) \nO47O48*1 (L_O47) \nO55O56*1 (L_O55) 
\nO55O57*1 (L_O55) \nO64O65*1 (L_O64) \nO64O66*1 (L_O64) \nO73O74*-1 (L_O74_n) 
\nO74O75*1 (L_O74) \nO76O77*-1 (L_O77_n) \nO77O78*1 (L_O77); \n\nOAdv BY 
\nO22O23*1 (L_O22) \nO22O24*1 (L_O22) \nO25O26*-1 (L_O26_n) \nO26O27*1 (L_O26) 
\nO31O33*-1 (L_O33_n) \nO32O33*-1 (L_O33_n) \nO49O50*1 (L_O49) \nO49O51*1 (L_O49) 
\nO52O54*-1 (L_O54_n) \nO53O54*-1 (L_O54_n) \nO58O59*-1 (L_O59_n) \nO59O60*1 (L_O59) 
\nO67O68*-1 (L_O68_n) \nO68O69*1 (L_O68) \nO76O78*-1 (L_O78_n) \nO77O78*-1 (L_O78_n) 
\nO79O81*-1 (L_O81_n) \nO80O81*-1 (L_O81_n) \nO82O83*1 (L_O82) \nO82O84*1 (L_O82); 
\n\nOEru BY \nO13O14*1 (L_O13) \nO13O15*1 (L_O13) \nO16O17*-1 (L_O17_n) 
\nO17O18*1 (L_O17) \nO25O26*1 (L_O25) \nO25O27*1 (L_O25) \nO37O38*-1 (L_O38_n) 
\nO38O39*1 (L_O38) \nO55O56*-1 (L_O56_n) \nO56O57*1 (L_O56) \nO58O60*-1 (L_O60_n) 
\nO59O60*-1 (L_O60_n) \nO73O74*1 (L_O73) \nO73O75*1 (L_O73) \nO82O84*-1 (L_O84_n) 
\nO83O84*-1 (L_O84_n) \nO85O86*-1 (L_O86_n) \nO86O87*1 (L_O86); \n
\n! means for all traits are set to 0
\n[OLea@0]; [OOrg@0]; [OAlt@0]; [OCre@0]; [OAna@0]; [OPro@0]; [OAdv@0]; [OEru@0];\n
\n! variances for all traits are set to 1
\nOLea@1; OOrg@1; OAlt@1; OCre@1; OAna@1; OPro@1; OAdv@1; OEru@1;\n
\n ! fixed correlations between traits
  OLea WITH OOrg@0.44;
  OLea WITH OAlt@0.31;
  OLea WITH OCre@0.23;
  OLea WITH OAna@0.31;
  OLea WITH OPro@.1;
  OLea WITH OAdv@0.39;
  OLea WITH OEru@0.33;

  OOrg WITH OAlt@0.16;
  OOrg WITH OCre@-0.04;
  OOrg WITH OAna@0.41;
  OOrg WITH OPro@0.15;
  OOrg WITH OAdv@0.2;
  OOrg WITH OEru@0.13;

  OAlt WITH OCre@0.36;
  OAlt WITH OAna@0.14;
  OAlt WITH OPro@0.17;
  OAlt WITH OAdv@0.08;
  OAlt WITH OEru@0.44;

  OCre WITH OAna@0.15;
  OCre WITH OPro@0.22;
  OCre WITH OAdv@0.05;
  OCre WITH OEru@0.54;

  OAna WITH OPro@0.38;
  OAna WITH OAdv@0.34;
  OAna WITH OEru@0.28;

  OPro WITH OAdv@0.5;
  OPro WITH OEru@0.18;

  OAdv WITH OEru@0.04;
\n\n! declare uniquenesses and set their starting values
\nO01O02*2 (eO01eO02); \nO01O03*2 (eO01eO03); \nO02O03*2 (eO02eO03); 
\nO04O05*2 (eO04eO05); \nO04O06*2 (eO04eO06); \nO05O06*2 (eO05eO06); 
\nO07O08*2 (eO07eO08); \nO07O09*2 (eO07eO09); \nO08O09*2 (eO08eO09); 
\nO10O11*2 (eO10eO11); \nO10O12*2 (eO10eO12); \nO11O12*2 (eO11eO12); 
\nO13O14*2 (eO13eO14); \nO13O15*2 (eO13eO15); \nO14O15*2 (eO14eO15); 
\nO16O17*2 (eO16eO17); \nO16O18*2 (eO16eO18); \nO17O18*2 (eO17eO18); 
\nO19O20*2 (eO19eO20); \nO19O21*2 (eO19eO21); \nO20O21*2 (eO20eO21); 
\nO22O23*2 (eO22eO23); \nO22O24*2 (eO22eO24); \nO23O24*2 (eO23eO24); 
\nO25O26*2 (eO25eO26); \nO25O27*2 (eO25eO27); \nO26O27*2 (eO26eO27); 
\nO28O29*2 (eO28eO29); \nO28O30*2 (eO28eO30); \nO29O30*2 (eO29eO30); 
\nO31O32*2 (eO31eO32); \nO31O33*2 (eO31eO33); \nO32O33*2 (eO32eO33); 
\nO34O35*2 (eO34eO35); \nO34O36*2 (eO34eO36); \nO35O36*2 (eO35eO36); 
\nO37O38*2 (eO37eO38); \nO37O39*2 (eO37eO39); \nO38O39*2 (eO38eO39); 
\nO40O41*2 (eO40eO41); \nO40O42*2 (eO40eO42); \nO41O42*2 (eO41eO42); 
\nO43O44*2 (eO43eO44); \nO43O45*2 (eO43eO45); \nO44O45*2 (eO44eO45); 
\nO46O47*2 (eO46eO47); \nO46O48*2 (eO46eO48); \nO47O48*2 (eO47eO48); 
\nO49O50*2 (eO49eO50); \nO49O51*2 (eO49eO51); \nO50O51*2 (eO50eO51); 
\nO52O53*2 (eO52eO53); \nO52O54*2 (eO52eO54); \nO53O54*2 (eO53eO54); 
\nO55O56*2 (eO55eO56); \nO55O57*2 (eO55eO57); \nO56O57*2 (eO56eO57); 
\nO58O59*2 (eO58eO59); \nO58O60*2 (eO58eO60); \nO59O60*2 (eO59eO60); 
\nO61O62*2 (eO61eO62); \nO61O63*2 (eO61eO63); \nO62O63*2 (eO62eO63); 
\nO64O65*2 (eO64eO65); \nO64O66*2 (eO64eO66); \nO65O66*2 (eO65eO66); 
\nO67O68*2 (eO67eO68); \nO67O69*2 (eO67eO69); \nO68O69*2 (eO68eO69); 
\nO70O71*2 (eO70eO71); \nO70O72*2 (eO70eO72); \nO71O72*2 (eO71eO72); 
\nO73O74*2 (eO73eO74); \nO73O75*2 (eO73eO75); \nO74O75*2 (eO74eO75); 
\nO76O77*2 (eO76eO77); \nO76O78*2 (eO76eO78); \nO77O78*2 (eO77eO78); 
\nO79O80*2 (eO79eO80); \nO79O81*2 (eO79eO81); \nO80O81*2 (eO80eO81); 
\nO82O83*2 (eO82eO83); \nO82O84*2 (eO82eO84); \nO83O84*2 (eO83eO84); 
\nO85O86*2 (eO85eO86); \nO85O87*2 (eO85eO87); \nO86O87*2 (eO86eO87); 
\n\n! declare correlated uniquenesses and set their starting values
\nO01O02 WITH O01O03*1 (eO01); \nO01O02 WITH O02O03*-1 (eO02_n); 
\nO01O03 WITH O02O03*1 (eO03); \nO04O05 WITH O04O06*1 (eO04); 
\nO04O05 WITH O05O06*-1 (eO05_n); \nO04O06 WITH O05O06*1 (eO06); 
\nO07O08 WITH O07O09*1 (eO07); \nO07O08 WITH O08O09*-1 (eO08_n); 
\nO07O09 WITH O08O09*1 (eO09); \nO10O11 WITH O10O12*1 (eO10); 
\nO10O11 WITH O11O12*-1 (eO11_n); \nO10O12 WITH O11O12*1 (eO12); 
\nO13O14 WITH O13O15*1 (eO13); \nO13O14 WITH O14O15*-1 (eO14_n); 
\nO13O15 WITH O14O15*1 (eO15); \nO16O17 WITH O16O18*1 (eO16); 
\nO16O17 WITH O17O18*-1 (eO17_n); \nO16O18 WITH O17O18*1 (eO18); 
\nO19O20 WITH O19O21*1 (eO19); \nO19O20 WITH O20O21*-1 (eO20_n);
\nO19O21 WITH O20O21*1 (eO21); \nO22O23 WITH O22O24*1 (eO22);
\nO22O23 WITH O23O24*-1 (eO23_n); \nO22O24 WITH O23O24*1 (eO24);
\nO25O26 WITH O25O27*1 (eO25); \nO25O26 WITH O26O27*-1 (eO26_n);
\nO25O27 WITH O26O27*1 (eO27); \nO28O29 WITH O28O30*1 (eO28);
\nO28O29 WITH O29O30*-1 (eO29_n); \nO28O30 WITH O29O30*1 (eO30);
\nO31O32 WITH O31O33*1 (eO31); \nO31O32 WITH O32O33*-1 (eO32_n);
\nO31O33 WITH O32O33*1 (eO33); \nO34O35 WITH O34O36*1 (eO34); 
\nO34O35 WITH O35O36*-1 (eO35_n); \nO34O36 WITH O35O36*1 (eO36);
\nO37O38 WITH O37O39*1 (eO37); \nO37O38 WITH O38O39*-1 (eO38_n); 
\nO37O39 WITH O38O39*1 (eO39); \nO40O41 WITH O40O42*1 (eO40); 
\nO40O41 WITH O41O42*-1 (eO41_n); \nO40O42 WITH O41O42*1 (eO42);
\nO43O44 WITH O43O45*1 (eO43); \nO43O44 WITH O44O45*-1 (eO44_n);
\nO43O45 WITH O44O45*1 (eO45); \nO46O47 WITH O46O48*1 (eO46); 
\nO46O47 WITH O47O48*-1 (eO47_n); \nO46O48 WITH O47O48*1 (eO48); 
\nO49O50 WITH O49O51*1 (eO49); \nO49O50 WITH O50O51*-1 (eO50_n); 
\nO49O51 WITH O50O51*1 (eO51); \nO52O53 WITH O52O54*1 (eO52); 
\nO52O53 WITH O53O54*-1 (eO53_n); \nO52O54 WITH O53O54*1 (eO54); 
\nO55O56 WITH O55O57*1 (eO55); \nO55O56 WITH O56O57*-1 (eO56_n); 
\nO55O57 WITH O56O57*1 (eO57); \nO58O59 WITH O58O60*1 (eO58); 
\nO58O59 WITH O59O60*-1 (eO59_n); \nO58O60 WITH O59O60*1 (eO60); 
\nO61O62 WITH O61O63*1 (eO61); \nO61O62 WITH O62O63*-1 (eO62_n); 
\nO61O63 WITH O62O63*1 (eO63); \nO64O65 WITH O64O66*1 (eO64); 
\nO64O65 WITH O65O66*-1 (eO65_n); \nO64O66 WITH O65O66*1 (eO66); 
\nO67O68 WITH O67O69*1 (eO67); \nO67O68 WITH O68O69*-1 (eO68_n); 
\nO67O69 WITH O68O69*1 (eO69); \nO70O71 WITH O70O72*1 (eO70); 
\nO70O71 WITH O71O72*-1 (eO71_n); \nO70O72 WITH O71O72*1 (eO72); 
\nO73O74 WITH O73O75*1 (eO73); \nO73O74 WITH O74O75*-1 (eO74_n); 
\nO73O75 WITH O74O75*1 (eO75); \nO76O77 WITH O76O78*1 (eO76); 
\nO76O77 WITH O77O78*-1 (eO77_n); \nO76O78 WITH O77O78*1 (eO78); 
\nO79O80 WITH O79O81*1 (eO79); \nO79O80 WITH O80O81*-1 (eO80_n); 
\nO79O81 WITH O80O81*1 (eO81); \nO82O83 WITH O82O84*1 (eO82); 
\nO82O83 WITH O83O84*-1 (eO83_n); \nO82O84 WITH O83O84*1 (eO84); 
\nO85O86 WITH O85O87*1 (eO85); \nO85O86 WITH O86O87*-1 (eO86_n); 
\nO85O87 WITH O86O87*1 (eO87); \n"

# fix eO63 (and not eO61)
tirt_os$model.constraint <- "MODEL CONSTRAINT: 
\n\n! factor loadings relating to the same item are equal in absolute magnitude 
\nL_O02_n = -L_O02; \nL_O05_n = -L_O05; \nL_O08_n = -L_O08; \nL_O11_n = -L_O11; 
\nL_O14_n = -L_O14; \nL_O17_n = -L_O17; \nL_O20_n = -L_O20; \nL_O23_n = -L_O23; 
\nL_O26_n = -L_O26; \nL_O29_n = -L_O29; \nL_O32_n = -L_O32; \nL_O35_n = -L_O35; 
\nL_O38_n = -L_O38; \nL_O41_n = -L_O41; \nL_O44_n = -L_O44; \nL_O47_n = -L_O47; 
\nL_O50_n = -L_O50; \nL_O53_n = -L_O53; \nL_O56_n = -L_O56; \nL_O59_n = -L_O59; 
\nL_O62_n = -L_O62; \nL_O65_n = -L_O65; \nL_O68_n = -L_O68; \nL_O71_n = -L_O71; 
\nL_O74_n = -L_O74; \nL_O77_n = -L_O77; \nL_O80_n = -L_O80; \nL_O83_n = -L_O83; 
\nL_O86_n = -L_O86; \n
\n! pair's uniqueness is equal to the sum of 2 utility uniquenesses\n
eO01eO02 = eO01 - eO02_n; \neO01eO03 = eO01 + eO03; \neO02eO03 = -eO02_n + eO03; 
\neO04eO05 = eO04 - eO05_n; \neO04eO06 = eO04 + eO06; \neO05eO06 = -eO05_n + eO06; 
\neO07eO08 = eO07 - eO08_n; \neO07eO09 = eO07 + eO09; \neO08eO09 = -eO08_n + eO09; 
\neO10eO11 = eO10 - eO11_n; \neO10eO12 = eO10 + eO12; \neO11eO12 = -eO11_n + eO12; 
\neO13eO14 = eO13 - eO14_n; \neO13eO15 = eO13 + eO15; \neO14eO15 = -eO14_n + eO15; 
\neO16eO17 = eO16 - eO17_n; \neO16eO18 = eO16 + eO18; \neO17eO18 = -eO17_n + eO18; 
\neO19eO20 = eO19 - eO20_n; \neO19eO21 = eO19 + eO21; \neO20eO21 = -eO20_n + eO21; 
\neO22eO23 = eO22 - eO23_n; \neO22eO24 = eO22 + eO24; \neO23eO24 = -eO23_n + eO24; 
\neO25eO26 = eO25 - eO26_n; \neO25eO27 = eO25 + eO27; \neO26eO27 = -eO26_n + eO27; 
\neO28eO29 = eO28 - eO29_n; \neO28eO30 = eO28 + eO30; \neO29eO30 = -eO29_n + eO30; 
\neO31eO32 = eO31 - eO32_n; \neO31eO33 = eO31 + eO33; \neO32eO33 = -eO32_n + eO33; 
\neO34eO35 = eO34 - eO35_n; \neO34eO36 = eO34 + eO36; \neO35eO36 = -eO35_n + eO36; 
\neO37eO38 = eO37 - eO38_n; \neO37eO39 = eO37 + eO39; \neO38eO39 = -eO38_n + eO39; 
\neO40eO41 = eO40 - eO41_n; \neO40eO42 = eO40 + eO42; \neO41eO42 = -eO41_n + eO42; 
\neO43eO44 = eO43 - eO44_n; \neO43eO45 = eO43 + eO45; \neO44eO45 = -eO44_n + eO45; 
\neO46eO47 = eO46 - eO47_n; \neO46eO48 = eO46 + eO48; \neO47eO48 = -eO47_n + eO48; 
\neO49eO50 = eO49 - eO50_n; \neO49eO51 = eO49 + eO51; \neO50eO51 = -eO50_n + eO51; 
\neO52eO53 = eO52 - eO53_n; \neO52eO54 = eO52 + eO54; \neO53eO54 = -eO53_n + eO54; 
\neO55eO56 = eO55 - eO56_n; \neO55eO57 = eO55 + eO57; \neO56eO57 = -eO56_n + eO57; 
\neO58eO59 = eO58 - eO59_n; \neO58eO60 = eO58 + eO60; \neO59eO60 = -eO59_n + eO60; 
\neO61eO62 = eO61 - eO62_n; \neO61eO63 = eO61 + eO63; \neO62eO63 = -eO62_n + eO63; 
\neO64eO65 = eO64 - eO65_n; \neO64eO66 = eO64 + eO66; \neO65eO66 = -eO65_n + eO66; 
\neO67eO68 = eO67 - eO68_n; \neO67eO69 = eO67 + eO69; \neO68eO69 = -eO68_n + eO69; 
\neO70eO71 = eO70 - eO71_n; \neO70eO72 = eO70 + eO72; \neO71eO72 = -eO71_n + eO72; 
\neO73eO74 = eO73 - eO74_n; \neO73eO75 = eO73 + eO75; \neO74eO75 = -eO74_n + eO75; 
\neO76eO77 = eO76 - eO77_n; \neO76eO78 = eO76 + eO78; \neO77eO78 = -eO77_n + eO78; 
\neO79eO80 = eO79 - eO80_n; \neO79eO81 = eO79 + eO81; \neO80eO81 = -eO80_n + eO81; 
\neO82eO83 = eO82 - eO83_n; \neO82eO84 = eO82 + eO84; \neO83eO84 = -eO83_n + eO84; 
\neO85eO86 = eO85 - eO86_n; \neO85eO87 = eO85 + eO87; \neO86eO87 = -eO86_n + eO87; 
\n
  ! fix one uniquness per block for identification
  eO01=1;
  eO04=1;
  eO07=1;
  eO10=1;
  eO13=1;
  eO16=1;
  eO19=1;
  eO22=1;
  eO25=1;
  eO28=1;
  eO31=1;
  eO34=1;
  eO37=1;
  eO40=1;
  eO43=1;
  eO46=1;
  eO49=1;
  eO52=1;
  eO55=1;
  eO58=1;
  !eO61=1;
  eO63=1;
  eO64=1;
  eO67=1;
  eO70=1;
  eO73=1;
  eO76=1;
  eO79=1;
  eO82=1;
  eO85=1;"

# save as input-file
cat(paste(tirt_os, collapse="\n\n"), file="3_Mplus/03_TIRT_OS.inp")








