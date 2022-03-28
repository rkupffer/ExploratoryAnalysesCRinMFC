# run TIRT models

library(MplusAutomation)
library(TirtAutomation)

runModels(target = "3_Mplus")


# Model fit

bt <- readModels("3_Mplus/03_TIRT_BT.out")
bt$warnings
bt$errors
bt$summaries
bt <- correctRmsea(n=3,p=20,bt$summaries["ChiSqM_Value"], bt$summaries["ChiSqM_DF"], bt$summaries["Observations"])
bt

bi <- readModels("3_Mplus/03_TIRT_BI.out")
bi$warnings
bi$errors
bi$summaries
bi <- correctRmsea(n=3,p=20,bi$summaries["ChiSqM_Value"], bi$summaries["ChiSqM_DF"], bi$summaries["Observations"])
bi

ho <- readModels("3_Mplus/03_TIRT_HO.out")
ho$warnings
ho$errors
ho$summaries
ho <- correctRmsea(n=3,p=20,ho$summaries["ChiSqM_Value"], ho$summaries["ChiSqM_DF"], ho$summaries["Observations"])
ho

ip <- readModels("3_Mplus/03_TIRT_IP.out")
ip$warnings
ip$errors
ip$summaries
ip <- correctRmsea(n=3,p=20,ip$summaries["ChiSqM_Value"], ip$summaries["ChiSqM_DF"], ip$summaries["Observations"])
ip

os <- readModels("3_Mplus/03_TIRT_OS.out")
os$warnings
os$errors
os$summaries
os <- correctRmsea(n=3,p=30,os$summaries["ChiSqM_Value"], os$summaries["ChiSqM_DF"], os$summaries["Observations"])
os
