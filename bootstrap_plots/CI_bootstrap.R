#this script will construct CI around the estimates of pin_pis and Kn_Ks using bootstrapping
#we will use sample() function in R which will sample for 'x' number of times from the already exisiting observations

#for bootstrapping we will write a function which will do the following:
# 1. take as input - the REGION-specific pi_n and pi_s
# 2. using this input - perform 500 bootstrapping cycles (replace = TRUE), and for every cycle obtain a mean of the pi_n and pi_s
# 3. get n bootstrapped pi_n/pi_s by taking the ratio of the pi_n and pi_s - where 'n' is the number of observations
# 4. return the 500 values

#importing relevant packages
library(ggplot2)
library(ggpubr)
library(cowplot)
library(tidyverse)

#### run this part only once to produce the bootstraps ####

# #writing the function which will perform the bootstrapping for pin_pis
# bootstrap_pinpis = function(pi_n, pi_s){
#   
#   #declaring an empty dataframe to store the bootstrapped values
#   bs_pinpis = data.frame()
#   
#   #using a for loop to perform bootstrapping 500 times
#   for(counter in seq(0,499, by = 1)){
#     
#     #obtaining a mean of every bootstrap for pin and pis
#     pin = mean(sample(x = pi_n, size = length(pi_n), replace = TRUE))
#     pis = mean(sample(x = pi_s, size = length(pi_s), replace = TRUE))
#     
#     #writing to a table
#     df = data.frame(pin = pin, pis = pis, pin_pis = pin/pis)
#     
#     #merging the means
#     bs_pinpis = rbind(bs_pinpis, df)
#     
#   }
#   
#   #returning
#   return(bs_pinpis[,"pin_pis"])
#   
# }
# 
# #writing the function which will perform the bootstrapping for Kn_Ks
# bootstrap_knks = function(k_n, k_s){
#   
#   #declaring an empty dataframe to store the bootstrapped values
#   bs_knks = data.frame()
#   
#   #using a for loop to perform bootstrapping 500 times
#   for(counter in seq(0,499, by = 1)){
#     
#     #obtaining a mean of every bootstrap for pin and pis
#     kn = mean(sample(x = k_n, size = length(k_n), replace = TRUE))
#     ks = mean(sample(x = k_s, size = length(k_s), replace = TRUE))
#     
#     #writing to a table
#     df = data.frame(kn = kn, ks = ks, kn_ks = kn/ks)
#     
#     #merging the means
#     bs_knks = rbind(bs_knks, df)
#     
#   }
#   
#   #returning
#   return(bs_knks[,"kn_ks"])
#   
# }
# 
# 
# ##### Section 1 - importing the outputs from Alag and post-processing #####
# 
# #### !! 1. H. sapiens !!
# 
# #### pop YRI ##### 
# 
# ### Part 1  - listing output files for both regions
# domain_stats_files_yri = list.files(path = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/TFBD/hsap/YRI/", pattern = "poly_div_dnabd_nondnabd_batch_*", full.names = TRUE)
# cds_stats_files_yri = list.files(path = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/wga_rerun/wga_hsap_YRI/", pattern = "poly_div_gene_stats_*", full.names = TRUE)
# 
# ### Part 2 - aggregating all batch files
# ##domain - YRI
# poly_div_stats_domain_yri = lapply(domain_stats_files_yri, function(x){read.delim(x, sep = "\t")})
# poly_div_stats_domain_yri = do.call(rbind, poly_div_stats_domain_yri)
# poly_div_stats_domain_yri = unique(poly_div_stats_domain_yri) ## -- taking unique due to the nature of DNABD naming file - scope for improvement
# 
# ##WGA - YRI
# poly_div_stats_yri = lapply(cds_stats_files_yri, function(x){read.delim(x, sep = "\t")})
# poly_div_stats_yri = do.call(rbind, poly_div_stats_yri)
# poly_div_stats_yri = unique(poly_div_stats_yri) ## -- taking unique due to the nature of DNABD naming file - scope for improvement
# 
# #### pop CEU ####
# 
# ### Part 1 - listing output files for both regions
# cds_stats_files_ceu = list.files(path = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/wga_rerun/wga_hsap_CEU/", pattern = "poly_div_gene_stats_*", full.names = TRUE)
# domain_stats_files_ceu = list.files(path = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/TFBD/hsap/CEU/", pattern = "poly_div_dnabd_nondnabd_batch_*", full.names = TRUE)
# 
# ### Part 2 - aggregating all batch files
# ## domain - CEU
# poly_div_stats_domain_ceu = lapply(domain_stats_files_ceu, function(x){read.delim(x, sep = "\t")})
# poly_div_stats_domain_ceu = do.call(rbind, poly_div_stats_domain_ceu)
# poly_div_stats_domain_ceu = unique(poly_div_stats_domain_ceu) ## -- taking unique due to the nature of DNABD naming file - scope for improvement
# 
# ##WGA - CEU
# poly_div_stats_ceu = lapply(cds_stats_files_ceu, function(x){read.delim(x, sep = "\t")})
# poly_div_stats_ceu = do.call(rbind, poly_div_stats_ceu)
# poly_div_stats_ceu = unique(poly_div_stats_ceu) ## -- taking unique due to the nature of DNABD naming file - scope for improvement
# 
# #### pop CHS ####
# 
# ### Part 1 - listing output files for both regions
# cds_stats_files_chs = list.files(path = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/wga_rerun/wga_hsap_CHS/", pattern = "poly_div_gene_stats_*", full.names = TRUE)
# domain_stats_files_chs = list.files(path = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/TFBD/hsap/CHS/", pattern = "poly_div_dnabd_nondnabd_batch_*", full.names = TRUE)
# 
# ### Part 2 - aggregating all batch files
# ## domain - CHS
# poly_div_stats_domain_chs = lapply(domain_stats_files_chs, function(x){read.delim(x, sep = "\t")})
# poly_div_stats_domain_chs = do.call(rbind, poly_div_stats_domain_chs)
# poly_div_stats_domain_chs = unique(poly_div_stats_domain_chs) ## -- taking unique due to the nature of DNABD naming file - scope for improvement
# 
# ##WGA - CHS
# poly_div_stats_chs = lapply(cds_stats_files_chs, function(x){read.delim(x, sep = "\t")})
# poly_div_stats_chs = do.call(rbind, poly_div_stats_chs)
# poly_div_stats_chs = unique(poly_div_stats_chs) ## -- taking unique due to the nature of DNABD naming file - scope for improveme
# 
# #cleaning up  
# rm(domain_stats_files_ceu, domain_stats_files_chs, domain_stats_files_yri, cds_stats_files_ceu, cds_stats_files_chs, cds_stats_files_yri)
# 
# #### !! 2. A. thaliana !! 
# 
# #### pop IB ##### 
# 
# ### Part 1  - listing output files for both regions
# domain_stats_files_ib = list.files(path = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/TFBD/aratha/IB/", pattern = "poly_div_dnabd_nondnabd_batch_*", full.names = TRUE)
# cds_stats_files_ib = list.files(path = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/wga_rerun/wga_aratha_IB/", pattern = "poly_div_gene_stats_*", full.names = TRUE)
# 
# ### Part 2 - aggregating all batch files
# ##domain - IB
# poly_div_stats_domain_ib = lapply(domain_stats_files_ib, function(x){read.delim(x, sep = "\t")})
# poly_div_stats_domain_ib = do.call(rbind, poly_div_stats_domain_ib)
# poly_div_stats_domain_ib = unique(poly_div_stats_domain_ib) ## -- taking unique due to the nature of DNABD naming file - scope for improvement
# 
# #cds - IB
# poly_div_stats_ib = lapply(cds_stats_files_ib, function(x){read.delim(x, sep = "\t")})
# poly_div_stats_ib = do.call(rbind, poly_div_stats_ib)
# poly_div_stats_ib = unique(poly_div_stats_ib) ## -- taking unique due to the nature of DNABD naming file - scope for improvement
# 
# #### pop NS ####
# 
# ### Part 1 - listing output files for both regions
# domain_stats_files_ns = list.files(path = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/TFBD/aratha/NS/", pattern = "poly_div_dnabd_nondnabd_batch_*", full.names = TRUE)
# cds_stats_files_ns = list.files(path = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/wga_rerun/wga_aratha_NS/", pattern = "poly_div_gene_stats_*", full.names = TRUE)
# 
# ### Part 2 - aggregating all batch files
# ## domain - NS
# poly_div_stats_domain_ns = lapply(domain_stats_files_ns, function(x){read.delim(x, sep = "\t")})
# poly_div_stats_domain_ns = do.call(rbind, poly_div_stats_domain_ns)
# poly_div_stats_domain_ns = unique(poly_div_stats_domain_ns) ## -- taking unique due to the nature of DNABD naming file - scope for improvement
# 
# #cds - NS
# poly_div_stats_ns = lapply(cds_stats_files_ns, function(x){read.delim(x, sep = "\t")})
# poly_div_stats_ns = do.call(rbind, poly_div_stats_ns)
# poly_div_stats_ns = unique(poly_div_stats_ns) ## -- taking unique due to the nature of DNABD naming file - scope for improvement
# 
# #### pop CA ####
# 
# ### Part 1 - listing output files for both regions
# domain_stats_files_ca = list.files(path = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/TFBD/aratha/CA/", pattern = "poly_div_dnabd_nondnabd_batch_*", full.names = TRUE)
# cds_stats_files_ca = list.files(path = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/wga_rerun/wga_aratha_CA/", pattern = "poly_div_gene_stats_*", full.names = TRUE)
# 
# ### Part 2 - aggregating all batch files
# ## domain - CA
# poly_div_stats_domain_ca = lapply(domain_stats_files_ca, function(x){read.delim(x, sep = "\t")})
# poly_div_stats_domain_ca = do.call(rbind, poly_div_stats_domain_ca)
# poly_div_stats_domain_ca = unique(poly_div_stats_domain_ca) ## -- taking unique due to the nature of DNABD naming file - scope for improvement
# 
# #cds - CA
# poly_div_stats_ca = lapply(cds_stats_files_ca, function(x){read.delim(x, sep = "\t")})
# poly_div_stats_ca = do.call(rbind, poly_div_stats_ca)
# poly_div_stats_ca = unique(poly_div_stats_ca) ## -- taking unique due to the nature of DNABD naming file - scope for improvement
# 
# #cleaning up
# rm(cds_stats_files_ib, cds_stats_files_ca, cds_stats_files_ns, domain_stats_files_ca, domain_stats_files_ns, domain_stats_files_ib)
# 
# 
# #### !! 3. D. melanogaster !!
# 
# #### pop ZAM ##### 
# 
# ### Part 1  - listing output files for both regions
# domain_stats_files_zam = list.files(path = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/TFBD/dmel/ZAM/", pattern = "poly_div_dnabd_nondnabd_batch_*", full.names = TRUE)
# cds_stats_files_zam = list.files(path = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/wga_rerun/wga_dmel_ZAM/", pattern = "poly_div_gene_*", full.names = TRUE)
# 
# ### Part 2 - aggregating all batch files
# ##domain - ZAM
# poly_div_stats_domain_zam = lapply(domain_stats_files_zam, function(x){read.delim(x, sep = "\t")})
# poly_div_stats_domain_zam = do.call(rbind, poly_div_stats_domain_zam)
# poly_div_stats_domain_zam = unique(poly_div_stats_domain_zam) ## -- taking unique due to the nature of DNABD naming file - scope for improvement
# 
# ##WGA - ZAM
# poly_div_stats_zam = lapply(cds_stats_files_zam, function(x){read.delim(x, sep = "\t")})
# poly_div_stats_zam = do.call(rbind, poly_div_stats_zam)
# poly_div_stats_zam = unique(poly_div_stats_zam) ## -- taking unique due to the nature of DNABD naming file - scope for improvement
# 
# #### pop SWE ####
# 
# ### Part 1 - listing output files for both regions
# domain_stats_files_swe = list.files(path = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/TFBD/dmel/SWE/", pattern = "poly_div_dnabd_nondnabd_batch_*", full.names = TRUE)
# cds_stats_files_swe = list.files(path = "/netscratch/dep_tsiantis/grp_laurent/joshi/TFBD_project_new/wga_rerun/wga_dmel_SWE/", pattern = "poly_div_gene_*", full.names = TRUE)
# 
# ### Part 2 - aggregating all batch files
# ## domain - NS
# poly_div_stats_domain_swe = lapply(domain_stats_files_swe, function(x){read.delim(x, sep = "\t")})
# poly_div_stats_domain_swe = do.call(rbind, poly_div_stats_domain_swe)
# poly_div_stats_domain_swe = unique(poly_div_stats_domain_swe) ## -- taking unique due to the nature of DNABD naming file - scope for improvement
# 
# ##WGA - SWE
# poly_div_stats_swe = lapply(cds_stats_files_swe, function(x){read.delim(x, sep = "\t")})
# poly_div_stats_swe = do.call(rbind, poly_div_stats_swe)
# poly_div_stats_swe = unique(poly_div_stats_swe) ## -- taking unique due to the nature of DNABD naming file - scope for improvement
# 
# #cleaning up
# rm(cds_stats_files_swe, cds_stats_files_zam, domain_stats_files_swe, domain_stats_files_zam)

# ##### Section 2 - calculating BS intervals and plotting #####
# 
# ###### pin_pis ######
# 
# ## 1. H. sapiens
# 
# ##pop - YRI
# yri_cds_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_yri$pi_nonsyn, pi_s = poly_div_stats_yri$pi_syn), region = "WGS", pop = "YRI", species = "H. sapiens")
# yri_dnabd_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_domain_yri$pi_nonsyn_dnabd, pi_s = poly_div_stats_domain_yri$pi_syn_dnabd), region = "DNABD", pop = "YRI", species = "H. sapiens")
# yri_nondnabd_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_domain_yri$pi_nonsyn_nondnabd, pi_s = poly_div_stats_domain_yri$pi_syn_nondnabd), region = "non-DNABD", pop = "YRI", species = "H. sapiens")
# 
# ##pop - CEU
# ceu_cds_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_ceu$pi_nonsyn, pi_s = poly_div_stats_ceu$pi_syn), region = "WGS", pop = "CEU", species = "H. sapiens")
# ceu_dnabd_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_domain_ceu$pi_nonsyn_dnabd, pi_s = poly_div_stats_domain_ceu$pi_syn_dnabd), region = "DNABD", pop = "CEU", species = "H. sapiens")
# ceu_nondnabd_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_domain_ceu$pi_nonsyn_nondnabd, pi_s = poly_div_stats_domain_ceu$pi_syn_nondnabd), region = "non-DNABD", pop = "CEU", species = "H. sapiens")
# 
# ##pop - CHS
# chs_cds_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_chs$pi_nonsyn, pi_s = poly_div_stats_chs$pi_syn), region = "WGS", pop = "CHS", species = "H. sapiens")
# chs_dnabd_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_domain_chs$pi_nonsyn_dnabd, pi_s = poly_div_stats_domain_chs$pi_syn_dnabd), region = "DNABD", pop = "CHS", species = "H. sapiens")
# chs_nondnabd_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_domain_chs$pi_nonsyn_nondnabd, pi_s = poly_div_stats_domain_chs$pi_syn_nondnabd), region = "non-DNABD", pop = "CHS", species = "H. sapiens")
# 
# 
# ## 2. A. thaliana
# 
# ##pop - IB
# ib_cds_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_ib$pi_nonsyn, pi_s = poly_div_stats_ib$pi_syn), region = "WGS", pop = "IB", species = "A. thaliana")
# ib_dnabd_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_domain_ib$pi_nonsyn_dnabd, pi_s = poly_div_stats_domain_ib$pi_syn_dnabd), region = "DNABD", pop = "IB", species = "A. thaliana")
# ib_nondnabd_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_domain_ib$pi_nonsyn_nondnabd, pi_s = poly_div_stats_domain_ib$pi_syn_nondnabd), region = "non-DNABD", pop = "IB", species = "A. thaliana")
# 
# ##pop - NS
# ns_cds_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_ns$pi_nonsyn, pi_s = poly_div_stats_ns$pi_syn), region = "WGS", pop = "NS", species = "A. thaliana")
# ns_dnabd_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_domain_ns$pi_nonsyn_dnabd, pi_s = poly_div_stats_domain_ns$pi_syn_dnabd), region = "DNABD", pop = "NS", species = "A. thaliana")
# ns_nondnabd_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_domain_ns$pi_nonsyn_nondnabd, pi_s = poly_div_stats_domain_ns$pi_syn_nondnabd), region = "non-DNABD", pop = "NS", species = "A. thaliana")
# 
# ##pop - CA
# ca_cds_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_ca$pi_nonsyn, pi_s = poly_div_stats_ca$pi_syn), region = "WGS", pop = "CA", species = "A. thaliana")
# ca_dnabd_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_domain_ca$pi_nonsyn_dnabd, pi_s = poly_div_stats_domain_ca$pi_syn_dnabd), region = "DNABD", pop = "CA", species = "A. thaliana")
# ca_nondnabd_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_domain_ca$pi_nonsyn_nondnabd, pi_s = poly_div_stats_domain_ca$pi_syn_nondnabd), region = "non-DNABD", pop = "CA", species = "A. thaliana")
# 
# 
# ## 3. D. melanogaster
# 
# ##pop - ZAM
# zam_cds_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_zam$pi_nonsyn, pi_s = poly_div_stats_zam$pi_syn), region = "WGS", pop = "ZAM", species = "D. melanogaster")
# zam_dnabd_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_domain_zam$pi_nonsyn_dnabd, pi_s = poly_div_stats_domain_zam$pi_syn_dnabd), region = "DNABD", pop = "ZAM", species = "D. melanogaster")
# zam_nondnabd_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_domain_zam$pi_nonsyn_nondnabd, pi_s = poly_div_stats_domain_zam$pi_syn_nondnabd), region = "non-DNABD", pop = "ZAM", species = "D. melanogaster")
# 
# ##pop - SWE
# swe_cds_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_swe$pi_nonsyn, pi_s = poly_div_stats_swe$pi_syn), region = "WGS", pop = "SWE", species = "D. melanogaster")
# swe_dnabd_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_domain_swe$pi_nonsyn_dnabd, pi_s = poly_div_stats_domain_swe$pi_syn_dnabd), region = "DNABD", pop = "SWE", species = "D. melanogaster")
# swe_nondnabd_pinpis = data.frame(pin_pis = bootstrap_pinpis(pi_n = poly_div_stats_domain_swe$pi_nonsyn_nondnabd, pi_s = poly_div_stats_domain_swe$pi_syn_nondnabd), region = "non-DNABD", pop = "SWE", species = "D. melanogaster")
# 
# 
# ### !! merging all stats together
# pin_pis_stats = rbind(yri_cds_pinpis, yri_dnabd_pinpis, yri_nondnabd_pinpis, ceu_cds_pinpis, ceu_dnabd_pinpis, ceu_nondnabd_pinpis, chs_cds_pinpis, chs_dnabd_pinpis, chs_nondnabd_pinpis,
#                       ib_cds_pinpis, ib_dnabd_pinpis, ib_nondnabd_pinpis, ns_cds_pinpis, ns_dnabd_pinpis, ns_nondnabd_pinpis, ca_cds_pinpis, ca_dnabd_pinpis, ca_nondnabd_pinpis, 
#                       zam_cds_pinpis, zam_dnabd_pinpis, zam_nondnabd_pinpis, swe_cds_pinpis, swe_dnabd_pinpis, swe_nondnabd_pinpis)
# 
# 
# #preparing to plot
# hsap = ggplot(data = pin_pis_stats[pin_pis_stats$species == "H. sapiens",], aes(x = pop, y = pin_pis, col = region, fill = region)) +
#   geom_boxplot(width = 0.2, show.legend = TRUE, alpha = 0.3, position = position_dodge(width = 0.6)) + ggtitle(expression(~italic("H. sapiens"))) + scale_y_continuous(breaks = seq(0,0.4,by=0.05), limits = c(0,0.4))+
#   xlab("Population") + ylab(bquote(~"\u03C0"[n]~"/"~"\u03C0"[s]))+ scale_colour_manual("Regions", values = c("DNABD" = "#1a1aff", "non-DNABD" = "#ffaa00", "WGS" = "#00cc66"))+
#   scale_fill_manual("Regions", values = c("DNABD" = "#1a1aff", "non-DNABD" = "#ffaa00", "WGS" = "#00cc66"))+
#   theme_bw()+theme(plot.title = element_text(hjust = 0.5), text = element_text(family = "times", size = 12))
# 
# aratha = ggplot(data = pin_pis_stats[pin_pis_stats$species == "A. thaliana",], aes(x = pop, y = pin_pis, col = region, fill = region)) +
#   geom_boxplot(width = 0.2, show.legend = TRUE, alpha = 0.3, position = position_dodge(width = 0.6)) + ggtitle(expression(~italic("A. thaliana"))) + scale_y_continuous(breaks = seq(0,0.4,by=0.05), limits = c(0,0.4))+
#   xlab("Population") + ylab(bquote(~"\u03C0"[n]~"/"~"\u03C0"[s]))+ scale_colour_manual("Regions", values = c("DNABD" = "#1a1aff", "non-DNABD" = "#ffaa00", "WGS" = "#00cc66"))+
#   scale_fill_manual("Regions", values = c("DNABD" = "#1a1aff", "non-DNABD" = "#ffaa00", "WGS" = "#00cc66"))+
#   theme_bw()+theme(plot.title = element_text(hjust = 0.5), text = element_text(family = "times", size = 12))
# 
# dmel = ggplot(data = pin_pis_stats[pin_pis_stats$species == "D. melanogaster",], aes(x = pop, y = pin_pis, col = region, fill = region)) +
#   geom_boxplot(width = 0.2, show.legend = TRUE, alpha = 0.3, position = position_dodge(width = 0.6)) + ggtitle(expression(~italic("D. melanogaster"))) + scale_y_continuous(breaks = seq(0,0.4,by=0.05), limits = c(0,0.4))+
#   xlab("Population") + ylab(bquote(~"\u03C0"[n]~"/"~"\u03C0"[s]))+ scale_colour_manual("Regions", values = c("DNABD" = "#1a1aff", "non-DNABD" = "#ffaa00", "WGS" = "#00cc66"))+
#   scale_fill_manual("Regions", values = c("DNABD" = "#1a1aff", "non-DNABD" = "#ffaa00", "WGS" = "#00cc66"))+
#   theme_bw()+theme(plot.title = element_text(hjust = 0.5), text = element_text(family = "times", size = 12))
# 
# 
# #arranging the plots and saving
# ggarrange(hsap, aratha, dmel, nrow = 1, common.legend = TRUE, legend = "bottom")
# ggsave(filename = "bootstrap_pinpis_allspecies.png", device = "png", dpi = 300, width = 8, height = 8)
# 
# 
# ###
# 
# 
# ###### Kn_Ks ######
# 
# ## 1. H. sapiens
# yri_cds_knks = data.frame(kn_ks = bootstrap_knks(k_n = poly_div_stats_yri$div_nonsyn, k_s = poly_div_stats_yri$div_syn), region = "WGS", pop = "YRI", species = "H. sapiens")
# yri_dnabd_knks = data.frame(kn_ks = bootstrap_knks(k_n = poly_div_stats_domain_yri$dnabd_div_nonsyn, k_s = poly_div_stats_domain_yri$dnabd_div_syn), region = "DNABD", pop = "YRI", species = "H. sapiens")
# yri_nondnabd_knks = data.frame(kn_ks = bootstrap_knks(k_n = poly_div_stats_domain_yri$nondnabd_div_nonsyn, k_s = poly_div_stats_domain_yri$nondnabd_div_syn), region = "non-DNABD", pop = "YRI", species = "H. sapiens")
# 
# ## 2. A. thaliana
# 
# ##pop - IB
# ib_cds_knks = data.frame(kn_ks = bootstrap_knks(k_n = poly_div_stats_ib$div_nonsyn, k_s = poly_div_stats_ib$div_syn), region = "WGS", pop = "IB", species = "A. thaliana")
# ib_dnabd_knks = data.frame(kn_ks = bootstrap_knks(k_n = poly_div_stats_domain_ib$dnabd_div_nonsyn, k_s = poly_div_stats_domain_ib$dnabd_div_syn), region = "DNABD", pop = "IB", species = "A. thaliana")
# ib_nondnabd_knks = data.frame(kn_ks = bootstrap_knks(k_n = poly_div_stats_domain_ib$nondnabd_div_nonsyn, k_s = poly_div_stats_domain_ib$nondnabd_div_syn), region = "non-DNABD", pop = "IB", species = "A. thaliana")
# 
# ## 3. D. melanogaster
# 
# ##pop - ZAM
# zam_cds_knks = data.frame(kn_ks = bootstrap_knks(k_n = poly_div_stats_zam$div_nonsyn, k_s = poly_div_stats_zam$div_syn), region = "WGS", pop = "ZAM", species = "D. melanogaster")
# zam_dnabd_knks = data.frame(kn_ks = bootstrap_knks(k_n = poly_div_stats_domain_zam$dnabd_div_nonsyn, k_s = poly_div_stats_domain_zam$dnabd_div_syn), region = "DNABD", pop = "ZAM", species = "D. melanogaster")
# zam_nondnabd_knks = data.frame(kn_ks = bootstrap_knks(k_n = poly_div_stats_domain_zam$nondnabd_div_nonsyn, k_s = poly_div_stats_domain_zam$nondnabd_div_syn), region = "non-DNABD", pop = "ZAM", species = "D. melanogaster")
# 
# 
# ##!! merging all stats together
# kn_ks_stats = rbind(yri_cds_knks, yri_dnabd_knks, yri_nondnabd_knks, ib_cds_knks, ib_dnabd_knks, ib_nondnabd_knks, zam_cds_knks, zam_dnabd_knks, zam_nondnabd_knks)
# 
# #preparing to plot
# hsap_knks = ggplot(data = kn_ks_stats[kn_ks_stats$species == "H. sapiens",], mapping = aes(x = pop, y = kn_ks, col = region, fill = region)) +
#   geom_boxplot(width = 0.2, show.legend = TRUE, alpha = 0.3, position = position_dodge(width = 0.6)) + ggtitle(expression(~italic("H. sapiens"))) + scale_y_continuous(breaks = seq(0,0.3,by=0.05), limits = c(0,0.3))+
#   ylab(bquote(~"K"[n]~"/"~"K"[s]))+ scale_colour_manual("Regions", values = c("DNABD" = "#1a1aff", "non-DNABD" = "#ffaa00", "WGS" = "#00cc66"))+
#   scale_fill_manual("Regions", values = c("DNABD" = "#1a1aff", "non-DNABD" = "#ffaa00", "WGS" = "#00cc66"))+
#   theme_bw()+theme(plot.title = element_text(hjust = 0.5), text = element_text(family = "times", size = 12), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
# 
# aratha_knks = ggplot(data = kn_ks_stats[kn_ks_stats$species == "A. thaliana",], mapping = aes(x = pop, y = kn_ks, col = region, fill = region)) +
#   geom_boxplot(width = 0.2, show.legend = TRUE, alpha = 0.3, position = position_dodge(width = 0.6)) + ggtitle(expression(~italic("A. thaliana"))) + scale_y_continuous(breaks = seq(0,0.3,by=0.05), limits = c(0,0.3))+
#   ylab(bquote(~"K"[n]~"/"~"K"[s]))+ scale_colour_manual("Regions", values = c("DNABD" = "#1a1aff", "non-DNABD" = "#ffaa00", "WGS" = "#00cc66"))+
#   scale_fill_manual("Regions",values = c("DNABD" = "#1a1aff", "non-DNABD" = "#ffaa00", "WGS" = "#00cc66"))+
#   theme_bw()+theme(plot.title = element_text(hjust = 0.5), text = element_text(family = "times", size = 12), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
# 
# dmel_knks = ggplot(data = kn_ks_stats[kn_ks_stats$species == "D. melanogaster",], mapping = aes(x = pop, y = kn_ks, col = region, fill = region)) +
#   geom_boxplot(width = 0.2, show.legend = TRUE, alpha = 0.3, position = position_dodge(width = 0.6)) + ggtitle(expression(~italic("D. melanogaster"))) + scale_y_continuous(breaks = seq(0,0.3,by=0.05), limits = c(0,0.3))+
#   ylab(bquote(~"K"[n]~"/"~"K"[s]))+ scale_colour_manual("Regions", values = c("DNABD" = "#1a1aff", "non-DNABD" = "#ffaa00", "WGS" = "#00cc66"))+
#   scale_fill_manual("Regions",values = c("DNABD" = "#1a1aff", "non-DNABD" = "#ffaa00", "WGS" = "#00cc66"))+
#   theme_bw()+theme(plot.title = element_text(hjust = 0.5), text = element_text(family = "times", size = 12), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
# 
# ggarrange(hsap_knks, aratha_knks, dmel_knks, nrow = 1, common.legend = TRUE, legend = "bottom")
# ggsave(filename = "bootstrap_knks_allspecies.png", device = "png", dpi = 300, width = 8, height = 8)

#### loading directly the bootstraps if exisiting ####

#loading the pre-written files
pinpis_stats = read.delim(file = "pinpis_stats_allspecies_500bs.txt", sep = "\t", header = TRUE)
knks_stats = read.delim(file = "knks_stats_allspecies_500bs.txt", sep = "\t", header = TRUE)

#!for the TFBD manuscript - skipping the CHS and CA populations from humans and thaliana respectively
pinpis_stats = pinpis_stats[!pinpis_stats$pop %in% c("CHS", "CA"),]

#setting region names - pinpis
pinpis_stats$region[pinpis_stats$region == "DNABD"] = "BD"
pinpis_stats$region[pinpis_stats$region == "non-DNABD"] = "non-BD"
pinpis_stats$region[pinpis_stats$region == "WGS"] = "All genes"

#setting region names - knks
knks_stats$region[knks_stats$region == "DNABD"] = "BD"
knks_stats$region[knks_stats$region == "non-DNABD"] = "non-BD"
knks_stats$region[knks_stats$region == "WGS"] = "All genes"

#### making plots with the pre-written data ####

###pinpis plots
#hsap
hsap = ggplot(data = pinpis_stats[pinpis_stats$species == "H. sapiens",], mapping = aes(x = pop, y = pin_pis, col = region, fill = region))+
  geom_boxplot(width = 0.2, show.legend = TRUE, alpha = 0.3, position = position_dodge(width = 0.6))+
  xlab("Population")+ylab(bquote(~"\u03C0"[n]~"/"~"\u03C0"[s]))+scale_colour_manual("Region", values = c("BD" = "#1a1aff", "non-BD" = "#ffaa00", "All genes" = "#00cc66"))+
  ggtitle(expression(~italic("H. sapiens"))) +scale_y_continuous(breaks = seq(0, 0.5, by = 0.1), limits = c(0,0.5))+ 
  scale_fill_manual("Region", values = c("BD" = "#1a1aff", "non-BD" = "#ffaa00", "All genes" = "#00cc66"))+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5), text = element_text(family = "times", size = 24), 
                        legend.text = element_text(family = "times", size = 20), legend.title = element_text(family = "times", size = 24), 
                        axis.title.x = element_text(margin = margin(t = 20)))

#aratha
aratha = ggplot(data = pinpis_stats[pinpis_stats$species == "A. thaliana",], mapping = aes(x = pop, y = pin_pis, col = region, fill = region))+
  geom_boxplot(width = 0.2, show.legend = TRUE, alpha = 0.3, position = position_dodge(width = 0.6))+
  xlab("Population")+ylab(bquote(~"\u03C0"[n]~"/"~"\u03C0"[s]))+scale_colour_manual("Region", values = c("BD" = "#1a1aff", "non-BD" = "#ffaa00", "All genes" = "#00cc66"))+
  ggtitle(expression(~italic("A. thaliana"))) +scale_y_continuous(breaks = seq(0, 0.5, by = 0.1), limits = c(0,0.5))+ 
  scale_fill_manual("Region", values = c("BD" = "#1a1aff", "non-BD" = "#ffaa00", "All genes" = "#00cc66"))+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5), text = element_text(family = "times", size = 24), 
                        legend.text = element_text(family = "times", size = 20), legend.title = element_text(family = "times", size = 24), 
                        axis.title.x = element_text(margin = margin(t = 20)))

#dmel
dmel = ggplot(data = pinpis_stats[pinpis_stats$species == "D. melanogaster",], mapping = aes(x = pop, y = pin_pis, col = region, fill = region))+
  geom_boxplot(width = 0.2, show.legend = TRUE, alpha = 0.3, position = position_dodge(width = 0.6))+
  xlab("Population")+ylab(bquote(~"\u03C0"[n]~"/"~"\u03C0"[s]))+scale_colour_manual("Region", values = c("BD" = "#1a1aff", "non-BD" = "#ffaa00", "All genes" = "#00cc66"))+
  ggtitle(expression(~italic("D. melanogaster"))) +scale_y_continuous(breaks = seq(0, 0.5, by = 0.1), limits = c(0,0.5))+ 
  scale_fill_manual("Region", values = c("BD" = "#1a1aff", "non-BD" = "#ffaa00", "All genes" = "#00cc66"))+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5), text = element_text(family = "times", size = 24), 
                        legend.text = element_text(family = "times", size = 20), legend.title = element_text(family = "times", size = 24), 
                        axis.title.x = element_text(margin = margin(t = 20)))

#arranging the plots
ggarrange(hsap, aratha, dmel, nrow = 1, common.legend = TRUE, legend = "right")

#saving the file
ggsave("pinpis_bootstraps_coding.tiff",device = "tiff", dpi = 300, width = 15, height = 12)

###knks plots
#hsap
hsap = ggplot(data = knks_stats[knks_stats$species == "H. sapiens",], mapping = aes(x = pop, y = kn_ks, col = region, fill = region))+
  geom_boxplot(width = 0.2, show.legend = TRUE, alpha = 0.3, position = position_dodge(width = 0.6))+
  xlab("Population")+ylab(bquote(~"K"[n]~"/"~"K"[s]))+scale_colour_manual("Region", values = c("BD" = "#1a1aff", "non-BD" = "#ffaa00", "All genes" = "#00cc66"))+
  ggtitle(expression(~italic("H. sap - P. tro"))) +scale_y_continuous(breaks = seq(0, 0.5, by = 0.1), limits = c(0,0.5))+ 
  scale_fill_manual("Region", values = c("BD" = "#1a1aff", "non-BD" = "#ffaa00", "All genes" = "#00cc66"))+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5), text = element_text(family = "times", size = 24), 
                        legend.text = element_text(family = "times", size = 20), legend.title = element_text(family = "times", size = 24), 
                        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

#aratha
aratha = ggplot(data = knks_stats[knks_stats$species == "A. thaliana",], mapping = aes(x = pop, y = kn_ks, col = region, fill = region))+
  geom_boxplot(width = 0.2, show.legend = TRUE, alpha = 0.3, position = position_dodge(width = 0.6))+
  xlab("Population")+ylab(bquote(~"K"[n]~"/"~"K"[s]))+scale_colour_manual("Region", values = c("BD" = "#1a1aff", "non-BD" = "#ffaa00", "All genes" = "#00cc66"))+
  ggtitle(expression(~italic("A. tha - A. lyr"))) +scale_y_continuous(breaks = seq(0, 0.5, by = 0.1), limits = c(0,0.5))+ 
  scale_fill_manual("Region", values = c("BD" = "#1a1aff", "non-BD" = "#ffaa00", "All genes" = "#00cc66"))+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5), text = element_text(family = "times", size = 24), 
                        legend.text = element_text(family = "times", size = 20), legend.title = element_text(family = "times", size = 24), 
                        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

#dmel
dmel =ggplot(data = knks_stats[knks_stats$species == "D. melanogaster",], mapping = aes(x = pop, y = kn_ks, col = region, fill = region))+
  geom_boxplot(width = 0.2, show.legend = TRUE, alpha = 0.3, position = position_dodge(width = 0.6))+
  xlab("Population")+ylab(bquote(~"K"[n]~"/"~"K"[s]))+scale_colour_manual("Region", values = c("BD" = "#1a1aff", "non-BD" = "#ffaa00", "All genes" = "#00cc66"))+
  ggtitle(expression(~italic("D. mel - D. sim"))) +scale_y_continuous(breaks = seq(0, 0.5, by = 0.1), limits = c(0,0.5))+ 
  scale_fill_manual("Region", values = c("BD" = "#1a1aff", "non-BD" = "#ffaa00", "All genes" = "#00cc66"))+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5), text = element_text(family = "times", size = 24), 
                        legend.text = element_text(family = "times", size = 20), legend.title = element_text(family = "times", size = 24), 
                        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

#arranging the plots
ggarrange(hsap, aratha, dmel, nrow = 1, common.legend = TRUE, legend = "right")

#saving the file
ggsave("knks_bootstraps_coding.tiff",device = "tiff", dpi = 300, width = 15, height = 12)



###### Section 3 - calculating the significance between sets using Welch's t-test ######

# #writing the files containing the pinpis and KnKs stats
# write.table(x = pin_pis_stats, file = "pinpis_stats_allspecies_500bs.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
# write.table(x = kn_ks_stats, file = "knks_stats_allspecies_500bs.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  

#reading the outputs from previous sections to perform significance tests 
knks_stats = read.delim(file = "knks_stats_allspecies_500bs.txt", header = TRUE, sep = "\t")
pin_pis_stats = read.delim(file = "pinpis_stats_allspecies_500bs.txt", header = TRUE, sep = "\t")

## here we will use Welch's t-test - which compares the differences in means of tow distributions, in contrast to student's t-test, Welch's t-test allows
## for differences in sample sizes and variances between two distributions
## the significance tests will always be DNABD vs non-DNABD and WGS for a species/population of species
## we will use the built-in t.test function in R - note - to perform standard t test use var.equal = TRUE and to perform Welch's t-test use var.equal = FALSE


#### pinpis significance tests

## !! 1. H. sapiens
hsap_pinpis = pin_pis_stats[pin_pis_stats$species == "H. sapiens",]

## pop-YRI
t.test(hsap_pinpis$pin_pis[hsap_pinpis$region == "DNABD" & hsap_pinpis$pop == "YRI"], 
       hsap_pinpis$pin_pis[hsap_pinpis$region == "non-DNABD" & hsap_pinpis$pop == "YRI"], var.equal = FALSE)
t.test(hsap_pinpis$pin_pis[hsap_pinpis$region == "DNABD" & hsap_pinpis$pop == "YRI"], 
       hsap_pinpis$pin_pis[hsap_pinpis$region == "WGS" & hsap_pinpis$pop == "YRI"], var.equal = FALSE)

## pop-CEU
t.test(hsap_pinpis$pin_pis[hsap_pinpis$region == "DNABD" & hsap_pinpis$pop == "CEU"], 
       hsap_pinpis$pin_pis[hsap_pinpis$region == "non-DNABD" & hsap_pinpis$pop == "CEU"], var.equal = FALSE)
t.test(hsap_pinpis$pin_pis[hsap_pinpis$region == "DNABD" & hsap_pinpis$pop == "CEU"], 
       hsap_pinpis$pin_pis[hsap_pinpis$region == "WGS" & hsap_pinpis$pop == "CEU"], var.equal = FALSE)

## pop-CHS
t.test(hsap_pinpis$pin_pis[hsap_pinpis$region == "DNABD" & hsap_pinpis$pop == "CHS"], 
       hsap_pinpis$pin_pis[hsap_pinpis$region == "non-DNABD" & hsap_pinpis$pop == "CHS"], var.equal = FALSE)
t.test(hsap_pinpis$pin_pis[hsap_pinpis$region == "DNABD" & hsap_pinpis$pop == "CHS"], 
       hsap_pinpis$pin_pis[hsap_pinpis$region == "WGS" & hsap_pinpis$pop == "CHS"], var.equal = FALSE)


## !! 2. A. thaliana
aratha_pinpis = pin_pis_stats[pin_pis_stats$species == "A. thaliana",]

## pop-IB
t.test(aratha_pinpis$pin_pis[aratha_pinpis$region == "DNABD" & aratha_pinpis$pop == "IB"], 
       aratha_pinpis$pin_pis[aratha_pinpis$region == "non-DNABD" & aratha_pinpis$pop == "IB"], var.equal = FALSE)
t.test(aratha_pinpis$pin_pis[aratha_pinpis$region == "DNABD" & aratha_pinpis$pop == "IB"], 
       aratha_pinpis$pin_pis[aratha_pinpis$region == "WGS" & aratha_pinpis$pop == "IB"], var.equal = FALSE)

## pop-NS
t.test(aratha_pinpis$pin_pis[aratha_pinpis$region == "DNABD" & aratha_pinpis$pop == "NS"], 
       aratha_pinpis$pin_pis[aratha_pinpis$region == "non-DNABD" & aratha_pinpis$pop == "NS"], var.equal = FALSE)
t.test(aratha_pinpis$pin_pis[aratha_pinpis$region == "DNABD" & aratha_pinpis$pop == "NS"], 
       aratha_pinpis$pin_pis[aratha_pinpis$region == "WGS" & aratha_pinpis$pop == "NS"], var.equal = FALSE)

## pop-CA
t.test(aratha_pinpis$pin_pis[aratha_pinpis$region == "DNABD" & aratha_pinpis$pop == "CA"], 
       aratha_pinpis$pin_pis[aratha_pinpis$region == "non-DNABD" & aratha_pinpis$pop == "CA"], var.equal = FALSE)
t.test(aratha_pinpis$pin_pis[aratha_pinpis$region == "DNABD" & aratha_pinpis$pop == "CA"], 
       aratha_pinpis$pin_pis[aratha_pinpis$region == "WGS" & aratha_pinpis$pop == "CA"], var.equal = FALSE)


## !! 3. D. melanogaster
dmel_pinpis = pin_pis_stats[pin_pis_stats$species == "D. melanogaster",]

## pop-ZAM
t.test(dmel_pinpis$pin_pis[dmel_pinpis$region == "DNABD" & dmel_pinpis$pop == "ZAM"], 
       dmel_pinpis$pin_pis[dmel_pinpis$region == "non-DNABD" & dmel_pinpis$pop == "ZAM"], var.equal = FALSE)
t.test(dmel_pinpis$pin_pis[dmel_pinpis$region == "DNABD" & dmel_pinpis$pop == "ZAM"], 
       dmel_pinpis$pin_pis[dmel_pinpis$region == "WGS" & dmel_pinpis$pop == "ZAM"], var.equal = FALSE)

## pop-SWE
t.test(dmel_pinpis$pin_pis[dmel_pinpis$region == "DNABD" & dmel_pinpis$pop == "SWE"], 
       dmel_pinpis$pin_pis[dmel_pinpis$region == "non-DNABD" & dmel_pinpis$pop == "SWE"], var.equal = FALSE)
t.test(dmel_pinpis$pin_pis[dmel_pinpis$region == "DNABD" & dmel_pinpis$pop == "SWE"], 
       dmel_pinpis$pin_pis[dmel_pinpis$region == "WGS" & dmel_pinpis$pop == "SWE"], var.equal = FALSE)


##### KnKs significance tests

## !! 1. H. sapiens
t.test(kn_ks_stats$kn_ks[kn_ks_stats$region == "DNABD" & kn_ks_stats$species == "H. sapiens"],
       kn_ks_stats$kn_ks[kn_ks_stats$region == "non-DNABD" & kn_ks_stats$species == "H. sapiens"], var.equal = FALSE)
t.test(kn_ks_stats$kn_ks[kn_ks_stats$region == "DNABD" & kn_ks_stats$species == "H. sapiens"],
       kn_ks_stats$kn_ks[kn_ks_stats$region == "WGS" & kn_ks_stats$species == "H. sapiens"], var.equal = FALSE)


## !! 2. A. thaliana
t.test(kn_ks_stats$kn_ks[kn_ks_stats$region == "DNABD" & kn_ks_stats$species == "A. thaliana"],
       kn_ks_stats$kn_ks[kn_ks_stats$region == "non-DNABD" & kn_ks_stats$species == "A. thaliana"], var.equal = FALSE)
t.test(kn_ks_stats$kn_ks[kn_ks_stats$region == "DNABD" & kn_ks_stats$species == "A. thaliana"],
       kn_ks_stats$kn_ks[kn_ks_stats$region == "WGS" & kn_ks_stats$species == "A. thaliana"], var.equal = FALSE)

## !! 3. D. melanogaster
t.test(kn_ks_stats$kn_ks[kn_ks_stats$region == "DNABD" & kn_ks_stats$species == "D. melanogaster"],
       kn_ks_stats$kn_ks[kn_ks_stats$region == "non-DNABD" & kn_ks_stats$species == "D. melanogaster"], var.equal = FALSE)
t.test(kn_ks_stats$kn_ks[kn_ks_stats$region == "DNABD" & kn_ks_stats$species == "D. melanogaster"],
       kn_ks_stats$kn_ks[kn_ks_stats$region == "WGS" & kn_ks_stats$species == "D. melanogaster"], var.equal = FALSE)

###### Section 4 - Plotting only the means and CI around pin_pis and Kn_Ks estimates ######

#### first for pin_pis
#the imported bootstrapping values would first need to be converted into a comptaible format - i.e extracting only means and CI
#for this we will use tidyverse
#the function used to calculate CI is mean_cl_normal --> comes from Hmisc
pin_pis_CI = pin_pis_stats %>% group_by(region, pop, species) %>% summarise(ci = list(mean_cl_normal(pin_pis) %>% rename(mean = y, LI = ymin, UI = ymax))) %>% unnest
 
#preparing to plot
hsap = ggplot(data = pin_pis_CI[pin_pis_CI$species == "H. sapiens",], aes(x = pop, y = mean, col = region, fill = region, ymin = LI, ymax = UI)) +
  geom_errorbar(width = 0.8, show.legend = TRUE, position = position_dodge(width = 0.6)) + geom_point(size = 2, alpha = 0.6, position = position_dodge(width = 0.6))+
  ggtitle(expression(~italic("H. sapiens"))) + scale_y_continuous(breaks = seq(0,0.4,by=0.05), limits = c(0,0.4))+
  xlab("Population") + ylab(bquote(~"\u03C0"[n]~"/"~"\u03C0"[s]))+ scale_colour_manual("Regions", values = c("DNABD" = "#1a1aff", "non-DNABD" = "#ffaa00", "WGS" = "#00cc66"))+
  scale_fill_manual("Regions", values = c("DNABD" = "#1a1aff", "non-DNABD" = "#ffaa00", "WGS" = "#00cc66"))+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5), text = element_text(family = "times", size = 12))

aratha = ggplot(data = pin_pis_CI[pin_pis_CI$species == "A. thaliana",], aes(x = pop, y = mean, col = region, fill = region, ymin = LI, ymax = UI)) +
  geom_errorbar(width = 0.8, show.legend = TRUE, position = position_dodge(width = 0.6)) + geom_point(size = 2, alpha = 0.6, position = position_dodge(width = 0.6))+
  ggtitle(expression(~italic("A. thaliana"))) + scale_y_continuous(breaks = seq(0,0.4,by=0.05), limits = c(0,0.4))+
  xlab("Population") + ylab(bquote(~"\u03C0"[n]~"/"~"\u03C0"[s]))+ scale_colour_manual("Regions", values = c("DNABD" = "#1a1aff", "non-DNABD" = "#ffaa00", "WGS" = "#00cc66"))+
  scale_fill_manual("Regions", values = c("DNABD" = "#1a1aff", "non-DNABD" = "#ffaa00", "WGS" = "#00cc66"))+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5), text = element_text(family = "times", size = 12))

dmel = ggplot(data = pin_pis_CI[pin_pis_CI$species == "D. melanogaster",], aes(x = pop, y = mean, col = region, fill = region, ymin = LI, ymax = UI)) +
  geom_errorbar(width = 0.8, show.legend = TRUE, position = position_dodge(width = 0.6)) + geom_point(size = 2, alpha = 0.6, position = position_dodge(width = 0.6))+
  ggtitle(expression(~italic("D. melanogaster"))) + scale_y_continuous(breaks = seq(0,0.4,by=0.05), limits = c(0,0.4))+
  xlab("Population") + ylab(bquote(~"\u03C0"[n]~"/"~"\u03C0"[s]))+ scale_colour_manual("Regions", values = c("DNABD" = "#1a1aff", "non-DNABD" = "#ffaa00", "WGS" = "#00cc66"))+
  scale_fill_manual("Regions", values = c("DNABD" = "#1a1aff", "non-DNABD" = "#ffaa00", "WGS" = "#00cc66"))+
  theme_bw()+theme(plot.title = element_text(hjust = 0.5), text = element_text(family = "times", size = 12))


#arranging the plots and saving
ggarrange(hsap, aratha, dmel, nrow = 1, common.legend = TRUE, legend = "bottom")
ggsave(filename = "bootstrap_pinpis_allspecies_onlyCI.png", device = "png", dpi = 300, width = 8, height = 8) 


