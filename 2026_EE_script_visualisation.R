
### INTRODUCTION #####

# Welcome. For the description of the project please visit: https://github.com/diekei/2026_EE_diekei_dicliptera_speciation
# Article is available at: 


## LIBRARY ####

library(tidyverse)
library(ggpubr)
library(ggplot2)
library(gghighlight)
library(ggnewscale)
library(ggside)
library(grid)
library(png)
library(scales)
library(Hmisc)
library(MASS)


## AESTHETICS ####

dic_m <- readPNG("icon/dic_m.png", native = FALSE)
dic_d <- readPNG("icon/dic_d.png", native = FALSE)
dic_mm <- readPNG("dic_mm_l.png", native = FALSE)
dic_dd <- readPNG("dic_dd_l.png", native = FALSE)
dic_md <- readPNG("dic_md_l.png", native = FALSE)
dic_dm <- readPNG("dic_dm_l.png", native = FALSE)


## HOST PREFERENCE ####

hp.pop <- c("BDG (6)", "PDL (17)", "TBR (14)", "PAT (42)", "RAN (48)")
hp.pref <- c("Mikania", "Mikania", "Mikania", "Dicliptera", "Dicliptera")
hp.prop <- c(1, 1, 0.857, -0.929, -0.833)
hp <- data.frame(hp.pop, hp.pref, hp.prop)

hp$hp.pop <- factor(hp$hp.pop, 
                    levels = c("BDG (6)", "PDL (17)", "TBR (14)", "PAT (42)", "RAN (48)"))

hp.labs <- c("BDG (6)" = "BDG [(6)]", "PDL (17)" = "PDL [(17)]", "TBR (14)" = "TBR [(14)]", 
             "PAT (42)" = "PAT [(42)]", "RAN (48)" = "RAN [(48)]")

ggplot(hp, aes(fill = hp.pref, y = hp.prop, x = hp.pop)) + 
  geom_bar(stat = "identity", width = 0.7) + 
  scale_fill_manual(values = c("#ABD1DC", "#BAC94A")) + 
  scale_color_manual(values = c("#ABD1DC", "#BAC94A")) + 
  scale_x_discrete(labels = parse(text = hp.labs)) +
  xlab("\nPopulation") + 
  ylab(expression(paste("Preference for ", italic("Mikania")))) +
  guides(fill = guide_legend(title = "pref")) +
  scale_y_continuous(
    limits = c(-1, 1.07),
    breaks = seq(-1, 1, 0.5),
    labels = function(x) gsub("\\.0$", "", x)
  ) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = "none",
        panel.grid = element_blank()) + 
  annotate("text", x = c(1,2,3,4,5), y = c(0.95, 0.95, 0.807, -0.869, -0.773), 
           label = c("a", "a", "a", "b", "b"), size = 2.5, color = "white") + 
  annotation_custom(rasterGrob(dic_m), xmin = "PAT (42)", xmax = "RAN (48)", ymin = 0.75, ymax = 1) + 
  annotation_custom(rasterGrob(dic_d), xmin = "BDG (6)", xmax = "PDL (17)", ymin = -1, ymax = -0.75)

ggsave(filename = "figures/diekei_dic_hp.png", width = 2.5, height = 5, device='png', dpi=1200)
ggsave(filename = "figures/diekei_dic_hp.pdf", width = 2.5, height = 5, device='pdf', dpi=1200)


## NON CHOICE MATING ####

#glm_ncm2.csv
ncm2 <- read.csv('data/01_2026_EE_vis_nonchoice_mating.csv')
ncm2
attach(ncm2)

ncm2$fem <- factor(ncm2$fem, levels = c("mrace", "drace"))
ncm2$mle <- factor(ncm2$mle, levels = c("mrace", "drace"))
ncm2$pair <- factor(ncm2$pair, levels = c("mm", "md", "dm", "dd"))

ncm.labs <- data.frame(fem = c("mrace", "mrace", "mrace", "mrace", 
                               "drace", "drace", "drace", "drace"), 
                       label = c("a", "a", "b", "c", "b", "bc", "ab", "ab"))

ncm.labs$fem <- factor(ncm.labs$fem, levels = c("mrace", "drace"))

fem.labs <- c("\u2640 M-race", "\u2640 D-race")
names(fem.labs) <- c("mrace", "drace")

mle.labs <- c(mrace = expression(italic(M) * "-race" ~ (20)), drace = expression(italic(D) * "-race" ~ (40)))
mle.labs <- c(
  mrace = bquote(atop(italic("\u2642 M-race"), "(20)")),
  drace = bquote(atop(italic("\u2642 D-race"), "(40)"))
)

ncm2 <- ncm2 %>%
  mutate(pair_exp = interaction(pair, exp, sep = "-"))

ncm2 |>
  dplyr::group_by(fem, mle, pair, pair_exp) |>
  dplyr::summarise(
    mean = mean(pro, na.rm = TRUE),
    sem = sd(pro, na.rm = TRUE) / sqrt(sum(!is.na(pro))),
    n = sum(!is.na(pro)),
    .groups = "drop"
  )

ncm.plot <- ggplot(ncm2, aes(x = mle, y = pro)) + 
  facet_wrap(~fem, labeller = labeller(fem = as_labeller(fem.labs))) + 
  stat_summary(aes(color = pair, shape = pair_exp, fill = pair, size = 2), 
               fun.data = "mean_se", fun.args = list(mult=1), 
               geom = "pointrange", size = 0.6, position = position_dodge(0.5)) + 
  scale_color_manual(values = c("#BAC94A", "#5D6E1E", "#3B5284", "#ABD1DC")) +
  scale_fill_manual(values = c("#BAC94A", "#5D6E1E", "#3B5284", "#ABD1DC")) +
  scale_shape_manual(values = c("mm-att" = 21, "mm-suc" = 1, 
                                "md-att" = 25, "md-suc" = 6, 
                                "dm-att" = 22, "dm-suc" = 0,
                                "dd-att" = 24, "dd-suc" = 2)) +
  xlab("\nMating pair") + 
  ylab("Proportion of trials (non-choice mating test)\n") + 
  scale_x_discrete(labels = (text = mle.labs)) + 
  scale_y_continuous(limits =  c(-0.35,1.25), breaks = seq(0,1,by=0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(face = "italic"),
        panel.grid = element_blank()) + 
  geom_text(x =c(0.9,1.125,1.9,2.125,0.88,1.125,1.875,2.125), 
            y = c(0.85,0.85,-0.05,-0.05,0.25,0.25,0.48,0.48), 
            aes(label = label), data = ncm.labs, size = 2.5, 
            color = c("#BAC94A","#BAC94A","#5D6E1E","#5D6E1E","#3B5284","#3B5284","#ABD1DC","#ABD1DC"))

ncm.plot

ncm.an <- function (grob, xmin = -Inf, xmax = Inf, 
                    ymin = -Inf, ymax = Inf, data){ layer(data = data, 
                                                          stat = StatIdentity, 
                                                          position = PositionIdentity, 
                                                          geom = ggplot2:::GeomCustomAnn, 
                                                          inherit.aes = TRUE, 
                                                          params = list(grob = grob, xmin = xmin, xmax = xmax, 
                                                                        ymin = ymin, ymax = ymax))}

ncm.an1 <- ncm.an(rasterGrob(dic_m, interpolate=TRUE), xmin=1, xmax=2, ymin=1, ymax=1.25, data = ncm2[1,])
ncm.an2 <- ncm.an(rasterGrob(dic_d, interpolate=TRUE), xmin=1, xmax=2, ymin=1, ymax=1.25, data = ncm2[9,])
ncm.an3 <- ncm.an(rasterGrob(dic_m, interpolate=TRUE), xmin=0.5, xmax=1.5, ymin=-0.35, ymax=-0.15, data = ncm2[1,])
ncm.an4 <- ncm.an(rasterGrob(dic_d, interpolate=TRUE), xmin=1.5, xmax=2.5, ymin=-0.35, ymax=-0.15, data = ncm2[1,])
ncm.an5 <- ncm.an(rasterGrob(dic_m, interpolate=TRUE), xmin=0.5, xmax=1.5, ymin=-0.35, ymax=-0.15, data = ncm2[9,])
ncm.an6 <- ncm.an(rasterGrob(dic_d, interpolate=TRUE), xmin=1.5, xmax=2.5, ymin=-0.35, ymax=-0.15, data = ncm2[9,])

ncm.plot + ncm.an1 + ncm.an2 + ncm.an3 + ncm.an4 + ncm.an5 + ncm.an6

png("figures/diekei_dic_ncm.png", 
    width = 4, height = 5, units = "in", res = 1200)
ncm.plot + ncm.an1 + ncm.an2 + ncm.an3 + ncm.an4 + ncm.an5 + ncm.an6
dev.off()

cairo_pdf(
  filename = "figures/diekei_dic_ncm.pdf",
  width = 4,
  height = 5,
  family = "Arial Unicode MS"
)

ncm.plot + ncm.an1 + ncm.an2 + ncm.an3 + ncm.an4 + ncm.an5 + ncm.an6

dev.off()


## EGG ####

# glm_egg2.csv
egg2 <- read.csv('data/02_2026_EE_vis_egg.csv', na.strings = "na")
egg2
attach(egg2)

egg2$pair <- factor(egg2$pair, levels = c("mm", "md", "dd", "dm"))

pair.labs <- as_labeller(c(mm = "MM [(18/280)]", 
                           md = "MD [(12/135)]", 
                           dd = "DD [(16/64)]", 
                           dm = "DM [(18/58)]"), 
                         default = label_parsed)

egg.labs <- data.frame(
  pair = c("mm", "mm", "mm", "md", "md", "md", "dd", "dd", "dd", "dm", "dm", "dm"),
  label = c(
    "italic('M') ~ x ~ italic('M') ~ ': 18(20)'", 
    "eb: '15.77'%+-%'0.37' ^a", 
    "hr: '0.434'%+-%'0.366' ^a",
    "italic('M') ~ x ~ italic('D') ~ ': 12(135)'",
    "eb: '15.91'%+-%'0.59' ^a", 
    "hr: '0.430'%+-%'0.471' ^a",
    "italic('D') ~ x ~ italic('D') ~ ': 16(64)'",
    "eb: '10.41'%+-%'0.53' ^b", 
    "hr: '0.311'%+-%'0.517' ^a",
    "italic('D') ~ x ~ italic('M') ~ ': 18(58)'",
    "eb: '9.64'%+-%'0.59' ^b", 
    "hr: '0.104'%+-%'0.338' ^b"
  ),
  stringsAsFactors = FALSE
)

egg.labs$pair <- factor(egg.labs$pair, levels = c("mm", "md", "dd", "dm"))

egg2_data <- read.csv('Statistics/glm_egg.csv', na.strings = "na")
egg2_data$pair <- factor(egg2_data$pair, levels = c("mm", "md", "dd", "dm"))

egg2_summary <- egg2_data %>%
  group_by(pair) %>%
  summarise(
    mean_prod = mean(prod, na.rm = TRUE),  # Mean of prod
    sem_prod = sd(prod, na.rm = TRUE) / sqrt(sum(!is.na(prod))),  # SEM of prod
    
    mean_htc = mean(htc, na.rm = TRUE),  # Mean of htc
    sem_htc = sd(htc, na.rm = TRUE) / sqrt(sum(!is.na(htc))),  # SEM of htc
    
    mean_hatch_rate = mean(htc / (htc + nhtc), na.rm = TRUE),  # Mean hatching rate
    sem_hatch_rate = sd(htc / (htc + nhtc), na.rm = TRUE) / sqrt(sum(!is.na(htc) & !is.na(nhtc)))  # SEM of hatching rate
  )

egg2_summary$pair <- factor(egg2_summary$pair, levels = c("mm", "md", "dd", "dm"))

egg.plot <- ggplot(egg2, aes(x = prod, y = htc)) + 
  geom_point(aes(shape = pair, color = pair, size = value, alpha = 0.5), 
             position = position_jitter(width = 0.4, height = 0.4, seed = 1830)) + 
  geom_point(data = egg2_summary, aes(x = mean_prod, y = mean_htc, color = pair, shape = pair, fill = pair), size = 2) +
  geom_hline(data = egg2_summary, aes(yintercept = mean_htc, color = pair), linetype = "dashed", size = 0.4) +
  geom_vline(data = egg2_summary, aes(xintercept = mean_prod, color = pair), linetype = "dashed", size = 0.4) +
  facet_wrap(~pair, labeller = labeller(pair = pair.labs)) +
  #gghighlight(unhighlighted_colour = "grey92") +
  #geom_line(aes(color = pair), stat = "smooth", method = lm, alpha = 0.5, size = 2) + 
  #geom_smooth(aes(fill = pair), method = lm, color = NA, alpha = 0.5) +
  scale_y_continuous(limits = c(0, max(egg2$htc, na.rm = TRUE)), expand = expansion(mult = c(0.08, 0.07), add = c(0, 0)), breaks = seq(0, 20, by = 10)) + 
  scale_x_continuous(limits = c(0, max(egg2$prod, na.rm = TRUE)), expand = expansion(mult = c(0, 0.05), add = c(0, 0))) + 
  scale_color_manual(values = c("#BAC94A", "#5D6E1E", "#ABD1DC", "#3B5284")) + 
  scale_fill_manual(values = c("#BAC94A", "#5D6E1E", "#ABD1DC", "#3B5284")) + 
  scale_shape_manual(values = c(21,25,24,22)) +
  xlab("\nEggs produced (per batch)") + 
  ylab("Eggs hatched (per batch)\n") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        #legend.position = "none",
        strip.background = element_rect(fill="white"), 
        strip.text = element_blank(),
        panel.grid = element_blank()) + 
  geom_text(x =c(21,21,21,21,21,21,21,21,21,21,21,21), 
            y = c(25,23,21,25,23,21,25,23,21,25,23,21), 
            aes(label = label), data = egg.labs, size = 2, hjust = 0, 
            color = c("black","black","black",
                      "black","black","black",
                      "black","black","black",
                      "black","black","black"), parse = TRUE)

egg.plot <- ggplot(egg2, aes(x = prod, y = htc)) + 
  geom_point(
    aes(shape = pair, color = pair, size = value),
    alpha = 0.5,
    position = position_jitter(width = 0.4, height = 0.4, seed = 1830)
  ) + 
  geom_point(
    data = egg2_summary,
    aes(x = mean_prod, y = mean_htc, color = pair, shape = pair, fill = pair),
    size = 2, color = "black"
  ) +
  geom_hline(data = egg2_summary, aes(yintercept = mean_htc, color = pair), linetype = "dashed", size = 0.4) +
  geom_vline(data = egg2_summary, aes(xintercept = mean_prod, color = pair), linetype = "dashed", size = 0.4) +
  facet_wrap(~pair, labeller = labeller(pair = pair.labs)) +
  scale_y_continuous(limits = c(0, max(egg2$htc, na.rm = TRUE)), 
                     expand = expansion(mult = c(0.08, 0.07), add = c(0, 0)), 
                     breaks = seq(0, 20, by = 10)) + 
  scale_x_continuous(limits = c(0, max(egg2$prod, na.rm = TRUE)), 
                     expand = expansion(mult = c(0, 0.05), add = c(0, 0))) + 
  scale_size_continuous(name = "n. obs.") +
  scale_color_manual(values = c("#BAC94A", "#5D6E1E", "#ABD1DC", "#3B5284"), guide = "none") + 
  scale_fill_manual(values = c("#BAC94A", "#5D6E1E", "#ABD1DC", "#3B5284"), guide = "none") + 
  scale_shape_manual(values = c(21,25,24,22)) +
  xlab("\nEggs produced (per batch)") + 
  ylab("Eggs hatched (per batch)\n") + 
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.background = element_rect(fill = "white"), 
    strip.text = element_blank(),
    panel.grid = element_blank()
  ) + 
  geom_text(
    x = c(21,21,21,21,21,21,21,21,21,21,21,21), 
    y = c(25,23,21,25,23,21,25,23,21,25,23,21), 
    aes(label = label), data = egg.labs, size = 2, hjust = 0, 
    color = "black", parse = TRUE
  ) +
  guides(
    shape = "none",
    size = guide_legend(order = 2)
  )

egg.plot

egg.an <- function (grob, xmin = -Inf, xmax = Inf, 
                    ymin = -Inf, ymax = Inf, data){ layer(data = data, 
                                                          stat = StatIdentity, 
                                                          position = PositionIdentity, 
                                                          geom = ggplot2:::GeomCustomAnn, 
                                                          inherit.aes = TRUE, 
                                                          params = list(grob = grob, xmin = xmin, xmax = xmax, 
                                                                        ymin = ymin, ymax = ymax))}

egg.an1 <- egg.an(rasterGrob(dic_mm, interpolate=TRUE), xmin=0, xmax=10, ymin=20, ymax=27, data = egg2[194,])
egg.an2 <- egg.an(rasterGrob(dic_dd, interpolate=TRUE), xmin=0, xmax=10, ymin=20, ymax=27, data = egg2[1,])
egg.an3 <- egg.an(rasterGrob(dic_md, interpolate=TRUE), xmin=0, xmax=10, ymin=20, ymax=27, data = egg2[82,])
egg.an4 <- egg.an(rasterGrob(dic_dm, interpolate=TRUE), xmin=0, xmax=10, ymin=20, ymax=27, data = egg2[50,])

egg.plot + egg.an1 + egg.an2 + egg.an3 + egg.an4

ggsave(filename = "figures/diekei_dic_egg.png", width = 5, height = 3.7, device='png', dpi=1200)
ggsave(filename = "figures/diekei_dic_egg.pdf", width = 5, height = 3.7, device='pdf', dpi=1200)


## LARVAL PERFORMANCE ####

lp2 <- read.csv('data/03_2026_EE_vis_larval_performance.csv')

lp2$pair <- factor(lp2$pair, levels = c("mm", "md", "dd", "dm"))
lp2$stg <- factor(lp2$stg, levels = c("Larval acceptance", "Survival to the second instar", "Reaching adulthood"))

lp2_summary <- lp2 %>%
  group_by(stg, pair) %>%
  summarise(
    mean_rsucc.mik = mean(rsucc.mik, na.rm = TRUE),
    sem_rsucc.mik = sd(rsucc.mik, na.rm = TRUE) / sqrt(sum(!is.na(rsucc.mik))), 
    
    mean_rsucc.dic = mean(rsucc.dic, na.rm = TRUE),
    sem_rsucc.dic = sd(rsucc.dic, na.rm = TRUE) / sqrt(sum(!is.na(rsucc.dic))),
  )

lp2_summary$pair <- factor(lp2_summary$pair, levels = c("mm", "md", "dd", "dm"))

pair_labs_expr <- c(
  mm = expression(italic(M) ~ "x" ~ italic(M)),
  md = expression(italic(M) ~ "x" ~ italic(D)),
  dd = expression(italic(D) ~ "x" ~ italic(D)),
  dm = expression(italic(D) ~ "x" ~ italic(M))
)

plot.lp <- ggplot(lp2, aes(x = rsucc.mik, y = rsucc.dic)) +
  geom_point(aes(colour = pair, shape = pair, fill = pair),
             size = 3, alpha = 0.5, color = "white") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey50", size = 0.5) +
  geom_segment(data = lp2_summary,
               aes(x = mean_rsucc.mik, xend = mean_rsucc.mik,
                   y = mean_rsucc.dic - sem_rsucc.dic,
                   yend = mean_rsucc.dic + sem_rsucc.dic, color = pair),
               size = 1) +
  geom_segment(data = lp2_summary,
               aes(y = mean_rsucc.dic, yend = mean_rsucc.dic,
                   x = mean_rsucc.mik - sem_rsucc.mik,
                   xend = mean_rsucc.mik + sem_rsucc.mik, color = pair),
               size = 1) +
  geom_point(data = lp2_summary,
             aes(x = mean_rsucc.mik, y = mean_rsucc.dic,
                 color = pair, shape = pair, fill = pair),
             size = 4, color = "black") +
  facet_wrap(~ stg, ncol = 3) +
  scale_color_manual(
    values = c("#BAC94A", "#5D6E1E", "#ABD1DC", "#3B5284"),
    labels = pair_labs_expr
  ) +
  scale_fill_manual(
    values = c("#BAC94A", "#5D6E1E", "#ABD1DC", "#3B5284"),
    labels = pair_labs_expr
  ) +
  scale_shape_manual(
    values = c(21, 25, 24, 22),
    labels = pair_labs_expr
  ) +
  scale_x_continuous(labels = c("0","0.25","0.5","0.75","1")) +
  scale_y_continuous(labels = c("0","0.25","0.5","0.75","1")) +
  xlab(expression(paste("Survival proportion on ", italic("Mikania")))) +
  ylab(expression(paste("Survival proportion on ", italic("Dicliptera")))) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.background = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.height = unit(0.8, "cm"),
    strip.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    legend.key.size = unit(0.5, "cm"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16),
    axis.title.x = element_text(margin = margin(t = 10, b = 0)),
    axis.title.y = element_text(margin = margin(r = 10, l = 0)),
    ggside.panel.scale = 0.2,
    ggside.axis.text = element_blank(),
    ggside.axis.ticks = element_blank(),
    ggside.panel.border = element_blank()
  )

plot.lp

png("figures/diekei_dic_lp.png", 
    width = 13, height = 5, units = "in", res = 1200)
plot.lp
dev.off()

pdf("figures/diekei_dic_lp.pdf", 
    width = 13, height = 5)
plot.lp
dev.off()


## CHOICE MATING ####

cm.mle <- c("M", "M", "M", "M", "M", "M", "D", "D", "D", "D", "D","D")
cm.fem <- c("mrace", "drace", "mrace", "drace", "mrace", "drace", "drace", "mrace", "drace", "mrace", "drace", "mrace")
cm.arr <- c("ab", "ab", "check", "check", "sep", "sep", "ab", "ab", "check", "check", "sep", "sep")
cm.rate <- c("0.6", "-0.4", "0.875", "-0.125", "0.91", "-0.09", "-0.667", "0.333", "-1", "0", "-0.92", "0.08")
cm.plot <- data.frame(cm.mle, cm.fem, cm.arr, cm.rate)

cm.plot$cm.mle <- factor(cm.plot$cm.mle, levels = c("M", "D"))
cm.plot$cm.arr <- factor(cm.plot$cm.arr, levels = c("ab", "check", "sep"))
cm.plot$cm.fem <- factor(cm.plot$cm.fem, levels = c("mrace", "drace"))

cm.arr.labs <- c("ab" = "Absence", "check" = "Checkerboard", "sep" = "Separate")

cm.mle.labs <- c("\u2642 M-race", "\u2642 D-race")
names(cm.mle.labs) <- c("M", "D")

cm.rate.labs <- data.frame(cm.mle = c("M", "M", "M", "M", "M", "M", "D", "D", "D", "D", "D", "D"), 
                           label = c("b","a","a","b","a","a","b","a","a","b","a","a"))

cm.rate.labs$cm.mle <- factor(cm.rate.labs$cm.mle, levels = c("M", "D"))

cm.hline <- data.frame(cm.mle = c("M", "D"),
                       hline = c(0.5, -0.5))

cm.hline$cm.mle <- factor(cm.hline$cm.mle, levels = c("D", "M"))

asm.plot <- ggplot(cm.plot, aes(fill = cm.fem, y = as.numeric(cm.rate), x = cm.arr)) +
  geom_bar(stat = "identity", width = 0.7) + 
  scale_fill_manual(values = c("#BAC94A", "#ABD1DC")) +
  scale_x_discrete(labels = (text = cm.arr.labs)) +
  xlab("\nHost-plant arrangements") +
  ylab("Rate of assortative mating\n") + 
  guides(fill = guide_legend(title = "fem")) +
  scale_y_continuous(expand=c(0,0),limits=c(-1.25,1.25), breaks = c(-1, -0.5, 0, 0.5, 1), labels = c("-1", "-0.5", "0", "0.5", "1")) +
  facet_wrap(~cm.mle, labeller = labeller(cm.mle = cm.mle.labs)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        strip.text = element_text(face = "italic"),
        panel.grid = element_blank()) +
  geom_text(x = c(0.93, 1.92, 2.94, 1.5, 1.5, 1.5, 0.91, 1.9, 2.91, 1.5, 1.5, 1.5), 
            y = c(0.52, 0.82, 0.855, -1.5, -1.5, -1.5, -0.57, -0.9, -0.83, 1.5, 1.5, 1.5), 
            aes(label = label), data = cm.rate.labs, size = 3.5, hjust = 0, 
            color = c("white","white","white","white","white","white",
                      "white","white","white","white","white","white")) +
  geom_hline(data = cm.hline, aes(yintercept = hline), linetype = "dashed")

asm.plot


asm.an <- function (grob, xmin = -Inf, xmax = Inf, 
                    ymin = -Inf, ymax = Inf, data){ layer(data = data, 
                                                          stat = StatIdentity, 
                                                          position = PositionIdentity, 
                                                          geom = ggplot2:::GeomCustomAnn, 
                                                          inherit.aes = TRUE, 
                                                          params = list(grob = grob, xmin = xmin, xmax = xmax, 
                                                                        ymin = ymin, ymax = ymax))}

asm.an1 <- asm.an(rasterGrob(dic_mm, interpolate=TRUE), xmin=1, xmax=3, ymin=1, ymax=1.25, data = cm.plot[1,])
asm.an2 <- asm.an(rasterGrob(dic_dm, interpolate=TRUE), xmin=1, xmax=3, ymin=-1.25, ymax=-1, data = cm.plot[1,])
asm.an3 <- asm.an(rasterGrob(dic_dd, interpolate=TRUE), xmin=1, xmax=3, ymin=-1.25, ymax=-1, data = cm.plot[7,])
asm.an4 <- asm.an(rasterGrob(dic_md, interpolate=TRUE), xmin=1, xmax=3, ymin=1, ymax=1.25, data = cm.plot[7,])

asm.plot + asm.an1 + asm.an2 + asm.an3 + asm.an4

png("figures/diekei_dic_cm.png", 
    width = 3, height = 5, units = "in", res = 1200)
asm.plot + asm.an1 + asm.an2 + asm.an3 + asm.an4
dev.off()

cairo_pdf(
  filename = "figures/diekei_dic_cm.pdf",
  width = 3,
  height = 5,
  family = "Arial Unicode MS"
)

asm.plot + asm.an1 + asm.an2 + asm.an3 + asm.an4

dev.off()


## SIGHTINGS ####

# glm_fid3.csv
fid3 <- read.csv('data/04_2026_EE_vis_fidelity.csv')
fid3
attach(fid3)

fid3$race <- factor(fid3$race, levels = c("mrace", "drace"))
fid3$arr <- factor(fid3$arr, levels = c("check", "sep"))
fid3$host <- factor(fid3$host,levels = c("mik", "dic"))
fid3$group <- factor(fid3$group, levels = c("mc", "ms", "dc", "ds"))


## SIGHTINGS (5th version)

fid.labs <- c("Separate: M-race", "Separate: D-race", "Checkerboard: M-race", "Checkerboard: D-race")
names(fid.labs) <- c("ms", "ds", "mc", "dc")

fid.plot <- ggplot(data = fid3[fid3$host != "dic",], aes(x = xcor, y = ycor)) + 
  stat_density_2d(aes(fill = ..count..), geom = "raster", contour = FALSE, 
                  h = c(ifelse(bandwidth.nrd(xcor) == 0, 0.1, bandwidth.nrd(xcor)),
                        ifelse(bandwidth.nrd(ycor) == 0, 0.1, bandwidth.nrd(ycor)))) +
  scale_fill_gradientn("mik", colours = c("white", "#BAC94A", "#5D6E1E"), 
                       values = rescale(x = c(0, 40, 86), from = c(0, 86), trans = "sqrt"), 
                       limits = c(0, 86), oob = scales::squish) +
  new_scale("fill") + 
  stat_density_2d(data = subset(fid3, host != "mik"), aes(fill = ..count.., alpha = ..count..), geom = "raster", contour = FALSE, 
                  h = c(ifelse(bandwidth.nrd(xcor) == 0, 0.1, bandwidth.nrd(xcor)),
                        ifelse(bandwidth.nrd(ycor) == 0, 0.1, bandwidth.nrd(ycor)))) +
  scale_alpha_continuous(range = c(0,1), limits = c(0, 1), guide = guide_none()) +
  scale_fill_gradientn("dic", colours = c("white", "#ABD1DC", "#3B5284"), 
                       values = rescale(x = c(0, 10, 40), from = c(0, 40), trans = "sqrt"), 
                       limits = c(0, 86), oob = scales::squish) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6)) +  
  guides(fill = guide_colorbar(reverse = TRUE)) + 
  facet_wrap(~group, labeller = labeller(group = fid.labs)) +
  theme_bw() +
  theme(strip.background = element_blank(), 
        strip.text = element_text(face = "italic"), 
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

fid.plot

fid.an <- function (grob, xmin = -Inf, xmax = Inf, 
                    ymin = -Inf, ymax = Inf, data){ layer(data = data, 
                                                          stat = StatIdentity, 
                                                          position = PositionIdentity, 
                                                          geom = ggplot2:::GeomCustomAnn, 
                                                          inherit.aes = TRUE, 
                                                          params = list(grob = grob, xmin = xmin, xmax = xmax, 
                                                                        ymin = ymin, ymax = ymax))}

fid.an1 <- fid.an(rasterGrob(dic_m, interpolate=TRUE), xmin=0, xmax=1, ymin=4.7, ymax=6, data = fid3[1,])
fid.an2 <- fid.an(rasterGrob(dic_m, interpolate=TRUE), xmin=0, xmax=1, ymin=4.7, ymax=6, data = fid3[341,])
fid.an3 <- fid.an(rasterGrob(dic_d, interpolate=TRUE), xmin=4, xmax=5, ymin=4.7, ymax=6, data = fid3[178,])
fid.an4 <- fid.an(rasterGrob(dic_d, interpolate=TRUE), xmin=4, xmax=5, ymin=4.7, ymax=6, data = fid3[450,])

fid.plot + fid.an1 + fid.an2 + fid.an3 + fid.an4


## DISPERSAL ####

dis2 <- read.csv('data/05_2026_EE_vis_dispersal.csv')
dis2
attach(dis2)

dis2$race <- factor(dis2$race, levels = c("mrace", "drace"))
dis2$group <- factor(dis2$group, levels = c("ms", "mc", "ds", "dc"))
dis2$arr <- factor(dis2$arr, levels = c("check", "sep"))

dis.labs <- c("M-race", "D-race")
names(dis.labs) <- c("mrace", "drace")

dis.plot <- ggplot(data = dis2, aes(x = dis, fill = group, color = group)) +
  geom_density(alpha = 0.6) + 
  scale_fill_manual(values = c("#5D6E1E", "#BAC94A", "#3B5284", "#ABD1DC")) +
  scale_color_manual(values = c("#BAC94A", "#BAC94A", "#ABD1DC", "#ABD1DC")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75), labels = c("0", "0.25", "0.5", "0.75")) +
  xlab(expression("\nDispersal distance (" * italic(m) / italic(d) * ")")) +
  ylab("Density\n") + 
  theme_bw() +
  facet_wrap(~race, labeller = labeller(race = dis.labs)) +
  theme(strip.background = element_blank(), 
        strip.text = element_text(face = "italic"),
        axis.title.y = element_blank(),
        panel.grid = element_blank())

dis.plot

dis.an <- function (grob, xmin = -Inf, xmax = Inf, 
                    ymin = -Inf, ymax = Inf, data){ layer(data = data, 
                                                          stat = StatIdentity, 
                                                          position = PositionIdentity, 
                                                          geom = ggplot2:::GeomCustomAnn, 
                                                          inherit.aes = TRUE, 
                                                          params = list(grob = grob, xmin = xmin, xmax = xmax, 
                                                                        ymin = ymin, ymax = ymax))}

dis.an1 <- dis.an(rasterGrob(dic_m, interpolate=TRUE), xmin=4, xmax=6, ymin=0.71, ymax=0.91, data = dis2[164,])
dis.an2 <- dis.an(rasterGrob(dic_d, interpolate=TRUE), xmin=4, xmax=6, ymin=0.71, ymax=0.91, data = dis2[1,])

dis.an1 <- dis.an(rasterGrob(dic_m, interpolate=TRUE), xmin=4, xmax=6, ymin=0.51, ymax=0.81, data = dis2[164,])
dis.an2 <- dis.an(rasterGrob(dic_d, interpolate=TRUE), xmin=4, xmax=6, ymin=0.51, ymax=0.81, data = dis2[1,])

dis.plot + dis.an1 + dis.an2




hfid1 <- ggarrange(fid.plot + fid.an1 + fid.an2 + fid.an3 + fid.an4,
                   dis.plot + dis.an1 + dis.an2,
                   labels = c("a", "b"),
                   heights = c(4.8, 2),
                   ncol = 1, nrow = 2, align = "v")

hfid1
ggsave(filename = "figures/diekei_dic_hfid.png", width = 5, height = 6, device = 'png', dpi = 1200)
ggsave(filename = "figures/diekei_dic_hfid.pdf", width = 5, height = 6, device = 'pdf', dpi = 1200)



