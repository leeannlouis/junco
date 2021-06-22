# Code used for generation of analyses and figures for the manuscript "Wing and leg bone microstructure reflects migratory demands in resident and migrant populations of the Dark-eyed Junco Junco hyemalis". Uses Junco hyemalis bone microstructure data table, which will be available as a supplement to the article and on Dryad (DOI to follow)

rm(list = ls()); cat("\014"); # Clear the workspace

#### LOAD PACKAGES AND CUSTOM FUNCTIONS ####
if (!require(dplyr)) install.packages('dplyr')
if (!require(ggplot2)) install.packages('ggplot2')
if (!require(lemon)) install.packages('lemon')
if (!require(gridExtra)) install.packages('gridExtra')
if (!require(grid)) install.packages('grid')
if (!require(colorblindr)) install.packages('colorblindr')
if (!require(RColorBrewer)) install.packages('RColorBrewer')
if (!require(ggnewscale)) install.packages('ggnewscale')
if (!require(dplyr)) install.packages('dplyr')
if (!require(tidyr)) install.packages('tidyr')
if (!require(stringr)) install.packages('stringr')
if (!require(nlme)) install.packages('nlme')
if (!require(rstatix)) install.packages('rstatix')
if (!require(Hmisc)) install.packages('Hmisc')
if (!require(corrplot)) install.packages('corrplot')
if (!require(knitr)) install.packages('knitr')
library(dplyr);   library(tidyr);     library(knitr)
library(Hmisc);   library(corrplot);  library(rstatix)
library(nlme);    library(stringr);   library(ggplot2); 
library(lemon);   library(gridExtra); library(grid); 
library(colorblindr); library(RColorBrewer); library(ggnewscale)

# Function for formatting numbers in graphs
scaleFUN <- function(x) sprintf("%.2f", x)

# Create a function that takes a P value and outputs asterisks for its significance
psig = function(pvec) {
  pval = character(length(pvec))
  for (i in 1:length(pvec)) {
    p=pvec[i]
    if (p<=0.05) {pval[i] = '*'}
  }
  return(pval)
}

#### LOAD AND SET UP DATA #####

raw = read.csv(file.choose())
colnames(raw)[1] = 'Museum' # fix import error
dir.export = choose.dir()

# To test, we first need to modify the data I bit. We will:
# -> make subspecies a factor, 
# -> derive a factor for migratory behavior, Behav: resident, Res or migrant, Mig
# -> derive a measurement Twist; see appendix for details. Uses the shear modulus from Reilly and Burstein for bovine fibrolamellar bone (5.1 GPa) and assumes that the TMD can be scaled from the density of bovine bone (2060 g cm-3, from Currey 1979)
data.full <- raw %>% 
  mutate(Subspp = factor(Subspp, levels = c('carolinensis', 'pontilis', 'aikeni', 
            'hyemalis', 'montanus')),
         Behav = factor(if_else(Subspp == 'carolinensis' | Subspp == 'pontilis', 
            'Res', 'Mig'), levels = c('Res', 'Mig'), ordered = FALSE),
         Twist = (Length * Mass * 9.81 * ( (TAr / pi)^(1/2))) / 
           (5.1 *(CtTMD/2060) * J * 1000))

# Reduce the data to the independent variables based on correlation work below
data <- data.full %>% select(c(SpecNum, Behav, Subspp, Location, Elev, Day, 
 Mass, Side, Bone, Length, BVTV, TbN, TbTh, TAr, CtTh, J, CtTMD))

#### AVERAGE VALUES (Tables 1 & 2) ####

# Collect average metadata (elevation, day, mass) by migratory behavior (Table 1)
avgs.meta.migr <- data %>% filter(Bone == 'Humerus') %>%
  select(Behav, Elev, Day, Mass) %>%
  pivot_longer(c('Elev', 'Day', 'Mass'), names_to='Metadata', values_to='Value') %>%
  group_by(Behav, Metadata) %>%
  summarise(Mean = round(mean(Value, na.rm = TRUE), 1), 
            StDev = round(sd(Value, na.rm = TRUE), 1), 
            .groups = 'drop') %>%
  unite(Mean_SD, c('Mean', 'StDev'))  %>%
  pivot_wider(names_from = Metadata, values_from=Mean_SD)
write.csv(avgs.meta.migr, 
          paste(dir.export, '\\MainT1_avg_metadata_by_migration.csv', sep=''))

# Collect average morphological data by migratory behavior (Table 2)
data.avgs <- data.full %>% 
  select(SpecNum, Subspp, Behav, Bone, Length, BVTV, TbN, TbTh, TAr, CtTh, J, CtTMD) %>%
  pivot_longer(c('Length', 'BVTV', 'TbN', 'TbTh', 'TAr', 'CtTh', 'J', 'CtTMD'),
               names_to = 'Morphology', values_to = 'Value') %>%
  group_by(Behav, Bone, Morphology) %>%
  summarise(Mean = signif(mean(Value), 4), StDev = signif(sd(Value), 3), 
            .groups = 'drop') %>%
  mutate(Morphology = factor(Morphology,
    levels = c('Length', 'BVTV', 'TbN', 'TbTh', 'TAr', 'CtTh', 'J', 'CtTMD'), 
    ordered = TRUE)) %>%
  unite(Mean_SD, c('Mean', 'StDev')) %>%
  pivot_wider(names_from = Bone, values_from = Mean_SD) %>%
  arrange(by_group = Morphology) %>%
  relocate(Humerus, .before = Femur)
write.csv(data.avgs, paste(dir.export, '\\MainT2_avg_morpho_by_migration.csv', sep=''))

#### MODEL (Table 3) ####


# Create datasest
dataLog <- data %>% 
  mutate_at(c('Length', 'TbN', 'TbN', 'TbTh', 'TAr', 'CtTh', 'J'), log10) %>%
  mutate(MassLog = log10(Mass)) 

# Chose analyses
analyses = c('Length', 'BVTV', 'TbN', 'TbTh', 'TAr', 'CtTh', 'J', 'CtTMD')
nAnaly = length(analyses)

# Set up matrix to store coefficients of the additive or interaction model. Use 4 rows for each for the three possible predictors (Mass, Behav, Mass*Behav) plus the intercept.
values = data.frame(analysis=c(rep(analyses,each=(4))), 
  effect=rep(c('Intercept', 'Mass', 'Behav', 'Mass*Behav'),nAnaly),
  humValue=numeric(nAnaly*(4)), humSE=numeric(nAnaly*(4)), 
  humP=numeric(nAnaly*(4)), humSig=character(nAnaly*(4)), 
  femValue=numeric(nAnaly*(4)), femSE=numeric(nAnaly*(4)), 
  femP=numeric(nAnaly*(4)), femSig=character(nAnaly*(4)),
  stringsAsFactors=FALSE)

for (i in 1:nAnaly) {

  row = 1 + 4 * (i-1)
  test = analyses[i]
  
  # Create the formula. Use MassLog if the morphological variable is a unit measurement
  if (!any(test == c('CtTMD', 'BVTV'))) {
    form = formula(paste0(test, '~MassLog+Behav'))
  } else {
    form = formula(paste0(test, '~Mass+Behav'))
  }

  # Perform the analyses: run the analysis, extract the coefficients and p-values,
  # and store the coefficients, their standard errors, and the p-values in the table
  
  # Humerus
  # Use an interaction model if the delta AICc for the interaction 
  # was lower for this morphological variable, and if the delta AICc for the additive
  # model was greater than 2 (nesting rule)
  if (any(test == c('Length', 'J', 'CtTMD'))) {
    intForm = formula(paste0(test, '~MassLog*Behav'))
    if (test == 'CtTMD') {intForm = formula(paste0(test, '~Mass*Behav'))}
    hAnaly = summary(lm(intForm, data = dataLog[dataLog$Bone == 'Humerus',]))
    hCoef = as.data.frame(hAnaly$coef)
    values[c(row:(row+3)), c('humValue', 'humSE')] = 
      c(signif(hCoef$Estimate, 4), signif(hCoef$`Std. Error`, 3))
    values[c(row:(row+3)), c('humP')] = signif(hCoef$`Pr(>|t|)`,3)
    values[c(row:(row+3)), c('humSig')] = c(psig(hCoef$`Pr(>|t|)`))
  
  } else {
    hAnaly = summary(lm(form, data = dataLog[dataLog$Bone == 'Humerus',]))
    hCoef = as.data.frame(hAnaly$coef)
    values[c(row:(row+2)), c('humValue', 'humSE')] = 
      c(signif(hCoef$Estimate, 4), signif(hCoef$`Std. Error`, 3))
    values[c(row:(row+2)), c('humP')] = signif(hCoef$`Pr(>|t|)`,3)
    values[c(row:(row+2)), c('humSig')] = c(psig(hCoef$`Pr(>|t|)`))
  }

  # Femur. Only BVTV was better fit by an interaction model
  if (test == 'BVTV') {
    intForm = formula(paste0(test, '~Mass*Behav'))
    fAnaly = summary(lm(intForm, data = dataLog[dataLog$Bone == 'Femur',]))
    fCoef = as.data.frame(fAnaly$coef)
    values[c(row:(row+3)), c('femValue', 'femSE')] = 
      c(signif(fCoef$Estimate, 4), signif(fCoef$`Std. Error`, 3))
    values[c(row:(row+3)), c('femP')] = signif(fCoef$`Pr(>|t|)`,3)
    values[c(row:(row+3)), c('femSig')] = c(psig(fCoef$`Pr(>|t|)`))
  
  } else {
    fAnaly = summary(lm(form, data = dataLog[dataLog$Bone == 'Femur',]))
    fCoef = as.data.frame(fAnaly$coef)
    values[c(row:(row+2)), c('femValue', 'femSE')] = 
      c(signif(fCoef$Estimate, 4), signif(fCoef$`Std. Error`, 3))
    values[c(row:(row+2)), c('femP')] = signif(fCoef$`Pr(>|t|)`,3)
    values[c(row:(row+2)), c('femSig')] = c(psig(fCoef$`Pr(>|t|)`))
  }

}

# Save to use 
write.csv(values, paste0(dir.export, '\\MainT3_coefficients.csv'))

#### PLOTTED DATA (Figures 4, 5, 6) #####

# To check if colorblind friendly:
#display.brewer.pal(n=11, name='PRGn')
#brewer.pal(n=11, name='PRGn')
#cvd_grid(myplot)

# Function to make point graphs
pointMaker = function(myVar, myBone) {
  myData <- data.full %>% filter(Bone == myBone) %>%
    select(Mass, Subspp, Behav, myVar) %>% rename(Value = myVar)
  myPlot <- ggplot(myData, aes(x = Mass, y = Value)) +
    geom_smooth(aes(x = Mass, y = Value, color = Behav), 
                method = 'lm', alpha = 0.35, data = myData) +
    geom_point(aes(fill = Subspp), shape = 21, size = 1.5, 
               color = 'black', stroke = 0.4) + 
    scale_fill_manual(values = c('#9970AB', '#762A83', '#001a0a', '#A6DBA0', '#D9F0D3'),
                      name='Subspecies') +
    scale_color_manual(labels=c('Resident', 'Migrant'), values=c('#C2A5CF', '#1B7837'))
    xlab('Mass (g)') + 
    theme(legend.position = 'none', axis.title.y=element_blank(),
          plot.title = element_text(hjust = 0.5, size=rel(1.1)), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black"), 
          axis.text.x = element_text(size=rel(1.35)), 
          axis.text.y = element_text(size=rel(1.35)))
  return(myPlot)
}

# Function to make bar plots
barMaker = function(myVar, myBone) {
  myData <- data.full %>% filter(Bone == myBone) %>%
    mutate(BVTV = BVTV * 100) %>%
    select(Mass, Subspp, Behav, myVar) %>% rename(Value = myVar)
  myAvgs <- myData %>% group_by(Behav) %>%
    dplyr::summarise(Mean = mean(Value), StDev = sd(Value), .groups = 'drop')
  myPlot <- ggplot(myAvgs, aes(x = Behav, y = Mean, fill = Behav)) + 
    geom_bar(stat = 'identity') + 
    geom_errorbar(aes(ymin=Mean-StDev, ymax=Mean+StDev), width=0.4, size=1.3) +
    scale_fill_manual(values=c('#C2A5CF', '#1B7837'), name = '', 
                      labels = c('Resident', 'Migrant')) + 
    new_scale('fill') +
    geom_jitter(aes(x=Behav, y=Value, fill=Subspp), data=myData, 
                width=0.1, color='black',
                shape = 21, size = 1.5, stroke = 0.4) +
    scale_fill_manual(values = c('#9970AB', '#762A83', '#001a0a', '#A6DBA0', '#D9F0D3'),
                      name = 'Subspecies') +
    scale_x_discrete(labels = c('Resident', 'Migrant')) + 
    scale_y_continuous(expand = c(0, 0)) +
    xlab('') + expand_limits(y = max(myData$Value) * 1.05) +
    theme(legend.position = 'none', axis.title.y=element_blank(),
          plot.title = element_text(hjust = 0.5, size=rel(1.1)), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size=rel(1.35)), 
          axis.text.y = element_text(size=rel(1.35)))
  return(myPlot)
}

# Grab the legend
my_legend = g_legend(barMaker('CtTh', 'Humerus') +
         guides(fill=guide_legend(label.theme=element_text(face='italic'),
               override.aes = list(size=3)))+
         theme(legend.position="right",
               legend.text=element_text(size=rel(1.1)), legend.key.size=unit(0.1, 'in'),
               legend.key.width=unit(0.3, 'in'), legend.key = element_blank(),
               legend.margin=margin(0.02, 0.02,0.02, 0.02,'in')))

# Length (Figure 4)
png(paste(dir.export, '\\MainF4_length.png', sep=''),
    width=6.5, height=2.25, res=1200, units='in')
grid.arrange(
  my_legend,
  arrangeGrob(
    pointMaker('Length', 'Humerus') + 
      labs(title = bquote('Humerus'), y = bquote('Length'~(mm))) + 
      theme(axis.title.y=element_text(), legend.position = 'none'),
    barMaker('Length', 'Femur') + 
      labs(title = bquote('Femur'), y = bquote('Length'~(mm))) + 
      theme(axis.title.y=element_text()), nrow=1),
  ncol=2, widths=c(1, 4))
dev.off(); dev.off()

# Trabecular Morphology (Figure 5)
png(paste(dir.export, '\\MainF5_trabecular.png', sep=''),
    width=6.5, height=3.25, res=1200, units='in')
grid.arrange(my_legend,
             arrangeGrob(
               arrangeGrob(
                 pointMaker('BVTV', 'Humerus') + labs(title = bquote('BV/TV'~('%'))) +
                   theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                         legend.position = 'none'),
                 barMaker('TbTh', 'Humerus') + 
                   labs(title = bquote('Tb.Th'~(mm)), y = bquote('Tb.Th'~(mm))) +
                   theme(axis.title.x=element_blank(), axis.text.x=element_blank()),
                 nrow=1, left=textGrob('Humerus', rot=90, gp=gpar(fontface='bold'))),
               arrangeGrob(
                 pointMaker('BVTV', 'Femur') + theme(legend.position = 'none'), 
                 barMaker('TbTh', 'Femur'),
                 nrow=1, left=textGrob('Femur', rot=90, gp=gpar(fontface='bold'))),
               nrow=2),
             ncol=2, widths=c(1, 4))
dev.off(); dev.off()

# Cortical Morphology (Figure 6)
png(paste(dir.export, '\\MainF6_cortical.png', sep=''),
    width=6.5, height=6.5, res=1200, units='in')
grid.arrange(
  arrangeGrob(
    my_legend,      
    barMaker('CtTh', 'Humerus') + labs(title = 'Humerus', y = bquote('Ct.Th'~(mm))) +
      theme(axis.title.x=element_blank(), axis.title.y=element_text()),      
    barMaker('CtTh', 'Femur') + labs(title='Femur', y = bquote('Ct.Th'~(mm))) + 
      theme(axis.title.x=element_blank(), axis.title.y=element_text()) + 
      scale_y_continuous(labels=scaleFUN), nrow=1),
  arrangeGrob(
    pointMaker('TAr', 'Humerus') + labs(title = bquote('T.Ar'~(mm^2))) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
            legend.position = 'none'),
    pointMaker('J', 'Humerus') + labs(title = bquote('J'~(mm^4))) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
            legend.position = 'none') + 
      scale_y_continuous(labels=scaleFUN),
    pointMaker('CtTMD', 'Humerus') + labs(title = bquote('TMD '(mg ~cm^-3))) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
            legend.position = 'none'),
    nrow=1, left=textGrob('Humerus', rot=90, gp=gpar(fontface='bold'))),
  arrangeGrob(
    pointMaker('TAr', 'Femur') + theme(legend.position = 'none'),
    pointMaker('J', 'Femur') + scale_y_continuous(labels=scaleFUN) + 
      theme(legend.position = 'none'), 
    pointMaker('CtTMD', 'Femur') + theme(legend.position = 'none'),
    nrow=1, left=textGrob('Femur', rot=90, gp=gpar(fontface='bold'))),
  nrow=3)
dev.off(); dev.off()

# Twist (Figure 7)
png(paste(dir.export, '\\MainF7_twist.png', sep=''),
    width=6.5, height = 3.25, res=1200, units='in')
grid.arrange(
  my_legend,
  arrangeGrob(
    pointMaker('Twist', 'Humerus') + labs(title = 'Humerus') + 
      theme(legend.position = 'none'),
    pointMaker('Twist', 'Femur') + labs(title = 'Femur') + 
      theme(legend.position = 'none'), nrow=1,
    left = textGrob(bquote('Twist Angle '~theta), rot=90, gp=gpar(fontface="bold"))),
  ncol=2, widths=c(1, 4))
dev.off(); dev.off()

#### SUPPLEMENT: CORRELATIONS (Figure S1) #####

# Which variables should we analyze? Only keep factors with < 50% correlation.
# Cortical: TAr, CtTh, CtTMD, but also J
# Trabecular: BVTV, TbN, TbTh. 

# Format data for correlation work
data.corr <- data.full %>% rename(BV = TbBV, 'BV/TV' = BVTV, 'Conn.D' = ConnDens, 
      'Tb.N' = TbN, 'Tb.Th' = TbTh, 'Tb.Sp' = TbSp,TMD = CtTMD, 'Tt.Ar' = TAr, 
      'Ct.Ar' = CtAr, 'Ma.Ar' = MaAr, 'Ct.Ar/Tt.Ar' = CtArTAr, 'Ct.Th' = CtTh)

# Get the correlations
ct.cor.h <- data.corr %>% filter(Bone == 'Humerus') %>% 
  select('Tt.Ar', 'Ct.Ar', 'Ma.Ar', 'Ct.Ar/Tt.Ar', 'Ct.Th', J, Imax, Imin, TMD) %>% 
  as.matrix() %>% rcorr()
ct.cor.f <- data.corr %>% filter(Bone == 'Femur') %>% 
  select('Tt.Ar', 'Ct.Ar', 'Ma.Ar', 'Ct.Ar/Tt.Ar', 'Ct.Th', J, Imax, Imin, TMD) %>% 
  as.matrix() %>% rcorr()
tb.cor.h <- data.corr %>% filter(Bone == 'Humerus') %>% 
  select(BV, 'BV/TV', 'Conn.D', SMI, 'Tb.N', 'Tb.Th', 'Tb.Sp') %>% as.matrix() %>% 
  rcorr()
tb.cor.f <- data.corr %>% filter(Bone == 'Femur') %>% 
  select(BV, 'BV/TV', 'Conn.D', SMI, 'Tb.N', 'Tb.Th', 'Tb.Sp') %>% as.matrix() %>% 
  rcorr()

# Send value to 0 if the correlation is < 50%
ct.cor.h$r[abs(ct.cor.h$r) < 0.5] = 0; ct.cor.f$r[abs(ct.cor.f$r) < 0.5] = 0;
tb.cor.h$r[abs(tb.cor.h$r) < 0.5] = 0; tb.cor.f$r[abs(tb.cor.f$r) < 0.5] = 0;

# Create plots
png(paste0(dir.export, '\\SuppF1_correlations.png'), width=5, 
    height=5, res=600, units='in', type = 'cairo') # export the image, if desired
par(oma = rep(0, 4), mfrow = c(2, 2))
corrplot(ct.cor.h$r[,2:9], type="upper", p.mat = ct.cor.h$P[,2:9], sig.level = 0.05, 
         cl.cex = 0.75, cl.align.text = 'r', cl.offset = -0.5, insig = 'n', 
         mar = c(0, 0, 2, 0), title = 'Humerus Cortical', tl.col = 'black')
corrplot(ct.cor.f$r[,2:9], type="upper", p.mat = ct.cor.f$P[,2:9], sig.level = 0.05, 
         cl.cex = 0.75, cl.align.text = 'r', cl.offset = -0.5, insig = 'n', 
         mar = c(0, 0, 2, 0), title = 'Femur Cortical', tl.col = 'black')
corrplot(tb.cor.h$r[,2:7], type="upper", p.mat = tb.cor.h$P[,2:7], sig.level = 0.05,
         insig = 'n', cl.cex = 0.75, cl.align.text = 'r', cl.offset = -0.5,
         title = 'Humerus Trabecular', tl.col = 'black', mar = c(0, 0, 2, 0))
corrplot(tb.cor.f$r[,2:7], type="upper", p.mat = tb.cor.f$P[,2:7], sig.level = 0.05,
         insig = 'n', cl.cex = 0.75, cl.align.text = 'r', cl.offset = -0.5, 
         mar = c(0, 0, 2, 0), title = 'Femur Trabecular', tl.col = 'black')

dev.off(); dev.off() # if exported the figure, close it

#### SUPPLEMENT: ADDITIONAL AVERAGE VALUES (Tables S3, S4, S5) ####

# Collect average metadata (elevation, day, mass) by subspecies behavior (Table S3)
avgs.meta.subspp <- data %>% filter(Bone == 'Humerus') %>%
  select(Subspp, Elev, Day, Mass) %>%
  pivot_longer(c('Elev', 'Day', 'Mass'), names_to='Metadata', values_to='Value') %>%
  group_by(Subspp, Metadata) %>%
  summarise(Mean = round(mean(Value, na.rm = TRUE), 1), 
            StDev = round(sd(Value, na.rm = TRUE), 1), 
            .groups = 'drop') %>%
  unite(Mean_SD, c('Mean', 'StDev'))  %>%
  pivot_wider(names_from = Metadata, values_from=Mean_SD)
write.csv(avgs.meta.subspp, 
          paste(dir.export, '\\SuppT3_avg_metadata_by_subspp.csv', sep=''))

# Collect additional morphological data by migratory behavior (Table S4)
data.avgs <- data.full %>% 
  select(SpecNum, Subspp, Behav, Bone, TbBV, ConnDens, 
         SMI, TbSp, CtAr, MaAr, CtArTAr) %>%
  pivot_longer(c('TbBV', 'ConnDens', 'SMI', 
                 'TbSp', 'CtAr', 'MaAr', 'CtArTAr'),
               names_to = 'Morphology', values_to = 'Value') %>%
  group_by(Behav, Bone, Morphology) %>%
  summarise(Mean = signif(mean(Value), 4), StDev = signif(sd(Value), 3), 
            .groups = 'drop') %>%
  mutate(Morphology = factor(Morphology,
         levels = c('TbBV', 'ConnDens', 'SMI', 
                    'TbSp', 'CtAr', 'MaAr', 'CtArTAr'), ordered = TRUE)) %>%
  unite(Mean_SD, c('Mean', 'StDev')) %>%
  pivot_wider(names_from = Bone, values_from = Mean_SD) %>%
  arrange(by_group = Morphology) %>%
  relocate(Humerus, .before = Femur)
write.csv(data.avgs, paste(dir.export, '\\SuppT4_avg_morpho_by_migration.csv', sep=''))

# Collect average morphological data by subspecies behavior (Table S5)
data.avgs <- data.full %>% 
  select(SpecNum, Subspp, Behav, Bone, Length, TbBV, BVTV, ConnDens, SMI, 
         TbN, TbTh, TbSp, TAr, CtAr, MaAr, CtArTAr, CtTh, J, CtTMD) %>%
  pivot_longer(c('Length', 'TbBV', 'BVTV', 'ConnDens', 'SMI', 
                 'TbN', 'TbTh', 'TbSp', 'TAr', 'CtAr', 'MaAr', 'CtArTAr', 
                 'CtTh', 'J', 'CtTMD'),
               names_to = 'Morphology', values_to = 'Value') %>%
  group_by(Subspp, Bone, Morphology) %>%
  summarise(Mean = signif(mean(Value), 4), StDev = signif(sd(Value), 3), 
            .groups = 'drop') %>%
  mutate(Morphology = factor(Morphology,
      levels = c('Length', 'TbBV', 'BVTV', 'ConnDens', 'SMI', 'TbN', 'TbTh', 'TbSp', 
        'TAr', 'CtAr', 'MaAr', 'CtArTAr', 'CtTh', 'J', 'CtTMD'),ordered = TRUE)) %>%
  unite(Mean_SD, c('Mean', 'StDev')) %>%
  pivot_wider(names_from = Bone, values_from = Mean_SD) %>%
  arrange(by_group = Morphology) %>%
  relocate(Humerus, .before = Femur)
write.csv(data.avgs, paste(dir.export, '\\SuppT5_avg_morpho_by_subspecies.csv', sep=''))

#### SUPPLEMENT: MODEL WITH AIC (Table S6) #####

# Set up the data for the AIC analyses
dataLog <- data.full %>% 
  mutate_at(c('Length', 'TbBV', 'ConnDens', 'TbN', 'TbTh', 
              'TbSp', 'TAr', 'CtAr', 'MaAr', 'CtTh', 'J'), log10) %>%
  mutate(MassLog = log10(Mass)) 

# Create a table to store AIC results
analyses = c('Length', 'TbBV', 'BVTV', 'ConnDens', 'SMI', 'TbN', 'TbTh', 
             'TbSp', 'TAr', 'CtAr', 'MaAr', 'CtArTAr', 'CtTh', 'J', 'CtTMD')
nAnaly = length(analyses)
scores = data.frame(bone=c(rep('Humerus', nAnaly), rep('Femur', nAnaly)),
                    analysis=rep(analyses, 2), additive=numeric(nAnaly*2), 
                    interaction=numeric(nAnaly*2), stringsAsFactors=FALSE)

# Calculate AIC scores for the additive and interaction models for each morphological response variable for each bone (humerus, femur)
for (i in 1:2) {
  for (j in 1:nAnaly) {
    bone = c('Humerus', 'Femur')[i]
    row = (i-1)*nAnaly + j
    test = analyses[j]
    
    # Use the log of mass for all tests except BMD and dimensionless variables (BVTV)
    massVar = 'MassLog'
    if (any(test == c('BVTV', 'SMI', 'CtArTAr', 'CtTMD'))) {massVar = 'Mass'}
    
    # Calculate AIC
    likAdd = AICcmodavg::AICc(lm(formula(paste(test, '~', massVar, '+ Behav')),
                                 data=dataLog[dataLog$Bone == bone,]))
    likCross = AICcmodavg::AICc(lm(formula(paste(test, '~', massVar, '* Behav')),
                                   data=dataLog[dataLog$Bone == bone,]))
    
    # Calculate delta AIC and add to matrix
    lowest = min(likAdd, likCross)
    scores[row, c(3)] = likAdd - lowest
    scores[row, c(4)] = likCross - lowest
  }
}

# Code to save as a .csv if desired
write.csv(scores, paste(dir.export, '\\SuppT6_AICc.csv', sep=''))


#### SUPLLEMENT: ADDITIONAL MODEL DATA (Table S7) ####

# Create datasest
dataLog <- data.full %>% 
  mutate_at(c('TbBV', 'ConnDens', 'TbSp', 'CtAr', 'MaAr'), log10) %>%
  mutate(MassLog = log10(Mass)) 

# Chose analyses
analyses = c('TbBV', 'ConnDens', 'SMI', 'TbSp', 'CtAr', 'MaAr', 'CtArTAr')
nAnaly = length(analyses)

# Set up matrix to store coefficients of the additive or interaction model. Use 4 rows for each for the three possible predictors (Mass, Behav, Mass*Behav) plus the intercept.
values = data.frame(analysis=c(rep(analyses,each=(4))), 
                    effect=rep(c('Intercept', 'Mass', 'Behav', 'Mass*Behav'),nAnaly),
                    humValue=numeric(nAnaly*(4)), humSE=numeric(nAnaly*(4)), 
                    humP=numeric(nAnaly*(4)), humSig=character(nAnaly*(4)), 
                    femValue=numeric(nAnaly*(4)), femSE=numeric(nAnaly*(4)), 
                    femP=numeric(nAnaly*(4)), femSig=character(nAnaly*(4)),
                    stringsAsFactors=FALSE)

for (i in 1:nAnaly) {
  
  row = 1 + 4 * (i-1)
  test = analyses[i]
  
  # Create the formula. Use MassLog if the morphological variable is a unit measurement
  if (!any(test == c('SMI', 'CtArTAr'))) {
    form = formula(paste0(test, '~MassLog+Behav'))
  } else {
    form = formula(paste0(test, '~Mass+Behav'))
  }
  
  # Perform the analyses: run the analysis, extract the coefficients and p-values,
  # and store the coefficients, their standard errors, and the p-values in the table
  
  # Humerus
  # Use an interaction model if the delta AICc for the interaction 
  # was lower for this morphological variable, and if the delta AICc for the additive
  # model was greater than 2 (nesting rule)
  if (any(test == c('CtAr'))) {
    intForm = formula(paste0(test, '~MassLog*Behav'))
    hAnaly = summary(lm(intForm, data = dataLog[dataLog$Bone == 'Humerus',]))
    hCoef = as.data.frame(hAnaly$coef)
    values[c(row:(row+3)), c('humValue', 'humSE')] = 
      c(signif(hCoef$Estimate, 4), signif(hCoef$`Std. Error`, 3))
    values[c(row:(row+3)), c('humP')] = signif(hCoef$`Pr(>|t|)`,3)
    values[c(row:(row+3)), c('humSig')] = c(psig(hCoef$`Pr(>|t|)`))
    
  } else {
    hAnaly = summary(lm(form, data = dataLog[dataLog$Bone == 'Humerus',]))
    hCoef = as.data.frame(hAnaly$coef)
    values[c(row:(row+2)), c('humValue', 'humSE')] = 
      c(signif(hCoef$Estimate, 4), signif(hCoef$`Std. Error`, 3))
    values[c(row:(row+2)), c('humP')] = signif(hCoef$`Pr(>|t|)`,3)
    values[c(row:(row+2)), c('humSig')] = c(psig(hCoef$`Pr(>|t|)`))
  }
  
  # Femur
  fAnaly = summary(lm(form, data = dataLog[dataLog$Bone == 'Femur',]))
  fCoef = as.data.frame(fAnaly$coef)
  values[c(row:(row+2)), c('femValue', 'femSE')] = 
    c(signif(fCoef$Estimate, 4), signif(fCoef$`Std. Error`, 3))
  values[c(row:(row+2)), c('femP')] = signif(fCoef$`Pr(>|t|)`,3)
  values[c(row:(row+2)), c('femSig')] = c(psig(fCoef$`Pr(>|t|)`))
  
}

# Save to use 
write.csv(values, paste0(dir.export, '\\SuppT7_coefficients.csv'))
