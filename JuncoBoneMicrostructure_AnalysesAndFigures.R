# Takes the junco bone microstructure .xlsx data file available on Dryad
# (DOI will be provided upon manuscript acceptance)
# and outputs spreadsheets of AICc values, effect sizes,
# and figures into user-provided folders. 
# These match the tables and figures in the article (DOI will be
# provided upon manuscript acceptance).

#### SETUP #####

# Clear workspace
rm(list = ls()); cat("\014"); dev.off()

library('ggplot2'); library('dplyr'); library('ggpubr'); library('scales'); library('lemon')
library('lme4'); library('nlme'); library(grid); library(gridExtra); library(openxlsx); library(viridis)
library('extrafont'); library('scales')

# Load data. Stored in JuncoBoneMicrostructure_Results.xlsx
print('Choose file to use for data.'); myFile = file.choose() 
raw = read.xlsx(myFile)

# Fix the importation of dates
raw$DateCollected = as.Date(as.numeric(raw$Date), origin=as.Date('1899-12-30'))

# Choose directory for data exports and figure exports
print('Choose directory to export data'); dir.export.data = choose.dir()
print('Choose directory to export figures'); dir.export.figs = choose.dir()

# Set up basic dataset. Only take males. 
# To only take complete cases: data = data[complete.cases(data), ]
data = raw[,c('Museum', 'SpecNum', 'Name', 'Group', 'Bone', 'Sex', 'Mass', 'DateCollected', 'Day', 'Location', 'Alti',
              'Length', 'BVTV', 'TbBV', 'TbBS', 'TbN', 'TbTh', 'TbSp', 'SMI', 'ConnDens', 'TA', 'BA', 'MA',
              'BATA', 'CtTMD', 'CtBV', 'CtTh', 'do', 'pMOI', 'Imax', 'Imin', 'Cmax')]
data = data[data$Sex=='M', -which(names(data)=='Sex')]  # only take males
data = data[-which(is.na(data$Mass)),]                  # only use animals with mass data

# Create treatment groups
data$Treat[data$Group %in% c('pinosus', 'carolinensis', 'pontilis')] = 'Res'
data$Treat[data$Group %in% c('hyemalis', 'montanus', 'aikeni')] = 'Mig'
data$Treat = factor(data$Treat, c('Res', 'Mig'))

# Separate the two bone types
data$Bone = factor(data$Bone, c('Humerus', 'Femur'))
hum = data[data$Bone=='Humerus',]; 
fem = data[data$Bone=='Femur',]

# We must have both humerus and femur data for all entries for the data to make sense. 
# Remove any humerus entries that are not in the femur data, and vice versa.
inFemNotHum = !fem$SpecNum %in% hum$SpecNum
if (any(inFemNotHum)) {
  print(sprintf('These specimens have femur data but not humerus data, and will be dropped: %s',
                paste(as.character(fem$SpecNum[inFemNotHum]), collapse=", ")))
  fem = fem[!inFemNotHum,]
}
inHumNotFem = !hum$SpecNum %in% fem$SpecNum
if (any(inHumNotFem)) {
  print(sprintf('These specimens have humerus data but not femur data, and will be dropped: %s ',
                paste(as.character(hum$SpecNum[inHumNotFem]), collapse=", ")))
  hum = hum[!inHumNotFem]
}

# Set up log-log data
unitMetrics = c('Alti', 'Length', 'TbBV', 'TbBS', 'TbN', 'TbTh', 'TbSp', 'ConnDens', 'TA', 'BA', 'MA',
                             'CtTMD', 'CtBV', 'CtTh', 'do', 'pMOI', 'Imax', 'Imin', 'Cmax')
dataLog = data; dataLog[,unitMetrics] = log(dataLog[,unitMetrics])
humLog = hum; humLog[,unitMetrics] = log(humLog[,unitMetrics])
femLog = fem; femLog[,unitMetrics] = log(femLog[,unitMetrics])
dataLog$MassLog = log(dataLog$Mass); humLog$MassLog = log(humLog$Mass); femLog$MassLog = log(femLog$Mass)
dataLog$AltiLog = log(dataLog$Alti); humLog$AltiLog = log(humLog$Alti); femLog$AltiLog = log(femLog$Alti)

# How many birds, and where are they from?
nrow(hum[(hum$Museum=='AMNH')&(hum$Group=='hyemalis'),])
nrow(hum[(hum$Museum=='AMNH')&(hum$Group=='carolinensis'),])
nrow(hum[(hum$Museum=='AMNH')&(hum$Group=='montanus'),])
nrow(hum[(hum$Museum=='AMNH')&(hum$Group=='aikeni'),])
nrow(hum[(hum$Museum=='AMNH')&(hum$Group=='pontilis'),])

# List the analyses of interest
analyses = c('Length', 'BVTV', 'TbBV', 'TbN', 'TbTh', 'TbSp', 'SMI', 'ConnDens',
             'BATA', 'CtTMD', 'CtBV', 'CtTh', 'do', 'pMOI', 'Imax', 'Imin')
nAnaly = length(analyses)

# Reorder levels
hum$Group = factor(hum$Group, levels = c('carolinensis', 'pontilis', 'aikeni', 'hyemalis', 'montanus'))
humLog$Group = factor(hum$Group, levels = c('carolinensis', 'pontilis', 'aikeni', 'hyemalis', 'montanus'))
fem$Group = factor(hum$Group, levels = c('carolinensis', 'pontilis', 'aikeni', 'hyemalis', 'montanus'))
femLog$Group = factor(hum$Group, levels = c('carolinensis', 'pontilis', 'aikeni', 'hyemalis', 'montanus'))
subspp_labels = c('J. h. carolinensis', 'J. h. pontilis', 'J. h. aikeni', 'J. h. hyemalis', 'J. h. montanus')

# Create plotting theme
transparent = rgb(0, 0, 0, max = 255, alpha = 0, names = "transparent")

mytheme = theme(axis.text = element_text(family="Tahoma", size=10, colour='black'),
                     axis.title = element_text(family="Tahoma", size=11, colour='black'), 
                     plot.title = element_blank(),
                     legend.position='none', panel.background = element_blank(),
                     axis.text.y=element_text(angle=90), axis.text.y.left = element_text(hjust = 0.5),
                     panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                     axis.ticks=element_line(color='black', size=2),
                     axis.line=element_line(color='black'))

# Function to print *'s for significance
psig = function(pvec) {
  pval = character(length(pvec))
  for (i in 1:length(pvec)) {
    p=pvec[i]
    if (p<=0.05 & p>0.01){
      pval[i] = '*'
    } else if (p<=0.01 & p>0.001) {
      pval[i] = '**'
    } else if (p<=0.001) {
      pval[i] = '***'
    }
  }
  return(pval)
}

# Function to print custom number of decimals in plots
scaleFUN1 <- function(x) sprintf("%.1f", x)
scaleFUN2 <- function(x) sprintf("%.2f", x)
scaleFUN3 <- function(x) sprintf("%.3f", x)

#### VALUES: AICC ANALYSES TO DEMONSTRATE WHICH TERMS ARE MEANINGFUL ####

getModel = function(theData, yval) {
  massVar = 'MassLog'
  if (any(yval == c('BVTV', 'SMI', 'BATA'))) {
    massVar = 'Mass'
  }
  form = formula(paste(test, '~', massVar, '+Treat', sep=''))
  
  df = theData[,c('Mass','MassLog', 'Treat', 'Group', 'Alti', 'AltiLog', 'Day', test)] 
  df = df[complete.cases(df),]

  form01 = formula(paste(yval, '~', massVar))
  form02 = formula(paste(yval, '~Group'))
  form03 = formula(paste(yval, '~Treat'))
  form04 = formula(paste(yval, '~', massVar, '+Group'))
  form05 = formula(paste(yval, '~', massVar, '+Treat'))
  form06 = formula(paste(yval, '~', massVar, '+Treat+Day'))
  form07 = formula(paste(yval, '~', massVar, '+Treat+Alti'))
  form08 = formula(paste(yval, '~', massVar, '+Treat+Day+Alti'))
  form09 = formula(paste(yval, '~', massVar, '*Treat'))
  
  lik01 = AICcmodavg::AICc(gls(form01, data=df))
  lik02 = AICcmodavg::AICc(gls(form02, data=df)); lik03 = AICcmodavg::AICc(gls(form03, data=df))
  lik04 = AICcmodavg::AICc(gls(form04, data=df)); lik05 = AICcmodavg::AICc(gls(form05, data=df))
  lik06 = AICcmodavg::AICc(gls(form06, data=df)); lik07 = AICcmodavg::AICc(gls(form07, data=df))
  lik08 = AICcmodavg::AICc(gls(form08, data=df)); lik09 = AICcmodavg::AICc(gls(form09, data=df))
  
  allmodels = c(lik01, lik02, lik03, lik04, lik05, lik06, lik07, lik08, lik09)
  
  lowest = min(allmodels)
  return(allmodels - lowest)
}

scores = data.frame(bone=c(rep('Humerus', nAnaly), rep('Femur', nAnaly)),
                    analysis=rep(analyses, 2), 
                    massOnly=numeric(nAnaly*2), groupOnly=numeric(nAnaly*2), 
                    treatOnly=numeric(nAnaly*2), massGroupAdd=numeric(nAnaly*2),
                    massTreatAdd=numeric(nAnaly*2), massTreatDayAdd=numeric(nAnaly*2),
                    massTreatAltiAdd=numeric(nAnaly*2), massTreatAltiDayAdd=numeric(nAnaly*2),
                    massTreatCross=numeric(nAnaly*2), stringsAsFactors=FALSE)
for (i in 1:2) {
  for (j in 1:nAnaly) {
    if (i==1) {dataset = humLog} else {dataset = femLog}
    row = (i-1)*nAnaly + j
    test = analyses[j]
    theScores = getModel(dataset, test)
    scores[row, c(3:11)] = theScores
    #scores$constant[row] = theScores[1]
  }
}

# Save to use 
write.csv(scores, paste(dir.export.data, '/AICc_Scores_Parameter_Selection.csv', sep=''))

#### VALUES: COMPARE MASS*TREAT AND MASS+TREAT MODELS ####

scores = data.frame(bone=c(rep('Humerus', nAnaly), rep('Femur', nAnaly)),
                    analysis=rep(analyses, 2), massTreatAdd=numeric(nAnaly*2), 
                    massTreatCross=numeric(nAnaly*2), stringsAsFactors=FALSE)
for (i in 1:2) {
  for (j in 1:nAnaly) {
    
    # Setup dataset, row, and which mass variable to use
    if (i==1) {dataset = humLog} else {dataset = femLog}
    row = (i-1)*nAnaly + j
    test = analyses[j]
    massVar = 'MassLog'
    if (any(test == c('BVTV', 'SMI', 'BATA'))) {
      massVar = 'Mass'
    }
    df = dataset[,c('Mass','MassLog', 'Treat', 'Group', 'Alti', 'AltiLog', 'Day', test)] 
    df = df[complete.cases(df),]
    
    formAdd = formula(paste(test, '~',  massVar, '+Treat'))
    likAdd = AICcmodavg::AICc(gls(formAdd, data=df))
    
    formCross = formula(paste(test, '~', massVar, '*Treat'))
    likCross = AICcmodavg::AICc(gls(formCross, data=df))
    
    lowest = min(likAdd, likCross)
    scores[row, c(3)] = likAdd - lowest
    scores[row, c(4)] = likCross - lowest
  }
}

# Save to use 
write.csv(scores, paste(dir.export.data, '/AICc_Scores_Add_vs_Cross.csv', sep=''))

#### VALUES: TABLE OF EFFECT SIZES ####

finalAnalyses = c(analyses) #, 'strBendLog', 'strTorsLog', 'relStr')
nFinalAnaly = length(finalAnalyses)
values = data.frame(analysis=c(rep(finalAnalyses,each=2)), 
                    effect=rep(c('Mass', 'Treatment'),nFinalAnaly),
                    humValue=numeric(nFinalAnaly*2), humSE=numeric(nFinalAnaly*2), 
                    humP=numeric(nFinalAnaly*2), humSig=character(nFinalAnaly*2), 
                    femValue=numeric(nFinalAnaly*2), femSE=numeric(nFinalAnaly*2), 
                    femP=numeric(nFinalAnaly*2), femSig=character(nFinalAnaly*2),
                    stringsAsFactors=FALSE)
for (i in 1:nAnaly) {
  row = 1 + 2*(i-1)
  test = finalAnalyses[i]
  massVar = 'MassLog'
  if (any(test == c('BVTV', 'SMI', 'BATA'))) {
    massVar = 'Mass'
  }
  form = formula(paste(test, '~', massVar, '+Treat', sep=''))
  
  humData = humLog[,c('Mass','MassLog', 'Treat', 'Group', 'Alti', 'AltiLog', 'Day', test)] 
  humData = humData[complete.cases(humData),]
  if (any(test == c('BVTV', 'SMI', 'BATA'))) {
    humAnaly = summary(lme(form, random=list(~1|Alti, ~1|Day), data=humData))$tTable
  } else {
    humAnaly = summary(lme(form, random=list(~1|AltiLog, ~1|Day), data=humData))$tTable
  }
  values[c(row:(row+1)),c(3,4)] = humAnaly[c(2:3), c(1,2)]
  values[c(row:(row+1)),5] = humAnaly[c(2:3), 5]
  values[c(row:(row+1)),6] = psig(humAnaly[c(2:3), 5])
  
  femData = femLog[,c('Mass','MassLog', 'Treat', 'Group', 'Alti', 'AltiLog', 'Day', test)] 
  femData = femData[complete.cases(femData),]
  if (any(test == c('BVTV', 'SMI', 'BATA'))) {
    femAnaly = summary(lme(form, random=list(~1|Alti, ~1|Day), data=femData))$tTable
  } else {
    femAnaly = summary(lme(form, random=list(~1|AltiLog, ~1|Day), data=femData))$tTable
  }
  values[c(row:(row+1)),c(7,8)] = femAnaly[c(2:3), c(1,2)]
  values[c(row:(row+1)),9] = femAnaly[c(2:3), 5]
  values[c(row:(row+1)),10] = psig(femAnaly[c(2:3), 5])
}

# Save to use 
write.csv(values, paste(dir.export.data, '/effect_sizes.csv', sep=''))

#### PLOTS: LENGTH ####

# Get x ticks (mass) for all plots
x.ticks = seq(min(hum$Mass, na.rm=T), max(hum$Mass, na.rm=T), length.out=3)
x.limits = c(x.ticks[1], tail(x.ticks, 1))

# Plot. Change y.val, hum vs. fem, formula Mass vs. MassLog, y= in ggplot, annote 'A', y axis label, trans=log_trans
y.val = 'Length'
y.values=c(hum[,y.val]); y.ticks=seq(min(y.values), max(y.values), length.out=4); y.limits=c(y.ticks[1], tail(y.ticks,1))
a = humLog[,c('Mass','MassLog', 'Group', 'Treat', 'Alti', 'Day', y.val)]; a = a[complete.cases(a),]
b = lme(formula(paste(y.val, '~MassLog+Treat')), random=list(~1|Alti, ~1|Day), data=a)$coefficients$fixed
int.res=b[1]; int.mig=int.res + b[3]; slope=b[2]
d = data.frame(s=c(slope), i=c(int.res))
h.len = ggplot(hum, aes(x=Mass, y=Length, fill=Group, shape=Group)) + geom_point(size=2, color='transparent', stroke=1.2) +
  scale_shape_manual(values=c(21, 21, 24, 24, 24), labels=subspp_labels) +
  scale_fill_manual(values=c("#E89047", "#8E4204", "#585EB1", "#090C54", "#C0C4F7"), labels=subspp_labels) +
  annotate("text", x=x.limits[1], y=y.limits[2], label="A", size=7, hjust=0, vjust=1) +
  scale_x_continuous(trans=log_trans(), labels=scaleFUN1, breaks=x.ticks, limits=x.limits) + 
  scale_y_continuous(trans=log_trans(), labels=scaleFUN1, breaks=y.ticks, limits=y.limits) + 
  xlab('Body mass (g)') + ylab('Humerus length (mm)') + mytheme +
  geom_abline(data=d, aes(slope=s, intercept=i, linetype=factor(i)), size=1.2) +
  scale_linetype_manual(values=c('solid'), labels=c('All Birds')) + 
  theme(legend.position=c(0.15, 1), legend.title=element_blank(), legend.text=element_text(size=8, face='italic'), 
        legend.key.size=unit(0.01, 'in'),
        legend.background=element_rect(linetype='solid', color='black', fill='white'), legend.box='horizontal',
        legend.justification=c(0,1), legend.spacing.y=unit(0,'in'), legend.box.margin=margin(0,0,0,0,'in'), 
        legend.box.spacing=unit(0, 'in'), legend.margin=margin(0.02, 0.02,0.02, 0.02,'in')) +
  guides(linetype=guide_legend(label.theme=element_text(size=8, face='plain')))

# Plot. Change y.val, hum vs. fem, formula Mass vs. MassLog, y= in ggplot, annote 'A', y axis label, trans=log_trans
y.val = 'Length'
y.values=c(fem[,y.val]); y.ticks=seq(min(y.values), max(y.values), length.out=4); y.limits=c(y.ticks[1], tail(y.ticks,1))
a = femLog[,c('Mass','MassLog', 'Group', 'Treat', 'Alti', 'Day', y.val)]; a = a[complete.cases(a),]
b = lme(formula(paste(y.val, '~MassLog+Treat')), random=list(~1|Alti, ~1|Day),data=a)$coefficients$fixed
int.res=b[1]; int.mig=int.res + b[3]; slope=b[2]
d = data.frame(s=c(slope, slope), i=c(int.res, int.mig))
f.len = ggplot(fem, aes(x=Mass, y=Length, fill=Group, shape=Group)) + geom_point(color='transparent', size=2, stroke=1.2) +
  scale_shape_manual(values=c(21, 21, 24, 24, 24), labels=subspp_labels) +
  scale_fill_manual(values=c("#E89047", "#8E4204", "#585EB1", "#090C54", "#C0C4F7")) +
  annotate("text", x=x.limits[1], y=y.limits[2], label="B", size=7, hjust=0, vjust=1) +
  scale_x_continuous(trans=log_trans(), labels=scaleFUN1, breaks=x.ticks, limits=x.limits) + 
  scale_y_continuous(trans=log_trans(), labels=scaleFUN1, breaks=y.ticks, limits=y.limits) + 
  xlab('Body mass (g)') + ylab('Femur length (mm)') + mytheme +
  geom_abline(data=d, aes(slope=s, intercept=i, color=factor(c(d$i[1], d$i[2]), levels=c(d$i[1], d$i[2]))), size=1.2) +
  scale_color_manual(values=c("#E89047", "#090C54"), labels=c('Resident', 'Migrant')) +
  theme(legend.position=c(0.15, 1), legend.title=element_blank(), legend.text=element_text(size=8), legend.key.size=unit(0.01, 'in'),
        legend.background=element_rect(linetype='solid', color='black', fill='white'), legend.box='horizontal',
        legend.justification=c(0,1), legend.spacing.y=unit(0,'in'), legend.box.margin=margin(0,0,0,0,'in'), 
        legend.box.spacing=unit(0, 'in'), legend.margin=margin(0.02, 0.02,0.02, 0.02,'in')) +
  guides(fill=FALSE, shape=FALSE)

# Save plots. To fit, cortical plots must be 3.5" wide (the Auk)
# (11-2-0.5)/4 = 2.125 tall
png(paste(dir.export.figs, '/length.png', sep=''), width=3.5*2, height=2.125*1, res=600, units='in')
grid.arrange(h.len, f.len, nrow=1)
dev.off()

#### PLOTS: TRABECULAR ####

# Get x ticks (mass) for all plots
x.ticks = seq(min(hum$Mass, na.rm=T), max(hum$Mass, na.rm=T), length.out=3)
x.limits = c(x.ticks[1], tail(x.ticks, 1))

# Plot. Change y.val, hum vs. fem, formula Mass vs. MassLog, y= in ggplot, annote 'A', y axis label, trans=log_trans
y.val = 'TbBV'
y.values=c(hum[,y.val]); y.ticks=seq(min(y.values), 0.30, length.out=3); y.limits=c(y.ticks[1], 0.30) # 
a = humLog[,c('Mass','MassLog', 'Treat', 'Alti', 'Day', y.val)]; a = a[complete.cases(a),]
b = lme(formula(paste(y.val, '~MassLog+Treat')), random=list(~1|Alti, ~1|Day), data=a)$coefficients$fixed
int.res=b[1]; int.mig=int.res + b[3]; slope=b[2]
d = data.frame(s=c(slope), i=c(int.res))
h.tb.bv = ggplot(hum, aes(x=Mass, y=TbBV, fill=Treat, shape=Treat)) + geom_point(size=2, color='transparent', stroke=1.2) +
  scale_shape_manual(values=c(21, 24), labels=c('Resident', 'Migrant')) +
  scale_fill_manual(values=c("#E89047", "#090C54"), labels=c('Resident', 'Migrant')) +
  annotate("text", x=x.limits[1], y=y.limits[2], label="A", size=7, hjust=0, vjust=1) +
  scale_x_continuous(trans=log_trans(), labels=scaleFUN1, breaks=x.ticks, limits=x.limits) + 
  scale_y_continuous(trans=log_trans(), labels=scaleFUN2, breaks=y.ticks, limits=y.limits) + 
  xlab('Body mass (g)') + ylab(bquote('Bone volume ('*~mm^3*')')) + mytheme +
  geom_abline(data=d, aes(intercept=i, slope=s, color=factor(i)), size=1.2) + 
  scale_color_manual(values=c('black'), labels=c('All Birds')) + 
  theme(legend.position=c(0.15 , 1), legend.title=element_blank(), legend.text=element_text(size=8), legend.key.size=unit(0.1, 'in'),
        legend.background=element_rect(linetype='solid', color='black', fill='white'), legend.justification=c(0,1), legend.spacing.y=unit(0,'in'),
        legend.box.margin=margin(0,0,0,0,'in'), legend.box.spacing=unit(0, 'in'), legend.margin=margin(0.02, 0.02,0.02, 0.02,'in'),
        legend.box='horizontal') + guides(fill=guide_legend(order=1), shape=guide_legend(order=1), color=guide_legend(order=2))

# Plot. Change y.val, hum vs. fem, formula Mass vs. MassLog, y= in ggplot, annote 'A', y axis label, trans=log_trans
y.val = 'TbTh'
y.values=c(hum[,y.val]); y.ticks=seq(min(y.values), max(y.values), length.out=4); y.limits=c(y.ticks[1], tail(y.ticks,1))
a = humLog[,c('Mass','MassLog', 'Treat', 'Alti', 'Day', y.val)]; a = a[complete.cases(a),]
b = lme(formula(paste(y.val, '~MassLog+Treat')), random=list(~1|Alti, ~1|Day), data=a)$coefficients$fixed
int.res=b[1]; int.mig=int.res + b[3]; slope=b[2]
h.tb.th = ggplot(hum, aes(x=Mass, y=TbTh, fill=Treat, shape=Treat)) +  geom_point(color='transparent', size=2, stroke=1.2) +
  scale_shape_manual(values=c(21, 24), labels=c('Resident', 'Migrant')) +
  scale_fill_manual(values=c("#E89047", "#090C54"), labels=c('Resident', 'Migrant')) + 
  annotate("text", x=x.limits[1], y=y.limits[2], label="B", size=7, hjust=0, vjust=1) +
  scale_x_continuous(trans=log_trans(), labels=scaleFUN1, breaks=x.ticks, limits=x.limits) + 
  scale_y_continuous(trans=log_trans(), labels=scaleFUN2, breaks=y.ticks, limits=y.limits) + 
  xlab('Body mass (g)') + ylab('Thickness (mm)') + mytheme

# Plot. Change y.val, hum vs. fem, formula Mass vs. MassLog, y= in ggplot, annote 'A', y axis label, trans=log_trans
y.val = 'TbBV'
y.values=c(fem[,y.val]); y.ticks=seq(min(y.values), max(y.values), length.out=3); y.limits=c(y.ticks[1], tail(y.ticks,1))
a = humLog[,c('Mass','MassLog', 'Treat', 'Alti', 'Day', y.val)]; a = a[complete.cases(a),]
b = lme(formula(paste(y.val, '~MassLog+Treat')), random=list(~1|Alti, ~1|Day), data=a)$coefficients$fixed
int.mig=int.res + b[3]; slope=b[2]
d = data.frame(s=c(slope), i=c(int.res))
f.tb.bv = ggplot(fem, aes(x=Mass, y=TbBV, fill=Treat, shape=Treat)) + geom_point(color='transparent', size=2, stroke=1.2) +
  scale_shape_manual(values=c(21, 24), labels=c('Resident', 'Migrant')) +
  scale_fill_manual(values=c("#E89047", "#090C54"), labels=c('Resident', 'Migrant')) + 
  annotate("text", x=x.limits[1], y=y.limits[2], label="C", size=7, hjust=0, vjust=1) +
  scale_x_continuous(trans=log_trans(), labels=scaleFUN1, breaks=x.ticks, limits=x.limits) + 
  scale_y_continuous(trans=log_trans(), labels=scaleFUN2, breaks=y.ticks, limits=y.limits) + 
  xlab('Body mass (g)') + ylab(bquote('Bone volume ('*~mm^3*')')) + mytheme +
  geom_abline(data=d, aes(intercept=i, slope=s, color=factor(i)), size=1.2) + 
  scale_color_manual(values=c('black'), labels=c('All Birds'))

# Plot. Change y.val, hum vs. fem, formula Mass vs. MassLog, y= in ggplot, annote 'A', y axis label, trans=log_trans
y.val = 'TbTh'
y.values=c(fem[,y.val]); y.ticks=seq(min(y.values), 0.05, length.out=4); y.limits=c(y.ticks[1], 0.05)
a = femLog[,c('Mass','MassLog', 'Treat', y.val)]; a = a[complete.cases(a),]
b = gls(formula(paste(y.val, '~MassLog+Treat')), data=a)$coef; int.res=b[1]; int.mig=int.res + b[3]; slope=b[2]
d = data.frame(s=c(slope, slope), i=c(int.res, int.mig))
f.tb.th = ggplot(fem, aes(x=Mass, y=TbTh, fill=Treat, shape=Treat)) + geom_point(color='transparent', size=2, stroke=1.2) +
  scale_shape_manual(values=c(21, 24), labels=c('Resident', 'Migrant')) +
  scale_fill_manual(values=c("#E89047", "#090C54"), labels=c('Resident', 'Migrant')) +
  annotate("text", x=x.limits[1], y=y.limits[2], label="D", size=7, hjust=0, vjust=1) +
  scale_x_continuous(trans=log_trans(), labels=scaleFUN1, breaks=x.ticks, limits=x.limits) + 
  scale_y_continuous(trans=log_trans(), labels=scaleFUN2, breaks=y.ticks, limits=y.limits) + 
  xlab('Body mass (g)') + ylab('Thickness (mm)') + mytheme +
  geom_abline(data=d, aes(slope=s, intercept=i, color=factor(c(d$i[1], d$i[2]), levels=c(d$i[1], d$i[2]))), size=1.2) + 
  scale_color_manual(values=c("#E89047", "#090C54"), labels=c('Resident', 'Migrant')) +
  theme(legend.position=c(0.15, 1), legend.title=element_blank(), legend.text=element_text(size=8), legend.key.size=unit(0.1, 'in'),
        legend.background=element_rect(linetype='solid', color='black', fill='white'), legend.justification=c(0,1),
        legend.box.margin=margin(0,0,0,0,'in'), legend.box.spacing=unit(0, 'in'), legend.margin=margin(0.02, 0.02,0.02, 0.02,'in')) +
  guides(fill=FALSE, shape=FALSE)

# Save plots. To fit, cortical plots must be (8.5-2)/2 = 3.25 wide
# (11-2-0.5)/4 = 2.125 tall
png(paste(dir.export.figs, '/trabecular.png', sep=''), width=3.5*2, height=2.125*2, res=600, units='in')
grid.arrange(h.tb.bv, f.tb.bv, h.tb.th, f.tb.th, nrow=2)
dev.off()

#### PLOTS: CORTICAL ####

# Get x ticks (mass) for all plots
x.ticks = seq(min(hum$Mass, na.rm=T), max(hum$Mass, na.rm=T), length.out=3)
x.limits = c(x.ticks[1], tail(x.ticks, 1))

# Plot. Add hjust=0, vjust=1 to annotate, change legend.position=c(0.15, 1), add the following:
# legend.box='horizontal', legend.justification=c(0,1), legend.spacing.y=unit(0,'in'), legend.box.margin=margin(0,0,0,0,'in'), 
# legend.box.spacing=unit(0, 'in'), legend.margin=margin(0.02, 0.02,0.02, 0.02,'in')) +
# guides(fill=FALSE, shape=FALSE)
y.val = 'CtBV'
y.values=c(hum[,y.val]); y.ticks=seq(min(y.values), max(y.values), length.out=3); y.limits=c(y.ticks[1], 0.632)
a = humLog[,c('Mass','MassLog', 'Treat', y.val)]; a = a[complete.cases(a),]
b = gls(formula(paste(y.val, '~MassLog+Treat')), data=a)$coef; int.res=b[1]; int.mig=int.res + b[3]; slope=b[2]
d = data.frame(s=c(slope), i=c(int.res))
h.ct.bv = ggplot(hum, aes(x=Mass, y=CtBV, fill=Treat, shape=Treat)) + geom_point(color='transparent', size=2, stroke=1.2) +
  scale_shape_manual(values=c(21, 24), labels=c('Resident', 'Migrant')) +
  scale_fill_manual(values=c("#E89047", "#090C54"), labels=c('Resident', 'Migrant')) +
  annotate("text", x=x.limits[1], y=y.limits[2], label="A", size=7, hjust=0, vjust=1) +
  scale_x_continuous(trans=log_trans(), labels=scaleFUN1, breaks=x.ticks, limits=x.limits) + 
  scale_y_continuous(trans=log_trans(), labels=scaleFUN3, breaks=y.ticks, limits=y.limits) + 
  xlab('Body mass (g)') + ylab(bquote('Bone volume ('*~mm^3*')')) + mytheme +
  geom_abline(data=d, aes(intercept=i, slope=s, color=c('black')), size=1.2) + 
  scale_color_manual(values=c('black'), labels=c('All Birds')) +
  theme(legend.position=c(0.15, 1), legend.title=element_blank(), legend.text=element_text(size=8), legend.key.size=unit(0.1, 'in'),
        legend.background=element_rect(linetype='solid', color='black', fill='white'),
        legend.box='horizontal', legend.justification=c(0,1), legend.spacing.y=unit(0,'in'), legend.box.margin=margin(0,0,0,0,'in'), 
        legend.box.spacing=unit(0, 'in'), legend.margin=margin(0.02, 0.02,0.02, 0.02,'in')) +
  guides(fill=guide_legend(order=1), shape=guide_legend(order=1), color=guide_legend(order=2))

# Plot. Add hjust=0, vjust=1 to annotate, change legend.position=c(0.15, 1), add the following:
# legend.box='horizontal', legend.justification=c(0,1), legend.spacing.y=unit(0,'in'), legend.box.margin=margin(0,0,0,0,'in'), 
# legend.box.spacing=unit(0, 'in'), legend.margin=margin(0.02, 0.02,0.02, 0.02,'in')) +
# guides(fill=FALSE, shape=FALSE)
y.val = 'CtTh'
y.values=c(hum[,y.val]); y.ticks=seq(min(y.values), max(y.values), length.out=3); y.limits=c(y.ticks[1], tail(y.ticks,1))
a = humLog[,c('Mass','MassLog', 'Treat', y.val)]; a = a[complete.cases(a),]
b = gls(formula(paste(y.val, '~MassLog+Treat')), data=a)$coef; int.res=b[1]; int.mig=int.res + b[3]; slope=b[2]
d = data.frame(s=c(slope, slope), i=c(int.res, int.mig))
h.ct.th = ggplot(hum, aes(x=Mass, y=CtTh, fill=Treat, shape=Treat)) + geom_point(color='transparent', size=2, stroke=1.2) +
  scale_shape_manual(values=c(21, 24), labels=c('Resident', 'Migrant')) +
  scale_fill_manual(values=c("#E89047", "#090C54"), labels=c('Resident', 'Migrant')) +
  annotate("text", x=x.limits[1], y=y.limits[2], label="B", size=7, hjust=0, vjust=1) +
  scale_x_continuous(trans=log_trans(), labels=scaleFUN1, breaks=x.ticks, limits=x.limits) + 
  scale_y_continuous(trans=log_trans(), labels=scaleFUN3, breaks=y.ticks, limits=y.limits) + 
  xlab('Body mass (g)') + ylab('Cortical thickness (mm)') + mytheme +
  geom_abline(data=d, aes(slope=s, intercept=i, color=factor(c(d$i[1], d$i[2]), levels=c(d$i[1], d$i[2]))), size=1.2) + 
  scale_color_manual(values=c("#E89047", "#090C54"), labels=c('Resident', 'Migrant')) +
  guides(fill=FALSE, shape=FALSE, color=FALSE)

# Plot. Change y.val, hum vs. fem, formula Mass vs. MassLog, y= in ggplot, annote 'A', y axis label, trans=log_trans
y.val = 'do'
y.values=c(hum[,y.val]); y.ticks=seq(min(y.values), max(y.values), length.out=3); y.limits=c(y.ticks[1], tail(y.ticks,1))
a = humLog[,c('Mass','MassLog', 'Treat', y.val)]; a = a[complete.cases(a),]
b = gls(formula(paste(y.val, '~MassLog+Treat')), data=a)$coef; int.res=b[1]; int.mig=int.res + b[3]; slope=b[2]
d = data.frame(s=c(slope, slope), i=c(int.res, int.mig))
h.ct.do = ggplot(hum, aes(x=Mass, y=do, fill=Treat, shape=Treat)) + geom_point(color='transparent', size=2, stroke=1.2) +
  scale_shape_manual(values=c(21, 24), labels=c('Resident', 'Migrant')) +
  scale_fill_manual(values=c("#E89047", "#090C54"), labels=c('Resident', 'Migrant')) + 
  annotate("text", x=x.limits[1], y=y.limits[2], label="C", size=7, hjust=0, vjust=1) +
  scale_x_continuous(trans=log_trans(), labels=scaleFUN1, breaks=x.ticks, limits=x.limits) + 
  scale_y_continuous(trans=log_trans(), labels=scaleFUN3, breaks=y.ticks, limits=y.limits) + 
  xlab('Body mass (g)') + ylab('Cortical diameter (mm)') + mytheme +
  geom_abline(data=d, aes(slope=s, intercept=i, color=factor(c(d$i[1], d$i[2]), levels=c(d$i[1], d$i[2]))), size=1.2) + 
  scale_color_manual(values=c("#E89047", "#090C54"), labels=c('Resident', 'Migrant')) +
  theme(legend.position=c(0.05, 0.7), legend.title=element_blank(), legend.text=element_text(size=8), legend.key.size=unit(0.1, 'in'),
      legend.background=element_rect(linetype='solid', color='black', fill='white')) +
  guides(fill=FALSE, shape=FALSE, color=FALSE)

# Plot. Add hjust=0, vjust=1 to annotate, change legend.position=c(0.15, 1), add the following:
# legend.box='horizontal', legend.justification=c(0,1), legend.spacing.y=unit(0,'in'), legend.box.margin=margin(0,0,0,0,'in'), 
# legend.box.spacing=unit(0, 'in'), legend.margin=margin(0.02, 0.02,0.02, 0.02,'in')) +
# guides(fill=FALSE, shape=FALSE)
y.val = 'pMOI'
y.values=c(hum[,y.val]); y.ticks=seq(min(y.values), max(y.values), length.out=3); y.limits=c(y.ticks[1], tail(y.ticks,1))
a = humLog[,c('Mass','MassLog', 'Treat', y.val)]; a = a[complete.cases(a),]
b = gls(formula(paste(y.val, '~MassLog+Treat')), data=a)$coef; int.res=b[1]; int.mig=int.res + b[3]; slope=b[2]
d = data.frame(s=c(slope), i=c(int.res))
h.ct.pMOI = ggplot(hum, aes(x=Mass, y=pMOI, fill=Treat, shape=Treat)) + geom_point(color='transparent', size=2, stroke=1.2) +
  scale_shape_manual(values=c(21, 24), labels=c('Resident', 'Migrant')) +
  scale_fill_manual(values=c("#E89047", "#090C54"), labels=c('Resident', 'Migrant')) +
  annotate("text", x=x.limits[1], y=y.limits[2], label="D", size=7, hjust=0, vjust=1) +
  scale_x_continuous(trans=log_trans(), labels=scaleFUN1, breaks=x.ticks, limits=x.limits) + 
  scale_y_continuous(trans=log_trans(), labels=scaleFUN3, breaks=y.ticks, limits=y.limits) + 
  xlab('Body mass (g)') + ylab(bquote('Polar moment of inertia ('*~mm^4*')')) + mytheme +
  geom_abline(data=d, aes(intercept=i, slope=s, color=c('black')), size=1.2) + 
  scale_color_manual(values=c('black'), labels=c('All Birds')) +
  guides(fill=FALSE, shape=FALSE, color=FALSE)

# Plot. Add hjust=0, vjust=1 to annotate, change legend.position=c(0.15, 1), add the following:
# legend.box='horizontal', legend.justification=c(0,1), legend.spacing.y=unit(0,'in'), legend.box.margin=margin(0,0,0,0,'in'), 
# legend.box.spacing=unit(0, 'in'), legend.margin=margin(0.02, 0.02,0.02, 0.02,'in')) +
# guides(fill=FALSE, shape=FALSE)
y.val = 'CtBV'
y.values=c(fem[,y.val]); y.ticks=seq(min(y.values), max(y.values), length.out=3); y.limits=c(y.ticks[1], 0.293)
a = femLog[,c('Mass','MassLog', 'Treat', y.val)]; a = a[complete.cases(a),]
b = gls(formula(paste(y.val, '~MassLog+Treat')), data=a)$coef; int.res=b[1]; int.mig=int.res + b[3]; slope=b[2]
d = data.frame(s=c(slope, slope), i=c(int.res, int.mig))
f.ct.bv = ggplot(fem, aes(x=Mass, y=CtBV, fill=Treat, shape=Treat)) + geom_point(color='transparent', size=2, stroke=1.2) +
  scale_shape_manual(values=c(21, 24), labels=c('Resident', 'Migrant')) +
  scale_fill_manual(values=c("#E89047", "#090C54"), labels=c('Resident', 'Migrant')) +
  annotate("text", x=x.limits[1], y=y.limits[2], label="E", size=7, hjust=0, vjust=1) +
  scale_x_continuous(trans=log_trans(), labels=scaleFUN1, breaks=x.ticks, limits=x.limits) + 
  scale_y_continuous(trans=log_trans(), labels=scaleFUN3, breaks=y.ticks, limits=y.limits) + 
  xlab('Body mass (g)') + ylab(bquote('Bone volume ('*~mm^3*')')) + mytheme +
  geom_abline(data=d, aes(slope=s, intercept=i, color=factor(c(d$i[1], d$i[2]), levels=c(d$i[1], d$i[2]))), size=1.2) + 
  scale_color_manual(values=c("#E89047", "#090C54"), labels=c('Resident', 'Migrant')) +
  theme(legend.position=c(0.15, 1), legend.title=element_blank(), legend.text=element_text(size=8), legend.key.size=unit(0.1, 'in'),
      legend.background=element_rect(linetype='solid', color='black', fill='white'),
      legend.box='horizontal', legend.justification=c(0,1), legend.spacing.y=unit(0,'in'), legend.box.margin=margin(0,0,0,0,'in'), 
      legend.box.spacing=unit(0, 'in'), legend.margin=margin(0.02, 0.02,0.02, 0.02,'in')) +
  guides(fill=FALSE, shape=FALSE)

# Plot. Change y.val, hum vs. fem, formula Mass vs. MassLog, y= in ggplot, annote 'A', y axis label, trans=log_trans
y.val = 'CtTh'
y.values=c(fem[,y.val]); y.ticks=seq(min(y.values), max(y.values), length.out=3); y.limits=c(y.ticks[1], tail(y.ticks,1))
a = femLog[,c('Mass','MassLog', 'Treat', y.val)]; a = a[complete.cases(a),]
b = gls(formula(paste(y.val, '~MassLog+Treat')), data=a)$coef; int.res=b[1]; int.mig=int.res + b[3]; slope=b[2]
d = data.frame(s=c(slope, slope), i=c(int.res, int.mig))
f.ct.th = ggplot(fem, aes(x=Mass, y=CtTh, fill=Treat, shape=Treat)) + geom_point(color='transparent', size=2, stroke=1.2) +
  scale_shape_manual(values=c(21, 24), labels=c('Resident', 'Migrant')) +
  scale_fill_manual(values=c("#E89047", "#090C54"), labels=c('Resident', 'Migrant')) +
  annotate("text", x=x.limits[1], y=y.limits[2], label="F", size=7, hjust=0, vjust=1) +
  scale_x_continuous(trans=log_trans(), labels=scaleFUN1, breaks=x.ticks, limits=x.limits) + 
  scale_y_continuous(trans=log_trans(), labels=scaleFUN3, breaks=y.ticks, limits=y.limits) + 
  xlab('Body mass (g)') + ylab('Cortical thickness (mm)') + mytheme +
  geom_abline(data=d, aes(slope=s, intercept=i, color=factor(c(d$i[1], d$i[2]), levels=c(d$i[1], d$i[2]))), size=1.2) + 
  scale_color_manual(values=c("#E89047", "#090C54"), labels=c('Resident', 'Migrant')) +
  theme(legend.position=c(0.05, 0.7), legend.title=element_blank(), legend.text=element_text(size=8), legend.key.size=unit(0.1, 'in'),
        legend.background=element_rect(linetype='solid', color='black', fill='white')) +
  guides(color=FALSE, fill=FALSE, shape=FALSE)

# Plot. Change y.val, hum vs. fem, formula Mass vs. MassLog, y= in ggplot, annote 'A', y axis label, trans=log_trans
y.val = 'do'
y.values=c(fem[,y.val]); y.ticks=seq(min(y.values), max(y.values), length.out=3); y.limits=c(y.ticks[1], tail(y.ticks,1))
a = femLog[,c('Mass','MassLog', 'Treat', y.val)]; a = a[complete.cases(a),]
b = gls(formula(paste(y.val, '~MassLog+Treat')), data=a)$coef; int.res=b[1]; int.mig=int.res + b[3]; slope=b[2]
d = data.frame(s=c(slope), i=c(int.res))
f.ct.do = ggplot(fem, aes(x=Mass, y=do, fill=Treat, shape=Treat)) + geom_point(color='transparent', size=2, stroke=1.2) +
  scale_shape_manual(values=c(21, 24), labels=c('Resident', 'Migrant')) +
  scale_fill_manual(values=c("#E89047", "#090C54"), labels=c('Resident', 'Migrant')) + 
  annotate("text", x=x.limits[1], y=y.limits[2], label="G", size=7, hjust=0, vjust=1) +
  scale_x_continuous(trans=log_trans(), labels=scaleFUN1, breaks=x.ticks, limits=x.limits) + 
  scale_y_continuous(trans=log_trans(), labels=scaleFUN3, breaks=y.ticks, limits=y.limits) + 
  xlab('Body mass (g)') + ylab('Cortical diameter (mm)') + mytheme +
  geom_abline(data=d, aes(intercept=i, slope=s, color=c('black')), size=1.2) + 
  scale_color_manual(values=c('black'), labels=c('All Birds')) +
  theme(legend.position=c(0.05, 0.7), legend.title=element_blank(), legend.text=element_text(size=8), legend.key.size=unit(0.1, 'in'),
        legend.background=element_rect(linetype='solid', color='black', fill='white')) +
  guides(color=FALSE, fill=FALSE, shape=FALSE)

# Plot. Change y.val, hum vs. fem, formula Mass vs. MassLog, y= in ggplot, annote 'A', y axis label, trans=log_trans
y.val = 'pMOI'
y.values=c(fem[,y.val]); y.ticks=seq(min(y.values), max(y.values), length.out=3); y.limits=c(y.ticks[1], tail(y.ticks,1))
a = femLog[,c('Mass','MassLog', 'Treat', y.val)]; a = a[complete.cases(a),]
b = gls(formula(paste(y.val, '~MassLog+Treat')), data=a)$coef; int.res=b[1]; int.mig=int.res + b[3]; slope=b[2]
d = data.frame(s=c(slope, slope), i=c(int.res, int.mig))
f.ct.pMOI = ggplot(fem, aes(x=Mass, y=pMOI, fill=Treat, shape=Treat)) + geom_point(color='transparent', size=2, stroke=1.2) +
  scale_shape_manual(values=c(21, 24), labels=c('Resident', 'Migrant')) +
  scale_fill_manual(values=c("#E89047", "#090C54"), labels=c('Resident', 'Migrant')) +
  annotate("text", x=x.limits[1], y=y.limits[2], label="H", size=7, hjust=0, vjust=1) +
  scale_x_continuous(trans=log_trans(), labels=scaleFUN1, breaks=x.ticks, limits=x.limits) + 
  scale_y_continuous(trans=log_trans(), labels=scaleFUN3, breaks=y.ticks, limits=y.limits) + 
  xlab('Body mass (g)') + ylab(bquote('Polar moment of inertia ('*~mm^4*')')) + mytheme +
  geom_abline(data=d, aes(slope=s, intercept=i, color=factor(c(d$i[1], d$i[2]), levels=c(d$i[1], d$i[2]))), size=1.2) + 
  scale_color_manual(values=c("#E89047", "#090C54"), labels=c('Resident', 'Migrant')) +
  theme(legend.position=c(0.05, 0.7), legend.title=element_blank(), legend.text=element_text(size=8), legend.key.size=unit(0.1, 'in'),
        legend.background=element_rect(linetype='solid', color='black', fill='white')) +
  guides(color=FALSE, fill=FALSE, shape=FALSE)

# Save plots. To fit, cortical plots must be (8.5-2)/2 = 3.25 wide
# (11-2-0.5)/4 = 2.125 tall
png(paste(dir.export.figs, '/cortical.png', sep=''), width=3.5*2, height=2*4, res=600, units='in')
grid.arrange(h.ct.bv, f.ct.bv, h.ct.th, f.ct.th, h.ct.do, f.ct.do, h.ct.pMOI, f.ct.pMOI, nrow=4)
dev.off()

#### PLOTS: SUPPLEMENTAL ####

# Plot. Change y.val, hum vs. fem, formula Mass vs. MassLog, y= in ggplot, annote 'A', y axis label, trans=log_trans
y.val = 'Alti'
y.ticks = c(0, 1000, 2000); y.limits=c(0, 2100)
a = hum[,c('Mass','Alti', 'Group', 'Treat')]; a = a[complete.cases(a),]
alti = ggplot(a, aes(x=Mass, y=Alti, fill=Group, shape=Group)) + geom_point(size=2, color='transparent', stroke=1.2) +
  scale_shape_manual(values=c(21, 21, 24, 24, 24), labels=subspp_labels) +
  scale_fill_manual(values=c("#E89047", "#8E4204", "#585EB1", "#090C54", "#C0C4F7"), labels=subspp_labels) +
  xlab('Body mass (g)') + ylab('Elevation (m)') + mytheme +
  scale_y_continuous(breaks=y.ticks, limits=y.limits) + 
  theme(legend.position=c(0.65, 0.45), legend.title=element_blank(), legend.text=element_text(size=8, face='italic'), 
        legend.key.size=unit(0.01, 'in'),
        legend.background=element_rect(linetype='solid', color='black', fill='white'), legend.box='horizontal',
        legend.justification=c(0,1), legend.spacing.y=unit(0,'in'), legend.box.margin=margin(0,0,0,0,'in'), 
        legend.box.spacing=unit(0, 'in'), legend.margin=margin(0.02, 0.02,0.02, 0.02,'in')); alti
png(paste(dir.export.figs, '/elev.png', sep=''), width=3.5, height=2, res=600, units='in'); alti; dev.off()

# Plot. Change y.val, hum vs. fem, formula Mass vs. MassLog, y= in ggplot, annote 'A', y axis label, trans=log_trans
y.val = 'Day'
a = hum[,c('Mass','Day', 'Group', 'Treat')]; a = a[complete.cases(a),]
day = ggplot(a, aes(x=Mass, y=Day, fill=Group, shape=Group)) + geom_point(size=2, color='transparent', stroke=1.2) +
  scale_shape_manual(values=c(21, 21, 24, 24, 24), labels=subspp_labels) +
  scale_fill_manual(values=c("#E89047", "#8E4204", "#585EB1", "#090C54", "#C0C4F7"), labels=subspp_labels) +
  xlab('Body mass (g)') + ylab('Day of the Year (1-365)') + mytheme +
  theme(legend.position=c(0.65, 0.8), legend.title=element_blank(), legend.text=element_text(size=8, face='italic'), legend.key.size=unit(0.01, 'in'),
        legend.background=element_rect(linetype='solid', color='black', fill='white'), legend.box='horizontal',
        legend.justification=c(0,1), legend.spacing.y=unit(0,'in'), legend.box.margin=margin(0,0,0,0,'in'), 
        legend.box.spacing=unit(0, 'in'), legend.margin=margin(0.02, 0.02,0.02, 0.02,'in')); day
png(paste(dir.export.figs, '/day.png', sep=''), width=3.5, height=2, res=600, units='in'); day; dev.off()
