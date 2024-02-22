# Analyses of SST and N-back behavioral data
# Written by Kevin Potter
# email: kevin.w.potter@gmail.com
# Please email me directly if you 
# have any questions or comments
# Last updated 2023-06-21

# Table of contents
# 1) Setup

#### 1) Setup ####

# Load in helpful packages
library(arfpam)
library(dplyr)
library(glmmTMB)

# Read in data
dtf_MMJ <- read.csv(
  "MMJ-Processed_data_for_Debbie_Burdinski-2021_11_08-18_12-a11d8a2.csv"
)
load( "MMJ-SST_summary_measures-2023_05_30-15_52.RData" )

chr_dir_working <- getwd()
setwd( 'MMJ-Nback_raw_data-2021_04_09-17_10' )

dtf_NBA <- readxl::read_excel(
  'nback_Accuracy_RTime.xlsx',
  sheet = 'data',
  col_names = TRUE
) %>% 
  data.frame()

setwd(chr_dir_working)

#### 2) Prep data ####

dtf_to_fit.SST <- SST %>% 
  filter(
    IDS.CHR.Subject %in% dtf_MMJ$IDS.CHR.Subject
  ) %>% 
  group_by(
    IDS.CHR.Subject,
    SSS.CHR.Time_point,
    IDX.INT.SST.Run
  ) %>% 
  summarise(
    SSS.CHR.Condition = 'Run ' %p% unique( IDX.INT.SST.Run ), 
    Y = unique( TSK.INT.SST.Stop.RT.Integration ),
    .groups = 'drop'
  ) %>% 
  data.frame() %>% 
  filter(
    Y > 0
  ) %>% 
  mutate(
    SSS.FCT.Time_point = factor(
      SSS.CHR.Time_point,
      levels = c( 'Baseline', 'One year' ), 
      labels = c( 'Baseline', 'One year' )
    ),
    PRDFollowHPHup = as.numeric(
      SSS.CHR.Time_point == 'One year'
    ),
    PRDSecondZZrun = as.numeric(
      IDX.INT.SST.Run %in% 2
    ), 
    PRDFollowHPHupZZbyZZsecondZZrun = 
      PRDFollowHPHup * PRDSecondZZrun, 
    Participant = factor(
      IDS.CHR.Subject
    )
  )

dtf_to_fit.NBA <- dtf_NBA %>% 
  filter(
    subject %in% dtf_MMJ$IDS.CHR.Subject
  ) %>% 
  group_by(
    IDS.CHR.Subject = subject,
    SSS.CHR.Time_point = timepoint
  ) %>% 
  summarise(
    Condition_0_back = unique( RT_0b_cor ),
    Condition_2_back = unique( RT_2b_cor ),
    .groups = 'drop'
  ) %>% 
  data.frame() %>% 
  tidyr::pivot_longer(
    cols = c( 'Condition_0_back', 'Condition_2_back' ),
    names_to = 'SSS.CHR.Condition',
    values_to = 'Y'
  ) %>% 
  data.frame() %>% 
  filter(
    !is.na( Y )
  ) %>% 
  mutate(
    SSS.CHR.Condition = recode(
      SSS.CHR.Condition,
      `Condition_0_back` = '0-back',
      `Condition_2_back` = '2-back',
    ),
    SSS.FCT.Time_point = factor(
      SSS.CHR.Time_point,
      levels = c( 'baseline', '1year' ), 
      labels = c( 'Baseline', 'One year' )
    ),
    PRDFollowHPHup = as.numeric(
      SSS.CHR.Time_point == '1year'
    ),
    PRD2HPHback = as.numeric(
      SSS.CHR.Condition %in% '2-back'
    ), 
    PRDFollowHPHupZZbyZZ2HPHback = 
      PRDFollowHPHup * PRD2HPHback,
    Participant = factor(
      IDS.CHR.Subject
    )
  )

#### 3) Analyses ####

fun_summarize_results <- function(
    obj_MLM,
    lgc_print = TRUE ) {
  
  lst_summary <- summary( obj_MLM )
  mat_coef <- lst_summary$coefficients$cond
  
  dtf_results <- data.frame(
    Effect = ffcar::ffcar_variables_to_labels(
      rownames(mat_coef)[-1]
    ),
    Estimate = '',
    PHPHvalue = ''
  )
  
  num_est <- format( round( mat_coef[-1, 1], 1 ), nsmall = 1 )
  num_lb <- format( round( 
      mat_coef[-1, 1] + qnorm(.025)*mat_coef[-1, 2], 
      1 ), nsmall = 1
    )
  num_ub <- format( round( 
    mat_coef[-1, 1] + qnorm(.975)*mat_coef[-1, 2], 
    1 ), nsmall = 1
  )
  
  dtf_results$Estimate <- 
    ffcar::ffcar_unicode_lookup()[['GLb']] %p% 
    ' = ' %p% 
    num_est %p% 
    ' (' %p% 
    num_lb %p% ' to ' %p% num_ub %p% ')'
  
  dtf_results$PHPHvalue <- 
    'p = ' %p% 
    format( round( mat_coef[-1, 4], 3), nsmall = 3 )
  dtf_results$PHPHvalue <- gsub(
    'p = 0.000', 'p < 0.001', 
    dtf_results$PHPHvalue, fixed = TRUE
  )
  
  if ( lgc_print ) {
    
    chr_header <- colnames(dtf_results)
    chr_header <- ffcar::ffcar_variables_to_labels(
      chr_header
    )
    
    int_chr <- apply(
      dtf_results, 2, function(x) max( nchar(x) )
    )
    int_chr <- sapply(
      seq_along(int_chr), function(i) {
        max( int_chr[i], nchar(chr_header) )
      }
    )
    
    for ( r in 1:nrow(dtf_results) ) {
      
      if ( r == 1 ) {
        
        chr_current <- squish(
          sapply( seq_along(chr_header), function(i) {
            stringr::str_pad( chr_header[i], width = int_chr[i] )
          } ),
          ' | '
        )
        cat( chr_current )
        cat('\n')
        
      }
      
      chr_current <- squish(
        sapply( seq_along(chr_header), function(i) {
          stringr::str_pad( dtf_results[r, i], width = int_chr[i] )
        } ),
        ' | '
      )
      cat( chr_current )
      cat('\n')
      
    }
    
  }
  
  return(dtf_results)
}

#### 3.1) SST ####

obj_MLM.SST.M0 <- glmmTMB(
  Y ~ 
    (1|Participant), 
  data = dtf_to_fit.SST
)

obj_MLM.SST.M1 <- glmmTMB(
  Y ~ 
    PRDFollowHPHup + 
    (1|Participant), 
  data = dtf_to_fit.SST
)

obj_MLM.SST.M2 <- glmmTMB(
  Y ~ 
    PRDSecondZZrun + 
    (1|Participant), 
  data = dtf_to_fit.SST
)

obj_MLM.SST.M3 <- glmmTMB(
  Y ~ 
    PRDFollowHPHup + 
    PRDSecondZZrun + 
    (1|Participant), 
  data = dtf_to_fit.SST
)

obj_MLM.SST.M4 <- glmmTMB(
  Y ~ 
    PRDFollowHPHup + 
    PRDSecondZZrun + 
    PRDFollowHPHupZZbyZZsecondZZrun + 
    (1|Participant), 
  data = dtf_to_fit.SST
)

dtf_comp.SST <- anova(
  obj_MLM.SST.M0,
  obj_MLM.SST.M1,
  obj_MLM.SST.M2,
  obj_MLM.SST.M3,
  obj_MLM.SST.M4
)
int_model <- which.min( dtf_comp.SST$BIC ) - 1
chr_best_model.SST <- paste0( 'obj_MLM.SST.M', int_model )

lst_summary.SST <- summary( get( chr_best_model.SST ) )

#### 3.2) N-back ####

obj_MLM.NBA.M0 <- glmmTMB(
  Y ~ 
    (1|Participant), 
  data = dtf_to_fit.NBA
)

obj_MLM.NBA.M1 <- glmmTMB(
  Y ~ 
    PRDFollowHPHup + 
    (1|Participant), 
  data = dtf_to_fit.NBA
)

obj_MLM.NBA.M2 <- glmmTMB(
  Y ~ 
    PRD2HPHback + 
    (1|Participant), 
  data = dtf_to_fit.NBA
)

obj_MLM.NBA.M3 <- glmmTMB(
  Y ~ 
    PRDFollowHPHup + 
    PRD2HPHback + 
    (1|Participant), 
  data = dtf_to_fit.NBA
)

obj_MLM.NBA.M4 <- glmmTMB(
  Y ~ 
    PRDFollowHPHup + 
    PRD2HPHback + 
    PRDFollowHPHupZZbyZZ2HPHback + 
    (1|Participant), 
  data = dtf_to_fit.NBA
)

dtf_comp.NBA <- anova(
  obj_MLM.NBA.M0,
  obj_MLM.NBA.M1,
  obj_MLM.NBA.M2,
  obj_MLM.NBA.M3,
  obj_MLM.NBA.M4
)
int_model <- which.min( dtf_comp.NBA$BIC ) - 1
chr_best_model.NBA <- paste0( 'obj_MLM.NBA.M', int_model )

lst_summary.NBA <- summary( get( chr_best_model.NBA ) )

#### 4) Figures ####

fun_plot_desc_pred <- function(
    dtf_to_fit,
    obj_MLM,
    chr_outcome,
    num_inc = 10, 
    lgc_model = TRUE, 
    lgc_new = TRUE ) {
  
  chr_predictors <- dtf_to_fit %>% column( 'PRD' )
  
  dtf_new <- dtf_to_fit %>%
    group_by(
      Time = SSS.FCT.Time_point,
      Condition = SSS.CHR.Condition
    ) %>%
    summarise_at(
      chr_predictors, mean
    ) %>% 
    data.frame
  dtf_new$Participant <- 'NEW'
  
  dtf_pred <- predict(
    obj_MLM, newdata = dtf_new,
    re.form = NA, 
    allow.new.levels = TRUE,
    se.fit = TRUE
  )
  dtf_pred <- data.frame( dtf_pred )
  
  dtf_obs <- dtf_to_fit %>% 
    group_by(
      Time = SSS.FCT.Time_point,
      Condition = SSS.CHR.Condition
    ) %>% 
    summarise(
      N = length(Y),
      M = mean(Y),
      SE = sd(Y)/sqrt(N),
      LB = M + SE*qt(.025, N - 1),
      UB = M + SE*qt(.975, N - 1),
      .groups = 'drop'
    )
  dtf_obs$col <- rep( palettes( index = 1:2 ), each = 2 )
  dtf_obs$pch <- rep( c( 21, 24 ), 2 )
  
  if (lgc_new) {
    x11( width = 5, height = 5)
  }
  
  xl <- c( .75, nrow(dtf_obs) + .25 )
  yl <- range( c( dtf_obs$LB, dtf_obs$UB ) )
  yl[1] <- floor(yl[1]/num_inc)*num_inc
  yl[2] <- ceiling(yl[2]/num_inc)*num_inc
  blank_plot(xl, yl)
  
  hv_line( h = seq( yl[1], yl[2], num_inc ), l = xl, lwd = 2, col = 'grey80' )
  
  error_bars(
    1:nrow(dtf_obs), lb = dtf_obs$LB, ub = dtf_obs$UB, lwd = 2,
    col = dtf_obs$col
  )
  points( 1:nrow(dtf_obs), dtf_obs$M,
          pch = dtf_obs$pch, cex = 1.35,
          bg = dtf_obs$col )
  
  if ( lgc_model ) {
    
    error_bars(
      1:nrow( dtf_obs ) + .2,
      lb = dtf_pred$fit + qnorm(.025)*dtf_pred$se.fit,
      ub = dtf_pred$fit + qnorm(.975)*dtf_pred$se.fit
    )
    points(
      1:nrow( dtf_obs ) + .2,
      dtf_pred$fit, pch = dtf_obs$pch, bg = 'white'
    )
    
  }
  
  hv_line( h = yl[1], l = xl, lwd = 2 )
  hv_line( v = xl[1], l = yl, lwd = 2 )
  
  add_axes(
    seq( yl[1], yl[2], num_inc ),
    line = -1.25, side = 2, cex = 1.15
  )
  mtext(
    chr_outcome,
    side = 2, line = 1.5, cex = 1.15
  )
  
  add_axes(
    c( 1.5, 3.5 ), c( 'Baseline', 'One year' ),
    line = -1.25, side = 1, cex = 1.15
  )
  
  legend(
    xl[1], yl[2] + diff(yl)*.2,
    unique( dtf_obs$Condition ),
    pch = c( 21, 24 ),
    pt.bg = 'black',
    bty = 'n',
    xpd = NA
  )
  
  if ( lgc_model ) {
    
    legend(
      xl[1] + diff(xl)*.5, yl[2] + diff(yl)*.2,
      c( 'Model predictions', '' ),
      pch = c( 21, 24 ),
      pt.bg = 'black',
      bty = 'n',
      xpd = NA
    )
    
  }
  
  
}


if ( FALSE ) {
  
  load( "MMJ-Processed_data-2022_05_31-11_00-eb9c29a.RData" )
  
  dtf_to_fit.SST <- SST %>% 
    filter(
      grepl( 'HC', IDS.CHR.Subject ) | 
        ( IDS.CHR.Subject %in% dtf_MMJ$IDS.CHR.Subject )
    ) %>% 
    group_by(
      IDS.CHR.Subject,
      SSS.CHR.Time_point,
      IDX.INT.SST.Run
    ) %>% 
    summarise(
      SSS.CHR.Condition = 'Run ' %p% unique( IDX.INT.SST.Run ), 
      Y = unique( TSK.INT.SST.Stop.RT.Integration ),
      .groups = 'drop'
    ) %>% 
    data.frame() %>% 
    filter(
      Y > 0
    ) %>% 
    mutate(
      SSS.FCT.Time_point = factor(
        SSS.CHR.Time_point,
        levels = c( 'Baseline', 'One year' ), 
        labels = c( 'Baseline', 'One year' )
      ),
      SSS.FCT.Group = factor(
        grepl( 'HC', IDS.CHR.Subject ),
        levels = c( TRUE, FALSE ),
        labels = c( 'HC', 'MM' )
      ), 
      PRDMMZZparticipants = as.numeric(
        !grepl( 'HC', IDS.CHR.Subject )
      ), 
      PRDHCZZparticipants = as.numeric(
        grepl( 'HC', IDS.CHR.Subject )
      ), 
      PRDFollowHPHup = as.numeric(
        SSS.CHR.Time_point == 'One year'
      ),
      PRDSecondZZrun = as.numeric(
        SSS.CHR.Condition %in% 'Run 2'
      ),
      PRDFollowHPHupZZbyZZsecondZZrun = 
        PRDFollowHPHup * PRDSecondZZrun, 
      PRDAge = NA, 
      Participant = factor(
        IDS.CHR.Subject
      )
    )
  
  for ( r in 1:nrow(dtf_to_fit.SST ) ) {
    
    lgc_match <- dat$IDS.CHR.Subject == dtf_to_fit.SST$IDS.CHR.Subject[r]
    
    dtf_to_fit.SST$PRDAge[r] <- 
      unique( strip_value( dat$SBJ.INT.Age[lgc_match] ) )[1]
    
  }
  dtf_to_fit.SST$PRDAge <- scale(
    dtf_to_fit.SST$PRDAge
  )[, 1]
  
  dtf_to_fit.SST <- dtf_to_fit.SST %>% 
    mutate(
      PRDMMZZparticipantsZZOPPBaselineCLP = 
        PRDMMZZparticipants * (1 - PRDFollowHPHup),
      PRDMMZZparticipantsZZOPPFollowHPHupCLP = 
        PRDMMZZparticipants * PRDFollowHPHup
    )
  
  dtf_to_fit.SST %>% 
    group_by(
      SSS.FCT.Group, 
      SSS.FCT.Time_point
    ) %>% 
    summarise(
      N = n_distinct( IDS.CHR.Subject ), 
      O = summa( Y, '[[M]] ([[SD]])', digits = 0 ),
      .groups = 'drop'
    )
  
  obj_MLM.SST.V1 <- glmmTMB(
    Y ~ PRDMMZZparticipantsZZOPPBaselineCLP + 
      PRDMMZZparticipantsZZOPPFollowHPHupCLP+ (1|Participant),
    data = dtf_to_fit.SST
  )
  obj_MLM.SST.V2 <- glmmTMB(
    Y ~ PRDHCZZparticipants + 
      PRDMMZZparticipantsZZOPPFollowHPHupCLP+ (1|Participant),
    data = dtf_to_fit.SST
  )
  
  obj_MLM.SST.Adjusted <- glmmTMB(
    Y ~ PRDMMZZparticipants + PRDFollowHPHup + PRDAge + (1|Participant),
    data = dtf_to_fit.SST
  )
  
  
  
  
}
