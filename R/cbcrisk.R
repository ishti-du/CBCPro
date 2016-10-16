cbcrisk <-
function(profile,start.age,pred.year=5, print.output=T)
{


abs_risk_cbc=function(profile,start.age,pred.year,h1star.seer,h2.seer)
{
  require(survival) 
  require(car)
  
  if(start.age <18|start.age >89)
    
  {stop("current age has to be within 18-89")}
  
  
  if (start.age+pred.year>89)
  {stop("current age + nyears.pred<89")}
  
  ### Preparing the data:
  data=data.frame(t(profile))
  ## Changing the variable names:
  colnames(data) <- c("age_1std","hormf","famhx_bc","lcis_atyph","brd_density","ER","bc_typ","cat_age1stb")
  
  ## Recoding the inputs:
  data$age_1std=recode(data$age_1std,"1='<30';2='30-40';else='40+'",as.factor.result=T)
  data$hormf=recode(data$hormf,"1='Yes';2='No';else='Unk'",as.factor.result=T)
  data$famhx_bc=recode(data$famhx_bc,"1='Yes';2='No';else='Unk'",as.factor.result=T)
  data$lcis_atyph=recode(data$lcis_atyph,"1='Yes';else='Unk'",as.factor.result=T)
  data$brd_density=recode(data$brd_density, "4='all_fat';3='scattered';2='heterog_dense';1='extrm_dense';else='Unk'",as.factor.result=T)
  data$ER=recode(data$ER,"1='Neg';2='Pos';else='Unk'",as.factor.result=T)
  data$bc_typ=recode(data$bc_typ,"1='Pure DCIS';2='Invasive_DCIS';else='Pure Invasive'",as.factor.result=T)
  data$cat_age1stb=recode(data$cat_age1stb,"1='<30/nulli';2='30-39';3='40+';else='Unk'",as.factor.result=T)
  
  data$grp=1010
  
  cbc_coeff <- c( 0.78496379 ,0.27152492,-0.05094703,-0.2511127, 0.18167327,0.44951127,0.45036138,0.70396452,0.53252559,0.42360416,0.39372522,0.12897412,-0.17737098,0.30621685, 0.50522831,0.26240513 ,
                  1.31195801,-0.02886171)
  
  names(cbc_coeff) <- c("age_1std<30", "age_1std30-40" ,"hormfUnk","hormfYes","famhx_bcUnk", "famhx_bcYes" ,"lcis_atyphYes","brd_densityextrm_dense","brd_densityheterog_dense","brd_densityscattered","brd_densityUnk","ERNeg",                   
                        "ERUnk","bc_typInvasive_DCIS","bc_typPure DCIS","cat_age1stb30-39", "cat_age1stb40+","cat_age1stbUnk" )
  
  #load("sysdata.rda", envir=environment())
  
  #model <- readRDS("cbc_relative_risk_model.rds")
  
  x1=model.matrix(fmod11.bcsc,data)
  RR=exp(x1%*%cbc_coeff)       ### Getting the RR#####  
  
  aa=c(18,30,35,40,45,50,55,60,65,70,75,80,85)
  bb=c(29,34,39,44,49,54,59,64,69,74,79,84,89)+1
  
  
  
  ##### Choosing the hazard rates:    
  
  if(h1star.seer==T){
    h1star=c(0.002642137,0.003639010,0.003480063,0.003355744,0.003450679,0.003212342,0.003505955,0.003559329,0.004015681,0.003960725,0.004364455,0.003938050,0.003553980)}
  else {h1star=c(0.005938242,0.004942339,0.004235975,0.003648188,0.003079530,0.003297890,0.003486713,0.003210002,0.004035436,0.00431467,0.003866516,
                 0.004041648,0.003070979)}
  
  if (h2.seer==T)
  {other.hazard=c(0.03465443,0.03191816,0.02610106,0.01507113,0.01299555,0.01317742,0.01413674,0.01558502,0.01776497,0.02320382,0.03280301,0.05083317,0.08456722)}
  else {other.hazard=c(0.01672640,0.02119907,0.01616915,0.01288945,0.01173593,0.01086192,0.01095524,0.01292321,0.01443535,0.02118252,0.02893834,0.04330877,0.05676192)}
  
  ar=c(0.7073222,0.7031257,0.6135916,0.5476586,0.4961982,0.4727390,0.4566321,0.4531747,0.4097558,0.4398585,0.4400653,0.4078442,0.4004540)
  # 
  cbc.hazard=h1star*(1-ar)
  
  
  
  st=findInterval(start.age,aa)
  
  if((start.age+pred.year) %in% aa){ 
    
    fn= (findInterval(start.age+pred.year,aa))-1
  } else {fn=findInterval(start.age+pred.year,aa)}
  
  
  surv.cbc=rep(0,length=fn+1)     ##### Survival function vector for cbc
  surv.other=rep(0,length=fn+1)     #####  Survival function vector for non-cbc
  surv.cbc[1]=1
  surv.other[1]=1
  surv.cbc_a=0
  surv.other_a=0
  
  abs_prob=rep(0,length=fn)
  
  ### Calculating the survival function and absolute risk###    
  for (i in 1:fn)
  {
    
    if (start.age>=aa[i] & start.age<=bb[i] & start.age+pred.year<=bb[i]){   #### this condition applies when the prediction ends in the same interval as the current age
      abs_prob[i] = ((cbc.hazard[i]*RR)/(cbc.hazard[i]*RR+other.hazard[i]))*(1-exp(-(cbc.hazard[i]*RR+other.hazard[i])*(start.age+pred.year-start.age)))
      
    }
    else if (start.age>=aa[i] & start.age<=bb[i] & start.age+pred.year>bb[i]){ ### This allows the prediction to go further
      abs_prob[i] = ((cbc.hazard[i]*RR)/(cbc.hazard[i]*RR+other.hazard[i]))*(1-exp(-(cbc.hazard[i]*RR+other.hazard[i])*(bb[i]-start.age)))
      
      surv.cbc[i+1]=surv.cbc[i]*exp(-cbc.hazard[i]*RR*(bb[i]-aa[i]))   ### Survival function till the previous interval
      surv.other[i+1]= surv.other[i]*exp(-other.hazard[i]*(bb[i]-aa[i]))
      
      surv.cbc_a=surv.cbc[i]*exp(-cbc.hazard[i]*RR*(start.age-aa[i])) #### Survival function till start.age
      surv.other_a=surv.other[i]*exp(-other.hazard[i]*(start.age-aa[i]))     #### Survival function till start.age
    } 
    else if (start.age<aa[i] & start.age+pred.year>aa[i] & start.age+pred.year<=bb[i]){   ### This is for the probability calculated in the last interval of the prediction
      abs_prob[i]=((cbc.hazard[i]*RR)/(cbc.hazard[i]*RR+other.hazard[i]))*(surv.cbc[i]/surv.cbc_a)*(surv.other[i]/surv.other_a)*(1-exp(-(cbc.hazard[i]*RR+other.hazard[i])*(start.age+pred.year-aa[i]))) 
      surv.cbc[i+1]=0
      surv.other[i+1]=0
    }
    else if (start.age<aa[i] & start.age+pred.year>bb[i]){      #### Probability for an interval within the prediction years
      abs_prob[i]=((cbc.hazard[i]*RR)/(cbc.hazard[i]*RR+other.hazard[i]))*(surv.cbc[i]/surv.cbc_a)*(surv.other[i]/surv.other_a)*(1-exp(-(cbc.hazard[i]*RR+other.hazard[i])*(bb[i]-aa[i])))
      
      surv.cbc[i+1]=surv.cbc[i]*exp(-cbc.hazard[i]*RR*(bb[i]-aa[i]))
      surv.other[i+1]= surv.other[i]*exp(-other.hazard[i]*(bb[i]-aa[i]))  
    }
    else 
    {abs_prob[i]=0
    surv.cbc[i+1]=surv.cbc[i]*exp(-cbc.hazard[i]*RR*(bb[i]-aa[i]))
    surv.other[i+1]= surv.other[i]*exp(-other.hazard[i]*(bb[i]-aa[i])) 
    }
    
  } 
  
  result <- c(start.age+pred.year, round(100*sum(abs_prob),2))
  #names(result) <- c("current.age", "nyears", "abs.risk(%)")
  return(result)
}
pred_year=seq(pred.year,89,pred.year)
output=data.frame(matrix(  ,nrow = 0, ncol = 2))
colnames(output)=c("by age", "CBC risk(%)")
for (i in 1:length(pred_year))
{  
  if(start.age+pred_year[i] <=89)
  {output[i,]=abs_risk_cbc(profile,start.age,pred.year=pred_year[i],h1star.seer=T,h2.seer=T)
  }
  else {break}
  
}
# prof=data.frame(t(profile))
# 
# colnames(prof) <- c("age_1std","hormf","famhx_bc","lcis_atyph","brd_density","ER","bc_typ","cat_age1stb")
# 
# ## Recoding the inputs:
# prof$age_1std=recode(prof$age_1std,"1='<30';2='30-40';else='40+'",as.factor.result=T)
# prof$hormf=recode(prof$hormf,"1='Yes';2='No';else='Unk'",as.factor.result=T)
# prof$famhx_bc=recode(prof$famhx_bc,"1='Yes';2='No';else='Unk'",as.factor.result=T)
# prof$lcis_atyph=recode(prof$lcis_atyph,"1='Yes';else='Unk'",as.factor.result=T)
# prof$brd_density=recode(prof$brd_density, "4='Almost ent. fat';3='Scattered';2='Heterog dense';1='Extreme dense';else='Unk'",as.factor.result=T)
# prof$ER=recode(prof$ER,"1='Neg';2='Pos';else='Unk'",as.factor.result=T)
# prof$bc_typ=recode(prof$bc_typ,"1='Pure DCIS';2='Invasive_DCIS';else='Pure Invasive'",as.factor.result=T)
# prof$cat_age1stb=recode(prof$cat_age1stb,"1='<30/nulli';2='30-39';3='40+';else='Unk'",as.factor.result=T)
# 
# colnames(prof)=c("Age at first BC diagnosis", "Anti-estrogen therapy","Family history of BC","High Risk Preneoplasia","Breast density","ER","First BC type","Age at first birth")
# 
# row.names(prof) <- NULL
# row.names(output) <- NULL

#out.list=list("profile" = prof, "current_age" = start.age, "risk" = output)

if (print.output==T)
{return(output)}
else {invisible(output)}


}
