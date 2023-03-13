#'Gives estimated values of SSI and its confidence bounds
#'@param x1 is sample from Exp Population 1.
#'@param x2 is sample from Exp Population 2.
#'@param x3 is sample from Exp Population 3.
#'@param x4 is sample from Exp Population 4.
#'@param alpha is the level of significance
#'@export
CIEXP4<-function(x1,x2,x3,x4,alpha){


  k<-4
  nsim<-rep(1,k)
  nsim[1]<- length(x1)
  nsim[2]<- length(x2)
  nsim[3]<- length(x3)
  nsim[4]<- length(x4)
  boot<-300

  reliabbayes_boot_Tulika<-matrix(1,boot,k)
  reliabmle_boot_Tulika<-matrix(1,boot,k)
  reliabumvue_boot_Tulika<-matrix(1,boot,k)
  bse_boot_Tulika<-matrix(1,boot,k)
  reliabumvue<-rep(1,k)



  nsim
  N=max(nsim)
  takematrix<-matrix(0,k,N+1)

  n<-rep(0,k)



  #------------------------------------------------------------------------------
  #                            MLE
  #-----------------------------------------------------------------------------

  tmle=rep(0,k)
  tmle[1]=(mean(x1))
  tmle[2]=(mean(x2))
  tmle[3]=(mean(x3))
  tmle[4]=(mean(x4))


  reliabmle=rep(0,k)
  for(i in 1:k)
  {
    sum1=0
    for(j in 1:k)
    {
      sum1<-sum1+(tmle[i]/(tmle[i]+tmle[j]))

    }
    reliabmle[i]=1-sum1/k
  }
  reliabmle
  #------------------------------------------------------------------------------
  # UMVUE
  #------------------------------------------------------------------------------
  te<-rep(0,k)
  part<-matrix(0,k,k)
  s<-rep(0,k)
  s[1]<-sum(x1)
  s[2]<-sum(x2)
  s[3]<-sum(x3)
  s[4]<-sum(x4)
 # library(pracma)

  v<-matrix(0,k,k)
  for(i in 1:k)
  {
    for(j in 1:k){
      v[i,j]<-s[i]/s[j]
    }
  }

  for(i in 1:k){
    for(j in 1:k){
      f <- function(x, y) {(nsim[j]-1)*(nsim[i]-1)*((1-x)^(nsim[i]-2))*((1-y)^(nsim[j]-2))}
      if(v[i,j]>1){

        xmin <- 0; xmax <- 1/v[i,j]
        ymin <- function(x) x*v[i,j]; ymax <- 1
        I <- integral2(f, xmin, xmax, ymin, ymax)
        te[j]<-I$Q
      }
      else{
        x2min <- 0; x2max <- 1
        y2min <- 0; y2max <- function(x) x*v[i,j]
        J <- integral2(f, x2min, x2max, y2min, y2max)
        te[j]<-1-J$Q
      }
      part[i,j]<-te[j]
    }

  }
  reliabumvue[1]<-(1-(1/(2*k)))-((1/k)*(part[2,1]+part[3,1]+part[4,1] )       )
  reliabumvue[2]<-(1-(1/(2*k)))-((1/k)*(part[1,2]+part[3,2]+part[4,2])        )
  reliabumvue[3]<-(1-(1/(2*k)))-((1/k)*(part[1,3]+part[2,3]+part[4,3])        )
  reliabumvue[4]<-(1-(1/(2*k)))-((1/k)*(part[1,4]+part[2,4]+part[3,4])        )

  #------------------------------------------------------------------------------

  # GB
  #------------------------------------------------------------------------------

  reliabbayes<-rep(0,k)
  y<-rep(0,k)
  y[1]<-sum(x1)
  y[2]<-sum(x2)
  y[3]<-sum(x3)
  y[4]<-sum(x4)

  Gauss2F1b <- function(a,b,c,x){
    if(x>=0 & x<1){
      hyperg_2F1(a,b,c,x)
    }else{
      hyperg_2F1(a,c-b,c,1-1/(1-x))/(1-x)^a
    }
  }

  term<-function(a,b,c,d){
    term1<-((b/a)^(d))*(d)/(c+d)
    term2<-Gauss2F1b(1+d,c+d,c+d+1,1-(b/a))
    term1*term2
  }

  bayes<-function(k,a,b,c,d)
  {
    (1/k)*(term(a,b,c,d))
  }
  reliabbayes[1]<-1-(1/(2*k))-(bayes(4,y[1],y[2],nsim[1],nsim[2])+bayes(4,y[1],y[3],nsim[1],nsim[3])+bayes(4,y[1],y[4],nsim[1],nsim[4]))
  reliabbayes[2]<-1-(1/(2*k))-(bayes(4,y[2],y[1],nsim[2],nsim[1])+bayes(4,y[2],y[3],nsim[2],nsim[3])+bayes(4,y[2],y[4],nsim[2],nsim[4]))
  reliabbayes[3]<-1-(1/(2*k))-(bayes(4,y[3],y[1],nsim[3],nsim[1])+bayes(4,y[3],y[2],nsim[3],nsim[2])+bayes(4,y[3],y[4],nsim[3],nsim[4]))
  reliabbayes[4]<-1-(1/(2*k))-(bayes(4,y[4],y[1],nsim[4],nsim[1])+bayes(4,y[4],y[2],nsim[4],nsim[2])+bayes(4,y[4],y[3],nsim[4],nsim[3]))


  #------------------------------------------------------------------------------

  # BSE
  #------------------------------------------------------------------------------

  s<-rep(0,k)
  s[1]<-sum(x1)
  s[2]<-sum(x2)
  s[3]<-sum(x3)
  s[4]<-sum(x4)
  bse<-rep(0,k)

  d<-matrix(0,k,k)
  for (i in 1:k) {
    for (j in 1:k) {
      d[i,j]<-((nsim[i]-2)*s[j])/((nsim[j]+1)*s[i])
    }
  }
  for (i in 1:k) {
    s2<-0
    for (j in 1:k) {
      s2<-s2+(1/(1+d[i,j]))
    }
    bse[i]<-1-(1/k)*(s2)
  }







  for (l in 1:boot) {
    x1_boot<- rexp(nsim[1],1/tmle[1])
    x2_boot<- rexp(nsim[2],1/tmle[2])
    x3_boot<- rexp(nsim[3],1/tmle[3])
    x4_boot<- rexp(nsim[4],1/tmle[4])

    tmle_boot=rep(0,k)
    tmle_boot[1]=(mean(x1_boot))
    tmle_boot[2]=(mean(x2_boot))
    tmle_boot[3]=(mean(x3_boot))
    tmle_boot[4]=(mean(x4_boot))


    reliabmle_boot=rep(0,k)
    for(i in 1:k)
    {
      sum1_boot=0
      for(j in 1:k)
      {
        sum1_boot<-sum1_boot+(tmle_boot[i]/(tmle_boot[i]+tmle_boot[j]))

      }
      reliabmle_boot[i]=1-sum1_boot/k
    }
    #------------------------------------------------------------------------------

    # UMVUE
    #------------------------------------------------------------------------------

    te_boot<-rep(0,k)
    part_boot<-matrix(0,k,k)
    s_boot<-rep(0,k)
    s_boot[1]<-sum(x1_boot)
    s_boot[2]<-sum(x2_boot)
    s_boot[3]<-sum(x3_boot)
    s_boot[4]<-sum(x4_boot)
 #   library(pracma)

    v_boot<-matrix(0,k,k)
    for(i in 1:k)
    {
      for(j in 1:k){
        v_boot[i,j]<-s_boot[i]/s_boot[j]
      }
    }

    for(i in 1:k){
      for(j in 1:k){
        f_boot <- function(x, y) {(nsim[j]-1)*(nsim[i]-1)*((1-x)^(nsim[i]-2))*((1-y)^(nsim[j]-2))}
        if(v_boot[i,j]>1){

          xmin_boot <- 0; xmax_boot <- 1/v_boot[i,j]
          ymin_boot <- function(x) x*v_boot[i,j]; ymax_boot <- 1
          I_boot <- integral2(f_boot, xmin_boot, xmax_boot, ymin_boot, ymax_boot)
          te_boot[j]<-I_boot$Q
        }
        if(v_boot[i,j]<1){
          x2min_boot <- 0; x2max_boot <- 1
          y2min_boot <- 0; y2max_boot <- function(x) x*v_boot[i,j]
          J_boot <- integral2(f_boot, x2min_boot, x2max_boot, y2min_boot, y2max_boot)
          te_boot[j]<-1-J_boot$Q
        }
        part_boot[i,j]<-te_boot[j]
      }

    }
    reliabumvue_boot<-rep(1,k)

    reliabumvue_boot[1]<-(1-(1/(2*k)))-((1/k)*(part_boot[2,1]+part_boot[3,1]+part_boot[4,1] )       )
    reliabumvue_boot[2]<-(1-(1/(2*k)))-((1/k)*(part_boot[1,2]+part_boot[3,2]+part_boot[4,2])        )
    reliabumvue_boot[3]<-(1-(1/(2*k)))-((1/k)*(part_boot[1,3]+part_boot[2,3]+part_boot[4,3])        )
    reliabumvue_boot[4]<-(1-(1/(2*k)))-((1/k)*(part_boot[1,4]+part_boot[2,4]+part_boot[3,4])        )

    #------------------------------------------------------------------------------

    # GB
    #------------------------------------------------------------------------------
    reliabbayes_boot<-rep(0,k)
    y_boot<-rep(0,k)
    y_boot[1]<-sum(x1_boot)
    y_boot[2]<-sum(x2_boot)
    y_boot[3]<-sum(x3_boot)
    y_boot[4]<-sum(x4_boot)


    reliabbayes_boot[1]<-1-(1/(2*k))-(bayes(4,y_boot[1],y_boot[2],nsim[1],nsim[2])+bayes(4,y_boot[1],y_boot[3],nsim[1],nsim[3])+bayes(4,y_boot[1],y_boot[4],nsim[1],nsim[4]))
    reliabbayes_boot[2]<-1-(1/(2*k))-(bayes(4,y_boot[2],y_boot[1],nsim[2],nsim[1])+bayes(4,y_boot[2],y_boot[3],nsim[2],nsim[3])+bayes(4,y_boot[2],y_boot[4],nsim[2],nsim[4]))
    reliabbayes_boot[3]<-1-(1/(2*k))-(bayes(4,y_boot[3],y_boot[1],nsim[3],nsim[1])+bayes(4,y_boot[3],y_boot[2],nsim[3],nsim[2])+bayes(4,y_boot[3],y_boot[4],nsim[3],nsim[4]))
    reliabbayes_boot[4]<-1-(1/(2*k))-(bayes(4,y_boot[4],y_boot[1],nsim[4],nsim[1])+bayes(4,y_boot[4],y_boot[2],nsim[4],nsim[2])+bayes(4,y_boot[4],y_boot[3],nsim[4],nsim[3]))

    #------------------------------------------------------------------------------

    # BSE
    #------------------------------------------------------------------------------

    bse_boot<-rep(0,k)

    d_boot<-matrix(0,k,k)
    for (i in 1:k) {
      for (j in 1:k) {
        d_boot[i,j]<-((nsim[i]-2)*s_boot[j])/((nsim[j]+1)*s_boot[i])
      }
    }
    for (i in 1:k) {
      s2<-0
      for (j in 1:k) {
        s2<-s2+(1/(1+d_boot[i,j]))
      }
      bse_boot[i]<-1-(1/k)*(s2)
    }

    for(d in 1:k)
    {

      reliabmle_boot_Tulika[l,d]<-reliabmle_boot[d]
      reliabumvue_boot_Tulika[l,d]<-reliabumvue_boot[d]
      reliabbayes_boot_Tulika[l,d]<-reliabbayes_boot[d]
      bse_boot_Tulika[l,d]<-bse_boot[d]

    }


  }
  reliabmle_boot_1<-reliabmle_boot_Tulika[ ,1]
  reliabmle_boot_2<-reliabmle_boot_Tulika[ ,2]
  reliabmle_boot_3<-reliabmle_boot_Tulika[ ,3]
  reliabmle_boot_4<-reliabmle_boot_Tulika[ ,4]
  #  reliabmle_boot_5<-reliabmle_boot_Tulika[ ,5]
  #   p1
  new_lambda_mle_1<-reliabmle_boot_1[order(reliabmle_boot_1)]
  C_star_mle1_1<-new_lambda_mle_1[floor((alpha/2)*boot)]
  C_star_mle2_1<-new_lambda_mle_1[floor((1-(alpha/2))*boot)]
  C_star_mle1_1
  C_star_mle2_1

  # length(lambda[lambda<c_star])/outer
  #   p2
  new_lambda_mle_2<-reliabmle_boot_2[order(reliabmle_boot_2)]
  C_star_mle1_2<-new_lambda_mle_2[floor((alpha/2)*boot)]
  C_star_mle2_2<-new_lambda_mle_2[floor((1-(alpha/2))*boot)]
  C_star_mle1_2
  C_star_mle2_2

  #   p3
  new_lambda_mle_3<-reliabmle_boot_3[order(reliabmle_boot_3)]
  C_star_mle1_3<-new_lambda_mle_3[floor((alpha/2)*boot)]
  C_star_mle2_3<-new_lambda_mle_3[floor((1-(alpha/2))*boot)]
  C_star_mle1_3
  C_star_mle2_3

  #   p4
  new_lambda_mle_4<-reliabmle_boot_4[order(reliabmle_boot_4)]
  C_star_mle1_4<-new_lambda_mle_4[floor((alpha/2)*boot)]
  C_star_mle2_4<-new_lambda_mle_4[floor((1-(alpha/2))*boot)]
  C_star_mle1_4
  C_star_mle2_4

  ################################################
  ##################### UMVUE
  ##############################################
  reliabumvue_boot_1<-reliabumvue_boot_Tulika[ ,1]
  reliabumvue_boot_2<-reliabumvue_boot_Tulika[ ,2]
  reliabumvue_boot_3<-reliabumvue_boot_Tulika[ ,3]
  reliabumvue_boot_4<-reliabumvue_boot_Tulika[ ,4]
  #   p1
  new_lambda_umvue_1<-reliabumvue_boot_1[order(reliabumvue_boot_1)]
  C_star_umvue1_1<-new_lambda_umvue_1[floor((alpha/2)*boot)]
  C_star_umvue2_1<-new_lambda_umvue_1[floor((1-(alpha/2))*boot)]

  ub_umvue<-c(C_star_umvue2_1,C_star_umvue2_1,C_star_umvue2_1,C_star_umvue2_1)
  # length(lambda[lambda<c_star])/outer
  #   p2
  new_lambda_umvue_2<-reliabumvue_boot_2[order(reliabumvue_boot_2)]
  C_star_umvue1_2<-new_lambda_umvue_2[floor((alpha/2)*boot)]
  C_star_umvue2_2<-new_lambda_umvue_2[floor((1-(alpha/2))*boot)]


  #   p3
  new_lambda_umvue_3<-reliabumvue_boot_3[order(reliabumvue_boot_3)]
  C_star_umvue1_3<-new_lambda_umvue_3[floor((alpha/2)*boot)]
  C_star_umvue2_3<-new_lambda_umvue_3[floor((1-(alpha/2))*boot)]


  #   p4
  new_lambda_umvue_4<-reliabumvue_boot_4[order(reliabumvue_boot_4)]
  C_star_umvue1_4<-new_lambda_umvue_4[floor((alpha/2)*boot)]
  C_star_umvue2_4<-new_lambda_umvue_4[floor((1-(alpha/2))*boot)]


  ########################### GB

  reliabbayes_boot_1<-reliabbayes_boot_Tulika[ ,1]
  reliabbayes_boot_2<-reliabbayes_boot_Tulika[ ,2]
  reliabbayes_boot_3<-reliabbayes_boot_Tulika[ ,3]
  reliabbayes_boot_4<-reliabbayes_boot_Tulika[ ,4]
  #  reliabbayes_boot_5<-reliabbayes_boot_Tulika[ ,5]
  #   p1
  new_lambda_bayes_1<-reliabbayes_boot_1[order(reliabbayes_boot_1)]
  C_star_bayes1_1<-new_lambda_bayes_1[floor((alpha/2)*boot)]
  C_star_bayes2_1<-new_lambda_bayes_1[floor((1-(alpha/2))*boot)]
  C_star_bayes1_1
  C_star_bayes2_1

  # length(lambda[lambda<c_star])/outer
  #   p2
  new_lambda_bayes_2<-reliabbayes_boot_2[order(reliabbayes_boot_2)]
  C_star_bayes1_2<-new_lambda_bayes_2[floor((alpha/2)*boot)]
  C_star_bayes2_2<-new_lambda_bayes_2[floor((1-(alpha/2))*boot)]
  C_star_bayes1_2
  C_star_bayes2_2

  #   p3
  new_lambda_bayes_3<-reliabbayes_boot_3[order(reliabbayes_boot_3)]
  C_star_bayes1_3<-new_lambda_bayes_3[floor((alpha/2)*boot)]
  C_star_bayes2_3<-new_lambda_bayes_3[floor((1-(alpha/2))*boot)]
  C_star_bayes1_3
  C_star_bayes2_3
  new_lambda_bayes_3

  #   p4
  new_lambda_bayes_4<-reliabbayes_boot_4[order(reliabbayes_boot_4)]
  C_star_bayes1_4<-new_lambda_bayes_4[floor((alpha/2)*boot)]
  C_star_bayes2_4<-new_lambda_bayes_4[floor((1-(alpha/2))*boot)]
  C_star_bayes1_4
  C_star_bayes2_4


  ############################## BSE

  bse_boot_1<-bse_boot_Tulika[ ,1]
  bse_boot_2<-bse_boot_Tulika[ ,2]
  bse_boot_3<-bse_boot_Tulika[ ,3]
  bse_boot_4<-bse_boot_Tulika[ ,4]
  #   p1
  new_lambda_bse_1<-bse_boot_1[order(bse_boot_1)]
  C_star_bse1_1<-new_lambda_bse_1[floor((alpha/2)*boot)]
  C_star_bse2_1<-new_lambda_bse_1[floor((1-(alpha/2))*boot)]
  C_star_bse1_1
  C_star_bse2_1

  # length(lambda[lambda<c_star])/outer
  #   p2
  new_lambda_bse_2<-bse_boot_2[order(bse_boot_2)]
  C_star_bse1_2<-new_lambda_bse_2[floor((alpha/2)*boot)]
  C_star_bse2_2<-new_lambda_bse_2[floor((1-(alpha/2))*boot)]
  C_star_bse1_2
  C_star_bse2_2

  #   p3
  new_lambda_bse_3<-bse_boot_3[order(bse_boot_3)]
  C_star_bse1_3<-new_lambda_bse_3[floor((alpha/2)*boot)]
  C_star_bse2_3<-new_lambda_bse_3[floor((1-(alpha/2))*boot)]
  C_star_bse1_3
  C_star_bse2_3

  #   p4
  new_lambda_bse_4<-bse_boot_4[order(bse_boot_4)]
  C_star_bse1_4<-new_lambda_bse_4[floor((alpha/2)*boot)]
  C_star_bse2_4<-new_lambda_bse_4[floor((1-(alpha/2))*boot)]
  C_star_bse1_4
  C_star_bse2_4

  lb_mle<-c(C_star_mle1_1,C_star_mle1_2,C_star_mle1_3,C_star_mle1_4)
  ub_mle<-c(C_star_mle2_1,C_star_mle2_2,C_star_mle2_3,C_star_mle2_4)
  len_mle<-ub_mle-lb_mle
  df_mle<-data.frame(estimator=reliabmle,lower.bound=lb_mle,upper.bound=ub_mle,lengths.interval=len_mle)

  lb_umvue<-c(C_star_umvue1_1,C_star_umvue1_2,C_star_umvue1_3,C_star_umvue1_4)
  ub_umvue<-c(C_star_umvue2_1,C_star_umvue2_2,C_star_umvue2_3,C_star_umvue2_4)
  len_umvue<-ub_umvue-lb_umvue
  df_umvue<-data.frame(estimator=reliabumvue,lower.bound=lb_umvue,upper.bound=ub_umvue,lengths.interval=len_umvue)

  lb_bayes<-c(C_star_bayes1_1,C_star_bayes1_2,C_star_bayes1_3,C_star_bayes1_4)
  ub_bayes<-c(C_star_bayes2_1,C_star_bayes2_2,C_star_bayes2_3,C_star_bayes2_4)
  len_bayes<-ub_bayes-lb_bayes
  df_bayes<-data.frame(estimator=reliabbayes,lower.bound=lb_bayes,upper.bound=ub_bayes,lengths.interval=len_bayes)


  lb_bse<-c(C_star_bse1_1,C_star_bse1_2,C_star_bse1_3,C_star_bse1_4)
  ub_bse<-c(C_star_bse2_1,C_star_bse2_2,C_star_bse2_3,C_star_bse2_4)
  len_bse<-ub_bse-lb_bse
  df_bse<-data.frame(estimator=bse,lower.bound=lb_bse,upper.bound=ub_bse,lengths.interval=len_bse)
  dflist<-list(df_mle,df_umvue,df_bayes,df_bse)
  dflist
}
