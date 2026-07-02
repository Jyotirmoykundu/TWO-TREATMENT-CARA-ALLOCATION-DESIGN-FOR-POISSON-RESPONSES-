#project_sem4
#different function of u ie e^u/1+e^u
u=rnorm(5000,mean=10,sd=3)
mean(exp(u)/(1+exp(u)))
myfunction=function(N,M,lambda_A,lambda_B)  #N denotes no. of patient
{
countA=NULL
propA=NULL
T1=NULL    #for test statistic 1
T2=NULL    #for test statistic 2
 for(z in 1:M)
 {
  X=NULL
  u=rnorm(N,10,3)
  u.trial=u[1:8]   
  n=length(u.trial)
  delta=NULL
  for(i in 1:n)
  {
    if(i<=n/2)
    {
     X[i]=rpois(n=1,lambda=lambda_A*(exp(u[i])/(1+exp(u[i]))))
     delta[i]=1
    }
    else
    {
     X[i]=rpois(n=1,lambda=lambda_B*(exp(u[i])/(1+exp(u[i]))))
     delta[i]=0
    }
  }
  Lam_nu.A=0
  Lam_de.A=0
  Lam_nu.B=0
  Lam_de.B=0
  for(i in 1:n)
  {
    if(delta[i]==1)
    {
     Lam_nu.A=Lam_nu.A+X[i]
     Lam_de.A=Lam_de.A+(exp(u[i])/(1+exp(u[i])))
    }
    else
    {
     Lam_nu.B=Lam_nu.B+X[i]
     Lam_de.B=Lam_de.B+(exp(u[i])/(1+exp(u[i])))
    }
  }
  lA.hat=(Lam_nu.A)/(Lam_de.A)    #lambdaA.hat upto 8 patients
  lB.hat=(Lam_nu.B)/(Lam_de.B)    #lambdaB.hat upto 8 patients
#------------function to define R_hat
  R.hat=function(m)
  {
    X_a=rpois(5000,lambda=lA.hat)
    X_b=rpois(5000,lambda=lB.hat)
    p1=0
    p2=0
    for(j in 1:5000)
    {
      if(X_a[j]<X_b[j])
       p1=p1+1
      else if(X_a[j]==X_b[j])
       p2=p2+1
      else
      {
       p1=p1
       p2=p2
      }
    }
    Ra.hat=(p1/5000)+0.5*(p2/5000)
    return(Ra.hat)
  }
  for(i in (n+1):N)
  {
    Ra.hat=R.hat(i-1)
    r=runif(1)
    if(r<=Ra.hat) 
    {
      Xa=rpois(1,lambda=lA.hat*(exp(u[i])/(1+exp(u[i]))))
      X=c(X,Xa)
      delta[i]=1
      Lam_nu.A=Lam_nu.A+X[i]
      Lam_de.A=Lam_de.A+(exp(u[i])/(1+exp(u[i])))
      lA.hat=Lam_nu.A/Lam_de.A
    }
    else
    {
      Xb=rpois(1,lambda=lB.hat*(exp(u[i])/(1+exp(u[i]))))
      X=c(X,Xb)
      delta[i]=0
      Lam_nu.B=Lam_nu.B+X[i]
      Lam_de.B=Lam_de.B+(exp(u[i])/(1+exp(u[i])))
      lB.hat=Lam_nu.B/Lam_de.B
    }
  }
  countA[z]=sum(delta)
  propA[z]=countA[z]/N
  lambda_n.hat=( lA.hat*sum(delta*(exp(u)/(1+exp(u))) + lB.hat*sum((1-delta)*(exp(u)/(1+exp(u)))) )/sum(exp(u)/(1+exp(u))))
  Shi_n.hat=(0.5)*0.9973
  T1[z]=(sqrt(N)*(lA.hat-lB.hat))/(sqrt((lambda_n.hat*2)/Shi_n.hat))
  T2[z]=(sqrt(N)*(lA.hat^(2/3)-lB.hat^(2/3)))/((2/3)*sqrt((2*lambda_n.hat^(1/3))/Shi_n.hat))
  }
  T0=sort(T1,decreasing=TRUE)
  T01=sort(T2,decreasing=TRUE)
  m=mean(propA)
  sd=sd(propA)
  crit_val_1=-1.35769190 # ---> T0[M-(M*0.05)]=-1.35769190(2,2)
  crit_val_2=-4.76108552  #---->T01[M-(M*0.05)]
  test=0
  for(t in 1:M)
  {
    if(T1[t]<crit_val_1)
    {
     test=test+1
    }
   else
    test=test
  }
  test1=0
  for(t in 1:M)
  {
    if(T2[t]<crit_val_2)
    {
     test1=test1+1
    }
   else
    test1=test1
  }
 P=test/M      #power corresponding to 1st test statistic
 P1=test1/M    #power corresponding to 2nd test statistic
  #---------for calculating bias:
  R.hat1=function(m)
  {
    X_a=rpois(m,lambda=lambda_A)
    X_b=rpois(m,lambda=lambda_B)
    p1=0
    p2=0
    for(j in 1:m)
    {
      if(X_a[j]<X_b[j])
       p1=p1+1
      else if(X_a[j]==X_b[j])
       p2=p2+1
      else
      {
       p1=p1
       p2=p2
      }
    }
    Ra.hat1=(p1/m)+0.5*(p2/m)
    return(Ra.hat1)
  }
  bias=mean(propA)-R.hat1(5000)
  mse=(sd(propA)^2)+(bias^2)
  return(c(m,sd,bias,mse,crit_val_1,crit_val_2,P,P1))
}
myfunction(100,1000,2,3)
#power calculation for test statistic 1
p1=myfunction(100,1000,2,2.4)[7]
p2=myfunction(100,1000,2,2.8)[7]
p3=myfunction(100,1000,2,3.2)[7]
p4=myfunction(100,1000,2,3.6)[7]
p5=myfunction(100,1000,2,4)[7]
p6=myfunction(100,1000,2,4.4)[7]
p7=myfunction(100,1000,2,4.8)[7]
p8=myfunction(100,1000,2,5.2)[7]
p9=myfunction(100,1000,2,5.6)[7]
p10=myfunction(100,1000,2,6)[7]
P1=c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)
#power calculation for test statistic 2
p11=myfunction(100,1000,2,2.4)[8]
p12=myfunction(100,1000,2,2.8)[8]
p13=myfunction(100,1000,2,3.2)[8]
p14=myfunction(100,1000,2,3.6)[8]
p15=myfunction(100,1000,2,4)[8]
p16=myfunction(100,1000,2,4.4)[8]
p17=myfunction(100,1000,2,4.8)[8]
p18=myfunction(100,1000,2,5.2)[8]
p19=myfunction(100,1000,2,5.6)[8]
p20=myfunction(100,1000,2,6)[8]
P_1=c(p11,p12,p13,p14,p15,p16,p17,p18,p19,p20)
x=seq(2.4,6,0.4)
#power curve using 1st test statistic
plot(x,P1,type="l",col="red",lwd=2,main="Comparison of power of the test for lambda A=2 and different values of lambda B",xlab="Lambda B",ylab="Power",ylim=c(0,max(max(P1),max(p))))
grid()
#same function and dist but different test statistics
plot(x,P1,type="l",col="green",lwd=2,main="Comparison of power for two seperate test statistics with lambdaA=2 and for different values of lambda B[e^u/1+e^u]",xlab="Lambda B",ylab="Power",ylim=c(0,max(max(p_1),max(p))))
lines(x,P_1,col="purple",lwd=2)
legend("topleft",col=c("green","purple"),legend=c("test_stat_1","test_stat_2"),lty=1:1,lwd=2:2,bty="n")
#consistency
c1=myfunction(100,1000,2,4)[7]
c2=myfunction(150,1000,2,4)[7]
c3=myfunction(200,1000,2,4)[7]
c4=myfunction(250,1000,2,4)[7]
c5=myfunction(300,1000,2,4)[7]
c6=myfunction(350,1000,2,4)[7]
c7=myfunction(400,1000,2,4)[7]
C=c(c1,c2,c3,c4,c5,c6,c7)

C1=myfunction(100,1000,2,4)[8]
C2=myfunction(150,1000,2,4)[8]
C3=myfunction(200,1000,2,4)[8]
C4=myfunction(250,1000,2,4)[8]
C5=myfunction(300,1000,2,4)[8]
C6=myfunction(350,1000,2,4)[8]
C7=myfunction(400,1000,2,4)[8]
C_1=c(C1,C2,C3,C4,C5,C6,C7)
s=seq(100,400,50)
plot(s,C,type="l",col="red",lwd=2,main="Checking consistency of the test based on two test statistics[for e^u/(1+e^u)]",xlab="sample size",ylab="Power",ylim=c(min(C,C_1),max(C,C_1)))
lines(s,C_1,col="green",lwd=2)
legend("topleft",col=c("red","green"),legend=c("test_stat_1","test_stat_2"),lty=1:1,lwd=2:2,bty="n")
grid()

