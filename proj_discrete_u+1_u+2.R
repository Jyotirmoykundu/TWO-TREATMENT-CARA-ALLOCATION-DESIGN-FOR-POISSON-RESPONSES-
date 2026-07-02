#project_sem4
#different function of u ie (u+1)/(u+2)
u=sample(c(0,1,2),replace=T,prob=c(0.25,0.5,0.25),size=5000)
mean((u+1)/(u+2))
myfunction=function(N,M,lambda_A,lambda_B)  #N denotes no. of patient
{
countA=NULL
propA=NULL
T1=NULL    #for test statistic 1
T2=NULL    #for test statistic 2
 for(z in 1:M)
 {
  X=NULL
  u=sample(c(0,1,2),replace=T,prob=c(0.25,0.5,0.25),size=N)
  u.trial=u[1:8]   
  n=length(u.trial)
  delta=NULL
  for(i in 1:n)
  {
    if(i<=n/2)
    {
     X[i]=rpois(n=1,lambda=lambda_A*((u[i]+1)/(u[i]+2)))
     delta[i]=1
    }
    else
    {
     X[i]=rpois(n=1,lambda=lambda_B*((u[i]+1)/(u[i]+2)))
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
     Lam_de.A=Lam_de.A+((u[i]+1)/(u[i]+2))
    }
    else
    {
     Lam_nu.B=Lam_nu.B+X[i]
     Lam_de.B=Lam_de.B+((u[i]+1)/(u[i]+2))
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
      Xa=rpois(1,lambda=lA.hat*((u[i]+1)/(u[i]+2)))
      X=c(X,Xa)
      delta[i]=1
      Lam_nu.A=Lam_nu.A+X[i]
      Lam_de.A=Lam_de.A+((u[i]+1)/(u[i]+2))
      lA.hat=Lam_nu.A/Lam_de.A
    }
    else
    {
      Xb=rpois(1,lambda=lB.hat*((u[i]+1)/(u[i]+2)))
      X=c(X,Xb)
      delta[i]=0
      Lam_nu.B=Lam_nu.B+X[i]
      Lam_de.B=Lam_de.B+((u[i]+1)/(u[i]+2))
      lB.hat=Lam_nu.B/Lam_de.B
    }
  }
  countA[z]=sum(delta)
  propA[z]=countA[z]/N
  lambda_n.hat=( lA.hat*sum(delta*((u+1)/(u+2))) + lB.hat*sum((1-delta)*((u+1)/(u+2))) )/sum(((u+1)/(u+2)))
  Shi_n.hat=(0.5)*0.64545
  T1[z]=(sqrt(N)*(lA.hat-lB.hat))/(sqrt((lambda_n.hat*2)/Shi_n.hat))
  T2[z]=(sqrt(N)*(lA.hat^(2/3)-lB.hat^(2/3)))/((2/3)*sqrt((2*lambda_n.hat^(1/3))/Shi_n.hat))
  }
  T0=sort(T1,decreasing=TRUE)
  T01=sort(T2,decreasing=TRUE)
  m=mean(propA)
  sd=sd(propA)
  crit_val_1=-10.35574649  # ---> T0[M-(M*0.05)]=-10.35574649(2,2)
  crit_val_2=-9.07531806   #---> T01[M-(M*0.05)]
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
p.1=myfunction(100,1000,2,2.4)[7]
p.2=myfunction(100,1000,2,2.8)[7]
p.3=myfunction(100,1000,2,3.2)[7]
p.4=myfunction(100,1000,2,3.6)[7]
p.5=myfunction(100,1000,2,4)[7]
p.6=myfunction(100,1000,2,4.4)[7]
p.7=myfunction(100,1000,2,4.8)[7]
p.8=myfunction(100,1000,2,5.2)[7]
p.9=myfunction(100,1000,2,5.6)[7]
p.10=myfunction(100,1000,2,6)[7]
P.2=c(p.1,p.2,p.3,p.4,p.5,p.6,p.7,p.8,p.9,p.10)
#power calculation for test statistic 2
p.11=myfunction(100,1000,2,2.4)[8]
p.12=myfunction(100,1000,2,2.8)[8]
p.13=myfunction(100,1000,2,3.2)[8]
p.14=myfunction(100,1000,2,3.6)[8]
p.15=myfunction(100,1000,2,4)[8]
p.16=myfunction(100,1000,2,4.4)[8]
p.17=myfunction(100,1000,2,4.8)[8]
p.18=myfunction(100,1000,2,5.2)[8]
p.19=myfunction(100,1000,2,5.6)[8]
p.20=myfunction(100,1000,2,6)[8]
P.22=c(p.11,p.12,p.13,p.14,p.15,p.16,p.17,p.18,p.19,p.20)
x=seq(2.4,6,0.4)
#same function and dist but different test statistics
plot(x,P.2,type="l",col="orange",lwd=2,main="Comparison of power for two seperate test statistics with lambdaA=2 and for different values of lambda B[(u+1)/(u+2)]",xlab="Lambda B",ylab="Power",ylim=c(0,max(max(p_1),max(p))))
lines(x,P.22,col="purple",lwd=2)
legend("topleft",col=c("orange","purple"),legend=c("test_stat_1","test_stat_2"),lty=1:1,lwd=2:2,bty="n")
#power curve using 1st test statistic
plot(x,p,type="l",col="red",lwd=2,main="Comparison of power of the test for lambda A and different values of lambda B",xlab="Lambda B",ylab="Power",ylim=c(0,max(max(p),max(P1))))
grid()
#consistency
c1=myfunction(100,1000,2,4)[c(7,8)]
c2=myfunction(150,1000,2,4)[c(7,8)]
c3=myfunction(200,1000,2,4)[c(7,8)]
c4=myfunction(250,1000,2,4)[c(7,8)]
c5=myfunction(300,1000,2,4)[c(7,8)]
c6=myfunction(350,1000,2,4)[c(7,8)]
c7=myfunction(400,1000,2,4)[c(7,8)]
C=c(c1,c2,c3,c4,c5,c6,c7)
odd=seq(1,14,2)
even=seq(2,14,2)
s=seq(100,400,50)
plot(s,C[odd],type="l",col="red",lwd=2,main="Checking consistency of the test based on two test statistics[for (u+1)/(u+2)]",xlab="sample size",ylab="Power",ylim=c(min(C),max(C)))
lines(s,C[even],col="green",lwd=2)
legend("topleft",col=c("red","green"),legend=c("test_stat_1","test_stat_2"),lty=1:1,lwd=2:2,bty="n")
grid()


