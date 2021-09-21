library(methods)

upperage=10000
lowerage=10

mitotic.age<-function(beta,b,c,d)
{
n=rep(500,ncol(beta)) #true parameter
            no.initial.n=5

### minimize the objective function "fr.j" for each patient in turn  #########
        for(j in 1:ncol(beta))
        {
            current.value=fr.j(n[j],b,c,d,beta[,j])

            columnoptim=vector("list",no.initial.n)
            val=rep(NA,no.initial.n)
            for(jj in 1:(no.initial.n-1)) 
            {
            temp=try(optim(par=lowerage+jj*(upperage-lowerage)/no.initial.n,fn=fr.j,gr=grr.j,b=b,c=c,d=d,betaj=beta[,j],method="L-BFGS-B",lower=lowerage,upper=upperage,control=list(factr=1)),silent=T)
            if(!is(temp,"try-error")){ columnoptim[[jj]]=temp; val[jj]=temp$value}
            } 
           temp=try(optim(par=n[j],fn=fr.j,gr=grr.j,b=b,c=c,d=d,betaj=beta[,j],method="L-BFGS-B",lower=lowerage,upper=upperage,control=list(factr=1)),silent=T)
            if(!is(temp,"try-error")) { columnoptim[[no.initial.n]]=temp; val[no.initial.n]=temp$value}
           temp=columnoptim[[ which(val==min(val,na.rm=T))[1] ]]


            if(!is(temp,"try-error"))
            {
                n[j]=temp$par
            } else {print(2);print(temp)}

            }

return(n)

}

fr.i <- function(x,n,betai)  ## objective function summed over all patients for CpG site i
{
    b <- x[1]
    c <- x[2]
    d <- x[3]

    return(sum((c+b^(n-1)*d-betai)^2,na.rm=T))
}

grr.i <- function(x,n,betai)  ## derivative of fr.i with respect to b,c,d
{
    b <- x[1]
    c <- x[2]
    d <- x[3]

    return(c( 2*sum( (c+b^(n-1)*d-betai)*(n-1)*b^(n-2)*d,na.rm=T ), 2*sum(c+b^(n-1)*d-betai,na.rm=T) , 2*sum( (c+b^(n-1)*d-betai)*b^(n-1),na.rm=T ))   )
}

fr.j <- function(x,b,c,d,betaj)  ## objective fuction summed over all CpG istes for patient j
{
    nj=x
    return(sum((c+b^(nj-1)*d-betaj)^2,na.rm=T))
}

grr.j <- function(x,b,c,d,betaj)  ## derivative of fr.j with respect to n_{j}
{
    nj=x
    return(2*sum((c+b^(nj-1)*d-betaj)*b^(nj-1)*log(b)*d,na.rm=T))
}





total <- function(n,b,c,d,beta)  ## objective function summed over all CpG sites and all patients
{
    result=0
    for(j in 1:length(n))
        result=result+sum((c+b^(n[j]-1)*d-beta[,j])^2,na.rm=T)

    return(result)
}


conversion.to.aqr<-function(x)  ### conversion of parameters b,c,d to a,q,r
{
#a=(delta.p+delta.d)/2
#q=(1+mu)/2
#r=E(X_1)

a=x[2]*(1-x[1])
q=a+x[1]
r=x[2]+x[3]
return(c(a,q,r))
}


conversion.to.bcd<-function(x)  ### conversion of parameters a,q,r to b,c,d
{
b=x[2]-x[1]
if(x[1]==0) c=0 else c=x[1]/(1-x[2]+x[1])
d=x[3]-c
return(c(b,c,d))
}

derivative.matrix<-function(x)
{
t1=x[1]/(1-x[2]+x[1])^2
t2=(1-x[2])/(1-x[2]+x[1])^2

return(rbind(c(-1,1,0), c(t2,t1,0), c(-t2,-t1,1)))
}

new.grr.i<-function(x,n,betai) ## derivative of fr.i with respect to a,q,r
{
return(rbind(grr.i( conversion.to.bcd(x),n,betai)) %*% derivative.matrix(x))
}
new.fr.i<-function(x,n,betai)  ## fr.i as a function of a,q,r
{
return(fr.i(conversion.to.bcd(x),n,betai))
}


## range of parameter values #########
q.min=0.5    
a.max=1
new.ui=rbind(c(1,0,0),c(0,1,0),c(0,0,1))
new.ui=rbind(new.ui,-new.ui)
new.ui=rbind(new.ui,c(-1,1,0))
new.ci=c(0, q.min , 0, -a.max , -1 , -1 , 0)-10^(-10)





ineqfun1=function(x,n,betai){ return(x[2]-x[1])}




age.estimate <- function(beta,temp)  ### main function to estimate tumor age
{
    initialb=0.99
    
### setting initial values for parameters a,q,r,b,c,d,n  ####
    a=rep(0.0001,nrow(beta))
    q=a+initialb
    r=rep(0.5,nrow(beta))

    b=rep(initialb,nrow(beta))
    c=a/(1-b)
    d=r-c

    largest=upperage*3/10
    smallest=lowerage*10
    n=(largest-smallest)/(max(temp)-min(temp))*(temp-min(temp))+smallest

    converge=0.00001

    count=0

    old.value=Inf
    new.value=10^10


### minimizing the objective function "total"  ##################
    while(abs(old.value-new.value)> converge & count < 30)
    {
        print(count)
        old.value=new.value

### minimize the objective function "new.fr.i" for each CpG site in turn  #########

        for(i in 1:nrow(beta))
        {#print(i)
            current.value=fr.i(c(b[i],c[i],d[i]),n,beta[i,])
            
            init=  conversion.to.aqr(c(b[i],c[i],d[i]))
            init=find.init(init,n,beta[i,])


            temp=try(constrOptim(theta=init, f=new.fr.i, grad=new.grr.i, ui=new.ui, ci=new.ci,n=n,betai=beta[i,]),silent=T)
            if(!is(temp,"try-error"))
            {
               bcd=conversion.to.bcd(temp$par)
            } else{
             temp=solnp(pars=init,fun=new.fr.i,ineqfun=ineqfun1,ineqLB=0,ineqUB=1,LB=c(0,q.min,0),UB=c(a.max,1,1),n=n,betai=beta[i,],control=list(trace=0))
             bcd=conversion.to.bcd(temp$pars)
            }
          
               b[i]=bcd[1]
               c[i]=bcd[2]
               d[i]=bcd[3]
           
        if( fr.i(c(b[i],c[i],d[i]),n,beta[i,])-current.value>10^(-9) ) 
            print(c("error1",i,count,fr.i(c(b[i],c[i],d[i]),n,beta[i,])-current.value))  ### check if the objective function gets reduced #####  


        if( fr.i(c(b[i],c[i],d[i]),n,beta[i,])-  fr.i(c(1,0,mean(beta[i,],na.rm=T)),n,beta[i,])   > - penaltycutoff ) 
            {b[i]=1;c[i]=0;d[i]=mean(beta[i,],na.rm=T);  }  ### filter out noninformative CpG sites ####


        }


        temp.new.value=total(n,b,c,d,beta)
        if( temp.new.value  >old.value) print("error2")  ### check if the objective function gets reduced ###

        
### minimize the objective function "fr.j" for each patient in turn  #########
        for(j in 1:ncol(beta))
        {
            print(j)
            current.value=fr.j(n[j],b,c,d,beta[,j])
            print(n[j])
            print(current.value)      
            no.initial.n=5
            columnoptim=vector("list",no.initial.n)
            val=rep(NA,no.initial.n)
            for(jj in 1:(no.initial.n-1)) 
            {
            temp=try(optim(par=lowerage+jj*(upperage-lowerage)/no.initial.n,fn=fr.j,gr=grr.j,b=b,c=c,d=d,betaj=beta[,j],method="L-BFGS-B",lower=lowerage,upper=upperage,control=list(factr=1)),silent=T)
            if(!is(temp,"try-error")){ columnoptim[[jj]]=temp; val[jj]=temp$value}
            } 
           temp=try(optim(par=n[j],fn=fr.j,gr=grr.j,b=b,c=c,d=d,betaj=beta[,j],method="L-BFGS-B",lower=lowerage,upper=upperage,control=list(factr=1)),silent=T)
            if(!is(temp,"try-error")) { columnoptim[[no.initial.n]]=temp; val[no.initial.n]=temp$value}
           temp=columnoptim[[ which(val==min(val,na.rm=T))[1] ]]


#print(temp)
            if(!is(temp,"try-error"))
            {
                n[j]=temp$par
            } else {print(2);print(temp)}

            if(fr.j(n[j],b,c,d,beta[,j])-current.value>10^(-9)) print(c("error3",j,count,fr.j(n[j],b,c,d,beta[,j])-current.value))  ### check if the objective function gets reduced ###
        }


        count=count+1
        new.value=   total(n,b,c,d,beta)
        if(new.value>temp.new.value) print("error4")
        print(new.value)
                
    }
    return(list(b,c,d,n))
}





find.init<-function(x,n,betai)  ### finding initial values for a,q,r out of 50 different values of b
{
b=c(runif(40,min=0,max=0.9),runif(10,min=0.9,max=1))
init=matrix(NA,nr=50,nc=3)
for(i in 1:50) init[i,]=find.init2(b[i],n,betai)
init=rbind(init,x)
y=apply(init,1,new.fr.i,n,betai)
id=which(y==min(y,na.rm=T))[1]
init=init[id,]
return(init)

}

find.init2<-function(b,n,betai)  ### finding initial values for a,q,r for the given b
{

min.a=max(0,q.min-b)
max.a=min(1-b,a.max)
min.r=0
max.r=1

c2= b^(n-1)
c1= (1-c2)/(1-b)
c3= -betai
d1=sum(c1^2)
d2=sum(c1*c2)
d3=sum(c1*c3)
d4=sum(c2^2)
d5=sum(c2*c3)


a=(d3*d4-d2*d5)/(d2^2-d1*d4)
r=(d1*d5-d2*d3)/(d2^2-d1*d4)
init=c(a,r)

r=0:1;a=(-d2*r-d3)/d1
init=rbind(init,cbind(a,r))

a=c(min.a,max.a);r= (-d2*a-d5)/d4
init=rbind(init,cbind(a,r))

init=init[rowSums(is.na(init))==0,]
if(!is.matrix(init)) init=matrix(init,nc=2)

sel= init[,1]>=min.a & init[,1]<=max.a & init[,2]>=0 & init[,2]<=1
init=init[sel,]

init=rbind(init,c(min.a,0),c(max.a,0),c(min.a,1),c(max.a,1))


x=rep(NA,nrow(init));for(i in 1:nrow(init)) x[i]=sum((c1*init[i,1]+c2*init[i,2]+c3)^2,na.rm=T)
#x=apply( init,1,fr.i2,b,n,betai)

id=which(x==min(x,na.rm=T))[1]
init=init[id,]


return(c(init[1],b+init[1],init[2]))
}

