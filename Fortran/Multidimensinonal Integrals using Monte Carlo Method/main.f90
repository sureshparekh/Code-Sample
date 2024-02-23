!Multidimentional Integrals using Monte Carlo sample mean method with Uniform Sampling

program MultInte
    implicit none
    real:: r,pi,sumfx,x(3),integral,f1,f,sumfx2,sigma,error
    integer::i,j,n,nmc

    n=3           					!total number of dimentions
    nmc=1000000   					!number of MC cycles
    pi=3.1415       
    sumfx=0.0  
    sumfx2=0.0
    do i=1,nmc    					!for MC Cycles

        do j=1,n    					!to generate array 
            r=rand()     				!Generating Random Number
            x(j)=pi*r      				!dx=a+(b-a)r here [a,b]=[0,pi]
        enddo
        
        f1=f(x(1),x(2),x(3))
        sumfx =sumfx+f1          			!summing values
        sumfx2=sumfx2+f1**2      			!summing squares
    enddo
    
    integral=(sumfx*pi**3)/nmc              		!Integral = {(f-e)(d-c)(b-a)/n}*sum(f(xi,x2,x3)), b,d,f=pi a,c,e=0
    
    sumfx=sumfx/nmc                			!taking mean of fx
    sumfx2=sumfx2/nmc              			!taking mean of fx**2
    sigma=sqrt(abs(sumfx**2 - sumfx2))  		!calculating sigma of fx
    error=sigma/sqrt(real(nmc))         		!calculating erroe in calculation of integral
    
    write(*,*)"Ingragral, I=",integral
    write(*,*)"Error =",error
end program MultInte

function f(x1,x2,x3)
    implicit none
    real,intent(in)::x1,x2,x3
    real::f

    f=sin(x1*x2*x3)
end function f
