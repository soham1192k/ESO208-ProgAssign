y=100;count=0;relerr=100;
%taking the three stopping criteria as input
end1=input("Enter the convergence criteria for relative approximation in percentage(upper limit): ");
end2=input("Enter the convergence criteria for the function value(upper limit): ");
end3=input("Enter the maximum iteration number: ");
x=input("Is the equation a polynomial?(Y/N): ",'s');
%checking if the function is a polynomial
if x=='N'||x=='n'
    y=input("Press 1 for Bisection, 2 for False-position, 3 for Fixed point, 4 for newton-raphson and 5 for Secant method: ");

%Bisection
    if y==1
        f=input("Enter the function: ",'s');
        x1=input("Enter the lower starting point: ");
        x2=input("Enter the higher starting point: ");
        ff = inline(f,'x') ;
        val=[1:1:101];
        x=[1:1:101];
        for i=1:1:101
            x(i)=0 + (i-1)*(5)/100;
            val(i)=feval(ff,x(i));
        end
        figure(1)
        plot(x,val,'b');
        xlabel('x'); ylabel('Functional Value');
        while(((abs(relerr))>=end1)&&(abs(y)>end2)&&(count<end3))
            p=(x1+x2)/2;
            y = feval(ff,p) ;
            y1=feval(ff,x1);
            y2=feval(ff,x2);
            if y>0
                if y1>0&&y2<0
                    relerr=((p-x1)/p)*100;
                    x1=p;
                else
                    relerr=((p-x2)/p)*100;
                    x2=p;
                end
            else if y<0
                if y1<0&&y2>0
                    relerr=((p-x1)/p)*100;
                    x1=p;
                else
                    relerr=((p-x2)/p)*100;
                    x2=p;
                end
            else
                break;
            end
        count=count+1;
        iteration(count)=count-1;
        rel(count)=relerr;
        end
end
fprintf('\nThe root is:%0.9ld\n',p);
if abs(p)<=end2
    fprintf('Terminated due to condition 2\n');
elseif abs(relerr)<=end1
    fprintf('Terminated due to condition 1\n');
elseif(count>=end3)
    fprintf('Terminated due to condition 3\n');
end
figure(2)
plot(iteration,rel,'r');
xlim([1 count])
xlabel('Iteration number'); ylabel('Relative Error');

%Regula Falsi
elseif y==2
f=input("Enter the function: ",'s');
xa=input("Enter the lower starting point: ");
xb=input("Enter the higher starting point: ");
ff = inline(f,'x') ;
val=[1:1:101];
x=[1:1:101];
for i=1:1:101
    x(i)=0 + (i-1)*(5)/100;
    val(i)=feval(ff,x(i));
end
figure(1)
plot(x,val,'b');
xlabel('x'); ylabel('Functional Value');
while(((abs(relerr))>=end1)&&(abs(y)>end2)&&(count<end3))
    y1=feval(ff,xa);
    y2=feval(ff,xb);
    p=xa-((xb-xa)*y1/(y2-y1));
    y=feval(ff,p);
    if y>0
        if y1>0&&y2<0
            relerr=((p-xa)/p)*100;
            xa=p;
        elseif y1<0&&y2>0
            relerr=((p-xb)/p)*100;
            xb=p;
        end
    elseif y<0
        if y2<0&&y1>0
            relerr=((p-xb)/p)*100;
            xb=p;
        elseif y1<0&&y2>0
            relerr=((p-xa)/p)*100;
            xa=p;
        end
    else
        break;
    end
    count=count+1;
    iteration(count)=count;
    rel(count)=relerr;
end
fprintf('\nThe root is:%0.9ld\n',p);
if abs(y)<=end2
    fprintf('Terminated due to condition 2\n');
elseif abs(relerr)<=end1
    fprintf('Terminated due to condition 1\n');
elseif (count>=end3)
    fprintf('Terminated due to condition 3\n');
end
figure(2)
plot(iteration,rel,'r');
xlim([1 count])
xlabel('Iteration number'); ylabel('Relative Error');

%Fixed Point
elseif y==3
f=input("Enter the function g(x), such that f(x)=0 is expressed as x=g(x):",'s');
x0=input("Enter the initial starting point: ");
ff = inline(f,'x') ;
val=[1:1:101];
x=[1:1:101];
for i=1:1:101
x(i)=0 + (i-1)*(5)/100;
val(i)=feval(ff,x(i));
end
figure(1)
plot(x,val,'b');
xlabel('x'); ylabel('Functional Value');
while((abs(relerr)>=end1)&&(abs(y)>end2)&&(count<end3))
y = feval(ff,x0);
if y==0
break;
else
relerr=((y-x0)/y)*100;
x0=y;
count=count+1;
iteration(count)=count;
rel(count)=relerr;
end
end
fprintf('\nThe root is:%0.9ld\n',x0);
if abs(y)<=end2
fprintf('Terminated due to condition 2\n');
elseif abs(relerr)<=end1
fprintf('Terminated due to condition 1\n');
elseif(count>=end3)
fprintf('Terminated due to condition 3\n');
end
figure(2)
plot(iteration,rel,'r');
xlim([1 count])
xlabel('Iteration number'); ylabel('Relative Error');

%Newton-Raphson
elseif y==4
f=input("Enter the function: ",'s');
f1=input("Enter the differentiation of the function: ",'s');
x0=input("Enter the initial starting point: ");
ff = inline(f,'x') ;
gg=inline(f1,'x');
val=[1:1:101];
x=[1:1:101];
for i=1:1:101
x(i)=0 + (i-1)*(5)/100;
val(i)=feval(ff,x(i));
end
figure(1)
plot(x,val,'b');
xlabel('x'); ylabel('Functional Value');
while(((abs(relerr))>=end1)&&(abs(y)>end2)&&(count<end3))
y1=feval(ff,x0);
y2=feval(gg,x0);
p=x0-(y1)/(y2);
relerr=((p-x0)/p)*100;
x0=p;
count=count+1;
iteration(count)=count;
rel(count)=relerr;
end
fprintf('\nThe root is:%0.9ld\n',p);
if abs(y1)<=end2
fprintf('Terminated due to condition 2\n');
elseif abs(relerr)<=end1
fprintf('Terminated due to condition 1\n');
elseif(count>=end3)
fprintf('Terminated due to condition 3\n');
end
figure(2)
plot(iteration,rel,'r');
xlim([1 count])
xlabel('Iteration number'); ylabel('Relative Error');

%Secant
else
f=input("Enter the function: ",'s');
xl=input("Enter the lower starting point: ");
xu=input("Enter the higher starting point: ");
ff = inline(f,'x') ;
val=[1:1:101];
x=[1:1:101];
for i=1:1:101
x(i)=0 + (i-1)*(5)/100;
val(i)=feval(ff,x(i));
end
figure(1)
plot(x,val,'b');
xlabel('x'); ylabel('Functional Value');
while(((abs(relerr))>=end1)&&(abs(y)>end2)&&(count<end3))
y1=feval(ff,xl);
y2=feval(ff,xu);
p=xu-((xu-xl)*(y2))/(y2-y1);
y=feval(ff,p);
relerr=((p-xu)/p)*100;
xl=xu;
xu=p;
count=count+1;
iteration(count)=count;
rel(count)=relerr;
end
fprintf('\nThe root is:%0.9ld\n',p);
if abs(y)<=end2
fprintf('Terminated due to condition 2\n');
elseif abs(relerr)<=end1
fprintf('Terminated due to condition 1\n');
elseif(count>=end3)
fprintf('Terminated due to condition 3\n');
end
figure(2)
plot(iteration,rel,'r');
xlim([1 count])
xlabel('Iteration number'); ylabel('Relative Error');
end
else
order=input("Enter the order of the polynomial: ");
for c=1:1:(order+1)
coeffs(c)=input('Enter the coefficients, starting from the lowest degree: ');
end
y=input("Press 1 for Muller, and 2 for Bairstow method: ");

%Muller
if y==1
y0=0;y1=0;y2=0;
x0=input("Enter the lowest starting point: ");
x1=input("Enter another starting point: ");
x2=input("Enter the highest starting point: ");
val=[1:1:101];
x=[1:1:101];
for i=1:1:101
val(i)=0;
end
for i=1:1:101
x(i)= 0+ (i-1)*(5)/100;
for j=1:1:(order+1)
val(i)=val(i)+(coeffs(j)*(x(i).^(j-1)));
end
end
figure(1)
plot(x,val,'b');
xlabel('x'); ylabel('Functional Value');
while(((abs(relerr))>=end1)&&(abs(y)>end2)&&(count<end3))
h0=x1-x0;
h1=x2-x1;
for c=1:1:(order+1)
y0=y0+(coeffs(c)*(x0.^(c-1)));
y1=y1+(coeffs(c)*(x1.^(c-1)));
y2=y2+(coeffs(c)*(x2.^(c-1)));
end
del0=(y1-y0)/(x1-x0);
del1=(y2-y1)/(x2-x1);
a=(del1-del0)/(h1+h0);
b=(a*h1)+del1;
c=y2;
p1=x2-((2*c)/(b-sqrt(b*b-4*a*c)));
p2=x2-((2*c)/(b+sqrt(b*b-4*a*c)));
if p1<=p2
p=p1;
else
p=p2;
end
relerr=((p-x2)/p)*100;
x0=x1;
x1=x2;
x2=p;
count=count+1;
iteration(count)=count;
rel(count)=relerr;
y0=0;y1=0;y2=0;
end
fprintf('\nThe root is:%0.9ld\n',p);
g=0;
for c=1:1:(order+1)
g=g+(coeffs(c)*(p.^(c-1)));
end
if abs(g)<=end2
fprintf('Terminated due to condition 2\n');
elseif abs(relerr)<=end1
fprintf('Terminated due to condition 1\n');
else
fprintf('Terminated due to condition 3\n');
end
figure(2)
plot(iteration,rel,'r');
xlim([1 count])
xlabel('Iteration number'); ylabel('Relative Error');

%Bairstow
else
delr=0;dels=0;plotcount=2;
r=input("Enter the value of r: ");
s=input("Enter the value of s: ");
val=[1:1:101];
x=[1:1:101];
for i=1:1:101
val(i)=0;
end
for i=1:1:101
x(i)= 0+ (i-1)*(5)/100;
for j=1:1:(order+1)
val(i)=val(i)+(coeffs(j)*(x(i).^(j-1)));
end
end
figure(1)
plot(x,val,'b');
xlabel('x'); ylabel('Functional Value');
relerrr=100; relerrs=100;
for i=order+1:-1:1
a(i)=coeffs(i);
end%forloop
count=0;
while(order>=3)
relerrr=100; relerrs=100;count=0;
while(abs(relerrr)>=0.01 && abs(relerrs)>=0.01 && count<=end3)
for i=order+1:-1:1
if i==order+1
b(i)=a(i);
elseif i==order
b(i)=a(i)+r*b(i+1);
else
b(i)=a(i)+r*b(i+1)+s*b(i+2);
end%ifelse
end%forloop
for i=order+1:-1:1
if i==order+1
c(i)=b(i);
elseif i==order
c(i)=b(i)+r*c(i+1);
else
c(i)=b(i)+r*c(i+1)+s*c(i+2);
end%if else
end%forloop
delr=((b(2)*c(3))-(b(1)*c(4)))/(c(2)*c(4) - c(3)*c(3));
dels=((b(2)*c(2))-(b(1)*c(3)))/(c(3)*c(3) - c(2)*c(4));
relerrr=(delr/(r+delr))*100;
relerrs=(dels/(s+dels))*100;
r=r+delr;
s=s+dels;
count=count+1;
if(relerrr==0 && relerrs==0)
    iteration(count)=count;
    a1(count)=0;
    a2(count)=0;
    iteration(count+1)=count+1;
    a1(count+1)=0;
    a2(count+1)=0;
    break;
end
iteration(count)=count;
a1(count)=relerrr;
a2(count)=relerrs;
end%whileloop
root1=(r+sqrt(r*r+4*s))/2;
root2=(r-sqrt(r*r+4*s))/2;
if((r*r+4*s)<0)
root1real=real(root1);
root2real=real(root2);
root1imag=imag(root1);
root2imag=imag(root2);
fprintf('\n\nA root is:%0.9ld%+0.9ldi\n',root1real,root1imag);
fprintf('A root is:%0.9ld%+0.9ldi\n',root2real,root2imag);
else
fprintf('A root is:%0.9ld\n',root1);
fprintf('A root is:%0.9ld\n',root2);
end%if
if(relerrr==end2 && relerrs==end2)
    fprintf('Terminated due to condition 2\n');
elseif abs(relerrr)<=end1 || abs(relerrs)<=end1
fprintf('Terminated due to condition 1\n');
elseif count>=end3
fprintf('Terminated due to condition 3\n');
end
b=b(3:order+1);
a=b;
order=order-2;
figure(plotcount)
plot(iteration,a1,'b');
xlim([1 count+1])
xlabel('Iteration number');
ylabel('Relative Error in r');
plotcount=plotcount+1;
figure(plotcount)
plot(iteration,a2,'r');
xlim([1 count+1])
xlabel('Iteration number');
ylabel('Relative Error in s');
plotcount=plotcount+1;
end%outerwhile
if order==2
root1=(-b(2)+sqrt(b(2)*b(2)-4*b(3)*b(1)))/(2*b(3));
root2=(-b(2)-sqrt(b(2)*b(2)-4*b(3)*b(1)))/(2*b(3));
if ((b(2)*b(2)-4*b(3)*b(1)))<0
root1real=real(root1);
root2real=real(root2);
root1imag=imag(root1);
root2imag=imag(root2);
fprintf('A root is:%0.9ld%+0.9ldi\n',root1real,root1imag);
fprintf('A root is:%0.9ld%+0.9ldi\n',root2real,root2imag);
else
fprintf('A root is:%0.9ld\n',root1);
fprintf('A root is:%0.9ld\n',root2);
end
elseif order==1
root1= -b(1)/b(2);
fprintf('A root is:%0.9ld\n',root1);
end
end%main ifelse
end


