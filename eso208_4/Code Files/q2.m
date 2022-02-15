lower=0;
upper=2;
h=input('Enter the grid size : ');
choice = input('Enter a for second order backward difference or b for ghost node : ','s');
x=lower:h:upper;
y=lower:h:upper;
len=length(x);
a=1:len;
b=1:len;
f=1:len;
for i=1:len
    a(i)=(-1*(x(i)+3))/(x(i)+1);
    b(i)=(x(i)+3)/((x(i)+1)*(x(i)+1));
    f(i)=(2*(x(i)+1))+3*b(i);
end
n=(upper-lower)/h;
if choice=='a'
    A=zeros(n-1,n-1);
    B=zeros(n-1,1);
    alpha=1:n-1;
    beta=1:n-1;
    gamma=1:n-1;
    for i=1:n-1
        alpha(i)= (2-h*a(i+1))/(2*h*h);
        beta(i)=(-2+(h*h*b(i+1)))/(h*h);
        gamma(i)= (2+h*a(i+1))/(2*h*h);
    end
    A(1,1)=beta(1);
    A(1,2)=gamma(1);
    index=2;
    B(1,1)=f(2)-(5*alpha(1));
    for i=2:n-2
        A(i,index-1)=alpha(i);
        A(i,index)=beta(i);
        A(i,index+1)=gamma(i);
        B(i,1)=f(i+1);
        index=index+1;
    end
    B(n-1,1)=f(n);    
    A(n-1,n-1)=beta(n-1)+((4/3)*gamma(n-1));
    A(n-1,n-2)=alpha(n-1) - (gamma(n-1)/3);
    soln=inv(A)*B;
    fprintf("***********************\n");
    fprintf("x           Temperature\n");
    fprintf("0             5\n");
    for i=2:n
        fprintf("%f      %f\n",x(i),soln(i-1));
    end
    fprintf("%f      %f\n",x(n+1),(4*soln(n-1)-soln(n-2))/3);
    fprintf("***********************\n");
    y(1)=5;
    for i=2:len-1
        y(i)=soln(i-1);
    end
    y(len)=(4*soln(n-1)-soln(n-2))/3;
    
    figure
    plot(x,y)
    xlabel('Length')
    ylabel('Temperature')
elseif choice=='b'
    A=zeros(n,n);
    B=zeros(n,1);
    alpha=1:n;
    beta=1:n;
    gamma=1:n;
    for i=1:n
        alpha(i)= (2-h*a(i+1))/(2*h*h);
        beta(i)=(-2+(h*h*b(i+1)))/(h*h);
        gamma(i)= (2+h*a(i+1))/(2*h*h);
    end
    A(1,1)=beta(1);
    A(1,2)=gamma(1);
    index=2;
    B(1,1)=f(2)-(5*alpha(1));
    for i=2:n-1
        A(i,index-1)=alpha(i);
        A(i,index)=beta(i);
        A(i,index+1)=gamma(i);
        B(i,1)=f(i+1);
        index=index+1;
    end
    B(n,1)=f(n+1);
    A(n,n)=beta(n);
    A(n,n-1)= alpha(n)+gamma(n);
    soln=inv(A)*B;
    fprintf("***********************\n");
    fprintf("x       Temperature\n");
    fprintf("0                5\n");
    for i=2:n+1
        fprintf("%f      %f\n",x(i),soln(i-1));
    end
    fprintf("***********************\n");
    y(1)=5;
    for i=2:len
        y(i)=soln(i-1);
    end
    figure
    plot(x,y)
    xlabel('Length')
    ylabel('Temperature')
end