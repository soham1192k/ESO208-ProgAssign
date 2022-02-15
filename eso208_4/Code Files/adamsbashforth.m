function [] = adamsbashforth(f, x0, y0, xf, h)
%ADAMSBASHFORTH Summary of this function goes here
%   Detailed explanation goes here
interval = (xf - x0)/h;
x(1:interval+3) = 0;
y(1:interval+3) = 0;
x(1) = x0;
y(1) = y0;
x(interval+2) = xf+1;

%EF
y(2)= y(1)+ h*f(x(1),y(1));

%AB2
x(2)=x(1)+h;
y(3)= y(2)+h*((1.5*f(x(2),y(2)) - 0.5*f(x(1),y(1)) ));

%AB3
x(3)=x(2)+h;
y(4)=y(3)+h*((23/12)*f(x(3),y(3)) - (16/12)*f(x(2),y(2)) +(5/12)*f(x(1),y(1)));

%AB4
x(4)=x(3)+h;
for i=5:(xf-x0)/h + 1
   y(i)=y(i-1)+h*((55/24)*f(x(i-1),y(i-1))-(59/24)*f(x(i-2),y(i-2))+(37/24)*f(x(i-3),y(i-3)) -(9/24)*f(x(i-4),y(i-4)));
   x(i) = x(i-1)+h;
end

f=fopen('outAdamsBashforth_4th_order.txt','w');
fprintf(f,"x\t\t\ty\n");
  for j=1:interval+1
    fprintf(f,'%f    %f\n', x(j),y(j));
  end
fclose(f);
X(1:interval+1) = x(1:interval+1);
Y(1:interval+1) = y(1:interval+1);
figure
plot(X,Y,'-o');
xlabel('x')
ylabel('y')
end
