function [] = rungekutta(f, x0, y0, xf, h)
%RUNGEKUTTA Summary of this function goes here
%   Detailed explanation goes here
interval = (xf - x0)/h;
x(1:interval+3) = 0;
y(1:interval+3) = 0;
x(1) = x0;
y(1) = y0;
x(interval+2) = xf+1;
i=1;
while x(i) <= xf
    phi0 = f(x(i),y(i));
    phi1 = f(x(i)+0.5*h,y(i)+0.5*h*phi0);
    phi2 = f(x(i)+0.5*h,y(i)+0.5*h*phi1);
    phi3 = f(x(i)+h,y(i)+h*phi2);
    y(i+1) = y(i) + h*( (1/6)*(phi0+phi3) + (1/3)*(phi1+phi2) ); %predictor
    
    x(i+1) = x(i)+h;
    
    i = i+1;
    
end
%x(i)=0;
%y(i)=0;
f=fopen('outRungeKutta.txt','w');
fprintf(f,"x\t\t\ty\n");
  for j=1:interval+1
    fprintf(f,'%f    %f\n', x(j),y(j));
  end
fclose(f);
X(1:interval+1) = x(1:interval+1);
Y(1:interval+1) = y(1:interval+1);
figure
plot(X,Y,'-o')
xlabel('x')
ylabel('y')
end
