function [] = adamsmoulton(f, x0, y0, xf, h)
%ADAMSMOULTON Summary of this function goes here
%   Detailed explanation goes here
interval = (xf - x0)/h;
x(1:interval+3) = 0;
y(1:interval+3) = 0;
x(1) = x0;
y(1) = y0;
x(interval+2) = xf+1;

%Trapezoidal
y(2) = y(1) + h* f(x(1),y(1)); %predictor- euler forward
    
    x(2) = x(1)+h;
    err = 1;
    while err > 0.0000001
        y_old = y(2);
        y(2) = y(1) + h/2*( f(x(2),y(2))+ f(x(1),y(1)) ); %corrector
        err  = abs((y_old-y(2))/y(2));
    end
    
%Adam-Moulton 3rd order
y(3) = y(2) + h* f(x(2),y(2)); %predictor- euler forward
    
    x(3) = x(2)+h;
    err = 1;
    while err > 0.0000001
        y_old = y(3);
        y(3) = y(2) + h*( 5/12*f(x(3),y(3))+ 2/3*f(x(2),y(2))-1/12*f(x(1),y(1)) ); %corrector
        err  = abs((y_old-y(3))/y(3));
    end
    

%Adam-Moulton 4th order
i=3;
while x(i) <= xf
    y(i+1) = y(i) + h* f(x(i),y(i)); %predictor- euler forward
    
    x(i+1) = x(i)+h;
    err = 1;
    while err > 0.0000001
        y_old = y(i+1);
        y(i+1) = y(i) + h*( 3/8*f(x(i+1),y(i+1))+19/24*f(x(i),y(i))-5/24*f(x(i-1),y(i-1))+1/24*f(x(i-2),y(i-2)) ); %corrector
        err  = abs((y_old-y(i+1))/y(i+1));
    end
    i = i+1;
    
end
%x(i)=0;
%y(i)=0;
f=fopen('outAdamsMoulton.txt','w');
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
