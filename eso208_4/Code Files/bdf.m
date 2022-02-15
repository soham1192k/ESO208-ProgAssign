function [] = bdf(f, x0, y0, xf, h)
%BDF Summary of this function goes here
%   Detailed explanation goes here
interval = (xf - x0)/h;
x(1:interval+3) = 0;
y(1:interval+3) = 0;
x(1) = x0;
y(1) = y0;
x(interval+2) = xf+1;

%bdf 1st order
y(2) = y(1) + h* f(x(1),y(1)); %predictor
x(2) = x(1)+h;

err = 1;
    while err > 0.0000001
        y_old = y(2);
        y(2) = y(1) + h*( f(x(2),y(2)) ); %corrector
        err  = abs((y_old-y(2))/y(2));
    end
    
%bdf 2nd order
y(3) = y(2) + h* f(x(2),y(2)); %predictor
x(3) = x(2)+h;

err = 1;
    while err > 0.0000001
        y_old = y(3);
        y(3) = 2/3*( 2*y(2)-1/2*y(1) + h*( f(x(3),y(3)) ) ); %corrector
        err  = abs((y_old-y(3))/y(3));
    end
   
%bdf 3rd order
y(4) = y(3) + h* f(x(3),y(3)); %predictor
x(4) = x(3)+h;

err = 1;
    while err > 0.0000001
        y_old = y(4);
        y(4) = 6/11*( 3*y(3)-3/2*y(2)+1/3*y(1) + h*( f(x(4),y(4)) ) ); %corrector
        err  = abs((y_old-y(4))/y(4));
    end

%bdf 4th order    
i=4;
while x(i) <= xf
    y(i+1) = y(i) + h* f(x(i),y(i)); %predictor
    
    x(i+1) = x(i)+h;
    err = 1;
    while err > 0.0000001
        y_old = y(i+1);
        y(i+1) = 12/25*( 4*y(i)-3*y(i-1)+4/3*y(i-2)-1/4*y(i-3)+h*f(x(i+1),y(i+1)) ); %corrector
        err  = abs((y_old-y(i+1))/y(i+1));
    end
    i = i+1;
    
end
%x(i)=0;
%y(i)=0;
f=fopen('outBDF.txt','w');
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
