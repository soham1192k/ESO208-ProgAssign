function [] = eulerbackward(f, x0, y0, xf, h)
%EULERBACKWARD Summary of this function goes here
%   Detailed explanation goes here
interval = (xf - x0)/h;
x(1:interval+3) = 0;
y(1:interval+3) = 0;
x(1) = x0;
y(1) = y0;
x(interval+2) = xf+1;
i=1;
while x(i) <= xf
    y(i+1) = y(i) + h* f(x(i),y(i)); %predictor-euler forward
    
    x(i+1) = x(i)+h;
    err = 1;
    while err > 0.0000001
        y_old = y(i+1);
        y(i+1) = y(i) + h*( f(x(i+1),y(i+1)) ); %corrector
        err  = abs((y_old-y(i+1))/y(i+1));
    end
    i = i+1;
    
end
%x(i)=0;
%y(i)=0;
f=fopen('outEulerBackward.txt','w');
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
