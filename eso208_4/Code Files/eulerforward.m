function [] = eulerforward(f, x0, y0, xf, h)
%EULERFORWARD Summary of this function goes here
%   Detailed explanation goes here
interval = (xf - x0)/h;
x(1:interval+3) = 0;
y(1:interval+3) = 0;
x(1) = x0;
y(1) = y0;

i=1;
while x(i) <= xf
    y(i+1) = y(i) + h* f(x(i),y(i));
    x(i+1) = x(i)+h;
    i = i+1;
    
end

  f = fopen("out_euler_forward.txt", 'w');
  fprintf(f,'x           y\n');
  
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
