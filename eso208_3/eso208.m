err = 'Error';
name = input('file name : ', 's');
data = importdata(name);
x = data(:, 1);
y = data(:, 2);
n = length(x);
disp('What do you want to do?');
disp('1. Fit a least square polynomial');
disp('2. Fit a Lngrange Interpolation polynomial');
disp('3. Fit a Newtons Interpolation Polynomial');
disp('4. Fit Cubic Splines');
choice = input('Make your choice : ');
switch choice
  case 1
    leastsquare(n, x, y);
  case 2
    lagrangeinterpolation(n, x, y);
  case 3
    newtonsinterpolation(n, x, y);
  case 4
    cubicsplines(n, x, y);
end

function leastsquare(n, x, y)
  m = input('order of the polynomial (<n): ');
  A = zeros(m);
  h(1:2*m-1) = 0; % our hero
  b(1:m) = 0; y_bar = 0;
  for i = 1:n
    y_bar = y_bar + y(i);
    for j = 1:m
      b(j) = b(j) + y(i)*x(i)^(j-1);
    end
    for j = 1:2*m-1
      h(j) = h(j) + x(i)^(j-1);
    end
  end
  for i = 1:m
    for j = 1:m
      A(i, j) = h(i+j-1);
    end
  end
  a = inv(A)*b';
  for i = 1:m
    b(i) = a(m-i+1);
  end
  q(1:n) = 0;  y_bar = y_bar / n;
  S = 0; S0 = 0;
  for i = 1:n
    q(i) = b(1);
    for j = 2:m;
      q(i) = q(i)*x(i) + b(j);
    end
    S = S + (y(i)-q(i))^2;
    S0 = S0 + (y(i)-y_bar)^2;
  end
  savePoly(b, "LeastSquare");
  plotLN(n, x, y, b, 'Least Square Polynomial');
  f = fopen('outLeastSquare.txt', 'a');
  fprintf(f,"\n\nR-sq = %f\n\n", 1 - S/S0);
  fclose(f);
end

function lagrangeinterpolation(n, x, y)
  a(1:n) = 0;
  for i = 1:n
    p = 1;
    q = 1;
    for j = 1:n
      if j ~= i
        p = conv(p, [-x(j), 1]);
        q = q * (x(i) - x(j));
      end
    end
    for j = 1:n
      a(n-j+1) = a(n-j+1) + p(j)*y(i)/q;
    end
  end
  savePoly(a, "Lagrange");
  plotLN(n, x, y, a, 'Lagrange Interpolation');
end
function newtonsinterpolation(n, x, y)
  c(1:n) = 0;
  a(1:n) = 0;
  for i = 1:n
    c(i) = nthDD(y, x, 1, i);
  end
  a(n) = c(1);
  for i = 2:n
    p = [-x(1);1];
    for j = 3:i
      p = conv(p,[-x(j-1);1]);
    end
    for j = 1:i
      a(n-j+1) = a(n-j+1) + c(i)*p(j);
    end
  end
  savePoly(a, 'Newton');
  plotLN(n, x, y, a, 'Newtons Interpolation');
end

function y = nthDD(f, x, i, j)
    if (i == j)
        y = f(i);
    else
        y = (nthDD(f,x,i+1,j) - nthDD(f,x,i,j-1)) / (x(j) - x(i));
    end
end

function savePoly(a, str)
  f = fopen("out"+str+".txt", 'w');
  fprintf(f, '%s Interpolation Polynomial.....\n\n', str);
  fprintf(f,'Coefficients of the Polynomial : \n\n');
  fprintf(f, '%f  ', a);
  fclose(f);
end

function plotLN(n, x, y, a, str)
  scatter(x, y, 'filled');
  hold on
  p = x(1):0.01:x(n);
  m = length(p);
  q(1:m) = 0;
  l = length(a);
  for i = 1:m
    q(i) = a(1);
    for j = 2:l;
      q(i) = q(i)*p(i) + a(j);
    end
  end
  plot(p, q);
  title(str);
  hold off
end
function cubicsplines(n, x, y)
  disp('Choose your spline : ');
  disp('1. Linear Spline');
  disp('2. Quadratic Spline');
  disp('3. Natural Spline');
  disp('4. Not-a-knot');
  disp('5. Periodic');
  disp('6. Clamped Spline');
  choice = input('Make your choice: ');
  switch choice
    case 1
      linearspline(n, x, y);
    case 2
      quadraticspline(n, x, y);
    case 3
      naturalspline(n, x, y);
    case 4
      notaknotspline(n, x, y);
    case 5
      periodicspline(n, x, y);
    case 6
      clampedspline(n, x, y);
  end
end

function linearspline(n, x, y)
  c0=1:n-1; c1=1:n-1;
  f = fopen('outLinearSpline.txt','w');
  fprintf(f, '+++ Output of Linear Spline Interpolation on given data set. +++');
  fprintf(f, "\n\nCoefficients of the equations in tabular form :\n\n");
  fprintf(f, "i         c1            c0 \n\n");
  scatter(x, y, 'filled');
  hold on
  for i = 1:n-1
    c1(i)=nthDD(y, x, i, i+1); c0(i)=y(i)-c1(i)*x(i);
    fprintf(f, "%d   %f    %f\n", i, c1(i), c0(i));
    p = x(i):0.01:x(i+1);
    q = c1(i).*p + c0(i);
    plot(p, q);
  end
  hold off
  fprintf(f, "\n\nInterpolated values y* at give x* :\n\n");
  fprintf(f, "   x*        y*\n\n");
  m = input('Values of x : ');
  t = 0;
  for i = 1:length(m)
    for j = 2:n
      if m(i)<x(j)
        t = j-1; break;
      end
    end
    fprintf(f, "%.3f     %f\n",m(i),m(i)*c1(t) + c0(t));
  end
  fclose(f);
end
function quadraticspline(n, x, y)
  u(1:n) = 0;
  u(1)= input('Enter the slope of the first node: ');
  c0=1:n-1; c1=1:n-1; c2=1:n-1;
  f = fopen('outQuadraticSpline.txt','w');
  fprintf(f, '+++ Output of Quadratic Spline Interpolation on given data set. +++');
  fprintf(f, "\n\nCoefficients of the equations in tabular form :\n\n");
  fprintf(f, "i         c2             c1            c0 \n\n");
  scatter(x, y, 'filled');
  hold on
  for i =1:n-1
    u(i+1) = 2*nthDD(y,x,i,i+1) - u(i);
    h = x(i+1) - x(i);
    c0(i) = u(i+1)*x(i)^2 - u(i)*x(i+1)^2 + y(i) + u(i)*h/2;
    c1(i) = u(i)*x(i+1)/h - u(i+1)*x(i)/h;
    c2(i) = (u(i+1)-u(i))/(2*h);
    fprintf(f,"%d     %f     %f      %f\n",i,c2(i),c1(i),c0(i));
    p = x(i):0.01:x(i+1);
    q = c2(i).*(p.^2) + c1(i).*p + c0(i);
    plot(p, q);
  end
  hold off
  fprintf(f, "\n\nInterpolated values y* at give x* :\n\n");
  fprintf(f, "   x*        y*\n\n");
  m = input('Values of x : ');
  t = 0;
  for i = 1:length(m)
    for j = 2:n
      if m(i)<x(j)
        t = j-1; break;
      end
    end
    fprintf(f, "%.3f     %f\n",m(i),m(i)^2*c2(t) + m(i)*c1(t) + c0(t));
  end
  fclose(f);
end
function naturalspline(n, x, y)
  h(1:n) = 0;
  g(1:n) = 0;
  a = zeros(n);
  b(1:n) = 0;
  A(1:n-1)=0; B(1:n-1)=0; C(1:n-1)=0; D(1:n-1)=0;
  c0(1:n-1)=0; c1(1:n-1)=0; c2(1:n-1)=0; c3(1:n-1)=0;
  for i = 2:n
    h(i) = x(i)-x(i-1);
    g(i) = (y(i)-y(i-1))/h(i);
  end
  b(1) = 0;
  b(n) = 0;
  for i = 2:n-1
    b(i) = 6*(g(i+1)-g(i));
    a(i, i-1) = h(i);
    a(i, i+1) = h(i+1);
    a(i, i) = 2*(h(i)+h(i+1));
  end
  a(1,1) = 1;
  a(n,n) = 1;
  alpha = inv(a)*b';
  f = fopen('outNaturalSpline.txt','w');
  fprintf(f, "Coefficients of the equations in tabular form :\n\n");
  fprintf(f, "i       c3         c2         c1         c0 \n\n");
  for i = 1:n-1
    A(i) = alpha(i+1)/(6*h(i+1));
    B(i) = alpha(i)/(6*h(i+1));
    C(i) = y(i+1)/h(i+1) - alpha(i+1)*h(i+1)/6;
    D(i) = y(i)/h(i+1) - alpha(i)*h(i+1)/6;
    c3(i) = A(i) - B(i);
    c2(i) = 3*x(i+1)*B(i) - 3*x(i)*A(i);
    c1(i) = 3*x(i)^2*A(i) - 3*x(i+1)^2*B(i) + C(i) - D(i);
    c0(i) = B(i)*x(i+1)^3 - A(i)*x(i)^3 - C(i)*x(i) + D(i)*x(i+1);
    fprintf(f, "%d   %f   %f   %f   %f\n",i,c3(i),c2(i),c1(i),c0(i));
  end
  fprintf(f, "\n\n1st derivative and 2nd derivative at each node :\n\n");
  fprintf(f, "i     1st der       2nd der\n\n");
  for i = 1:n-1
    fprintf(f, "%d   %f   %f\n",i-1,3*c3(i)*x(i)^2+2*c2(i)*x(i)+c1(i), 6*c3(i)*x(i)+2*c2(i));
  end
  fprintf(f, "%d   %f   %f\n",n-1,3*c3(n-1)*x(n)^2+2*c2(n-1)*x(n)+c1(n-1), 6*c3(n-1)*x(n)+2*c2(n-1));
  fprintf(f, "\n\nInterpolated values y* at give x* :\n\n");
  fprintf(f, "   x*        y*\n\n");
  m = input('Values of x : ');
  t = 0;
  for i = 1:length(m)
    for j = 2:n
      if m(i)<x(j)
        t = j-1; break;
      end
    end
    fprintf(f, "%.3f     %f\n",m(i),m(i)^3*c3(t) + m(i)^2*c2(t) + m(i)*c1(t) + c0(t));
  end
  fclose(f);
  scatter(x, y, 'filled');
  hold on
  for i = 1:n-1
    p = x(i):0.01:x(i+1);
    q = p.^3*c3(i) + p.^2*c2(i) + p*c1(i) + c0(i);
    plot(p, q);
  end
  hold off
end
function notaknotspline(n, x, y)
  h(1:n) = 0;
  g(1:n) = 0;
  a = zeros(n);
  b(1:n) = 0;
  A(1:n-1)=0; B(1:n-1)=0; C(1:n-1)=0; D(1:n-1)=0;
  c0(1:n-1)=0; c1(1:n-1)=0; c2(1:n-1)=0; c3(1:n-1)=0;
  for i = 2:n
    h(i) = x(i)-x(i-1);
    g(i) = (y(i)-y(i-1))/h(i);
  end
  b(1) = 0;
  b(n) = 0;
  for i = 2:n-1
    b(i) = 6*(g(i+1)-g(i));
    a(i, i-1) = h(i);
    a(i, i+1) = h(i+1);
    a(i, i) = 2*(h(i)+h(i+1));
  end
  a(1,1) = h(3)/(h(2)+h(3)); a(1,2) = -1; a(1,3) = h(2)/(h(2)+h(3));
  a(n,n) = h(n-1)/(h(n)+h(n-1)); a(n,n-1) = -1; a(n,n-2) = h(n)/(h(n)+h(n-1));
  alpha = inv(a)*b';
  f = fopen('outnotaknotSpline.txt','w');
  fprintf(f, "Coefficients of the equations in tabular form :\n\n");
  fprintf(f, "i       c3         c2         c1         c0 \n\n");
  for i = 1:n-1
    A(i) = alpha(i+1)/(6*h(i+1));
    B(i) = alpha(i)/(6*h(i+1));
    C(i) = y(i+1)/h(i+1) - alpha(i+1)*h(i+1)/6;
    D(i) = y(i)/h(i+1) - alpha(i)*h(i+1)/6;
    c3(i) = A(i) - B(i);
    c2(i) = 3*x(i+1)*B(i) - 3*x(i)*A(i);
    c1(i) = 3*x(i)^2*A(i) - 3*x(i+1)^2*B(i) + C(i) - D(i);
    c0(i) = B(i)*x(i+1)^3 - A(i)*x(i)^3 - C(i)*x(i) + D(i)*x(i+1);
    fprintf(f, "%d   %f   %f   %f   %f\n",i,c3(i),c2(i),c1(i),c0(i));
  end
  fprintf(f, "\n\n1st derivative and 2nd derivative at each node :\n\n");
  fprintf(f, "i     1st der       2nd der\n\n");
  for i = 1:n-1
    fprintf(f, "%d   %f   %f\n",i-1,3*c3(i)*x(i)^2+2*c2(i)*x(i)+c1(i), 6*c3(i)*x(i)+2*c2(i));
  end
  fprintf(f, "%d   %f   %f\n",n-1,3*c3(n-1)*x(n)^2+2*c2(n-1)*x(n)+c1(n-1), 6*c3(n-1)*x(n)+2*c2(n-1));
  fprintf(f, "\n\nInterpolated values y* at give x* :\n\n");
  fprintf(f, "   x*        y*\n\n");
  m = input('Values of x : ');
  t = 0;
  for i = 1:length(m)
    for j = 2:n
      if m(i)<x(j)
        t = j-1; break;
      end
    end
    fprintf(f, "%.3f     %f\n",m(i),m(i)^3*c3(t) + m(i)^2*c2(t) + m(i)*c1(t) + c0(t));
  end
  fclose(f);
  scatter(x, y, 'filled');
  hold on
  for i = 1:n-1
    p = x(i):0.01:x(i+1);
    q = p.^3*c3(i) + p.^2*c2(i) + p*c1(i) + c0(i);
    plot(p, q);
  end
  hold off
end
function clampedspline(n, x, y)
  slope1=input('Enter slope at first point: ');
  slope2=input('Enter slope at last point: ');
  h(1:n) = 0;
  g(1:n) = 0;
  a = zeros(n);
  b(1:n) = 0;
  A(1:n-1)=0; B(1:n-1)=0; C(1:n-1)=0; D(1:n-1)=0;
  c0(1:n-1)=0; c1(1:n-1)=0; c2(1:n-1)=0; c3(1:n-1)=0;
  for i = 2:n
    h(i) = x(i)-x(i-1);
    g(i) = (y(i)-y(i-1))/h(i);
  end
  b(1)=((g(2)-slope1)*6)/(x(2)-x(1));
  b(n)=((slope2-g(n))*6)/(h(n));
  for i = 2:n-1
    b(i) = 6*(g(i+1)-g(i));
    a(i, i-1) = h(i);
    a(i, i+1) = h(i+1);
    a(i, i) = 2*(h(i)+h(i+1));
  end
  a(1,1) = 2; a(1,2) =1;
  a(n,n) = 2; a(n,n-1) =1;
  alpha = inv(a)*b';
  f = fopen('outclampedSpline.txt','w');
  fprintf(f, "Coefficients of the equations in tabular form :\n\n");
  fprintf(f, "i       c3         c2         c1         c0 \n\n");
  for i = 1:n-1
    A(i) = alpha(i+1)/(6*h(i+1));
    B(i) = alpha(i)/(6*h(i+1));
    C(i) = y(i+1)/h(i+1) - alpha(i+1)*h(i+1)/6;
    D(i) = y(i)/h(i+1) - alpha(i)*h(i+1)/6;
    c3(i) = A(i) - B(i);
    c2(i) = 3*x(i+1)*B(i) - 3*x(i)*A(i);
    c1(i) = 3*x(i)^2*A(i) - 3*x(i+1)^2*B(i) + C(i) - D(i);
    c0(i) = B(i)*x(i+1)^3 - A(i)*x(i)^3 - C(i)*x(i) + D(i)*x(i+1);
    fprintf(f, "%d   %f   %f   %f   %f\n",i,c3(i),c2(i),c1(i),c0(i));
  end
  fprintf(f, "\n\n1st derivative and 2nd derivative at each node :\n\n");
  fprintf(f, "i     1st der       2nd der\n\n");
  for i = 1:n-1
    fprintf(f, "%d   %f   %f\n",i-1,3*c3(i)*x(i)^2+2*c2(i)*x(i)+c1(i), 6*c3(i)*x(i)+2*c2(i));
  end
  fprintf(f, "%d   %f   %f\n",n-1,3*c3(n-1)*x(n)^2+2*c2(n-1)*x(n)+c1(n-1), 6*c3(n-1)*x(n)+2*c2(n-1));
  fprintf(f, "\n\nInterpolated values y* at give x* :\n\n");
  fprintf(f, "   x*        y*\n\n");
  m = input('Values of x : ');
  t = 0;
  for i = 1:length(m)
    for j = 2:n
      if m(i)<x(j)
        t = j-1; break;
      end
    end
    fprintf(f, "%.3f     %f\n",m(i),m(i)^3*c3(t) + m(i)^2*c2(t) + m(i)*c1(t) + c0(t));
  end
  fclose(f);
  scatter(x, y, 'filled');
  hold on
  for i = 1:n-1
    p = x(i):0.01:x(i+1);
    q = p.^3*c3(i) + p.^2*c2(i) + p*c1(i) + c0(i);
    plot(p, q);
  end
  hold off
end
function periodicspline(n, x, y)
  h(1:n) = 0;
  g(1:n) = 0;
  a = zeros(n);
  b(1:n) = 0;
  A(1:n-1)=0; B(1:n-1)=0; C(1:n-1)=0; D(1:n-1)=0;
  c0(1:n-1)=0; c1(1:n-1)=0; c2(1:n-1)=0; c3(1:n-1)=0;
  for i = 2:n
    h(i) = x(i)-x(i-1);
    g(i) = (y(i)-y(i-1))/h(i);
  end
  b(1)=0;
  b(n)=0;
  for i = 2:n-1
    b(i) = 6*(g(i+1)-g(i));
    a(i, i-1) = h(i);
    a(i, i+1) = h(i+1);
    a(i, i) = 2*(h(i)+h(i+1));
  end
  a(1,1)=2*h(1);
  a(1,2)=h(1);
  a(1,n-1)=h(n-1);
  a(1,n)=2*h(n-1);
  a(n,1)=1;
  a(n,n)=-1;
  alpha = inv(a)*b';
  f = fopen('outperiodicSpline.txt','w');
  fprintf(f, "Coefficients of the equations in tabular form :\n\n");
  fprintf(f, "i       c3         c2         c1         c0 \n\n");
  for i = 1:n-1
    A(i) = alpha(i+1)/(6*h(i+1));
    B(i) = alpha(i)/(6*h(i+1));
    C(i) = y(i+1)/h(i+1) - alpha(i+1)*h(i+1)/6;
    D(i) = y(i)/h(i+1) - alpha(i)*h(i+1)/6;
    c3(i) = A(i) - B(i);
    c2(i) = 3*x(i+1)*B(i) - 3*x(i)*A(i);
    c1(i) = 3*x(i)^2*A(i) - 3*x(i+1)^2*B(i) + C(i) - D(i);
    c0(i) = B(i)*x(i+1)^3 - A(i)*x(i)^3 - C(i)*x(i) + D(i)*x(i+1);
    fprintf(f, "%d   %f   %f   %f   %f\n",i,c3(i),c2(i),c1(i),c0(i));
  end
  fprintf(f, "\n\n1st derivative and 2nd derivative at each node :\n\n");
  fprintf(f, "i     1st der       2nd der\n\n");
  for i = 1:n-1
    fprintf(f, "%d   %f   %f\n",i-1,3*c3(i)*x(i)^2+2*c2(i)*x(i)+c1(i), 6*c3(i)*x(i)+2*c2(i));
  end
  fprintf(f, "%d   %f   %f\n",n-1,3*c3(n-1)*x(n)^2+2*c2(n-1)*x(n)+c1(n-1), 6*c3(n-1)*x(n)+2*c2(n-1));
  fprintf(f, "\n\nInterpolated values y* at give x* :\n\n");
  fprintf(f, "   x*        y*\n\n");
  m = input('Values of x : ');
  t = 0;
  for i = 1:length(m)
    for j = 2:n
      if m(i)<x(j)
        t = j-1; break;
      end
    end
    fprintf(f, "%.3f     %f\n",m(i),m(i)^3*c3(t) + m(i)^2*c2(t) + m(i)*c1(t) + c0(t));
  end
  fclose(f);
  scatter(x, y, 'filled');
  hold on
  for i = 1:n-1
    p = x(i):0.01:x(i+1);
    q = p.^3*c3(i) + p.^2*c2(i) + p*c1(i) + c0(i);
    plot(p, q);
  end
  hold off
end