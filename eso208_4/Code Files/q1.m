file= input('Enter the file name where the data is: ','s');
data=importdata(file);
syms f(x,y)
f(x,y)=str2sym(data.textdata{1,1});
x0=data.data(1);
y0=data.data(2);
xf=data.data(3);
h=data.data(4);

disp('What do you want to do?');
disp('1. Euler Forward');
disp('2. Euler Backward');
disp('3. Trapezoidal');
disp('4. 4th-order Adams-Bashforth');
disp('5. 4th-order Adams-Moulton');
disp('6. 4th-order Backward Difference Formulation (BDF)');
disp('7. 4th Order Runge-Kutta');
choice = input('Make your choice : ');


if choice == 1
    eulerforward(f, x0, y0, xf, h);
elseif choice == 2
    eulerbackward(f, x0, y0, xf, h);
elseif choice == 3
    trapezoidal(f, x0, y0, xf, h);
elseif choice == 4
    adamsbashforth(f, x0, y0, xf, h);
elseif choice == 5
    adamsmoulton(f, x0, y0, xf, h);
elseif choice == 6
    bdf(f, x0, y0, xf, h);
elseif choice == 7
    rungekutta(f, x0, y0, xf, h);
end

disp('Output Sent');