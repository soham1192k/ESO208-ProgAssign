
disp('This program is created to estimate eigenvalues.');
file = input('Enter file name for n, matrix A and % tolerance : ', 's');
data = importdata(file);
n = data(1);
A = zeros(n);
for i = 1:n*n
    A(i) = data(i+1);
end
A = A';
tol = data(n*n+2);
disp('Enter "L" for only the largest eigenvalue');
disp('Enter "A" for All eigenvalues');
choice = input('Your Response : ', 's');
if choice == 'L' || choice == 'l'
    powerM(n,A,tol);
elseif choice == 'A' || choice == 'a'
    QRdecomp(n,A,tol);
else
    disp(err);
end

function powerM(n,A,t)
    f = fopen('outPower.txt','w');
    fprintf(f,"Power Method\n\n");
    z = 1:n;
    z = z';
    y = 1:n;
    for i = 2:n
        y(i) = 0;
    end
    y = y';
    lbd1 = 0;
    lbd2 = 0;
    itr = 0;
    while true
        itr = itr + 1;
        z = A * y;
        y1 = z./norm(z);
        lbd1 = y' * z;
        if abs((lbd1-lbd2)/lbd1)*100 <= t
            break;
        end
        lbd2 = lbd1;
        y = y1;
    end
    fprintf(f,"The Largest eigenvalue : \n%f\n\n",lbd1);
    fprintf(f,"Eigenvector : \n");
    for i = 1:n
        fprintf(f,"%f\n",y(i));
    end
    fprintf(f, "\nIterations : \n%d\n", itr);
    disp("The largest eigenvalue : ");
    disp(lbd1);
    disp("Iterations : ");
    disp(itr);
    fclose(f);
    disp("--> The result has also been sent to file outPower.txt.");
end

function QRdecomp(n,A,t)
    f = fopen('outQRd.txt','w');
    fprintf(f,"QR Decomposition Method\n\n");
    q = zeros(n);
    z = 1:n;
    r = zeros(n);
    rold = r;
    itr = 0;
    while true
        itr= itr + 1;
        for i = 1:n
            for j = 1:n
                z(j) = A(j,i);
            end
            for j = 1:i-1
                p = 0;
                for k = 1:n
                    p = p + A(k,i)*q(k,j);
                end
                for k = 1:n
                    z(k) = z(k) - p*q(k,j);
                end
            end
            zz = norm(z);
            for j = 1:n
                q(j,i) = z(j) ./ zz;
            end
        end
        r = q' * A;
        A = r * q;
        % stopping criteria
        max = 0;
        p = 0;
        for i = 1:n
            if abs(r(i,i) - rold(i,i)) > max
                max = abs(r(i,i) - rold(i,i));
            end
            p = p + r(i,i)*r(i,i);
        end
        p = sqrt(p);
        if (max/p)*100 <= t
            break;
        end
        rold = r;
    end
    fprintf(f,"Eigenvalues : \n");
    for i = 1:n
        fprintf(f,"%f\n",r(i,i));
    end
    fprintf(f,"\nIterations : \n%d\n",itr);
    fclose(f);
    disp("--> The result has been sent to file outQRd.txt");
end
