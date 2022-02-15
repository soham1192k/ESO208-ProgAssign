disp('This program is created to solve system of linear equations.');
disp('Choose one of these choices');
disp('1.) Solve a System of Equations');
disp('2.) Perform a LU decomposition');
disp('3.) Perform a Matrix Inversion');
choice = input('Enter a number 1-3 (default 1): ');
if choice < 1 || choice > 3
    choice = 1;
end
switch choice
    case 1
        str = input('Is the system tri-diagonal (Y/N) default N :', 's');
        if isempty(str)
            str = 'N';
        end
        if str=='y' || str=='Y'
            file1 = input('Enter the file name where the data is : ', 's');
            thomas(file1);
        elseif str=='N' || str=='n'
            file1 = input('Enter the file name where the data is : ', 's');
            gauss(file1);
        else 
            disp(err);
        end
    case 2
        str = input('Is the Matrix symmetric and positive definite (Y/N) default N : ', 's');
        if isempty(str)
            str = 'N';
        end
        if str=='y' || str=='Y'
            file1 = input('Enter the file name where the data is : ', 's');
            cholesky(file1);
        elseif str=='N' || str=='n'
            file1 = input('Enter the file name where the data is : ', 's');
            doocr(file1);
        else 
            disp(err);
        end
    case 3
        file1 = input('Enter the file name where the matrix is : ', 's');
        gaussjordon(file1);
    otherwise
        disp(err);
end

function  thomas(file)
    data = importdata(file);
    n = data(1);
    % ---------- Preallocate all required space for speed -----------
    l = 1:n-1;
    u = 1:n-1;
    d = 1:n;
    b = 1:n;
    alpha = 1:n;
    beta = 1:n;
    x = 1:n;
    % ------------- preallocation complete --------------------------
    for i = 1:n-1
        l(i) = data(1+i);
        u(i) = data(2*n+i);
    end
    for i = 1:n
        d(i) = data(n+i);
        b(i) = data(3*n-1+i);
    end
    % Algorithm 
    alpha(1) = d(1);
    beta(1) = b(1);
    for i = 2:n
        alpha(i) = d(i) - (l(i-1)/alpha(i-1))*u(i-1);
        beta(i) = b(i) - (l(i-1)/alpha(i-1))*beta(i-1);
    end
    x(n) = beta(n)/alpha(n);
    for i = 1:n-1
        x(n-i) = (beta(n-i) - u(n-i)*x(n-i+1))/alpha(n-i);
    end
    f = fopen('outThomas.txt', 'w');
    for i = 1:n
        fprintf(f, "X%d = %f\n", i, x(i));
    end
    fclose(f);
    disp("--> The result have been sent to file out_p1_1.txt.");
end

function gauss(file)
    data = importdata(file);
    n = data(1);
    % --------- Preallocate all required space for apeed ------------------
    A = zeros(n);
    b = 1:n;
    x = 1:n;
    % ----------------- Preallocation Complete -------------------------
    for i = 1:n*n
        A(i) = data(i+1);
    end
    for i = 1:n
        b(i) = data(1+n*n+i);
    end
    % A = A';
    % Algorithm 
    for i = 1:n-1
        % --------- partial pivoting starts ----------
        max = A(i,i);
        maxindex = i;
        for j = i+1:n
            if A(j,i) > max
                max = A(j,i);
                maxindex = j;
            end
        end
        if maxindex ~= i
            % swap row i with row j
            for k = 1:n
                A(i,k) = A(i,k) + A(j,k);
            end
            b(i) = b(i) + b(j);
            for k = 1:n
                A(j,k) = A(i,k) - A(j,k);
            end
            b(j) = b(i) - b(j);
            for k = 1:n
                A(i,k) = A(i,k) - A(j,k);
            end
            b(i) = b(i) - b(j);
        end
        for j = i+1:n
            l = A(j,i) / A(i,i);
            for k = i:n
                A(j,k) = A(j,k) - l*A(i, k);
            end
            b(j) = b(j) - l*b(i);
        end
    end
    % --------- partial pivoting ends ------------
    x(n) = b(n) / A(n, n);
    for i = 1:n-1
        x(n-i) = b(n-i);
        for j = 1:i
            x(n-i) = x(n-i) - A(n-i,n-j+1)*x(n-j+1);
        end
        x(n-i) = x(n-i) / A(n-i,n-i);
    end
    f = fopen('outGauss.txt', 'w');
    for i = 1:n
        fprintf(f, "X%d = %f\n", i, x(i));
    end
    fclose(f);
    disp("--> The result have been sent to file outGauss.txt.");
end

function cholesky(file)
    data = importdata(file);
    n = data(1);
    % ----------- Preallocate all required space for speed ------------
    A = zeros(n);
    b = 1:n;
    L = zeros(n);
    % ------------------- Preallocation Complete ---------------------
    for i = 1:n*n
        A(i) = data(i+1);
    end
    for i = 1:n
        b(i) = data(1+n*n+i);
    end
    A = A';
    % Algorithm
    f = fopen('outCholesky.txt','w');
    % ------------- Full Pivoting Starts --------------
    fprintf(f,"A (before pivoting) :\n");
    for i = 1:n
        for j = 1:n
            fprintf(f,"%f   ",A(i,j));
        end
        fprintf(f,"\n");
    end
    fprintf(f,"\nLog of all row and column swaps :\n");
    for i = 1:n-1
        % find max
        max = A(i,i);
        mx = i;
        my = i;
        for j = i:n
            if A(i,j) > max
                max = A(i,j);
                mx = i;
                my = j;
            end
            if A(j,n) > max
                max = A(j,n);
                mx = j;
                my = n;
            end
            if A(n,j) > max
                max = A(n,j);
                mx = n;
                my = j;
            end
            if A(j,i) > max
                max = A(j,i);
                mx = j;
                my = i;
            end
        end
        % swap column i with my
        if my ~= i
            for j = 1:n
                A(j,i) = A(j,i) + A(j,my);
            end
            for j = 1:n
                A(j,my) = A(j,i) - A(j,my);
            end
            for j = 1:n
                A(j,i) = A(j,i) - A(j,my);
            end
            fprintf(f,"Swap Column %d with Column %d\n", i, my);
        end
        % swap row i with mx
        if mx ~= i
            for j = 1:n
                A(i,j) = A(i,j) + A(mx,j);
            end
            for j = 1:n
                A(mx,j) = A(i,j) - A(mx,j);
            end
            for j = 1:n
                A(i,j) = A(i,j) - A(mx,j);
            end
            fprintf(f,"Swap Row %d with Row %d\n", i, mx);
        end
    end
    % ------------- Full Pivoting Ends --------------
    fprintf(f, "\nL\n");
    for i = 1:n
        for j = 1:i
            L(i,j) = A(i,j);
            if i==j
                for k = 1:j-1
                    L(i,j) = L(i,j) - L(j,k)*L(j,k);
                end
                L(i,j) = sqrt(L(i,j));
            else
                for k = i:j-1
                    L(i,j) = L(i,j) - L(i,k)*L(j,k);
                end
                L(i,j) = L(i,j) / L(j,j);
            end
            fprintf(f,"%f  ",L(i,j));
        end
        for j = i+1:n
            fprintf(f,"%f  ",L(i,j));
        end
        fprintf(f,"\n");
    end
    fclose(f);
    disp("--> The result have been sent to file outCholesky.txt");
end

function doocr(file)
    data = importdata(file);
    n = data(1);
    % ----------- Preallocate all required space for speed ------------
    A = zeros(n);
    b = 1:n;
    L  = zeros(n);
    U = zeros(n);
    % ------------------- Preallocation Complete ---------------------
    for i = 1:n*n
        A(i) = data(i+1);
    end
    for i = 1:n
        b(i) = data(1+n*n+i);
    end
    A = A';
    % Algorithm
    f = fopen('outDooCr.txt','w');
    % ------------- Full Pivoting Starts --------------
    fprintf(f,"A (before pivoting) :\n");
    for i = 1:n
        for j = 1:n
            fprintf(f,"%f   ",A(i,j));
        end
        fprintf(f,"\n");
    end
    fprintf(f,"\nLog of all row and column swaps :\n");
    for i = 1:n-1
        % find max
        max = A(i,i);
        mx = i;
        my = i;
        for j = i:n
            if A(i,j) > max
                max = A(i,j);
                mx = i;
                my = j;
            end
            if A(j,n) > max
                max = A(j,n);
                mx = j;
                my = n;
            end
            if A(n,j) > max
                max = A(n,j);
                mx = n;
                my = j;
            end
            if A(j,i) > max
                max = A(j,i);
                mx = j;
                my = i;
            end
        end
        % swap column i with my
        if my ~= i
            for j = 1:n
                A(j,i) = A(j,i) + A(j,my);
            end
            for j = 1:n
                A(j,my) = A(j,i) - A(j,my);
            end
            for j = 1:n
                A(j,i) = A(j,i) - A(j,my);
            end
            fprintf(f,"Swap Column %d with Column %d\n", i, my);
        end
        % swap row i with mx
        if mx ~= i
            for j = 1:n
                A(i,j) = A(i,j) + A(mx,j);
            end
            for j = 1:n
                A(mx,j) = A(i,j) - A(mx,j);
            end
            for j = 1:n
                A(i,j) = A(i,j) - A(mx,j);
            end
            fprintf(f,"Swap Row %d with Row %d\n", i, mx);
        end
    end
    % ------------- Full Pivoting Ends --------------
    disp("Choose one of these methods for LU Decomposition");
    disp('1.) Doolittle Decomposition');
    disp('2.) Crout Decomposition');
    choice = input('Your choice (default 1) : ');
    if isempty(choice)
        choice = 1;
    end
    if choice < 1 || choice > 2
        choice = 1;
    end
    if choice==1
        % doolittle(n,A,L,U);
        for i = 1:n
            L(i,i) = 1.0;
        end
        for i = 1:n
            for j = 1:i
                U(j,i) = A(j,i);
                for k = 1:j-1
                    U(j,i) = U(j,i) - L(j,k)*U(k,i);
                end
                if j ~= i
                    L(i,j) = A(i,j);
                    for k = 1:j-1
                        L(i,j) = L(i,j) - L(i,k)*U(k,j);
                    end
                    L(i,j) = L(i,j) / U(j,j);
                end
            end
        end
    elseif choice==2
        % crout(n,A,L,U);
        for i = 1:n
            U(i,i) = 1.0;
        end
        for i = 1:n
            for j = 1:i
                L(i,j) = A(i,j);
                for k = 1:j-1
                    L(i,j) = L(i,j) - L(i,k)*U(k,j);
                end
                if j ~= i
                    U(j,i) = A(j,i);
                    for k = 1:j-1
                        U(j,i) = U(j,i) - L(j,k)*U(k,i);
                    end
                    U(j,i) = U(j,i) / L(j,j);
                end
            end
        end
    else
        disp(err);
    end
    fprintf(f,"\nL\n");
    for i = 1:n
        for j = 1:n
            fprintf(f,"%f  ",L(i,j));
        end
        fprintf(f,"\n");
    end
    fprintf(f,"\nU\n");
    for i = 1:n
        for j = 1:n
            fprintf(f,"%f  ",U(i,j));
        end
        fprintf(f,"\n");
    end
    fclose(f);
    disp("--> The result have been sent to file outDooCr.txt");
end

function gaussjordon(file)
    data = importdata(file);
    n = data(1);
    % ----------- Preallocate all required space for speed ------------
    A = zeros(n);
    d  = zeros(n);
    % ------------------- Preallocation Complete ---------------------
    for i = 1:n*n
        A(i) = data(i+1);
    end
    A = A';
    for i = 1:n
        d(i,i) = 1.0;
    end
    for i = 1:n
        t = A(i,i);
        for j = 1:n
            d(i,j) = d(i,j) / t;
        end
        for j = 1:n
            A(i,j) = A(i,j) / t;
        end
        for j = 1:n
            if j ~= i
                t = A(j,i);
                for k = 1:n
                    A(j,k) = A(j,k) - A(i,k) * t;
                    d(j,k) = d(j,k) - d(i,k) * t;
                end
            end
        end
    end
    f = fopen('outGaussJordon.txt', 'w');
    fprintf(f, "Inverse of Given Matrix :\n\n");
    for i = 1:n
        for j = 1:n
            fprintf(f,"%f   ",d(i,j));
        end
        fprintf(f,"\n");
    end
    fclose(f);
    disp("--> The result have been sent to file outGaussJordon.txt.");
end