function [iter err] = testRGNEP(n,type,param,kappa,fname)
% test one step of randomized divide-and-conquer
% for the generalized eigenproblem

    if nargin < 5
        fname = 'noeps';
    end

    % create test problems
    if strcmp(type,'Sym'),
        % create random real eigenvalues between -1.5 and 1.5
        w = 3*(rand(1,n)-.5);
        B = diag(w);
    elseif strcmp(type,'Nonsym'),
        % create random matrix with eigenvalues inside the 
        % circle of radius 1.5
        B = randn(n);
        B = B / max(abs(eig(B))) * 2;
    elseif strcmp(type,'SmallReal'),
        % create block diagonal matrix with eigenvalues having 
        % small real part (equal to param) and random imag part
        B = zeros(n);
        for k = 1:2:n-1,
            imaginary = randn(1);
            sign = randn(1);
            if sign < 0,
                sign = -1;
            else
                sign = 1;
            end
            B(k,k) = sign*param;
            B(k+1,k) = -imaginary;
            B(k,k+1) = imaginary;
            B(k+1,k+1) = sign*param;				
        end
        if mod(n,2) == 1,
            B(end,end) = param;
        end    
    elseif strcmp(type,'Cond1'),
        % create example with all eigenvalues having specified real
        % part (param) and one evec with specified angle (kappa)
        B = zeros(n);
        for k = 1:2:n-1,
            imaginary = randn(1);
            sign = randn(1);
            if sign < 0,
                sign = -1;
            else
                sign = 1;
            end
            B(k,k) = sign*param;
            B(k+1,k) = -imaginary;
            B(k,k+1) = imaginary;
            B(k+1,k+1) = sign*param;				
        end
        if mod(n,2) == 1,
            B(end,end) = param;
        end
        angle = kappa;
        T = eye(n); T(n,n) = sin(angle);
        T(1:n-1,n) = sqrt(1/(n-1))*cos(angle)*ones(n-1,1);
        B = T*B/T;     
        % set kappa to 1 for random orthogonal similarity
        kappa = 1;
    elseif strcmp(type,'Cond2')
        % create example with one eigenvalue having specified real
        % part (param) and specified condition number (kappa) and all
        % other eigenvalues with real part +1 or -1
        B = zeros(n);
        for k = 1:2:n-1,
            imaginary = .5*randn(1);
            sign = randn(1);
            if sign < 0,
                sign = -1;
            else
                sign = 1;
            end
            B(k,k) = sign;
            B(k+1,k) = -imaginary;
            B(k,k+1) = imaginary;
            B(k+1,k+1) = sign;				
        end
        if mod(n,2) == 1,
            B(end,end) = param;
        else
            B(end-1,end-1) = -1;
            B(end,end-1) = 0; B(end-1,end) = 0;
            B(end,end) = param;
        end
        angle = kappa;
        T = eye(n); T(n,n) = sin(angle);
        T(1:n-1,n) = sqrt(1/(n-1))*cos(angle)*ones(n-1,1);
        B = T*B/T;     
        % set kappa to 1 for random orthogonal similarity
        kappa = 1;
    elseif strcmp(type,'Jordan'),
        % create example with a Jordan block centered at param
        % containing half the eigenvalues (other eigenvalues at -1)
        W = param*eye(n/2) + diag(ones(n/2-1,1),1);
        B = [W zeros(n/2); zeros(n/2) -eye(n/2)];
        kappa = 1;
    elseif strcmp(type,'special')
        for j=0:2:14,
            param = 10^(-j);
            B = zeros(n);
            for k = 1:2:n-1,
                imaginary = randn(1);
                sign = randn(1);
                if sign < 0,
                    sign = -1;
                else
                    sign = 1;
                end
                B(k,k) = sign*param;
                B(k+1,k) = -imaginary;
                B(k,k+1) = imaginary;
                B(k+1,k+1) = sign*param;				
            end
            if mod(n,2) == 1,
                B(end,end) = param;
            end
            [U R] = qr(randn(n));
            A=U*B*U';
            [Ap,Bp,bckerr,Rconv,splt,numit] = irs(A,eye(n),1,1,1,-1,60);
            semilogy(1:numit,bckerr(1:numit),'--or',1:numit,Rconv(1:numit),'-+b'), hold on
            title('Convergence')
            xlabel('IRS Iterations')
            ylabel('Backward Error')
            axis([1 60 10^(-17) 1])
        end
        hold off
        if ~strcmp(fname,'noeps')
            saveas(gcf,sprintf('%s.eps',fname),'epsc');
        end
        E=eig(A);
        neg_count = 0;
        for i=1:length(E),
            if real(E(i)) < 0,
                neg_count = neg_count + 1;
            end
        end
        if neg_count ~= n-splt(numit-1),
            splt(numit-1)
        end
        iter = Inf; err = Inf;
        return
    else
        error('invalid type');
    end
   % perform similarity transformation with eigenvector matrix having
   % condition number kappa
   [U R] = qr(randn(n));
   [V R] = qr(randn(n));
   s=ones(n,1);
   s(end) = 1/kappa;
   S = diag(s);
   iS = diag(1./s);
   A=U*(S*(V'*B*V)*iS)*U';
   %eA = sort(eig(A));
   %eB = sort(eig(B));
   %matlabs_max_eval_err_after_trans = max(abs(eA-eB)./abs(eA))
   %cond(U*S*V'*T)

    
    % perform implicit repeated squaring with imag axis as splitting line
    % (use identity for initial B matrix)
    [Ap,Bp,bckerr,Rconv,splt,numit] = irs(A,eye(n),1,1,1,-1,60);
    
    % complete the divide-and-conquer with GRURV
	[U,R1,V] = rurv(Ap);
	[R2,Qr] = rq(U'*(Ap+Bp));
	A2 = Qr*A*Qr';
	[E21norm,l] = split(A2);
    
    E = eig(A);
    E1 = eig(A2(1:l,1:l));
    E2 =  eig(A2(l+1:n,l+1:n));
    
    figure, plot(real(E),imag(E),'b*'), hold on
    t = linspace(-1.5,1.5);
    plot(1i*t,'--k'), hold off
    axis([-1.5 1.5 -1.5 1.5])
    axis square
    
    neg_count = 0;
    for i=1:length(E),
        if real(E(i)) < 0,
            neg_count = neg_count + 1;
        end
    end
    if neg_count ~= length(E2),
        err = Inf;
    end
    
    % Note: full algorithm would recursively work on 
    %       A(1:l,1:l) and A(l+1:n,l+1:n)
    
    % eigenvalue plot
%     figure, subplot(1,2,1)
%         % plot splitting line
%         t = linspace(-1.5,1.5);
%         plot(1i*t,'--k'), hold on
%         % plot original eigenvalues            
%         plot(real(E),imag(E),'b*'), hold off
%         % plot eigenvalues of transformed submatrices                  
%         %plot(real(E1),imag(E1),'go');
%         %plot(real(E2),imag(E2),'bo'), hold off
%         if strcmp(type,'Cond2')
%             format short G
%             title(['Eigenvalues, angle = ' num2str(angle,'%1.1e')])
%         else
%             title('Eigenvalues')
%         end
%         xlabel('Real')
%         ylabel('Imag')
%         axis([-1.5 1.5 -1.5 1.5])
%         axis square
    figure, subplot(1,2,1)
        plot(1:numit-1,splt(1:numit-1),'-+b',1:numit,(n-neg_count)*ones(1,numit),'k--')
        axis([1 numit 1 n])
        axis square
        xlabel('IRS Iterations')
        ylabel('Split size')
        title('Best Split')
        
        
    % convergence plot
    subplot(1,2,2),
        semilogy(1:numit,bckerr(1:numit),'--or',1:numit,Rconv(1:numit),'-+b')
        %semilogy(1:numit,bckerr(1:numit),'-or')
        title('Convergence')
        xlabel('IRS Iterations')
        ylabel('Backward Error')
        axis([1 numit 10^(-17) 1])
        axis square
        %legend('BckErr','Rconv')
    % Note: convergence plot shows backward error of divide and conquer 
    %       and NOT the criterion within the IRS iteration (convergence    
    %       of the R matrix), though that's interesting to plot as well
    %subplot(2,2,3),
    %plot(s)
    %title('Singular Values of Eigenvector Matrix')
    
    if ~strcmp(fname,'noeps')
        saveas(gcf,sprintf('%s.eps',fname),'epsc');
    end
    iter = numit-1;
    err = bckerr(numit-1);
      

    
end

function [Ap,Bp,bckerr,Rconv,splt,numit] = irs(A,B,alpha,beta,gamma,delta,maxit)
% implicit repeated squaring (with Mobius transformation)

    n = size(A,1);
    normA = norm(A);

    % initialize with Mobius transformation
    Ap = alpha*A + beta*B;
    Bp = gamma*A + delta*B;
    R = zeros(n);
    Rold = ones(n);
    bckerr = zeros(1,maxit);
    Rconv = zeros(1,maxit);
    splt = zeros(1,maxit);
    
    % perform repeated squaring
    numit = 1; converged = false;
    while norm(R-Rold,1)/norm(Rold,1) > 10*n*eps && numit <= maxit,
    %while numit <= 7 || (norm(R-Rold,1)/norm(Rold,1) > 10*n*eps && ~converged && numit <= maxit),
        
        Rold = R;
       
        % QR decomposition of size 2n x n
        [Q,R] = qr([Bp;-Ap]);

        % enforce positive diagonal entries of R
        D = diag(sign(diag(R(1:n,:))));
        R(1:n,:)=D*R(1:n,:);
        Q(:,1:n)=Q(:,1:n)*D;

        % two n x n matrix multiplies
        Ap = Q(1:n,n+1:2*n)'*Ap;
        Bp = Q(n+1:2*n,n+1:2*n)'*Bp;

        % extract square R
        R = R(1:n,:);

        % complete the divide and conquer with GRURV to compute 
        % the backward error (for the purposes of plotting)
        [U,R1,V] = rurv(Ap);
        [R2,Qr] = rq(U'*(Ap+Bp));
        A2 = Qr*A*Qr';
        [E21norm,l] = split(A2);
        bckerr(numit) = E21norm/normA;
        Rconv(numit) = norm(R-Rold,1)/norm(Rold,1);
        splt(numit) = l;

        if (numit >= 3 && bckerr(numit-2) < 10*n*eps && bckerr(numit-1) < 10*n*eps)
        	converged = true;
        end
        
        numit = numit+1;

    end
    
    if numit > maxit,
        numit = maxit;
    end
    
    
end

function [U,R,V] = rurv(A)
% compute rank-revealing decomposition of A (A=URV)
% U is orthogonal
% R is upper triangular and rank-revealing (whp)
% V is Haar matrix (random orthogonal)

    % create Haar matrix
    [V,R1] = qr(randn(size(A,1)));

    % compute rank-revealing decomposition of A
    [U,R] = qr(A*V');
    
end

function [minnorm,l] = split(A)
% find SW block of A with minimum norm

    n = size(A,1);
    val = ones(1,n-1);
    
    % compute norm of SW block for each dimension
    for j = 1:n-1, 
        val(j) = norm(A(j+1:n,1:j),2);
    end

    [minnorm,l] = min(val);
    
end

function [R,Q] = rq(A)
% compute RQ factorization using QR
% R is upper triangular
% Q is orthogonal

    [Q,R] = qr(flipud(A)');
    R = fliplr(flipud(R'));
    Q = flipud(Q');
    
end

