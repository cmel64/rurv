function singularNTest()
0
%singularNTest2(1000,400,10000,[20,40,80,160,320],10000,150,4,1);
1
%singularNTest2(1000,400,1,[1000,100000,10000000,1000000000,100000000000,10000000000000],10000,400,4,2);
2
%singularNTest2([100,200,400,600,800,1000,1200],...
%    [50,100,200,300,400,500,600],10000,40,10000,150,3,1);
3
singularNTest2([100,200,400,600,800,1000,1200],...
     [50,100,200,300,400,500,600],1,1000,10000,100,3,2);
4
%singularNTest2([800,800,800,800,800,800,800,800],...
%     [50,150,250,350,450,550,650,750],10000,40,10000,150,2,1);
5
%singularNTest2([800,800,800,800,800,800,800,800],...
%     [50,150,250,350,450,550,650,750],1,40,10000,150,2,2);

 
end


function [riiSigMax, riiSigMin, c1, c2, c3] = singularNTest2(n,k,initSig,C,endSig,p,mode,eigt)
% produces test statistics for matrix with singular values in initSig, kSig
% and endSig for the rurv algorithm

n
kSig = C * n

if (max(size(C)) == 1)
    C = C *ones(size(kSig));
end

delta = 3.3/C(1)*ones(size(kSig));

riiSigMax = zeros(max(size(initSig)), max(size(kSig)), max(size(endSig)), p);
riiSigMin = zeros(max(size(initSig)), max(size(kSig)), max(size(endSig)), p);
%c1 = zeros(max(size(initSig)), max(size(kSig)), max(size(endSig)), p,max(size(n)));
%c2 = zeros(max(size(initSig)), max(size(kSig)), max(size(endSig)), p,max(size(n)));
%c3 = zeros(max(size(initSig)), max(size(kSig)), max(size(endSig)), p,max(size(n)));
minR = zeros(max(size(initSig)), max(size(kSig)), max(size(endSig)), p,max(size(n)));
maxR = zeros(max(size(initSig)), max(size(kSig)), max(size(endSig)), p,max(size(n)));
RinvR = zeros(max(size(initSig)), max(size(kSig)), max(size(endSig)), p,max(size(n)));

for o = 1:p
for j = 1:max(size(initSig))
    for l = 1:max(size(kSig))
        if (C(l)*delta < 1.01)
            error(strcat('ERROR C * delta = ', num2str(C(l)*delta), ' < 1.01 for kSig =', num2str(kSig(l))));
        end
        for m = 1:max(size(endSig))
            for q = 1:max(size(n))
            [tMinR, tMaxR, tRinvR] = singTestStat(n(q),k(q),initSig(j),kSig(l),endSig(m),eigt);
            minR(j,l,m,o,q) = tMinR;
            maxR(j,l,m,o,q) = tMaxR;
            RinvR(j,l,m,o,q) = tRinvR;
%            c1(j,l,m,o,q) = temp2;
%            c2(j,l,m,o,q) = temp3;
%            c3(j,l,m,o,q) = temp4;
            
            end
        end
    end    
end
o
end

avgMinR = mean(minR,4);
avgMaxR = mean(maxR,4);
avgRinvR = mean(RinvR,4);

q1MinR = quantile(minR,0.25,4);
q1MaxR = quantile(maxR,0.25,4);
q1RinvR = quantile(RinvR,0.25,4);

q2MinR = quantile(minR,0.5,4);
q2MaxR = quantile(maxR,0.5,4);
q2RinvR = quantile(RinvR,0.5,4);

q3MinR = quantile(minR,0.75,4);
q3MaxR = quantile(maxR,0.75,4);
q3RinvR = quantile(RinvR,0.75,4);

h = figure;

midn = floor(max(size(n))/2)+1;
midk = floor(max(size(k))/2)+1;
midks = floor(max(size(kSig))/2)+1;
diagA = ones(n(midn),1);
if(eigt == 1)
    dEndSig = endSig^(1/(n(midn)-k(midk)-1));
    dInitSig = initSig^(1/(k(midk)-1));
    diagA(k(midk)+1:n(midn)) = dEndSig.^(n(midn)-k(midk)-1:-1:0);
    diagA(k(midk)) = diagA(k(midk)+1) * kSig(midks);
    diagA(1:k(midk)-1) = diagA(k(midk)) * dInitSig.^(k(midk)-1:-1:1);
else
    diagA(k(midk)+1:n) = 1;
    diagA(k(midk)) = kSig(midks);
    diagA(1:k(midk)-1) = kSig(midks);
end
    
for j = 1:max(size(initSig))
    for l = 1:max(size(kSig))
        for m = 1:max(size(endSig))
            c = initSig(j)
            dleta = delta(l)
            d1 = ((2.02*initSig(j)/delta(l))^3)
            d2 = (1.01/(delta(l)*C(l)))
            d3 = (1-(1.01/(delta(l)*C(l)))^2)
            const1 = (((2.02*initSig(j)/delta(l))^3)*(1.01/(delta(l)*C(l)))/(1-(1.01/(delta(l)*C(l)))^2))
            const2 = ((2.02*initSig(j))^2)/(1-(1.01/(delta(l)*C(l)))^2)
            if (mode == 1)
            subplot(2,2,1), semilogy(1:n(midn), diagA);
            title('Singular value distribution of test matrices');
            ylabel('\sigma');
            xlabel('singular values in descending order');
                
            bound2 = (2.02/delta(l)) * sqrt(k.*(n-k));
            n_tot = reshape(transpose(repmat(reshape(n,[max(size(n)),1]),[1,p])),1,[]);
            minR_tot = reshape(minR(j,l,m,:,:), 1, []);
            subplot(2,2,2), semilogy(n_tot, minR_tot, '.r',...
                n, reshape(q1MinR(j,l,m,:),[1,max(size(n))]), '-o',...
                n, reshape(q2MinR(j,l,m,:),[1,max(size(n))]), '-og',...
                n, reshape(q3MinR(j,l,m,:),[1,max(size(n))]), '-om',...
                n, bound2, '-oc');
            title(sprintf('Testing Inequality 1 with r = %d and delta = %d',k(1),delta(l)));
            ylabel('\sigma_r / \sigma_{min} (R_{11})');
            xlabel('n');
            
            bound3 = const1*(k.*(n-k)).^(3/2) + const2*(k.*(n-k)) + 1;
            maxR_tot = reshape(maxR(j,l,m,:,:), 1, []);
            %print(h,'-djpeg',sprintf('nmin_%d_%d_%d.jpg',initSig(j),kSig(l),endSig(m)));
            subplot(2,2,3), semilogy(n_tot, maxR_tot, '.r',...
                n, reshape(q1MaxR(j,l,m,:),[1,max(size(n))]), '-o',...
                n, reshape(q2MaxR(j,l,m,:),[1,max(size(n))]), '-og',...
                n, reshape(q3MaxR(j,l,m,:),[1,max(size(n))]), '-om',...
                n, bound3, '-oc');
            title(sprintf('Testing Inequality 2 with r = %d and delta = ',k(1),delta(l)));
            ylabel('\sigma_{max} (R_{22}) / \sigma_{r+1}');
            xlabel('n');
            %print(h,'-djpeg',sprintf('nmax_%d_%d_%d.jpg',initSig(j),kSig(l),endSig(m)));
            
            bound4 = (1/(1-(1.01/(delta(l)*C(l)))^2))*((2.02/delta(l))*sqrt(k.*(n-k)) + (1.01^2)/((delta(l)*C(l))^2));
            RinvR_tot = reshape(RinvR(j,l,m,:,:), 1, []);
            subplot(2,2,4), semilogy(n_tot, RinvR_tot, '.r',...
                n, reshape(q1RinvR(j,l,m,:),[1,max(size(n))]), '-o',...
                n, reshape(q2RinvR(j,l,m,:),[1,max(size(n))]), '-og',...
                n, reshape(q3RinvR(j,l,m,:),[1,max(size(n))]), '-om',...
                n, bound4, '-oc');
            title(sprintf('Testing Inequality 3 with r = %d and delta = %d',k(1), delta(l)));
            ylabel('|| R_{11}^{-1} R_{12} ||_2');
            xlabel('n');
            
            
            print(h,'-depsc',sprintf('n_%d_%d_%d_%d.eps',initSig(j),kSig(l),endSig(m),eigt));

            elseif (mode == 2)
            subplot(2,2,1), semilogy(1:n(midn), diagA);
            title('Singular value distribution of test matrices');
            ylabel('\sigma');
            xlabel('singular values in descending order');
            
            bound2 = (2.02/delta(l)) * sqrt(k.*(n-k));
            k_tot = reshape(transpose(repmat(reshape(k,[max(size(k)),1]),[1,p])),1,[]);
            minR_tot = reshape(minR(j,l,m,:,:), 1, []);
            subplot(2,2,2), semilogy(k_tot, minR_tot, '.r',...
                k, reshape(q1MinR(j,l,m,:),[1,max(size(n))]), '-o',...
                k, reshape(q2MinR(j,l,m,:),[1,max(size(n))]), '-og',...
                k, reshape(q3MinR(j,l,m,:),[1,max(size(n))]), '-om',...
                k, bound2, '-oc');
            title(sprintf('Testing Inequality 1 with n = %d and delta = %d',n(1),delta(l)));
            ylabel('\sigma_r / \sigma_{min} (R_{11})');
            xlabel('r');
            
            bound3 = const1*(k.*(n-k)).^(3/2) + const2*(k.*(n-k)) + 1;
            maxR_tot = reshape(maxR(j,l,m,:,:), 1, []);
            %print(h,'-djpeg',sprintf('kmin_%d_%d_%d.jpg',initSig(j),kSig(l),endSig(m)));
            subplot(2,2,3), semilogy(k_tot, maxR_tot, '.r',...
                k, reshape(q1MaxR(j,l,m,:),[1,max(size(n))]), '-o',...
                k, reshape(q2MaxR(j,l,m,:),[1,max(size(n))]), '-og',...
                k, reshape(q3MaxR(j,l,m,:),[1,max(size(n))]), '-om',...
                k, bound3 , '-oc');
            title(sprintf('Testing Inequality 2 with n = %d and delta = %d',n(1),delta(l)));
            ylabel('\sigma_{max} (R_{22}) / \sigma_{r+1}');
            xlabel('r');
            %print(h,'-djpeg',sprintf('kmax_%d_%d_%d.jpg',initSig(j),kSig(l),endSig(m)));
            
            bound4 = (1/(1-(1.01/(delta(l)*C(l)))^2))*((2.02/delta(l))*sqrt(k.*(n-k)) + (1.01^2)/((delta(l)*C(l))^2));
            RinvR_tot = reshape(RinvR(j,l,m,:,:), 1, []);
            subplot(2,2,4), semilogy(k_tot, RinvR_tot, '.r',...
                k, reshape(q1RinvR(j,l,m,:),[1,max(size(n))]), '-o',...
                k, reshape(q2RinvR(j,l,m,:),[1,max(size(n))]), '-og',...
                k, reshape(q3RinvR(j,l,m,:),[1,max(size(n))]), '-om',...
                k, bound4, '-oc');
            title(sprintf('Testing Inequality 3 with n = %d and delta = %d',n(1),delta(l)));
            ylabel('|| R_{11}^{-1} R_{12} ||_2');
            xlabel('r');
            
            print(h,'-depsc',sprintf('k_%d_%d_%d_%d.eps',initSig(j),kSig(l),endSig(m),eigt));
            
            elseif (mode == 3)    
            subplot(2,2,1), semilogy(1:n(midn), diagA);
            title('Singular value distribution of test matrices');
            ylabel('\sigma');
            xlabel('singular values in descending order');
                
            bound2 = (2.02/delta(l)) * sqrt(k.*(n-k));
            n_tot = reshape(transpose(repmat(reshape(n,[max(size(n)),1]),[1,p])),1,[]);
            minR_tot = reshape(minR(j,l,m,:,:), 1, []);
            subplot(2,2,2), plot(n_tot, minR_tot, '.r',...
                n, reshape(q1MinR(j,l,m,:),[1,max(size(n))]), '-o',...
                n, reshape(q2MinR(j,l,m,:),[1,max(size(n))]), '-og',...
                n, reshape(q3MinR(j,l,m,:),[1,max(size(n))]), '-om',...
                n, bound2, '-oc');
            title(sprintf('Testing Inequality 1 with r = %f * n and delta = %d',k(1)/n(1),delta(l)));
            ylabel('\sigma_r / \sigma_{min} (R_{11})');
            xlabel('n');
            
            bound3 = const1*(k.*(n-k)).^(1/2) + const2*(k.*(n-k)).^(1/2) + 1;
            maxR_tot = reshape(maxR(j,l,m,:,:), 1, []);
            %print(h,'-djpeg',sprintf('nkmin_%d_%d_%d.jpg',initSig(j),kSig(l),endSig(m)));
            subplot(2,2,3), plot(n_tot, maxR_tot, '.r',...
                n, reshape(q1MaxR(j,l,m,:),[1,max(size(n))]), '-o',...
                n, reshape(q2MaxR(j,l,m,:),[1,max(size(n))]), '-og',...
                n, reshape(q3MaxR(j,l,m,:),[1,max(size(n))]), '-om',...
                n, bound3, '-oc');
            title(sprintf('Testing Inequality 2 with r = %f * n and delta = %d',k(1)/n(1),delta(l)));
            ylabel('\sigma_{max} (R_{22}) / \sigma_{r+1}');
            xlabel('n');
            %print(h,'-djpeg',sprintf('nkmax_%d_%d_%d.jpg',initSig(j),kSig(l),endSig(m)));
            
            bound4 = (1/(1-(1.01/(delta(l)*C(l)))^2))*((2.02/delta(l))*sqrt(k.*(n-k)) + (1.01^2)/((delta(l)*C(l))^2));
            RinvR_tot = reshape(RinvR(j,l,m,:,:), 1, []);
            subplot(2,2,4), plot(n_tot, RinvR_tot, '.r',...
                n, reshape(q1RinvR(j,l,m,:),[1,max(size(n))]), '-o',...
                n, reshape(q2RinvR(j,l,m,:),[1,max(size(n))]), '-og',...
                n, reshape(q3RinvR(j,l,m,:),[1,max(size(n))]), '-om',...
                n, bound4, '-oc');
            title(sprintf('Testing Inequality 3 with r = %f * n and delta = %d',k(1)/n(1),delta(l)));
            ylabel('|| R_{11}^{-1} R_{12} ||_2');
            xlabel('n');
            print(h,'-depsc',sprintf('nk2_%d_%d_%d_%d.eps',initSig(j),kSig(l),endSig(m),eigt));
            end
        end
    end
end

for j = 1:max(size(initSig))
    for l = 1:max(size(n))
        for m = 1:max(size(endSig))
            const1 = ((2.02*initSig(j)./delta).^3).*(1.01./(delta.*C))./(1-(1.01./(delta.*C)).^2);
            const2 = ((2.02*initSig(j))^2)./(1-(1.01./(delta.*C)).^2);
            if (mode == 4)
            subplot(2,2,1), semilogy(1:n(midn), diagA);
            title('Singular value distribution of test matrices');
            ylabel('\sigma');
            xlabel('singular values in descending order');
            
            bound2 = (2.02./delta) .* sqrt(k(1).*(n(1)-k(1)))
            kSig_tot = reshape(transpose(repmat(reshape(kSig,[max(size(kSig)),1]),[1,p])),1,[]);
            minR_tot = reshape(reshape(minR(j,:,m,:,l),[size(minR,2),size(minR,4)])', 1, []);
            subplot(2,2,2), loglog(kSig_tot, minR_tot, '.r',...
                kSig, reshape(q1MinR(j,:,m,l),[1,max(size(kSig))]), '-o',...
                kSig, reshape(q2MinR(j,:,m,l),[1,max(size(kSig))]), '-og',...
                kSig, reshape(q3MinR(j,:,m,l),[1,max(size(kSig))]), '-om',...
                kSig, bound2, '-oc');
            title(sprintf('Testing Inequality 1 with n = %d, r = %d, delta = %d',n(1),k(1),delta(1)));
            ylabel(' \sigma_r / \sigma_{min} (R_{11})');
            xlabel('\sigma_r / \sigma_{r+1}');
            
            bound3 = const1.*(k(1)*(n(1)-k(1)))^(3/2) + const2.*(k(1)*(n(1)-k(1))) + 1;
            maxR_tot = reshape(reshape(maxR(j,:,m,:,l),[size(maxR,2),size(maxR,4)])', 1, []);
            %print(h,'-djpeg',sprintf('kSmin_n%d_%d_%d.jpg',n(1),initSig(j),endSig(m)));
            subplot(2,2,3), loglog(kSig_tot, maxR_tot, '.r',...
                kSig, reshape(q1MaxR(j,:,m,l),[1,max(size(kSig))]), '-o',...
                kSig, reshape(q2MaxR(j,:,m,l),[1,max(size(kSig))]), '-og',...
                kSig, reshape(q3MaxR(j,:,m,l),[1,max(size(kSig))]), '-om',...
                kSig, bound3, '-oc');
            title(sprintf('Testing Inequality 2 with n = %d, r = %d, delta = %d',n(1),k(1),delta(1)));
            ylabel('\sigma_{max} (R_{22}) / \sigma_{r+1}');
            xlabel('\sigma_r / \sigma_{r+1}');
            %print(h,'-djpeg',sprintf('kSmax_n%d_%d_%d.jpg',n(1),initSig(j),endSig(m)));
            
            bound4 = (1./(1-(1.01./(delta.*C)).^2)).*((2.02./delta)*sqrt(k(1)*(n(1)-k(1))) + (1.01^2)./((delta.*C).^2));
            RinvR_tot = reshape(reshape(RinvR(j,:,m,:,l),[size(RinvR,2),size(RinvR,4)])', 1, []);
            subplot(2,2,4), loglog(kSig_tot, RinvR_tot, '.r',...
                kSig, reshape(q1RinvR(j,:,m,l),[1,max(size(kSig))]), '-o',...
                kSig, reshape(q2RinvR(j,:,m,l),[1,max(size(kSig))]), '-og',...
                kSig, reshape(q3RinvR(j,:,m,l),[1,max(size(kSig))]), '-om',...
                kSig, bound4, '-oc');
            title(sprintf('Testing Inequality 3 with n = %d, r = %d, delta = %d',n(1),k(1),delta(1)));
            ylabel('|| R_{11}^{-1} R_{12} ||_2');
            xlabel('\sigma_r / \sigma_{r+1}');
            
            print(h,'-depsc',sprintf('kS2_%d_%d_%d_%d.eps',n(1),initSig(j),endSig(m),eigt));
            end
        end
    end
end

end

function [minR, maxR, RinvR] = singTestStat(n,k,initSig, kSig, endSig,eigt)
        %A = eye(n,n);
        
        A = zeros(n,n);
        
        %for i = (n-1):-1:k+1
        %    A(i,i) = A(i+1,i+1) * endSig^(1/(n-k-1));
        %end
        
        %for i = (k-1):-1:1
        %    A(i,i) = A(i+1,i+1) * initSig^(1/(k-1));
        %end
        
        diagInd = linspace(1,numel(A),length(A));
        
        if (eigt == 1)
            dEndSig = endSig^(1/(n-k-1));
            dInitSig = initSig^(1/(k-1));
            
            A(diagInd(k+1:n)) = dEndSig.^(n-k-1:-1:0);
            A(k,k) = A(k+1,k+1) * kSig;
            A(diagInd(1:k-1)) = A(k,k) * dInitSig.^(k-1:-1:1);
        else
            A(diagInd(k+1:n)) = 1;
            A(k,k) = kSig;
            A(diagInd(1:k-1)) = kSig;
        end
         
        B = randn(n,n);
        [Q,~] = qr(B);
        A = Q*A;
        [~,R,~] = rurv(A);
        Sigma = sort(svd(A),'descend');
        
        %riiSig = sort(abs(diag(R)),'descend');
        %riiSig = riiSig ./ Sigma;
        
        R11 = R(1:k,1:k);
        R12 = R(1:k,(k+1):n);
        R22 = R((k+1):n, (k+1):n);
        R11min = min(svd(R11));
        R22max = max(svd(R22));
        minR = Sigma(k)/R11min;
        maxR = R22max/Sigma(k+1);
        RinvR = norm(R11\R12);
%        c1 = (sqrt(k*(n-k)) * R11min) / Sigma(k);
%        c2 = 4*(initSig^3)*(k^2)*((n-k)^2)*Sigma(k+1)/R22max;
%        c3 = norm((R11^(-1))*R12)/sqrt(k*(n-k));

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

