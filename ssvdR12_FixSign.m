function [U1 U2 V] = ssvdR12_FixSign(A,B,alpha) 
% Usage: [U1 U2 V] = ssvdR12_FixSign(A,B,alpha) 
% Input: A = P1xN data matrix
%        B = P2xN data matrix
%        alpha = l1 "sparsity" parameter.
% Output: U1 = vector of loadings of rows in A
%         U2 = vector of loadings of rows in B
%         V = vector of scores for dominant signal in the row-space
%             of matrix stack [A;B].
% Description: This version of l1 SVD fixes the sign problem inherent
%              to the SVD.  

% Get dimensions of input data matrices A and B
[m1,n1]=size(A);
[m2,n2]=size(B);

% Check that A and B have the same number of columns
if n1==n2
    n=n1;
else
    fprintf('The data should have the same number of columns \n');
    U1=[];
    U2=[];
    V=[];
    return;
end;

% Form stacked matrix
X=[A;B];

% Set initial state vector
[t1 t2 t3] = svd(X,'econ');

% initial value of left singular vector
u0 = t1(:,1);
% ud is used to to decide convergence
ud = 1;
% number of iterations
iter_cnt = 0; 

% Fix sign of initial state
mu = mean(X)';
t31 = t3(:,1);
corrmut31 = corr(mu,t31);
if corrmut31 <= 0
    u0 = -u0;
end

% Compute sparse, rank-1 approximation
while (ud > 0.0001) 
    iter_cnt = iter_cnt+1;    
    v =  X'*u0/sqrt(sum((X'*u0).^2)); %updating v    
    u = sign(X*v).*max(abs(X*v)-alpha,0);%updating u
    s =  sqrt(sum(u.^2)); %singular value
    u = u/s;% normalizing u    
    ud = sqrt(sum((u0-u).^2));%ud=||u-u0||
    u0 = u;
    if iter_cnt > 5000 %5000 is the maximum number of iterations
        disp('Failed to converge! Increase the limit on the maximum number of iterations')
        break
    end    
end

% Get dominant signal in row-space of stacked matrix
V=sqrt(sum((X'*u).^2))*v;
% Get loadings on rows of matrix A
U1=u(1:m1);
% Get loadings on rows of matrix B
U2=u(m1+1:m1+m2);

