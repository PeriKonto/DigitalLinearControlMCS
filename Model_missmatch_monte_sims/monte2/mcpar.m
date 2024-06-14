function [BB,A,B]=mcpar(TH,nt,Pfac);
% MCPAR - Monte Carlo realisation of a transfer function model
%
% function [BB,A,B]=mc_par(th,nt); 
% 
% th:    theta matrix as obtained from e.g. riv
% nt:    required number of realisations
% Pfac:  if required, a factor to multiply the P matrix
%        - normally 1
%
% BB:    nt row vectors of random parameters normally distributed
%        according to the parametric information in theta
% AA,BB: for SISO systems contain nt rows of denominator
%        and numerator coefficients respectively
%
%  See the example run in mcparex.m

% W. Tych 1993/2000

if nargin==1, nt=50; Pfac=1; end
if nargin==2, Pfac=1; end

[m1,m2]=size(TH);
if m1<=3
 error('MCPAR: No covariance information given in TH');
end
m=m1-3;
randn('seed',213184);
P=TH(4:m1,1:m);  
if any(eig(P)<=0), 
   disp('Warning  MCPAR: P matrix is not positive-definite');
   BB=[]; A=[]; B=[]; return,
end

P=chol(P)*Pfac;  
par=TH(3,1:m);

%BB=zeros(nt+1,m);
%for k=1:nt
%   BB(k+1,:)=par+randn(1,m)*P;
%end

BB=randn(nt+1,m)*P+par(ones(nt+1,1),:);  % large gain on time

BB(1,:)=par;

na=TH(1,4);
nu=TH(1,3);
if nu>1; return;end

if nu>0
   nb=TH(1,5:4+nu);
else
   nb=0;
end
nc=TH(1,5+nu);
nk=TH(1,7+2*nu);

A=[ones(nt+1,1) BB(:,1:na)];
B=[zeros(nt+1,nk) BB(:,(na+1):(na+nb))];
