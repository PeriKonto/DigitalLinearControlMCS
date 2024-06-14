%second case: keeping the model parameters constant and changing the gains
clear all
clc
startup
sriv
cd c:\Matlab~1\work\control\monte2
[BB,A,B]=mcpar(TH,1000,1);
Alength=length(A);
Blength=length(B);
if (Alength < Blength) & (Alength > Blength)
   error('Invalid Size!!!')
   return
end

A1=A(:,1); B1=B(:,1);
A2=A(:,2); B2=B(:,2); 
fs =[]; gs =[]; ks =[]; Ys=[]; Us=[];
figure(2)
   
for i = 1:1000
   at=ad;
   bt=bd;
   [f,g,k]=con_gains(at(2),bt(2),at,bt);
   fs=[fs f];
   gs=[gs g];
   ks=[ks k];
   sim('pipcontrol');
   Ys=[Ys y];
 end
 
 subplot(331),plot(fs);
 subplot(331),title('  Proportional Gain Variation f ');
 subplot(332),plot(gs);
 subplot(332),title('Feed Forward Gain G ');
 subplot(333),plot(ks);
 subplot(333),title('Integral Gain Variation Ki ');
 
 subplot(312),plot(Ys);
 subplot(312),title('Monte Carlo Simulations x1000:  Systems response');
 subplot(313),plot(Us);
 subplot(313),title('Monte Carlo Simulations x1000:  Systems Input');
