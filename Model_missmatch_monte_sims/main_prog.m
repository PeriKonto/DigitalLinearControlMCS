%first case: keeping the model gains constant and changing the parameters 

clear all
clc
startup
sriv
cd c:\Matlab~1\work\control\monte

[BB,A,B]=mcpar(TH,1000,1); %monte carlo simulation function 
Alength=length(A);
Blength=length(B);
if (Alength < Blength) & (Alength > Blength)
   error('Invalid Size!!!')
   return
end


A1=A(:,1); B1=B(:,1);
A2=A(:,2); B2=B(:,2); 
fs =[]; gs =[]; ks =[]; Ys=[]; Us=[];

%introducing 10% uncertainty to the 2nd parameter of the numerator
A2=(A2.*0.1)+A2;

at=[A1(1) A2(1)];
bt=[B1(1) B2(1)];
[f,g,k]=con_gains(A2(1), B2(1),at,bt);
figure(3)
sim('pipcontrol');
subplot(211),plot(y);
subplot(211),title('Systems response '); 
subplot(211),xlabel('Time intervals ');
subplot(212),plot(u,'r');
subplot(212),title('Systems input' );
subplot(212),xlabel('Time intervals ');


for i = 2:1000
        at=[A1(i) A2(i)];
        bt=[B1(i) B2(i)];
        sim('pipcontrol');
        Ys=[Ys y];
        Us=[Us u];
end
figure(4)
subplot(211),plot(Ys);
subplot(211),title('Monte Carlo Simulations x1000:  Systems response');
subplot(212),plot(Us);
subplot(212),title('Monte Carlo Simulations x1000:  Systems Input');

