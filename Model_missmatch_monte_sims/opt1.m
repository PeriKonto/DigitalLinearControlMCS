%Optimal Control Program
%By Peri Kontoroupis
%Effect of changing the diagonal variables in a Q matrix
%as seen on K values

A = [0.5870 32.8684 0; 0 0 0; -0.5870 -32.8684 1;];
B = [0 1 0]';
D = [0 0 1];
R = 1;


%x(k)=Ax(k-1)+Bu(k-1)+Dyd(k)
%Ks=[];
h = [1 0 0];
%for i = 1:-0.001:0.01
%for i = 1:1:100
Q = [1 0 0;0 1 0; 0 0 0.01;];
[K,S,E] = dlqr(A,B,Q,R);
%Ks=[Ks K];   
%end

F_closed = (A-B*K);
d_closed = zeros(size(h,1),size(D,2));
sys = ss(F_closed,D',h,0,-1);
rlocus(sys);




%----------------------------------------------------
%Values1=[];
%Values2=[];
%Values3=[];
%l = length(Ks);
%Extract every third value from Ks 

%for j = 3:3:l
%Value1=Ks(j-2);   
%Value2=Ks(j-1);   
%Value3=Ks(j); 
%Values1 = [Values1 Value1];
%Values2 = [Values2 Value2];
%Values3 = [Values3 Value3];
%end
%subplot(221),plot(Values1);
%subplot(221),title('Yk) variation');
%subplot(221),grid
%subplot(221),xlabel('Number of iteration');
%subplot(221),ylabel('K value');
%subplot(222),plot(Values2,'r');
%subplot(222),title('U(k-1) variation');
%subplot(222),grid
%subplot(222),xlabel('Number of iteration');
%subplot(222),ylabel('K value');
%subplot(223),plot(Values3,'c');
%subplot(223),title('Z(k) variation');
%subplot(223),grid
%subplot(223),xlabel('Number of iteration');
%subplot(223),ylabel('K value');

%fprintf('Integral of error state is %d\n', i); 
%fprintf('Matrix Weights: %0.2f\n', 0.01, 1, i);
%fprintf('K values: %0.2f\n', Ks(i), Ks(i+1),Ks(i+2));

