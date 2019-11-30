clc;
clear all;
%Creating graph
N=50;  %number of vertices
M=10;  %bandwidth
S=10;  %number of samples
G=gsp_sensor(N);   %graph
figure(1);
gsp_plot_graph(G);
title('Sensor Graph');


%Initialisation
mu=0.5;     %learning rate
Cv=diag(0.01 .* rand(1,N));   %covariance matrix


%Laplace Transform
G = gsp_compute_fourier_basis(G);    %laplace transform


%Bandlimited signal
s=zeros(N,1);
s(1:M)= -2 + 4.*rand(M,1);
f = gsp_igft(G,s);    
figure(2);
gsp_plot_signal(G,f);
title('Sensor graph with bandlimited signal');


%B matrix
sigma=zeros(N,N);
for i=1:10
    sigma(i,i)=1;
end
B= G.U * sigma * (G.U)';


%max_det sampling
[s_md,~]=maxdet(M,S,N,G);
disp('Sampled set of vertices for Max-Det Algorithm');
disp(s_md);
param.vertex_highlight=s_md;
figure(3);
gsp_plot_signal(G,f,param);
title('Sensor graph with Sampled vertices (Max-Det)');


%max_mineig sampling
[s_me,~]=max_mineig(M,S,N,G);
disp('Sampled set of vertices for Max-mineig Algorithm');
disp(s_me);
param.vertex_highlight=s_me;
figure(4);
gsp_plot_signal(G,f,param);
title('Sensor graph with Sampled vertices (Max-mineig)');


%minmsd sampling
[s_mm,~]=minmsd(M,S,N,G,mu,Cv);
disp('Sampled set of vertices for Min-MSD Algorithm');
disp(s_mm);
param.vertex_highlight=s_mm;
figure(5);
gsp_plot_signal(G,f,param);
title('Sensor graph with Sampled vertices (Min-MSD)');



%LMS algorithm
itr=100;
lms(f,mu,M,N,S,G,Cv,B,itr);


%LMS algorithm transient MSD error
lms_msd(f,mu,N,M,G,B,itr);


%LMS algorithm steady state error
lms_steadymsd(f,mu,N,M,G,B,itr);


%LMS with unknown bandwidth
lambda=0.1;
s0=zeros(N,1);
s0(1:M)= ones(M,1);
lms_withbdwidth(N,M,G,mu,s,lambda);