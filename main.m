clc;
clear all;
%Creating graph
N=50;
M=10;
G=gsp_sensor(N);
figure(1);
gsp_plot_graph(G);
title('Sensor Graph');

%Initialisation
mu=0.5;
Cv=0.01.*rand(N,N);

%Laplace Transform
G = gsp_compute_fourier_basis(G);

%Bandlimited signal
b=10;
s=zeros(N,1);
s(1:10)= -2 + 4.*rand(b,1);
f = gsp_igft(G,s);
figure(2);
gsp_plot_signal(G,f);
title('Sensor graph with bandlimited signal');

%max_mineig sampling
s_md=maxdet(M,N,G);
disp('Sampled set of vertices for Max-Det Algorithm');
disp(s_md);
param.vertex_highlight=s_md;
figure(3);
gsp_plot_signal(G,f,param);
title('Sensor graph with Sampled vertices (Max-Det)');

%maxdet sampling
s_me=max_mineig(M,N,G);
disp('Sampled set of vertices for Max-mineig Algorithm');
disp(s_me);
param.vertex_highlight=s_me;
figure(4);
gsp_plot_signal(G,f,param);
title('Sensor graph with Sampled vertices (Max-mineig)');

%minmsd sampling
s_mm=minmsd(M,N,G,mu,Cv);
disp('Sampled set of vertices for Min-MSD Algorithm');
disp(s_mm);
param.vertex_highlight=s_mm;
figure(5);
gsp_plot_signal(G,f,param);
title('Sensor graph with Sampled vertices (Min-MSD)');