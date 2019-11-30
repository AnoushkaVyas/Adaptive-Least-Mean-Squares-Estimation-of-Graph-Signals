%LMS algorithm for graph signals
function lms(x0,mu,M,N,S,G,Cv,B,itr)
    %noise
    mean=zeros(1,N);
    v=mvnrnd(mean,Cv,1)';
    
    %signal initialization
    s=zeros(N,1);
    s(1:M)= -1 + 2 .*rand(M,1);
    f = gsp_igft(G,s);
    
    %D matrix
    [~,D]=maxdet(M,S,N,G);
   
    %Lms algorithm
    j=1;
    while j<=itr
        y= D* B * x0 + D * v;
        f=f+ mu * B * D * (y-f);
        j=j+1;
    end
    
    %plot
    figure(6);
    plot(f,'-o','LineWidth',2,'MarkerSize',10);
    hold on;
    plot(x0,'-o','LineWidth',2,'MarkerSize',10);
    title('Signal Reconstruction using Max-Det strategy with 10 samples');
    xlabel('Node Index');
    ylabel('Graph Signal');
    legend('Estimated Signal','Original Signal');
    grid on;
end