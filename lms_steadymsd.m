%LMS algorithm MSD analysis
function lms_steadymsd(x0,mu,N,t,G,B,itr)
    M=[10,15,20,30,40,50];
    Uf=G.U(:,1:t);
    steady_msd=zeros(3,6);
    for i=1:200
        s_msd=zeros(3,6);
        
        %noise
        mean=zeros(1,N);
        Cv=diag(0.01 .* rand(1,N));
        v=mvnrnd(mean,Cv,1)';

        %signal initialization
        s=zeros(N,1);
        s(1:t)= -1 + 2 .*rand(t,1);
        f1 = gsp_igft(G,s);
        f2 = gsp_igft(G,s);
        f3 = gsp_igft(G,s);
        
        
        for l=1:6
            %D matrix
            [~,D1]=maxdet(t,M(l),N,G);
            [~,D2]=max_mineig(t,M(l),N,G);
            r=randperm(N,M(l));
            D3=zeros(N,N);
            for z=1:M(l)
                D3(r(z),r(z))=1;
            end 
            
            
            %phi matrix
           ph1=(eye(t)-(mu .* Uf'*D1*Uf))*(eye(t)-(mu .* Uf'*D1*Uf));
           ph2=(eye(t)-(mu .* Uf'*D2*Uf))*(eye(t)-(mu .* Uf'*D2*Uf));
           ph3=(eye(t)-(mu .* Uf'*D3*Uf))*(eye(t)-(mu .* Uf'*D3*Uf));



            %Lms algorithm
            j=1;
            while j<=itr
                y1= D1* B * x0 + D1 * v;
                f1=f1+ mu * B * D1 * (y1-f1);
                y2= D2* B * x0 + D2 * v;
                f2=f2+ mu * B * D2 * (y2-f2);
                y3= D3* B * x0 + D3 * v;
                f3=f3+ mu * B * D3* (y3-f3);
                j=j+1;
            end
            s1=Uf'*(f1-x0);
            s2=Uf'*(f2-x0);
            s3=Uf'*(f3-x0);
            
            s_msd(1,l)= s1'* ph1 * s1;
            s_msd(2,l)=s2' * ph2 * s2;
            s_msd(3,l)=s3' * ph3 * s3;
         end
        steady_msd=steady_msd+s_msd;
    end
    steady_msd=steady_msd/200;
    
    %plot
    figure(8);
    plot(M,10*log10(steady_msd(1,:)),'-o','LineWidth',2,'MarkerSize',10);
    hold on;
    plot(M,10*log10(steady_msd(2,:)),'-o','LineWidth',2,'MarkerSize',10);
    hold on;
    plot(M,10*log10(steady_msd(3,:)),'-o','LineWidth',2,'MarkerSize',10);
    title('Steady state MSD');
    xlabel('Number of samples');
    ylabel('Steady state-MSD MSD (db)');
    legend('Max-det','Max-Mineig','Random Sampling strategy');
    grid on;
end