%LMS algorithm transient MSD analysis
function lms_msd(x0,mu,N,M,G,B,itr)
    Uf=G.U(:,1:M);
    [~,D1]=maxdet(M,10,N,G);
    [~,D2]=maxdet(M,14,N,G);
    [~,D3]=maxdet(M,18,N,G);
    final_msd=zeros(3,itr);
    for i=1:200
        %error initialization
        msd=zeros(3,itr);
        
        %noise
        mean=zeros(1,N);
        Cv=diag(0.01 .* rand(1,N));
        v=mvnrnd(mean,Cv,1)';

        %signal initialization
        s=zeros(N,1);
        s(1:M)= -1 + 2 .*rand(M,1);
        f = gsp_igft(G,s);
        s0= Uf'*(f-x0);
        
        for l=1:3
            if l==1
                D=D1;
            end
             if l==2
                D=D2;
             end
             if l==3
                D=D3;
             end
            TT=(eye(M)-(mu .* Uf'*D*Uf));
            ph=TT * TT;
            msd(l,1)=s0' * ph * s0;

            %Lms algorithm
            j=2;
            while j<=itr
                y= D* B * x0 + D * v;
                f=f+ mu * B * D * (y-f);
                s0= Uf'* (f-x0);
                msd(l,j)= s0' * ph * s0;
                j=j+1;
           end
         end
        final_msd=final_msd+msd;
    end
    final_msd=final_msd/200;
    
    %plot
    figure(7);
    plot(10*log10(final_msd(1,:)),'-.','LineWidth',2,'MarkerSize',10);
    hold on;
    plot(10*log10(final_msd(2,:)),'-','LineWidth',2,'MarkerSize',10);
    hold on;
    plot(10*log10(final_msd(3,:)),'--','LineWidth',2,'MarkerSize',10);
    title('Transient MSD using Max-Det strategy');
    xlabel('Iterations');
    ylabel('Transient MSD (db)');
    legend('|S|=10','|S|=14', '|S|=18');
    grid on;
end