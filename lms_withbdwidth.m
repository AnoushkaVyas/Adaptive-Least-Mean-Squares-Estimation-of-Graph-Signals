%LMS without Bandwidth knowledge
function lms_withbdwidth(N,M,G,mu,s0,lambda)
    itr=100;
    final_card=zeros(4,itr);
    final_nmsd=zeros(3,itr);
    for p=1:10
        card_f=zeros(4,itr);
        card_f(1,:)= M * ones(1,itr);
        
        %noise
        mean=zeros(1,N);
        Cv=diag(0.0004*ones(1,N));
        v=mvnrnd(mean,Cv,1)';
        nmsd=zeros(3,itr);
        
        %Lasso
        s=rand(N,1);
        D=eye(N);
        F=[1:N];
        j=1;
        gamma= mu*lambda;
        while j<=itr
            T=zeros(N,1);
            y=D * G.U * s0 + D * v;
            sm=s+ mu* (G.U)' * D * (y-G.U*s);
            for i=1:N
                if sm(i) > gamma
                    T(i)= sm(i)-gamma;
                end
                if sm(i) < (-1*gamma)
                    T(i)=sm(i)+gamma;
                end
                if (-1*gamma) <= sm(i) <= gamma
                    T(i)=0;
                end
            end
            s=T;
            nmsd(1,j)=(norm(s-s0)/norm(s0))^2;
            for i=1:N
                if s(i) == 0
                    F(i)=0;
                end
            end
            card_f(2,j)=nnz(F);
            D=maxdetext(nnz(F),N,F,G);
            j=j+1;
        end


        %Garotte
        s=rand(N,1);
        D=eye(N);
        F=[1:N];
        j=1;
        gamma= mu*lambda;
        while j<=itr
            T=zeros(N,1);
            y=D * (G.U) * s0 + D * v;
            sm=s+ mu* (G.U)' * D * (y-(G.U)*s);
            for i=1:N
                if abs(sm(i)) > gamma
                    T(i)= sm(i)*(1-(gamma^2/sm(i)^2));
                end
                if abs(sm(i)) <= gamma
                    T(i)=0;
                end
            end
            s=T;
            nmsd(2,j)=(norm(s-s0)/norm(s0))^2;
            for i=1:N
                if s(i) == 0
                    F(i)=0;
                end
            end
            card_f(3,j)=nnz(F);
            D=maxdetext(nnz(F),N,F,G);
            j=j+1;
        end

        %Hard Thresholding
        s=rand(N,1);
        D=eye(N);
        F=[1:N];
        j=1;
        gamma= mu*lambda;
        while j<=itr
            T=zeros(N,1);
            y=D * (G.U) * s0 + D * v;
            sm=s+ mu* (G.U)' * D * (y-(G.U)*s);
            for i=1:N
                if abs(sm(i)) > gamma
                    T(i)= sm(i);
                end
                if abs(sm(i)) <= gamma
                    T(i)=0;
                end
            end
            s=T;
            nmsd(3,j)=(norm(s-s0)/norm(s0))^2;
            for i=1:N
                if s(i) == 0
                    F(i)=0;
                end
            end
            card_f(4,j)=nnz(F);
            D=maxdetext(nnz(F),N,F,G);
            j=j+1;
        end
        final_card=final_card+card_f;
        final_nmsd=final_nmsd+nmsd;
    end
    
    %plot
    figure(9);
    plot(final_card(1,:)/10,'-','LineWidth',2,'MarkerSize',10);
    hold on;
    plot(final_card(2,:)/10,'--','LineWidth',2,'MarkerSize',10);
    hold on;
    plot(final_card(3,:)/10,'-.','LineWidth',2,'MarkerSize',10);
    hold on;
    plot(final_card(4,:)/10,'.','LineWidth',2,'MarkerSize',10);
    xlabel('Iterations');
    ylabel('Cardinality of F for Max-Det sampling');
    legend('True Bandwidth','Lasso','Garotte','Hard Thresholding');
    grid on;
    
    figure(10);
    plot(10*log10(final_nmsd(1,:)/10),'-','LineWidth',2,'MarkerSize',10);
    hold on;
    plot(10*log10(final_nmsd(2,:)/10),'--','LineWidth',2,'MarkerSize',10);
    hold on;
    plot(10*log10(final_nmsd(3,:)/10),'-.','LineWidth',2,'MarkerSize',10)
    xlabel('Iterations');
    ylabel('Transient NMSD (db) for Max-Det sampling');
    legend('Lasso','Garotte','Hard Thresholding');
    grid on;
    
end
    
            
        
    