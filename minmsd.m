function S=minmsd(M,T,G,mu,Cv)
    S=zeros(M,1);
    Uf=G.U(:,1:M);
    D=zeros(T,T);
    C=zeros(T,T);
    visited=zeros(T,1);
    j=1;
    while nnz(S)<M
        min=10000000000;
        index=0;
        for i=1:T
            if visited(i) ==0
                C(i,i)=1;
                G=Uf'*C * Cv * C * Uf;
                TT=(eye(M)-(mu*Uf'*C*Uf));
                Q=kron(TT,TT);
                vG=reshape(G,[M*M,1]);
                vI=reshape(eye(M),[M*M,1]);
                temp=vG' * pinv(Q) * vI;
                if (temp < min)
                    min=temp;
                    index=i;
                end
                C(i,i)=0;
            end
        end
        visited(index)=1;
        D(index,index)=1;
        C=D;
        S(j)=index;
        j=j+1;
    end
end