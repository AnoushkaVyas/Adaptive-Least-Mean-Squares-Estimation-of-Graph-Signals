%Maxdet sampling without known bandwidth
function D=maxdetext(f,T,F,G)
    S=zeros(f,1);
    Uf=zeros(T,f);
    k=1;
    for i=1:T
        if F(i) ~= 0
            Uf(:,k)=G.U(:,F(i));
            k=k+1;
        end
    end
    D=zeros(T,T);
    C=zeros(T,T);
    visited=zeros(T,1);
    j=1;
    while nnz(S)<f
        max=-9223372036854775808;
        index=0;
        for i=1:T
            if visited(i) ==0
                C(i,i)=1;
                e=eig(Uf'*C * Uf);
                o=1;
                u=zeros(1,nnz(e));
                for p=1:f
                    if e(p) ~= 0
                        u(o)=e(p);
                        o=o+1;
                    end
                end     
                temp=prod(u);
                if (temp > max)
                    max=temp;
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
    
            
                
                