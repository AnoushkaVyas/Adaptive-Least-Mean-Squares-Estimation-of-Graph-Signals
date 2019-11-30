%Maxdet sampling
function [S,D]=maxdet(M,F,T,G)
    S=zeros(M,1);  %sampling set
    Uf=G.U(:,1:M);
    D=zeros(T,T);
    C=zeros(T,T);
    visited=zeros(T,1);
    j=1;
    while nnz(S)<F
        max=-9223372036854775808;
        index=0;
        for i=1:T
            if visited(i) ==0
                C(i,i)=1;
                e=eig(Uf'*C * Uf);
                o=1;
                u=zeros(1,nnz(e));
                for p=1:M
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
    
            
                
                