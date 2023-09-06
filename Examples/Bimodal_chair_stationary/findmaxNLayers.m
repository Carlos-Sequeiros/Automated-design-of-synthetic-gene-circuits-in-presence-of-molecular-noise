function [vals,crdx,crdy]=findmaxNLayers(M,n)
%%King because checks inmediate cells also in diagonal
vals=[0];
crdx=[0];
crdy=[0];
sz=size(M);
i=0;
    for j=n+1:(sz(1)-(n+1))
        for k=n+1:(sz(2)-(n+1))
            if(M(j,k)==max(max(M(j-n:j+n,k-n:k+n))) && M(j,k)~=0 && prod(prod(M(j-n:j+n,k-n:k+n)))~=0)
                if i==0
                    i=i+1;
                    vals(i)=M(j,k);
                    crdx(i)=j;
                    crdy(i)=k;
                else if i>0 && (j~=crdx(i)||k~=crdy(i))
                    i=i+1;
                    vals(i)=M(j,k);
                    crdx(i)=j;
                    crdy(i)=k;
                    end
                end
            end
        end
    end