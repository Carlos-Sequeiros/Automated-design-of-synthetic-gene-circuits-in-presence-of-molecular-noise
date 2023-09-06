function cx=f_cx_v2(x,H,K,epsilon,ci)
    n=size(H,2);
    indreg=cell(n,1);
    for i=1:n
        indreg{i}=find(H(i,:));
    end
    RHO=cell(n,n);
    for i=1:n
        for j=1:n
            if H(i,j)>0
                RHO{j,i}= x(j).^H(i,j)./(x(j).^H(i,j)+K(i,j)^H(i,j));
            elseif H(i,j)==0
                RHO{j,i}=1;
            else
                RHO{j,i}=K(i,j)^(-H(i,j))./(x(j).^(-H(i,j))+K(i,j)^(-H(i,j)));
            end
        end
    end
    ON_RHO=cell(n,n);
    for i=1:numel(RHO)
        ON_RHO{i}=1-RHO{i};
    end
    S_OFF_ON=cell(n,2*n);
    S_OFF_ON(:,1:2:2*n-1)=RHO;
    S_OFF_ON(:,2:2:2*n)=ON_RHO;
    % Calculate the c_i(x) for i=ci
    ind_J=length(indreg{ci});
    CX=cell(1,2^ind_J);
    CX{1}=S_OFF_ON{indreg{ci}(1),2*ci-1}; 
    CX{2}=S_OFF_ON{indreg{ci}(1),2*ci};
    if ind_J>1
        for i=2:ind_J
            AAA=CX(1:2^(i-1));
            BBB=S_OFF_ON(indreg{ci}(i),2*ci-1:2*ci);
            [ja,jb] = meshgrid(1:2^(i-1),[1 2]);
            for jj=1:2^i
                CX{jj}=AAA{ja(jj)}.*BBB{jb(jj)};
            end
        end
    end
    cx=0;
    for i=1:2^ind_J
        try
        cx=cx+CX{i}*epsilon{ci}(i);
        catch
            keyboard
        end
    end
end