function [nu,propensity]=rate_ctes(j,n,H,K,reaction_constants)
 % reactions asociated 12 reactions 6 species
        % 0     -> mRNA1      : rc=c1(P3)*k1_1
        % mRNA1 -> mRNA1 + P1 : rc=k2_1
        % mRNA1 -> 0          : rc=gamma1_1
        % P1    -> 0          : rc=gamma2_1
        % 0     -> mRNA2      : rc=c2(P1)*k1_2
        % mRNA2 -> mRNA2 + P2 : rc=k2_2
        % mRNA2 -> 0          : rc=gamma1_2
        % P2    -> 0          : rc=gamma2_2
        % 0     -> mRNA3      : rc=c3(P2)k1_3
        % mRNA3 -> mRNA3 + P3 : rc=k2_3
        % mRNA3 -> 0          : rc=gamma1_3
        % P3    -> 0          : rc=gamma2_3
        
        %u_tempo = 0.25e4; % dimensionless
        %gamma1_1  = 0.01*u_tempo; 
        %gamma2_1  = 0.0004*u_tempo;
        %a_param_1 = 10;
        %b_param_1 = 10;
        %k1_1      = a_param_1*gamma2_1;
        %k2_1      = b_param_1*gamma1_1;
        %gamma1_2  = 0.01*u_tempo; 
        %gamma2_2  = 0.0004*u_tempo;
        %a_param_2 = 10;
        %b_param_2 = 5;
        %k1_2      = a_param_2*gamma2_2;
        %k2_2      = b_param_2*gamma1_2;
        %H_param_2 = -4;
        %k_param_2 = 70;
        %epsilon_2 = 0.1;
        nonzero_elements=sum(H~=0);
        disp('entra')
        c0=zeros(1,j);
        c0=c0+reaction_constants;
        basicStoichiometricMatrix=[1 0;
            0 1;
            -1 0;
            0 -1];
        S_mat=zeros(j,n);
        for i=1:n/2
            S_mat(4*(i-1)+1:4*i,2*(i-1)+1:2*i)=basicStoichiometricMatrix;
        end
        nu=S_mat;
        epsilon=cell(1,3);
        for i=1:size(H,1)
            epsilon{i}=ones(1,2^numel(find(H(i,:))))*0.15;
            epsilon{i}(end)=1;
        end
        propensity=@(x)[f_cx_v2(x(2:2:end),H,K,epsilon,1)*c0(1),...
            x(1)*c0(2),x(1)*c0(3),x(2)*c0(4),...
            f_cx_v2(x(2:2:end),H,K,epsilon,2)*c0(5),x(3)*c0(6),...
            x(3)*c0(7),x(4)*c0(8),...
            f_cx_v2(x(2:2:end),H,K,epsilon,3)*c0(9),x(5)*c0(10),...
            x(5)*c0(11),x(6)*c0(12)];
end