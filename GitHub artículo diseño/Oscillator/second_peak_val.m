function  val = second_peak_val(x)
%SECOND_PEAK_VAL Summary of this function goes here
%   Detailed explanation goes here
	%K=reshape(x(1:9),3,3);
	%reaction_constants=reshape(x(10:21),3,4);
	%H=reshape(x(22:30),3,3);
    K=[0 x(1) 0; 0 0 x(1); x(1) 0 0];
    reaction_constants=repmat([x(2:4) 1],1,3);
    H=[0 x(5) 0; 0 0 x(5); x(5) 0 0];
% 	if(sum(sum(H~=0))==0) | (sum(sum((H~=0)==eye(2)))==4) | (((x(11)==0)&(x(12)==0))|((x(13)==0)&(x(14)==0)))
% 		fprintf('Invalid Topology, assuming high cost function')
% 		val = 10;
% 		return;
% 	end
    [nu,propensity]=rate_ctes(12,6,H,K,reaction_constants);
    Tgrid=0:0.1:1e2;
    simulation=SSA_simulation(propensity,nu,5*ones(1,6),Tgrid,1);
    correlation=zeros(1,1000);
    for i=1:size(correlation,2)
        correlation(i)=trapz(simulation{1}(2,:).*circshift(simulation{1}(2,:),[0 i-1]));
    end
    correlation=correlation-mean(correlation);
    correlation=correlation/max(correlation);
    End_condition=false;
    index=2;
    while(~End_condition && index<size(correlation,2))
        if(correlation(index)>correlation(index-1) &&...
                correlation(index)>correlation(index+1))
            val=-correlation(index);
            End_condition=true;
            fprintf('K=%f, H=%i, km=%f, kx=%f, gamma_m=%f, valor=%f\n',x(1),x(5),x(2),x(3),x(4),-val)
        end
        index=index+1;
    end
end