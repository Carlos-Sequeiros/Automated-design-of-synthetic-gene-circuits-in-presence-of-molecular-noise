function simulation=SSA_simulation(propensity,nu,x0,Tgrid,nsimula)
    simulation=cell(nsimula,1);
    for nsim=1:nsimula
        simulation{nsim}=zeros(numel(x0),numel(Tgrid));
        x=x0';
        simulation{nsim}(:,1)=x;
        t=Tgrid(1);
        for i=2:numel(Tgrid)
            while t<Tgrid(i)
                props=propensity(x);
                a0=sum(props);
                r1=rand();
                r2=rand();
                delta_t=1/a0*log(1/r1);
                indicator=0;
                reaction_index=0;
                while indicator<r2*a0
                    reaction_index=reaction_index+1;
                    indicator=indicator+props(reaction_index);
                end
                if reaction_index==0
                    t=Tgrid(i);
                elseif t+delta_t>=Tgrid(i)
                    t=Tgrid(i);
                else
                    t=t+delta_t;
                    x=x+nu(reaction_index,:)';
                end
            end
            simulation{nsim}(:,i)=x;
        end
    end
end