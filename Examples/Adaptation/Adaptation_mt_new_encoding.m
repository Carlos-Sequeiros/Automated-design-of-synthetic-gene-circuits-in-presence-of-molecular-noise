function cost_function = Adaptation_mt_new_encoding( x )
%RESPONSE_DIFERENCE Summary of this function goes here
%   Detailed explanation goes here
Discretization_intervals=128;
Domain_limit=350;
Snapshots_num=40;
H_mat=[1 0 0; x(23:25); x(26:28)];
indreg=cell(1,3);
for i=1:numel(indreg)
    indreg{i}=find(H_mat(i,:));
end
mat=[3,[1 0 0 x(23:28)],[1 0 0 x(1:6)]];
%Epsilon coefficients of input gene
auxmat=0.250*ones(1,2^numel(indreg{1})+1);
auxmat(1)=2^numel(indreg{1});
auxmat(end)=0;
mat=[mat auxmat];
%Epsilon coefficients of intermediary gene
auxmat=zeros(1,2^numel(indreg{2})+1);
auxmat(1)=2^numel(indreg{2});
auxmat(2:end)=x(7:14);
mat=[mat auxmat];
%Epsilon coefficients of output gene
auxmat=zeros(1,2^numel(indreg{3})+1);
auxmat(1)=2^numel(indreg{3});
auxmat(2:end)=x(15:22);
mat=[mat auxmat];

mat=[mat 8 400 25 1 8 400 25 1 8 400 25 1];
file_content={mat, 0, 'Hill'};
writecell(file_content,'dato.txt','Delimiter','tab')
fileID=fopen('malla.txt','w');
    fprintf(fileID,['3 0 ' num2str(Domain_limit) ' ' num2str(Discretization_intervals)...
        ' 0 ' num2str(Domain_limit) ' ' num2str(Discretization_intervals)...
        ' 0 ' num2str(Domain_limit) ' ' num2str(Discretization_intervals)...
        ' 0 80 400 ' num2str(Snapshots_num)]);
fclose(fileID);
[flag,output]=system('"../parSELANSI-C.exe" dato.txt malla.txt opciones.txt');
fileID=fopen('test_results.txt','r');
marginal_dist=fread(fileID,'double');
fclose(fileID);
marginal_dist=reshape(marginal_dist,numel(marginal_dist)/Snapshots_num,Snapshots_num);
expected_proteins=zeros(1,Snapshots_num);
for i=1:40
    expected_proteins(i)=trapz(marginal_dist(1:2:end,i),marginal_dist(1:2:end,i).*marginal_dist(2:2:end,i));
end
derivative_of_expected_proteins=diff(expected_proteins);
if(abs(derivative_of_expected_proteins(20))>1e-3 || abs(derivative_of_expected_proteins(39))>1e-3)
    fprintf('\nWarning, simulation has not converged to stationary state, assuming extremely high cost (1e10)\n')
    cost_function=1e10;
    return
end
O1=expected_proteins(21);
O2=expected_proteins(40);
Opeak=0;
for i=20:size(derivative_of_expected_proteins,2)
    prev_value=derivative_of_expected_proteins(i-1);
    act_value=derivative_of_expected_proteins(i);
    if(abs(prev_value)<1e-3)
        prev_value=0;
    end
    if(abs(act_value)<1e-3)
        act_value=0;
    end
    if(act_value*prev_value<0)
        Opeak=expected_proteins(i);
    end
end
if (~Opeak)
    fprintf('\nNo peak found, assuming extremely high cost (1e10)\n')
    cost_function=1e10;
    return;
end
sensitivity=abs((Opeak-O1)/O1/(0.75-0.25)*0.25);
precision=1/abs((O2-O1)/O1/(0.75-0.25)*0.25);
cost_function=sensitivity^(-2)+100*precision^(-2)+0.01*O2^(-2);
fprintf('\n Cost function value: %f \n',cost_function)
end