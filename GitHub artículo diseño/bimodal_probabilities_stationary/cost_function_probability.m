function f = cost_function_probability(x,domains)
Discretization_intervals=400;
Domain_limit=200;
Snapshots_num=1;
discriminationFactor=1/12;
H_mat=[x(5) x(6); x(7) x(8)].*[x(9) x(10); x(11) x(12)];%x(13:16) holds -1 0 or 1 in order to decide network topology
indreg=cell(1,2);
for i=1:numel(indreg)
    indreg{i}=find(H_mat(i,:));
end
mat=[2,x(5:8).*x(9:12),x(1:4).*abs(x(9:12))];
for i=1:2
    auxmat=0.1*ones(1,2^numel(indreg{i})+1);
    auxmat(1)=2^numel(indreg{i});
    auxmat(end)=1;
    mat=[mat auxmat];
end
mat=[mat 10 200 25 1 10 200 25 1];
file_content={mat, 0, 'Hill'};
writecell(file_content,'dato.txt','Delimiter','tab')
fileID=fopen('malla.txt','w');
    fprintf(fileID,['2 0 ' num2str(Domain_limit) ' ' num2str(Discretization_intervals)...
        ' 0 ' num2str(Domain_limit) ' ' num2str(Discretization_intervals)...
        ' 0 5 1000 ' num2str(Snapshots_num)]);
fclose(fileID);
[flag,output]=system('"../parSELANSI-C.exe" dato.txt malla.txt opciones.txt');
fileID=fopen('test_results.txt');
datos_SELANSI=fread(fileID,'double');
fclose(fileID);
datos_SELANSI=reshape(datos_SELANSI,numel(datos_SELANSI)/Snapshots_num,Snapshots_num);
pdf=reshape(datos_SELANSI(3:3:end,1),Discretization_intervals+1,Discretization_intervals+1);
f=0;
domain_probs=[0.50 0.30 0.20];
for i=1:numel(domains)
    f=f+(domain_probs(i)-0.25*trapz(trapz(double(domains{i}).*pdf))).^2;
end
end