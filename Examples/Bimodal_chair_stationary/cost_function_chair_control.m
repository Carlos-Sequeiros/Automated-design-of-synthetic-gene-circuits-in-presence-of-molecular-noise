function f = cost_function_chair_control(x)
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

[value,cordx,cordy]=findmaxNLayers(pdf,2);
    if size(cordx,2)<2
        f=1e10;
        return
    end
        X=reshape(datos_SELANSI(1:3:end,1),Discretization_intervals+1,Discretization_intervals+1);
        Y=reshape(datos_SELANSI(2:3:end,1),Discretization_intervals+1,Discretization_intervals+1);
        cordxTrue=zeros(size(cordx));
        cordyTrue=zeros(size(cordy));
        try
            for i=1:size(cordx,2)
                cordxTrue(i)=X(cordx(i),cordy(i));
                cordyTrue(i)=Y(cordx(i),cordy(i));
            end
        catch
            keyboard
        end
        chairHolder=zeros(Discretization_intervals+1,Discretization_intervals+1);
        [gx,gy]=gradient(pdf.*(pdf>max(pdf(:))*discriminationFactor),0.5,0.5);
        [gxx,gxy]=gradient(gx,0.5,0.5);
        [gyx,gyy]=gradient(gy,0.5,0.5);
        for i=1:Discretization_intervals+1
            for j=1:Discretization_intervals+1
                if sign(prod(eig([gxx(i,j),gxy(i,j);gyx(i,j),gyy(i,j)])))==-1
                    chairHolder(i,j)=-0.5*log(gx(i,j)^2+gy(i,j)^2);
                end
            end
        end
        chairHolder(isinf(chairHolder))=0;
        [~,I]=max(chairHolder(:));
        %[~,Xch,Ych]=findmaxNLayers(chairHolder,2);
%         if (Xch(1)==0||Ych(1)==0)
%             f=10^5;
%             return
%         end
        Vch=pdf(I());
        leftMaximum=[cordxTrue(find(cordxTrue==min(cordxTrue))) cordyTrue(find(cordxTrue==min(cordxTrue)))];
        rightMaximum=[cordxTrue(find(cordxTrue==max(cordxTrue))) cordyTrue(find(cordxTrue==max(cordxTrue)))];
        valores=[value(find(cordxTrue==min(cordxTrue))) value(find(cordxTrue==max(cordxTrue)))];
        fprintf('The number of maximums of the function is: %d at (%3.1f, %3.1f) and (%3.1f, %3.1f)\n',size(cordxTrue,2),leftMaximum(1),leftMaximum(2),rightMaximum(1),rightMaximum(2))
        f=1/norm(leftMaximum-rightMaximum)+((3*valores(1)-valores(2))/valores(2))^2+((10*Vch(1)-valores(2))/valores(2))^2;
end