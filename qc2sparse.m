function H=qc2sparse(qcHFileName)
fid=fopen(qcHFileName,'r');
n=fscanf(fid,'%d',[1 1]);
m=fscanf(fid,'%d',[1 1]);
z=fscanf(fid,'%d',[1 1]);
I=sparse(eye(z));
Z=sparse(zeros(z));
H=sparse([]);
for i=1:m
    IH=sparse([]);
    for j=1:n
        shift=fscanf(fid,'%d',[1 1]);
        if shift==-1
            IH=[IH Z];
        else
            IH=[IH circshift(I, [0 shift])];
        end
    end
    H=[H;IH];
end
fclose(fid);
end