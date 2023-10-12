
fid = fopen('/home/siavash/Workspace/geoid/INPUT/RefGeo/eigen5c2');
s=fgets(fid);

n=0;
m=0;

while s>0
    zwi = str2num(s);
    if m==0 
        cnm(n+1, m+1) = zwi(1);
        snm(n+1, m+1) = 0.0;
    else
        cnm(n+1, m+1) = zwi(1);
        snm(n+1, m+1) = zwi(2);
    end
    if m==n
        n=n+1;
        m=0;
    else
        m=m+1;
    end
    s=fgets(fid);
end
fclose(fid);

fid = fopen('/home/siavash/Workspace/geoid/INPUT/RefGeo/eigen5c2_new','w');
l_max=size(cnm,1)-1;

fprintf(fid,'%6s %6s %16s %16s \n','n','m','cnm', ...
        'snm');
for i=0:l_max
    for j=0:i
        fprintf(fid,'%6i %6i %16.5e %16.5e \n',i,j,cnm(i+1,j+1), ...
        snm(i+1,j+1) );
    end
end

fclose(fid);

