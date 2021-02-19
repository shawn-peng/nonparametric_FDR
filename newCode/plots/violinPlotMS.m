function sigma = violinPlotMS(s1,s2,e)
[I1,ix1]=sort(s1);
%[I2,ix2]=sort(s2);
IX1={};
x= round(min(s1)):3:round(max(s1));
s_I2 = s1-s2; 
s_I2=s_I2(ix1);
sigma=nan(length(x),1);
%bp={};
bp=[];
X=[];
for i = 1:length(x)
   ix=find(I1<x(i)+e & I1>x(i)-e);
   IX1=[IX1,ix];
   sigma(i)= sqrt(var([s_I2(ix);-s_I2(ix)])); 
   %bp=[bp,s_I2(ix)];
   bp=[bp;s_I2(ix)];
   X = [X;repmat(x(i),length(s_I2(ix)),1)];
end
S.bp=bp;
S.X=X;
S.s1=s1;
jsonStr = jsonencode(S);
fid = fopen(['test_search/results/', 'vp' ,'.json'], 'w');
if fid == -1, error('Cannot create JSON file'); end
fwrite(fid, jsonStr, 'char');
fclose(fid);
violinplot(bp,X);
xlabel('S_1')
ylabel('S_1-S_2')
end