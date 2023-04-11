function y=Hammersley(d,k)
%d is the dimension of variable; k is the index of sample
ytemp=[];
for i=1:1:d
    for j=1:1:k
if i==1
    ytemp(i,j)=(j-1)/k;
else
        ytemp(i,j)=Hammersleybasis(i,j);
end    
end
end
y=ytemp;