function [M2, down, across] = clusterarray(a)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
M = []; 
M2 = []; 
Y = pdist(a,'correlation');
%Y = pdist(a,'euclidean');

Y = squareform(Y);
Z = linkage(Y,'single');
[H,T,perm] = dendrogram(Z,0);
down = perm; 
for i = 1:length(a(:,1))
   
     M(i,:) = a(perm(i),:);

end

Y = pdist(M','correlation');


%Y = pdist(M','euclidean');
Y = squareform(Y);
Z = linkage(Y,'single');
[H,T,perm] = dendrogram(Z,0);
across = perm; 

for i = 1:length(a(1,:))
   
     M2(:,i) = M(:,perm(i));

end
