load('heinzmaterials'); 
hall = csvread('HeintzAllGenes.csv');
load('heinzcells'); 

%%

 %%
halll = heinzd; 
 %hcells = heinzcells; 
%hallf1 = hall(:,hcells);
hallf1 = halll(:,1:76); 
 hi = find(max(hallf1')>800& min(hallf1')<200); 

hallf= hallf1(hi,:);
j = find(hallf>5000); 
hallf(j) = 5000; 
%%
for j = 1:76
    hallf(:,j)= hallf(:,j)/mean(hallf(:,j)); 
end; 
%%
for j= 1:length(hallf)
     hallf(j,:)= hallf(j,:)/mean(hallf(j,:)); 
end;
%%


%%
h = []; 
Do=[]; 
w=[]; 
for n = 2:35
      
  [w1,h1,D]=nnmf(hallf,n,'replicates',10);
 % [w2,h2]=nmf(hallf,m,1);
  h = [h;h1];
   w = [w,w1];
  Do(n)=D; 
end; 
   %%
   
   dim = length(h); 
dimb = length(hallf);

%%
r =corrcoef([normc(hallf);h]');
%%
rm = r(dimb+1:dimb+dim,:);
%%
rmc = rm(:,1:dimb); 
%%
rmcs = sum(rmc'>.6);
   %%
   k = find(rmcs>(.01*length(hallf))); 
   %%
    hn = []; 
%allaviv = allaviv1(:,ai);
for j = 1:length(h)
    hn(j,:) = h(j,:)/max(h(j,:)); 
end;
   %%
 rh = corrcoef(w(:,k));
 rt =(rh>.9).*(rh<1);  
    imagesc(rt(downr,downr))
   
   %%
     dim = length(h)+1; 
  statusc = []; 
    statuscv = []; 
    hcorr=[]; 
    genes = {}; 
    genel =[]; 
    'starting'
   for k = 1:(dim-1)
       
      r = corrcoef([hallf',h(k,:)']); 
    hcorr(k) = length(find(r(17187,1:17186)>.6)); 
    genes{k} = find(r(17187,1:17186)>.6); 
    genel = [genel,genes{k}];
  end;
  'done'