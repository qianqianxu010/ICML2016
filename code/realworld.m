tic
load universitydata.mat;
method=2;
knockoff_plus=0;
q =0.1;

data=univer_qqxu;
[compare,~]=size(data);
user = data(:,1);
y = data(:,4);
data = data(:,2:3);

user_id = unique(user);
A = zeros(compare,length(user_id));
p = length(user_id);
for i = 1:p
    A(:,i) = (user==user_id(i));
end
n=max(data(:));
i = [1:compare,1:compare];
j = data(:);
k = [ones(1,compare),-ones(1,compare)];
d = sparse(i,j,k,compare,n);
X = (eye(compare) - d*((d'*d)\d'));
Y = X*y;
X = X*A;

if method == 1
    myStatistic = @(X, X_ko, y) knockoff.stats.lassoSignedMax(X, X_ko, y);
elseif method == 2
    warning off
    myStatistic = @(X, X_ko, y) IssStatistic(X, X_ko, y,p);
elseif method == 3
    myStatistic = @(X, X_ko, y) knockoff.stats.forwardSelection(X, X_ko, y);
elseif method == 4
    myStatistic = @(X, X_ko, y) knockoff.stats.forwardSelectionOMP(X, X_ko, y);
end
[S,W,info] = knockoff.filter(X, Y, q, 'Statistic', myStatistic,'Knockoffs', 'equi','Normalize',true);
if knockoff_plus == 1;
    S = knockoff.selectVars(W, q, 'knockoff+');
end
a=user_id(S);

  res = y - A(:,S)*((A(:,S)'*A(:,S))\(A(:,S)'*y));
score = (d'*d)\(d'*res);
score = score - mean(score);


b=a';
len=length(a)
c=zeros(len,3);
dd=[];
aa=[];
for i=1:len
    ind=find(univer_qqxu(:,1)==b(i));
    com=length(ind);
    ind2=find(univer_qqxu(ind,4)==1);
    com2=length(ind2);
    ind3=find(univer_qqxu(ind,4)==-1);
    com3=length(ind3);
    c(i,1)=com2;
    c(i,2)=com3;
    c(i,3)=com-com2-com3;
    if com2==0 || com3==0
        dd=[dd
            a(i)];
    else
        aa=[aa
            a(i)];
    end
end

path = [zeros(p,1),info.hist_u(1:p,:)];
t = [0,info.hist_t];
for  i=1:340
    if ismember(i,dd)
        type = 'r';
    else if ismember(i,aa)
            type = 'b';
        else
            type = 'g';
        end
    end
    stairs(t,path(i,:),type);
    hold on;
end
xlabel('t','fontsize',24)
ylabel('\gamma','fontsize',24)
 

