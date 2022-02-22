clear all
close all
clc

Ans = inf;
% scalars
fnum = input('please enter number of test function : ');
Ni = 50;
Nj = 5;
D = fminmax(fnum);
[~,xmin,~] = fminmax(fnum);
[~,~,xmax] = fminmax(fnum);
iter = 2000;
l1 = 0.05;
l2 = 100;
l3 = 0.01;
t = 1;
e = 0.05;
bta = 1;
alfa = 1;
K = 1;
M1 = 0.1;
M2 = 0.2;
Ttta = 298.15;
gama = 0;


% vectors
c = l3 * rand(1 , Nj);
T = zeros(1 , iter);
Fbest = zeros(1 , iter);

% matrixes
h = zeros(iter , Nj);
h(1,:) = l1 * rand(1 , Nj);
p = l2 * rand(Ni , Nj);
s = zeros(Ni , Nj);
x = zeros(Ni , Nj , D);
F = zeros(Ni , Nj);
Xbest = zeros(iter , D);
Xi_best = zeros(Ni , D);
Fi_best = zeros(Ni);
x = zeros(Ni , Nj , D);
for i=1:Ni
    for j=1:Nj
        for d=1:D
            x(i,j,d) = xmin + rand * (xmax - xmin);
        end
    end
end

for i=1:Ni
    for j=1:Nj
        F(i,j) = f(x(i,j,:),fnum);
    end
end
% finde best solution
minvar = F(1,1);
tempbest = zeros(1,1,D);
for i=1:Ni
   for j=1:Nj
       if F(i,j) < minvar
           tempbest = x(i,j,:);
           minvar = F(i,j);
       end
   end
end
Xbest(t,:) = tempbest;
Fbest(t) = minvar;
minvar = F(1,1);
jsave = 1;
for i=1:Ni
    minvar = F(i,1);
   for j=1:Nj
       if F(i,j) < minvar
           minvar = F(i,j);
           jsave = j;
       end
   end
   Xi_best(i,:) = x(i,jsave,:);
   Fi_best(i) = f(Xi_best(i,:),fnum);
end
   
while t<iter
%      for each agnet update position by eq10;
    gama = bta * exp( (Fbest(t)+e) / (F(i,j)+e) );
    if rand < 0.5
        flag = -1;
    else
        flag = 1;
    end
    r = rand;
    for i=1:Ni
        for j=1:Nj
            for d=1:D
                x(i,j,d) = x(i,j,d) + ...
                flag * r * gama * (Xi_best(i,d)-x(i,j,d)) + ...
                flag * r * alfa * ( s(i,j) * Xbest(t,d) - x(i,j,d) );
            end
        end 
    end
%      Update Henry’s coefficient of each gas type using eq8
    T(t) = exp(-t/iter);
    for j=1:Nj
        h(t+1,j) = h(t,j)*exp(c(j)*(1/T(t) - 1/Ttta));
    end
%      Update solubility of each gas using eq9
    for i=1:Ni
        for j=1:Nj
            s(i,j) = K * h(t+1,j) * p(i,j);
        end
    end
%      Rank and select the number of worst agents using eq11
    Nw = Ni * (rand * (M1-M2) + M1);
    Nw = fix(Nw);
    [~,worstfinder] = sort(Fi_best);
    Nw = Ni - Nw;
    Nw = worstfinder(Nw);
%      Update the position of the worst agents using eq12
    for j=1:Nj
        for d=1:D
            x(Nw,j,d) = xmin + rand * (xmax - xmin);
        end
    end
%     update Fittnes of each agnet
    for i=1:Ni
        for j=1:Nj
            F(i,j) = f(x(i,j,:),fnum);
        end
    end
    
%      Update the best gas Xi,best, and the best search agent Xbest.
   minvar = F(1,1);
   tempbest = zeros(1,1,D);
   for i=1:Ni
       for j=1:Nj
           if F(i,j) < minvar
               tempbest = x(i,j,:);
               minvar = F(i,j);
           end
       end
   end
   Xbest(t,:) = tempbest;
   Fbest(t) = minvar;
%    Fbest(t)
   minvar = F(1,1);
   jsave = 1;
   for i=1:Ni
       minvar = F(i,1);
       for j=1:Nj
           if F(i,j) < minvar
               minvar = F(i,j);
               jsave = j;
           end
       end
       Xi_best(i,:) = x(i,jsave,:);
   end
    disp(['t = ' num2str(t) ' Fbest = ' num2str(Fbest(t))]);
   Ans = Fbest(t);
    t = t + 1;
    
end
    disp(['t = ' num2str(t) ' Fbest = ' num2str(Ans)]);
