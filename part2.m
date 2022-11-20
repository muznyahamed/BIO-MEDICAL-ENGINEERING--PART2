%% Q 2

clear;
close all;
clc;

figure(1) 
of = randi([1 10],5,5);
subplot(2,5,1);imshow(of,[0 maximum(of,[],'all')]);
title('Orginal image');


proj_angles = [0,-45,45,90];
proj = radon(of,proj_angles);
subplot(2,5,2);imshow(proj',[0 maximum(proj,[],'all')]);
title('4 Forward projections');


rec = ones(size(of));

subplot(2,5,4);imshow(rec,[0 maximum(rec,[],'all')]);
title('initial image with 0 iterations');
image_ones = ones(size(proj));
sense = iradon(image_ones,proj_angles,'none',5);
col_o_f = of(:);
for i=1:6
    proj_r = radon(rec,proj_angles);
    ratio = proj./(proj_r + 0.0001);
    subplot(2,5,3);imshow(ratio',[0 maximum(ratio,[],'all')]);
    title('ratio of real vs reconstructed');
    bp_ratio = iradon(ratio,proj_angles,'none',5);
    rec = (rec.*bp_ratio)./sense;
    subplot(2,5,i+4);imshow(rec,[0 maximum(rec,[],'all')]);
    title("Iteration no."+i+"i")
    col_rec = rec(:);
    error = sum ((col_o_f-col_rec).^2);
    rmse1(i) = sqrt(error/25);
end  
figure(2)
i = 1:1:6;
plot(i,rmse1,'-o','MarkerEdgeColor','red')
ylim([1.5 2.6])
xlim([0 7])
title('NRMS Error of 6 iterations');

ylabel('NRMSE');
xlabel('Iterations');


%% Projections for Linear Algebra
azi_angles =[0,-45,45,90];
Projection1=radon(of,azi_angles);
imshow(Projection1');



of = randi([1 10],5,5);
imagesc(of);
P1 =zeros(1,11);
P2 =zeros(1,11);
P1 = sum(of,1);
P2 = sum(of,2)';
% P3=zeros(1,9);
summ = zeros(1,1);
for i=1:5
    for j=1:i
        k=j;
        m=5-i+j;
        summ = summ + of(k,m);
    end
    P3(i)=summ;
end    

for i=1:5
    for j=1:i
        m=j;
        k=5-i+j;
        summ = summ + of(k,m);
    end
    P3(10-i)=summ;
    summ=0;
end 

P4=zeros(1,9);
for i=1:5
    for j=1:i
        k=j;
        m=i+1-j;
        summ = summ + of(k,m);
    end
    P4(i)=summ;
    summ=0;
end    

for i=1:5
    for j=1:i
        k=j;
        m=5-i+j;
        summ = summ + of(k,m);
    end
    P4(10-i)=summ;
    summ=0;
end 

%% Back Projection using Linear Algebra

syms X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 X14 X15 X16 X17 X18 X19 X20 X21 X22 X23 X24 X25 ;
X = [X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 X14 X15 X16 X17 X18 X19 X20 X21 X22 X23 X24 X25 ];
clear X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 X14 X15 X16 X17 X18 X19 X20 X21 X22 X23 X24 X25;


for i= 1:5
    eqn(i) = X(i)+X(i+5)+X(i+10)+X(i+15)+X(i+20)-P1(i);
end


for i= 1:5
    eqn(i+5) = X(5*(i-1)+1)+X(5*(i-1)+2)+X(5*(i-1)+3)+X(5*(i-1)+4)+X(5*(i-1)+5)-P2(i);
end

for i= 1:5
    maximum=i*5;
    maximum=6-i;
    now=maximum;
    summ=X(now);
    for j=1:i-1
    if i>1    
    now = maximum +j*(maximum-maximum)/(i-1);   
    summ = summ+X(now);
    end
    end
    eqn3(1,i) = summ -P3(i) ;
end
for i= 1:4
    maximum=25-i;
    maximum=5*i+1;
    now=maximum;
    summ=X(now);
    for j=1:4-i
    if 4-i>0    
    now = maximum +j*(maximum-maximum)/(4-i);   
    summ = summ+X(now);
    end
    end
    eqn3(1,i+5) = summ-P4(i) ;
end



for i= 1:5
    maximum=5*(i-1)+1;
    maximum=i;
    now=maximum;
    summ=X(now);
    for j=1:i-1
    if i>1    
    now = maximum +j*(maximum-maximum)/(i-1);   
    summ = summ+X(now);
    end
    end
    eqn4(1,i) = summ -P4(i) ;
end
for i= 1:4
    maximum=21+i;
    maximum=5*(i+1);
    now=maximum;
    summ=X(now);
    for j=1:4-i
    if 4-i>0    
    now = maximum +j*(maximum-maximum)/(4-i);   
    summ = summ+X(now);
    end
    end
    eqn4(1,i+5) = summ -P4(i) ;
end



equations = horzcat(eqn,eqn3,eqn4);


F=@(X)[double(equations(:))];
fun=@(u) F(u(1),u(2),u(3));
x0=0;
[val1,fval1]=fsolve(F,x0),
[u2,~,fval2]=lsqnonlin(fun,x0)
    
[A,B] = equationsToMatrix(equations, X);
Xvalues = mldivide(A,B);

sol = solve(equations,X);
X1val = sol.X1;