clc;
clear;
prob=6.21;
%% Initialize polynomial equations for prob 6.6
if prob==6.6
    fprintf('Problem 6.6\n');
    m=5;                % define the number of equations
    x=[-2 -1 0 1 2];    % input x
    y=[37 7 5 13 61];   % input y=f(x)
    a=zeros(m,m);
    b=zeros(m,1);
    r=b;
    for i=1:m
        for j=1:m
            a(i,j)=x(i)^(j-1);
        end
        b(i,1)=y(i);
    end
    fprintf('Solved by A\\B\n');
    c=[a,b]
    r=a\b
end
%% Initialize coefficient matrix for prob 6.17
if prob==6.17
    fprintf('\nProblem 6.17\n');
    m=4;
    a=zeros(m,m);
    b=zeros(m,1);
    r=zeros(m,1);
    a=[ 2 1 0 6
        5 2 0 0
        0 7 2 2 
        0 0 8 9];
    b=[ 64
        37
        66
        104];
    % after pivoting
%     a=[ 5 2 0 0
%         0 7 2 2 
%         0 0 8 9
%         2 1 0 6];
%     b=[ 37
%         66
%         104
%         64];
    fprintf('Solved by A\\B\n');
    c=[a,b]
    r=a\b
end
%% Initialize coefficient matrix for prob 6.21
if prob==6.21
    fprintf('\nProblem 6.21\n');
    m=3;
    a=zeros(m,m);
    b=zeros(m,1);
    r=zeros(m,1);
    a=[ 8 1 2 
        1 9 1
        2 3 7];
    b=[ 29
        34
        48];
    fprintf('Solved by A\\B\n');
    c=[a,b]
    r=a\b
end
%% Initialize coefficient matrix for prob 6.24
if prob==6.24
    fprintf('\nProblem 6.24\n');
    m=9;
    tnorth=40;
    tsouth=100;
    twest=50;
    teast=20;
    a=zeros(m,m);
    b=zeros(m,1);
    r=zeros(m,1);
    a=[ -4 1 0 1 0 0 0 0 0
        1 -4 1 0 1 0 0 0 0
        0 1 -4 0 0 1 0 0 0
        1 0 0 -4 1 0 1 0 0
        0 1 0 1 -4 1 0 1 0
        0 0 1 0 1 -4 0 0 1
        0 0 0 1 0 0 -4 1 0
        0 0 0 0 1 0 1 -4 1
        0 0 0 0 0 1 0 1 -4];
    b=-[40+50;40;20+40;50;0;20;50+100;100;100+20];
    fprintf('Solved by A\\B\n');
    c=[a,b]
    r=a\b
end
%% Initialize coefficient matrix for Example 6.2
if prob==6.2e3
    fprintf('\nProblem Example 6.2\n');
    x0=0;
    x1=1;
    T0=1;
    T1=0.5;
    const_G=5;
    m=4;
    dx=(x1-x0)/(m+1);
    %r=zeros(m,1);
    %a=zeros(m,m);
    %b=zeros(m,1);
    p=-(2+const_G*dx^2);
    a=[ p 1 0 0 
        1 p 1 0 
        0 1 p 1 
        0 0 1 p];
    b=-[T0;0;0;T1];
    fprintf('Solved by A\\B\n');
    c=[a,b]
    r=a\b
end
%% Gaussian
fprintf('\nSolved by Gauseian Elimination\n');
c=[a,b]
r=zeros(m,1);
for i=2:m
    for j=i:m
        c(j,:)=c(j,:)-c(i-1,:)/c(i-1,i-1)*c(j,i-1);
    end
end
c
r(m)=c(m,m+1)/c(m,m);
for i=m-1:-1:1
    sum=c(i,m+1);
    if i~=m
        for j=i+1:m
            sum=sum-r(j)*c(i,j);
        end
        
    end
    r(i)=sum/c(i,i);    
end
r
%% Gauss Jordan
fprintf('\nSolved by Gauss Jordan\n');
c=[a,b]
r=zeros(m,1);
for i=1:m
    if i~=1
        for j=i:m
            c(j,:)=c(j,:)-c(i-1,:)/c(i-1,i-1)*c(j,i-1);
        end
    end
    c(i,:)=c(i,:)/c(i,i);
end
c
for i=m:-1:1
    if i~=m
        for j=i:-1:1
            c(j,:)=c(j,:)-c(i+1,:)/c(i+1,i+1)*c(j,i+1);
        end
    end
    c(i,:)=c(i,:)/c(i,i);
end
c
%% Crout
fprintf('\nSolved by Crout''s Method\n');
c=[a,b]
u=zeros(m,m+1);
l=zeros(m,m);
r=zeros(m,1);
[ai, aj] = size(c);
       for i = 1:ai
           l(i, 1) = c(i, 1);
           u(i, i) = 1;
       end
       for j = 2:aj
       u(1, j) = c(1, j) / l(1, 1);
       end
       for i = 2:ai
           for j = 2:i
               l(i, j) = c(i, j) - l(i, 1:j - 1) * u(1:j - 1, j);
           end
           
           for j = i + 1:aj
               u(i, j) = (c(i, j) - l(i, 1:i - 1) * u(1:i - 1, j)) / l(i, i);
           end
       end
l
u
%% LU decomposition
fprintf('\nSolved by LU decomposition\n');
c=[a,b]
[l,u,p]=lu(a)
r1=l\(p*b)
r2=u\r1

%% Jacobi
fprintf('\nSolved by Jacobi\n');
c=[a,b]
r=zeros(m,1);
rold=zeros(m,1);
per1=1;
itr=0;
while(per1>=1)
	itr=itr+1;
	per1=0;
    for i=1:m
        r1=r;
        sum=0;
        for j=1:m
            if (j~=i)  
            sum=sum+a(i,j)*rold(j,1);
            end
        end
        b(i,1);
        sum;
        r(i,1)=(b(i,1)-sum)/a(i,i);    
        eps=abs((r(i,1)-r1(i,1)));
        if(eps>0.001)
            per1=per1+1; 
        end    
    end 
    rold=r1;
end
itr
r


%% Gauss Seidel
fprintf('\nSolved by Gauss Seidel\n');
% a=[4 3 0 0;
%    2 6 1 0;
%    0 3 8 5;
%    0 0 2 9];
% c=[2;5;7;1];
% [n,m]=size(a);
% %[n1,1]=size(c)
% x=[0;0;0;0];
c=[a,b]
r=zeros(m,1);
rold=zeros(m,1);
per1=1;
itr=0;
while(per1>=1)
    itr=itr+1;
    per1=0;
    for i=1:m
        r1=r;
        sum=0;
        for j=1:m
            if (j~=i)  
            sum=sum+a(i,j)*r1(j,1);        
            end
        end
        r(i,1)=(b(i,1)-sum)/a(i,i);    
        eps=abs((r(i,1)-r1(i,1))/r(i,1));
        if(eps>0.001)
            per1=per1+1;    
        end    
    end    
end
itr
r
