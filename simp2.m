function s=simp2(f,a,b,n)

for i=1:10
    
h=(b-a)/n;
m=n/2;

s1=0;
s2=0;
for k=1:m
    x=a+h*(2*k-1);
    s1=s1+feval(f,x);
end
for k=1:(m-1)
    x=a+2*h*k;
    s2=s2+feval(f,x);
end
s = h*(feval(f,a)+feval(f,b)+4*s1+2*s2)/3;

fprintf('n=%4g        Approx. Value of Integral=%.4f\n',n,s);

n = 2*n;

end

    


