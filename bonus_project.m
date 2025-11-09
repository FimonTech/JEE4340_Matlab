syms b
R = 5;
L = 5e-3;
f=60;
w = 2*pi*f;
X = w*L;
th = atan(X/R);
a = pi/3;
eq = sin(b-th) == sin(a-th)*exp(R/L*(a-b)/w);
sol = solve(eq,b);
disp(rad2deg(sol));
b = pi;
for i=1:1:1000
    b1 = b + 0.001;
    it_1 = sin(b1-th) -sin(a-th)*exp(R/L*(a-b1)/w);
    b2 = b - 0.001;
    it_2 = sin(b2-th) -sin(a-th)*exp(R/L*(a-b2)/w);
    if (it_1 - it_2 > 0)
        b = b2;
    else
        b = b1;
    end
    if (it_1 < 0.0001 || it_2 < 0.0001)
        break;
    end
end
disp(rad2deg(abs(b)));
