% initial starting file for modeling the Alzheimer's Disease system of
% equations at 4.1

% rate of healthy tissue production
beta = 15;
lambda(t) = 2*(1-a(t))*(1+beta*a(t));

% functions scaled by N for monomers, oligomers, and immobile aggregates
N = 10^11;
X(t)= w1(t)/N;
Y(t)= w2(t)/N;
Z(t)= w3(t)/N;

% equations at 4.1
dXdt = -k*X^2 - k*X*Y - kstar*X*Z - M(1)*X + lambda(t);
dYdt = 0.5*k*X^2 - k*X*Y - k*Y^2 - kstar*Y*Z - M(2)*Y;
dZdt = 0.5*k*Y^2 + k*X*Y - M(3)*Z;

% rates of polymer productions
function f = polymers(t,w)
% length of shortest immobile aggregates
n = 3;
% summation of polymers of lengths 1:n-1
for j = 1:n-1
    polysum = polysum + w(j);
end
f(1) = -K*w(1)*polysum -Kstar*w(1)*w(n) + lambda - M(1)*w(1);
% oligomer rates can be generated generically using this loop
for s = 2:n-1
    % sum of aggregates of length s created from smaller polymers chaining
    for i = 1:s
        olisum = olisum + w(i)*w(s-i);
    end
    f(s) = K/2*olisum - K*w(s)*polysum - Kstar*w(s)w(n) - M(s)*w(s);
end
% sum of immobile aggregates of length >= n from polymer chaining
for i = 1:n-1
    for j = i:n-1
        if i+j >= n
            aggresum = aggresum + w(i)*w(j);
        end
    end
end
f(n) = K/2*aggresum - M(n)*w(n);
f = f(:);
end





