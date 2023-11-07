clc;
clear all;

ts = 0.5e-3;
k = 0;
t1 = 5*k;
t2 = 5*k+2.5;
t3 = 5*k+5;

%% r(t) when t:[5k,5k+2.5]
t_vec1 = t1:ts:t2;
r_vec1 = (2/5)*(t_vec1-5*k)-(1/(2*pi))*sin((4/5)*pi*(t_vec1-5*k));
r_vec1 = r_vec1';

%% r(t) when t:[5k+2.5,5k+5]
t_vec2 = t2:ts:t3;
r_vec2 = 2-(2/5)*(t_vec2-5*k)+(1/(2*pi))*sin((4/5)*pi*(t_vec2-5*k));
r_vec2 = r_vec2';

%% r(t)
r_vec = [r_vec1;r_vec2(2:end)];
t_vec = [t_vec1';t_vec2(2:end)';];

plot(t_vec,r_vec);

%% time domain expression of basis function phi1:= 1
phi1 = r_vec;

%% time domain expression of basis function phi2:= (z-1)/ts, time advance shift
phi2 = (1/ts)*(r_vec([2:end,end])-r_vec);

%% time domain expression of basis function phi3 
% initialize phi3[1]=[y1];
y1 = (2/ts^2)*(r_vec(2)-2*r_vec(1));
phi3 = [y1];

%% perform recursion to express phi3
yk_prev = y1;
for k = 2:size(r_vec,1)-1
    yk = (2/ts^2)*(r_vec(k+1)-2*r_vec(k)+r_vec(k-1)) -yk_prev;
    yk_prev = yk;
    phi3 = [phi3;yk];
end

% adjust value of k := size(r_vec,1) to express the last element in phi3, assuming r_vec(k+1)=0
k = k+1;
p = (2/ts^2)*(-r_vec(k)+r_vec(k-1)) -yk_prev;
phi3 = [phi3;p];

%% lifted PHI=[phi1, phi2, phi3]
PHI = [phi1, phi2, phi3];
display(PHI);

%% define system parameters
m = 2;
d = 4;
k = 200;
zeta = 0.01;
wr = 3000;
wa = 4000;
wq = 605;

s = tf('s');
C = 20000*(1+s/50)/(1+s/800);
Pn = 1/(m*s^2+d*s+k); 
P = Pn*(wr^2*(s^2+2*zeta*wa*s+wa^2))/(wa^2*(s^2+2*zeta*wr*s+wr^2)); % plant 

mn = 1.6;
dn = 3.2;
kn = 160;
Pn = 1/(mn*s^2+dn*s+kn); % nominal plant 
% discretization 
Cd = c2d(C,ts,'tustin');
Pnd = c2d(Pn,ts,'zoh');
Pd = c2d(P,ts,'zoh');

% H 
H = -Pd/(1+Pd*Cd );
Hn = -Pnd/(1+Pnd*Cd );
Hss = impulse(H,0:ts:5);
Hss = tril(toeplitz(Hss));
Hnss = impulse(Hn,0:ts:5);
Hnss = tril(toeplitz(Hnss));

% S
S = 1/(1+Pd*Cd);
Sss = impulse(S,0:ts:5);
Sss = tril(toeplitz(Sss));

L = -pinv(Hnss*PHI);

conv = norm(eye(size(L*Hss*PHI))+L*Hss*PHI,2);


iter = 10; % number of iteration

dim = (t2-t1)/ts;
ej = zeros(dim+1,iter);
fj = zeros(dim+1,iter+1);
theta_j = zeros(3,iter);

for j = 1:iter
ej(:,j) = Hss*fj(:,j)+lsim(S,r_vec,t_vec);
theta_j(:,j+1) = theta_j(:,j)+L*ej(:,j);
fj(:,j+1) = PHI*theta_j(:,j+1);
end

