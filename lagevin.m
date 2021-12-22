clc;
clear;
alpha=1.8;
m=0.01;
%k1=0.0004; 
k2=0.2;
sigma=0.2;
F=0.1;
c=alpha*gamma((1+alpha)/2)/(2^(1-alpha)*pi^(1/2)*gamma(1-alpha/2));
zz=1/(2-alpha);

x0=10;
v0=1;

h=0.0001;
T=10;
nT=T/h;

xf=zeros(1,nT);
qfm=zeros(1,nT);
vf=zeros(1,nT);
qff=zeros(1,nT);
xf(:,1)=x0;
vf(:,1)=v0;
qfm(:,1)=x0;
qff(:,1)=x0;

L=h^(1/alpha)*stblrnd(alpha,0,1,0,1,nT);
for i=1:nT
    v0=vf(:,i);
    x0=xf(:,i);
    xf(:,i+1)=x0+v0*h;
    vf(:,i+1)=v0+(F/m-k2*v0*x0/m)*h+(sigma/m)*x0*L(:,i);
end

for i=1:nT
    q0m=qfm(:,i);
    qfm(:,i+1)=q0m+(F/(k2*q0m))*h+(sigma/k2)*L(:,i)-(sigma^2/(2*k2^2*q0m))*L(:,i)-c*sigma^2/(2*k2^2*q0m)*zz*h;
end

for i=1:nT
    q0f=qff(:,i);
    qff(:,i+1)=q0f+(F/(k2*q0f))*h+(sigma/k2)*L(:,i);
end


t=0:h:T;
figure(1)
plot(t,xf)
hold on
plot(t,qfm)
hold on
plot(t,qff)
xlabel('t'); ylabel('q');

figure(2)
var1=xf-qfm;
var2=xf-qff;
plot(t,var1,'r')
hold on
plot(t,var2,'c')
xlabel('t'); ylabel('\delta q')