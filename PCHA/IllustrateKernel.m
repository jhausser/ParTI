% This script was used to generate Figure 5 in the paper 
% Morten Mørup and Lars K. Hansen, "Archetypal
% Analysis for Machine Learning and Data Mining", submitted NeuroComputing
% 2011

clear all;
close all;

N1=250;      % Number of observations in outer circle
R1=1;        % Radius of outer circle
widthR1=0.3; % Width of outer circle
N2=250;      % Number of observations in inner circle
R2=0;        % Radius of inner circle
widthR2=0.1; % Width of inner circle

% Generate synetic data
theta1=rand(1,N1)
X1=(R1+widthR1*(ones(2,1)*rand(1,N1))).*[cos(2*pi*theta1); sin(2*pi*theta1)];
theta2=rand(1,N2)
X2=(R2+widthR2*(ones(2,1)*rand(1,N1))).*[cos(2*pi*theta2); sin(2*pi*theta2)];
X=[X1 X2];

% Plot the generated data
figure;
hold on;
plot(X(1,:),X(2,:),'.');

% Estimate the AA/PCH model using PCHA.m and PCHAkernel.m
noc=10;     % Number of components
delta=0;    % Value of \delta
colors={'black','red'};
for k=[1 2]
    if k==1
        [XC,S,C,SSE,varexpl]=PCHA(X,noc,1:size(X,2),1:size(X,2));   
    else
        D = squareform(pdist(X'));
        sigma=0.1;        
        [S,C,SSE,varexpl]=PCHAkernel(exp(-D.^2/sigma),noc,1:size(X,2),1:size(X,2));   
        XC=X*C;
    end
    for n=1:noc
        plot(XC(1,n),XC(2,n),'o','color',colors{k});        
        plot(XC(1,n),XC(2,n),'x','color',colors{k});        
    end
    axis off;
    axis equal;
end

