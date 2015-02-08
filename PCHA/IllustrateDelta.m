% This script was used to generate Figure 4 in the paper 
% Morten Mørup and Lars K. Hansen, "Archetypal
% Analysis for Machine Learning and Data Mining", submitted NeuroComputing
% 2011

clear all;
N=1000;     % Number of observations
tresh=0.8;  % Level of truncation of simplex
DD=3;       % Dimensionality of simplex

% Generate synethetic data
XC=[cos(0) cos(2*pi/3)  cos(2*pi/3*2); sin(0) sin(2*pi/3) sin(2*pi/3*2)];
S=-log(rand(DD,N));
S=S./repmat(sum(S),DD,1);
[I,J]=find(S>tresh);
S(:,J)=[];
NN=size(S,2);
% Add noise with standard deviation sigma
sigma=0.0;
X=XC*S+sigma*randn(2,NN);

% Plot the generated data
figure;
hold on;
plot(X(1,:),X(2,:),'.');

% Estimate the AA/PCH model using PCHA.m
noc=3;  % Number of components
ms2=1;  % Widht of line in generated plot
delta=[0 0.25 0.5]; % values of \delta

colors={'black','red','green'}; % Mark in seperate colors the 3 different PCHA solutions
for k=1:3
    [XC,S,C,SSE,varexpl]=PCHA(X,noc,1:size(X,2),1:size(X,2),delta(k));   
    line(XC(1,[1 2 3 1]),XC(2,[1 2 3 1]),'color',colors{k},'linewidth',1.5*ms2,'LineStyle','-');
    axis off;
    axis equal;
end
legend({'observations','\delta=0','\delta=0.25','\delta=0.5'})