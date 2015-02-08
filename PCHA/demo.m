%% Example of how to use PCHA.m, PCHAkernel.m and PCHAsparse.m
X=rand(15,100); % random data

noc=10; % Number of archetypes

U=1:size(X,2); % Entries in X used that is modelled by the AA model
I=1:size(X,2); % Entries in X used to define archetypes
% if two expensive to useall entries for I find N relevant observations by
% the following procedure:
% N=100;
% I=FurthestSum(X,N,ceil(rand*size(X,2)));

opts.maxiter=1000;
opts.conv_crit=1e-6;

% Use PCHA.m
[XC,S,C,SSE,varexpl]=PCHA(X,noc,I,U,0.1,opts);

% Use PCHAkernel.m
[S,C,SSE,varexpl]=PCHAkernel(X'*X,noc,I,U,0.1,opts);

% Use PCHAsparse.m such that everything that is zero will be treated as missing
X=sprand(100,1000,0.1); 
[XC,S,C,SSE,varexpl]=PCHAsparse(X,noc,opts);

