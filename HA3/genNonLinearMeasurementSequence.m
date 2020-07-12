function Y = genNonLinearMeasurementSequence(X, h, R)
%GENNONLINEARMEASUREMENTSEQUENCE generates ovservations of the states 
% sequence X using a non-linear measurement model.
%
%Input:
%   X           [n x N+1] State vector sequence
%   h           Measurement model function handle
%   h           Measurement model function handle
%               [hx,Hx]=h(x) 
%               Takes as input x (state) 
%               Returns hx and Hx, measurement model and Jacobian evaluated at x
%   R           [m x m] Measurement noise covariance
%
%Output:
%   Y           [m x N] Measurement sequence
% Your code here
%Size of sequence vector
M=length(R);N=size(X,2)-1;

%space for measurement vector is allocated
Y=zeros(M,N);

%Sampling
r=mvnrnd(zeros(M,1),R,N)';

%generate sequence
for i=1:N
    Y(:,i)=h(X(:,i+1))+r(:,i);
end
end