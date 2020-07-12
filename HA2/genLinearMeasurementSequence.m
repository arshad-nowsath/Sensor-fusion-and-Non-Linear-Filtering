function Y = genLinearMeasurementSequence(X, H, R)
%GENLINEARMEASUREMENTSEQUENCE generates a sequence of observations of the state
% sequence X using a linear measurement model. Measurement noise is assumed to be
% zero mean and Gaussian.
%
%Input:
% X [n x N+1] State vector sequence. The k:th state vector is X(:,k+1)
% H [m x n] Measurement matrix
% R [m x m] Measurement noise covariance
%
%Output:
% Y [m x N] Measurement sequence
% your code here
% Extracting the size of our measurement sequence vector
m = length(R);
N = size(X,2)-1;
% Allocating for our measurement sequence vector
Y = zeros(m,N);
% Drawing samples for our measurement noise vector
r = mvnrnd(zeros(m,1), R, N)';
% Generating the measurement sequence using the linear measurement model
for i=1:N
 Y(:,i) = H*X(:,i+1)+r(:,i);
end
end
