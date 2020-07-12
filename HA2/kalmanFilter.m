function [xHat, P] = kalmanFilter(Y, x_0, P_0, A, Q, H, R)
%KALMANFILTER Filters measurements sequence Y using a Kalman filter.
%
%Input:
% Y [m x N] Measurement sequence
% x_0 [n x 1] Prior mean
% P_0 [n x n] Prior covariance
% A [n x n] State transition matrix
% Q [n x n] Process noise covariance
% H [m x n] Measurement model matrix
% R [m x m] Measurement noise covariance
%
%Output:
% x [n x N] Estimated state vector sequence
% P [n x n x N] Filter error convariance
%
%% Parameters
N = size(Y,2);
n = length(x_0);
m = size(Y,1);
%% Data allocation
xHat = zeros(n,N);
P = zeros(n,n,N);
for k=1:N
 if(k==1)
 [xHat(:,k),P(:,:,k)] = linearPrediction(x_0, P_0, A, Q);
 [xHat(:,k),P(:,:,k)] = linearUpdate(xHat(:,k), P(:,:,k), Y(:,k), H, R);
 else
 [xHat(:,k),P(:,:,k)] = linearPrediction(xHat(:,k-1), P(:,:,k-1), A, Q);
 [xHat(:,k),P(:,:,k)] = linearUpdate(xHat(:,k), P(:,:,k), Y(:,k), H, R);
 end
end
end