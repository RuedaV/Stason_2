function [ loss ] = QLIKE(h2_proxy, h2_pred)
%QLIKE Summary of this function goes here
%   Detailed explanation goes here

if (h2_pred > 0)
    loss = exp(log(h2_pred) + h2_proxy/h2_pred);
else
    fprintf('%50s\n','-----------------------------------------------------');
    loss = 0;    
end
end

