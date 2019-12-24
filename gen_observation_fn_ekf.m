function [Z2 Z1]= gen_observation_fn_ekf(model,X,W)

%r/t observation equation
if model.outliers_flg == 1
    if rand<0.1
        DD = model.D_outliers;
    else
        DD = model.D;
    end
else
    DD = model.D;
end


if ~isnumeric(W)
    if strcmp(W,'noise')
        W= DD*randn(size(DD,2),size(X,2));
    elseif strcmp(W,'noiseless')
        W= zeros(size(model.D,1),size(X,2));
    end
end

if isempty(X)
    Z1= [];
    Z2= [];
else %modify below here for user specified measurement model
    P= X([1 3],:);
    Z2(1,:)= atan2(P(1,:),P(2,:)) + W(1,:);   
    Z2(2,:)= sqrt(sum(P.^2)) + W(2,:);  
    Z1(1,:)= (sqrt(sum(P.^2)) + W(2,:)).*cos(-atan2(P(1,:),P(2,:)) + W(1,:)+pi/2);
    Z1(2,:)= (sqrt(sum(P.^2)) + W(2,:)).*sin(-atan2(P(1,:),P(2,:)) + W(1,:)+pi/2);
end
