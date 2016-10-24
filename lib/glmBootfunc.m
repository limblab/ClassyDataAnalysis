function [coeffs] = glmBootfunc(X,y,propName,propVal)
% GLMBOOTFUNC function used for bootstrapping the GLM, returning only
% the coefficients to save memory
%   Argument order is same as expected by fitglm

model = fitglm(X,y,propName,propVal);

coeffs = model.Coefficients.Estimate;

end