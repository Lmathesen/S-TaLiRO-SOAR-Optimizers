function [run, history] = BayesOpt_Taliro(inpRanges,opt)
% UR_Taliro - Performs random sampling in the state and input spaces.
%
% USAGE:
%   [run, history] = UR_Taliro(inpRanges,opt)
%
% INPUTS:
%
%   inpRanges: n-by-2 lower and upper bounds on initial conditions and
%       input ranges, e.g.,
%           inpRanges(i,1) <= x(i) <= inpRanges(i,2)
%       where n = dimension of the initial conditions vector +
%           the dimension of the input signal vector * # of control points
%
%   opt : staliro options object
%
% OUTPUTS:
%   run: a structure array that contains the results of each run of
%       the stochastic optimization algorithm. The structure has the
%       following fields:
%
%           bestRob : The best (min or max) robustness value found
%
%           bestSample : The sample in the search space that generated
%               the trace with the best robustness value.
%
%           nTests: number of tests performed (this is needed if
%               falsification rather than optimization is performed)
%
%           bestCost: Best cost value. bestCost and bestRob are the
%               same for falsification problems. bestCost and bestRob
%               are different for parameter estimation problems. The
%               best robustness found is always stored in bestRob.
%
%           paramVal: Best parameter value. This is used only in
%               parameter querry problems. This is valid if only if
%               bestRob is negative.
%
%           falsified: Indicates whether a falsification occured. This
%               is used if a stochastic optimization algorithm does not
%               return the minimum robustness value found.
%
%           time: The total running time of each run. This value is set by
%               the calling function.
%
%   history: array of structures containing the following fields
%
%       rob: all the robustness values computed for each test
%
%       samples: all the samples generated for each test
%
%       cost: all the cost function values computed for each test.
%           This is the same with robustness values only in the case
%           of falsification.
%
% See also: staliro, staliro_options, UR_Taliro_parameters

% (C) 2010, Sriram Sankaranarayanan, University of Colorado
% (C) 2010, Georgios Fainekos, Arizona State University
params = opt.optim_params;
nSamples = params.n_tests;
StopCond = opt.falsification;

[nInputs, ~] = size(inpRanges); 

% Create sampling grid across solution space
n0_grid = ceil(100*nInputs^(7/4)); 
design_grid = lhsdesign(n0_grid, nInputs, 'criterion', 'maximin').*(inpRanges(:,2) - inpRanges(:,1))' + inpRanges(:,1)';

% Initialize outputs
run = struct('bestRob',[],'bestSample',[],'nTests',[],'bestCost',[],'paramVal',[],'falsified',[],'time',[]);
history = struct('rob',[],'samples',[],'cost',[]);

% %initialize curSample vector
% curSample = repmat({0}, 1, opt.n_workers);

% get polarity and set the fcn_cmp
if isequal(opt.parameterEstimation,1)
    if isequal(opt.optimization,'min')
        fcn_cmp = @le;
        minmax = @min;
    elseif isequal(opt.optimization,'max')
        fcn_cmp = @ge;
        minmax = @max;
    end
else
    fcn_cmp = @le;
    minmax = @min;
end

% if rem(nSamples/opt.n_workers,1) ~= 0
%     error('The number of tests (opt.ur_params.n_tests) should be divisible by the number of workers.')
% end

% Start BO Here 
crowded_EI_flag = params.crowded_EI_flag; %if equal to 1 then use crowded EI, if 0 use EI
n_0 = 5*nInputs;   %number of initial design points
B = 1;     %number of replications per design point, set to 1 for deterministic problem
B_n0_setting=1;

%Crowded EI level set threshold
alpha_lvl_set = params.crowding_threshold; %i.e. EIs within 5% of maxEI

% Instantiate and scale initial design
sim_count = 0;
x_0 = lhsdesign(n_0,nInputs,'criterion','maximin')';
x_0 = x_0.*(inpRanges(:,2) - inpRanges(:,1)) + inpRanges(:,1);


%take first samples, check falsification
for i = 1:n_0  
    curSample{i} = x_0(:,i);
    curVal{i} = Compute_Robustness(curSample{i});
    sim_count = sim_count + 1;
    
    %instantiate storage/history if first sample
    if nargout>1 && i == 1
        if isa(curVal{1},'hydis')
            history.cost = hydis(zeros(nSamples,1));
            history.rob = hydis(zeros(nSamples,1));
        else
            history.cost = zeros(nSamples,1);
            history.rob = zeros(nSamples,1);   
        end
        history.samples = zeros(nSamples,nInputs);
    end
    
    %store as necessary
    if nargout>1
        if isa(curVal{i},'hydis')
            history.cost(i) = hydisc2m(curVal{i})';
            history.rob(i) = hydisc2m(curVal{i})';
        else
            history.cost(i) = curVal{i}';
            history.rob(i) = curVal{i}'; 
        end
        history.samples(i,:) = curSample{i}'; 
    end
    
    %find and store the best value seen so far
    if isa(curVal{1},'hydis')
        [minmax_val, minmax_idx] = minmax(hydisc2m(curVal));
    else
        [minmax_val, minmax_idx] = minmax(cell2mat(curVal));
    end
    bestCost = minmax_val;
    run.bestCost = minmax_val;
    run.bestSample = curSample{minmax_idx};
    run.bestRob = minmax_val;
    run.falsified = minmax_val<=0;
    run.nTests = sim_count;
    
    %check if best value is falsifying, if so, exit as necessary
    if (fcn_cmp(minmax_val,0) && StopCond)
        if nargout>1
            if isa(minmax_val,'hydis')
                history.cost(i+1:end) = hydis([],[]);
                history.rob(i+1:end) = hydis([],[]);
            else
                history.cost(i+1:end) = [];
                history.rob(i+1:end) = [];
            end
            history.samples(i+1:end,:) = [];
        end
        disp(' BayesOpt_Taliro: FALSIFIED by initializing samples!');
        return;
    end
end

%set up for surrogate modeling
xTrain = cell2mat(curSample)';
yTrain = cell2mat(curVal)';
all_x = xTrain;
all_y = yTrain;

clear curSample curVal;

while(sim_count < nSamples)
    %Fit Gaussian Process Meta Model
    GPmod = OK_Rmodel_kd_nugget(xTrain, yTrain, 0, 2);
     
    %Evaluate expected improvement acquisition function over grid
    [pred, mse] = OK_Rpredict(GPmod,design_grid,0,yTrain);    
    EI_values = EIcalc_kd_data(design_grid, xTrain, pred, mse, yTrain);
    
    EI_range = max(EI_values) - min(EI_values);
    relative_top_percent_threshold = EI_range*alpha_lvl_set;
    
   [maxEI, Index] = max(EI_values);
   if crowded_EI_flag == 1 %use crowded EI
       all_maxEI_pointer = 1;
       clear all_maxEI_locations;
       clear crowding_distance;

       %Find all locations in top EI alpha level se %%
       for i=1:size(EI_values,1)
           if EI_values(i) >= (maxEI - relative_top_percent_threshold)
               all_maxEI_locations(all_maxEI_pointer,:)= design_grid(i,:);
               maxEI_value(all_maxEI_pointer,1) = EI_values(i);
               all_maxEI_pointer = all_maxEI_pointer+1;
           end
       end

       %Determine crowding distance at all locations in alpha level set
       for i=1:size(all_maxEI_locations,1)
           for j = 1:nInputs
               DimWise_Crowd(j) = min(abs((ones(size(all_x,1),1).*all_maxEI_locations(i,j))-all_x(:,j)));
           end
           crowding_distance(i) = sum(DimWise_Crowd);
           clear DimWise_Crowd;
       end

       %Find the candidate with "max EI" that has the largest crowding distance
       [max_crowd, Index] = max(crowding_distance);
       x0 = all_maxEI_locations(Index,:); 
   else %use standard EI
       x0 = design_grid(Index,:);
   end
       
   %store acquisition function sample appropriately and sample it
   curSample{1} = x0';
   curVal = Compute_Robustness(curSample);
   f0 = curVal{1};
   sim_count = sim_count + 1;
   
   %store as necessary
    if nargout>1
        if isa(curVal{1},'hydis')
            history.cost(sim_count) = hydisc2m(curVal)';
            history.rob(sim_count) = hydisc2m(curVal)';
        else
            history.cost(sim_count) = curVal{1};
            history.rob(sim_count) = curVal{1}; 
        end
        history.samples(sim_count,:) = curSample{1}'; 
    end
    
    %find and store the best value seen so far
    if isa(curVal{1},'hydis')
        [minmax_val, minmax_idx] = minmax(hydisc2m(curVal));
    else
        [minmax_val, minmax_idx] = minmax(cell2mat(curVal));
    end
    
    if (fcn_cmp(minmax_val,bestCost))
        bestCost = minmax_val;
        run.bestCost = minmax_val;
        run.bestRob = minmax_val;
        run.bestSample = curSample{minmax_idx};
        if opt.dispinfo>0
            if isa(minmax_val,'hydis')
                disp(['Best ==> <',num2str(get(minmax_val,1)),',',num2str(get(minmax_val,2)),'>']);
            else
                disp(['Best ==> ' num2str(minmax_val)]);
            end
        end
    end
    %check if best value is falsifying or if , if so, exit as necessary
    if (fcn_cmp(bestCost,0) && StopCond)
        run.falsified = 1;
        run.nTests = sim_count;
        if nargout>1
            if isa(minmax_val,'hydis')
                history.cost(sim_count+1:end) = hydis([],[]);
                history.rob(sim_count+1:end) = hydis([],[]);
            else
                history.cost(sim_count+1:end) = [];
                history.rob(sim_count+1:end) = [];
            end
            history.samples(sim_count+1:end,:) = [];
        end
        disp(' BayesOpt_Taliro: FALSIFIED!');
        return;
    end
    %check if budget has been exhausted
    if sim_count >= nSamples
        run.nTests = sim_count;
        disp(' BayesOpt_Talro: Samples Exhausted!');
        return;
    end

    all_x = [all_x; x0];
    all_y = [all_y; curVal{1}];
   
    %%% Update the global model training data w/ last centroid sampled %%%
    xTrain = [xTrain; x0];
    yTrain = [yTrain; f0];      
end
disp(' BayesOpt_Talro: Samples Exhausted!');
run.nTests = nSamples; 
end

 
function f = OK_Rlh_kd_nugget(params,k,d,D_X,Y,regr,corr_model, delta)
% Likelihood function for the parameter estimation of modified nugget effect model 
% params - parameters to be estimated, sigma_z and theta
% k,d - size of the input locations
% D_X - distance matrix for the input locations
% Y - observed simulation output values, size [k, 1], k points
% regr - the regression function
% sigma_e - the noise matrix
% corr_model - the correlation model used for the spatial correlation
% corr_model = 0: linear correlation function
% corr_model = 1: exponential correlation function
% corr_model = 2: gaussian correlation function
% corr_model = 3: cubic spline correlation function
% Modified Nugget Effect Kriging toolbox. By YIN Jun, QUAN Ning, NG Szu
% Hui, 2011-2012.

if(min(params(1:d)) <= 0.001)
    f = inf;
    return;
end

theta = params(1:d);
% get correlation matrix given theta
R = OK_corr(corr_model,theta,D_X);
% addition of nugget to increase stability of correlation matrix inversion
R = R + delta.*eye(size(R,1),size(R,2));   

% sum of determinisitc correlation matrix and noise matrix
CR  = R;
[U,pd] = chol(CR);
if(pd>0)
%     if pd1==0
%         U=U1;
%     elseif pd2==0
%         U=U2;
%     else
save data;
error('covariance matrix is nearly singular');
%     end
end

%
L = U';
Linv = inv(L);
Sinv = Linv'*Linv;

% the optimal beta given sigma_z and theta
beta = inv(regr'*Sinv*regr)*(regr'*(Sinv*Y)); 
sigma_z = (1/k)*(Y-regr*beta)'*Sinv*(Y-regr*beta);
%Z = L\(Y-regr*beta);

% negative log likelihood function
f = k*log(sigma_z)+log(det(R));
%f = (log(det(L)) + 0.5*Z'*Z + 0.5*k*log(2*pi));
% else
%     save data;
%     error('covariance matrix is nearly singular');
% end

% Calculate the inverse matrix
end



function [ei_0] = EIcalc_kd_data(x_0,x,predictions,variance,y) % EI calculator

[curr_best, curr_best_ind] =  min(y);

b_0 = ones(size(x_0,1),1); %design for constant mean regression?(vector of 1's as b)

%%%% PLOT THE PREDICTED SURFACE %%%%
%[sortedX, sortIndex] = sort(x_0); %sort x_0 in ascending order
%sortedY = y_0(sortIndex); %sort y_0 maintaing original x_0/y_0 index pairs
%figure();
%plot(sortedX, sortedY) %plot of the global prediction over design_grid (x_0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1;
found=0;
while (i<=size(x_0,1) && found==0)
    if x_0(i,:) == x(curr_best_ind,:)
        curr_best = predictions(i);
        found=1;
    else
        i=i+1;
    end
end
counts = size(x_0,1);
ei_0 = zeros(size(x_0,1),1);
variance = sqrt(variance);

for i = 1:counts
    if variance(i)>0
        ei_0(i) = (curr_best-predictions(i)) * normcdf((curr_best-predictions(i))/variance(i),0,1) + variance(i) * normpdf((curr_best-predictions(i))/variance(i),0,1);
    end
    if variance(i)<=0
        ei_0(i)=0;
    end
    
end

end



function R = OK_corr(corr_model,theta,D_X)
% calculate the correlation matrix for the modified nugget effect model
% corr_model - the correlation model used for the spatial correlation
% corr_model = 0: linear correlation function
% corr_model = 1: exponential correlation function
% corr_model = 2: gaussian correlation function
% corr_model = 3: cubic spline correlation function
% theta - the sensitivity parameters
% D_X - the distance matrix for the input locations
% Modified Nugget Effect Kriging toolbox. By YIN Jun, QUAN Ning, NG Szu
% Hui, 2011-2012.

d1 = size(D_X,1);
d2 = size(D_X,2);
d = size(theta,1);

switch corr_model
    case 0%linear correlation function
        R = prod(max(1 - abs(D_X) .* repmat(reshape(theta,[1 1 d]),[d1 d2]), 0));
    case 1%exponential correlation function
        R = exp(sum((-abs(D_X).^corr_model).*repmat(reshape(theta,[1 1 d]),[d1 d2]),3));
    case 2%Gaussian correlation function  
        R = exp(sum((-abs(D_X).^corr_model).*repmat(reshape(theta,[1 1 d]),[d1 d2]),3));
    case 3%Cubic correlation function
        T = repmat(reshape(theta,[1 1 d]),[d1 d2]);
        R = prod(((D_X<=(T./2)).*(1-6*(D_X./T).^2+6*(D_X./T).^3)+((T./2)<D_X & D_X<=T).*(2*(1-D_X./T).^3)),3);
end


end

function M_model = OK_Rmodel_kd_nugget(Xtrain, Ytrain, regr_model, corr_model) %, sigma_e)
% Build a Ordinary Kriging model for estimation
% X - design locations for the simulation inputs, size [k, d], k points with d
% dimensions 
% Y - observed simulation output values, size [k, 1], k points
% regr_model - the underlying regression model for the mean function:
% regr_model = 0: constant mean function;
% regr_model = 1: linear mean function;
% regr_model = 2: quadratic mean function;
% corr_model - the correlation model used for the spatial correlation
% corr_model = 0: linear correlation function
% corr_model = 1: exponential correlation function
% corr_model = 2: gaussian correlation function
% corr_model = 3: cubic spline correlation function

% Example
%       M_model = OK_model(X,Y, 0,1);
% This function uses OK model with a constant mean function and the 
% exponential correlation function to fit the data, (X,Y)
% where X is the design points and Y are the outputs at design points X 

% Modified Nugget Effect Kriging toolbox. By YIN Jun, QUAN Ning, NG Szu
% Hui, 2011-2012.

%k is the number of records in Xtrain, d is the problem dimension
[k,d] = size(Xtrain);

%error checking
if (size(Ytrain,1) ~= k)
    error('Simulation outputs and inputs should have the same number of rows.')
elseif(size(Ytrain,2) ~= 1)
    error('Simulation outputs should only have one value at each location.')
end

if ((regr_model ~= 0)&(regr_model ~= 1)&(regr_model ~= 2)) %#ok<*AND2>
    error('please type of regression model: 0 = constant, 1 = linear, 2 = quadratic.')
end

if ((corr_model ~= 0)&(corr_model ~= 1)&(corr_model ~= 2)&(corr_model ~= 3))
    error('please type of correlation function: 0 = linear, 1 = exponential, 2 = Gaussian, 3 = cubic spline.')
end

% normalize regression locations
min_X = min(Xtrain);  
max_X = max(Xtrain);
Xtrain = (Xtrain - repmat(min_X,k,1)) ./ repmat(max_X-min_X,k,1);

% calculate the distance between all normalized regression locations
tmp = 0;
D_X = zeros(k, k, d);
tempD_X = zeros(k*k,d);
for h = 1:d
    hh=1;
    for i = 1:k
        for l = 1:k 
%D_X is a matrix of distance between locations, one page per dimension
%(makes use of all loops, l=1:k tracks across each matrix row, i=1:k tracks down matrix tow)
            D_X(i,l,h) = (Xtrain(i,h) - Xtrain(l,h));
%tempD_X has vector per dimension of serialized distances between locations
%(doesn't make use of l=1:k or i=1:k, instead usses hh to track down vector)
            tempD_X(hh,h) = (Xtrain(i,h) - Xtrain(l,h)); 
            hh=hh+1;
        end
    end
end

%call the regression model
regr = OK_regr(Xtrain,regr_model);

% initialize parameters, least square estimation for beta and sample
% variance for sigma_z, 
beta_0 = (regr'*regr)\(regr'*Ytrain); %basic mean function when regr = 1 (constant)
sigma_z0 = var(Ytrain-regr*beta_0);   %variance of residuals (res = observ-mean)
theta_0 = zeros(d,1);                 %initialize hyperparameters, 1 per dimension

if (corr_model == 0 || corr_model ==3)  
    theta_0(1:d) = 0.5;               %if using linear or cubic correlation, initalize all theta to 0.5
else                                      %o.w. theta is: log(2)/d * (mean absolute value of all distances)^(-2 or -1) 
    theta_0(1:d) = (log(2)/d)*(mean(abs(tempD_X)).^(-corr_model));
end

%initialize correlation matrix
% calculate the correlation matrix R based on the initial theta selected and correlation model selected
R = OK_corr(corr_model,theta_0,D_X);           %R=correlation matrix\

% initialize the lower bounds 
a=13;
delta_lb = max((max(eig(R))*(cond(R)-exp(a)))/(cond(R)*(exp(a)-1)),0);

lob_sigma_z = 0.00001*sigma_z0;    
lob_theta = 0.001*ones(d,1); 

lob = [lob_theta];
% maximize profile log-likelihood function (-"logPL")
% subject to lower bounds on theta
myopt = optimset('Display','off','MaxFunEvals',1000000,'MaxIter',500);

for iterations = 1:40
try
%estimate hyperparameters, x represents the paramters [sigma_z; theta1; theta2..]
parms = fmincon(@(x) OK_Rlh_kd_nugget(x,k,d,D_X,Ytrain,regr,corr_model,delta_lb),[theta_0],[],[],[],[],lob,[],[],myopt); 
break
catch
    lob_sigma_z = 0.00001;    
    lob_theta = 0.0001*ones(d,1); 
    lob = [lob_sigma_z;lob_theta];
    
    sigma_z0 = rand();
    theta_0 = rand(d,1);
    parms = fmincon(@(x) OK_Rlh_kd_nugget(x,k,d,D_X,Ytrain,regr,corr_model,delta_lb),[sigma_z0;theta_0],[],[],[],[],lob,[],[],myopt); 
end
end
% record MLEs for theta 
theta = parms(1:length(parms));

% calculate the correlation matrix R based on the theta and correlation model selected
R = OK_corr(corr_model,theta,D_X);           %R=correlation matrix\

CR = (R+delta_lb.*eye(size(R,1),size(R,2)));
[U0,pd0] = chol(CR);
%%[U1,pd1] = chol(sigma_z*R);
%%[U2,pd2] = chol(sigma_e); 
if(pd0>0)

    save data;
 %   error('covariance matrix is nearly singular');
    %end
end
CR=U0;
%cholesky decomposition
L = U0';
D_L = U0';
L_inv = inv(L);
R_inv = L_inv'*L_inv;
beta = inv(regr'*R_inv*regr)*(regr'*(R_inv*Ytrain));
beta_v = inv(regr'*R_inv*regr)*(regr'*(R_inv));
sigma_z = (1/k)*(Ytrain-regr*beta)'*R_inv*(Ytrain-regr*beta);

% output MLEs and other things useful in prediction
M_model.sigma_z =  sigma_z;
M_model.min_X = min_X;
M_model.max_X = max_X;
M_model.regr =  regr;
M_model.beta = beta;
M_model.beta_v = beta_v;
M_model.theta = theta;
M_model.X = Xtrain;
M_model.corr = corr_model;
M_model.L = L;
M_model.D_L = D_L;
M_model.Z = L\(Ytrain-regr*beta);
M_model.Z_v = L\(eye(max(size(Ytrain)))-regr*beta_v);
M_model.Z_m = inv(L);
M_model.DZ_m = inv(D_L);
M_model.Rinv = R_inv;
M_model.nugget = delta_lb;
end



function regr = OK_regr(X,regr_model)
% Call the regression function for the MNEK model
% X - design locations for the simulation inputs, size [k, d], k points with d
% dimensions 
% regr_model - the underlying regression model for the mean function:
% regr_model = 0: constant mean function;
% regr_model = 1: linear mean function;
% regr_model = 2: quadratic mean function;
% Modified Nugget Effect Kriging toolbox. By YIN Jun, QUAN Ning, NG Szu
% Hui, 2011-2012.

%call regression function for the MNEK model
[length,dim] = size(X);
switch regr_model
    case 0
        regr = ones(length,1);
    case 1
        regr = [ones(length,1),X];
    case 2
        mid = (dim+1)*(dim+2)/2;
        regr = [ones(length,1),X,zeros(length,mid-dim-1)];
        j = dim+1;
        q = dim;
        for i = 1:dim
            regr(:,j+(1:q)) = repmat(X(:,i),1,q).*X(:,i:n)
        end
   
end
end



function [f,mse] = OK_Rpredict(model,X_pred,regr_model,Y)
% Build the modified nugget effect predictor based on the model given
% model - modified nugget effect kriging model, given by MNEK_model function 
% X_pred - locations to be predicted 
% regr_model - the underlying regression model for the mean function:
% regr_model = 0: constant mean function;
% regr_model = 1: linear mean function;
% regr_model = 2: quadratic mean function;

% Exmaple
%      M_predict  = MNEK_predict(model,X0,0);
% Using the parameter estimates of the MNEK model obtained from MNEK_model.m,
% the function predicts the response values at prediction points X0 with 
% a constant mean function

% Modified Nugget Effect Kriging toolbox. By YIN Jun, QUAN Ning, NG Szu
% Hui, 2011-2012.

% Obtain model parameters from MNEK_model
X = model.X;
min_X = model.min_X;
max_X = model.max_X;
[k,d] = size(X);
theta = model.theta;
beta = model.beta;
Z = model.Z;
L = model.L;
Rinv = model.Rinv;
sigma_z = model.sigma_z;
corr_model = model.corr;
F = ones(k,1);

% get the size of the locations to be predicted
K = size(X_pred,1);
% get the regression model for the locations to be predicted
regr_pred = OK_regr(X_pred,regr_model);

% normalize distance matrix for prediction points and training points
X_pred = (X_pred - repmat(min_X,K,1)) ./ repmat(max_X-min_X,K,1);

distXpred = zeros(k,K,d);
for h = 1:d 
    for i = 1:k
        for j = 1:K
            distXpred(i,j,h) = (X_pred(j,h)-X(i,h));
        end
    end
end

R_pred = OK_corr(2,theta,distXpred);

% calculate prediction responses and MSE at prediction points 
mse = NaN(K,1);
f = regr_pred*beta + R_pred'*(Rinv*(Y-ones(k,1)*beta));
for r = 1:K
    mse(r,1) = sigma_z * (1 - R_pred(:,r)'*Rinv*R_pred(:,r) + (1-F'*Rinv*R_pred(:,r))'*inv(F'*Rinv*F)*(1-F'*Rinv*R_pred(:,r)));
end
end





