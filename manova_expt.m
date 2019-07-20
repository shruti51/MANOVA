clear
load manova_data.mat
cat_ = reshape(repmat(cat,8,1),32,1);
t = table(cat_,data(:,1),data(:,2),data(:,3),'VariableNames',{'soil_type ','yield','water','herbicide'});

%% Problem Statement:

% Problem: Want to determine whether there was a significant difference
% between the types of soils (loam, sandy, salty, clay) based on the yield of the
% crop, amount of water required and amount of herbicide needed. Eight
% fields are chosen for analysis. 

%% Means Calculations
% Total Mean
total_mean = [mean(t.yield), mean(t.water), mean(t.herbicide)];

% Group Mean (corresponding to 4 groups namely loam, sandy, salty, clay)
group_means = zeros(4,3);
group_means_minus_total_mean = zeros(4,3);
data_minus_total_mean = bsxfun(@minus,data,total_mean);
t_minus_means = table(cat_,data_minus_total_mean(:,1),data_minus_total_mean(:,2),data_minus_total_mean(:,3),...
    'VariableNames',{'soil_type ','yield','water','herbicide'});
for group = 1:4
   group_=8*(group-1)+1:8*group;
   group_means(group,:) = [mean(t.yield(group_)), mean(t.water(group_)), mean(t.herbicide(group_))]; 
   group_means_minus_total_mean(group,:)=[mean(t_minus_means.yield(group_)), mean(t_minus_means.water(group_)),...
       mean(t_minus_means.herbicide(group_))];
end

% Total Mean and Group Mean table:
cat__= cat';
cat__{5,1} = 'total';
t_means = table(cat__, [group_means(:,1);total_mean(:,1)],[group_means(:,2);total_mean(:,2)],...
    [group_means(:,3);total_mean(:,3)],'VariableNames',{'soil_type ','yield','water','herbicide'});

%% Various sum of squares parameters calculations
G = 4; % Number of groups.
n_g = 8;  % Number of elements in each group (assumed same in this example).
n = n_g * G; 
D = 3;  % Number of measurements (dimensionality of every mean).

T = zeros(D,D);
H = zeros(D,D);
E = zeros(D,D);
data_minus_group_mean = zeros(size(data));

for g = 1: G
    group_=8*(g-1)+1:8*g;
    
    % Total Cross Products for p and q:
    % Degrees of freedom in T:  n-1 (total data -1) where total data is in
    % one coordinate.
    T = T + data_minus_total_mean(group_,:)'*data_minus_total_mean(group_,:) ;
    
    % Hypothesis sum of squares and cross products:
    % dof in H: g-1. (B)
    H = H + n_g * group_means_minus_total_mean(g,:)'* group_means_minus_total_mean(g,:);
    
    % Error (or residual) sum of squares and cross products:
    % dof in E: n-g (W)
    data_minus_group_mean(group_,:) = bsxfun(@minus, data(group_,:),group_means(g,:));
    E = E +  data_minus_group_mean(group_,:)'*data_minus_group_mean(group_,:);
    
    % Note that dof(T) = dof(H)+dof(E)
end

% Observation T= H+E
assert(norm(T-(H+E))/9 < 1e-5)

%% Create an F-test (in multivariate setting)
% By dividing H by E or HE^{-1}. In manova, we reject the null hypothesis if
% H is 'large' compared to E. 

% 1) Wilk's Lambda: \lambda = |E|/|H+E| = |E|/|T|.
% We reject the null hypothesis when Wilk's Lambda is close to 0.
lambda = det(E)/det(T); 
a = n - g - (D-g+2)/2;
b = sqrt((D^2*(G-1)^2 - 4)/(D^2+(G-1)^2 - 5));
c = (D*(G-1)-2)/2;
df1 = D*(G-1);
df2 = a*b - c;
F = (1-lambda^(1/b))/(lambda^(1/b)) * df2/df1;

% 2) Hotelling-Lawley Trace: T_0^2 = trace(HE^{-1})
% H is large compared to E when T_0^2 is large and we thus reject the null
% hypothesis.
T0_sq = trace(H/E);

% 3) Pillai-Barlett Trace: V = trace(H(H+E)^{-1}=trace(HT^{-1}
% If H is large compared to E then this statistics will be large and we
% will reject the null hypothesis.
V = trace(H/T);

% 4) Roy's Largest root: \theta = largest eigenvalue of HE^{-1}
% Reject the null hypothesis if the statistic is large. Sometimes we use
% \lambda_p / (1+\lambda_p) where \lambda_p is the largest eigenvalue of
% HE^{-1}
lambda_p = eigs(H/E,1);
theta = lambda_p/(1+lambda_p);







