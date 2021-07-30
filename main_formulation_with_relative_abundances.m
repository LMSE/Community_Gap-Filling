%% Initialize cobratoolbox and choose solver (CobraToolbox v.3, CPLEX 12.8, MATLAB R2017b)
initCobraToolbox(0);
changeCobraSolver('ibm_cplex','all');

%% Load input models and files
DB = readCbModel('BiGG.mat');  % database
model1 = readCbModel('model1.mat');  % model of organism 1
model2 = readCbModel('model2.mat');  % model of organism 2
community_media = textscan(fopen('community_media.txt'),'%s %s %d %d');  % media for the community
f1 = 0.5;  % relative abundance of organism 1 in the community
f2 = 0.5;  % relative abundance of organism 2 in the community

%% Create the community model

% Make the organism compartments and combine them.
member1 = make_member(1,DB,model1);
member2 = make_member(2,DB,model2);
delete(gcp('nocreate'))
community = mergeTwoModels(member1,member2);

% Make sure that the biomass production is the objective function.
[~,ibio1] = ismember('BIOMASS_m1',community.rxns);
[~,ibio2] = ismember('BIOMASS_m2',community.rxns);
community.c(ibio1) = 1;
community.c(ibio2) = 1;
community = changeRxnBounds(community,community.rxns(ibio1),0.9,'l');
community = changeRxnBounds(community,community.rxns(ibio2),0.09,'l');

%% Make the common exchange compartment
EX_name = community_media{1};
EX_lb = double(community_media{3});
EX_ub = double(community_media{4});

% Add irreversible exchange reactions for all the exchange metabolites of the community model.
ie = (~cellfun('isempty',regexp(community.mets,'.*\[e\]$','match')));
metse = community.mets(ie);
lb(1:length(metse),1) = 0;
ub(1:length(metse),1) = 1000;
community = addExchangeRxn(community,metse,lb,ub);

% Add the exchange reactions from the media to the community model.
rxnsEx = strrep(EX_name,'(e)','[e]');
community = changeRxnBounds(community,rxnsEx,EX_lb,'l');
community = changeRxnBounds(community,rxnsEx,EX_ub,'u');

% Rename the exchange reactions of the community model.
[~,iEx] = ismember(rxnsEx,community.rxns);
community.rxns(iEx) = strrep(rxnsEx,'[e]','(e)_e');
iEX = (~cellfun('isempty',regexp(community.rxns,'.*\[e\]$','match')));
community.rxns(iEX) = strrep(community.rxns(iEX),'[e]','(e)_d3');

%% Check the community model for discrepancies

% Find in the community model the reactions of each group
% and check if each reaction exists only one time.
m1match = regexp(community.rxns,'.*_m1$','match');
d1match = regexp(community.rxns,'.*_d1$','match');
m2match = regexp(community.rxns,'.*_m2$','match');
d2match = regexp(community.rxns,'.*_d2$','match');
ematch = regexp(community.rxns,'.*_e$','match');
d3match = regexp(community.rxns,'.*_d3$','match');
m1 = (~cellfun('isempty',m1match));
d1 = (~cellfun('isempty',d1match));
m2 = (~cellfun('isempty',m2match));
d2 = (~cellfun('isempty',d2match));
e = (~cellfun('isempty',ematch));
d3 = (~cellfun('isempty',d3match));
test = m1 + d1 + m2 + d2 + e + d3;
if any(test==1) && sum(test)==length(community.rxns)
    fprintf('The reactions are properly organized. Proceed.\n');
else
    error('The reactions are not organized.');
end

%% Save community model
%save('community.mat','community');

%% Solve the community gap-filling problem
[nmets,nrxns] = size(community.S);

% MILP formulation.
% x = [v y] (2*nrxns x 1)
ctype = char(['C' * ones(1,nrxns) 'B' * ones(1,nrxns)]);
f = [sparse(nrxns,1); d1+d2+d3];

lb_all = spdiags(f1*m1+f1*d1+f2*m2+f2*d2+e+d3,0,sparse(nrxns,nrxns))*community.lb;
ub_all = spdiags(f1*m1+f1*d1+f2*m2+f2*d2+e+d3,0,sparse(nrxns,nrxns))*community.ub;

Aineq = [-spdiags(m1+d1+m2+d2+e+d3,0,sparse(nrxns,nrxns)) spdiags(lb_all,0,sparse(nrxns,nrxns));
          spdiags(m1+d1+m2+d2+e+d3,0,sparse(nrxns,nrxns)) -spdiags(ub_all,0,sparse(nrxns,nrxns))];
bineq = sparse(2*nrxns,1);
Aeq = [community.S sparse(nmets,nrxns)];
beq = sparse(nmets,1);
lb = [(f1*m1+f1*d1+f2*m2+f2*d2+e+d3).*community.lb; m1+m2+e];
ub = [(f1*m1+f1*d1+f2*m2+f2*d2+e+d3).*community.ub; ones(nrxns,1)];

prob.ctype = ctype;
prob.f = f;
prob.Aineq = Aineq;
prob.bineq = bineq;
prob.Aeq = Aeq;
prob.beq = beq;
prob.lb = lb;
prob.ub = ub;
cplex = Cplex(prob);

% Solve the MILP and find alternative solutions.
CPX_PARAM_SOLNPOOLINTENSITY = 2;
CPX_PARAM_EPINT = 1e-10;
CPX_PARAM_EPRHS = 1e-9;
CPX_PARAM_MIPEMPHASIS = 1;
CPX_PARAM_FPHEUR = 1;
CPX_PARAM_NODESEL = 2;
CPX_PARAM_EPMRK = 0.99;
CPX_PARAM_NUMERICALEMPHASIS = 1;
CPX_PARAM_TILIM = 7200;
solutions = populate(cplex);
save('solutions.mat','solutions');

% Save the fluxes of the reactions for the best solution.
x = solutions.pool.solution(1).x;
v = x(1:nrxns);
y = x(nrxns+1:end);
data = [community.rxns(y>0.5) num2cell(v(y>0.5))];
T = cell2table(data,'VariableNames',{'Reaction','Flux'});
writetable(T,'Solution_Fluxes.xlsx');

% Save the fluxes of the exchange and the added database reactions for the 10 best solutions.
sols = solutions.pool.solution;
objval_i = (1:length(sols))';
for i = 1:length(sols)
    objval(i,1) = sols(i).objval;
    xs(:,i) = sols(i).x;
end
objval_table = sortrows([objval,objval_i],1);
vs = xs(1:nrxns,:);
ys = xs(nrxns+1:end,:);
Vs = vs(:,objval_table(1:10,2));
Ys = ys(:,objval_table(1:10,2));

ex_rxns = [];
d_rxns = [];
for i = 1:10
    vrxns = community.rxns(abs(Vs(:,i))>=1e-9);
    yrxns = community.rxns(Ys(:,i)>0.5);
    rxns = intersect(yrxns,vrxns);
    is = (~cellfun('isempty',regexp(rxns,'^EX_.*\d$','match')));
    ex_rxns = [ex_rxns;rxns(is)];
    is = (~cellfun('isempty',regexp(rxns,'.*_d[0-9]$','match')));
    d_rxns = [d_rxns;rxns(is)];
end

ex_rxns = unique(ex_rxns);
[~,index] = ismember(ex_rxns,community.rxns);
fluxes = Vs(index,:);
data = [ex_rxns num2cell(fluxes)];
T = cell2table(data,'VariableNames',{'ExchangeReaction','Solution1','Solution2','Solution3','Solution4','Solution5','Solution6','Solution7','Solution8','Solution9','Solution10'});
writetable(T,'Solution_ExchangeReaction_Fluxes.xlsx');

d_rxns = unique(d_rxns);
[~,index] = ismember(d_rxns,community.rxns);
fluxes = Vs(index,:);
data = [d_rxns num2cell(fluxes)];
T = cell2table(data,'VariableNames',{'Reaction','Solution1','Solution2','Solution3','Solution4','Solution5','Solution6','Solution7','Solution8','Solution9','Solution10'});
writetable(T,'Solution_AddedDatabaseReaction_Fluxes.xlsx');