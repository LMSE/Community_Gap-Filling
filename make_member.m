function [community_member] = make_member(n,DB,model)
% This function creates a compartment for the given microorganism in 
% the given medium, that contains a reduced database and can exchange 
% metabolites with the common exchange compartment.
%
% INPUTS
%   n:                 identifier for the microorganism (1 or 2)
%   DB:                database with biochemical reactions in cobra model structure
%   model:             genome-scale metabolic model of the microorganism in cobra model structure
%
% OUTPUTS
%   community_member:  a cobra model structure for the microorganism compartment in the community
%%
DB0 = DB;
model0 = model;

%% Modify the exchange reactions in the database and in the model

% Make the exchange reactions to the database irreversible.
iEX = findExcRxns(DB);
DB = changeRxnBounds(DB,DB.rxns(iEX),0,'l');
DB = changeRxnBounds(DB,DB.rxns(iEX),1000,'u');

% Add blocked exchange reactions for all the exchange metabolites of the model
% that don't participate in an exchange reaction.
iEx = findExcRxns(model);
metsEx = findMetsFromRxns(model,model.rxns(iEx));
ie = (~cellfun('isempty',regexp(model.mets,'.*\[e\]$','match')));
metse = model.mets(ie);
[~,i] = setdiff(metse,metsEx);
lb = zeros(length(i),1);
ub = zeros(length(i),1);
model = addExchangeRxn(model,metse(i),lb,ub);
iEX = findExcRxns(model);
model.rxns(iEX) = regexprep(model.rxns(iEX),'\[(?<compartment>[a-z])\]$','($<compartment>)');

%% Make the organism compartment

% Add or replace in the database the reactions of the model.
for i = 1:length(model.rxns)
    rxn = char(model.rxns(i));
    [row,~,value] = find(model.S(:,i));
    metlist = model.mets(row);
    DB = addReaction(DB,rxn,metlist,value,true,model.lb(i),model.ub(i));
    [~,i_old] = ismember(rxn,DB.rxns);
    DB = moveRxn(DB,i_old,i);
end

% Remove redundant reactions.
community_member = checkDuplicateRxn(DB);

% Make sure that the biomass production is the objective function of the organism compartment.
[~,ibio] = ismember('BIOMASS',community_member.rxns);
community_member = changeObjective(community_member,community_member.rxns(ibio),1);
sol = optimizeCbModel(community_member);

% Rename the metabolites in the organism compartment.
tag = strcat('_',num2str(n));
community_member.mets = strcat(community_member.mets,tag);

% Find and rename the reactions of the database and the model in the organism compartment.
[~,index] = ismember(model.rxns,community_member.rxns);
tagm = strcat('_m',num2str(n));
community_member.rxns(nonzeros(index)) = strcat(community_member.rxns(nonzeros(index)),tagm);
mask = ones(length(community_member.rxns),1);
mask(nonzeros(index)) = 0;
tagd = strcat('_d',num2str(n));
community_member.rxns(logical(mask)) = strcat(community_member.rxns(logical(mask)),tagd);

%% Perform FVA in the organism compartment

% Perform FVA for all the reactions that don't belong to the model.
indices = (cellfun('isempty',regexp(community_member.rxns,'.*_m[0-9]$','match')));
rxns = community_member.rxns(indices);
setWorkerCount(8);
[minFlux,maxFlux] = fastFVA(community_member,0,'max','ibm_cplex',rxns);

% Set new reaction bounds and eliminate the blocked reactions of the database.
community_member = changeRxnBounds(community_member,rxns,minFlux,'l');
community_member = changeRxnBounds(community_member,rxns,maxFlux,'u');
indices_remove = arrayfun(@(x,y) x == 0 && y == 0, minFlux, maxFlux);
community_member = removeRxns(community_member,rxns(indices_remove));

%% Modify the exchange reactions in the microorganism compartment

% Make the exchange reactions between the organism compartment
% and the common exchange compartment.
iEX = findExcRxns(community_member);
rxnsEX = community_member.rxns(iEX);
iEX = full(sum(community_member.S(:,iEX)) == -1)';
rxnsEX = rxnsEX(iEX);
for i = 1:length(rxnsEX)
    metsextra(i,1) = regexprep(regexprep(findMetsFromRxns(community_member,rxnsEX(i)),'_[0-9]$',''),'\[[a-z]\]$','[e]');
end

community_member = addMultipleMetabolites(community_member,unique(metsextra));
[~,i] = ismember(rxnsEX,community_member.rxns);
[~,is] = ismember(metsextra,community_member.mets);
for indexaki = 1:length(i)
    community_member.S(is(indexaki),i(indexaki)) = 1;
end

%% Save outputs
%save(strcat('community_member',tag,'.mat'),'community_member');


end