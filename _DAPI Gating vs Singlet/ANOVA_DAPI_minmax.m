% Setting up ANOVA for the minmax results 
DAPI_minmax = readtable("DAPI_minmax_export.csv");


%% 
y = DAPI_minmax.range;

% Use experiment as factor and ligand concenctration as a continuous fct
g1 = categorical(DAPI_minmax.Experiment); % ANOVA requires categorical time, rather than datetime
g2 = log10(DAPI_minmax.Ligand_conc);


p = anovan(y, {g1,g2}, "continuous", 2); % "continuous", 2 to treat g2 as continuous. 

