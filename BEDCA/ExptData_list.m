
%
% Place your data entries here. 
%
% ExptDataLoc is location relative to MCMC directory. Can be relative or absolute.
% Give file name and in ExptDats give type of fluorophore
%
% Add more as needed.
%

% 
% This is a treatment name list, first is name used in MCMC files, 2nd entry is search string in your Experiment file names - your KiT file names
% must include it and uniquely identify the files (it is vital not to use search names that are contained within other search terms). 
% Ensure unique. Code does not distinguish between capitals v small letters.
% Note DMSO == untreated so appears twice as two different search terms.
  % Later entries are aka; can be used as runtreatment in EuclDistCorrection_Driver_fn
%
% We used standard concentrations (Noc, Taxol) so allowed names that didnt specify conc.

% Add new treatments to array
TreatmentLib={{'DMSO','DMSO'},{'DMSO','untreated'},{'Taxol15min','tax15min','tax'},{'Nocod2hr','noc2hr','noc'},{'Nocod45min','noc45min'},{'Nocod30min','noc30min'},{'Nocod15min','noc15min'},{'DMSOMG','DMSOMG'},{'3uMnocMG','3uMnocMG'},{'3uMnocRevMG','3uMnocRevMG'},{'1uMTaxolMG','1uMtaxMG','taxMG'},{'1uMTaxolRevMG','1uMtaxRevMG','taxRevMG'},{'CTRrnai','CTRi'},{'CenpCrnai','CenpCi'},{'CenpTrnai','CenpTi'}};

% Directory with data (relative location to this file)
ExptDataLoc= '../BEDCA/ExptData/';

% FILE NAMES IN ABOVE DIRECTORY
% Compilation sets over days (in ExptLabel)
% the exp26/exp28/exp30 pool and exp27/29/30 pool are DMSO, exp28/exp30 pool are taxol)
%FileNames{1}='iM_Pool_exp26_exp28_exp30_NGSka9G3CenpC_SCenpCSka_9G3CenpC.mat'; DMSO,example

% Create an associated ExptDat. Use concatenated markers name.
%ExptDats{1}.name='HeLASka9G3CenpC_m23.mat'; example


% This is groupings by paired markers. .name should match across group
%ExptDats{1}.group=[1 2];ExptDats{2}.group=[1 2];  % Bub1

%
% Treatments
%
%ExptDats{1}.treatment='DMSO'; example

%Fluorophore

%ExptDats{1}.fluorophore='Ab'; example

