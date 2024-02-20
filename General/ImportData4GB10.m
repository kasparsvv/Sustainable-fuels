function [Data]=ImportData4GB10(filename,cO)
% [Data]=ImportData4GB10(filename,cO)
%  WARNING: sometimes datafiles use comma as decimal separator. Than the
%           function FAILS. We have a fast matlab tool to convert these files
%           batchwise (see ConvertFiles.m)
%  Input 
%  [mandatory]  filename = full name to file to be read
%  [optional]   cO = cell array containing column order {'time','Encoder','Sensor'} but in
%           the order used in your data files. Above order is the default
%  Reads datafile and stores it in a struct.
%           Data.t          = time array
%           Data.pulse      = pulse array
%           Data.Volt       = sensor voltage array
%           Data.RevEnds    = Indices to revolution end point.
%           Data.NRevs      = Number of complete revolutions
% Note: function trims measurement to full cycles (so even number of revs).
%
delim='\t';
D=importdata(filename,delim, 0);
%     whos D
if (nargin < 2)
    iT=1;iEncoder=2;iSensor=3;
else
    ii = myfind(cO,{'time','Encoder','Sensor'});
    iT          =ii(1);
    iEncoder    = ii(2);
    iSensor     = ii(3);
end
t=D(:,iT);p=D(:,iEncoder);V=D(:,iSensor);
endloc = length(t);
% Find double tooth
[apks,alocs]    = findpeaks(p);
id              = find(apks > 0);
pks             = apks(id);
locs            = alocs(id);
tlocs           = t(locs);
dtlocs          = diff(tlocs);  % compute dt between maxima. Reason, double tooth has bigger dt!
[dlocs]         = find(dtlocs > max(dtlocs)*0.5); % Pick out double tooth
% Take full cycles (is two revolutions! this is a four-stroke engine)
CycleIndices    = locs(dlocs);
NRevs           = length(CycleIndices);
Vc              = V(CycleIndices);
if mod(NRevs,2)==0 % To make sure that only full cycles are taken.
    fprintf('%20s\n %3i revs (will be patched\n',filename,NRevs);
    if (Vc(1)<Vc(2))
        fprintf('IF condition:  %9.3f < %9.3f apparently!\n',Vc(1),Vc(2));
        CycleIndices=CycleIndices(1:end-1);
    else
        fprintf('ELSE condition:%9.3f > %9.3f apparently!\n',Vc(1),Vc(2));
        CycleIndices=CycleIndices(2:end);
    end
else
    fprintf('%20s\n %3i revs (will NOT be patched)\n',filename,NRevs);
    if (Vc(1)<Vc(2))
        fprintf('IF condition:  %9.3f < %9.3f apparently!\n',Vc(1),Vc(2));
        CycleIndices=CycleIndices(1:end-1);
    else
        fprintf('ELSE condition:%9.3f > %9.3f apparently!\n',Vc(1),Vc(2));
        CycleIndices=CycleIndices(2:end);
    end
end
%     FullCycles   = [CycleIndices(1)+1:CycleIndices(end)];
FullCycles   = [CycleIndices(1)+1:endloc];
% Store it in the struct
Data.t       = t(FullCycles)-t(FullCycles(1));  % Start at t=0. Shift it.
Data.pulse   = p(FullCycles);
Data.Volt    = V(FullCycles);
Data.RevEnds = CycleIndices(2:end)-FullCycles(1);
Data.NRevs   = length(Data.RevEnds);
end

