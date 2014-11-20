function [result, why] = structeq(struct1, struct2, funh2string, ignorenan)
% STRUCTEQ performs an equality comparison between two structures by
% recursively comparing the elements of the struct array, their fields and
% subfields. This function requires companion function CELLEQ to compare
% two cell arrays.
%
% USAGE:
%
% structeq(struct1, struct2)
%       Performs a comparison and returns true if all the subfields and
%       properties of the structures are identical. It will fail if
%       subfields include function handles or other objects which don't
%       have a defined eq method. 
%
% [iseq, info] = structeq(struct1, struct2)
%       This syntax returns a logical iseq and a second output info which
%       is a structure that contains a field "Reason" which gives you a
%       text stack of why the difference occurred as well as a field
%       "Where" which contains the indices and subfields of the structure
%       where the comparison failed. If iseq is true, info contains empty
%       strings in its fields.
%
% [...] = structeq(struct1, struct2, funh2string, ignorenan)
%       Illustrates an alternate syntax for the function with additional
%       input arguments. See the help for CELLEQ for more information on the
%       meaning of the arguments
%
% METHOD:
% 1. Compare sizes of struct arrays
% 2. Compare numbers of fields
% 3. Compare field names of the arrays
% 4. For every element of the struct arrays, convert the field values into
% a cell array and do a cell array comparison recursively (this can result
% in multiple recursive calls to CELLEQ and STRUCTEQ)
%
% EXAMPLE:
% % Compare two handle graphics hierarchies
% figure;
% g = surf(peaks(50));
% rotate3d
% hg1 = handle2struct(gcf);
% set(g,'XDataMode', 'manual');
% hg2 = handle2struct(gcf);
% 
% structeq(hg1, hg2)
% [iseq, info] = structeq(hg1, hg2)
% [iseq, info] = structeq(hg1, hg2, true)


if nargin < 3
    funh2string = false;
end
if nargin < 4
    ignorenan = false;
end

why = struct('Reason','','Where','');

if any(size(struct1) ~= size(struct2))
    result = false;
    why = struct('Reason','Sizes are different',Where,'');
    return
end

fields1 = fieldnames(struct1);
fields2 = fieldnames(struct2);

% Check field lengths
if length(fields1) ~= length(fields2)
    result = false;
    why = struct('Reason','Number of fields are different','Where','');
    return
end

% Check field names
result = celleq(fields1,fields2);
result = all(result);
if ~result
    why = struct('Reason','Field names are different','Where','');
    return
end

for i = 1:numel(struct1)
    props1 = struct2cell(struct1(i));
    props2 = struct2cell(struct2(i));
    [result, subwhy] = celleq(props1,props2,funh2string,ignorenan);
    result = all(result);
    if ~result
        
        [fieldidx, subwhy.Where] = strtok(subwhy.Where, '}');
        fieldidx = str2double(fieldidx(2:end));
        %str2double(regexp(subwhy.Where,'{([0-9]+)}','tokens','once'));
        where = sprintf('(%d).%s%s',i,fields1{fieldidx},subwhy.Where(2:end));
        why = struct('Reason',sprintf('Properties are different <- %s',subwhy.Reason),'Where',where);
        return
    end
end
