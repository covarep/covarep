function [LENG,numb] = beginend(filename)
% This function imports the times of beginning and end of the intervals over
% which the electroglottographic signal is to be analyzed.
% It can handle text files of two types: either a text file containing columns
% of figures (in which case it reads the last two columns), or a file
% produced by copying the Regions List written under the software SoundForge and
% pasting it into a .txt file. Annotations produced in SoundForge contains a heading
% and, for each token, its name, beginning, end, and length according to an
% implicit format.
% If the input is structured in neither of these two different way, the function
% stops. (This function can be expanded and adapted to read other input
% formats.)
% In the output matrix, times are converted into milliseconds.  

% finding out whether there are any nonnumeric characters in the file

% % Previous version: looking for any nonnumeric character; but in fact
% % individual nonnumeric characters do not prevent the correct loading of a
% % file: e.g. the two lines
% 24a	61929	62207
% 24b	65404	65742
% % load properly. Test to find out whether the input file is a SoundForge
% % Regions List file or a matrix: if there is more than one nonnumeric
% % character in the first line: treat as SoundForge Regions List.
numer = 0;
	% fid = fopen(filename, 'rt');
	% while and ( feof(fid) == 0, numer == 0)
	%    tline = '';
	%    numline = '';
	%    tline = fgetl(fid);
	%    for n = 1:length(tline)
	%        % removing nonnumeric characters from the line, retaining spaces
	%        if tline(n) ~= ' '
	%            % special case for i, which if it becomes a number is an imaginary
	%            % number and not one of the 1..9 signs.
	%            if tline(n) ~= 'i'
	%                numline = [numline str2num(tline(n))];
	%            end
	%        else
	%            numline = [numline,' '];
	%        end
	%    end
	%    if length(numline) ~= length(tline)
	%        % in case there is a nonnumeric character, <numline> and <tline> are
	%        % not identical; in this case: set numer at 1.
	%        % The file cannot be loaded as a matrix.
	%        numer = 1;
	%    end
	% end
    
fid = fopen(filename, 'rt');
tline = ''; numline = '';
tline = fgetl(fid);
fclose(fid);
nbletters = length(nonzeros(isletter(tline)));
if nbletters > 1
   % in case there is more than one nonnumeric character in the first line:
   % set numer at 1; the file will be handled as Regions List.
   % The file cannot be loaded as a matrix.
   % Note: the END OF LINE character is not transformed into a numerical
   % character; the length of numline will therefore be less than that of
   % tline by one character even when there are no letters in the line
   numer = 1;
end
    
	% for n = 1:length(tline)
	%    % removing nonnumeric characters from the line, retaining spaces
	%    if tline(n) ~= ' '
	%        % special case for i, which if it becomes a number is an imaginary
	%        % number and not one of the 1..9 signs.
	%        if tline(n) ~= 'i'
	%            numline = [numline num2str(str2num(tline(n)))];
	%        end
	%    else
	%        numline = [numline,' '];
	%    end
	% end

  
% in case it has been detected that the file is basically made up of figures: 
% it can be loaded as a matrix:
if numer == 0
    numb = [];
    LENG = load(filename);
    [m,n] = size(LENG);
    % if there are more than 2 columns: retrieving the last two, and placing the
    % first in <numb> vector
    if n > 2
        % if the number of columns is higher than 2: selecting the last two columns.
        % This can happen in case the first column contains numbers associated with the
        % items, for instance.
        numb = rot90(LENG(:,1));
        LENG = LENG(:,n - 1:n);
        disp('There are more than two columns in .txt file containing time boundaries.')
        disp('Last two columns selected.')
    elseif n < 2
        error('Input text file must contain two values for each item: beginning and end.')
    end
    
% if nonnumeric characters were found: attempt at reading the file as SoundForge
% Regions list.
elseif numer == 1
    fid = fopen(filename, 'rt');
    linecount = 0;
    % variable containing beginning and end of sample in ms
    LENG = [];
    % variable containing number-label attached by user to item in SoundForge, if any
    numb = [];
    % variable counting excluded lines
    badlines = 0;

    while feof(fid) == 0
       linecount = linecount + 1;
       tline = fgetl(fid);
       % Test that works for files created by SoundForge's "Regions Lists" functions:
       % if the last character of line is a number then this is a line that
       % contains relevant data.
       if length(tline) > 0
           if ~isempty(str2num(tline(length(tline))))
               begpt = tline(length(tline) - 39:length(tline) - 28);
               % calculating beginning of item, in ms, converting it from the
               % <hours:minutes:seconds:milliseconds> format used by SoundForge
               LENG(linecount - (badlines+1),1) = str2num(begpt(1:2)) * 3600000 + ...
                   str2num(begpt(4:5)) * 60000 + str2num(begpt(7:8)) * 1000 + ...
                   str2num(begpt(10:12));
               % calculating end of item, in ms
               endpt = tline(length(tline) - 25:length(tline) - 14);
               LENG(linecount - (badlines+1),2) = str2num(begpt(1:2)) * 3600000 + ...
                   str2num(endpt(4:5)) * 60000 + str2num(endpt(7:8)) * 1000 + ...
                   str2num(endpt(10:12));
           else
               badlines = badlines + 1;
           end
           % retrieving number given to the item in the Regions List, if any
           NUMB = [];
           i = 1;
           while ismember(str2num(tline(i)),0:9)
               NUMB = [NUMB tline(i)];
               i = i + 1; 
           end
           if ~isempty(NUMB)
               numb(linecount - (badlines+1)) = str2num(NUMB);
           end
       end
    end
    fclose(fid);
end