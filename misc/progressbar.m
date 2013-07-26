% Display ASCII waiting bar
%
% Inputs
%  max_index : maximum index
%  -- or --
%  pb        : progress bar structure
%  i         : current
%  title     : text shown after the remaining time
%
% Outputs
%  pb        : progress bar structure
%
% Example
%  pb = progressbar(max_index);
%  for index=1:max_index;
%    [...]
%    pb = progressbar(pb, index);
%  end
%
% Copyright (c) 2008 Ircam-CNRS UMR9912-STMS
% 
% License
%  This file is under the LGPL license,  you can
%  redistribute it and/or modify it under the terms of the GNU Lesser General 
%  Public License as published by the Free Software Foundation, either version 3 
%  of the License, or (at your option) any later version. This file is
%  distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
%  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
%  PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
%  details.
%
% This function is part of the Covarep project: http://covarep.github.io/covarep
%
% Author
%  Gilles Degottex <gilles.degottex@ircam.fr>
%

function pb = progressbar(pbin, i, text)

	if nargin==1
		pb.max_n = pbin;
		pb.history = [];
		pb.i_old = 0;
        pb.finished = false;
        pb.nbp = 0;
		tic;
		return;
	else
		pb = pbin;
	end

	if toc<1/16 && i<pb.max_n; return; end

%  	if nargin~=2
%  		return;
%  	end

    if ~pb.finished

	    pb.history = [pb.history, toc];
	    pb.iteration_time = mean(pb.history);

	    pos = i / pb.max_n;
	    str='|=';
	    for c = 1:floor(pos*50)
		    str = [str, '='];
	    end
	    for c = floor(pos*50):50-1
		    str = [str, ' '];
	    end

	    di = i - pb.i_old;
	    nrest = (pb.max_n-i)/di;
	    minutes = floor(sum(pb.history)/60);
	    seconds = rem(sum(pb.history), 60);
	    str_times = ['  elapsed: ', num2str(minutes), ':', num2str(floor(seconds))];
	    if nrest>0
		    rem_minutes = floor(pb.iteration_time*nrest/60);
		    rem_seconds = rem(pb.iteration_time*nrest, 60);
		    str_times = [str_times, ' remaining: ', num2str(rem_minutes), ':', num2str(round(rem_seconds))];
	    else
            str_times = [str_times, '               '];
        end

		str_full = [num2str(floor(pos * 100)),'%% ', str_times];
		if 1
			str_full = [repmat('\b',1,pb.nbp),str_full];
		else
			str_full = ['\r',str_full];
		end

        if nargin>2
        	str_full = [str_full, ': ', text];
        end
       	str_full = [str_full, '       '];
        nbp = fprintf(str_full);
        pb.nbp = nbp - pb.nbp;

	    pb.i_old = i;
	    tic;

        if i>=pb.max_n;
            pb.finished = true;
            disp(' ');
        end
    end

return
