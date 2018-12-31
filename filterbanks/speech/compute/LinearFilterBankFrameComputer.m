% Frame computers whose features are derived from linear filter banks
% 
% Description
% Computers based on linear filter banks have a predictable number of
% coefficients and organization. Like the banks, the features with
% lower indices correspond to filters with lower bandwidths.
% `num_coeffs` will be simply `bank.num_filts + int(include_energy)`.
% 
% Properties
% ----------
% bank : [LinearFilterBank]
%     Each filter in the bank corresponds to a coefficient in a
%     frame vector. Can be a LinearFilterBank or something compatible
%     with `pydrobert.speech.alias_factory_subclass_from_arg`
% num_coeffs : [int]
% include_energy : [bool], optional
%     Whether to include a coefficient based on the energy of the
%     signal within the frame. If ``True``, the energy coefficient
%     will be inserted at index 0.
%     (def. false)
% 
%
% Copyright (c) 2018 Department of Computer Science,
%                    University of Toronto, Canada,
%                    Vector Institute, Canada
%
% License
% This file is under the LGPL license,  you can
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
%  Yingxue Wang <yingxue@cs.toronto.edu>
%  Sean Robertson <sdrobert@cs.toronto.edu>
%


classdef LinearFilterBankFrameComputer <  FrameComputer
    properties (Access = private)
        bank %= LinearFilterBank()
        num_coeffs
        include_energy
    end
    
    methods
        function obj = LinearFilterBankFrameComputer(bank, include_energy)
            % Constructor of the LinearFilterBankFrameComputer
            %
            % Inputs
            % bank : [LinearFilterBank]
            % includes_energy : [bool], optional, (def. false)
            % 
            % Outputs
            % obj : LinearFilterBankFrameComputer object
            % 
            % Example
            % >> bank = LinearFilterBank();
            % >> include_energy = false;
            % >> LinearFilterBankFrameComputer(bank, include_energy);
            %
            
            if nargin == 1
                include_energy = false;
            end
            
            % Call superclass constructor before accessing object
            obj = obj@FrameComputer();
            
            obj.bank = bank;
            obj.include_energy = include_energy;
            obj.num_coeffs = obj.bank.num_filts + cast(include_energy, 'int8');
        end
        
        % No other methods needed for linear Filters
    
    end
end