function [velocityMap, optOut] = runParallelSTICS(opt, varargin)
%runParallelSTICS Allows parallel computing of STIC(C)S. 
%
%   Each run can be either a different file using the same parameters or
%   the same file using different parameters, or both. By defining a second
%   stack, one can use this same code to perform either STICS
%   (auto-correlation, default) or STICCS (cross-correlation). 
%   
%   
%   Rodrigo Migueles, 2021 under the supervision of Paul Wiseman and Arnold
%   Hayer at McGill University.

%     if isempty(varargin) == 0
%         for i = 1:length(opt)
%             if opt(i).STICCS
%                 stack1(i) = readFileToStack(opt(i).filePath{1});
%                 stack2(i) = readFileToStack(opt(i).filePath{2});
%                 disp('Starting STICCS analysis');
%             else
%                 stack1{i} = readFileToStack(opt(i).filePath{1});  
%                 disp('Starting STICS analysis...');
%             end
%         end
%     else
    if length(varargin) == 1
        stack1 = varargin{1};
        disp('Starting STICS analysis...');
    elseif length(varargin) == 2
%         disp('Starting STICCS analysis');
%         stack1 = varargin{1};
%         stack2 = varargin{2};       
    end
    
    
    parfor i = 1:length(opt)
        disp(['Run #', num2str(i)])
        if opt(i).STICCS
            [velocityMap{i,1},  ~, ~, ~, optOut{i}] = sticcs_vectormapping(stack1, stack2, opt(i));
%             varFile(i) = load([opt.path 'VelocityMap' opt(i).outputName '.mat'], 'velocityMap');
%             positionsMap{i} = {position_x, position_y, position_t};
        elseif ~opt(i).STICCS
            [velocityMap{i,1},  ~, ~, ~, optOut{i}] = stics_vectormapping(stack1, opt(i));
%             positionsMap{i} = {position_x, position_y, position_t};
        end
    end

end
