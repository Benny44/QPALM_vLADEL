classdef ladel < handle
    %Matlab interface for the LADEL C library
    %
    %ladel methods
    %
    %
    %
    %
    
    properties (SetAccess = private, Hidden = true)
        ncol;
    end
    
    methods
        
        %% Constructor - Create a new solver instance
        function this = ladel(ncol)
            ladel_mex('init', ncol);
        end
        
        function delete(~)
            ladel_mex('delete');
        end
        
        function varargout = factorize(~, M, varargin)
            M = triu(M);
            if nargin == 3
                ordering = varargin{1};
                ladel_mex('factorize', M, ordering);
            else
                ladel_mex('factorize', M);
            end
            
            if nargout > 0
                if nargout == 2
                    [L, D] = ladel_mex('return');
                    varargout{1} = L;
                    varargout{2} = D;
                elseif nargout == 3
                    [L, D, p] = ladel_mex('return');
                    varargout{1} = L;
                    varargout{2} = D;
                    varargout{3} = p;
                else
                    error('Wrong number of output arguments for factorize');
                end
            end
        end
        
        function y = dense_solve(~, x)
            y = ladel_mex('solve', x);
        end
        
        function varargout = factorize_advanced(~, M, Mbasis, varargin)
            M = triu(M);
            Mbasis = triu(Mbasis);
            if nargin == 4
                ordering = varargin{1};
                ladel_mex('factorize_advanced', M, Mbasis, ordering);
            else
                ladel_mex('factorize_advanced', M, Mbasis);
            end
            
            if nargout > 0
                if nargout == 2
                    [L, D] = ladel_mex('return');
                    varargout{1} = L;
                    varargout{2} = D;
                elseif nargout == 3
                    [L, D, p] = ladel_mex('return');
                    varargout{1} = L;
                    varargout{2} = D;
                    varargout{3} = p;
                else
                    error('Wrong number of output arguments for factorize_advanced');
                end
            end
                    
        end
        
        function varargout = factorize_with_prior_basis(~, M)
            M = triu(M);         
            ladel_mex('factorize_with_prior_basis', M);
            if nargout > 0
                if nargout == 2
                    [L, D] = ladel_mex('return');
                    varargout{1} = L;
                    varargout{2} = D;
                elseif nargout == 3
                    [L, D, p] = ladel_mex('return');
                    varargout{1} = L;
                    varargout{2} = D;
                    varargout{3} = p;
                else
                    error('Wrong number of output arguments for factorize_with_prior_basis');
                end
            end 
        end
        
        function varargout = row_mod(~, row, varargin)
            
            if nargin ~= 2 && nargin ~= 4
                error('Wrong number of input arguments for row_mod.\n Use .row_mod(row) to delete a row (with index row) or .row_mod(row, w, diag_elem) to add a row w and with on the diagonal diag_elem.\n');
            end
            if nargin == 2
                ladel_mex('rowmod', int64(row), length(row));
            else
                w = varargin{1};
                diag_elem = varargin{2};
                ladel_mex('rowmod', int64(row), length(row), w, diag_elem);
            end
            
            if nargout > 0
                if nargout == 2
                    [L, D] = ladel_mex('return');
                    varargout{1} = L;
                    varargout{2} = D;
                elseif nargout == 3
                    [L, D, p] = ladel_mex('return');
                    varargout{1} = L;
                    varargout{2} = D;
                    varargout{3} = p;
                else
                    error('Wrong number of output arguments for factorized_advanced');
                end
            end
        end
    end
    
end

