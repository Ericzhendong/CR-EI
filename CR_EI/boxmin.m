% function  [t, f, fit, perf] = boxmin(t0, lo, up, par)
function  [t, f, perf] = boxmin(X0,fname, lo, up, varargin)
%BOXMIN  Minimize with positive box constraints
   global Z;
   global varInputs;
   Z  = X0;
   t0 = unwrap(X0);
   varInputs = varargin;

    % Initialize
    [t, f, itpar] = start(fname,t0, lo, up);
    if  ~isinf(f)  % Iterate
        p = length(t);
        if  p <= 2,  kmax = 2; else,  kmax = min(p,4); end
        
        for  k = 1 : kmax
            
            th = t;        
            [t, f, itpar] = explore(fname,t, f, itpar);
            [t, f, itpar] = move(fname,th, t, f, itpar);
        end
    end
    perf = struct('nv',itpar.nv, 'perf',itpar.perf(:,1:itpar.nv));
    t    = rewrap(Z,t);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function  [obj, fit] = objfunc(theta, par)
% [obj, fit] = feval(par.post, theta', par); % row theta should be supplied
% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function  [t, f, fit, itpar] = start(t0, lo, up, par)
function  [t, f,itpar] = start(fname,t0, lo, up)

% Get starting point and iteration parameters
    global Z;
    global varInputs;
% Initialize
    t = t0(:);  lo = lo(:);   up = up(:);   p = length(t);
    D = 2 .^ ([1:p]'/(p+2));
    ee = find(up == lo);  % Equality constraints
    if  ~isempty(ee)
        D(ee) = ones(length(ee),1);   t(ee) = up(ee);
    end
    ng = find(t < lo | up < t);  % Free starting values
    if  ~isempty(ng)
        %t(ng) = (lo(ng) .* up(ng).^7).^(1/8);  % Starting point
        t(ng) = (lo(ng)+up(ng))/2;  % Starting point
    end
    ne = find(D ~= 1);
    
    % Check starting point and initialize performance info
    % [f  fit] = objfunc(t,par);   nv = 1;
    t = rewrap(Z,t);
    [f,tmp] = feval(fname, t, varInputs{:});
    t = unwrap(t);
    nv = 1;
    
    itpar = struct('D',D, 'ne',ne, 'lo',lo, 'up',up, ...
        'perf',zeros(p+2,200*p), 'nv',1);
    itpar.perf(:,1) = [t; f; 1];
    
    if  isinf(f)    % Bad parameter region
        return
    end
    
    
    if  length(ng) > 1  % Try to improve starting guess
        d0 = 16;  d1 = 2;   q = length(ng);
        th = t;   fh = f;   jdom = ng(1);
        
        for  k = 1 : q
            
            j = ng(k);    fk = fh;  tk = th;
            DD = ones(p,1);  DD(ng) = repmat(1/d1,q,1);  DD(j) = 1/d0;
            alpha = min(log(lo(ng) ./ th(ng)) ./ log(DD(ng))) / 5;
            
            v = DD .^ alpha;   tk = th;
            
            for  rept = 1 : 4
                tt = tk .* v;
                tt = rewrap(Z,tt);
                ff = feval(fname, tt, varInputs{:});
                tt = unwrap(tt);
                nv = nv+1;
                
                itpar.perf(:,nv) = [tt; ff; 1];
                
                if  ff <= fk
                    tk = tt;  fk = ff;
                    if  ff <= f
                        t = tt;  f = ff;  jdom = j;
                    end
                else
                    itpar.perf(end,nv) = -1;   break
                end
            end
        end % improve
        
        % Update Delta
        if  jdom > 1
            D([1 jdom]) = D([jdom 1]);
            itpar.D = D;
        end
    end % free variables
    
    itpar.nv = nv;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function  [t, f, fit, itpar] = explore(t, f, fit, itpar, par)
function  [t, f, itpar] = explore(fname,t, f, itpar)
% Explore step
    global Z;
    global varInputs;
    nv = itpar.nv;   ne = itpar.ne;
    for  k = 1 : length(ne)
        j = ne(k);   tt = t;   DD = itpar.D(j);
        
        if  t(j) == itpar.up(j)
            atbd = 1;   tt(j) = t(j) / sqrt(DD);
        elseif  t(j) == itpar.lo(j)
            atbd = 1;  tt(j) = t(j) * sqrt(DD);
        else
            atbd = 0;  tt(j) = min(itpar.up(j), t(j)*DD);
        end
        tt = rewrap(Z,tt);
        [ff,tmp] = feval(fname, tt, varInputs{:});
        tt = unwrap(tt);
        nv = nv+1;
        itpar.perf(:,nv) = [tt; ff; 2];
        
        if  ff < f
            t = tt;  f = ff; 
        else
            itpar.perf(end,nv) = -2;
            if  ~atbd  % try decrease
                tt(j) = max(itpar.lo(j), t(j)/DD);
                tt = rewrap(Z,tt);
                [ff,tmp] = feval(fname, tt, varInputs{:});  
                tt = unwrap(tt);
                nv = nv+1;
            
                itpar.perf(:,nv) = [tt; ff; 2];
                
                if  ff < f
                    t = tt;  f = ff; 
                else
                    itpar.perf(end,nv) = -2;
                end
            end
        end
    end % k
    
    itpar.nv = nv;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function  [t, f, fit, itpar] = move(th, t, f, fit, itpar, par)
function  [t, f, itpar] = move(fname,th, t, f, itpar)
    % Pattern move
    global Z;
    global varInputs;
    nv = itpar.nv;   ne = itpar.ne;   p = length(t);
    v = t ./ th;
    if  all(v == 1)
        itpar.D = itpar.D([2:p 1]).^.2;
        return
    end
    
    % Proper move
    rept = 1;
    while  rept
        tt = min(itpar.up, max(itpar.lo, t .* v));
        
        tt = rewrap(Z,tt);
        [ff,tmp] = feval(fname, tt, varInputs{:});  
        tt = unwrap(tt);
        nv = nv+1;
        
        
        itpar.perf(:,nv) = [tt; ff; 3];
        if  ff < f
            t = tt;  f = ff; 
            v = v .^ 2;
        else
            itpar.perf(end,nv) = -3;
            rept = 0;
        end
        
        if  any(tt == itpar.lo | tt == itpar.up), rept = 0; end
    end
    
    itpar.nv = nv;
    itpar.D = itpar.D([2:p 1]).^.25;

end



function v = unwrap(s)
% Extract the numerical values from "s" into the column vector "v". The
% variable "s" can be of any type, including struct and cell array.
% Non-numerical elements are ignored. See also the reverse rewrap.m. 

    v = [];   
    if isnumeric(s)
        v = s(:);                        % numeric values are recast to column vector
    elseif isstruct(s)
        v = unwrap(struct2cell(orderfields(s))); % alphabetize, conv to cell, recurse
    elseif iscell(s)
        for i = 1:numel(s)             % cell array elements are handled sequentially
            v = [v; unwrap(s{i})];
        end
    end                                                   % other types are ignored
end

function [s v] = rewrap(s, v)
% Map the numerical elements in the vector "v" onto the variables "s" which can
% be of any type. The number of numerical elements must match; on exit "v"
% should be empty. Non-numerical entries are just copied. See also unwrap.m.
    if isnumeric(s)
        if numel(v) < numel(s)
            error('The vector for conversion contains too few elements')
        end
        
        s = reshape(v(1:numel(s)), size(s));            % numeric values are reshaped
        v = v(numel(s)+1:end);                        % remaining arguments passed on
    elseif isstruct(s) 
        [s p] = orderfields(s); p(p) = 1:numel(p);      % alphabetize, store ordering
        [t v] = rewrap(struct2cell(s), v);                 % convert to cell, recurse
        s = orderfields(cell2struct(t,fieldnames(s),1),p);  % conv to struct, reorder
    elseif iscell(s)
        for i = 1:numel(s)             % cell array elements are handled sequentially
            [s{i} v] = rewrap(s{i}, v);
        end
    end
end
