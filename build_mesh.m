function [TR, Db, nC, d, nE, dNodes, fNodes, s,m,vol_T,mp_T, r_c4n, theta_c4n] = build_mesh(c4n, n4e, R)

% Sets up all the bookkeeping for the FEM mesh defined by [c4n, n4e].

    cla,patch('vertices',c4n,'faces',n4e,'edgecol','k','facecol',[.8,.9,1]);
    xlim([-R R])
    ylim([-R R])
    axis off
    drawnow
    % pause

    TR = triangulation(n4e,c4n);
    Db = freeBoundary(TR);
    [nC,d]  = size(c4n);            % number of nodes
    nE      = size(n4e,1);          % number of elements
    dNodes  = unique(Db);           % Dirichlet boundary
    fNodes  = setdiff(1:nC,dNodes); % free nodes
    [s,m,vol_T,mp_T] = fe_matrices(c4n,n4e);
    r_c4n     = sqrt(c4n(:,1).^2+c4n(:,2).^2);
    theta_c4n = atan2(c4n(:,2),c4n(:,1));
end