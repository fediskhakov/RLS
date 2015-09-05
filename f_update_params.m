function [par mp] = f_update_params(par,mp)

    % dt
    % if dt is large and close to 1, mp.df and par.nC can be controlled directly
    % if mp.dt<1
    %     mp.df=exp(-0.05*mp.dt);
    %     par.nC=ceil(1/mp.dt);
    % end
    
    % tpm
    % NB: treat second column as dependent
    mp.tpm(:,2)=[1 1]'-mp.tpm(:,1);

    % pti
    par.pti=f_pti(mp, par);

    % add nC to mp
    mp.nC=par.nC;

end %function
