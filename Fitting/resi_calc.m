function [val_cell_viab val_viral_proj] = resi_calc(time, U, V, cell_time, viral_time, p, cell_data, viral_data, initialconds)

    Usol = interp1(time,U,cell_time);
    Vsol = interp1(time,V,viral_time);

    if isempty(find(Usol>1e10))==0
       Usol  =repmat(1e10, 1, length(cell_time))';
       Vsol = repmat(1e10, 1, length(viral_time));
    elseif sum(isnan(Usol))>0
       Usol  =repmat(1e10, 1, length(cell_time))';
       Vsol = repmat(1e10, 1, length(viral_time));
    elseif sum(isnan(Vsol))>0
       Usol  =repmat(1e10, 1, length(cell_time))';
       Vsol = repmat(1e10, 1, length(viral_time));
    end

    Usol = interp1(time,U,cell_time);
    cell_viab = (Usol)/(initialconds(1))*100;

    Vsol = interp1(time,V,viral_time);
    val_cell_viab = [cell_viab'- cell_data]./cell_data;%./max(cell_data);

    val_viral_proj = [Vsol - viral_data./p.alpha]./viral_data.*[1 1 10 10 1 1];%/max(viral_data);
    
    val = [val_cell_viab val_viral_proj];

    end