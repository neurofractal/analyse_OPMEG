function dataout = fixsens(data1)
    chans2keep = ismember(data1.grad.label,data1.label);

    grad_mod         = data1.grad;

    % For doubles
    pot_grad_fields = {'chanori','chanpos','coilori','coilpos'};

    for p = 1:length(pot_grad_fields)
        if isfield(grad_mod,pot_grad_fields{p})
            eval(['grad_mod.' pot_grad_fields{p} ' = grad_mod.' pot_grad_fields{p} '(chans2keep,:)']);
        end
    end

    % For cells
    pot_grad_fields = {'chantype','chanunit','label'};

    for p = 1:length(pot_grad_fields)
        if isfield(grad_mod,pot_grad_fields{p})
            eval(['grad_mod.' pot_grad_fields{p} ...
                ' = grad_mod.' pot_grad_fields{p} '(chans2keep)']);
        end
    end

    % For tra matrix
    if isfield(grad_mod,'tra')
        grad_mod.tra = grad_mod.tra(chans2keep,chans2keep);
    end

    % Replace
    dataout=data1;
    dataout.grad = grad_mod;
end
