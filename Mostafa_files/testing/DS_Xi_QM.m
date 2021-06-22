function Xi = DS_Xi_QM(DSs, names)

    n = size(DSs,1);
    m = max(names(:,1));
    md = max(names(:,2));

    I = names(:,1);
    J = names(:,2);

    Xi = zeros(n,m,md);
    for k = 1:size(DSs,2)
        Xi(:,I(k),J(k)) = DSs(:,k);
    end

end