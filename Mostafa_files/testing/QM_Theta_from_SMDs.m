function Theta = QM_Theta_from_SMDs(SMDs, names)

    n = size(SMDs,1);
    m = max(names(:));
    I = names(:,1);
    J = names(:,2);

    Theta = zeros(n,m,m);
    for k = 1:size(SMDs,2)
        Theta(:,I(k),J(k)) = SMDs(:,k);
        Theta(:,J(k),I(k)) = SMDs(:,k);
    end

end