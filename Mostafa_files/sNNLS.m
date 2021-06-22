function [wei]=sNNLS(G,b,tole)
EE=[];
Z=[1:size(G,2)];
wei=sparse(zeros(size(G,2),1));
zeta=sparse(zeros(size(G,2),1));
fullset=[1:size(G,2)];

 cntt=0;
 neta=0;
while norm(G*wei-b,'fro')>tole*norm(b,'fro')
    muu=G'*(b-G*wei);
    [vv,ee]=max(muu);
    vv=abs(vv);
    EE=union(EE,ee);
   
    Z=setdiff(Z,ee);
    while true
        GE=full(G(:,EE));
        zeta(EE)=pinv(GE)*b;
        zeta(Z)=0;
        if zeta(EE)>0
            wei=zeta;
            break
        end
        
        for ii=EE
            netaMin1=min(wei(ii)/(wei(ii)-zeta(ii)));
            if cntt==0
                neta=netaMin1;
                cntt=cntt+1;
            elseif netaMin1<neta
                neta=netaMin1  ;
            end
        end
        wei=wei+neta*(zeta-wei);
        Z=find(~wei);
         EE=setdiff(fullset,Z);
    end
        
end
end
