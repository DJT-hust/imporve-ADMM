function res= find_connect(n,connect_mat)
    
    M_set=find(connect_mat(n,:));
    num=numel(M_set);
    res=zeros(num,2);
    res=res(:,1)+n;
    res(:,2)=M_set';
    
