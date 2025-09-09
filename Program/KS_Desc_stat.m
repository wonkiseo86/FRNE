function Sum_stat = KS_Desc_stat(x_t,Fn)

s=x_t;
sumstats = [] ;
for i=1:size(Fn,1)

    f = Fn(i,:)';
    ff= s.*f ;
    mm= mrsum(s,ff);
    
    ss=s;
    for k=2:4
       ss = [ss (s-mm).^k] ; 
    end
    ff2 = ss.*repmat(f,1,size(ss,2)) ;
    a0  = [(1+i-1) mrsum(s,ff2)'] ;
    sumstats = [sumstats; a0];
end
sumstats(:,1) = cumsum(ones(size(Fn,1),1));

c1 = sumstats(:,2); % mean
c2 = sumstats(:,3); % variance
c3 = sumstats(:,4)./(c2.^(3/2)) ; % skewness
c4 = sumstats(:,5)./(c2.^2) ; % kurtosis

Sum_stat = [c1 c2 c3 c4];