    n=8
    c = zeros(n + 1, 1);
    for k = 0 : 1 : (n + 1) / 2
        c(2 * k + 1)=(-1)^k * factorial(2 * (n + 1) - 2 * k) / (factorial((n + 1) - 2 * k)*factorial(n + 1 - k) * factorial(k) * 2^(n + 1));
        if k ~= (n + 1) / 2 
            c(2 * k + 2)=0;
        end
    end
    
    str="";
    for i = 1 : n + 1
        str = strcat(str,num2str(c(i)),", ");
        
    end  
    str=str(1:length(str)-1);
    str
    
    t2=roots(c); 

    t2    
    t1=sort(t2);

    options = optimset('TolFun',1e-100, 'Display','off');
    t = zeros(n + 1, 1);
    for i = 1 : n + 1
        t(i) = fsolve(@(x)LegendrePolinom(n + 1, 0, x), t1(i), options);
    end  
    
    str="";
    for i = 1 : n + 1
        str = strcat(str,num2str(t(i)),", ");
        
    end  
    str=str(1:length(str)-1);
    str
    
