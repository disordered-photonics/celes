%======================================================================
%> @brief Double factorial
%>
%> @param n (int)
%> @retval f (int): f=n!!
%======================================================================
function f=factd2(n)
% f=factd2(n)
% Return the double factorial of n, i.e. n!!

f=n;
for j1=n-2:-2:1
    f=f*j1;
end

if n==0
    f=1;
end
    
    
    