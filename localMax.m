function [f,X,T] = localMax(p,k)
%vyberie pri kmitaní body, kde má funkcia
%lokálne maximum
%k - poèet okolitých bodov na porovnávanie
%za lokálne maximum berie bod pre ktorý platí,
%že všetkých k*2 okolitých bodov je väèších ako 
%tento bod
%pre k=1 to isté ako funkcia ampCheck
%p - vektor, harmonický priebeh
%f - hodnoty vystupnych bodov, X - indexy lokalnych maxim
%T - perioda maxim
n = length(p);
j = 1;
pR = p(k+2:2*k+1);%body napravo
pL = p(1:k);%body na¾avo
m = 0;%index posledného maxima
T = [];
for i = k+1:n-k
    if all(pL < p(i)) && all(pR < p(i))
        f(j) = p(i);  
        pL = p(i-k+1:i);
        if i<n-k
            pR = p(i+2:i+k+1);
        end
        T(j) = i-m;
        X(j) = i;
        j = j+1;
        m=i;
    else
        if i<n-k
            pL = p(i-k+1:i);
            pR = p(i+2:i+k+1);
        end
    end
end
T = round(mean(T));%priemerný rozostup maxím
end