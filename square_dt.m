function sq = square_dt(n,period)
    % Square wave in discrete time
    repeats = n/period;
    sq = square(2*pi*linspace(0,repeats-1/n,n));
end
