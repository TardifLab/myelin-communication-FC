function names = canonical_base_predictors(tests)
% Both stages and interaction builders use this canonical base order.
% Either input ordering is accepted; other base compositions fail explicitly.
requested = getopt(tests,'base_predictors',{'binary','caliber'});
requested = lower(string(requested(:)));
assert(numel(requested)==2 && isequal(sort(requested),["binary";"caliber"]), ...
    'bicsdual:BasePredictors', ...
    'Dual BICS requires binary and caliber once each; ED is added internally.');
names = {'binary','caliber'};
end
