function [models, labels] = models_from_tests(tests)
if isfield(tests,'models')
    models = tests.models; labels = tests.labels;
else
    models = { {tests.X1}, {tests.X2}, {tests.X3} };
    labels = ensure_labels(tests,'labels_tests',{'X1','X2','X3'});
end
end