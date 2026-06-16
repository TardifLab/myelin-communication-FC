function U = unpack_useCM(ALLX)
% Expect ALLX = {caliber, MTsat, gratio, delay, binary, ED}, each 1x6 cell, except ED N x N
U = struct();
U.caliber = ALLX{1};
U.MTsat   = ALLX{2};
U.gratio  = ALLX{3};
U.delay   = ALLX{4};
U.binary  = ALLX{5};
U.ED      = ALLX{6};
end