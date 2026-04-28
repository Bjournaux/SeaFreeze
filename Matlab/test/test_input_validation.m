function test_input_validation()
% Smoke-test that bad inputs to SeaFreeze and SF_WhichPhase produce the
% expected SeaFreeze:* / SF_WhichPhase:* errors instead of cryptic crashes.
% Baptiste Journaux - 2026

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'LocalBasisFunction'));
addpath(fileparts(here));

cases = {
    % --- Bad material -------------------------------------------------------
    {'unknown material',          @() SF_getprop([100 280], 'IceX'),                         'SeaFreeze:unknownMaterial'}
    {'non-string material',       @() SF_getprop([100 280], 5),                              'SeaFreeze:badInput'}
    % --- Bad PT shape -------------------------------------------------------
    {'wrong cell length (water)', @() SF_getprop({100, 280, 0.5}, 'water1'),                 'SeaFreeze:badInput'}
    {'wrong cell length (NaCl)',  @() SF_getprop({100, 280}, 'NaClaq'),                      'SeaFreeze:badInput'}
    {'scatter wrong cols (water)',@() SF_getprop([100 280 0.5], 'water1'),                   'SeaFreeze:badInput'}
    {'scatter wrong cols (NaCl)', @() SF_getprop([100 280], 'NaClaq'),                       'SeaFreeze:badInput'}
    {'PT not numeric',            @() SF_getprop('hello', 'water1'),                         'SeaFreeze:badInput'}
    {'PT empty',                  @() SF_getprop(zeros(0,2), 'water1'),                      'SeaFreeze:badInput'}
    {'PT contains NaN',           @() SF_getprop([100 NaN], 'water1'),                       'SeaFreeze:badInput'}
    {'PT cell contains NaN',      @() SF_getprop({[100 200], [NaN 280]}, 'water1'),          'SeaFreeze:badInput'}
    % --- Bad property names ------------------------------------------------
    {'unknown property',          @() SF_getprop([100 280], 'water1', 'banana'),             'SeaFreeze:unknownProperty'}
    {'shear on liquid',           @() SF_getprop([100 280], 'water1', 'Vp'),                 'SeaFreeze:unknownProperty'}
    {'mixing on ice',             @() SF_getprop([100 250], 'Ih', 'mus'),                    'SeaFreeze:unknownProperty'}
    {'props wrong type',          @() SF_getprop([100 280], 'water1', 42),                   'SeaFreeze:badInput'}
    % --- SF_WhichPhase ------------------------------------------------------
    {'WhichPhase bad solute',     @() SF_WhichPhase({0.1,280}, 'solute', 'KCl'),            'SF_WhichPhase:badInput'}
    {'WhichPhase NaCl no m',      @() SF_WhichPhase({0.1,280}, 'solute','NaCl'),            'SeaFreeze:badInput'}
};

n_pass = 0; n_fail = 0;
for i = 1:size(cases,1)
    name = cases{i}{1}; fn = cases{i}{2}; want_id = cases{i}{3};
    try
        fn();
        fprintf('  [FAIL] %-32s  (no error thrown)\n', name);
        n_fail = n_fail + 1;
    catch err
        if strcmp(err.identifier, want_id)
            fprintf('  [pass] %-32s  -> %s\n', name, err.identifier);
            n_pass = n_pass + 1;
        else
            fprintf('  [FAIL] %-32s  got %s, want %s\n', ...
                    name, err.identifier, want_id);
            n_fail = n_fail + 1;
        end
    end
end

% --- Sanity: valid inputs still work --------------------------------------
try
    SF_getprop([100 280], 'water1', 'rho');
    SF_getprop({0.1:50:200, 273:5:300, [0.1 0.5]}, 'NaClaq', {'rho','Cp'});
    SF_WhichPhase({0.1, 280});
    SF_WhichPhase({0.1, 280, 1.0}, 'solute', 'NaCl');
    fprintf('  [pass] valid inputs still work\n');
    n_pass = n_pass + 1;
catch err
    fprintf('  [FAIL] valid inputs raised: %s (%s)\n', err.message, err.identifier);
    n_fail = n_fail + 1;
end

% --- Deprecation: SeaFreeze() should warn once and return same result -----
try
    clear SeaFreeze   % reset the persistent `warned` flag
    s = warning('on', 'SeaFreeze:deprecated');
    cleanup = onCleanup(@() warning(s));
    lastwarn('', '');
    a = SeaFreeze([100 280], 'water1', 'rho');
    [msg, id] = lastwarn;
    if ~strcmp(id, 'SeaFreeze:deprecated')
        error('expected SeaFreeze:deprecated warning, got id=''%s''', id);
    end
    b = SF_getprop([100 280], 'water1', 'rho');
    if abs(a.rho - b.rho) > 1e-12 * abs(b.rho)
        error('SeaFreeze and SF_getprop returned different rho values');
    end
    fprintf('  [pass] SeaFreeze deprecation warning + result match\n');
    n_pass = n_pass + 1;
catch err
    fprintf('  [FAIL] deprecation alias: %s\n', err.message);
    n_fail = n_fail + 1;
end

fprintf('\n%d passed, %d failed\n', n_pass, n_fail);
if n_fail > 0
    error('test_input_validation:fail', '%d validation cases failed', n_fail);
end
end
