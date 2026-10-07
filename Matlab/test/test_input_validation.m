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
    {'wrong cell length (water)', @() SF_getprop({100, 280, 0.5}, 'water_Bollengier2019'),                 'SeaFreeze:badInput'}
    {'wrong cell length (NaCl)',  @() SF_getprop({100, 280}, 'NaClaq'),                      'SeaFreeze:badInput'}
    {'wrong cell length (NaCl_LP)',@() SF_getprop({100, 280}, 'NaClaq_LP'),                  'SeaFreeze:badInput'}
    {'scatter wrong cols (water)',@() SF_getprop([100 280 0.5], 'water_Bollengier2019'),                   'SeaFreeze:badInput'}
    {'scatter wrong cols (NaCl)', @() SF_getprop([100 280], 'NaClaq'),                       'SeaFreeze:badInput'}
    {'scatter wrong cols (NaCl_HP)',@() SF_getprop([100 280], 'NaClaq_HP'),                  'SeaFreeze:badInput'}
    {'PT not numeric',            @() SF_getprop('hello', 'water_Bollengier2019'),                         'SeaFreeze:badInput'}
    {'PT empty',                  @() SF_getprop(zeros(0,2), 'water_Bollengier2019'),                      'SeaFreeze:badInput'}
    {'PT contains NaN',           @() SF_getprop([100 NaN], 'water_Bollengier2019'),                       'SeaFreeze:badInput'}
    {'PT cell contains NaN',      @() SF_getprop({[100 200], [NaN 280]}, 'water_Bollengier2019'),          'SeaFreeze:badInput'}
    % --- Bad property names ------------------------------------------------
    {'unknown property',          @() SF_getprop([100 280], 'water_Bollengier2019', 'banana'),             'SeaFreeze:unknownProperty'}
    {'shear on liquid',           @() SF_getprop([100 280], 'water_Bollengier2019', 'Vp'),                 'SeaFreeze:unknownProperty'}
    {'mixing on ice',             @() SF_getprop([100 250], 'Ih', 'mus'),                    'SeaFreeze:unknownProperty'}
    {'props wrong type',          @() SF_getprop([100 280], 'water_Bollengier2019', 42),                   'SeaFreeze:badInput'}
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
    SF_getprop([100 280], 'water_Bollengier2019', 'rho');
    SF_getprop({0.1:50:200, 273:5:300, [0.1 0.5]}, 'NaClaq',          {'rho','Cp'});
    SF_getprop({0.1:50:200, 273:5:300, [0.1 0.5]}, 'NaClaq_LP',       {'rho','Cp'});
    SF_getprop({1000:500:3000, 300:100:500, [0.1 0.5]}, 'NaClaq_HP',  {'rho','Cp'});
    SF_getprop({0.1:50:200, 273:5:300, [0.1 0.5]}, 'NaClaq_5GPa_2024',{'rho','Cp'});
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
    % Escalate the warning to an error to detect it: lastwarn is unreliable under
    % Octave, which also records disabled warnings (e.g. from its own fullfile).
    s = warning('error', 'SeaFreeze:deprecated');
    cleanup = onCleanup(@() warning(s));
    id = '';
    try
        SeaFreeze([100 280], 'water_Bollengier2019', 'rho');
    catch w
        id = w.identifier;
    end
    if ~strcmp(id, 'SeaFreeze:deprecated')
        error('expected SeaFreeze:deprecated warning, got id=''%s''', id);
    end
    warning('off', 'SeaFreeze:deprecated');
    a = SeaFreeze([100 280], 'water_Bollengier2019', 'rho');
    b = SF_getprop([100 280], 'water_Bollengier2019', 'rho');
    if abs(a.rho - b.rho) > 1e-12 * abs(b.rho)
        error('SeaFreeze and SF_getprop returned different rho values');
    end
    fprintf('  [pass] SeaFreeze deprecation warning + result match\n');
    n_pass = n_pass + 1;
catch err
    fprintf('  [FAIL] deprecation alias: %s\n', err.message);
    n_fail = n_fail + 1;
end

% --- Renamed materials (1.2): old names warn once per session, same results --
DM = 'SeaFreeze:deprecatedMaterial';
renamed = {'water1', 'water_Bollengier2019', [100 280]; ...
           'water2', 'water_Brown2018',      [1000 400]; ...
           'water3', 'water_Brown2026',      [0.1 300]};
for k = 1:size(renamed, 1)
    old = renamed{k,1}; new = renamed{k,2}; PT = renamed{k,3};
    try
        clear sf_material_name                       % a fresh session for the warning
        w1 = issues_warning(@() SF_getprop(PT, old, {'rho','G'}), DM);
        w2 = issues_warning(@() SF_getprop(PT, old, {'rho','G'}), DM);   % second use: silent
        w3 = issues_warning(@() SF_getprop(PT, new, {'rho','G'}), DM);   % new name: never
        s = warning('off', DM);
        a = SF_getprop(PT, old, {'rho','G'}); b = SF_getprop(PT, new, {'rho','G'});
        warning(s);
        if ~(w1 && ~w2 && ~w3), error('warning pattern %d%d%d, expected 100', w1, w2, w3); end
        if ~isequal(a.rho, b.rho) || ~isequal(a.G, b.G), error('results differ from %s', new); end
        fprintf('  [pass] material %s -> %s: warns once, identical results\n', old, new);
        n_pass = n_pass + 1;
    catch err
        fprintf('  [FAIL] material %s -> %s: %s\n', old, new, err.message);
        n_fail = n_fail + 1;
    end
end
try
    clear sf_material_name
    s = warning('off', DM);
    c1 = SF_PhaseLines('Ih', 'water1', 'segment', 'stable');
    c2 = SF_PhaseLines('Ih', 'water_Bollengier2019', 'segment', 'stable');
    r1 = SF_rho2P(1000, 300, 'water1'); r2 = SF_rho2P(1000, 300, 'water_Bollengier2019');
    q1 = SF_coexistence('saturation', 400, 'fluid', 'water3'); q2 = SF_coexistence('saturation', 400);
    g1 = SF_phase_range('water2'); g2 = SF_phase_range('water_Brown2018');
    h1 = SF_WhichPhase({[0.1 300], [250 300]}, 'liquid', 'water3');
    h2 = SF_WhichPhase({[0.1 300], [250 300]}, 'liquid', 'water_Brown2026');
    warning(s);
    ok = isequal(c1.T, c2.T) && isequal(r1, r2) && isequal(q1.P, q2.P) && isequal(g1, g2) && isequal(h1, h2);
    if ~ok, error('an entry point gave different results for the old name'); end
    fprintf('  [pass] old names in SF_PhaseLines / SF_rho2P / SF_coexistence / SF_phase_range / SF_WhichPhase\n');
    n_pass = n_pass + 1;
catch err
    fprintf('  [FAIL] old names in the other entry points: %s\n', err.message);
    n_fail = n_fail + 1;
end
clear sf_material_name

fprintf('\n%d passed, %d failed\n', n_pass, n_fail);
if n_fail > 0
    error('test_input_validation:fail', '%d validation cases failed', n_fail);
end
end


function warned = issues_warning(fn, id)
% True if fn issues warning `id` (escalated to an error for the call; lastwarn
% is unreliable under Octave, which also records disabled warnings).
s0 = warning('query', id);   % restore only this id: Octave keeps a per-id
                             % 'error' state through a whole-state restore
warning('error', id);
try
    fn();
    warned = false;
catch e
    warned = strcmp(e.identifier, id);
    if ~warned, warning(s0.state, id); rethrow(e); end
end
warning(s0.state, id);
end
