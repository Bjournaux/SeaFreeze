% strip_NaClaq_splines.m  —  remove fitting metadata from NaClaq .mat files
%
% Keeps only the fields needed at runtime (knots, coefs, order, etc.),
% saves originals as *_full.mat, and re-saves stripped versions in -v7
% format for faster loading.

here = fileparts(mfilename('fullpath'));
spdir = fullfile(here, '..', 'splines');

runtime_fields = {'knots','coefs','number','order','dim','form', ...
                  'MW','nu','Tc','cutoff','mask','PTm_mask'};

required_fields = {'knots','coefs','number','order','dim'};

targets = {
    fullfile('NaCl_aq_LP_2026',    'NaCl_aq_LP_2026.mat')
    fullfile('NaCl_aq_HP_2026',    'NaCl_aq_HP_2026.mat')
    fullfile('NaCl_aq_HP_2026_v1', 'NaCl_aq_HP_2026_v1.mat')
    fullfile('NaCl_aq_HP_2026_v2', 'NaCl_aq_HP_2026_v2.mat')
    fullfile('NaCl_aq_HP_2026_v3', 'NaCl_aq_HP_2026_v3.mat')
};

fprintf('=== Stripping NaClaq spline files ===\n\n');
fprintf('%-45s  %10s  %10s  %8s\n', 'File', 'Before', 'After', 'Ratio');
fprintf('%s\n', repmat('-', 1, 78));

for k = 1:numel(targets)
    fpath = fullfile(spdir, targets{k});
    if ~isfile(fpath)
        fprintf('  SKIP (not found): %s\n', targets{k});
        continue
    end

    d = dir(fpath);
    size_before = d.bytes;

    S = load(fpath, 'sp');
    sp_orig = S.sp;

    sp = struct();
    for j = 1:numel(runtime_fields)
        fn = runtime_fields{j};
        if isfield(sp_orig, fn)
            sp.(fn) = sp_orig.(fn);
        end
    end

    for j = 1:numel(required_fields)
        fn = required_fields{j};
        assert(isfield(sp, fn) && ~isempty(sp.(fn)), ...
            'Missing required field ''%s'' in %s', fn, targets{k});
    end

    [fdir, fname, ~] = fileparts(fpath);
    backup = fullfile(fdir, [fname '_full.mat']);
    copyfile(fpath, backup);

    save(fpath, 'sp', '-v7');

    V = load(fpath, 'sp');
    assert(isequal(sp.coefs, V.sp.coefs), ...
        'Verification failed for %s: coefs mismatch after reload', targets{k});

    d2 = dir(fpath);
    size_after = d2.bytes;

    fprintf('%-45s  %8.1f MB  %8.1f KB  %6.0fx\n', ...
        targets{k}, size_before/1e6, size_after/1e3, size_before/size_after);
end

fprintf('\nDone. Originals saved as *_full.mat alongside stripped files.\n');
