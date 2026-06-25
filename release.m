%release.m  DynaProg release script.
% Run from the repository root in MATLAB.

% Pre-requisites: 
% - Create and push a tag with the version number (e.g. "v1.7.0") on origin/main
% - Ensure the changelog is up-to-date

%% Read version
% Read version from latest tag on main remote branch
[status, out] = system("git describe --abbrev=0 origin/main");
if status ~= 0
    error("git describe failed: %s", strtrim(out));
end
tag_fromRepo = strtrim(out);

% Read version from source
src = fileread(fullfile('src', '@DynaProg', 'DynaProg.m'));
tok = regexp(src, 'Version\s*=\s*"([^"]+)"', 'tokens', 'once');
assert(~isempty(tok), 'Could not parse Version from DynaProg.m.');
version = tok{1};
tag_fromSource = ['v' version];

if compareVersions(tag_fromSource, tag_fromRepo)  ~= 0
    error("Tag from DynaProg.m (%s) must be equal to the latest tag " + ...
        "on GitHub (%s). Ensure both are up-to-date.\n", ...
        tag_fromSource, tag_fromRepo);
else
    tag = tag_fromSource;
end

[status, ~] = system(sprintf("gh release view %s", tag));
if status == 0
    error("Release %s already exists on the GitHub repository", tag)
end
% Check the changelog
notes = readChangelog(tag);
release_notes = sprintf('%s', notes);

%% 1. Export documentation
fprintf('[1/4] Creating html documentation...\n');
export_docs();
fprintf('Done.\n\n');

%% 2. Run tests
fprintf('[2/4] Running tests...\n');
results  = runtests('tests');
n_failed = sum([results.Failed]);
if n_failed > 0
    error('DynaProg:release', '%d test(s) failed - aborting.', n_failed);
end
fprintf('%d tests passed.\n\n', numel(results));

%% 3. Package toolbox
fprintf('[3/4] Packaging toolbox...\n');

% Stable identifier for the File Exchange / Add-On Manager.
tbx_ID = '3df3f3df-feb6-4cfb-84f3-df0529d46b38';
toolboxFolder = "src";

% Initialize empty options
opts = matlab.addons.toolbox.ToolboxOptions('src', tbx_ID);

% Set options
opts.ToolboxName = "DynaProg";
opts.ToolboxVersion = version;

opts.Description = "DynaProg provides a flexible tool to solve a " + ...
    "finite horizon multi-stage deterministic decision problem, " + ...
    "which is a problem where a decision must be made at each stage " + ...
    "for a system that evolves through a finite number of stages, " + ...
    "minimizing the total cost incurred. The toolbox allows you to " + ...
    "define your own optimization problem and attempts to solve it " + ...
    "using Dynamic Programming." + ...
    "This toolbox comes with html documentation integrated in the " + ...
    "help browser.";
opts.Summary = "Solve finite horizon multi-stage deterministic " + ...
    "decision problems using Dynamic Programming.";

opts.AuthorName = "Federico Miretti";
opts.AuthorEmail = "federico.miretti@polito.it";
opts.ToolboxImageFile = fullfile("doc", "card.png");

opts.SupportedPlatforms.Win64 = true;
opts.SupportedPlatforms.Maci64 = true;
opts.SupportedPlatforms.Glnxa64 = true;
opts.SupportedPlatforms.MatlabOnline = true;

opts.MinimumMatlabRelease = "R2020b";

mltbx_name = sprintf('DynaProg_%s.mltbx', tag);
opts.OutputFile = mltbx_name;

% Package the toolbox
matlab.addons.toolbox.packageToolbox(opts);
fprintf('Created: %s\n\n', mltbx_name);

%% 4. Create GitHub release
fprintf('[4/4] Creating GitHub release...\n');

% Write notes to a temp file to avoid shell-quoting issues with
% multi-line notes or embedded quotes/special characters.
notes_file = [tempname '.md'];
fid = fopen(notes_file, 'w');
fwrite(fid, release_notes);
fclose(fid);
cleanup = onCleanup(@() delete(notes_file));

cmd = sprintf('gh release create %s "%s" --title "%s" --notes-file "%s"', ...
    tag, mltbx_name, tag, notes_file);

[status, out] = system(cmd);

if status ~= 0
    error('DynaProg:release', 'gh release create failed:\n%s', out);
end
fprintf('Done: %s\n', strtrim(out));

%% -----------------------------------------------------------------------
function notes = readChangelog(versionTag)
text = fileread('CHANGELOG.md');
lines = splitlines(text);
header = ['# ' versionTag];

start_line = find(strcmp(lines, header), 1);
if isempty(start_line)
    error('No changelog entry found for "%s".', header);
end

rest = lines(start_line+1 : end);
next_h1 = find(startsWith(rest, '# '), 1);
if isempty(next_h1)
    section = rest;
else
    section = rest(1 : next_h1-1);
end

while ~isempty(section) && strtrim(section{1}) == ""
    section(1) = [];
end
while ~isempty(section) && strtrim(section{end}) == ""
    section(end) = [];
end

section = strrep(section, "##", "#");
notes = strjoin(section, newline);

end

function v = parseVersion(tag)
    % Extract all digit groups, ignoring the 'v' prefix and separators
    v = str2double(regexp(tag, '\d+', 'match'));
end

function c = compareVersions(a, b)
    % Returns  1 if a > b,  -1 if a < b,  0 if equal
    va = parseVersion(a);
    vb = parseVersion(b);
    n = max(numel(va), numel(vb));
    va(end+1:n) = 0;          % pad shorter one with zeros
    vb(end+1:n) = 0;          % so 'v1.6' == 'v1.6.0'
    d = va - vb;
    idx = find(d ~= 0, 1);    % first component that differs
    if isempty(idx)
        c = 0;
    else
        c = sign(d(idx));
    end
end