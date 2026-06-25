function export_docs()
    % Export all doc source files to HTML.
    % Run from the repository root.
    source_docfiles = dir(fullfile('doc', '*.m'));
    source_docfiles = source_docfiles(~strcmp({source_docfiles.name}, 'export_docs.m'));
    [~, docnames, ~] = cellfun(@fileparts, {source_docfiles.name}, 'UniformOutput', false);
    
    for n = 1:numel(docnames)
        src  = fullfile('doc', [docnames{n} '.m']);
        dest = fullfile('src', 'html', [docnames{n} '.html']);
        export(src, dest, Format="html");
        fprintf('%d/%d  %s\n', n, numel(docnames), docnames{n});
    end
end