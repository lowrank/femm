function logger(c, s)
    global logger_level
    %
    %  '*[0 ,0.8, 0]': g
    %  '*[1, 0.6, 0.2]': o
    %  '*[0, 1., 1.]' : c
    %  ''*[1, 1, 0.2]'': m
    %  '*[1, 0 ,0.]' : r
    if logger_level > 0
        out = strcat({datestr(datetime('now'))}, {' '}, {s});
        cprintf(c, out{1});
        cprintf('\n');
    end
end