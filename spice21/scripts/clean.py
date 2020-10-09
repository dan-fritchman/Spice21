
"""
Our home-grown `cargo fix`, that does what we want. 
"""

import glob 
import re 

skip = ["Length", "Width", "Type", "Vth", "MAX", "MIN"]

fh = open('test.log', 'r')
replacements = {}
line = fh.readline()
while line:
    if 'warning' in line and 'should have a snake case name' in line:
        splits = line.split('`')
        if len(splits) != 3: 
            print(f"Invalid snake-case target: {line}")
            assert False, line
            line = fh.readline()
            continue 
        target = splits[1]
        line = fh.readline()
        while True:
            if 'warning' in line: 
                print(f"Invalid snake-case warning not closed: {line}")
                assert False, line
                break 
            if 'snake case' in line and '`' in line:
                splits = line.split('`')
                if len(splits) != 3: 
                    print(f"Invalid snake-case replacement detected: {line}")
                    assert False, line
                    line = fh.readline()
                    continue 
                suggestion = splits[1]
                if target not in replacements and target not in skip:
                    if len(target) > 2:
                        print(f"{target} -> {suggestion}")
                        replacements[target] = suggestion 
                    else:
                        print(f"Target short enough for concern: {target}")
                break 
            line = fh.readline()
    line = fh.readline()
    
# print(replacements)
print(len(replacements))

files = glob.glob('src/**/analysis.rs', recursive=True)
print(files)

for f in files:
    txt = open(f, 'r').read()
    # Sort by length, longest to shortest, to avoid collisions 
    for k in sorted(list(replacements.keys()), key=len, reverse=True):
        # r = re.compile(f"([^A-Za-z0-9_]){k}([^A-Za-z0-9_])")
        # txt = r.sub(rf"{\1}{k}{\2}", txt)
        txt = txt.replace(k, replacements[k])
    open(f,'w').write(txt)

