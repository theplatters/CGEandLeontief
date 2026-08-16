"""Fix common notebook issues after round-trip sync from Mac."""
import json, sys

def fix_notebook(path):
    with open(path) as f:
        nb = json.load(f)

    changed = False
    for cell in nb['cells']:
        if cell['cell_type'] == 'code':
            for i, line in enumerate(cell['source']):
                # Fix double-backslash in Julia operator: \\ -> \
                if 'data.Ω)' in line and '\\\\' in line:
                    # Python str after JSON decode: \\\\ in JSON -> \\ in Python
                    # Replace \\ (2 chars) -> \ (1 char)
                    new_line = line.replace('\\\\', '\\')
                    print(f"  FIXED: {line.strip()!r} -> {new_line.strip()!r}")
                    cell['source'][i] = new_line
                    changed = True

        elif cell['cell_type'] == 'markdown':
            src_joined = ''.join(cell['source'])
            if '5-10 min' in src_joined:
                src_joined = src_joined.replace(
                    'It takes ~5-10 min on a Mac; use the container wrapper for smaller tests.',
                    'It runs in ~2-3 sec on a Mac (only 76 cheap equilibrium solves for the Hessian; MC is matrix ops).'
                )
                cell['source'] = [src_joined]
                print(f"  FIXED timing comment")
                changed = True

    if changed:
        with open(path, 'w') as f:
            json.dump(nb, f, indent=1, ensure_ascii=False)
        print(f"Wrote {path}")
    else:
        print(f"No changes needed in {path}")

if __name__ == '__main__':
    fix_notebook(sys.argv[1])