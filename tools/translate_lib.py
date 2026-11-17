"""Translates the user-facing string literals of the helper library by rewriting the SOURCE (positions come from the AST), so the
translated library keeps exactly the same code. Fragments of f-strings are translated one by one (keep the word order compatible)."""
import ast


def _literal(text, quote, in_fstring):
    esc = text.replace("\\", "\\\\").replace("\n", "\\n").replace("'", "\\'").replace('"', '\\"')
    if in_fstring:
        esc = esc.replace("{", "{{").replace("}", "}}")
        return esc
    return repr(text)


def translate_source(src, table):
    """Returns (new_source, used_keys). table: {english literal: spanish literal}."""
    tree = ast.parse(src)
    used_extra = set()
    lines = src.split("\n")
    offs = [0]
    for ln in lines:
        offs.append(offs[-1] + len(ln) + 1)
    pos = lambda l, c: offs[l - 1] + len(lines[l - 1].encode()[:c].decode())
    fstring_consts, jobs, skip = set(), [], set()
    for n in ast.walk(tree):
        if isinstance(n, ast.JoinedStr):
            multi = n.end_lineno != n.lineno                      # several adjacent literals: the positions of its pieces are unreliable
            rebuilt = False
            for v in n.values:
                if isinstance(v, ast.Constant):
                    fstring_consts.add(id(v))
                    if multi and v.value in table:
                        rebuilt = True
            if rebuilt:
                parts = []
                for v in n.values:
                    if isinstance(v, ast.Constant):
                        parts.append(_literal(table.get(v.value, v.value), None, True))
                        if v.value in table:
                            used_extra.add(v.value)
                    else:
                        seg = ast.get_source_segment(src, v)
                        base = pos(v.lineno, v.col_offset)
                        inner = [(pos(c.lineno, c.col_offset) - base, pos(c.end_lineno, c.end_col_offset) - base, c.value) for c in ast.walk(v)
                                 if isinstance(c, ast.Constant) and isinstance(c.value, str) and c.value in table]
                        for a2, b2, key in sorted(inner, reverse=True):          # strings inside the {expression}: translate them in place
                            seg = seg[:a2] + repr(table[key]) + seg[b2:]
                            used_extra.add(key)
                            skip.add(("inner", key, a2 + base))
                        parts.append(seg)
                    skip.add(id(v))
                    for c in ast.walk(v):
                        skip.add(id(c))
                jobs.append((pos(n.lineno, n.col_offset), pos(n.end_lineno, n.end_col_offset), "f'''" + "".join(parts) + "'''", "whole"))
    for n in ast.walk(tree):
        if isinstance(n, ast.Constant) and isinstance(n.value, str) and n.value in table and id(n) not in skip:
            jobs.append((pos(n.lineno, n.col_offset), pos(n.end_lineno, n.end_col_offset), n.value, id(n) in fstring_consts))
    used = set(used_extra)
    out = src
    for a, b, key, infs in sorted(jobs, key=lambda j: (j[0], j[1]), reverse=True):
        if infs == "whole":
            out = out[:a] + key + out[b:]
            continue
        used.add(key)
        if infs:
            raw = out[a:b]
            expected = key.replace("{", "{{").replace("}", "}}")
            # a fragment inside an f-string: replace only its text (the slice is the raw source of the fragment)
            out = out[:a] + _literal(table[key], None, True) + out[b:]
        else:
            out = out[:a] + _literal(table[key], None, False) + out[b:]
    return out, used
