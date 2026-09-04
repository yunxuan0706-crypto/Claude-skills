"""Turn the matplotlib mathtext used in the figure into PowerPoint text runs.

Output is a list of {t, sub, sup, i, b} dicts: PowerPoint has no stacked
super+subscript, so x^{(a)}_{sir} becomes three sequential runs.
"""
SYM = {
    "Phi": "Φ", "Psi": "Ψ", "psi": "ψ", "phi": "φ",
    "beta": "β", "lambda": "λ", "theta": "θ", "delta": "δ",
    "rho": "ρ", "mu": "μ", "alpha": "α", "gamma": "γ",
    "epsilon": "ε", "Rightarrow": "⇒", "rightarrow": "→",
    "to": "→", "infty": "∞", "geq": "≥", "leq": "≤",
    "times": "×", "pm": "±", "langle": "⟨", "rangle": "⟩",
    "ldots": "…", "cdot": "·", "approx": "≈", "neq": "≠",
}
UPRIGHT = set("ΦΨ∞⇒→≥≤×±⟨⟩…·≈≠")
SPACE = {",": " ", ";": " ", "!": "", " ": " ", ":": " "}


def _emit(out, t, sub, sup, i=False, b=False):
    if not t:
        return
    WS = (" ", "\u2009")          # \, renders as a thin space, not a plain one
    if out and t.startswith(WS) and out[-1]["t"].endswith(WS):
        t = t.lstrip("".join(WS))
        if not t:
            return
    if out and out[-1]["sub"] == sub and out[-1]["sup"] == sup \
            and out[-1]["i"] == i and out[-1]["b"] == b:
        out[-1]["t"] += t
    else:
        out.append({"t": t, "sub": sub, "sup": sup, "i": i, "b": b})


def _group(src, k):
    """Read the argument after _ ^ or a command: {..} or one character."""
    if k < len(src) and src[k] == "{":
        depth, j = 1, k + 1
        while j < len(src) and depth:
            depth += {"{": 1, "}": -1}.get(src[j], 0)
            j += 1
        return src[k + 1:j - 1], j
    return (src[k], k + 1) if k < len(src) else ("", k)


def _math(src, out, sub=False, sup=False):
    k = 0
    while k < len(src):
        c = src[k]
        if c == "\\":
            j = k + 1
            while j < len(src) and src[j].isalpha():
                j += 1
            name = src[k + 1:j]
            if not name:                                  # \, \; \! \ escapes
                _emit(out, SPACE.get(src[j] if j < len(src) else " ", " "), sub, sup)
                k = j + 1
            elif name == "mathbf":
                g, k = _group(src, j)
                _emit(out, g, sub, sup, i=False, b=True)
            elif name in SYM:
                ch = SYM[name]
                _emit(out, ch, sub, sup, i=ch not in UPRIGHT)
                k = j
            else:
                _emit(out, name, sub, sup)
                k = j
        elif c in "_^":
            g, k = _group(src, k + 1)
            _math(g, out, sub or c == "_", sup or c == "^")
        elif c.isalpha():
            _emit(out, c, sub, sup, i=True)
            k += 1
        elif c == "-":
            _emit(out, "−", sub, sup)                # proper minus sign
            k += 1
        elif c == "=" and not (sub or sup):
            _emit(out, " = ", sub, sup)
            k += 1
        elif c == "," and not (sub or sup):
            _emit(out, ", ", sub, sup)
            k += 1
        elif c == " ":
            _emit(out, " ", sub, sup)
            k += 1
        else:
            _emit(out, c, sub, sup)
            k += 1


def runs(s):
    out, parts = [], s.split("$")
    for n, part in enumerate(parts):
        if n % 2:
            _math(part, out)
        else:
            _emit(out, part, False, False)
    return [r for r in out if r["t"]]


if __name__ == "__main__":
    for probe in ["$x^{(a)}_{sir}$", "$k^{(a)}(u)=2$", "Within-hyperedge factor $C$",
                  "$P(\\mathbf{k})\\rightarrow\\Psi,\\,\\psi_a$",
                  "$\\Rightarrow\;S(t),\\,I(t),\\,R(t),\;R(\\infty)$",
                  "$C\;\\geq\;(m-1)\\,T$", "$-\\,\\delta_{ab}$", "$K_{ab}=B_{ab}\\,C_b$"]:
        flat = "".join(("[%s]" if r["sub"] else "{%s}" if r["sup"] else "%s") % r["t"]
                       for r in runs(probe))
        print("%-46s -> %s" % (probe, flat))
