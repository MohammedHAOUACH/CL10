"""One-pass fix of the manuscript's section numbering and reference list.

Sections: the manuscript jumps 14 -> 16; the hole is closed by renumbering
16..20 down by one, and every internal cross-reference follows.

References: 't Hooft 1976 (a pseudoparticle/instanton paper) sat inside the
causal-set / lattice-gauge / QCA citation cluster and was never cited again;
it is replaced by Wilson's lattice gauge theory, which the sentence names but
never cited. Two references the text uses without citing are added: Lieb and
Robinson (the bound the Wilson factor's tail is measured against) and Furey
(Clifford-algebra unification, the cluster the text claims exists). Hisano's
initial is corrected. All citation markers in the body are remapped through the
two insertions, so no marker is left pointing at the wrong entry.
"""
import re
import sys

PATH = "doc/final_publication.md"

SEC_MAP = {16: 15, 17: 16, 18: 17, 19: 18, 20: 19}

# old number -> new number: 1-5 unchanged, 6-12 shift by one (Lieb-Robinson
# inserted after 5), 13-53 shift by two (Furey inserted after old 12)
def ref_new(n: int) -> int:
    if n <= 5:
        return n
    if n <= 12:
        return n + 1
    return n + 2


WILSON = ('K. G. Wilson, "Confinement of quarks," *Phys. Rev. D* **10**, '
          '2445 (1974).')
LIEB = ('E. H. Lieb and D. W. Robinson, "The finite group velocity of quantum '
        'spin systems," *Commun. Math. Phys.* **28**, 251 (1972).')
FUREY = ('N. Furey, "Generations: Three Prints, in Colour," *J. High Energy '
         'Phys.* **10**, 046 (2014), arXiv:1405.4601.')

CIT = re.compile(r"\[(\d{1,2}(?:[,\u2013]\s?\d{1,2})*)\]")


def remap_citations(text: str) -> str:
    def sub(m):
        parts = re.split(r"([,\u2013])", m.group(1))
        out = [str(ref_new(int(p))) if p.strip().isdigit() else p
               for p in parts]
        return "[" + "".join(out) + "]"
    return CIT.sub(sub, text)


def main():
    src = open(PATH, encoding="utf-8").read()
    head, _, refs_block = src.partition("\n## References\n")
    if not refs_block:
        sys.exit("no References section")

    # ---- sections ----------------------------------------------------------
    head = re.sub(r"(?m)^(## )(\d{1,2})\. ",
                  lambda m: m.group(1)
                  + str(SEC_MAP.get(int(m.group(2)), int(m.group(2))))
                  + ". ", head)

    def sec_sub(m):
        nums = [n for n in (m.group(2),
                            (m.group(3) or "").lstrip("\u2013")) if n]
        return m.group(1) + "\u2013".join(
            str(SEC_MAP.get(int(n), int(n))) for n in nums)

    head = re.sub(r"(Section[s]? )(\d{1,2})(\u2013\d{1,2})?(?![.\d])",
                  sec_sub, head)

    # ---- citation markers --------------------------------------------------
    head = remap_citations(head)

    # the Clifford cluster gains Furey; the Lieb-Robinson bound gains its cite
    head = head.replace("Clifford-algebra formulations of this unification "
                        "exist [11\u201313].",
                        "Clifford-algebra formulations of this unification "
                        "exist [11\u201314].")
    head = head.replace("Lieb\u2013Robinson bound.",
                        "Lieb\u2013Robinson bound [6].")

    # ---- reference list ----------------------------------------------------
    entries, num = {}, None
    for line in refs_block.splitlines():
        m = re.match(r"^(\d{1,2})\. (.*)$", line)
        if m:
            num = int(m.group(1))
            entries[num] = m.group(2)
        elif num is not None and line.strip():
            entries[num] += "\n" + line
    if len(entries) != 53:
        sys.exit(f"expected 53 references, found {len(entries)}")

    entries[2] = WILSON
    entries[46] = entries[46].replace("H. Hisano", "Y. Hisano")

    new = {ref_new(n): text for n, text in entries.items()}
    new[6] = LIEB
    new[14] = FUREY

    out = "\n".join(f"{n}. {new[n]}" for n in sorted(new))
    open(PATH, "w", encoding="utf-8").write(head + "\n## References\n\n" + out
                                            + "\n")
    print(f"{len(new)} references; sections renumbered; markers remapped")


if __name__ == "__main__":
    main()
