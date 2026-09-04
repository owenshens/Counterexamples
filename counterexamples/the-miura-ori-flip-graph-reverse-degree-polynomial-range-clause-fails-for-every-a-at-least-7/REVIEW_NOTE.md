# Referee note

**The Reverse Degree-Polynomial Range Clause for the $2\times n$ Miura-ori Flip Graph Fails for Every $a\ge 7$**

Files a referee has here: `paper.tex` and `paper.pdf` (the paper), `verify.py` (the program), `verify.output.txt` (the recorded run of that program), and this note. There is nothing else, and nothing outside this folder is needed to read the paper.

## 1. What the paper claims

Let $v_n^d$ be the number of vertices of degree $d$ in the origami flip graph $\mathrm{OFG}(M_{2,n})$ of the $2\times n$ Miura-ori, and put $f_a(n):=v_n^{2n-a}$. Section 1 ("The statement, with its locator") quotes verbatim, from lines 1256--1258 of the LaTeX source of the cited e-print arXiv:2506.19700v2 of Christensen, Hull, O'Neil, Pappano, Ter-Saakov and Yang, a theorem printed there with no proof: $v_n^{2n-a}$ is a polynomial in $n$ of degree $\lfloor a/2\rfloor$ for $n\ge\lceil a/2\rceil+1$, for all $a\in\mathbb N_0$. The paper calls this the *reverse-diagonal theorem* and settles its **range clause** in the negative.

* **Theorem 1** (Section 2): the reverse-diagonal theorem is false at $a=7$. The witness is five cells of the source's own degree-distribution table, $f_7(5),\dots,f_7(9)=12,88,296,680,1288$ (equation (3)), whose fourth forward difference is $4\ne0$ although the window $n=5,\dots,9$ lies inside the asserted range $n\ge5$ and the source's exclusion clause is vacuous at odd $a$.
* **Proposition 2** and **Corollary 3** (Section 3): with $k_a=\lfloor a/2\rfloor$ and $G_a=(1-x)^{k_a+1}\sum_{n\ge1}f_a(n)x^n$, one has $\deg G_{2m}=3m-1$, $\deg G_{2m+1}=3m$, leading coefficient exactly $4(-1)^{\lfloor a/2\rfloor+1}$, and $G_a(1)>0$; consequently $f_a$ agrees with a polynomial of degree exactly $\lfloor a/2\rfloor$ on $[a-1,\infty)$ and not before, the last nonvanishing value of $\Delta^{k_a+1}f_a$ sitting at base $n=a-2$ with value $4(-1)^{\lfloor a/2\rfloor+1}$.
* **Theorem 4** (Section 3): hence the reverse-diagonal theorem is false for **every** $a\ge7$; the failing base indices inside the asserted range lie in $\{s(a),\dots,a-2\}$ with $s(a)=\lfloor a/2\rfloor+2$ (equation (2)), and the largest is $n=a-2$. For $2\le a\le6$, granting the source's own exclusion of the degree-$2$ cell, the statement does hold on $n\ge s(a)$.
* **What survives** (Section 4, with **Corollary 5**): polynomiality and degree exactly $\lfloor a/2\rfloor$ are *proved* for every $a\ge2$ — the source asserts them without proof — and the repaired range clause is $n\ge\max(\lfloor a/2\rfloor+2,\ a-1)$ (equation (9)), which is $6$ at $a=7$. The weaker repair $n\ge\lfloor a/2\rfloor+2$ is explicitly not correct, being $5$ at $a=7$.

The paper also observes in Section 1 that, read literally, the theorem's first admissible $n$ for even $a$ is the very cell the source excludes, so on that reading it already fails at $a=6$; it then sets that reading aside and grants the exclusion clause throughout, so the $a=7$ witness turns on no quantifier quibble. It asserts nothing at the degenerate index $a=1$.

## 2. What the program checks

`verify.output.txt` records **45 checks, all passing**, and closes `VERDICT: ALL 45 CHECKS PASS`. By block, with the claim each supports:

* **5 checks — calibration of the recurrence before it is used.** 39 distinct cells published by the source are reproduced from the four seeds $v_1^2=2$, $v_2^2=4$, $v_2^3=0$, $v_2^4=2$ and the pure-$v$ recurrence, equation (1): the 30 cells of Table 1 ($d=3,5,7,9,11$, $n=2,\dots,9$), the column $v_3^2,\dots,v_3^6=4,4,8,0,2$ (of which $v_3^4=8$ and $v_3^6=2$ need the $v_1^2=2$ anomaly), the $a=6$ diagonal $36,128,292,544,900$, and $v_4^2=4$. Zero mismatches.
* **6 checks — Theorem 1.** The witness values as printed cells; the same five values forced by (1) out of *other* published cells, exactly as the remark following Theorem 1 prints the sums; $\Delta^4=4\ne0$; the difference tower of Section 2 reproduced line for line as $[76,208,384,608]$, $[132,176,224]$, $[44,48]$, $[4]$; the window lying inside $n\ge5$ with $2n-7=2$ insoluble; and the cubic through $(6,88),(7,296),(8,680),(9,1288)$ predicting $8$ at $n=5$ against $f_7(5)=12$.
* **4 checks — the seed rows and the admissible start.** $f_2(n)=4(n-1)$ for $n\le59$ and $f_3(n)=4(n-2)$ for $2\le n\le59$ (consistent with, not proving, the closed forms of Section 3); the least defined non-excluded index is $\lfloor a/2\rfloor+2$ for $a=0,\dots,80$, giving row lengths $8,8,7,7,6,6$ over $n\le9$ as the source prints them; and equation (2) versus the literal bound for $a=0,\dots,80$.
* **3 checks — the ladder and the dictionary.** $G_2,\dots,G_8$ rebuilt from the recursion (6) agree coefficient for coefficient with the polynomials printed in Section 3; $G_4(1),\dots,G_7(1)=8,16,16,48$ are computed from those coefficient lists (the paper states only $G_a(1)>0$); and $[x^n]G_a=(\Delta^{k_a+1}f_a)(n-k_a-1)$, equation (8), holds in 1414 instances over $a=2,\dots,40$, tying the generating function to the recurrence.
* **4 checks — Proposition 2 over $a=2,\dots,80$:** the two degree formulas; $L_a=4(-1)^{\lfloor a/2\rfloor+1}$, so $|L_a|=4$ throughout; $G_a(1)>0$ (minimum $4$); and the $x=1$ specialisation $G_a(1)=2G_{a-2}(1)$ (even $a$), $G_{a-1}(1)+2G_{a-2}(1)$ (odd $a$) for $a=4,\dots,80$.
* **3 checks — Corollary 3:** $\deg G_a-\lfloor a/2\rfloor=a-1$ for $a=2,\dots,80$; $(\Delta^{k_a+1}f_a)(a-2)=4(-1)^{\lfloor a/2\rfloor+1}\ne0$ with vanishing at $n=a-1,a,a+1,a+7$; and $\Delta^{k_a}f_a$ equal to the constant $G_a(1)$ on $n=a-1,\dots,a+24$ for $a=2,\dots,40$.
* **7 checks — the census and the repair, $a=2,\dots,80$:** the statement holds on its whole admissible range for exactly $a=2,3,4,5,6$ and fails for the 74 values $a=7,\dots,80$; $f_0(n)=2$ and $f_1(n)=0$ for $n\le39$; the failing indices form exactly the block $s(a),\dots,a-2$; their number is $\lceil a/2\rceil-3$; the counts at eight sampled $a$ (e.g. $1$ at $a=7$, $37$ at $a=80$); the least correct start is exactly (9), namely $6$ at $a=7$; and the two terms of (9) coincide for $a\ge5$.
* **1 check — an independent census** computed from raw difference tables of the recurrence with no generating function anywhere, agreeing with the $G$-based census for all 41 values $a=0,\dots,40$.
* **2 checks — shortest in-range windows.** For $a=7,\dots,16$ the final difference on $n=s(a),\dots,s(a)+k_a+1$ is $4,-4,36,-44,172,-260,488,-996,344,-2164$, all nonzero and each equal to the corresponding coefficient of $G_a$; and first obstruction equals last exactly at $a=7,8$, the only $a$ with a single failing $n$.
* **3 checks — checksums on published integers alone:** $1288-4\cdot680+6\cdot296-4\cdot88+12=4=[x^9]G_7$, recovering the witness with no recurrence at all; the same combination on the published $a=6$ diagonal giving $0$; and $v_{10}^{13}=2168$ by two independent routes (the recurrence, and expanding $G_7/(1-x)^4$).
* **5 controls, both polarities:** $a=6$ must stay silent and does, $\Delta^4\equiv0$ on $n=5,\dots,60$; $a=6$ must fire once the excluded cell $v_4^2=4$ is prepended, and does, at $n=4=a-2$ where Corollary 3 puts it; $a=4$ against the wrong hypothesis "degree 1" fires, and against its correct degree 2 stays silent; and corrupting the witness $12\to8$ makes $\Delta^4$ vanish, so the detector reads the data rather than always reporting failure.
* **2 checks — the companion null** on the source's *proved* sibling theorem (fixed degree $d$): corroborated at $m=2$ on the sampled window $n=\lceil d/2\rceil+1,\dots,85$ with degree exactly $d-2$ for all 39 values $d=2,\dots,40$, zero obstructions; and the same test does fire one step below that range at $d=8$, so the null is not vacuous.

## 3. What the program does **not** check

**The paper's universal results are hand proofs and the program is a control.** Proposition 2, Corollary 3 and the falsity assertion of Theorem 4 are proved for all $a$ in Section 3; computation only corroborates them. Theorem 1 needs no computer at all — five printed integers and ten subtractions. The run states its own scope in a closing `SCOPE` block; carried over faithfully:

1. **The source's recurrence proposition is not re-proved.** Equation (1) is taken as given, and the conventions around it ($v_m^e=0$ for $e<2$ or $e>2m$, $v_1^2=2$, $v_m^2=4$ for $m\ge2$) are pinned against the source's own printed integers only. Section 4 says the same, and adds that if that proposition failed at a boundary not exercised by the cells of Table 1, the reverse-diagonal recurrence (4) would need redoing.
2. **Flat-foldability itself is NOT RE-RUN.** No mountain--valley assignment of $M_{2,n}$ is enumerated anywhere, so the model behind (1) is trusted. Theorem 1 does not depend on it, its five values being printed cells; every value at $n\ge10$ does.
3. **The census stops at $a=80$ and $n=130$.** The exact failing set of Theorem 4, and hence attainment of the upper bound $\lceil a/2\rceil-3$ on the number of failures, rests on that census alone and is not claimed beyond $a=80$. The paper asserts neither: Theorem 4 says only that the failing indices lie in $\{s(a),\dots,a-2\}$ with $n=a-2$ the largest, and Section 4 states that which other indices of that set fail "is not asserted here". The counts and the ten first-obstruction values above are quantities the program computes and the paper does not print.
4. **The line numbers, byte count and md5 quoted for the cited e-print are not checked**: they are properties of an external file the program does not fetch, and no such file is in this folder.

Two further limits, both consistent with the paper's own Section 5. First, several blocks are **sampled, not quantified**: the dictionary (8) over $a\le40$, the constant-difference check over $a\le40$ and $n\le a+24$, the independent census over $a\le40$, and the sibling-theorem null over $d\le40$, $n\le85$ — none of these is a proof for all $n$ or all $a$. Second, the witness values and every other cell of Table 1 are **transcribed from the cited source**, not recomputed from the flip graph; the program checks that they are mutually consistent under (1), which is the paper's stated ground for trusting them, but nothing here checks the source's printed table against the graph it describes. The bibliographic comparison with the two cited Gupta e-prints in Section 4 is likewise not verified by any computation.

## 4. How to check it

```sh
python3 verify.py            # Python 3.9+, standard library only, exact integer arithmetic
shasum -a 256 verify.py
```

The program prints one `PASS` line per check, then the `SCOPE` block, then the verdict, and exits 0 only if every check passes. The header of `verify.output.txt` carries the SHA-256 of the program beside its output, so transcript and program can be paired. For the `verify.py` shipped here that digest is

    3bdd93f279c976a61b30fe898876887a46e96cdefe5b25dc6c80d6cc50678386

and it is the digest recorded in the transcript header; the recorded run was made under Python 3.9.25 and ended with status 0. Theorem 1, Proposition 2, Corollary 3, Theorem 4 and the repair (9) can all be checked by hand from the paper alone, with $G_2,\dots,G_8$ printed out in Section 3.
