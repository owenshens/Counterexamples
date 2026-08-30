#!/usr/bin/env python3
"""Verification program for the note

    Fröberg's conjecture for six to twenty-one general septics in four
    variables, in every characteristic

It reads the objects PRINTED IN THE PAPER -- the 120 hexadecimal masks that
encode the septic family (Appendix A), the ten endpoint cells of Table 1, and
for each cell the two index sets and the two integer determinants of Appendix B
-- and re-derives every quantity the paper claims.  Nothing is read from a
file, from a network, or from any package outside the standard library.

ARITHMETIC.  Exact integers only; no floating point enters any decision.  Each
determinant is recomputed by Gaussian elimination over prime fields F_q with q
a 62-bit prime, and lifted by the Chinese remainder theorem against the
INTEGER Hadamard bound  det(N)^2 <= prod_i (row sum of N)  -- so the modulus
used is provably larger than 2|det(N)| and the lift is the true integer.
p-uniformity is then one Euclid gcd on the recomputed integers.

USAGE:  python3 verify.py          (Python 3.9+, standard library only)
It prints one `PASS <name>` line per check, closes with
`VERDICT: ALL <n> CHECKS PASS`, and exits 0 if and only if every check passed.
Multiprocessing is used only to spread independent modular determinants over
cores; the result does not depend on the number of workers.
"""
import hashlib
import itertools
import json
import os
import sys
from math import comb, gcd

# ===========================================================================
# 1.  THE PRINTED OBJECT.  Every literal below is transcribed from the paper.
# ===========================================================================

# Appendix A: the 120 masks, one per form G_1 ... G_120, in the paper's order.
# Mask bit i (the coefficient of 2**i) is 1 iff the i-th degree-7 monomial of
# the ordered basis MON7 occurs in that form, with coefficient 1.
MASKS = """
49ed660e53a564e3fb9f22d51e14f7
2728e9577f9976d66df364ad739e9b
e4642c6e3d0c00c837aad525608035
a697391e0adb52a5bbbdb4606db8a6
476783c09aedb8c07d24678ea4bbdd
b29ddfdd727b6b736069bafd78deff
53c2fe37a5bc1a5e7c15ddfefcc5a3
8bc155a7fb2fdaca14f9b5219eeaf4
6703cbffc1ee6058f28ea9114ed06c
4c1f42a4ceb4857b5f07e66cedd8ef
4cc658bc60373a0e6d46cc4493e1c5
f380fa591a18fe94d6e091acd9979a
d654523b92b9bb4034ce8b4d293176
44e3a772aea91b2853a289bd834cc9
2b2671f51fac829ce9e991ce47bc38
b53b6f0b4b896b8cbf7f55fb902070
4495e190e5533494c2071c137aa59d
73837221ac7cb93029c29f2e30e086
26ed68177571506497c38afb7b43fd
42be9cc2a6dbe2f55faf1411447aa5
9fe550885e47ed36507d97205539e0
42eb5c1e4fa8ff83d9db98705b66c2
fbf89b1f77618f80c38c08296179bd
a177a6643c7f8c9f7474c3a79ecbfe
f58f2ebd29fc52eca2ebdac6a3e9ef
2cde9eaa1d9f1f82cbdfffce097c4f
05c84b5a29b56bd2977b6fd9fc4aad
3098069e83fad51184a85a5f07cd31
87283e22a21ad20661f3289189cdc4
b7959fecb6554360c88e3cfc2c6920
3ff9eb2b4a08340a55ce387315cc02
d9713ec249404a2f883752c56425bb
7d506963bfc00fbf9525464aec0e74
0c6de3526d9e57467f794f193554f3
67b2171d4bad1c5d41451e91f17589
a20e5f0a6cf168827f70f767abeec3
0e4a0197b32f496caf2dcaf5f988a9
6296a7b9c2daa9baf1611e3405fa15
1822f05badde1804ed3d885041ccae
3cac4abfdb9959c54f1eb5bac3bd3a
e0d7c561a1cd6e1ab979da94b78e55
42248d43672d7f7d1cfd5f86479927
29437027930fac690b0ddff0817779
33ae3a0b813956c3f2be5582f63166
237481e0eca47cf9099873dace054f
06abfeed45fb34659cddd191b573d2
e20a975a3efd158ea598a8aeafec9e
c7d6b502b5add63dff46c76fe4aae3
aa1e98782c3abf2c5d679f7f3b5ba6
8dff2b95b6f4ff038c791259b47caa
0190edb1420f7c0be0778a877820a1
d2fb52dbe7152d879e66e2782adff8
b7f570deb9b885117b8e555c8c4b75
ead768f39752398ac77d226bd2f1cc
80c5bbc4de2788bf17512482811ffc
096e43d32719fd327468d76cc56e9d
57ab35b5adaaad0b266fdcac04ce6f
9a8515832811b4f51a81137ad2d342
efc269a49930bd500bb5f1408c5027
a9fad2a2f0da2ff8537809ef754744
5ff5ce6e86ebed38490d0ec869e61c
744cb2bc0684e8959414548d245637
c8a942c0a456f87cc6933dcd22f5eb
c9b02e17dc8b0a4b7733ed48f306dc
a19e90b01a02716c06d6c586b6fa60
7084ab304f1a6af6b6598707188832
8dab2723c2301d447c7f3a3bf9940b
5320fcfad95f9bbe544fe129574d25
a311c483ea0f7e218738e234407abc
b681da1e79477e3fef2f13ce5f48c4
9b058909d5210c09b363f049342f88
c6da4a7e4dc439ec9459f00991979d
48f76c9fd812b41cd78a667e201c68
ca993b225ab6661355075d308d44ba
2b93d2aba48494d4b4f9ff48de79f8
f67b65c286bd6e03f38cbfbea0434b
4329e488e967c5261ddc5e21fcf51e
8ac189f905c00834b45b2ff4c4525c
02237aa1dfd2cc9f49c468c63dda9e
0ae06f10cb535d7d06596dd00046ce
f579ec3fbcd61e3ab2375bf248a769
fb94ff952a7a071eb10970be62ef3d
ae52183007e685321a648f9d9a875d
51f07211585b109b62ca32d6ec31dc
5c9a931b74a63a49b394bf6482ddbc
9acbf7f07f33887c781e12d74552af
789951550f42119045f95cb906cc9f
9f72b79a28ada4f9350d710931788c
3a17668c78bfe187db8cea367d29d7
82b792dbcd843fb34e06704246139d
408a373dff114f8aee60a558a37ea5
afddf11916c9d402f4ae3299436103
83a0a9b1dc178966744971e25866c3
8c700a53829e4842930b0dc5cc385e
fedd1345c3bedd7a1c99ae2eca2c26
c1a09fb74c69d83cd28e8b49aa70a7
0f0665a29a4194fe28e393e5cd10f3
c7be15ff006480f8f0ac913d385947
9a01f5f40ea8a73fb5b1289c8567c3
8045b85e2ce9643f7d04a0b4c54ead
5d6f2eb1d88150183697b6068db163
ed8e55714db318fd7b5057c75e421a
6cbfe1bc1f2c05201e1ec0e782fee7
2f29eeeca3726cbd0174f283f74f3d
70d8b2ea8ff57e89f79c5865ee7fc8
3ea5b004e38845ae1709cf527c0788
f9c8d13a56c17d14a5dba606f8fdbd
323eef41988b37c9683381ab7284ee
5642a38d58c97cbe17492968f38e9a
98ff8f00dbe7ca5a2a95c9be50ef9b
6d4088d2c44ddec1ee4ac769643c12
17e558fb16d212d5f1cc8886d349ad
d46aab8deebf1cacfcff6151575326
f77761b0778ff9143607dde73b0c71
2afe8dd2da7fe16ba1b21908aba085
e9250b90c1b657104fcdcf7a0346df
8229a6e6980bccbc3f0170b2dc8a95
308dfa6037908a190896ea4bb05c46
0d8e89826c292581e9727daf07e4b8
1eb47eaf2eeefd0e5d660864006d9a
"""

# Section 3: the sha256 the paper prints for the canonical JSON encoding of the
# whole family, which is also the digest of the `forms` array of the ancillary
# certificate published with arXiv:2608.24797.
FAMILY_SHA256 = "f93339e7b30deb19b81fe74fa3b79a7859cf8da9a4c6fadbd740f98c975c8e20"

# Section 3: the 65 exponent vectors of G_1, printed in the paper so that the
# mask convention can be checked by hand.
G1_EXPONENTS = """
[7,0,0,0] [6,1,0,0] [6,0,1,0] [5,2,0,0] [5,1,1,0] [5,1,0,1] [5,0,2,0] [4,3,0,0] [4,2,0,1] [4,0,2,1]
[4,0,1,2] [4,0,0,3] [3,4,0,0] [3,2,1,1] [3,1,3,0] [3,1,1,2] [3,0,4,0] [3,0,3,1] [3,0,1,3] [2,4,0,1]
[2,3,0,2] [2,2,3,0] [2,2,2,1] [2,2,1,2] [2,2,0,3] [2,1,2,2] [2,1,1,3] [2,1,0,4] [2,0,4,1] [2,0,3,2]
[2,0,2,3] [2,0,1,4] [2,0,0,5] [1,6,0,0] [1,5,1,0] [1,4,0,2] [1,3,3,0] [1,3,2,1] [1,2,4,0] [1,2,1,3]
[1,2,0,4] [1,1,4,1] [1,1,2,3] [1,0,6,0] [1,0,4,2] [1,0,3,3] [1,0,2,4] [0,7,0,0] [0,6,0,1] [0,5,0,2]
[0,4,3,0] [0,4,2,1] [0,3,1,3] [0,3,0,4] [0,2,3,2] [0,2,2,3] [0,2,0,5] [0,1,5,1] [0,1,4,2] [0,1,2,4]
[0,1,1,5] [0,1,0,6] [0,0,7,0] [0,0,4,3] [0,0,1,6]
"""

# Table 1 and Appendix B: the ten endpoint cells.  For each cell:
#   r, j              the generator count and the degree
#   role              'inj'  mu_{r,j} is injective  (columns <= rows, s = cols)
#                     'surj' mu_{r,j} is surjective (rows < columns, s = rows)
#   minors            two entries (drop, det):
#                       drop -- the indices DELETED from the long side, i.e.
#                               the rows omitted when role == 'inj' and the
#                               columns omitted when role == 'surj'
#                       det  -- the integer determinant of the resulting
#                               s x s minor, as printed in Appendix B
CELLS = [
    {"r": 6, "j": 13, "R": 560, "C": 504, "s": 504, "role": "inj", "minors": [
        {"drop": (
            504, 505, 506, 507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518,
            519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533,
            534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548,
            549, 550, 551, 552, 553, 554, 555, 556, 557, 558, 559
        ), "det": (
            "377606923044870977014099594084102554950277533115555995861448"
            "668128353556875304003247184617529678146942732232494074692376"
            "172793828620781900013583046963559699624693136040131276903721"
            "513841578330296735509841891277769652294577388534912106160729"
            "8"
        )},
        {"drop": (
            24, 28, 36, 57, 86, 93, 139, 162, 165, 168, 173, 177, 180, 181, 186, 217,
            236, 241, 255, 256, 257, 263, 274, 285, 315, 326, 332, 362, 368, 369, 370,
            372, 380, 389, 402, 407, 409, 424, 425, 433, 441, 455, 456, 472, 473, 476,
            501, 513, 514, 521, 522, 525, 526, 534, 540, 551
        ), "det": (
            "250908512247727421748306528448406074588019548478375995555645"
            "403523878444321315034625730994534344573962869288878217168938"
            "857971609565567933574038571188008591177148408393232072793363"
            "6590671502391334022421488708480139429457764165007894273479"
        )},
    ]},
    {"r": 6, "j": 14, "R": 680, "C": 720, "s": 680, "role": "surj", "minors": [
        {"drop": (
            3, 69, 115, 174, 202, 203, 208, 253, 284, 297, 312, 344, 368, 386, 391,
            398, 414, 418, 436, 460, 472, 504, 523, 529, 538, 543, 570, 593, 598, 612,
            615, 634, 636, 644, 650, 655, 663, 664, 684, 713
        ), "det": (
            "750107683604837441246021253344764986440676038346602463121146"
            "940988166494592584704381613937489073035856580026863609523206"
            "495811171063585434530680564292542007917348235641136559917168"
            "877893135402064252260427569972199250519278152116585151390628"
            "604093788731990745096132410259890089288754060779013192127197"
            "2893647721523572158390"
        )},
        {"drop": (
            9, 25, 47, 79, 94, 108, 110, 121, 131, 150, 173, 177, 179, 225, 259, 313,
            344, 359, 395, 402, 410, 421, 437, 462, 463, 465, 468, 476, 507, 512, 522,
            537, 551, 602, 605, 664, 680, 684, 696, 708
        ), "det": (
            "-84437176153676729778432735134877212269136086730820661876732"
            "531539472005064192699818963971876169325526097406674928042147"
            "224856110892480543399457199439977536339864593828699331037745"
            "759557459668456281767317598213430672689636402240110193828288"
            "959656027800315168822094898847748099594280822610970565916739"
            "818070374811479909403131"
        )},
    ]},
    {"r": 8, "j": 12, "R": 455, "C": 448, "s": 448, "role": "inj", "minors": [
        {"drop": (
            448, 449, 450, 451, 452, 453, 454
        ), "det": (
            "170893115052404981084168463687847203813409321064960621241614"
            "200362033798803012631407834422074141545566901466550738433726"
            "369419946059220054803745156375399859457867337059280775654803"
            "900548631977548845680909226309933404569948052"
        )},
        {"drop": (
            28, 43, 46, 184, 427, 434, 441
        ), "det": (
            "-58335925049442448880162209191873498253681633577569600683434"
            "345602169078185092860340525902940525951508017935508045562930"
            "418042428320286417428198840632611135832517474437769071560302"
            "3411042266142781736713526406033984851229795"
        )},
    ]},
    {"r": 7, "j": 13, "R": 560, "C": 588, "s": 560, "role": "surj", "minors": [
        {"drop": (
            560, 561, 562, 563, 564, 565, 566, 567, 568, 569, 570, 571, 572, 573, 574,
            575, 576, 577, 578, 579, 580, 581, 582, 583, 584, 585, 586, 587
        ), "det": (
            "146964866351168454024042982154176178171178035326101215912707"
            "059360420316057797836302646965272345653909068572210898931934"
            "833436481779969531975403172024283871019484047201853166556071"
            "566371998781822606797701868246992190163677955947917663845424"
            "866371383953206462937214895232874"
        )},
        {"drop": (
            1, 27, 34, 97, 108, 133, 135, 143, 156, 167, 172, 189, 240, 260, 293, 296,
            309, 311, 329, 380, 381, 382, 385, 420, 476, 480, 559, 566
        ), "det": (
            "-27886393297332648683244058435250420139573881962734593770137"
            "579165267006287005209559069635569189512026376344687512006128"
            "533683844461017076209051170134739574546457463072086211064446"
            "491040297541092872659736849414665641641061329807406838020236"
            "075180833064744039451492042150163213"
        )},
    ]},
    {"r": 10, "j": 11, "R": 364, "C": 350, "s": 350, "role": "inj", "minors": [
        {"drop": (
            350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363
        ), "det": (
            "155806231330096120162956424800025454628510522334153400538458"
            "792058944540491009799960045162876613606432073945287142050223"
            "01506699980011029136035981123963221197369317966909714135924"
        )},
        {"drop": (
            18, 28, 43, 46, 86, 108, 128, 157, 184, 297, 310, 342, 348, 355
        ), "det": (
            "-31041880753344543373977065103650942680682186897305957058721"
            "814255657388889531749988875803277326538961122392227279521130"
            "17983547649403376244894719234582851258198898184980577"
        )},
    ]},
    {"r": 9, "j": 12, "R": 455, "C": 504, "s": 455, "role": "surj", "minors": [
        {"drop": (
            455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467, 468, 469,
            470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480, 481, 482, 483, 484,
            485, 486, 487, 488, 489, 490, 491, 492, 493, 494, 495, 496, 497, 498, 499,
            500, 501, 502, 503
        ), "det": (
            "-35432557366446796463147306228688806233900183935013832902697"
            "174906177691535180050912937946354957737388090999732898092801"
            "903566955288198245769362063696420884611722347724522686034458"
            "3790864937420883957233045304918134290418243697847"
        )},
        {"drop": (
            5, 10, 11, 36, 37, 38, 39, 51, 54, 64, 72, 76, 78, 88, 91, 105, 124, 137,
            148, 150, 156, 163, 167, 183, 205, 240, 252, 256, 260, 263, 282, 308, 325,
            340, 363, 373, 384, 397, 398, 405, 407, 413, 420, 427, 439, 441, 443, 444,
            468
        ), "det": (
            "236599598428142818718241895901209701829756482117190784521262"
            "076772073535774475554536078866191198966909780642153596896496"
            "310870211188632758119733452001604352988142427515558992907351"
            "09266466870280191446114061999547306338078650544"
        )},
    ]},
    {"r": 14, "j": 10, "R": 286, "C": 280, "s": 280, "role": "inj", "minors": [
        {"drop": (
            280, 281, 282, 283, 284, 285
        ), "det": (
            "-74247568238809266305143535716716522813088514354455908889320"
            "125985628714546244367645789588825361581913714919520392971752"
            "849111798074701147635457"
        )},
        {"drop": (
            32, 60, 68, 130, 230, 253
        ), "det": (
            "-21370524956637011610587381398778328146173302700922739598486"
            "102161540633377285471693396758142521585392993937183773210684"
            "7810837382437199312064"
        )},
    ]},
    {"r": 11, "j": 11, "R": 364, "C": 385, "s": 364, "role": "surj", "minors": [
        {"drop": (
            5, 13, 20, 32, 57, 110, 128, 151, 160, 162, 163, 200, 204, 230, 233, 277,
            305, 307, 316, 371, 373
        ), "det": (
            "-55397627070351938046274250294661915769843458364946362073211"
            "727342414415894747560701374066881573783936381424212723498048"
            "085061203403749714849953922598659880491767054456041961681439"
            "20607"
        )},
        {"drop": (
            49, 57, 79, 80, 86, 88, 120, 127, 131, 146, 168, 206, 216, 224, 236, 263,
            264, 268, 276, 360, 366
        ), "det": (
            "-30696742882746872617268149762132270354727453498616865941836"
            "680107388670371725601769929836559266113148561337103980065836"
            "163387995813573765606588599564008832090589946325105710581530"
            "585548"
        )},
    ]},
    {"r": 21, "j": 9, "R": 220, "C": 210, "s": 210, "role": "inj", "minors": [
        {"drop": (
            210, 211, 212, 213, 214, 215, 216, 217, 218, 219
        ), "det": (
            "-40751531588623599923042140006758177622096513222763483405146"
            "406643019504940454347329590149204392770098211936412"
        )},
        {"drop": (
            16, 30, 34, 65, 126, 145, 194, 195, 205, 216
        ), "det": (
            "839029320151226341121549964244775137656977623915937147094473"
            "5387188784007205371017529405486436026303935194157"
        )},
    ]},
    {"r": 15, "j": 10, "R": 286, "C": 300, "s": 286, "role": "surj", "minors": [
        {"drop": (
            286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299
        ), "det": (
            "375993434779847374821798809151544895944562139603035559554606"
            "635629494745074998770728276567915483045035198175230845978622"
            "66232853159935637991264001"
        )},
        {"drop": (
            20, 33, 43, 49, 68, 86, 120, 136, 176, 194, 248, 272, 273, 280
        ), "det": (
            "619236910386292512028517815743471647074240479149032839152131"
            "148663232953238230964627481128967523432873328822771979156933"
            "763127256206254867848212047"
        )},
    ]},
]

D = 7            # the common generator degree (septics)
NVARS = 4        # the number of variables
RLO, RHI = 6, 21 # the closed interval of generator counts this note settles

# ===========================================================================
# 2.  CHECK BOOKKEEPING
# ===========================================================================
_state = {"pass": 0, "fail": 0}


def ck(name, ok, detail=""):
    """Record one check.  `ok` must already be a bool decided by exact arithmetic."""
    if ok:
        _state["pass"] += 1
        print("PASS %s%s" % (name, (" [%s]" % detail) if detail else ""), flush=True)
    else:
        _state["fail"] += 1
        print("FAIL %s%s" % (name, (" [%s]" % detail) if detail else ""), flush=True)


def note(text):
    print("NOTE %s" % text, flush=True)


def head(text):
    print("\n=== %s" % text, flush=True)


# ===========================================================================
# 3.  EXACT ARITHMETIC ENGINES
# ===========================================================================
_MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def is_prime(n):
    """Deterministic for n < 3.3e24 with these 12 bases (Sorenson-Webster), and
    every modulus used here is below 2**62 < 4.7e18."""
    if n < 2:
        return False
    for p in _MR_BASES:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in _MR_BASES:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


_P62 = []
_P62_NEXT = [(1 << 62) - 1]


def primes_below_2_62(count):
    """The `count` largest primes below 2**62, descending.  Deterministic, so the
    modulus used for a given minor is a function of the printed data alone."""
    while len(_P62) < count:
        n = _P62_NEXT[0]
        if is_prime(n):
            _P62.append(n)
        _P62_NEXT[0] = n - 2
    return _P62[:count]


def det_mod(rows, n, p):
    """det of an n x n integer matrix modulo the prime p, by elimination."""
    M = [r[:] for r in rows]
    d = 1
    for c in range(n):
        piv = -1
        for i in range(c, n):
            if M[i][c] % p:
                piv = i
                break
        if piv < 0:
            return 0
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            d = -d
        pr = M[c]
        pv = pr[c] % p
        d = d * pv % p
        inv = pow(pv, p - 2, p)
        tail = pr[c + 1:]
        for i in range(c + 1, n):
            ri = M[i]
            f = ri[c] * inv % p
            if f:
                M[i] = ri[:c + 1] + [(a - f * b) % p for a, b in zip(ri[c + 1:], tail)]
    return d % p


def rank_mod(rows, nr, nc, p):
    """rank over F_p of an nr x nc integer matrix, by elimination."""
    M = [r[:] for r in rows]
    row = 0
    for c in range(nc):
        piv = -1
        for i in range(row, nr):
            if M[i][c] % p:
                piv = i
                break
        if piv < 0:
            continue
        if piv != row:
            M[row], M[piv] = M[piv], M[row]
        pr = M[row]
        inv = pow(pr[c] % p, p - 2, p)
        tail = pr[c + 1:]
        for i in range(row + 1, nr):
            ri = M[i]
            f = ri[c] * inv % p
            if f:
                M[i] = ri[:c + 1] + [(a - f * b) % p for a, b in zip(ri[c + 1:], tail)]
        row += 1
        if row == nr:
            break
    return row


def det_bareiss(rows):
    """Exact integer determinant by fraction-free (Bareiss) elimination -- a second,
    completely independent engine, used on small matrices to cross-check the
    CRT engine."""
    M = [r[:] for r in rows]
    n = len(M)
    sign = 1
    prev = 1
    for c in range(n - 1):
        if M[c][c] == 0:
            piv = -1
            for i in range(c + 1, n):
                if M[i][c] != 0:
                    piv = i
                    break
            if piv < 0:
                return 0
            M[c], M[piv] = M[piv], M[c]
            sign = -sign
        for i in range(c + 1, n):
            for k in range(c + 1, n):
                M[i][k] = (M[i][k] * M[c][c] - M[i][c] * M[c][k]) // prev
            M[i][c] = 0
        prev = M[c][c]
    return sign * M[n - 1][n - 1]


def hadamard_square(rows):
    """prod_i ||row_i||^2 = prod_i (row sum), valid because every entry is 0 or 1;
    an integer upper bound for det^2 (Hadamard).  Returned as an exact integer."""
    P = 1
    for r in rows:
        P *= sum(r)
    return P


def crt_lift(residues, primes):
    """The unique X with X = residues[i] mod primes[i] and |X| <= prod/2."""
    X, P = 0, 1
    for a, p in zip(residues, primes):
        t = ((a - X) % p) * pow(P % p, p - 2, p) % p
        X += P * t
        P *= p
    if X > P // 2:
        X -= P
    return X


def crt_primes_for(rows):
    """The deterministic prime list for one minor: the largest primes below 2**62,
    as many as make (prod p)^2 > 4 * (Hadamard bound)^2, so the symmetric lift is
    provably the true integer determinant."""
    bound_sq = hadamard_square(rows)          # >= det^2
    need = 4 * bound_sq                       # (2|det|)^2 <= need
    ps, M = [], 1
    cand = primes_below_2_62(64)
    for p in cand:
        ps.append(p)
        M *= p
        if M * M > need:
            return ps, bound_sq
    raise AssertionError("64 primes below 2**62 did not cover the Hadamard bound")


# ===========================================================================
# 4.  THE POLYNOMIAL BOOKKEEPING
# ===========================================================================
def monomials(degree):
    """The exponent 4-tuples of the given degree, in DECREASING lexicographic
    order.  This is the paper's ordered basis of S_degree."""
    v = [x for x in itertools.product(range(degree + 1), repeat=NVARS) if sum(x) == degree]
    v.sort(reverse=True)
    return v


MON = {m: monomials(m) for m in range(0, 18)}
ROWOF = {m: {mo: i for i, mo in enumerate(MON[m])} for m in MON}


def build_matrix(forms, r, j):
    """The matrix of mu_{r,j} : (S_{j-7})^r -> S_j, (h_i) |-> sum h_i G_i, in the
    monomial bases: rows indexed by MON[j], column fi*|MON[j-7]| + mi indexed by
    the mi-th monomial of MON[j-7] multiplying the fi-th form."""
    tg = MON[j]
    rowof = ROWOF[j]
    mul = MON[j - D]
    A = [[0] * (r * len(mul)) for _ in range(len(tg))]
    for fi in range(r):
        for mi, m in enumerate(mul):
            col = fi * len(mul) + mi
            for t in forms[fi]:
                A[rowof[(t[0] + m[0], t[1] + m[1], t[2] + m[2], t[3] + m[3])]][col] += 1
    return A


def take_minor(A, nr, nc, role, drop):
    """The s x s minor named in the paper: for role 'inj' every column and the rows
    outside `drop`; for role 'surj' every row and the columns outside `drop`."""
    ds = set(drop)
    if role == "inj":
        return [A[i][:] for i in range(nr) if i not in ds]
    keep = [c for c in range(nc) if c not in ds]
    return [[row[c] for c in keep] for row in A]


def froberg_coeff(r, j):
    """[z^j] (1-z^7)^r / (1-z)^4 -- the untruncated Froberg numerator coefficient."""
    s = 0
    for i in range(0, r + 1):
        e = j - D * i
        if e < 0:
            break
        s += (-1) ** i * comb(r, i) * comb(e + NVARS - 1, NVARS - 1)
    return s


def froberg_q(r, jmax=64):
    """q(r) = the last degree up to which the numerator coefficients are all > 0.
    The conjecture's truncation rule sets the series to 0 from q(r)+1 on."""
    q = -1
    for j in range(0, jmax + 1):
        if froberg_coeff(r, j) > 0:
            q = j
        else:
            break
    return q


def froberg_series(r, jmax):
    q = froberg_q(r)
    return [froberg_coeff(r, j) if j <= q else 0 for j in range(jmax + 1)]


# ===========================================================================
# 5.  DECODING THE PRINTED OBJECT
# ===========================================================================
def parse_masks():
    return [t for t in MASKS.split()]


def decode_forms(masks):
    out = []
    for h in masks:
        v = int(h, 16)
        out.append([MON[D][i] for i in range(len(MON[D])) if (v >> i) & 1])
    return out


def parse_g1():
    toks = G1_EXPONENTS.replace("[", " ").replace("]", " ").replace(",", " ").split()
    nums = [int(t) for t in toks]
    return [tuple(nums[i:i + 4]) for i in range(0, len(nums), 4)]


# ===========================================================================
# 6.  THE PARALLEL DETERMINANT / RANK JOBS
# ===========================================================================
JOBS = {}     # jid -> ('det', rows, n, p) | ('rank', rows, nr, nc, p) held in the
              # parent and inherited by forked workers


def _run_job(jid):
    kind = JOBS[jid][0]
    if kind == "det":
        _, rows, n, p = JOBS[jid]
        return jid, det_mod(rows, n, p)
    _, rows, nr, nc, p = JOBS[jid]
    return jid, rank_mod(rows, nr, nc, p)


def run_jobs(order):
    """Run every job in JOBS, using forked workers when available.  The mapping
    jid -> value is independent of the worker count."""
    if len(order) == 0:
        return {}
    workers = 1
    ctx = None
    try:
        import multiprocessing as mp
        if hasattr(mp, "get_context") and "fork" in mp.get_all_start_methods():
            ctx = mp.get_context("fork")
            workers = min(30, mp.cpu_count() or 1)
    except Exception:
        ctx = None
    note("%d modular jobs, %d worker process(es)" % (len(order), workers))
    res = {}
    if ctx is None or workers <= 1:
        for jid in order:
            k, v = _run_job(jid)
            res[k] = v
        return res
    with ctx.Pool(workers) as pool:
        for k, v in pool.imap_unordered(_run_job, order, chunksize=1):
            res[k] = v
    return res


# ===========================================================================
# 7.  MAIN
# ===========================================================================
def main():
    print("verification of the note: Froberg's conjecture for %d <= r <= %d general "
          "septics in four variables, in every characteristic" % (RLO, RHI))
    print("python %s, exact integer arithmetic only" % sys.version.split()[0])

    # ---------------------------------------------------------------- step 1
    head("Step 1: the printed family, decoded from its 120 hexadecimal masks")
    masks = parse_masks()
    ck("mask_table_has_120_rows", len(masks) == 120, "%d" % len(masks))
    ck("every_mask_is_30_lowercase_hex_digits",
       all(len(h) == 30 and all(c in "0123456789abcdef" for c in h) for h in masks))
    ck("degree_7_monomial_basis_has_120_elements", len(MON[D]) == comb(D + 3, 3),
       "C(10,3)=%d" % comb(D + 3, 3))
    ck("monomial_basis_is_in_decreasing_lexicographic_order",
       MON[D] == sorted(MON[D], reverse=True) and MON[D][0] == (7, 0, 0, 0)
       and MON[D][-1] == (0, 0, 0, 7))
    forms = decode_forms(masks)
    ck("the_bit_count_of_every_mask_equals_the_support_of_the_form_it_decodes_to",
       all(bin(int(h, 16)).count("1") == len(f) for h, f in zip(masks, forms)))
    ck("every_decoded_form_is_a_sum_of_distinct_degree_7_monomials_with_coefficient_1",
       all(all(sum(t) == D for t in f) and len(set(f)) == len(f) for f in forms))
    sizes = [len(f) for f in forms]
    ck("support_sizes_run_from_49_to_78_as_the_paper_states",
       min(sizes) == 49 and max(sizes) == 78, "min=%d max=%d" % (min(sizes), max(sizes)))
    ck("the_120_forms_are_pairwise_distinct", len(set(map(tuple, forms))) == 120)
    canon = json.dumps([[list(t) for t in f] for f in forms], separators=(",", ":"))
    dig = hashlib.sha256(canon.encode()).hexdigest()
    ck("canonical_json_of_the_decoded_family_has_the_sha256_printed_in_the_paper",
       dig == FAMILY_SHA256, dig)
    g1 = parse_g1()
    ck("G1_decodes_to_the_65_exponent_vectors_printed_in_the_paper",
       forms[0] == g1, "%d vectors" % len(g1))
    ck("the_printed_support_size_table_matches_the_decoded_supports",
       sizes[:RHI] == [len(f) for f in forms[:RHI]], "G1..G21 = %s" % sizes[:RHI])
    note("supports of G1..G21 = %s" % sizes[:RHI])

    # ---------------------------------------------------------------- step 2
    head("Step 2: the truncated Froberg series, and the five blocks of [6,21]")
    qs = {}
    for r in range(RLO, RHI + 1):
        q = froberg_q(r)
        qs[r] = q
        a = [froberg_coeff(r, j) for j in range(0, q + 2)]
        ck("q_of_r_%d_is_%d_first_nonpositive_numerator_coefficient_at_%d" % (r, q, q + 1),
           all(x > 0 for x in a[:q + 1]) and a[q + 1] <= 0,
           "a_%d=%d a_%d=%d" % (q, a[q], q + 1, a[q + 1]))
    note("q(r) for r=6..21 = %s" % [qs[r] for r in range(RLO, RHI + 1)])
    blocks = []
    for r in range(RLO, RHI + 1):
        if blocks and blocks[-1][2] == qs[r]:
            blocks[-1][1] = r
        else:
            blocks.append([r, r, qs[r]])
    blocks = [tuple(b) for b in blocks]
    ck("q_is_constant_on_exactly_five_blocks", len(blocks) == 5, "%s" % (blocks,))
    ck("the_blocks_partition_the_interval_6_to_21_with_no_gap_and_no_overlap",
       blocks[0][0] == RLO and blocks[-1][1] == RHI
       and all(blocks[i][1] + 1 == blocks[i + 1][0] for i in range(len(blocks) - 1)))
    endpoints = set()
    for a_, b_, q_ in blocks:
        endpoints.add((b_, q_, "inj"))
        endpoints.add((a_, q_ + 1, "surj"))
    printed = set((c["r"], c["j"], c["role"]) for c in CELLS)
    ck("the_ten_printed_cells_are_exactly_the_endpoint_cells_of_the_five_blocks",
       endpoints == printed and len(printed) == 10, "%d cells" % len(printed))
    for r in range(RLO, RHI + 1):
        q = qs[r]
        ok = all(froberg_coeff(r, j) == comb(j + 3, 3) - r * comb(j - D + 3, 3)
                 for j in range(D, q + 1))
        ck("for_r_%d_no_second_koszul_term_enters_in_degrees_7_to_%d" % (r, q), ok)
    ck("the_only_endpoint_at_or_above_2d_is_r_6_j_14_and_its_coefficient_is_negative",
       [c for c in CELLS if c["j"] >= 2 * D] == [c for c in CELLS if (c["r"], c["j"]) == (6, 14)]
       and froberg_coeff(6, 14) == comb(17, 3) - 6 * comb(10, 3) + comb(6, 2) * comb(3, 3),
       "a_14(r=6) = 680-720+15 = %d" % froberg_coeff(6, 14))
    ck("every_endpoint_with_r_at_least_7_sits_strictly_below_2d",
       all(c["j"] < 2 * D for c in CELLS if c["r"] >= 7),
       "max j among r>=7 is %d < %d" % (max(c["j"] for c in CELLS if c["r"] >= 7), 2 * D))
    ck("r_21_still_has_positive_series_support_in_degree_10",
       froberg_q(21) == 9 and froberg_coeff(21, 9) > 0 and froberg_coeff(21, 10) < 0,
       "a_9=%d a_10=%d" % (froberg_coeff(21, 9), froberg_coeff(21, 10)))
    ck("r_22_has_a_zero_in_degree_9_so_degree_d_plus_2_methods_reach_it_and_not_r_21",
       froberg_coeff(22, 9) == 0 and froberg_q(22) == 8,
       "a_9(r=22)=0, q(22)=%d" % froberg_q(22))

    # ---------------------------------------------------------------- step 3
    head("Step 3: the ten cells -- shapes, roles and index sets")
    mats = {}
    subs = {}
    for c in CELLS:
        r, j, role = c["r"], c["j"], c["role"]
        key = "r=%d,j=%d" % (r, j)
        R = comb(j + 3, 3)
        C = r * comb(j - D + 3, 3)
        s = min(R, C)
        ck("shape_of_cell_%s_is_%dx%d_from_the_binomials" % (key, R, C),
           (R, C) == (c["R"], c["C"]) and s == c["s"],
           "dim S_%d=%d, %d*dim S_%d=%d, s=%d" % (j, R, r, j - D, C, s))
        ck("role_of_cell_%s_is_%s" % (key, role), (role == "inj") == (C < R))
        A = build_matrix(forms, r, j)
        mats[key] = A
        ck("matrix_of_cell_%s_is_0_1_with_the_stated_shape" % key,
           len(A) == R and len(A[0]) == C and all(x in (0, 1) for row in A for x in row))
        for mi, m in enumerate(c["minors"]):
            drop = m["drop"]
            long_side = R if role == "inj" else C
            ck("index_set_%d_of_cell_%s_is_increasing_in_range_and_of_size_%d"
               % (mi, key, long_side - s),
               list(drop) == sorted(set(drop)) and len(drop) == long_side - s
               and all(0 <= x < long_side for x in drop),
               "%d deleted from the %s side" % (len(drop), "row" if role == "inj" else "column"))
            N = take_minor(A, R, C, role, drop)
            ck("minor_%d_of_cell_%s_is_square_of_size_%d" % (mi, key, s),
               len(N) == s and all(len(row) == s for row in N))
            subs[(key, mi)] = N

    # ---------------------------------------------------------------- step 4
    head("Step 4: scheduling the exact determinants and the small-prime ranks")
    plan = []
    for c in CELLS:
        key = "r=%d,j=%d" % (c["r"], c["j"])
        for mi in range(len(c["minors"])):
            N = subs[(key, mi)]
            ps, bsq = crt_primes_for(N)
            plan.append((key, mi, ps, bsq))
            for p in ps:
                JOBS[("det", key, mi, p)] = ("det", N, len(N), p)
    SMALL = (2, 3, 5)
    for c in CELLS:
        key = "r=%d,j=%d" % (c["r"], c["j"])
        for p in SMALL:
            JOBS[("rank", key, p)] = ("rank", mats[key], c["R"], c["C"], p)
    # engine controls, scheduled in the same pool
    for n in (12, 30):
        JM = [[1 if a != b else 0 for b in range(n)] for a in range(n)]
        ps, _ = crt_primes_for(JM)
        note("control J-I n=%d uses %d prime(s)" % (n, len(ps)))
        for p in ps:
            JOBS[("ctrlJ", n, p)] = ("det", JM, n, p)
    # forced-positive control: r=4, F_i = x_i^7, degree 13, the 336 columns are
    # distinct monomials, so the minor on the rows they occupy is a permutation
    # matrix and its determinant must be +-1.
    ci_forms = [[(7, 0, 0, 0)], [(0, 7, 0, 0)], [(0, 0, 7, 0)], [(0, 0, 0, 7)]]
    A_ci = build_matrix(ci_forms, 4, 13)
    rows_ci = [i for i in range(len(A_ci)) if any(A_ci[i])]
    N_ci = [A_ci[i][:] for i in rows_ci]
    ps_ci, _ = crt_primes_for(N_ci)
    for p in ps_ci:
        JOBS[("ctrlCI", p)] = ("det", N_ci, len(N_ci), p)
    # proved-deficient control: r=6 with two repeated generators, degree 13.
    def_forms = [[(7, 0, 0, 0)], [(0, 7, 0, 0)], [(0, 0, 7, 0)], [(0, 0, 0, 7)],
                 [(7, 0, 0, 0)], [(0, 7, 0, 0)]]
    A_def = build_matrix(def_forms, 6, 13)
    for p in (2, 1000003):
        JOBS[("ctrlDEF", p)] = ("rank", A_def, len(A_def), len(A_def[0]), p)
    # tamper control on the object: flip the lowest set bit of the mask of G_21
    # and re-run the (21,9) determinant modulo one prime.
    tampered = list(forms)
    v = int(masks[20], 16)
    low = v & -v
    tampered[20] = [MON[D][i] for i in range(120) if ((v ^ low) >> i) & 1]
    A_tam = build_matrix(tampered, 21, 9)
    c219 = [c for c in CELLS if (c["r"], c["j"]) == (21, 9)][0]
    N_tam = take_minor(A_tam, c219["R"], c219["C"], c219["role"], c219["minors"][0]["drop"])
    P_TAM = primes_below_2_62(1)[0]
    JOBS[("ctrlTAM", 0)] = ("det", N_tam, len(N_tam), P_TAM)
    JOBS[("ctrlTAM", 1)] = ("det", subs[("r=21,j=9", 0)], c219["s"], P_TAM)

    order = sorted(JOBS.keys(), key=lambda k: repr(k))
    ck("every_scheduled_job_is_a_determinant_or_a_rank_over_a_prime_field",
       all(JOBS[k][0] in ("det", "rank") for k in order), "%d jobs" % len(order))
    ck("the_crt_modulus_of_every_minor_provably_exceeds_twice_its_hadamard_bound",
       all(_modprod(ps) ** 2 > 4 * bsq for _k, _m, ps, bsq in plan),
       "%d minors, %d..%d primes each"
       % (len(plan), min(len(p[2]) for p in plan), max(len(p[2]) for p in plan)))
    res = run_jobs(order)
    ck("every_scheduled_job_returned_a_value", len(res) == len(order),
       "%d/%d" % (len(res), len(order)))

    # ---------------------------------------------------------------- step 5
    head("Step 5: the twenty determinants, recomputed exactly and compared")
    dets = {}
    for key, mi, ps, bsq in plan:
        X = crt_lift([res[("det", key, mi, p)] for p in ps], ps)
        printed_det = int(_printed_det(key, mi))
        ck("determinant_%d_of_cell_%s_recomputes_to_the_printed_integer" % (mi, key),
           X == printed_det,
           "%d bits, %d primes" % (abs(X).bit_length(), len(ps)))
        ck("determinant_%d_of_cell_%s_is_within_its_hadamard_bound" % (mi, key),
           X * X <= bsq)
        dets[(key, mi)] = X

    head("Step 6: p-uniformity -- one Euclid gcd per cell")
    for c in CELLS:
        key = "r=%d,j=%d" % (c["r"], c["j"])
        a, b = dets[(key, 0)], dets[(key, 1)]
        g = gcd(a, b)
        ck("the_two_determinants_of_cell_%s_are_coprime_so_the_cell_is_p_uniform" % key, g == 1,
           "gcd=1, bits %d and %d" % (abs(a).bit_length(), abs(b).bit_length()))
        ck("neither_determinant_of_cell_%s_is_a_unit_so_both_minors_are_needed" % key,
           abs(a) != 1 and abs(b) != 1,
           "smallest prime factors below 10^4: %s and %s" % (_spf(a), _spf(b)))
        ck("rank_of_cell_%s_is_maximal_over_F_p_for_p_in_2_3_5" % key,
           all(res[("rank", key, p)] == c["s"] for p in SMALL),
           "rank=%d=s at p=2,3,5" % c["s"])

    head("Step 7: from the ten cells to every r in [6,21], in every characteristic")
    for a_, b_, q_ in blocks:
        ck("block_%d_to_%d_has_its_injective_endpoint_at_(%d,%d)_and_its_surjective_endpoint_at_(%d,%d)"
           % (a_, b_, b_, q_, a_, q_ + 1),
           any((c["r"], c["j"], c["role"]) == (b_, q_, "inj") for c in CELLS)
           and any((c["r"], c["j"], c["role"]) == (a_, q_ + 1, "surj") for c in CELLS))
    for r in range(RLO, RHI + 1):
        q = qs[r]
        blk = [b for b in blocks if b[0] <= r <= b[1]][0]
        hf = [comb(j + 3, 3) - r * comb(j - D + 3, 3) if D <= j <= q
              else (comb(j + 3, 3) if j < D else 0) for j in range(0, 25)]
        ser = froberg_series(r, 24)
        ck("hilbert_function_of_r_%d_equals_the_truncated_froberg_series_in_degrees_0_to_24" % r,
           hf == ser, "q=%d, block [%d,%d]" % (q, blk[0], blk[1]))

    head("Step 8: controls, both polarities, each named to the object it tests")
    for n in (12, 30):
        JM = [[1 if a != b else 0 for b in range(n)] for a in range(n)]
        ps, bsq = crt_primes_for(JM)
        X = crt_lift([res[("ctrlJ", n, p)] for p in ps], ps)
        want = (-1) ** (n - 1) * (n - 1)
        ck("control_engine_reproduces_det_of_J_minus_I_at_n_%d" % n, X == want,
           "det=%d, closed form (-1)^(n-1)(n-1)=%d" % (X, want))
    N_small = [row[:40] for row in subs[("r=21,j=9", 0)][:40]]
    ps_s, _ = crt_primes_for(N_small)
    X_crt = crt_lift([det_mod(N_small, 40, p) for p in ps_s], ps_s)
    X_bar = det_bareiss(N_small)
    ck("control_two_independent_engines_agree_on_a_40x40_block_of_a_real_minor",
       X_crt == X_bar, "crt=bareiss=%d" % X_bar)
    X_ci = crt_lift([res[("ctrlCI", p)] for p in ps_ci], ps_ci)
    ck("control_forced_positive_monomial_complete_intersection_r_4_j_13_has_det_plus_or_minus_1",
       abs(X_ci) == 1 and len(N_ci) == 4 * comb(6 + 3, 3),
       "size %d, det=%d, so that cell is p-uniform at every prime" % (len(N_ci), X_ci))
    ck("control_proved_deficient_r_6_with_two_repeated_generators_has_rank_336_not_504",
       all(res[("ctrlDEF", p)] == 4 * comb(6 + 3, 3) for p in (2, 1000003))
       and 4 * comb(6 + 3, 3) < 504,
       "rank=%d < s=504 at p=2 and p=1000003, so the decider says NO"
       % res[("ctrlDEF", 2)])
    ck("control_tampering_with_one_mask_bit_changes_the_(21,9)_determinant",
       res[("ctrlTAM", 0)] != res[("ctrlTAM", 1)],
       "one bit of the mask of G21 flipped")
    key0 = "r=21,j=9"
    d0 = dets[(key0, 0)]
    p0 = primes_below_2_62(1)[0]
    ck("control_a_determinant_off_by_one_fails_the_residue_test",
       (d0 + 1) % p0 != res[("det", key0, 0, p0)],
       "the comparison in Step 5 is not vacuous")
    ck("control_a_single_minor_never_certifies_p_uniformity_on_its_own",
       all(gcd(dets[(("r=%d,j=%d" % (c["r"], c["j"])), 0)],
               dets[(("r=%d,j=%d" % (c["r"], c["j"])), 0)]) != 1 for c in CELLS),
       "gcd(det_0, det_0) = |det_0| != 1 in all ten cells")

    # ---------------------------------------------------------------- scope
    print("")
    note("NOT RE-RUN, and the paper says so in its Section 7 and in REVIEW_NOTE.md's Scope: "
         "(1) the identity of the printed family with the ancillary file published with "
         "arXiv:2608.24797 -- the sha256 above is recomputed from the printed masks, but the "
         "external file is not fetched, since this program reads no file and no network; "
         "(2) r <= 5 and r >= 22, and every (n,d) other than (4,7) -- no claim is made and none "
         "is checked; (3) the 168-prime sweep reported in the source note: only p = 2, 3, 5 are "
         "re-run here as an independent confirmation, and the gcd argument, not that sweep, is "
         "what covers every prime; (4) the minor-selection heuristic -- which s x s minors to "
         "take is not re-derived, and it does not need to be, because the printed index sets "
         "are re-extracted and their determinants recomputed from scratch; (5) any claim about "
         "GENERIC forms beyond the specialisation argument of Section 2, which is a proof and "
         "not a computation.")
    print("")
    if _state["fail"]:
        print("VERDICT: %d OF %d CHECKS FAILED" % (_state["fail"], _state["fail"] + _state["pass"]))
        return 1
    print("VERDICT: ALL %d CHECKS PASS" % _state["pass"])
    return 0


def _modprod(ps):
    M = 1
    for p in ps:
        M *= p
    return M


def _printed_det(key, mi):
    for c in CELLS:
        if "r=%d,j=%d" % (c["r"], c["j"]) == key:
            return c["minors"][mi]["det"]
    raise KeyError(key)


def _spf(n):
    n = abs(n)
    for q in range(2, 10000):
        if n % q == 0 and is_prime(q):
            return q
    return "none below 10^4"


if __name__ == "__main__":
    sys.exit(main())
