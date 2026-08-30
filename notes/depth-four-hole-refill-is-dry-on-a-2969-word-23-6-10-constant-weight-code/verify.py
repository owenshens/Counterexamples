#!/usr/bin/env python3
"""Verification program for the paper

    "Depth-four hole refill is dry on a 2969-word (23,6,10) constant-weight code"

Python 3.9+, STANDARD LIBRARY ONLY (hashlib, itertools, math, sys, time).  No third-party
package, no network, no external data file: the code C_0 is the listing printed in Appendix A of
the paper, reproduced verbatim in C0_LISTING below, and every number the paper states is
re-derived from those 2969 words.

Runtime about 1-2 minutes on one core; peak memory about 1.5 GB.

WHAT IT DOES

  A  the object            C_0 parses, is a (23,6,10) constant-weight code, sha256 identity
  B  the blocker census    built TWICE by disjoint methods (2-swap neighbourhoods; 8-subset
                           containment index) and the two are compared word for word
  C  calibration           the eight integers the source paper publishes about C_0, reproduced,
                           including the two DIFFERENT pair-counting conventions it uses
  D  the deletion space    the 5,094,818 sets, with soundness and saturation checks
  E  the sweep             exact maximum clique in Cand(D) for every one of the 5,094,818 sets
  F  attainment            the explicit break-even exchanges at k = 1, 2, 3, 4

WHAT IT DOES NOT DO -- printed again at the end of the run as `NOT RE-RUN:` lines.
"""
import hashlib
import itertools
import sys
import time
from math import comb

# ---------------------------------------------------------------------------
# THE OBJECT.  Appendix A of the paper, verbatim: 2969 six-hex-digit words, 16 to a line, in the
# order of the source listing file, 13 to a line.  Bit i of a word (i = 0 .. 22, least significant first) means
# "coordinate i+1 belongs to the word".
# ---------------------------------------------------------------------------
C0_LISTING = """\
2822BD 283BB0 2843B6 284F89 287AC1 287C98 28863E 288AAB 288BC6 28935C 289CE2 28A14F 28A5D1
28B52A 28C6D2 28CC39 28E4E4 28EA70 2901AF 2906CD 29196A 291C33 29319C 293361 293626 293C4C
29459A 294CE8 295259 295534 29588B 29630D 298E58 299703 29A694 29A82E 29D065 29D2B0 29DB08
29E815 2A0697 2A074E 2A09B3 2A1C56 2A1E49 2A21E9 2A2A47 2A2EE0 2A5163 2A5A94 2A615A 2A6A8A
2A7506 2A8955 2A9E30 2AA534 2AAB48 2AB4C8 2ACE44 2AD590 2AE423 2B0176 2B083D 2B14D1 2B16A8
2B23A4 2B246A 2B28D4 2B4750 2B59C0 2B8E81 2B9149 2B9898 2BC60A 2BE8A0 2C04F3 2C2137 2C2C4B
2C3554 2C5938 2C609B 2C720E 2C7C22 2C859C 2C9287 2C964A 2CA362 2CA929 2CC0AE 2CD0D8 2CEE80
2D03D2 2D0F0A 2D15A2 2D1916 2D1A8C 2D2A32 2D30E8 2D6151 2D7601 2D9870 2DAA44 2DC168 2DC704
2DC892 2E0CA9 2E1F80 2E21C6 2E22A3 2E2E18 2E4515 2E46B0 2E6249 2E8C0E 2EA06C 2EC132 2ED228
2EE181 2F04E4 2F2C05 2F3430 2F4921 2F601C 2F7082 2F8942 2F9124 301AC7 3026D3 302B2B 3030DE
303374 30351B 304D55 304ECC 305467 305AF0 305D92 30628F 308787 309C3C 30AD25 30B546 30C3E2
30D0EC 30F603 30F988 3107B4 311F11 31274A 3132A3 3149C6 31538C 3160BA 316523 316B60 316E82
317630 318B2C 319076 31921B 31A097 31A271 31AE09 31B581 31D6C0 31D905 320277 3209AD 320CF2
32155C 321F0A 322BD0 324B16 3264CA 32A0E6 32A11E 32A469 32C709 32DA0C 330593 33061D 33085E
3311AA 3312D2 331A89 3322C5 333035 3340F4 336213 336849 33890B 3392A4 33B043 33D118 3400EF
340F38 3411BC 34160F 341A6C 3423C9 343A19 3445B1 34581E 346654 3468E8 3489A3 3490D5 34AC1A
34C24D 34E12A 34FA40 35089B 350B54 352C34 3534D0 353528 35384A 35446A 3550C9 3585C4 359192
35A226 35A318 35C03C 35C053 35CA88 3606AA 3642D8 365381 365434 367608 368839 3698E0 36A20B
36A88C 36C164 36D04A 370D22 370E90 371238 371C41 373106 376130 37A640 3805E6 380CD9 382A5C
382B85 3834AC 383C83 3854D4 386C61 387191 38918D 389453 38984E 38B0F0 38C1B8 38D28A 38E648
38F026 39015B 391D88 394037 39489C 395C06 39606C 3962D0 396584 398A07 398CA4 39960C 39A4C2
39B920 39E00B 3A0678 3A0A1B 3A0C27 3A1871 3A2598 3A3909 3A4724 3A5348 3A5612 3A5CA0 3A816A
3A82CC 3A85A1 3A8896 3A9464 3AAE02 3AB214 3AC02D 3ACB80 3AD803 3B0AA2 3B3540 3B4C11 3B7A00
3B81D0 3C08B5 3C09CA 3C2D0C 3C302B 3C3252 3C415C 3C423A 3C448D 3C8549 3CD121 3CE211 3D022D
3D0C52 3D1505 3D4845 3D50A4 3D6302 3D82E0 3E1322 3E1944 3E2A24 3E3890 3E410B 3E45C0 3E8A41
3E8D10 3F0388 3F0868 3F8085 3FD400 40197E 402EEC 4034D7 403B99 404677 4047BC 4049DD 405AE6
4061B7 40655B 406AF1 406B3A 40724F 407768 407C2E 408AD7 408BF8 408CBE 4095EC 4099CB 409D27
409E4E 40AE59 40B70B 40BB2C 40BCA9 40C7CA 40D1D6 40DC1D 40E0EE 40E52D 40E687 40F331 40F958
41053F 410DF4 410E9D 411B2B 411C6D 4127D8 412C7A 412F13 4138E3 4139CC 413E8A 414B6C 414CCE
4152DA 4154B3 415719 4163E2 41670E 41698B 416E29 41703D 417654 4186E3 4189AD 41925D 419372
41AB4A 41B0BA 41B195 41C07B 41C55C 41CA35 41D383 41D6A4 41D8C5 41D90E 41E29C 41E8D2 41F164
41F443 41F588 41FA22 4203F6 420B3D 421CEA 4226B3 4229DA 422D36 4230AF 423D43 4242EB 42456E
424CAD 424F07 4255D8 425696 425A55 426365 4268C7 426DE0 426E4A 427471 4278B2 427F10 428DC6
4291B3 42999C 429A36 429D70 42A2CD 42A467 42AAE2 42AF21 42B354 42B8D1 42C25E 42C395 42C4D3
42C8F4 42D06D 42D3E0 42D623 42D70C 42DD81 42E0B9 42E172 42E634 42EA13 42F107 42F28A 42F4C4
43094F 430AE5 430DA3 430E2E 4310F9 431653 4318D6 432387 432BA8 432CC9 432F44 43321E 4334B4
433A31 4340BE 434333 4346F0 434D52 43528D 435346 435827 4361D4 436259 436417 436538 43706A
4371A1 437AC0 43829B 438C39 439325 439487 439F02 43A07C 43AD90 43B1C2 43B660 43B909 43C50B
43C645 43C9C8 43CA86 43D151 43D238 43E4A2 43E861 43F814 4409E7 440F8E 44131F 441BD4 441CB5
4426CB 442D69 442E27 4431F8 443933 4444FA 444B93 44516B 445B25 445DA8 44661D 446974 44738C
44817D 44863B 44968D 449BA2 449E61 44A2B5 44A39A 44A4DC 44A81F 44B46A 44B586 44B945 44C32E
44C44F 44CCA3 44CF09 44D0F1 44DE90 44E1C9 44E263 44E6A8 44EE44 44F216 4505D3 45065E 450B71
450EB2 4511AE 4516E8 45233C 4524B9 452996 452C55 453266 453349 45340F 4548E9 454E43 455057
45549C 455570 4562C5 456CA4 4572B0 45781A 458B07 459934 459C13 45A12B 45A9E0 45AE81 45B618
45C1B8 45D429 45DB40 45E911 46049F 460AD9 4611CD 461E83 461F48 462257 462678 4632E1 46332A
463A0D 463C64 4646CC 46513C 465272 4662A6 46642B 466781 466D0C 467095 46834B 4684E9 468B64
468D2A 469B11 46A1AC 46A60E 46C037 46DC42 46E2D0 46F920 470BC2 471398 47143A 471515 472925
473950 474B14 474C91 4750E4 477203 4780D5 47811E 478E50 479063 47944C 479CA0 47A232 47B284
47C2A1 47C524 47D092 47E046 47E308 4804EF 4811BD 4816DC 481D59 481E17 48232F 482BC3 4832B6
48335A 48343B 483CF0 48469B 48562D 485C63 485F0A 48647C 486B15 4870E5 4874CA 487E84 488B33
488C75 489369 4896AA 489B85 48A517 48A9D4 48AE92 48B44D 48B88E 48B962 48C1AB 48C947 48CE26
48D274 48D398 48E34C 48EBA0 48ECC1 48F093 4903BA 490679 490B56 491BE0 491CA6 49215D 49249E
4929B1 492E34 493907 493E41 4942A7 49434B 4945E1 495878 495994 4961AC 496D48 49816E 49860F
498ACC 498D8A 49911B 49943C 49A265 49A389 49A4E8 49AC23 49BB10 49C2D1 49CCB0 49D648 49D921
49E036 49E419 4A017B 4A07D1 4A1267 4A189B 4A1CC5 4A239C 4A2B70 4A316C 4A3D88 4A3E22 4A4BC4
4A52B1 4A5C34 4A7056 4A866C 4A89E1 4A8E43 4A92D2 4A9619 4AA48B 4AAB06 4AACA4 4AB035 4ACD18
4AD0A6 4AD14A 4AF241 4AF428 4B06CA 4B0F09 4B1932 4B1E90 4B20E6 4B242D 4B3458 4B4075 4B4199
4B461C 4B50C3 4B6292 4B7224 4B82B4 4B8358 4B8C16 4B91A8 4BA053 4BC262 4BC82C 4BE085 4BEE00
4C07A9 4C0E4D 4C0F62 4C1476 4C1D92 4C1E38 4C23E4 4C254E 4C2CAA 4C2D85 4C3691 4C44D5 4C4AB4
4C4C1E 4C4D31 4C53C1 4C5887 4C5D44 4C6632 4C6853 4C7126 4C7A60 4C7C09 4C81B6 4C86F0 4C895A
4C8C99 4C902F 4C9346 4C94C3 4CA351 4CB1A1 4CB624 4CBA03 4CC98C 4CCAC2 4CCC68 4CD115 4CD832
4CE825 4CEA18 4CED02 4CF450 4D00FC 4D0395 4D0837 4D1165 4D1623 4D1989 4D1C4A 4D2172 4D221B
4D2AD0 4D3382 4D34C4 4D3D20 4D41C6 4D450D 4D4664 4D492A 4D4E88 4D5A11 4D6321 4D6483 4D6A06
4D84A5 4D8D41 4D91D0 4DA642 4DAC0C 4DC312 4DC854 4DF880 4E0635 4E09B8 4E11E2 4E12AC 4E190E
4E203E 4E24D2 4E2519 4E304B 4E41A5 4E420F 4E458A 4E4866 4E4949 4E5413 4E58D0 4E60E8 4E6342
4E80CE 4E9078 4E9494 4EA229 4EB112 4ED089 4EDA04 4F00AB 4F0326 4F0447 4F0568 4F21C1 4F224C
4F26A0 4F405A 4F6830 4F8131 4F920A 4FC4C0 4FC803 5005DE 500B5B 5012ED 50172E 5023D5 502739
502E1E 5035C9 503663 50492F 504DC3 504EAA 5055A5 505A33 505B1C 505C5A 506BC8 506C4D 506F24
5078D4 50856B 508E2D 5093B4 509693 509758 5098F2 509939 50A2AB 50A78C 50A8E5 50ACCA 50B153
50B474 50BA15 50C6D4 50C99A 50D26A 50D345 50E433 50F08D 50F1A2 50FE08 510367 51078B 5109EA
51109F 5115B8 5117C4 511AD1 5124A7 51256C 512AC6 512B0D 513278 513C19 5142B9 514396 514C3C
515607 515949 51605E 516171 5174E0 51823E 518457 5191E1 51950D 519B88 519C2A 51A0D9 51A51A
51A6B0 51B306 51CAE0 51CC89 51DD10 51E807 51E928 5206C7 5207E8 520D65 5212BA 521371 521437
5221E3 52226E 52250F 522F82 52305D 523196 5236D0 524157 524978 524AD2 524D94 524E19 5250CE
525664 5258E1 5263B0 526835 527229 527483 527A06 52853C 528672 528AB1 528D13 529166 52920F
52958A 529AC4 529C49 52A495 52A856 52AB18 52B0E8 52B381 52B823 52BD04 52C2AC 52CC0E 52CF40
52D312 52D4B0 52E04B 5301B5 5302DC 530559 530A17 530C9A 5316A1 531A62 53203B 532352 5328F0
533446 533708 533885 5343C1 534463 534B0A 536126 536684 536D01 537098 538269 5383A2 5384E4
5388C3 53905A 539614 53A08E 53A145 53A603 53AA24 53B130 53C01D 53C832 53D184 53DA01 53E450
540E74 5413CA 5414AB 541C47 5422F2 5428BC 542C93 54312D 543287 54364C 543D0A 543EA0 5441EC
544669 544B46 544E85 5450B6 545199 54648E 5469A1 546A0B 546C62 547251 547438 54818F 548BC1
54941E 54A368 54A926 54B990 54CB30 54CC51 54D4C8 54D903 54DC24 54E11C 5504CD 550BA4 55115C
551B12 55262A 552943 55441B 5549D0 555328 555C82 556E10 557504 5580B3 558478 558A19 558C94
5590C6 55A254 55B221 55BC40 55C065 55C14A 55C60C 55E282 560B23 560CA6 560D89 5610D3 5629C4
562E41 563234 56540D 565522 565A88 566912 567148 568296 5690A5 569628 56A071 56B242 56C0E2
56C129 56D054 56E205 56EC80 570076 570730 57110B 571245 57182C 572291 57241C 5730A2 574087
5744A8 574C44 575031 576260 58037C 580F98 58114F 5819A3 581EC2 582077 58219B 5825B4 5833A8
583B44 583C25 58438D 58464E 584B61 584CE4 585C91 58616A 5864A9 586703 586986 587930 5882E6
5886C9 588C87 589916 58A42E 58AD09 58AE60 58B498 58C0B5 58C159 58C638 58D12C 58D40B 58D9C0
58E252 58EC14 5906AC 5909C5 591235 59161A 591C54 59244B 592B22 5930D2 5944D8 594913 595066
59518A 595E20 596318 597281 59780C 59819C 598711 598970 599243 59B0A4 59B148 59CB04 59CC42
5A01AE 5A062B 5A0D4A 5A10F4 5A22A5 5A2349 5A24CC 5A3213 5A4236 5A5780 5A582A 5A5905 5A60D1
5A641A 5A7C40 5A8255 5A8CD0 5A9B20 5AA0B2 5AC283 5ACC21 5AE064 5AE188 5B0A38 5B1286 5B1524
5B18C8 5B1C03 5B3061 5B414C 5B49A0 5B5250 5B8027 5B9091 5BA2C0 5BF002 5C01F1 5C053A 5C0653
5C08D6 5C091D 5C1968 5C1A26 5C1C8C 5C2316 5C2A31 5C2C58 5C38C1 5C40CB 5C43A2 5C4427 5C5152
5C5614 5C6145 5C622C 5C66C0 5C8325 5C8554 5C8863 5C8DA0 5C8E0A 5C9431 5C9A50 5CA1C2 5CAA84
5CC492 5CF300 5D04E2 5D0A83 5D0C29 5D0D06 5D18B0 5D2588 5D2605 5D2864 5D3111 5D4134 5D5441
5D9122 5D9680 5D9805 5DA812 5DC181 5DD018 5DE420 5E006D 5E070C 5E0AE0 5E228A 5E2807 5E409C
5E4311 5E4470 5E801B 5E8834 5E9141 5E9882 5EA610 5EC248 5F0192 5F0851 5F4602 5F6009 5F8908
60077A 6012F3 601575 601B4D 6028BB 603AA5 6047A3 60523E 6056C5 6059B1 60694E 6070D9 60859B
60938E 60B51C 60CCD8 60D472 60D8AA 60E61A 60E83C 61154E 61199A 611D85 6120F5 6122AE 614355
614897 614FC0 6151E8 615C51 616466 616591 617186 617B01 6183E4 618729 618867 618959 6194D4
61A64C 61B431 61ED04 62231B 6225C5 622E94 6233C8 623607 62485B 624E32 62518B 626929 6282A7
628716 629057 62942E 629E88 62A54A 62B2B0 62C531 62C926 62E096 630B91 631274 631961 631C1C
632263 633093 633522 6352A2 637045 6380F2 6392C1 63A215 6406E6 641659 64219D 64226D 643742
64385C 644E2C 6450AD 645516 646307 646358 647521 647A82 648C56 649338 64AA2A 64AB14 64AD88
64B0E4 64C2B2 64C961 64D14C 64D20B 64E055 651724 6521CA 654389 654D18 655261 65680D 657448
659099 65B052 65C48A 65C611 24D266 6605B2 661127 66124E 662983 66348A 664BA0 665868 666898
66A849 66AC30 66B806 66CC05 66D418 6704D8 6708B4 670A0B 673029 674162 676C02 67822C 67A501
6801D7 68071D 680BAC 680D2B 680EB1 681E64 682965 68328B 6848F2 685B50 6862C6 68640F 686651
686A23 6898B4 68A0E3 68AE05 68C11E 68C585 68C760 68F109 68F844 690786 690CC3 694A1A 694C25
697412 698552 699A82 69D08C 6A04BC 6A0A6A 6A1B03 6A288D 6A2C13 6A3151 6A34A1 6A3838 6A38C2
6A4338 6A4469 6A4C86 6A501D 6A61A2 6A8A1C 6A91C4 6AA780 6AC891 6B130C 6B290A 6B4A41 6B5501
6B6488 6B8183 6B8C60 6BB404 6BE140 6C0879 6C10BA 6C268C 6C3B08 6C4543 6C54E0 6C69C0 6C83C8
6C88C5 6CC434 6CD182 6D042E 6D094C 6D2461 6D2710 6D40B1 6D804B 6D9214 6D9508 6DA106 6DCA20
6E0354 6E0A92 6E0D24 6E6444 6E7210 6E810D 6E8451 6E8622 6E9821 6EE00A 6F1640 6F4284 6FA090
700BB2 701656 701DE0 7025AA 702917 703869 70407D 706CB0 707415 7082DA 7084F1 708D4C 709B42
70B720 70C217 70D2A1 70D486 70E1D0 70F812 710C0F 710E61 711133 712D50 713390 713826 714272
715A44 7160C3 716225 71720A 71828D 718526 718E12 71A06A 71A982 71C541 720A8E 721295 721529
721934 721A58 722689 723072 7240B3 725243 72630C 7270A4 728B05 72A960 72C0C5 73058C 732C28
734648 735016 738246 73A0A1 73C290 73D060 7407D0 740AA9 741951 741986 741C32 74205B 7430B1
747046 748172 7484AC 749C81 74C384 74C642 74E409 750169 75030E 751E08 754823 756094 758115
75B00C 760546 7624E0 763510 76402E 767801 768198 770425 7711C0 772B00 779810 78083E 7811D8
78146A 781741 782C46 783E10 7842E8 785306 786819 78830B 789232 78A131 7900CE 790495 7910A9
7921E0 792A88 794528 79A841 79C0A2 7A03C2 7A300E 7A4854 7A88A8 7AA058 7AC502 7B0321 7B0432
7B0E04 7B2114 7B8409 7C1017 7C12C4 7C20A6 7C7088 7CB402 7CC806 7D0258 7D5900 7E0483 7E1209
7E2128 7E4C08 7F2042 705698 1A31CA 141E96 113D45 11593A 111769 1C4CB8 6641D1 09D952 1E4645
60EA89 1E492C 1A6F08 10AD96 123F60 241D2D 02DE58 0BC478 110E3B 1A4E62 180F0F 190EF0 1C2E89
1D0E1C 1D4B09 1C5E48 1D4D60 1CCD05 1DC831 11097D 1C7A05 280D7C 124FA1 1E42A9 1C0C6E 1CCA64
134A71 1146E5 1D852C 1C9D42 2812EE 1CC788 1D5390 1C4AD1 1E5A30 1A5998 192B51 1D6AA0 1D3324
1936C8 182DE8 189965 0C51F4 1C9B0C 1DAF00 1CD894 1F5881 18DB11 1DDA02 1A1BC1 189E23 1D6C0A
194CA3 198BA1 1889D3 308B99 18A3D8 0B0DD8 1803EB 18AB34 1884FC 1961C9 1C6836 1C6370 10F2B8
1D21D4 1D2938 1A22D6 1B32B0 199338 19A1AA 19330B 19E990 1BC941 1AB950 10EB0E 1B04E9 1B90E2
19D4A8 1352E8 1C58E2 1D7850 1167A8 18762A 18669C 0CC51B 18E329 1AA263 19F260 1E3268 1D12AA
19866A 17066C 17AC11 1CBC28 193C62 12ED30 1B6072 1B7128 1F1170 1E15A8 1F9029 1D1A61 0BD40D
19142F 19B033 0C1B4B 1950F1 126D46 18447B 1C7464 0D55C8 1C550E 1E64A2 1E1463 1AF0A1 1B70C4
1670F0 1662C3 1EE038 1073E1 1A40E7 1A6C85 19FC01 1D34A1 191E85 1C5683 1C4197 18E2C5 1CE883
1CA855 1CB107 11C1D5 1AE017 1BA119 1E221D 1A01DD 19A48D 1E108F 1D11C3 1E7114 1FE804 1AC31C
1BD034 18139E 1A943A 1A233A 1B149C 1F301A 1CA09E 1DB414 1EA344 1FC320 1B8B12 1B910E 1D8386
1A844F 1F805C 1E822E 1AA6A8 1E03B4 1FC08A 1AC499 00B0FD 1E00FA 1EC852 0F50B8 0DF130 1F20AC
1D3886 1253D4 13A394 17F210 19A946 1B2983 112B9A 1E2153 1EA50A 1331D1 10796C 1423AE 142B65
14EE21 1CC4E1 1C24C7 1EA425 17B064 13E1E0 0B62E1 176055 0304F7 1364B1 1250BD 1340DB 1F024B
13623C 13E289 1742B2 1EA0C9 0BE254 09E62C 13E50C 0926AB 1541AB 1D2692 1C31B2 0A30F3 1222F9
0BC0B3 0B41EA 1C8695 0AD0D5 1390CD 103CBA 1029CF 0730C7 15211F 09613B 0821FE 1461DA 12E0DC
09A31E 11626B 13B22A 0AF232 12F462 126627 11D863 13C06E 1302AF 1042FE 16342E 1A132D 08733C
15634C 172329 11A2EC 152659 16217C 0C29D9 0D2E68 106E78 0B6C43 12F324 123B8C 183179 1620B7
08F1E8 097469 15642D 1459C5 09B12D 0C2C3D 0F3123 0953A9 14A63C 1456AC 178235 14932B 1586A9
1709E1 1318B3 179283 157093 19D087 11D22D 14433D 0A4CF1 140EE3 14C2A7 0C66A5 0A6639 0E8738
0DA631 1CA472 16F086 17C2C4 0BB096 0B2237 132867 10D05F 00E2DB 1502D7 12C23B 08D2E3 1CB641
0B832B 0BA649 0B9271 16D261 12C475 14946D 0B82C7 0BA1CC 0F2079 0F3289 096E91 136896 0944BD
0993CA 1503F8 12927C 101E5D 1E9298 0E80BD 18A61B 0EE28C 0FC219 0E924D 0981F9 09568E 0DCA85
0FA00F 0F106E 0A147D 0ED01E 09361D 08C65D 0C176C 11E615 14D619 15F109 14D178 179348 0D9E09
0F1B28 0A13F8 112AB5 09A2F2 0DE0A9 06E0E5 0D64F0 10F4D1 256239 0D3453 0DA193 09665A 05D233
09714E 0C707A 0DB90A 0979A2 11F0CA 1560E6 0E3A51 14F035 03E135 0DE445 08D46E 1171B4 155274
077036 155C15 080BF5 1011F7 0150EF 0851DB 091B93 1047D9 0CC3B1 1283E5 0833CD 0D2347 131336
151927 111B4E 1D1646 14744B 14B0E3 0EE660 111CCB 09C4CB 0914FA 1504BE 15107B 0D16B4 1354A6
141BB1 0406FD 1616C9 0AC3C9 12D685 1194B5 0A7287 0E73A0 0A69B4 0CAAE1 0AB2E4 1234E5 03326D
0A764C 12725A 123347 0EB446 1612E6 08A6CE 113E2C 08B678 05347C 15329C 0E30DC 0722F4 11B91C
188E56 127C1C 131C78 0D4C59 0752D1 0D18D5 0F5454 07590D 0F8865 0B19A5 12F845 0ABC61 125738
125876 11E874 093974 1294D6 1E86C2 0E27C8 0212DF 024AB7 0CD6C4 151CE4 0F13C4 0592E5 095B45
0D98A3 0DA8B4 0391F4 03F0AC 11DC4C 0FF048 15A44E 17A0D2 1D8CC8 1484DB 03A0EB 0C85EA 14B68A
15E498 10BC0F 01B28F 11FA84 0EB885 0CB394 05F0D4 0D84D6 158E45 099C47 05A077 0D285E 04666E
08949F 104C9F 084EC7 0D0C8F 0DC24E 0D0AE6 05C4EC 05417E 0C36E2 14D2D2 0D7185 0E427C 10BA66
12AA87 0AAB91 15A985 09B9C1 087C55 0A1F1C 0DE862 0E52CA 05564D 17C449 0F40CD 16618D 124B4D
1701CE 153BC0 0664BC 0266D5 23E0D1 2262F2 156CC1 0C6CCC 0E74C1 1076C6 10770D 08D30F 03903F
0A311F 0C690F 02F21D 037370 0B2595 0115DD 089ED1 121D87 10D59C 02D0FA 06D2B4 1055EA 08E4BA
10E58B 11A3C3 02A1D7 12E351 02A2BE 07289D 120EBC 13BC88 171994 0F350C 1A8D8C 0AE859 03E91A
0C8EAC 0EA25A 03B2D8 03B455 0682F3 118E8E 154D8C 12B631 06B227 07B0B1 02B6C3 0E96A1 04A4AF
06C62D 0E052F 0714AD 0574AA 0C358B 0A9D0B 11358E 081DCE 092F8C 0CAD64 0A2D4D 0DE3C0 01EEC8
01E94D 057C46 168947 0EA833 169932 00F665 0C9C5C 0E99C8 0D1F50 143C71 0136F1 143595 04E479
0DA558 16B458 04F34A 0C9173 083F29 0C1DE1 141DD8 26344D 058CF1 16C991 16039B 0E8317 0F4383
122F15 148C37 12C1B6 0271E6 0743AC 0F8589 02B58D 0E0B8D 03A6A5 01C78D 170785 07035D 0107EE
050FC9 2702E9 0EC54C 27A154 05D31C 06A699 12E692 17861A 12C55A 16549A 06F198 0F6988 0B4D0E
134656 01D999 0CDC8A 14CCC6 00DBCC 0E8AD4 105B8B 006BAD 2483AD 1078A7 03F883 017C8D 0D9D84
137382 187943 034E8B 077419 085937 04543F 03F606 0A369A 01B4E6 05BC25 08A9A7 099B26 017327
02E3A3 03A366 0AE5C2 05EB24 132DA4 05B3A8 117558 035197 03170F 0313E3 12196B 0298E7 0233B5
2817A5 0195AB 11D70A 130F43 0B3B42 13EA42 09EB03 0E9A62 02AA75 09CA69 158B62 159E30 05962E
16436A 079256 07D688 03D4E1 0396B2 129F90 0C8F83 078E23 033F81 00279F 0075B9 0166B6 035635
177C20 04FB81 06DC31 0D1539 061975 161E25 179C06 153E03 057961 0D5D03 21595C 22F074 005CFC
027EA8 421FA4 011BBC 0307B9 1087BA 089DB8 023759 0335E8 03A78A 0721BA 03255E 07364A 07614B
162E32 023676 06EC26 16CE14 055EA1 06CEC1 060DEC 07AB41 07ACC4 15ACA2 03AD07 1489F4 0287DC
028DB5 00BE9C 01875B 03C857 012E4F 014B1F 135F04 027D25 02A53B 03856D 27053C 03D24B 027A63
04237B 062F0B 082E73 0A4B53 2360A7 238A4D 1816B9 001EAF 0EB0AA 26BE01 0CBD11 06F053 129A53
00BAB3 20B26B 009C7B 011A77 01FA51 11B652 212AD9 01395B 25E858 24F42C 04FCE0 163CC2 0374D2
023DD4 030AFA 22743A 0A5D68 02976A 268E68 27A4A8 071D8A 133D12 057F08 111FA2 01373A 0133D6
003BEA 061EF0 214AE3 217073 21B178 04BB70 228379 14AB13 22CAB8 03C3D2 236398 11C4F2 0C0D57
030F96 160D1E 2783B0 056B52 0C0EDA 05D45A 01FC38 06465B 05C761 059D68 108F71 147B22 0BDC22
031D66 0BAD28 035DB0 03CF11 052BA3 2527E0 032D71 0B84AE 01CBAA 01CC2F 11E332 23C236 01D8B6
09DE14 10D636 005F95 0CDF20 154E26 052D2E 0D6D14 15C916 1547C2 1C4F12 056C33 0D8C3A 0D2DC2
05C9C3 0CCDD0 0519F2 13CD82 074CE2 2329E2 21F01E 233714 239B50 053B15 00575E 04CE72 01E56A
0C4E2B 10C76C 024F74 04CB55 064F98 04CD3C 040DBB 0169F8 054ED4 146B94 047E34 00EC5E 02DD16
00799E 04DE07 08BC36 0C3F06 0C2A97 049997 0C7496 039A95 09BEA0 21BA34 02BDA2 213C96 2510B7
20B1B6 20BB07 042CF6 045CD3 24B632 0F0672 02AF52 20EF28 0C6F41 086B66 0C5A56 216936 066B31
272853 0CB8D2 007D72 273A60 05AE16 05914F 26169C 252687 063B92 034DC5 222D8E 0C4DA6 025FC2
055B86 175942 0649D6 2618DA 257990 03FD40 0B7911 026D93 226743 25F502 06631E 02C6E6 075E12
26C4D4 274B48 07D12A 01A5BC 01AB39 03B313 042F5C 0443CF 0549B5 0789A6 08FE42 0A0FAA 0A44DE
0A54AB 0B531A 0C4BE8 0D5C2C 0DC6A2 0E5626 0E5709 0ED350 0FC194 0FE412 10323F 1099AE 10A55D
11C539 11CA93 12478E 127A91 12ADC1 12D1C3 12DAA2 135165 135589 13700F 137641 138173 1386D1
138D54 13D8D0 143956 14835E 14A1B9 14AED0 14B1CC 146D19 160CD5 165217 16911D 17E023 182EA6
185BA4 1862B3 18BCC4 18DC70 1909B6 1911EC 192C17 19432E 19511D 196925 197216 1992D4 199C92
19C41E 1A0D39 1A9C15 1AC2F0 1B1057 1B2E21 1B4295 1B5223 1B544A 1B9E40 1C0FC4 1C1C1B 1C1D34
1C8279 1C8AB2 1C94A6 1CC343 1CE606 1D5432 1DC650 1E0B58 1E3684 1E4A86 1E604E 1F4518 1F520C
1FB180 208EF4 20EAA6 276624 460C73 4DD406 00B5DA 04715D 05E60B 078F0C 07C970 08C37A 08F689
0943DC 094F38 09C9E4 09EC86 0A6E16 0B2CB2 0B4EA4 0DC127 0DD491 10CE4B 112CDC 11CB58 126AE4
127133 12AE4C 133A54 139D21 13A436 13C307 13C8A5 18266D 184E35 18BB82 18F11A 194ACA 196E44
19834D 1A2C74 1AD246 1B0B64 1B0CC6 1C1355 1C3718 1D304D 1DD144 1F84B0 218DB2 244CE5 251CB8
2693C2 08B257 094C76 0B0E55 0D4617 0EA1F0 0F091B 0F1A07 128F26 1326E2 13D413 144973 18ACB1
18C1CE 18CD2A 1A3926 1D0475 1F1611 24D589 267C50 26E710 276A81 01AC9B 01CF46 02ACF8 115353
1269AA 128CAB 16DD08 17248B 172D48 18389D 1CE948 219798 21CB94 2454CE 250567 25511B 270E46
274471 0960D7 122C5B 186CD2 254732 260B72 0CE156 106A57 10E167 2023E7 2065D6 20E953 21AE62
22B931 00E3F4 26790A 0E296A 0183B7 602A76 0189DE 23819D 031ECC 2178E4 217295 0215BE 0E18B6
2248EE 2C1867 083567 24316E 074267 028E1F 04EC95 0C78B1 033AA6 03B872 02B94E 08FD0C 033C2B
03DA64 22D293 2198D3 059EC2 049F1A 26CB03 346C07 1C7D80 03CC9C 275C84 2635A4 06B169 063D38
0EA91C 0A887E 068576 05169B 242F91 1D0D91 05B4C9 0521ED 092CE5 0659A3 2E48C3 0EC06B 06E8CA
0362CE 21E247 20C0F7 21C517 27440F 2F9013 27D205 23D720 13CE28 04E9B2 26D4A2 068D59 199859
085EB2 0453BA 00F92B 2436A9 01EDA1 0A9786 069A8E 23D0C6 0AD139 029F45 28554D 0279C9 23D829
0D0F25 159551 029BA9 128BCA 141772 181A7A 0C0B3E 2294B9 12B09B 0D8B98 0F0AB1 16ABA0 72E222
23E268 18813F 008B6F 0C8667 043EC5 2672C4 266945 20A9EC 19A561 15A869 10C9E9 15D9A0 025B2E
230B27 00A7E9 20E6B1 22272D 26282F 242ACE 0E702D 2A6C2C 23546C 216E1C 10653E 22413F 20247F
1638A9 2139A9 23B984 259B21 2431D3 08CA9E 04C796 065547 0A584F 0EC487 076586 25E18C 227817
05C09F 0F029E 3E2416 268C93 2605CB 2155C3 11454F 394643 3D8423 20DB32 02C98F 24E4C3 006CEB
0C61E3 166561 22C5E8 02CD63 201F63 648745 0CA70D 02EF84 20778A 2C65A8 0493D9 2E1199 2A51AC
20F1C5 1D0333 28133B 061A3B 165859 1A5551 20D751 005B79 225D19 0C5671 0C529D 264A1D 0ACB25
2C4365 1389B8 30C87A 382972 1C2D23 3E2701 26B30C 349714 00973D 0C9A35 261731 271552 6A151A
42BC1A 23AA92 2DA28A 0D8A53 2C025F 2A80DB 28A99A 24B88B 038FE0 1897E0 2CB01D 2CDC41 207F44
40BFC0 047AD8 0877D0 203E5A 047713 087A1B 20DDA4 08D8AD 10ECAC 1CE1A4 24CD4A 239C4A 208CCF
20BD49 10BAC9 05B744 098774 058A7C 24A2F8 2588EA 619A68 69A238 392439 098D1D 25845D 27483A
036F22 0B6564 07686C 172A0E 1B2706 385A29 21562B 055A6A 095762 08E19D 0AE12E 0BBA0C 0C39AC
0F622A 12B8B4 204BDA 2495F0 369503
"""

PUBLISHED_SHA256 = '2cfc5d1bd8e2e61dc73a2ab31e3ec985d268d67cc635468912828c8a21b5cbfd'

# The single deletion set exhibited in the paper (0-based indices into the listing).
EXHIBIT_D = (27, 31, 585, 1333)
EXHIBIT_D_WORDS = ('294CE8', '29630D', '463A0D', '623607')
EXHIBIT_CAND = ('0B630D', '214D69', '423E0D', '4A3705', '4E3305')
EXHIBIT_BLOCKERS = {'0B630D': (31,), '214D69': (27,), '423E0D': (585, 1333),
                    '4A3705': (1333,), '4E3305': (585,)}
EXHIBIT_CONFLICTS = ((('423E0D', '4A3705'), 8), (('4A3705', '4E3305'), 9))
EXHIBIT_CLIQUE = ('0B630D', '214D69', '423E0D', '4E3305')
EXHIBIT_CLIQUE_INTS = (5, 6, 7, 4, 2, 7)

T0 = time.time()
_np = [0]
_nf = [0]


def ck(name, ok, detail=''):
    """One check.  `PASS <name>  <detail>` or `FAIL <name>  <detail>`."""
    tag = 'PASS' if ok else 'FAIL'
    (_np if ok else _nf)[0] += 1
    sys.stdout.write('%s %s%s\n' % (tag, name, ('  ' + detail) if detail else ''))
    sys.stdout.flush()
    return ok


def say(*a):
    sys.stdout.write('    ' + ' '.join(str(x) for x in a) + '\n')
    sys.stdout.flush()


def stamp(msg):
    sys.stdout.write('--- [%6.1fs] %s\n' % (time.time() - T0, msg))
    sys.stdout.flush()


POP = bytes(bin(i).count('1') for i in range(1 << 16))


def popc(x):
    return POP[x & 0xFFFF] + POP[x >> 16]


def hx(v):
    return '%06X' % v


# ===========================================================================
# A.  THE OBJECT
# ===========================================================================
stamp('A. the object')
TOK = C0_LISTING.split()
WORDS_HEX = [t.upper() for t in TOK]
ck('A1-listing-parses',
   len(WORDS_HEX) == 2969
   and all(len(t) == 6 for t in WORDS_HEX)
   and all(c in '0123456789ABCDEF' for t in WORDS_HEX for c in t)
   and all(int(t, 16) < (1 << 23) for t in WORDS_HEX),
   '%d six-hex-digit tokens, every value < 2^23' % len(WORDS_HEX))
W = [int(t, 16) for t in WORDS_HEX]
N = len(W)
WSET = set(W)
ck('A2-words-distinct', len(WSET) == N, '%d distinct words' % len(WSET))
ck('A3-constant-weight-10', set(popc(v) for v in W) == {10},
   'popcount set %s' % sorted(set(popc(v) for v in W)))
mx = 0
for i in range(N):
    wi = W[i]
    for j in range(i + 1, N):
        p = POP[(wi & W[j]) & 0xFFFF] + POP[(wi & W[j]) >> 16]
        if p > mx:
            mx = p
ck('A4-min-distance-6-and-tight', mx == 7,
   'max pairwise intersection over all %d pairs = %d, so min distance = 2*(10-%d) = %d'
   % (N * (N - 1) // 2, mx, mx, 2 * (10 - mx)))
sha = hashlib.sha256(''.join(WORDS_HEX).encode()).hexdigest()
ck('A5-identity-sha256', sha == PUBLISHED_SHA256, sha)

# ===========================================================================
# B.  THE BLOCKER CENSUS, TWICE
# ===========================================================================
# A CANDIDATE is a weight-10 word on 23 points that is not in C_0.  Its BLOCKERS are the words of
# C_0 it meets in >= 8 points, i.e. the words it is NOT compatible with.
stamp('B. blocker census, method 1: 2-swap neighbourhoods')
cell = {}
put = cell.setdefault
for i in range(N):
    d = W[i]
    ib = [b for b in range(23) if (d >> b) & 1]
    ob = [b for b in range(23) if not (d >> b) & 1]
    a1 = [1 << b for b in ob]
    a2 = [(1 << x) | (1 << y) for x, y in itertools.combinations(ob, 2)]
    for x in ib:                                  # |s cap d| = 9
        base = d ^ (1 << x)
        for a in a1:
            put(base | a, []).append(i)
    for x, y in itertools.combinations(ib, 2):    # |s cap d| = 8
        base = d ^ ((1 << x) | (1 << y))
        for a in a2:
            put(base | a, []).append(i)
NCAND = comb(23, 10) - N
hist1 = {}
for v in cell.values():
    k = len(v)
    if k <= 4:
        hist1[k] = hist1.get(k, 0) + 1
ck('B1-candidate-count', NCAND == 1141097,
   'C(23,10) - 2969 = %d - %d = %d' % (comb(23, 10), N, NCAND))
ck('B2-no-zero-blocker-candidate', len(cell) == NCAND,
   'every one of the %d candidates has at least one blocker, so C_0 is a MAXIMAL (23,6,10) code; '
   'zero-blocker candidates = %d' % (NCAND, NCAND - len(cell)))
ck('B3-exactly-k-histogram-method-1',
   [hist1.get(k, 0) for k in (1, 2, 3, 4)] == [70, 1178, 6503, 22496],
   'exactly 1/2/3/4 blockers: %s' % [hist1.get(k, 0) for k in (1, 2, 3, 4)])
cum = []
t = 0
for k in (1, 2, 3, 4):
    t += hist1.get(k, 0)
    cum.append(t)
ck('B4-cumulative-pool-sizes', cum == [70, 1248, 7751, 30247],
   '#{|B| <= 1,2,3,4} = %s  (the source publishes 1248, 7751, 30247)' % cum)
# blocker set (as a sorted tuple) -> the candidates it belongs to, for |B| <= 4
bs = {}
bmap = {}
for m, v in cell.items():
    if len(v) <= 4:
        tup = tuple(sorted(v))
        bs.setdefault(tup, []).append(m)
        bmap[m] = tup
del cell
ck('B5-pool-partitions-into-blocker-sets',
   sum(len(v) for v in bs.values()) == 30247 and len(bmap) == 30247,
   '%d candidates with at most 4 blockers, carried by %d distinct blocker sets'
   % (len(bmap), len(bs)))

stamp('B. blocker census, method 2: 8-subset containment index (independent)')
# |s cap c| >= 8 with both of weight 10 holds iff c contains one of the C(10,8) = 45 8-subsets of
# s.  So index every 8-subset of every codeword, then union.  This method shares no line of
# reasoning with the 2-swap enumeration above.
idx8 = {}
for i, c in enumerate(W):
    bits = [b for b in range(23) if (c >> b) & 1]
    bi = 1 << i
    for T in itertools.combinations(bits, 8):
        m = 0
        for b in T:
            m |= 1 << b
        idx8[m] = idx8.get(m, 0) | bi
get8 = idx8.get
hist2 = {}
bmap2 = {}
for cbits in itertools.combinations(range(23), 10):
    m = 0
    for b in cbits:
        m |= 1 << b
    if m in WSET:
        continue
    u = 0
    for T in itertools.combinations(cbits, 8):
        t8 = 0
        for b in T:
            t8 |= 1 << b
        u |= get8(t8, 0)
    k = bin(u).count('1')
    if k <= 4:
        hist2[k] = hist2.get(k, 0) + 1
        blk = []
        uu = u
        while uu:
            b = uu & -uu
            blk.append(b.bit_length() - 1)
            uu ^= b
        bmap2[m] = tuple(blk)
del idx8
ck('B6-exactly-k-histogram-method-2',
   [hist2.get(k, 0) for k in (0, 1, 2, 3, 4)] == [0, 70, 1178, 6503, 22496],
   'exactly 0/1/2/3/4 blockers: %s  (the 0 is an exhaustive sweep of all C(23,10) masks and is '
   'the independent proof that C_0 is maximal)' % [hist2.get(k, 0) for k in (0, 1, 2, 3, 4)])
ck('B7-the-two-censuses-agree', bmap == bmap2,
   'the %d (candidate, blocker set) pairs are identical under both methods' % len(bmap))
del bmap2

# ===========================================================================
# C.  CALIBRATION AGAINST THE SOURCE'S PUBLISHED INTEGERS
# ===========================================================================
stamp('C. calibration: the source paper\'s own pair counts')
# Two words CONFLICT iff they meet in >= 8 points.  Counting all conflicting pairs in a pool of
# weight-10 words exactly, without touching all C(pool,2) pairs: with A_j = sum over j-subsets T
# of C(n_T, 2), where n_T counts pool words containing T, one has A_9 = N_9 and
# A_8 = N_8 + 9 N_9, where N_k is the number of pairs meeting in exactly k points (N_10 = 0
# because the pool words are distinct).  Hence #conflicting = N_8 + N_9 = A_8 - 8 A_9.
POOLS = {}
for k in (2, 3, 4):
    POOLS[k] = sorted(m for m in bmap if len(bmap[m]) <= k)


def all_conflicting(pool):
    c8 = {}
    c9 = {}
    for m in pool:
        bits = [b for b in range(23) if (m >> b) & 1]
        for T in itertools.combinations(bits, 8):
            c8[T] = c8.get(T, 0) + 1
        for T in itertools.combinations(bits, 9):
            c9[T] = c9.get(T, 0) + 1
    a8 = sum(v * (v - 1) // 2 for v in c8.values())
    a9 = sum(v * (v - 1) // 2 for v in c9.values())
    return a8 - 8 * a9


def share_and_conflict(pool):
    """(#pairs sharing at least one blocker, #pairs that BOTH conflict AND share a blocker)."""
    n = len(pool)
    byb = {}
    for j, m in enumerate(pool):
        for e in bmap[m]:
            byb.setdefault(e, []).append(j)
    share = set()
    for L in byb.values():
        for a in range(len(L)):
            ia = L[a]
            for b in range(a + 1, len(L)):
                jb = L[b]
                share.add(ia * n + jb if ia < jb else jb * n + ia)
    conf = 0
    for key in share:
        i, j = divmod(key, n)
        if popc(pool[i] & pool[j]) >= 8:
            conf += 1
    return len(share), conf


ac = {k: all_conflicting(POOLS[k]) for k in (2, 3, 4)}
sc = {k: share_and_conflict(POOLS[k]) for k in (2, 3, 4)}
ck('C1-all-conflicting-pairs-in-the-1248-pool', ac[2] == 5261,
   'the source publishes 5261 conflicting pairs among the 1248 candidates with at most two '
   'blockers; recomputed %d' % ac[2])
ck('C2-conflict-and-share-pairs-in-the-7751-pool', sc[3][1] == 43105,
   'the source publishes 43105; recomputed %d' % sc[3][1])
ck('C3-conflict-and-share-pairs-in-the-30247-pool', sc[4][1] == 887964,
   'the source publishes "about 888k"; recomputed %d' % sc[4][1])
ck('C4-the-two-conventions-are-different',
   (ac[3], ac[4]) == (200311, 2870892) and (sc[3][0], sc[4][0]) == (135753, 3475687),
   'in the 7751 and 30247 pools: ALL conflicting pairs = %d and %d; pairs sharing a blocker = %d '
   'and %d; the conjunction = %d and %d. Reading 887964 as "all conflicting pairs among the '
   '30247" is wrong by a factor 3.2 -- the source is right, the trap is in the reading.'
   % (ac[3], ac[4], sc[3][0], sc[4][0], sc[3][1], sc[4][1]))
ck('C5-conflict-and-share-in-the-1248-pool-is-not-5261', sc[2][1] == 803,
   'so the 5261 of C1 is unambiguously the ALL-conflicting convention, and 803 is the '
   'conjunction; the two conventions cannot be confused on this pool')

# ===========================================================================
# D.  THE DELETION SPACE
# ===========================================================================
stamp('D. the deletion space')
SEED = sorted(bs.keys())
BY_SIZE = {k: [t for t in SEED if len(t) == k] for k in (1, 2, 3, 4)}
say('distinct blocker sets by size: %s' % {k: len(v) for k, v in BY_SIZE.items()})
idx1 = {}
idx2 = {}
idx4 = {1: {}, 2: {}, 3: {}}
SMALL = BY_SIZE[1] + BY_SIZE[2] + BY_SIZE[3]
for tup in SMALL:
    for e in tup:
        idx1.setdefault(e, []).append(tup)
    for pr in itertools.combinations(tup, 2):
        idx2.setdefault(pr, []).append(tup)
for tup in BY_SIZE[4]:
    for r in (1, 2, 3):
        for sub in itertools.combinations(tup, r):
            idx4[r].setdefault(sub, []).append(tup)


def partners(A):
    """Every seed blocker set B with |A u B| <= 4, for a sorted tuple A of size 1..3.

    |A u B| <= 4 iff |B \\ A| <= 4 - |A|.  So: any B of size <= 3 qualifies when |A| = 1; when
    |A| = 2 a size-3 B needs one element inside A; when |A| = 3 a size-2 B needs one and a size-3
    B needs two; and a size-4 B qualifies exactly when A is a subset of it.  The index makes each
    of those a lookup instead of a scan over all 30117 seeds -- check D1 verifies it is complete."""
    a = len(A)
    if a == 1:
        return SMALL + idx4[1].get(A, [])
    if a == 2:
        out = list(BY_SIZE[1]) + list(BY_SIZE[2])     # 2+1 and 2+2 are unconstrained
        seen = set()
        for e in A:
            for tup in idx1.get(e, ()):
                if len(tup) == 3 and tup not in seen:
                    seen.add(tup)
                    out.append(tup)
        return out + idx4[2].get(A, [])
    out = list(BY_SIZE[1])
    seen = set()
    for e in A:
        for tup in idx1.get(e, ()):
            if len(tup) == 2 and tup not in seen:
                seen.add(tup)
                out.append(tup)
    seen3 = set()
    for pr in itertools.combinations(A, 2):
        for tup in idx2.get(pr, ()):
            if len(tup) == 3 and tup not in seen3:
                seen3.add(tup)
                out.append(tup)
    return out + idx4[3].get(A, [])


# D1: the indexed partner enumeration is complete -- checked against a full scan of all 30117
# seeds on 100 members of each size.
bad = 0
for k in (1, 2, 3):
    for A in BY_SIZE[k][:100]:
        want = set(t for t in SEED if len(set(A) | set(t)) <= 4)
        got = set(partners(A))
        if want != got:
            bad += 1
ck('D1-partner-index-is-complete', bad == 0,
   'on 300 sampled blocker sets the indexed enumeration equals a brute-force scan over all %d '
   'seeds' % len(SEED))

DSPACE = set(SEED)
frontier = list(SMALL)
rounds = 0
while frontier:
    rounds += 1
    nxt = []
    for A in frontier:
        sA = set(A)
        for B in partners(A):
            U = sA | set(B)
            if len(U) <= 4:
                tu = tuple(sorted(U))
                if tu not in DSPACE:
                    DSPACE.add(tu)
                    if len(tu) <= 3:
                        nxt.append(tu)
    frontier = nxt
    say('closure round %d -> |D-space| = %d, new small sets = %d' % (rounds, len(DSPACE), len(nxt)))
layer = {}
for D in DSPACE:
    layer[len(D)] = layer.get(len(D), 0) + 1
ck('D2-deletion-space-size',
   [layer.get(k, 0) for k in (1, 2, 3, 4)] == [70, 3567, 143112, 4948069]
   and len(DSPACE) == 5094818,
   '|D| = 1,2,3,4: %s ; 70 + 3567 + 143112 + 4948069 = %d'
   % ([layer.get(k, 0) for k in (1, 2, 3, 4)], len(DSPACE)))
ck('D3-reduction-factor', comb(2969, 4) == 3231108527626,
   'C(2969,4) = %d unrestricted 4-subsets against %d in the D-space, a factor %.3e'
   % (comb(2969, 4), layer[4], comb(2969, 4) / layer[4]))
# D4 saturation: one further closure round over EVERY member of size <= 3 adds nothing.  Together
# with Lemma 2 of the paper this is what makes the sweep exhaustive.
addl = 0
for D in DSPACE:
    if len(D) > 3:
        continue
    sD = set(D)
    for B in partners(D):
        U = sD | set(B)
        if len(U) <= 4 and tuple(sorted(U)) not in DSPACE:
            addl += 1
ck('D4-deletion-space-is-saturated', addl == 0,
   'adjoining any one blocker set to any of the %d members of size <= 3 produces nothing outside '
   'the family, so it is the least fixed point' % (layer[1] + layer[2] + layer[3]))

# ===========================================================================
# E.  THE SWEEP
# ===========================================================================
stamp('E. the sweep: exact maximum clique in Cand(D) for every D')


def maxclique(masks):
    """Exact maximum clique size in the compatibility graph (u ~ v iff |u cap v| <= 7).
    Greedy-bound branch and bound; the inputs here have at most 6 vertices."""
    m = len(masks)
    adj = [0] * m
    for i in range(m):
        for j in range(i + 1, m):
            if popc(masks[i] & masks[j]) <= 7:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
    best = [0]

    def ext(sz, cand):
        if sz > best[0]:
            best[0] = sz
        if sz + bin(cand).count('1') <= best[0]:
            return
        c = cand
        while c:
            b = c & -c
            v = b.bit_length() - 1
            c ^= b
            ext(sz + 1, cand & adj[v] & ~((b << 1) - 1))

    ext(0, (1 << m) - 1)
    return best[0]


getbs = bs.get
comb_it = itertools.combinations
hist_cand = {}
hist_clique = {}
hist_margin = {}
tested = 0
live = 0
witness = None
for D in DSPACE:
    k = len(D)
    pool = []
    for r in range(1, k + 1):
        for T in comb_it(D, r):
            v = getbs(T)
            if v:
                pool += v
    tested += 1
    if k == 4:
        hist_cand[len(pool)] = hist_cand.get(len(pool), 0) + 1
    if len(pool) > k:
        live += 1
        cs = maxclique(pool)
        hist_clique[(k, cs)] = hist_clique.get((k, cs), 0) + 1
        hist_margin[cs - k] = hist_margin.get(cs - k, 0) + 1
        if cs > k:
            witness = (D, pool, cs)
            break
    else:
        # pruned: max clique <= |Cand(D)| <= |D|, so no witness here.  The margin recorded is the
        # UPPER BOUND min(|Cand(D)|,|D|) - |D|, not a measured clique size.
        hist_margin[min(len(pool), k) - k] = hist_margin.get(min(len(pool), k) - k, 0) + 1
ck('E1-every-set-was-tested', tested == 5094818 and witness is None,
   'swept %d of %d deletion sets, no shard, no cap, no early exit; %d had |Cand(D)| > |D| and '
   'received a full exact maximum-clique search' % (tested, len(DSPACE), live))
ck('E2-live-count', live == 14447, '%d deletion sets passed the necessary count filter' % live)
ck('E3-candidate-histogram',
   [hist_cand.get(j, 0) for j in range(1, 7)] == [17027, 936347, 2791429, 1189193, 13911, 162]
   and sum(hist_cand.values()) == 4948069 and max(hist_cand) == 6,
   '|Cand(D)| = 1..6 over the 4948069 four-sets: %s ; sums to %d ; max |Cand(D)| = %d ; only '
   '13911 + 162 = %d four-sets carry as many as 5 candidates'
   % ([hist_cand.get(j, 0) for j in range(1, 7)], sum(hist_cand.values()), max(hist_cand),
      hist_cand[5] + hist_cand[6]))
WANT_CLIQUE = {(2, 1): 2, (2, 2): 3, (3, 1): 7, (3, 2): 154, (3, 3): 208,
               (4, 1): 16, (4, 2): 753, (4, 3): 6179, (4, 4): 7125}
ck('E4-clique-histogram', hist_clique == WANT_CLIQUE and sum(hist_clique.values()) == 14447,
   '(|D|, max clique) over the live sets: %s ; sums to %d'
   % (sorted(hist_clique.items()), sum(hist_clique.values())))
ck('E5-the-two-tallies-cross-agree',
   sum(v for (k, c), v in hist_clique.items() if k == 4) == hist_cand[5] + hist_cand[6],
   'live four-sets counted by the clique tally = %d = %d = four-sets with |Cand(D)| >= 5 counted '
   'by the candidate tally; two independently accumulated histograms agree'
   % (sum(v for (k, c), v in hist_clique.items() if k == 4), hist_cand[5] + hist_cand[6]))
ck('E6-THEOREM-no-clique-ever-exceeds-D',
   all(c <= k for (k, c) in hist_clique) and witness is None,
   'for every one of the %d deletion sets the maximum set of pairwise-compatible insertable '
   'candidates has size at most |D|; there is no depth-<=4 exchange with |S| > |D|' % tested)
ck('E7-margin-histogram',
   [hist_margin.get(j, 0) for j in (-3, -2, -1, 0)] == [17043, 942405, 2874164, 1261206]
   and sum(hist_margin.values()) == 5094818
   and all(j <= 0 for j in hist_margin),
   'margin -3/-2/-1/0 = %s, sums to %d, and the "+1 or more" bucket is EMPTY'
   % ([hist_margin.get(j, 0) for j in (-3, -2, -1, 0)], sum(hist_margin.values())))

# ===========================================================================
# F.  ATTAINMENT: THE BREAK-EVEN EXCHANGES
# ===========================================================================
stamp('F. attainment')
D4 = EXHIBIT_D
ck('F1-exhibited-deletion-set-is-in-the-space',
   tuple(sorted(D4)) in DSPACE
   and tuple(hx(W[i]) for i in D4) == EXHIBIT_D_WORDS,
   'D = %s = %s' % (list(D4), list(EXHIBIT_D_WORDS)))
pool = []
for r in range(1, 5):
    for T in comb_it(D4, r):
        pool += bs.get(T, [])
ck('F2-Cand-of-the-exhibited-set',
   tuple(sorted(hx(m) for m in pool)) == tuple(sorted(EXHIBIT_CAND))
   and all(bmap[int(w, 16)] == EXHIBIT_BLOCKERS[w] for w in EXHIBIT_CAND)
   and all(popc(int(w, 16)) == 10 and int(w, 16) not in WSET for w in EXHIBIT_CAND),
   '|Cand(D)| = 5 > 4 = |D|, so this D is one of the 14447 live sets; each of %s has weight 10, '
   'is not in C_0, and has all its blockers inside D: %s'
   % (list(EXHIBIT_CAND), {w: list(EXHIBIT_BLOCKERS[w]) for w in EXHIBIT_CAND}))
pairs = []
for a, b in itertools.combinations(sorted(EXHIBIT_CAND), 2):
    pairs.append(((a, b), popc(int(a, 16) & int(b, 16))))
for (a, b), p in pairs:
    say('|0x%s cap 0x%s| = %d%s' % (a, b, p, '   CONFLICT' if p >= 8 else ''))
ck('F3-exactly-two-conflicting-pairs',
   tuple((pr, p) for pr, p in pairs if p >= 8) == EXHIBIT_CONFLICTS,
   'the only conflicts inside Cand(D) are |0x423E0D cap 0x4A3705| = 8 and '
   '|0x4A3705 cap 0x4E3305| = 9')
best = 0
argmax = []
CL = sorted(EXHIBIT_CAND)
for r in range(1, 6):
    for sel in itertools.combinations(CL, r):
        if all(popc(int(u, 16) & int(v, 16)) <= 7 for u, v in itertools.combinations(sel, 2)):
            if r > best:
                best = r
                argmax = [sel]
            elif r == best:
                argmax.append(sel)
ck('F4-maximum-clique-is-4-and-unique',
   best == 4 and len(argmax) == 1 and tuple(sorted(argmax[0])) == tuple(sorted(EXHIBIT_CLIQUE)),
   'brute force over all 31 non-empty subsets of Cand(D): the largest pairwise-compatible set has '
   'size 4, is unique, and equals %s; there is NO 5-clique, so the margin at this D is exactly 0'
   % list(EXHIBIT_CLIQUE))
ck('F5-clique-pairwise-intersections',
   tuple(popc(int(u, 16) & int(v, 16))
         for u, v in itertools.combinations(EXHIBIT_CLIQUE, 2)) == EXHIBIT_CLIQUE_INTS,
   'pairwise intersections inside the maximum clique, in order 01 02 03 12 13 23: %s'
   % list(EXHIBIT_CLIQUE_INTS))


def exchange_is_valid(delete_idx, insert_hex):
    """Full O(n^2) recheck that (C_0 \\ D) u S is a (23,6,10) constant-weight code."""
    drop = set(delete_idx)
    ws = [W[i] for i in range(N) if i not in drop] + [int(h, 16) for h in insert_hex]
    if len(set(ws)) != len(ws):
        return None
    if set(popc(v) for v in ws) != {10}:
        return None
    worst = 0
    for i in range(len(ws)):
        wi = ws[i]
        for j in range(i + 1, len(ws)):
            p = POP[(wi & ws[j]) & 0xFFFF] + POP[(wi & ws[j]) >> 16]
            if p > worst:
                worst = p
    return (len(ws), worst, ws)


EXCH = (
    (4, (27, 31, 585, 1333), ('214D69', '0B630D', '423E0D', '4E3305')),
    (3, (27, 585, 1333), ('214D69', '423E0D', '4E3305')),
    (2, (585, 1333), ('423E0D', '4E3305')),
    (1, (27,), ('214D69',)),
)
new_codes = []
for tag, Dk, Sk in EXCH:
    r = exchange_is_valid(Dk, Sk)
    ok = r is not None and r[0] == 2969 and r[1] <= 7
    new_codes.append(r[2] if r else None)
    ck('F6-break-even-exchange-at-k-%d' % tag, ok,
       'delete %s, insert %s -> a valid (23,6,10) code with %d words and worst pairwise '
       'intersection %d; gain = |S| - |D| = %d exactly'
       % (list(Dk), list(Sk), r[0] if r else -1, r[1] if r else -1, len(Sk) - len(Dk)))
ck('F7-k-1-attainment-is-universal',
   len(BY_SIZE[1]) == 70 and all(len(bs[t]) == 1 for t in BY_SIZE[1]),
   'all %d singleton deletion sets in the D-space have |Cand(D)| = 1, so every one of them gives a '
   '1-for-1 break-even exchange and none gives a 2-for-1' % len(BY_SIZE[1]))
ck('F8-the-refilled-code-is-a-different-code',
   new_codes[0] is not None and set(new_codes[0]) != WSET
   and len(set(new_codes[0]) & WSET) == 2965,
   'the k = 4 exchange yields a 2969-word (23,6,10) code sharing 2965 words with C_0, so it is a '
   'genuinely different code of the same size, not a relabelling')
# ===========================================================================
# SCOPE AND VERDICT
# ===========================================================================
sys.stdout.write('\n')
stamp('scope')
for line in (
    'NOT RE-RUN: Section XII of Brouwer-Shearer-Sloane-Smith 1990, which the source paper '
    'describes as certifying k-optimality of this move for k = 2..5, is paywalled and was never '
    'read. If it certifies k-optimality for the 2969-word code specifically, the |D| = 4 layer of '
    'Theorem 1 is not new. Nothing in this program can settle that.',
    'NOT RE-RUN: the |D| <= 3 layer of the sweep (70 + 3567 + 143112 = 146749 sets) RE-PROVES a '
    'theorem already published in the source paper; only the |D| = 4 layer is new, and within it '
    'the source\'s own depth-unbounded theorem already excludes every S all of whose members carry '
    'at most two blockers.',
    'NOT RE-RUN: in check E7 the -2, -1 and 0 buckets are UPPER-BOUND-DERIVED for the 5080371 '
    'sets with |Cand(D)| <= |D|: the recorded margin there is min(|Cand(D)|,|D|) - |D| and no '
    'clique search was run, so a set counted in one of those buckets may truly sit lower. Only the '
    '-3 bucket and the 14447 live margins are measured. The conclusion "+1 or more is empty" is '
    'untouched, a clique being unable to exceed |Cand(D)|.',
    'NOT RE-RUN: depth 5 and above; any code other than C_0; and any lower bound on A(23,6,10). '
    'C_0 is superseded -- the source paper itself proves A(23,6,10) >= 2979 and Brouwer\'s table '
    'reads 2992 -- so nothing here moves a bound.',
    'NOT RE-RUN: the two forced-positive controls on modified codes (2968 and 2961 words) that '
    'were used during the search to show the decision procedure can return YES. Their role is '
    'replaced here by check F6, which exhibits genuine break-even exchanges at every k in 1..4, '
    'and by check B7, which compares two disjoint censuses.',
):
    sys.stdout.write(line + '\n')
sys.stdout.write('\n')
if _nf[0]:
    sys.stdout.write('VERDICT: %d CHECK(S) FAILED\n' % _nf[0])
    sys.exit(1)
sys.stdout.write('VERDICT: ALL %d CHECKS PASS\n' % _np[0])
sys.exit(0)
