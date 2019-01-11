################################################################################################################
targets, factors
ESE2, (!NFkB & !Snai2 & !ESE2) | (!NFkB & !Snai2 & ESE2) | (!NFkB & Snai2 & ESE2) | (NFkB & !Snai2 & !ESE2) | (NFkB & !Snai2 & ESE2)
E2F, (!Rb & !p53 & !Snai2 & !Cyclin) | (!Rb & !p53 & !Snai2 & Cyclin)
Cyclin, (!ESE2 & !E2F & !p16 & !NFkB & !Snai2) | (!ESE2 & !E2F & !p16 & NFkB & !Snai2) | (!ESE2 & !E2F & !p16 & NFkB & Snai2) | (!ESE2 & E2F & !p16 & !NFkB & !Snai2) | (!ESE2 & E2F & !p16 & NFkB & !Snai2) | (!ESE2 & E2F & !p16 & NFkB & Snai2) | (ESE2 & !E2F & !p16 & !NFkB & !Snai2) | (ESE2 & !E2F & !p16 & NFkB & !Snai2) | (ESE2 & E2F & !p16 & !NFkB & !Snai2) | (ESE2 & E2F & !p16 & NFkB & !Snai2)
Snai2, (!ESE2 & !NFkB & !Snai2) | (!ESE2 & !NFkB & Snai2) | (!ESE2 & NFkB & !Snai2) | (!ESE2 & NFkB & Snai2) | (ESE2 & NFkB & Snai2)
NFkB, !(!ESE2 & !p16 & !Snai2 & !NFkB)
TELasa, (!Snai2 & !ESE2) | (Snai2 & !ESE2)
Rb, (!Cyclin & !p16 & p53) | (!Cyclin & p16 & !p53) | (!Cyclin & p16 & p53) | (Cyclin & !p16 & p53) | (Cyclin & p16 & !p53) | (Cyclin & p16 & p53)
p16, (!p16 & !E2F & p53 & !TELasa & !Snai2) | (!p16 & !E2F & p53 & !TELasa & Snai2) | (!p16 & !E2F & p53 & TELasa & !Snai2) | (!p16 & E2F & p53 & !TELasa & !Snai2) | (!p16 & E2F & p53 & !TELasa & Snai2) | (!p16 & E2F & p53 & TELasa & !Snai2) | (p16 & !E2F & !p53 & !TELasa & !Snai2) | (p16 & !E2F & p53 & !TELasa & !Snai2) | (p16 & !E2F & p53 & !TELasa & Snai2) | (p16 & !E2F & p53 & TELasa & !Snai2) | (p16 & E2F & !p53 & !TELasa & !Snai2) | (p16 & E2F & !p53 & !TELasa & Snai2) | (p16 & E2F & !p53 & TELasa & !Snai2) | (p16 & E2F & !p53 & TELasa & Snai2) | (p16 & E2F & p53 & !TELasa & !Snai2) | (p16 & E2F & p53 & !TELasa & Snai2) | (p16 & E2F & p53 & TELasa & !Snai2) | (p16 & E2F & p53 & TELasa & Snai2) | (p16 & !E2F & !p53 & TELasa & !Snai2)
p53, (!p53 & !NFkB & !TELasa & !p16 & !Snai2) | (!p53 & !NFkB & !TELasa & p16 & !Snai2) | (!p53 & NFkB & !TELasa & p16 & !Snai2) | (p53 & !NFkB & !TELasa & !p16 & !Snai2) | (p53 & !NFkB & !TELasa & p16 & !Snai2) | (p53 & NFkB & !TELasa & p16 & !Snai2)
################################################################################################################
targets, factors
AG, (! EMF1 & ! AP2 & ! TFL1) | (! EMF1 & ! AP1 & LFY) | (! EMF1 & ! AP2 & LFY) | (! EMF1 & ! TFL1 & LFY & (AG & SEP)) | (! EMF1 & (LFY & WUS))
AP1, (! AG & ! TFL1) | (FT & LFY & ! AG) | (FT & ! AG & ! PI) | (LFY & ! AG & ! PI) | (FT & ! AG & ! AP3) | (LFY & ! AG & ! AP3)
AP2, ! TFL1
AP3, (LFY & UFO) | (PI & SEP & AP3 & (AG | AP1))
EMF1, ! LFY
FT, ! EMF1
FUL, ! AP1 & ! TFL1
LFY, ! EMF1 | ! TFL1
PI, (LFY & (AG | AP3)) | (PI & SEP & AP3 & (AG | AP1))
SEP, LFY
TFL1, ! AP1 & (EMF1 & ! LFY)
WUS, WUS & (! AG | ! SEP)
################################################################################################################
targets, factors
Gene1, (!Gene1 & Gene8) | (Gene1 & !Gene8)
Gene2, (!Gene5 & !Gene3) | (Gene5 & Gene3)
Gene3, (!Gene3) | (!Gene8)
Gene4, (!Gene1 & Gene4)
Gene5, (Gene4 & !Gene5)
Gene6, (Gene4) | (Gene2)
Gene7, (!Gene8) | (!Gene7)
Gene8, (Gene5 & Gene1)
################################################################################################################

