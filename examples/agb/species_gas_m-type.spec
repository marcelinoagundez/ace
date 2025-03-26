!
!elements
!--------
34  2            ! No of elements, (1/2) abundances (relative/log relative to H=12)
 H  12.00        ! element, abundance
 He 10.93        ! element, abundance
 Li -0.30        ! element, abundance
 Be  1.38        ! element, abundance
 B   2.70        ! element, abundance
 C   8.43        ! element, abundance
 N   7.83        ! element, abundance
 O   8.69        ! element, abundance
 F   4.48        ! element, abundance
 Ne  7.93        ! element, abundance
 Na  6.24        ! element, abundance
 Mg  7.60        ! element, abundance
 Al  6.45        ! element, abundance
 Si  7.51        ! element, abundance
 P   5.41        ! element, abundance
 S   7.12        ! element, abundance
 Cl  5.50        ! element, abundance
 Ar  6.40        ! element, abundance
 K   5.03        ! element, abundance
 Ca  6.34        ! element, abundance
 Sc  3.15        ! element, abundance
 Ti  4.95        ! element, abundance
 V   3.93        ! element, abundance
 Cr  5.64        ! element, abundance
 Mn  5.43        ! element, abundance
 Fe  7.50        ! element, abundance
 Co  4.99        ! element, abundance
 Ni  6.22        ! element, abundance
 Cu  4.19        ! element, abundance
 Zn  4.56        ! element, abundance
 Rb  2.52        ! element, abundance
 Sr  2.87        ! element, abundance
 Zr  2.58        ! element, abundance
 Ba  2.18        ! element, abundance
!
! .................... Gas-phase neutral species
!
! ... Hydrogen : 1
H
H2
! ... Helium : 2
He
! ... Carbon : 6
C 
C2
C2(3)              ! triplet
C3
C4
C5
C6
C7
C8
C9
C10
C11
C12
! ... C,H          ! still various isomers to implement 
!                  ! to be done when building chemical network
!                  ! for the moment only the most stable is included
CH
CH2
CH2(1)             ! singlet
CH3
CH4
C2H
C2H2
C2H3
C2H4
C2H5
C2H6
C3H
c-C3H2             ! cyclopropenyldiene
H2C3               ! propadienylidene
HC3H               ! propynylidene
C3H3               ! propargyl
P-C3H4             ! propyne
A-C3H4             ! allene
C3H5               ! allyl
C3H6               ! propylene
n-C3H7             ! n-propyl
i-C3H7             ! isopropyl
C3H8               ! propane
C4H
C4H2               ! butadiyne
C4H3               ! 1,2,3-butatriene-4-yl; i-1-butene-3yne-2-yl (same energy)
C4H4               ! 1-butene-3yne CH2=CH-CCH
C4H5               ! 2-Butayn-1-yl CH3-CC-CH2*
C4H6               ! 1,3-butadiene
C4H7               ! trans-1-methylallyl CH2=CHC*HCH3
C4H8               ! isobuten CH3C(CH3)=CH2
C4H9               ! t-butyl
C4H10              ! isobutane
C5H
C5H2               ! *HC=C=C=C=CH*
C5H3               ! 1,4-pentadiyne-3-yl HCCCH*CCH
C5H4               ! 1,3 pentadiyne HCC-CC-CH3
C6H
C6H2               ! hexatriyne
C7H
C8H
C8H2
C9H
C10H
C10H2
C11H
C12H
C12H2
! ... Oxygen : 8
O 
O(1D) 
O2
O3
! ... H,O
OH
H2O
HO2
H2O2
! ... C,H,O
CO
HCO
H2CO
CH2OH
CH3O
CH3OH
C2O
HCCO               ! ketenyl
H2CCO              ! ketene
CH3CO              ! acetyl
CH2CHO             ! vinyloxy
CH3CHO             ! acetaldehyde
CH2CHOH            ! vinyl alcohol
c-C2H4O            ! oxirane, ethylene oxide
C2H5O              ! ethoxy
CH2CH2OH           ! 
CH3CHOH            !
CH3OCH2            ! methyl methoxy
C2H5OH             ! ethanol
CH3OCH3            ! dimethyl ether
C3H3O              ! acrolein CH2=CH*CO
CO2
HOCO               ! trans-hydroxyformyl
HCOOH              ! cis-formic acid
CH3O2              ! methyl peroxy
CH2(OH)2           ! methanediol
CH3OOH             ! peroxymethane
(CHO)2             ! glyoxal O=CHCH=O
CH3COO             ! acetyloxyl (acetic acid radical)
COOCH3             ! methyl formate radical
CH3COOH            ! acetic acid
HCOOCH3            ! methyl formate
C3O2
! ... Nitrogen : 7
N 
N2
N3
! ... N,H
NH
NH2
NH3
N2H                ! diazenyl
N2H2               ! diazine (trans/cis-HNNH equilibrium)
N2H3               ! hydrazine radical
N2H4               ! hydrazine
! ... C,N,H
CN
HCN
HNC
H2CN               ! H2C=N*
HCNH               ! trans H*C=NH
CH2NH              ! methaneimine
CH2NH2             ! methylene amine
CH3NH2             ! methyl amine
C2N
HCCN
CH2CN              ! cyanomethyl
CH3CN              ! methyl cyanide
CH3NC              ! methyl isocyanide
!C3N               ! <----- no therm data found
!HC3N              ! <----- no therm data found
!HC5N              ! <----- no therm data found
!HC7N              ! <----- no therm data found
NCN
CNN
HNCN
HCNN
C2N2               ! cyanogen
CNCN               ! isocyanogen
C4N2
! ... C,N,H,O
NO
HNO
HNOH
NH2OH              ! hydroxylamine
NO2
HNO2               ! nitrous acid (trans/cis-HONO equilibrium)
NO3
HNO3               ! nitric acide
N2O
N2O3
N2O4
NCO                ! isocyanate
HNCO               ! isocyanic acid
HOCN               ! cyanic acid
HCNO               ! fulminic acid
! ... Sulfur : 16
S 
S2
S3
S4
S5
S6
S7
S8
! ... H,S
SH
H2S
HS2
H2S2
! ... H,C,S
CS
HCS
H2CS               ! thioformaldehyde
CH3S
CH3SH              ! methyl mercaptane
HCCS               ! ethyne thiol radical
CH2CHSH            ! vinyl mercaptan
C2H5S              ! ethyl thio radical
C2H5SH             ! ethane thiol
CH3SCH3            ! dimethylsulfide
CS2
C2S2
C3S2
! ... H,C,S,O
SO
HSO
SOH
SO2
HSO2
SO3
HSO3
H2SO4              ! sulfuric acid
S2O
OCS
C3SO
! ... H,C,S,O,N
NS
! ... Silicon : 14
Si
Si2
Si3
! ... Si,H
SiH
SiH2
SiH3
SiH4
Si2H4
Si2H5
Si2H6
! ... Si,C,H
SiC
SiC2
SiC3
SiC4
SiC5
Si2C
Si2C2
Si2C3
Si2C4
Si3C
Si3C2
Si3C3
Si4C
Si4C2
Si5C
SiCH
SiCH2
CH2SiH2
CH3SiH2
CH3SiH3
! ... Si,C,H,O
SiO
SiO2
! ... Si,C,H,O,N
SiN
SiNH
Si2N
! ... Si,C,H,O,N,S
SiS
SiS2
! ...Phosphorus : 15
P 
P2
P3
P4
! ...P,H
PH
PH2
PH3
P2H
P2H2
P2H4               ! byphosphine
! ...P,H,C
CP
HCP
CH2PH              ! methylenephosphine
CH3PH2             ! methylphosphine
! ...P,H,C,O
PO
HPO
H2PO               ! HPOH
H3PO               ! H3PO
PO2
HPO2               ! HOPO
PO3
HPO3               ! HOPO2
H3PO3              ! phosphonic acid
H3PO4              ! phosphoric acid
P2O3
P2O4
P2O5
P3O6
P4O6
P4O7
P4O8
P4O9
P4O10
! ...P,H,C,O,N
PN
! ...P,H,C,O,S
PS
! ...Chlorine : 17 ! larger organohalogens in Burcat's files
Cl
Cl2
HCl
! ...Cl,H,C
CCl
CCl2
CCl3
CCl4
CHCl
CHCl2
CHCl3              ! chloroform
CH2Cl
CH2Cl2
CH3Cl
C2Cl
C2Cl2
C2Cl3
C2Cl4
C2Cl6
C2HCl
C2HCl3
C2H2Cl
C2H2Cl2
C2H3Cl
C2H4Cl
C2H5Cl
! ...Cl,H,C,O
ClO
ClO2
ClO3
ClO4
Cl2O
Cl2O2
HOCl
HClO2
HClO3
HClO4
COCl
COCl2
CHClO
CCl3OH
CH3COCl
C2H3ClO2
! ...Cl,H,C,O,N
NClH2
NCl2H
NCl3
ClCN
NOCl
NO2Cl
! ...Cl,H,C,O,S
SCl
SCl2
S2Cl
S2Cl2
SO2Cl2
! ...Cl,H,C,O,Si
SiCl
SiCl2
SiCl3
SiCl4
SiHCl
SiHCl3
SiH2Cl2
SiH3Cl
! ...Cl,H,C,O,P
PCl
PCl2
PCl3
PCl5
POCl3
! ...Fluorine : 9
F 
F2
! ...F,H
HF
H2F2
H3F3
H4F4
H5F5
H6F6
H7F7
! ...F,H,C
CF
CHF
CH2F
CH3F
CF2
CHF2
CH2F2
CF3
CHF3
CF4
C2F
C2HF
C2F2
C2F3
C2F4
C2HF3
C2H2F2             ! 1,1 difluoroethylene FC-1132a
C2H3F
C2F6
C3F
! ...F,H,C,O
FO
HOF
FO2                ! FOO
F2O
F2O2               ! FOOF
FCO
HFCO
COF2
! ...F,H,C,O,N
NF
NHF
NH2F
NF2
NHF2
NF3
N2F2
N2F4
FCN
NOF
NOF3
NO2F
NO3F
! ...F,H,C,O,S
SF
SF2
SF3
SF4
SF5
SF6
S2F2               ! thiothionyl fluoride SSF2
SOF2
SO2F2
HSO3F
! ...F,H,C,O,Si
SiF
SiF2
SiF3
SiF4
SiHF
SiH3F
SiH2F2
SiHF3
Si2F6
! ...F,H,C,O,P
PF
PF2
PF3
PF5
POF3
! ...F,H,C,O,S,Si,P,Cl
ClF
ClF3
ClF5
CFCl
CHFCl
CH2FCl
CFCl2
CHFCl2
CFCl3
CF2Cl
CHF2Cl
CF2Cl2
CF3Cl
C2FCl
C2H2FCl
C2HFCl2
C2FCl3
C2HClF2            ! cis-CHF=CFCl
C2F2Cl2
C2F3Cl
COFCl
SO2FCl
SiFCl
PFCl
PFCl2
PFCl4
PF2Cl
PF2Cl3
PF3Cl2
PF4Cl
POFCl2
POF2Cl
! ... Sodium : 11
Na
Na2
NaH
NaO
Na2O
Na2O2
NaOH
Na2O2H2
NaCN
Na2C2N2
NaNO2
NaNO3
Na2SO4
NaCl
Na2Cl2
Na3Cl3
NaF
Na2F2
Na3F3
! ... Magnesium : 12
Mg
Mg2
MgH
MgO
MgOH
Mg(OH)2
MgN
MgS
MgCl
MgCl2
Mg2Cl4
MgF
MgF2
MgClF
Mg2F4
! ... Aluminium : 13
Al
Al2
AlH
AlH2
AlH3
AlC
AlC2
Al2C2
AlO
AlO2
Al2O
Al2O2
Al2O3
AlOH
Al(OH)2
Al(OH)3
HAlO
HAlO2
AlN
AlS
AlS2
Al2S
Al2S2
AlCl
AlCl2
AlCl3
AlHCl
AlHCl2
Al2Cl6
AlOCl
AlOCl2
AlOHCl
AlOHCl2
AlF
AlF2
AlF3
Al2F6
AlHF
AlHF2
AlH2F
AlOF
AlOF2
AlOHF
AlOHF2
Al(OH)2F
Al(OH)2Cl
AlFCl
AlFCl2
AlF2Cl
AlHFCl
AlH2Cl
LiAlF4
NaAlF4
! ... Potassium : 19
K
K2
KH
KO
K2O
K2O2
KOH
K2(OH)2
K2CO3
KCN
K2C2N2
KNO2
KNO3
K2SO4
KCl
K2Cl2
KF
K2F2
NaK
KAlF4
KLi
! ... Calcium : 20
Ca
Ca2
CaH
CaO
CaOH
Ca(OH)2
CaS
CaCl
CaCl2
CaF
CaF2
! ... Titanium : 22
Ti
TiH
TiC
TiC2
TiC3
TiC4
Ti2C
Ti2C2
Ti2C3
Ti2C4
Ti3C
Ti3C2
Ti3C3
Ti3C4
Ti4C
Ti4C2
Ti4C3
Ti4C4
Ti3C8
Ti4C8
Ti6C13
Ti7C13
Ti8C12
Ti9C15
Ti13C22
TiO
TiO2
TiN
TiS
TiCl
TiCl2
TiCl3
TiCl4
TiOCl
TiOCl2
TiF
TiF2
TiF3
TiF4
TiOF
TiOF2
! ... Chromium : 24
Cr
Cr2
CrH
CrO
CrO2
CrO3
Cr(OH)6
CrN
CrS
CrCl
CrCl2
CrCl3
CrCl4
CrCl5
CrCl6
CrO2Cl2
! ... Manganese : 25
Mn
MnH
MnO
MnF
MnCl
! ... Iron : 26
Fe
FeH
FeO
FeO2
Fe(OH)2
Fe(CO)5
FeS
FeCl
FeCl2
FeCl3
Fe2Cl4
Fe2Cl6
FeF
FeF2
FeF3
! ... Cobalt : 27
Co
CoH
CoCl
CoCl2
CoCl3
Co2Cl4
CoF2
! ... Nickel : 28
Ni
NiH
NiO
NiS
NiCl
NiCl2
NiF
! ... Copper : 29
Cu
Cu2
CuH
CuO
CuOH
CuS
CuCl
Cu3Cl3
CuF
CuF2
! ... Zinc : 30
Zn
ZnH
CH3Zn
C2H6Zn
ZnO
ZnS
ZnF
ZnCl
ZnCl2
! ... Lithium : 3
Li
Li2
LiH
LiO
LiOH
Li2O
Li2O2
Li2(OH)2
LiN
LiON
LiNO2
LiNO3
Li2SO4
LiCl
Li2Cl2
Li3Cl3
LiOCl
LiF
Li2F2
Li3F3
Li2FCl
LiOF
LiNa
LiONa
! ... Berilium : 4
Be
Be2
BeH
BeH2
BeO
Be2O
Be2O2
Be3O3
Be4O4
BeOH
Be(OH)2
BeN
BeS
BeCl
BeCl2
Be2Cl4
BeF
BeF2
Be2F4
Be2OF2
! ... Boron : 5
B
B2
BH
BH2
BH3
BH4
BH5
B2H
B2H2
B2H3
B2H4
B2H5
B2H6
B3H7
B3H9
B4H4
B4H10
B4H12
B5H9
BC
BC2
B2C
BO
BO2
B2O
B2O2
B2O3
BOH
HBO
HBOH
H2BOH
HBO2
B(OH)2
HB(OH)2
H3BO3
B2(OH)4
H3B3O3
H3B3O6
BN
BH3NH3
B3N3H6
BS
BS2
B2S
B2S2
B2S3
HBS
BCl
BCl2
BCl3
BHCl
BHCl2
BH2Cl
BOCl
BOCl2
B2Cl4
BClOH
BCl2OH
BCl(OH)2
BF
BF2
BF3
B2F4
BHF
BHF2
BH2F
BOF
BOF2
B3O3F3
BFOH
BF2OH
BF(OH)2
BFCl
BFCl2
BF2Cl
BHFCl
B3O3F2Cl
B3O3FCl2
B3O3Cl3
LiBO2
NaBO2
KBO2
MgB2
AlB2
! ... Scandium : 21
Sc
ScO
ScO2
Sc2O
Sc2O2
ScS
ScF
ScCl
! ... Vanadium : 23
V
VO
VO2
V2O4
V4O10
VN
VCl4
! ... Rubidium : 37
Rb
RbH
RbO
Rb2O
Rb2O2
RbOH
Rb2O2H2
RbNO2
RbNO3
Rb2SO4
RbCl
Rb2Cl2
RbF
Rb2F2
RbNa
RbK
RbLi
RbBO2
! ... Srontium : 38
Sr
Sr2
SrH
SrO
SrOH
Sr(OH)2
SrS
SrCl
SrCl2
SrF
SrF2
! ... Zirconium : 40
Zr
ZrN
ZrO
ZrO2
ZrCl2
ZrCl4
ZrF
ZrF2
ZrF4
! ... Barium : 56
Ba
Ba2
BaH
BaO
BaOH
Ba(OH)2
BaS
BaCl
BaCl2
BaF
BaF2
! ... Neon : 10
Ne
NeH
! ... Argon : 18
Ar
ArH
!
! .................... Gas-phase ion species
!
e-
H+
H-
H2+
H3+
He+
HeH+
C+
C-
O+
O-
N+
N-
S+
S-
Si+
Si-
P+
P-
Cl+
Cl-
F+
F-
Na+
Na-
Mg+
Al+
Al-
K+
K-
Ca+
Ti+
Ti-
Cr+
Cr-
Mn+
Fe+
Fe-
Co+
Co-
Ni+
Ni-
Cu+
Cu-
Zn+
Li+
Li-
Be+
B+
B-
Sc+
Sc-
V+
V-
Rb+
Rb-
Sr+
Zr+
Zr-
Ba+
Ne+
NeH+
Ar+
ArH+
