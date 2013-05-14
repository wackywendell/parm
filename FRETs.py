from numpy import sqrt, mean, std, array, zeros
from os.path import expanduser
mydir = expanduser('~/idp/')

Forsterdist = 54

aSijs = [(54, 72), (72, 92), (9, 33), (54, 92), (92, 130), (33, 72),
            (9, 54), (72, 130), (9, 72), (54, 130), (33, 130), (9, 130)]
aSETs = [0.85, 0.85, 0.86, 0.66, 0.65, 0.72, 0.7, 0.55, 0.64, 0.5, 0.53, 0.36]
aSerrs = [0.0437, 0.0451, 0.0503, 0.0674, 0.0624, 0.0643, 0.0608,
                    0.0621, 0.0658, 0.0627, 0.0608, 0.0529]
aSseq = "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA"
aSpdb = mydir + 'pdb/aS.pdb'
aS = aSijs, aSETs, aSerrs

tauijs = [(291, 322), (291, 354), (354, 433), (103, 184),
                        (17, 103), (184, 291), (244, 354), (322, 433),
                        (291, 433), (103, 291), (17, 291), (17, 433)]
tauETs = [0.81, 0.61, 0.62, 0.7, 0.17, 0.36,
                                    0.37, 0.51, 0.42, 0.33, 0.42, 0.22]
tauerrs = [0.054, 0.069, 0.068, 0.067, 0.032, 0.057, 
                                0.058, 0.067, 0.062, 0.055, 0.062, 0.04]
tauseq = "MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPGSETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDKKAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL"
taupdb = mydir + 'pdb/tau_phyre2.pdb'
tau = tauijs, tauETs, tauerrs

bSijs = [(9,33), (35,59), (59,83), (102, 126), (18, 126)]
bSETs = [0.831, 0.851, 0.863, 0.607, 0.234]
bSerrs = [0.055, 0.052, 0.048, 0.074, 0.047]

bSseq = "MDVFMKGLSMAKEGVVAAAEKTKQGVTEAAEKTKEGVLYVGSKTREGVVQGVASVAEKTKEQASHLGGAVFSGAGNIAAATGLVKREEFPTDLKPEEVAQEAAEEPLIEPLMEPEGESYEDPPQEEYQEYEPEA"
bSpdb = mydir + 'pdb/beta-synuclein-phyre2.pdb'
bS = bSijs, bSETs, bSerrs

gSijs = [(9,33), (35,59), (59,83), (92, 117), (9, 117)]
gSETs = [0.88, 0.825, 0.791, 0.799, 0.343]
gSerrs = [0.044,0.055,0.059,0.060,0.057]

gSseq = "MDVFKKGFSIAKEGVVGAVEKTKQGVTEAAEKTKEGVMYVGAKTKENVVQSVTSVAEKTKEQANAVSEAVVSSVNTVATKTVEEAENIAVTSGVVRKEDLRPSAPQQEGEASKEKEEVAEEAQSGGD"
gSpdb = mydir + 'pdb/gamma-synuclein-phyre2.pdb'
gS = gSijs, gSETs, gSerrs

def ETeffrsq(ETs, expETs):
    ETs, expETs = array(ETs), array(expETs)
    dists = (ETs -expETS)**2
    return sqrt(mean(dists))
    
def ETeffrsqerrs(ETs, expETs, simETerrs=None, expETerrs=None):
    """\Delta=\sqrt{\left\langle \left(x_i-y_i\right)^2\right\rangle}
    \sigma_{\Delta}^2=\sum_{i=1}^{N}\left(\frac{\partial\Delta}{\partial x_{i}}\right)^{2}\sigma_{x_{i}}^{2}+\left(\frac{\partial\Delta}{\partial y_{i}}\right)^{2}\sigma_{y_{i}}^{2}
    \frac{\partial\Delta}{\partial x_i} = \frac{x_i-y_i}{N\Delta}"""
    ETs, expETs = array(ETs), array(expETs)
    N = len(ETs)
    dxs = ETs - expETs
    sigxs = array(simETerrs) if simETerrs is not None else zeros((N,))
    sigys = array(expETerrs) if expETerrs is not None else zeros((N,))
    
    D = sqrt(mean(dxs**2))
    
    dDdxi = dxs/D/N
    Derr = sqrt(sum((dDdxi * sigys)**2 + (dDdxi * sigxs)**2))
    
    return D, Derr

aSseq3 = "\
MET ASP VAL PHE MET LYS GLY LEU SER LYS ALA LYS GLU \
GLY VAL VAL ALA ALA ALA GLU LYS THR LYS GLN GLY VAL \
ALA GLU ALA ALA GLY LYS THR LYS GLU GLY VAL LEU TYR \
VAL GLY SER LYS THR LYS GLU GLY VAL VAL HIS GLY VAL \
ALA THR VAL ALA GLU LYS THR LYS GLU GLN VAL THR ASN \
VAL GLY GLY ALA VAL VAL THR GLY VAL THR ALA VAL ALA \
GLN LYS THR VAL GLU GLY ALA GLY SER ILE ALA ALA ALA \
THR GLY PHE VAL LYS LYS ASP GLN LEU GLY LYS ASN GLU \
GLU GLY ALA PRO GLN GLU GLY ILE LEU GLU ASP MET PRO \
VAL ASP PRO ASP ASN GLU ALA TYR GLU MET PRO SER GLU \
GLU GLY TYR GLN ASP TYR GLU PRO GLU ALA".split(' ')
