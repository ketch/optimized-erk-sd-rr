# python module erk
#
# Description:
#
# Module to extract the ERK 3S* coefficients from the textfiles, as python lists,
# which can be used to configure the ERK schemes implemented in coolfluid
#
#
# Usage:
#
# 1) as module to be imported
# import erk
# erk_coeffs = erk.parse(filepath)
#
# 2) standalone script to test. It #prints the python lists
# python erk.py filepath



# import regular expressions module
import re
import os

# class that contains all coefficients needed for coolfluid
class CF3_ERK_Coeffs:
    order=0
    nb_stages=0
    c=[]*nb_stages
    beta=[0]*nb_stages
    gamma1=[0]*nb_stages
    gamma2=[0]*nb_stages
    gamma3=[0]*nb_stages
    delta=[0]*nb_stages
    
    def __str__(self):
        out =  "order     = "+str(self.order)+"\n"
        out += "nb_stages = "+str(self.nb_stages)+"\n"
        out += "c         = "+str(self.c)+"\n"
        out += "beta      = "+str(self.beta)+"\n"
        out += "gamma1    = "+str(self.gamma1)+"\n"
        out += "gamma2    = "+str(self.gamma2)+"\n"
        out += "gamma3    = "+str(self.gamma3)+"\n"
        out += "delta     = "+str(self.delta)+"\n"
        out += "cfl       = "+str(self.cfl())+"\n"
        out += "cfl/stage = "+str(self.cfl()/self.nb_stages)
        return out
    
    def tex(self):
        out = '\\resizebox{\\linewidth}{!}{\\begin{tabular}{r|r|r|r|r|r}\n'
        out += '\\multicolumn{1}{c|}{c}  &  \\multicolumn{1}{c|}{$\\beta$}  &  \\multicolumn{1}{c|}{$\\gamma_1$}  &  \\multicolumn{1}{c|}{$\\gamma_2$}  &  \\multicolumn{1}{c|}{$\\gamma_3$} &  \\multicolumn{1}{c}{$\\delta$} \\\\\n'
        out += '\\hline\n'
        for i in range(self.nb_stages) :
            out += '{c:23.16e}  & {b:23.16e}  & {g1:23.16e}  & {g2:23.16e}  & {g3:23.16e} & {d:23.16e} \\\\\n'.format(
                    c=self.c[i],
                    b=self.beta[i],
                    g1=self.gamma1[i],
                    g2=self.gamma2[i],
                    g3=self.gamma3[i],
                    d=self.delta[i])
        out += '\\hline \n \\end{tabular}}'
        return out
    
    # cfl number was extracted from Unconstrained SD coefficient files
    def cfl(self):
        if (self.order == 2):
            if (self.nb_stages == 3):
                m_cfl = 0.5876247771084309
            elif (self.nb_stages == 4):
                m_cfl = 0.8171781152486801
            elif (self.nb_stages == 5):
                m_cfl = 1.0345799848437309
            elif (self.nb_stages == 6):
                m_cfl = 1.2511042691767216
            elif (self.nb_stages == 7):
                m_cfl = 1.4647344406694174
            elif (self.nb_stages == 8):
                m_cfl = 1.6774455457925797
            elif (self.nb_stages == 9):
                m_cfl = 1.8892222177237272
            elif (self.nb_stages == 10):
                m_cfl = 2.1005510352551937
            elif (self.nb_stages == 11):
                m_cfl = 2.3115957109257579
            elif (self.nb_stages == 12):
                m_cfl = 2.5223528221249580
            elif (self.nb_stages == 13):
                m_cfl = 2.7329694945365191
            elif (self.nb_stages == 14):
                m_cfl = 2.9434285033494234
            elif (self.nb_stages == 15):
                m_cfl = 3.1539227580651641
            elif (self.nb_stages == 16):
                m_cfl = 3.3642986416816711
            elif (self.nb_stages == 17):
                m_cfl = 3.5746607976034284
            elif (self.nb_stages == 18):
                m_cfl = 3.7849829904735088
            elif (self.nb_stages == 19):
                m_cfl = 3.9953222870826721
            elif (self.nb_stages == 20):
                m_cfl = 4.2056231759488583
            else:
                print "no cfl number defined for ",self.nb_stages," rk stages"
        elif (self.order == 3):
            if (self.nb_stages == 4):
                m_cfl = 0.3345983475446701
            elif (self.nb_stages == 5):
                m_cfl = 0.4535968415439129
            elif (self.nb_stages == 6):
                m_cfl = 0.5728160217404366
            elif (self.nb_stages == 7):
                m_cfl = 0.6915322504937649
            elif (self.nb_stages == 8):
                m_cfl = 0.8083578944206238
            elif (self.nb_stages == 9):
                m_cfl = 0.9242995362728834
            elif (self.nb_stages == 10):
                m_cfl = 1.0390471667051315
            elif (self.nb_stages == 11):
                m_cfl = 1.1524990852922201
            elif (self.nb_stages == 12):
                m_cfl = 1.2652794085443020
            elif (self.nb_stages == 13):
                m_cfl = 1.3778314273804426
            elif (self.nb_stages == 14):
                m_cfl = 1.4892721455544233
            elif (self.nb_stages == 15):
                m_cfl = 1.6004448523744941
            elif (self.nb_stages == 16):
                m_cfl = 1.7116628587245941
            elif (self.nb_stages == 17):
                m_cfl = 1.8220776738598943
            elif (self.nb_stages == 18):
                m_cfl = 1.9325973931699991
            elif (self.nb_stages == 19):
                m_cfl = 2.0428037177771330
            elif (self.nb_stages == 20):
                m_cfl = 2.1530115976929665
            else:
                print "no cfl number defined for ",self.nb_stages," rk stages"
        elif (self.order == 4):
            if (self.nb_stages == 5):
                m_cfl = 0.2366644330322742
            elif (self.nb_stages == 6):
                m_cfl = 0.3020527213811874
            elif (self.nb_stages == 7):
                m_cfl = 0.3698835335671902
            elif (self.nb_stages == 8):
                m_cfl = 0.4403182119131088
            elif (self.nb_stages == 9):
                m_cfl = 0.5127951130270958
            elif (self.nb_stages == 10):
                m_cfl = 0.5859489552676678
            elif (self.nb_stages == 11):
                m_cfl = 0.6599807459861040
            elif (self.nb_stages == 12):
                m_cfl = 0.7341530174016953
            elif (self.nb_stages == 13):
                m_cfl = 0.8082503871992230
            elif (self.nb_stages == 14):
                m_cfl = 0.8820748422294855
            elif (self.nb_stages == 15):
                m_cfl = 0.9557751473039389
            elif (self.nb_stages == 16):
                m_cfl = 1.0291188955307007
            elif (self.nb_stages == 17):
                m_cfl = 1.1015691235661507
            elif (self.nb_stages == 18):
                m_cfl = 1.1741869524121284
            elif (self.nb_stages == 19):
                m_cfl = 1.2464623991400003
            elif (self.nb_stages == 20):
                m_cfl = 1.3186175376176834
            else:
                print "no m_cfl number defined for ",self.nb_stages," rk stages"
        elif (self.order == 5):
            if (self.nb_stages == 7):
                m_cfl = 0.2321297861635685
            elif (self.nb_stages == 8):
                m_cfl = 0.2740361541509628
            elif (self.nb_stages == 9):
                m_cfl = 0.3175020497292280
            elif (self.nb_stages == 10):
                m_cfl = 0.3616421483457088
            elif (self.nb_stages == 11):
                m_cfl = 0.4063922166824341
            elif (self.nb_stages == 12):
                m_cfl = 0.4514714889228344
            elif (self.nb_stages == 13):
                m_cfl = 0.4975992953404784
            elif (self.nb_stages == 14):
                m_cfl = 0.5458587780594826
            elif (self.nb_stages == 15):
                m_cfl = 0.5958906607702374
            elif (self.nb_stages == 16):
                m_cfl = 0.6462623924016953
            elif (self.nb_stages == 17):
                m_cfl = 0.6967312470078468
            elif (self.nb_stages == 18):
                m_cfl = 0.7464754488319159
            elif (self.nb_stages == 19):
                m_cfl = 0.7955208839848638
            elif (self.nb_stages == 20):
                m_cfl = 0.8438967168331146
            elif (self.nb_stages == 21):
                m_cfl = 0.8920353557914495
            elif (self.nb_stages == 22):
                m_cfl = 0.9395490400493145
            else:
                print "no cfl number defined for ",self.nb_stages," rk stages"
        elif (self.order == 6):
            if (self.nb_stages == 7):
                m_cfl = 0.1538090000000000
            elif (self.nb_stages == 8):
                m_cfl = 0.2004983276128769
            elif (self.nb_stages == 9):
                m_cfl = 0.2385877724736929
            elif (self.nb_stages == 10):
                m_cfl = 0.2732296474277973
            elif (self.nb_stages == 11):
                m_cfl = 0.3076814347878098
            elif (self.nb_stages == 12):
                m_cfl = 0.3423537500202656
            elif (self.nb_stages == 13):
                m_cfl = 0.3773931460455060
            elif (self.nb_stages == 14):
                m_cfl = 0.4127601720392704
            elif (self.nb_stages == 15):
                m_cfl = 0.4482912132516503
            elif (self.nb_stages == 16):
                m_cfl = 0.4838886857032776
            elif (self.nb_stages == 17):
                m_cfl = 0.5195324262604117
            elif (self.nb_stages == 18):
                m_cfl = 0.5552181880921125
            elif (self.nb_stages == 19):
                m_cfl = 0.5909251142293215
            elif (self.nb_stages == 20):
                m_cfl = 0.5859374068677425
            else:
                print "no cfl number defined for ",self.nb_stages," rk stages"
        
        else:
            print "order ",self.order," not supported"
        return m_cfl

# function that parses the file, and returns a valid CF3_ERK_Coeffs object
def parse(filepath):
    mydata = open(filepath).read()
    #print "Parsing file "+filepath
    regex = re.compile(r"#stage.*order.*\n(?P<nb_stages>[0-9]+?)\s+(?P<order>[0-9]+?)\s+A",re.DOTALL)
    m = re.match(regex, mydata).groupdict()
    
    class parsed_strings:
        order = m['order']
        #print "order = ",order
        nb_stages = m['nb_stages']
        #print "nb_stages = ",nb_stages
        A = re.compile(r".*\nA\n(.*?)\n\n",re.DOTALL).match(mydata).group(1)
        #print "A=\n",A
        b = re.compile(r".*\nb\n(.*?)\n\n",re.DOTALL).match(mydata).group(1)
        #print "b=\n",b
        c = re.compile(r".*\nc\n(.*?)\n\n",re.DOTALL).match(mydata).group(1)
        #print "c=\n",c
        alpha = re.compile(r".*\nalpha\n(.*?)\n\n",re.DOTALL).match(mydata).group(1)
        #print "alpha=\n",alpha
        beta = re.compile(r".*\nbeta\n(.*?)\n\n",re.DOTALL).match(mydata).group(1)
        #print "beta=\n",beta
        gamma1 = re.compile(r".*\ngamma1\n(.*?)\n\n",re.DOTALL).match(mydata).group(1)
        #print "gamma1=\n",gamma1
        gamma2 = re.compile(r".*\ngamma2\n(.*?)\n\n",re.DOTALL).match(mydata).group(1)
        #print "gamma2=\n",gamma2
        gamma3 = re.compile(r".*\ngamma3\n(.*?)\n\n",re.DOTALL).match(mydata).group(1)
        #print "gamma3=\n",gamma3
        delta = re.compile(r".*\ndelta\n(.*?)\n\n",re.DOTALL).match(mydata).group(1)
        #print "delta=\n",delta
    
    class parsed_lists:
        order = int(parsed_strings.order)
        #print "order = ",order

        nb_stages = int(parsed_strings.nb_stages)
        #print "nb_stages = ",nb_stages

        A = parsed_strings.A.split("\n")
        for k in range(nb_stages):
            A[k] = [float(s) for s in A[k].split()]
        #print "A=\n",A

        b = [float(s) for s in parsed_strings.b.split()]
        #print "b=\n",b

        c = [float(s) for s in parsed_strings.c.split()]
        #print "c=\n",c

        alpha = parsed_strings.alpha.split("\n")
        for k in range(nb_stages+1):
            alpha[k] = [float(s) for s in alpha[k].split()]
        #print "alpha = \n",alpha

        beta = parsed_strings.beta.split("\n")
        for k in range(nb_stages+1):
            beta[k] = [float(s) for s in beta[k].split()]
        #print "beta = \n",beta
        
        gamma1 = [float(s) for s in parsed_strings.gamma1.split()]
        #print "gamma1 = \n",gamma1

        gamma2 = [float(s) for s in parsed_strings.gamma2.split()]
        #print "gamma2 = \n",gamma2
        
        gamma3 = [float(s) for s in parsed_strings.gamma3.split()]
        #print "gamma3 = \n",gamma3
        
        delta = [float(s) for s in parsed_strings.delta.split()]
        #print "delta = \n",delta
        
    #------------------------------------------------------------
    # coolfluid needs only c, beta, gamma1, gamma2, gamma3, delta
    # from beta it only needs beta(i,i+1)
    # from gamma's it doesn't need the first element
    erk_coeffs = CF3_ERK_Coeffs()
    erk_coeffs.order = parsed_lists.order;
    erk_coeffs.nb_stages = parsed_lists.nb_stages;
    erk_coeffs.c = parsed_lists.c;
    erk_coeffs.beta = [0]*parsed_lists.nb_stages;
    for k in range(parsed_lists.nb_stages):
        erk_coeffs.beta[k] = parsed_lists.beta[k+1][k];
    erk_coeffs.gamma1 = parsed_lists.gamma1[1:len(parsed_lists.gamma1)];
    erk_coeffs.gamma2 = parsed_lists.gamma2[1:len(parsed_lists.gamma2)];
    erk_coeffs.gamma3 = parsed_lists.gamma3[1:len(parsed_lists.gamma3)];
    erk_coeffs.delta  = parsed_lists.delta;
    return erk_coeffs

def get_coeffs(order,stages) :
    this_dir=os.path.dirname(__file__)
    if (this_dir == ""):
        this_dir = "."
    rk = parse(this_dir+'/erk-x'+str(order)+'/3sstar/ERK-'+str(order)+'-'+str(stages)+'.txt')
    return rk

##########################################################################
# STANDALONE EXECUTION
##########################################################################

if __name__ == "__main__":
    import sys
    erk_coeffs = parse(sys.argv[1])
    print erk_coeffs
    print erk_coeffs.tex()
