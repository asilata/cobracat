from sage.combinat.root_system.weyl_group import WeylGroup
def reduced_word(root,ct):
    r"""
    Return the reduced word for a given root.
    
    Every root 'alpha' can be written as a word 'w' acting on a simple root alpha_i, i.e. alpha = w(alpha_i), for some i. 
    Here roots are written in terms of vectors, where alpha_i = e_i in R^n, n=rank(ct).
    
    Raises a type error if the inputted vector is not a root of the associated WeylGroup of the cartan type
    
    INPUT:
    
    root -- root vector
    ct -- Cartan type
    
    OUTPUT:
    
    reduced expression for the root and a vector corresponding to the reduced expression
    
    EXAMPLES:
    
        sage: reduced_word_new((1,1,1),A3)
        ('s3*s2*alpha1', [2, 1,0])
        
        sage: reduced_word_new((1,0),A2)
        ('1*alpha1', [0])
    
    """
    W=WeylGroup(ct,prefix='s',implementation='permutation')
    alpha=W.roots()
    is_root=False
    for x in alpha:
        if str(root)==str(x):
            is_root=True
            break
    if is_root==False:
        raise TypeError("Inputed vector is not a root of the root system of the associated Weyl Group.")
    max_len=float('inf')
    for j in W.index_set():
        for w in W:
            if str(w.action(alpha[j-1],side="left"))==str(root):
                if w.length()<max_len:
                    max_len=w.length()
                    word=str(w)+'*alpha'+str(j)
                    red_word=w._reduced_word
                    red_word.append(j-1)
    return word,red_word

def central_charge(charge,ct):
    r"""
    Assigns a central charge to each root of the root system, given the values of central charge for the simple roots.
    
    INPUT:
    
    charge -- a list (of length = no of simple roots) containing the central charge of the simple roots (basis elements)
    ct -- Cartan Type
    
    OUTPUT:
    
    returns a dictionary with keys:values as root:central charge respectively.
    
    EXAMPLES:
        sage: central_charge([1,1],'A2')
        {(1, 0): 1, (0, 1): 1, (1, 1): 2, (-1, 0): -1, (0, -1): -1, (-1, -1): -2}
        
    """
    W = WeylGroup(ct,prefix="s",implementation="permutation")
    if len(charge)<len(W.index_set()):
        raise TypeError("Inputted central charge hasn't been defined for all basis vectors")
    if len(charge)>len(W.index_set()):
        raise TypeError("Inputted central charge should only be defined for the basis vectors")
    alpha=W.roots()
    cent_charge={} # initialize a dictionary
    for x in alpha:
        cent_charge[str(x)] = sum(charge[j] * x[j] for j in range(len(W.index_set())))
    return cent_charge

def root_sequence(root,ct):
    r"""
    Given a root, it calculates the root sequence of the associated word.
    
    INPUT:
    
    root -- root vector
    ct -- Cartan type
    
    OUTPUT:
    
    the root sequence as a list of roots
    
    EXAMPLES:
    
        sage:root_sequence((1,1,1),'A3')
        [(0, 0, 1), (0, 1, 1), (1, 1, 0)]
    
    """
    word,red_word=reduced_word(root,ct)
    W = WeylGroup(ct,prefix="s",implementation="permutation")
    alpha_simp=W.simple_roots()
    simp_refl=W.gens()
    id=W.list()[0]
    alpha_seq=[]
    refl_seq=[id]
    for x in red_word[:-1]:
        alpha_seq.append(alpha_simp[x+1])
        refl_seq.append(refl_seq[-1]*simp_refl[x])
    alpha_seq.append(alpha_simp[red_word[-1]+1])
    root_seq=[]
    for i in range(len(red_word)):
        root_seq.append(refl_seq[i].action(alpha_seq[i]))
    return root_seq

def braid_lift(root,charge,ct):
    r"""
    Given a root, it gives the associated braid lift to the root.
    This code currently only works for Cartan Types A_n
    
    INPUT:
    
    root -- root vector
    charge -- a list (of length = no of simple roots) containing the central charge of the simple roots (basis elements)
    ct -- Cartan Type of Type A_n
    
    OUTPUT:
    
    a dictionary where the key:value pari correspond to the i^th generator of the braid group and its respective power
    
    EXAMPLES:
    
        sage:braid_lift((1,1,0),[I,3+I,1-I],'A3')
        {1: a1^-1*a0*a1}
        
        sage:braid_lift((1,1,1,0),[I,3+I,1-I,-2],'A4')
        {1: a2^-1*a1^-1*a0*a1^-1*a2}
        
        sage:braid_lift((1,1,1),[1,1,1],'A3')
        {1: a2*a1^-1*a0*a1^-1*a2,
         2: a2^-1*a1*a0*a1*a2,
         3: a2*a1*a0*a1*a2,
         4: a2^-1*a1^-1*a0*a1^-1*a2}
         
    """
    W = WeylGroup(ct,prefix="s",implementation="permutation")
    simp_refl=W.gens()
    word,red_word=reduced_word(root,ct)
    cent_charge=central_charge(charge,ct)
    root_seq=root_sequence(root,ct)
    
    tracker=[0]*len(root_seq)
    
    for i, x in enumerate(root_seq):
        arg_x = arg(cent_charge[str(x)])
        arg_root = arg(cent_charge[str(root)])
        
        if arg_x < arg_root:
            tracker[i] = -1
        elif arg_x > arg_root:
            tracker[i] = 1
        else:
            tracker[i] = 0
            
    tracker[-1]=1
    
    B=BraidGroup(int(ct[1])+1,'a')
    sigma=B.gens()
    perm = [sigma[i] for i in red_word]
    counter = sum(1 for x in tracker if x == 0)
    
    dictor={}
    list=['']*len(perm)
    for j in range(1,2**counter+1):
        dictor[j]=list[:]
        
    k=-1
    for i in range(len(perm)):
        if tracker[i]==0:
            for j in range(1,2**counter+1):
                if bin(j)[k]=='b':
                    (dictor[j])[i]=0
                else:
                    (dictor[j])[i]=int(bin(j)[k])
        k=k-1
        if tracker[i]==1: 
            for j in range(1,2**counter+1):
                (dictor[j])[i]=1
        if tracker[i]==-1: 
            for j in range(1,2**counter+1):
                (dictor[j])[i]=-1
    for j in range(1, 2 ** counter + 1):
        for i in range(len(perm)):
            dictor[j][i] = -1 if dictor[j][i] == 0 else dictor[j][i]
            
    perm_copy=perm[:]
    perm_copy.reverse()
    for x in perm_copy[1:]:
        perm.append(x)
    for j in range(1,2**counter+1):
        c=dictor[j][:]
        c.reverse
        for x in c[1:]:
            dictor[j].append(x)
            
    output={}
    for j in range(1, 2 ** counter + 1):
        result = perm[0] ** -1 * perm[0]
        for i, perm_elem in enumerate(perm):
            result *= perm_elem ** dictor[j][i]
        output[j] = result
    return output

def stab_obj(root,charge,ct,minimize=False):
    r"""
    Given a root, it gives the associated braid lift to the root.
    This code currently only works for Cartan Types A_n
    
    INPUT:
    
    root -- root vector
    charge -- a list (of length = no of simple roots) containing the central charge of the simple roots (basis elements)
    ct -- Cartan Type
    minimize -- takes T/F values. determines whether a not the complex is to be minimized using Gaussian elimination.
    
    OUTPUT:
    
    a dictionary where the values are the stable objects under the specified stability conditions
    
    NOTE: there may be more than one stable object since we are considering a generic stability condition
    
    EXAMPLES:
    
        sage:stab_obj((1,1,0),[I,3+I,1-I],'A3',False)
        {1: [0]: P1<0>   P1<0>+P1<2> :[1]}
        
        sage:stab_obj((1,1,0),[I,3+I,1-I],'A3',True)
        {1: [1]: P1<2> :[1]}
        
        sage:stab_obj((1,1,1,0),[I,3+I,1-I,-2],'A4',False)
        {1: [0]: P1<0>   P1<0>+P1<2>   P1<2>+P1<4> :[2]}
        
        sage:stab_obj((1,1,1,0),[I,3+I,1-I,-2],'A4',True)
        {1: [2]: P1<4> :[2]}
        
        sage:stab_obj((1,1,1),[1,1,1],'A3',False)
        {1: [-1]: 0   P1<2>+P1<0>+P1<0>   P1<0>+P1<2> :[1],
         2: [-1]: P1<0>+P1<-2>   P1<0>+P1<-2>+P1<0> :[0],
         3: [-2]: P1<-2>+P1<-4>   P1<0>+P1<-2>   P1<0> :[0],
         4: [0]: P1<0>   P1<0>+P1<2>   P1<2>+P1<4> :[2]}
         
         sage:stab_obj((1,1,1),[1,1,1],'A3',True)
         {1: [0]: P1<0> :[0],
          2: [0]: P1<0> :[0],
          3: [-2]: P1<-4> :[-2],
          4: [2]: P1<4> :[2]}
         
    """
    W = WeylGroup(ct,prefix="s",implementation="permutation")
    simp_refl=W.gens()
    word,red_word=reduced_word(root,ct)
    cent_charge=central_charge(charge,ct)
    root_seq=root_sequence(root,ct)
    p,s=braid_actions_setup(ct,QQ,minimize)
    
    tracker=[0]*len(root_seq)
    
    for i, x in enumerate(root_seq):
        arg_x = arg(cent_charge[str(x)])
        arg_root = arg(cent_charge[str(root)])
        
        if arg_x < arg_root:
            tracker[i] = -1
        elif arg_x > arg_root:
            tracker[i] = 1
        else:
            tracker[i] = 0
            
    tracker[-1]=1
    counter = sum(1 for x in tracker if x == 0)
    leng=len(red_word)
    
    dictor={}
    list=['']*leng
    for j in range(1,2**counter+1):
        dictor[j]=list[:]
        
    k=-1
    for i in range(leng):
        if tracker[i]==0:
            for j in range(1,2**counter+1):
                if bin(j)[k]=='b':
                    (dictor[j])[i]=0
                else:
                    (dictor[j])[i]=int(bin(j)[k])
        k=k-1
        if tracker[i]==1: 
            for j in range(1,2**counter+1):
                (dictor[j])[i]=1
        if tracker[i]==-1: 
            for j in range(1,2**counter+1):
                (dictor[j])[i]=-1
    for j in range(1, 2 ** counter + 1):
        for i in range(leng):
            dictor[j][i] = -1 if dictor[j][i] == 0 else dictor[j][i]
            
    obj=p[red_word[-1]+1];obj
    leng=len(dictor[1])
    stab_obj={}
    for j in dictor.keys():
        stab_obj[j]=obj
    for i in range(leng-1):
        for j in dictor.keys():
            stab_obj[j]=s[dictor[j][-i-2]](stab_obj[j])
    stab_obj
    return stab_obj
